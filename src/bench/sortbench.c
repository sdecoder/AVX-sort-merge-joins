/**
 * @file    sortbench.c
 * @author  Cagri Balkesen <cagri.balkesen@inf.ethz.ch>
 * @date    Tue Dec 18 23:20:10 2012
 * @version $Id $
 *
 * @brief   A microbenchmark for single-threaded sorting routines.
 *
 * (c) 2012-2014, ETH Zurich, Systems Group
 *
 * \ingroup MicroBenchmarks
 */
#include <stdio.h>
#include <sys/time.h>           /* gettimeofday */
#include <stdlib.h>             /* qsort */
#include <string.h>             /* memcpy */

#if defined(__cplusplus)
#include <algorithm>            /* sort() */
#include <iostream>
 using namespace std;
#endif

#include "testutil.h" /* is_sorted_tuples() */
#include "affinity.h"
#include "types.h" /* relation_t, tuple_t */
#include "generator.h"
#include "avxsort.h"
#include "avxsort_multiway.h"

#ifdef PERF_COUNTERS
#include "perf_counters.h"      /* PCM_x */
#endif

/** A size greater than the L3 cache for cooling the L3 cache */
#define MAX_L3_CACHE_SIZE (24*1024*1024)

#ifndef L2_CACHE_SIZE
#define L2_CACHE_SIZE (256*1024)
#endif

/** Number of tuples that can fit into L2 cache divided by 2 */
#ifndef BLOCKSIZE
#define BLOCKSIZE (L2_CACHE_SIZE / (2 * sizeof(tuple_t)))
#endif

#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif

void
print_timing2(uint64_t numtuples, struct timeval * start, struct timeval * end,
             FILE * out)
{
    double diff_usec = (((*end).tv_sec*1000000L + (*end).tv_usec)
                        - ((*start).tv_sec*1000000L+(*start).tv_usec));

    fprintf(out, "NUM-TUPLES = %lld TOTAL-TIME-USECS = %.4lf ", numtuples, diff_usec);
    fprintf(out, "TUPLES-PER-SECOND = ");
    fflush(out);

    fprintf(out, "%.4lf ", (numtuples/(diff_usec/1000000L)));
    fflush(out);
}

static inline int __attribute__((always_inline))
tuplekeycmp(const void * k1, const void * k2)
{
    int val = ((tuple_t *)k1)->key - ((tuple_t *)k2)->key;
    return val;
}

/**
 * @addtogroup MicroBenchmarks
 * Microbenchmarks for single-thread AVX sort routines.
 * @{
 */

/**
 * A microbenchmark for different single-threaded AVX sort routines.
 * @param ntuples
 * @param what
 * @return
 */
int
bench_avxsort(uint32_t ntuples,
              int what /* 0: avxsort() , 1: avxsort_multiwaymerge(),
                          2: avxsort_aligned() */)
{
    struct timeval start, end;
    relation_t rel, outrel, relcpy, outrel2;
    tuple_t * dummybuffer = (tuple_t *) malloc_aligned(MAX_L3_CACHE_SIZE);

    outrel.num_tuples = ntuples;
    posix_memalign((void**)&outrel.tuples, CACHE_LINE_SIZE, ntuples * sizeof(tuple_t));
    outrel2.num_tuples = ntuples;
    posix_memalign((void**)&outrel2.tuples, CACHE_LINE_SIZE, ntuples * sizeof(tuple_t));

    fprintf(stderr,
            "[INFO ] Creating relation with size = %.3lf MiB, #tuples = %d : ",
            (double) sizeof(tuple_t) * ntuples/1024.0/1024.0,
            ntuples);
    seed_generator(12345);
    posix_memalign((void**)&rel.tuples, CACHE_LINE_SIZE, ntuples * sizeof(tuple_t));
    create_relation_pk(&rel, ntuples);
    /* parallel_create_relation(&rel, ntuples, 1, ntuples); */
    fprintf(stderr, " OK\n");
    fflush(stderr);

    /* What about qsort() and std::sort() ? */
    relcpy.num_tuples = ntuples;
    posix_memalign((void**)&relcpy.tuples, CACHE_LINE_SIZE, ntuples * sizeof(tuple_t));

#if 0
    memcpy(relcpy.tuples, rel.tuples, ntuples * sizeof(tuple_t));
    gettimeofday(&start, NULL);
    qsort (relcpy.tuples, ntuples, sizeof (int64_t), tuplekeycmp);
    gettimeofday(&end, NULL);
    fprintf(stdout, "[INFO ] C stdlib qsort() ==> ");
    print_timing2(ntuples, &start, &end, stdout);
#endif

#if defined(__cplusplus)
    /* C++ std::sort() */
    memcpy(relcpy.tuples, rel.tuples, ntuples * sizeof(tuple_t));
    gettimeofday(&start, NULL);
    std::sort((int64_t *) relcpy.tuples, (int64_t *)(relcpy.tuples + ntuples));
    gettimeofday(&end, NULL);
    fprintf(stdout, "[INFO ] C++ std::sort()  ==> ");
    print_timing2(ntuples, &start, &end, stdout);
    fprintf(stderr, ", \n");
#endif

    free(relcpy.tuples);

    /* cool-down the caches */
    intkey_t garbage = 0xdeadbeef;
    for(unsigned int i = 0; i < (MAX_L3_CACHE_SIZE/sizeof(tuple_t)); i++){
        garbage += dummybuffer[i].key;
        dummybuffer[i].key = garbage;
    }

    /* AVX sort */
    fprintf(stderr, "[INFO ] Running single core AVX sort...\n");
#ifdef PERF_COUNTERS
    PCM_initPerformanceMonitor("pcm.cfg", NULL);
    PCM_start();
#endif

    gettimeofday(&start, NULL);

    if(what == 0 || what == 2){
        /* chooses aligned version depending on input alignment */
        avxsort_tuples(&rel.tuples, &outrel.tuples, ntuples);
    }
    else if(what == 1){
        /* if size is larger than L3-size; then uses multi-way merging */
        avxsortmultiway_tuples(&rel.tuples, &outrel.tuples, ntuples);
    }

    gettimeofday(&end, NULL);

#ifdef PERF_COUNTERS
    PCM_stop();
    PCM_log("========== Profiling results ==========\n");
    PCM_printResults();
    PCM_cleanup();
#endif

    char impl[256];

    if(what == 0)
        sprintf(impl, "%s", "avxsort()");
    else if(what == 1)
        sprintf(impl, "%s", "avxsort_multiwaymerge()");
    else if(what == 2)
        sprintf(impl, "%s", "avxsort_aligned()");
    else
        sprintf(impl, "%s", "none()");

    fprintf(stdout, "[INFO ] %s ==> ", impl);
    print_timing2(ntuples, &start, &end, stdout);
    fprintf(stdout, "\n");

    int rv = 1;

    if(is_sorted_tuples(outrel.tuples, outrel.num_tuples)) {
        fprintf(stderr, "[INFO ] Output relation is now sorted. (garbage=%d)\n", garbage%4);
    }
    else {
        fprintf(stderr, "[ERROR] Output relation is not sorted.\n");
        rv = 0;
    }


    /* clean-up temporary space */
    free(rel.tuples);
    free(outrel.tuples);
    free(dummybuffer);

    return rv;
}

/** @} */

/** Benchmark sorting over numa spread data by using numactl [options] */
int
bench_sort_over_numa()
{
    struct timeval start, end;
    relation_t rel, relcpy, outrel;
    int32_t ntuples = 50000000;
    int rv = 1;

    outrel.num_tuples = ntuples;
    posix_memalign((void**)&outrel.tuples, CACHE_LINE_SIZE, ntuples * sizeof(tuple_t));

    fprintf(stdout,
            "[INFO ] Creating relation with size = %.3lf MiB, #tuples = %d : ",
            (double) sizeof(tuple_t) * ntuples/1024.0/1024.0,
            ntuples);
    fflush(stdout);
    seed_generator(12345);
    create_relation_pk(&rel, ntuples);
    fprintf(stdout, " OK\n");

    /* What about qsort() and std::sort() ? */
    relcpy.num_tuples = ntuples;
    posix_memalign((void**)&relcpy.tuples, 64, ntuples * sizeof(tuple_t));
    memcpy(relcpy.tuples, rel.tuples, ntuples * sizeof(tuple_t));
    gettimeofday(&start, NULL);
    qsort (relcpy.tuples, ntuples, sizeof (tuple_t), tuplekeycmp);
    gettimeofday(&end, NULL);
    fprintf(stdout, "[INFO ] C stdlib qsort() ==> ");
    print_timing2(ntuples, &start, &end, stdout);


#if defined(__cplusplus)
    /* C++ std::sort() */
    memcpy(relcpy.tuples, rel.tuples, ntuples * sizeof(tuple_t));
    printf("[INFO ] Press enter to start C++ sort()...\n");getchar();
    gettimeofday(&start, NULL);
    std::sort((int64_t *) relcpy.tuples, (int64_t *)(relcpy.tuples + ntuples));
    gettimeofday(&end, NULL);
    fprintf(stdout, "[INFO ] C++ std::sort()  ==> ");
    print_timing2(ntuples, &start, &end, stdout);
#endif
    free(relcpy.tuples);


    /* clean-up temporary space */
    free(rel.tuples);
    free(outrel.tuples);

    return rv;
}

/**
 * Benchmark Multi-Way Merge
 *
 * Arguments:
 * [1]: number of tuples in 2^20;
 * [2]: 0=>avxsort(), 1=>avxsort_multiwaymerge(), 2:
 * [3]: is number of tuples a power of 2 ? 0 or 1
 */

int gen_random_int()
{
  const int value = rand();
  int sign = rand();
  return sign > RAND_MAX / 2 ? value: -value; 
  
}

inline void SetPtr(int64_t& carrier, int ptr){
    int _ptr = ptr & 0x000FFFFF; 
    carrier =  carrier | _ptr;  
}

inline int GetPtr(int64_t& carrier){
    return (int) carrier & 0x000FFFFF;  
}

inline void SetKeyInt(int64_t& carrier, int64_t key){   
    if( key >= 0 )  carrier =  carrier | (key << 20);   
    else {
        carrier = carrier | (((int64_t)1) << 63);
        int64_t reg = - key; 
        carrier = carrier | (reg << 20);
    }
}


inline int GetKeyInt(int64_t& carrier){ 
    int64_t sig = carrier & (((int64_t) 1) << 63);
    if(sig < 0) {
        int reg_ = (int) ((carrier &  0x000FFFFFFFF00000 ) >> 20); 
        return -reg_;       
    }
    else {
        return (int) ((carrier &  0x000FFFFFFFF00000 ) >> 20);
    }
}

void _hybridsort(int64_t* data, int nitems){
  int64_t* iptr = new int64_t[nitems];
  const int start = 0;
  int bound_;   
  {//inline int partition_(float* array, int length){
    int astart_ = 0; 
    int aend_ = nitems -1; 
    while( astart_ < aend_){
      while ( data[astart_] < 0 && astart_ < nitems) ++astart_; 
      while ( data[aend_] > 0 && aend_ >= 0 ) --aend_; 
      if( astart_ < aend_){
        //swap 
        int tmp = data[astart_]; 
        data[astart_] = data[aend_];
        data[aend_]   = tmp; 
      }
    }  
      bound_ = aend_;
  }

  int64_t * neg_input; 
  int64_t * neg_output;  
  int64_t * pos_input;
  int64_t * pos_output;  

  if( bound_ > 0 ) {
 
    const int neg_part_length =  bound_ + 1;
    posix_memalign((void**)&neg_input, CACHE_LINE_SIZE, neg_part_length* sizeof(int64_t));
    posix_memalign((void**)&neg_output, CACHE_LINE_SIZE, neg_part_length* sizeof(int64_t));
    for (int i = 0; i < bound_ + 1; ++i)
    {
      //neg_input[i] = fdata[i];
      neg_input[i] = data[i] & 0x7FFFFFFFFFFFFFFF;//0;
    }
    avxsort_int64(&neg_input, &neg_output, neg_part_length);  
    for (int i = 0; i < neg_part_length; ++i){
      iptr[start + i] = GetPtr(neg_output[ neg_part_length - 1 - i]);
    }
    free(neg_input);
    free(neg_output);

    const int pos_part_length = (nitems -  bound_ - 1);
    posix_memalign((void**)&pos_input, CACHE_LINE_SIZE, pos_part_length* sizeof(tuple_t));
    posix_memalign((void**)&pos_output, CACHE_LINE_SIZE, pos_part_length* sizeof(tuple_t));
    memcpy(pos_input, data + (bound_ + 1), pos_part_length * sizeof(tuple_t));
    avxsort_int64(&pos_input, &pos_output, pos_part_length);   
    for (int i = 0; i < pos_part_length; ++i){
      iptr[start + neg_part_length + i] = GetPtr(pos_output[i]);
    }

    free(pos_input); 
    free(pos_output);
    
  } else {
    //LOG(INFO) << "[dbg] SortArray: all positive/1 negetive";
    int64_t* output;
    posix_memalign((void**)&output, CACHE_LINE_SIZE, nitems* sizeof(tuple_t));
    avxsort_int64(&data, &output, nitems);   
    for (int i = 0; i < nitems; ++i){   
     iptr[start + i] = GetPtr(output[i]);      
    }  
   // free(output);
  }

}


void bench_hybridsort(){
  cout<< "[dbg] benchmark hybrid AVX sort:" << endl;
  const int ntuples = 100000000; 
  struct timeval start, end;
  int64_t* inputdata;
  int64_t* outputdata; 
  posix_memalign((void**)&inputdata, CACHE_LINE_SIZE, ntuples * sizeof(int64_t));
  posix_memalign((void**)&outputdata, CACHE_LINE_SIZE, ntuples * sizeof(int64_t));
  for (int i = 0; i < ntuples; ++i)
  {
    inputdata[i] = 0; 
    SetPtr(inputdata[i], i);
    SetKeyInt(inputdata[i], gen_random_int());
  }

  puts("benchmark _hybridsort: ");
  gettimeofday(&start, NULL);
  _hybridsort(inputdata, ntuples); 
  gettimeofday(&end, NULL);
  long _time_consumed = (end.tv_sec-start.tv_sec) * 1000000 +(end.tv_usec - start.tv_usec); 
  cout<<"[dbg] time consumed using _hybridsort: " << _time_consumed << " ms" << endl;
  print_timing2(ntuples, &start, &end, stdout); puts("");

  puts("benchmark direct sort: ");
  gettimeofday(&start, NULL);
  avxsort_int64(&inputdata, &outputdata, ntuples); 
  gettimeofday(&end, NULL);
  _time_consumed = (end.tv_sec-start.tv_sec) * 1000000 +(end.tv_usec - start.tv_usec); 
  cout<<"[dbg] time consumed using direct sort: " << _time_consumed << " ms" << endl;
  print_timing2(ntuples, &start, &end, stdout); puts("");

}


void debug_pos_neg_sort(){
  const int ntuples = 100; 
  cout<< "[dbg] debug positive/negetive sort, array length: " << ntuples << endl;
  int64_t* inputdata;
  int64_t* outputdata;
  posix_memalign((void**)&inputdata, CACHE_LINE_SIZE, ntuples * sizeof(int64_t));
  posix_memalign((void**)&outputdata, CACHE_LINE_SIZE, ntuples * sizeof(int64_t));
  for (int i = 0; i < ntuples; ++i)
  {
    inputdata[i] = 0; 
    SetPtr(inputdata[i], i);
    SetKeyInt(inputdata[i], gen_random_int());
    //SetKeyInt(inputdata[i], -i);
    
  }

  avxsort_int64(&inputdata, &outputdata, ntuples); 
	cout<< "[dbg] dumping out sorted array:" << endl;
	for (int i = 0; i < ntuples; ++i)
  {
  	int ptr = GetPtr(outputdata[i]);
		int key = GetKeyInt(outputdata[i]);
  	cout << "\tINDEX: " << i << " PTR: " << ptr << " KEY: " << key << endl;
  }

  cout << "[dbg] verifying the array: " << endl;
  bool verify_result = true; 
	for (int i = 1; i < ntuples; ++i)
  {

  	int pre_ptr = GetPtr(outputdata[i-1]);
		int pre_key = GetKeyInt(outputdata[i-1]);

  	int cur_ptr = GetPtr(outputdata[i]);
		int cur_key = GetKeyInt(outputdata[i]);
		if(pre_key <= cur_key) continue; 
		else {
      verify_result = false;
			cout << "[dbg] Out of order at: " << endl; 
			cout << "\tINDEX: " << i -1  << " PTR: " << pre_ptr << " KEY: " << pre_key << endl;
			cout << "\tINDEX: " << i  << " PTR: " << cur_ptr << " KEY: " << cur_key << endl;
			break; 
		}
	
  }

  if( verify_result ) cout<< "[dbg] PASS" << endl;
  else cout << "[dbg] FAIL" << endl;

}


int
main(int argc, char ** argv)
{
    /* start initially on CPU-0 */
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(0, &set);
    if (sched_setaffinity(0, sizeof(set), &set) <0) {
        perror("sched_setaffinity");
    }
    if(argc != 4){
        printf("[INFO ] Usage: %s [numtuples in 2^20] [0=>avxsort(), 1=>avxsort_multiwaymerge(), 2=>aligned] [numtuples pow2?]\n",
               argv[0]);
        return 0;
    }

    int seed = 20121984;/*time(NULL);*/
    fprintf(stderr, "[INFO ] seed = %d\n", seed);
    srand(seed);

    int32_t ntuples = 1024*1024;
    int what = 0; /* 0: avxsort() , 1: avxsort_multiwaymerge(), 2: */

    ntuples *= atoi(argv[1]);
    what = atoi(argv[2]);
    if(atoi(argv[3]) == 0)
        ntuples +=  rand() % BLOCKSIZE;

    //bench_avxsort(ntuples, what);
    //bench_hybridsort();
    debug_pos_neg_sort();
    /* numa-bench */
#if 0
    fprintf(stderr, "[INFO ] Sorting over NUMA-spread data benchmark...\n");
    bench_sort_over_numa();
    fprintf(stderr, "[INFO ] ----------------------\n");
#endif

    return 0;
}
