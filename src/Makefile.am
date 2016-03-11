AUTOMAKE_OPTIONS = foreign

INCLUDE  = -I. -I$(top_srcdir)/src/util/ -I$(top_srcdir)/src/datagen/ 
INCLUDE += -I$(top_srcdir)/src/merge/ -I$(top_srcdir)/src/avxsort/ -I$(top_srcdir)/src/scalarsort/ 
INCLUDE += -I$(top_srcdir)/src/partition/ -I$(top_srcdir)/src/joins/

#
# Note: Best way to globally set system cache parameters is editing the following:
#
GLOBAL_SYSTEM_PARAMS= #-DCACHE_LINE_SIZE=64 -DL2_CACHE_SIZE=262144 -DL3_CACHE_SIZE=20971520

DEFINES = $(GLOBAL_SYSTEM_PARAMS) -D_GNU_SOURCE
MYLDFLAGS =
MISC = 
MISC_SRC = 

if HAVE_NUMA
MYLDFLAGS += -lnuma
DEFINES += -DHAVE_NUMA
endif

if KEY_8B
DEFINES += -DKEY_8B
endif

if PERF_COUNTERS
DEFINES += -DPERF_COUNTERS
MYLDFLAGS += -lperf
MISC += -I$(top_srcdir)/lib/intel-pcm-2.35/ -L$(top_srcdir)/lib/intel-pcm-2.35/lib/
MISC_SRC += $(top_srcdir)/src/util/perf_counters.c $(top_srcdir)/src/util/perf_counters.h
LIBS += -lperf
endif

if DEBUG
DEFINES += -DDEBUG
else
DEFINES += -DNDEBUG
endif

if NOTIMING
DEFINES += -DNO_TIMING
endif

if JOIN_MATERIALIZE
DEFINES += -DJOIN_MATERIALIZE
endif

if PERSIST_RELATIONS
DEFINES += -DPERSIST_RELATIONS
endif

AM_CFLAGS  = -Wall -Wno-format -fno-strict-aliasing -std=c99 #-ansi 
if HAVE_AVX
DEFINES += -DHAVE_AVX 
AM_CFLAGS += -mavx $(DEFINES)
else
AM_CFLAGS += -msse4.2 $(DEFINES)
endif

# for compiling with Intel Compiler (icc)
#ICC_FLAGS =
#-opt_report

AM_LDFLAGS = #-pthread
CFLAGS = -O3 -g

bin_PROGRAMS = sortmergejoins bench_partitioning bench_merge bench_multiwaymerge bench_sort tputbench


### 
# Main program for the sort-merge join algorithms
###
sortmergejoins_LDFLAGS = $(AM_LDFLAGS) $(MYLDFLAGS) -lpthread #-pthread 
sortmergejoins_CFLAGS  = $(AM_CFLAGS) $(MISC) $(INCLUDE)
sortmergejoins_SOURCES = $(MISC_SRC) \
	$(top_srcdir)/config.h \
	$(top_srcdir)/src/types.h \
	$(top_srcdir)/src/util/affinity.h \
	$(top_srcdir)/src/util/cpu_mapping.h \
	$(top_srcdir)/src/util/cpu_mapping.c \
	$(top_srcdir)/src/util/numa_shuffle.h \
	$(top_srcdir)/src/util/numa_shuffle.c \
	$(top_srcdir)/src/util/memalloc.h \
	$(top_srcdir)/src/util/memalloc.c \
	$(top_srcdir)/src/util/barrier.h \
	$(top_srcdir)/src/util/rdtsc.h \
	$(top_srcdir)/src/avxsort/avxcommon.h \
	$(top_srcdir)/src/avxsort/avxsort_core.h \
	$(top_srcdir)/src/avxsort/avxsort.h \
	$(top_srcdir)/src/avxsort/avxsort.c \
	$(top_srcdir)/src/merge/merge.h \
	$(top_srcdir)/src/merge/merge.c \
	$(top_srcdir)/src/scalarsort/scalarsort.h \
	$(top_srcdir)/src/scalarsort/scalarsort.c \
	$(top_srcdir)/src/datagen/genzipf.h \
	$(top_srcdir)/src/datagen/genzipf.c \
	$(top_srcdir)/src/datagen/generator.h \
	$(top_srcdir)/src/datagen/generator.c \
	$(top_srcdir)/src/partition/partition.h \
	$(top_srcdir)/src/partition/partition.c \
	$(top_srcdir)/src/joins/joincommon.h \
	$(top_srcdir)/src/joins/joincommon.c \
	$(top_srcdir)/src/joins/sortmergejoin_multipass.h \
	$(top_srcdir)/src/joins/sortmergejoin_multipass.c \
	$(top_srcdir)/src/merge/scalar_multiwaymerge.h \
	$(top_srcdir)/src/merge/scalar_multiwaymerge.c \
	$(top_srcdir)/src/merge/avx_multiwaymerge.h \
	$(top_srcdir)/src/merge/avx_multiwaymerge.c \
	$(top_srcdir)/src/joins/sortmergejoin_multiway.h \
	$(top_srcdir)/src/joins/sortmergejoin_multiway.c \
	$(top_srcdir)/src/joins/sortmergejoin_mpsm.h \
	$(top_srcdir)/src/joins/sortmergejoin_mpsm.c \
	$(top_srcdir)/src/main.c

###
# Partitioning microbenchmark
###
bench_partitioning_LDFLAGS = $(AM_LDFLAGS) $(MYLDFLAGS)
bench_partitioning_CFLAGS  = $(AM_CFLAGS) $(MISC) $(INCLUDE)
bench_partitioning_SOURCES =	$(MISC_SRC) $(top_srcdir)/src/types.h	\
				$(top_srcdir)/src/util/rdtsc.h		\
				$(top_srcdir)/src/util/cpu_mapping.h	\
				$(top_srcdir)/src/util/cpu_mapping.c	\
				$(top_srcdir)/src/util/affinity.h	\
				$(top_srcdir)/src/util/lock.h		\
				$(top_srcdir)/src/util/memalloc.h	\
				$(top_srcdir)/src/util/memalloc.c	\
				$(top_srcdir)/src/util/barrier.h	\
				$(top_srcdir)/src/datagen/genzipf.h	\
				$(top_srcdir)/src/datagen/genzipf.c	\
				$(top_srcdir)/src/datagen/generator.h	\
				$(top_srcdir)/src/datagen/generator.c	\
				$(top_srcdir)/src/partition/partition.h	\
				$(top_srcdir)/src/partition/partition.c	\
				$(top_srcdir)/src/bench/partitioningbench.c			

###
# Merge microbenchmarks
###
bench_merge_LDFLAGS = $(AM_LDFLAGS) $(MYLDFLAGS)
bench_merge_CFLAGS  = $(AM_CFLAGS) $(MISC) $(INCLUDE) -I$(top_srcdir)/tests/
bench_merge_SOURCES =	$(MISC_SRC) $(top_srcdir)/src/util/affinity.h \
			$(top_srcdir)/src/avxsort/avxsort_core.h \
			$(top_srcdir)/tests/testutil.h \
			$(top_srcdir)/tests/testutil.c \
			$(top_srcdir)/src/bench/mergebench.c

###
# Multi-way merge microbenchmarks
###
bench_multiwaymerge_LDFLAGS = $(AM_LDFLAGS) $(MYLDFLAGS)
bench_multiwaymerge_CFLAGS  = $(AM_CFLAGS) $(MISC) $(INCLUDE) -I$(top_srcdir)/tests/
bench_multiwaymerge_SOURCES = $(MISC_SRC) $(top_srcdir)/src/util/affinity.h \
			  $(top_srcdir)/src/merge/avx_multiwaymerge.h \
			  $(top_srcdir)/src/merge/avx_multiwaymerge.c \
			  $(top_srcdir)/src/merge/scalar_multiwaymerge.h \
			  $(top_srcdir)/src/merge/scalar_multiwaymerge.c \
			  $(top_srcdir)/tests/testutil.h \
			  $(top_srcdir)/tests/testutil.c \
			  $(top_srcdir)/src/bench/multiwaymergebench.c

###
# Sorting microbenchmarks
### 
bench_sort_LDFLAGS = $(AM_LDFLAGS) $(MYLDFLAGS)
bench_sort_CFLAGS  = $(AM_CFLAGS) $(MISC) $(INCLUDE) -I$(top_srcdir)/tests/
bench_sort_SOURCES =	$(MISC_SRC) $(top_srcdir)/src/util/affinity.h	\
			$(top_srcdir)/src/util/cpu_mapping.h		\
			$(top_srcdir)/src/util/cpu_mapping.c		\
			$(top_srcdir)/tests/testutil.h			\
			$(top_srcdir)/tests/testutil.c			\
			$(top_srcdir)/src/avxsort/avxcommon.h		\
			$(top_srcdir)/src/avxsort/avxsort_core.h	\
			$(top_srcdir)/src/avxsort/avxsort.h		\
			$(top_srcdir)/src/avxsort/avxsort.c		\
			$(top_srcdir)/src/merge/merge.h			\
		        $(top_srcdir)/src/merge/merge.c			\
			$(top_srcdir)/src/merge/avx_multiwaymerge.h	\
			$(top_srcdir)/src/merge/avx_multiwaymerge.c	\
			$(top_srcdir)/src/avxsort/avxsort_multiway.h	\
			$(top_srcdir)/src/avxsort/avxsort_multiway.c	\
			$(top_srcdir)/src/datagen/genzipf.h		\
			$(top_srcdir)/src/datagen/genzipf.c		\
			$(top_srcdir)/src/datagen/generator.h		\
			$(top_srcdir)/src/datagen/generator.c		\
			$(top_srcdir)/src/bench/sortbench.c							


### 
# Throughput benchmarks for partitioning/merging and NUMA memcpy/read using m-way join
###
tputbench_LDFLAGS = $(AM_LDFLAGS) $(MYLDFLAGS) -lpthread #-pthread 
tputbench_CFLAGS  = $(AM_CFLAGS) $(MISC) $(INCLUDE)
tputbench_SOURCES = $(MISC_SRC) \
	$(top_srcdir)/config.h \
	$(top_srcdir)/src/types.h \
	$(top_srcdir)/src/util/affinity.h \
	$(top_srcdir)/src/util/cpu_mapping.h \
	$(top_srcdir)/src/util/cpu_mapping.c \
	$(top_srcdir)/src/util/barrier.h \
	$(top_srcdir)/src/util/memalloc.h \
	$(top_srcdir)/src/util/memalloc.c \
	$(top_srcdir)/src/util/rdtsc.h \
	$(top_srcdir)/src/avxsort/avxcommon.h \
	$(top_srcdir)/src/avxsort/avxsort_core.h \
	$(top_srcdir)/src/avxsort/avxsort.h \
	$(top_srcdir)/src/avxsort/avxsort.c \
	$(top_srcdir)/src/merge/merge.h \
	$(top_srcdir)/src/merge/merge.c \
	$(top_srcdir)/src/scalarsort/scalarsort.h \
	$(top_srcdir)/src/scalarsort/scalarsort.c \
	$(top_srcdir)/src/datagen/genzipf.h \
	$(top_srcdir)/src/datagen/genzipf.c \
	$(top_srcdir)/src/datagen/generator.h \
	$(top_srcdir)/src/datagen/generator.c \
	$(top_srcdir)/src/partition/partition.h \
	$(top_srcdir)/src/partition/partition.c \
	$(top_srcdir)/src/joins/joincommon.h \
	$(top_srcdir)/src/joins/joincommon.c \
	$(top_srcdir)/src/merge/scalar_multiwaymerge.h \
	$(top_srcdir)/src/merge/scalar_multiwaymerge.c \
	$(top_srcdir)/src/merge/avx_multiwaymerge.h \
	$(top_srcdir)/src/merge/avx_multiwaymerge.c \
	$(top_srcdir)/src/bench/tputbench.c
