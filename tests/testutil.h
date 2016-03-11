#ifndef TESTUTIL_H
#define TESTUTIL_H

#include "types.h"

int
is_sorted_int64(int64_t * items, uint64_t nitems);

int
is_sorted_int32(int32_t * items, uint64_t nitems);

int
is_sorted_tuples(tuple_t * items, uint64_t nitems);

int
is_sorted_tuples_noassert(tuple_t * items, uint64_t nitems);

int
is_array_equal(int64_t *arr1, int64_t *arr2, uint64_t sz1, uint64_t sz2);

int64_t *
generate_rand_int64(int num);

int32_t *
generate_rand_int32(int num);

tuple_t *
generate_rand_tuples(int num);

int64_t *
generate_rand_ordered_int64(int num);

int32_t *
generate_rand_ordered_int32(int num);

tuple_t *
generate_rand_ordered_tuples(int num);

/**
 * Cache-line aligned memory allocation with posix_memalign()
 *
 * @param size
 */
void *
malloc_aligned(size_t size);

#endif
