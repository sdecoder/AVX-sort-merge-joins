/**
 * @file   numa_shuffle.h
 * @author Cagri Balkesen <cagri.balkesen@inf.ethz.ch>
 * @date   Wed May 23 16:43:30 2012
 *
 * @brief  Provides a machine specific implementation for NUMA-shuffling strategies.
 *
 * @note   The implementation needs to be customized for different hardware
 *         based on the NUMA topology.
 *
 * (c) 2014, ETH Zurich, Systems Group
 *
 */
#ifndef NUMA_SHUFFLE_H_
#define NUMA_SHUFFLE_H_

#include "types.h" /* enum numa_strategy_t */

/**
 * \ingroup numa
 *
 * Initialize the NUMA shuffling strategy with one of the following:
 *
 * RANDOM, RING, NEXT
 */
void
numa_shuffle_init(enum numa_strategy_t numastrategy, int nthreads);

/**
 * \ingroup numa
 *
 * Machine specific implementation of NUMA-shuffling strategies.
 *
 * @note The implementation needs to be customized for different hardware
 *       based on the NUMA topology.
 *
 * @param my_tid logical thread id of the calling thread.
 * @param i next thread index for shuffling (between 0 and nthreads)
 * @param nthreads number of total threads
 * @return the logical thread id of the destination thread for data shuffling
 */
int
get_numa_shuffle_strategy(int my_tid, int i, int nthreads);


#endif /* NUMA_SHUFFLE_H_ */
