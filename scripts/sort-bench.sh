#!/bin/bash
#
# Sorting benchmarks; Run from the scripts/ folder.
#
BIN=../src/bench_sort
NAME=sortbench-mway-buf1.25MB
OUTFILE=../measurements/$NAME.txt
LOGFILE=../measurements/$NAME.log

BUFNTUPLES=163840 #1.25MiB

# With AVX
echo "# Sorting with different avxsort routines (unaligned, mway, aligned); data (pow2, not-pow2)" >> $OUTFILE
for ALGO in 0 1 2 #0->avxsort() 1->avxsortmultiway() 2->avxsort_aligned()
do
    for ISPOW2 in 1 0 # numtuples-pow2 not-pow2
    do
		for NTUPLES in 1 2 4 8 16 32 64 128 256
		do
	    	for REP in 1 2 3 # Repeat 3 times
	    	do
				EXE="$BIN $NTUPLES $ALGO $ISPOW2" #$BUFNTUPLES"
				echo "========================================" >> $LOGFILE
				echo $EXE >> $LOGFILE
				$EXE 2>> $LOGFILE >> $OUTFILE
	    	done
		done
    done
done
