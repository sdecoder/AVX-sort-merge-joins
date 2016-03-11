#!/bin/bash
# Run Sort merge joins with different number of threads
#
# Binary executable options:
BIN=../src/sortmergejoins
RESULTS=../measurements/
#
# Configure the data set from following:
#
# Configuration for Workload B (Kim et al., ORACLE/Intel Dataset):
#NTUPLES=128000000
#OUTPUTFILE=tput-WorkloadB
#
#
# Configuration for Workload A (Albutiu et al., TUM Dataset):
NTUPLES=1600000000
OUTPUTFILE=tput-WorkloadA
#

OUTFILE="$RESULTS/$OUTPUTFILE.txt"
LOGFILE="$RESULTS/$OUTPUTFILE.log"

# Algorithms implemented with AVX
echo "========== Algorithms implemented with AVX ==========" >> $OUTFILE
echo "#ALGO NTHREADS RUNNO PARTCYC SORTCYC MERGE1CYC MERGERESTCYC MJOINCYC NUMTUP USECS TPUT" >> $OUTFILE
for ALGO in m-pass m-way mpsm
do
    for NTHREADS in 64 32 16 8 4 2 1
    do
		for REP in 1 2 3 # Repeat 3 times
		do
	    	EXE="$BIN --algo=$ALGO --nthreads=$NTHREADS --r-size=$NTUPLES --s-size=$NTUPLES"
	    	echo "========================================" >> $LOGFILE
	    	echo $EXE >> $LOGFILE
	    	echo -n "$ALGO $NTHREADS $REP " >> $OUTFILE
	    	$EXE >> $LOGFILE 2>> $OUTFILE
	    	echo "" >> $OUTFILE
		done
    done
done



# Scalar algorithms without AVX
echo "========== Scalar algorithms without AVX ==========" >> $OUTFILE
echo "#ALGO NTHREADS RUNNO PARTCYC SORTCYC MERGE1CYC MERGERESTCYC MJOINCYC NUMTUP USECS TPUT" >> $OUTFILE
for ALGO in m-pass m-way mpsm
do
    for NTHREADS in 64 32 16 8 4 2 1
    do
		for REP in 1 2 3 # Repeat 3 times
		do
		    EXE="$BIN --algo=$ALGO --nthreads=$NTHREADS --r-size=$NTUPLES --s-size=$NTUPLES --scalarmerge --scalarsort"
		    echo "========================================" >> $LOGFILE
		    echo $EXE >> $LOGFILE
		    echo -n "$ALGO-scalar $NTHREADS $REP " >> $OUTFILE
		    $EXE >> $LOGFILE 2>> $OUTFILE
		    echo "" >> $OUTFILE
		done
    done
done
