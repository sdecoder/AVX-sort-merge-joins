#!/bin/bash
###### 
#
# Cardinality benchmarks varies the outer table size with various ratios to inner.
#    
######
# Binary executable options:
BIN=../src/sortmergejoins
RESULTS=../measurements/

##### Just scalar without AVX; CODE MUST BE CONFIGURED & COMPILED WITH --enable-key8B
#
# OUTPUTFILE=cardinality-scalar
# RTUPLES=1600000000
# JUSTSCALAR=YES
#
##### AVX-based implementations #####
#
OUTPUTFILE=cardinality-128M #-1600M
RTUPLES=134217728 #536870912 #500000000 #1600000000 
JUSTSCALAR=NO
#

OUTFILE="$RESULTS/$OUTPUTFILE.txt"
LOGFILE="$RESULTS/$OUTPUTFILE.log"

# Algorithms implemented with AVX
echo "========== Cardinality benchmarks: Algorithms implemented with AVX ==========" >> $OUTFILE
if [ "$JUSTSCALAR" == "NO" ] 
then
    echo "#ALGO NTHREADS NUMR NUMS RUNNO PARTCYC SORTCYC MERGE1CYC MERGERESTCYC MJOINCYC NUMTUP USECS TPUT" >> $OUTFILE
    for ALGO in m-pass m-way mpsm
    do
		for NTHREADS in 64 #32 16 8 4 2 1 # NOTE: Enable other nr. of threads as needed
		do
		    for SM in 1 2 4 8 16 # Outer to inner relation size ratio; NOTE: only up to 8 for R > 1B
		    do
				STUPLES=$(($RTUPLES*$SM))
				for REP in 1 2 3 # Repeat 3 times
				do
				    EXE="$BIN --algo=$ALGO --nthreads=$NTHREADS --r-size=$RTUPLES --s-size=$STUPLES"
			    	echo "========================================" >> $LOGFILE
			    	echo $EXE >> $LOGFILE
			    	echo -n "$ALGO $NTHREADS $RTUPLES $STUPLES $REP " >> $OUTFILE
				$EXE >> $LOGFILE 2>> $OUTFILE
			    	echo "" >> $OUTFILE
				done
		    done
		done
    done
fi


# Scalar algorithms without AVX
echo "========== Cardinality benchmarks: Scalar algorithms without AVX ==========" >> $OUTFILE
echo "#ALGO NTHREADS NUMR NUMS RUNNO PARTCYC SORTCYC MERGE1CYC MERGERESTCYC MJOINCYC NUMTUP USECS TPUT" >> $OUTFILE
for ALGO in m-pass m-way mpsm
do
    for NTHREADS in 64 #32 16 8 4 2 1 # NOTE: Enable other nr. of threads as needed
    do
		for SM in 1 2 4 8 16 # Outer to inner relation size ratio; NOTE: only up to 8 for R > 1B
		do
		    STUPLES=$(($RTUPLES*$SM))
		    for REP in 1 2 3 # Repeat 3 times
		    do
				EXE="$BIN --algo=$ALGO --nthreads=$NTHREADS --r-size=$RTUPLES --s-size=$STUPLES --scalarmerge --scalarsort"
		    	echo "========================================" >> $LOGFILE
		    	echo $EXE >> $LOGFILE
		    	echo -n "$ALGO $NTHREADS $RTUPLES $STUPLES $REP " >> $OUTFILE
			    $EXE >> $LOGFILE 2>> $OUTFILE
		    	echo "" >> $OUTFILE
		    done
		done
    done
done
