#!/bin/bash
#
# This is an experiment which mimics the input size scale experiment
# from the original Sort vs. Hash paper by Graefe et. al.
# Tables of equal size are scaled up to 2Billion tuples.
#
# Binary executable options:
BIN=../src/sortmergejoins
RESULTS=../measurements/
MINRELSIZE=33554432

# Just scalar with 16-Byte tuples; CODE MUST BE CONFIGURED & COMPILED WITH --enable-key8B
#
# OUTPUTFILE=scalesize-scalar
# JUSTSCALAR=YES
#
# With avx:
OUTPUTFILE=scalesize
JUSTSCALAR=NO
#

OUTFILE="$RESULTS/$OUTPUTFILE.txt"
LOGFILE="$RESULTS/$OUTPUTFILE.log"

# Algorithms implemented with AVX
echo "========== Input size scale experiment with equal sized tables (with AVX) ==========" >> $OUTFILE
if [ "$JUSTSCALAR" == "NO" ] 
then
    echo "# Input size scale experiment with equal sized tables" >> $OUTFILE
    echo "#ALGO NTHREADS NUMR NUMS RUNNO PARTCYC SORTCYC MERGE1CYC MERGERESTCYC MJOINCYC NTUP USECS TPUT" >> $OUTFILE
    for ALGO in m-pass m-way mpsm
    do
		for NTHREADS in 64 32 16 8 4 2 1
		do
		    for SM in 1 2 4 8 16 24 32 40 48 56 60 #64
		    do
				RTUPLES=$(($MINRELSIZE*$SM))
				STUPLES=$(($MINRELSIZE*$SM))
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
echo "========== Input size scale experiment with equal sized tables (without AVX) ==========" >> $OUTFILE
echo "#ALGO NTHREADS NUMR NUMS RUNNO PARTCYC SORTCYC MERGE1CYC MERGERESTCYC MJOINCYC NTUP USECS TPUT" >> $OUTFILE
for ALGO in m-pass m-way mpsm
do
    for NTHREADS in 64 32 16 8 4 2 1
    do
		for SM in 1 2 4 8 16 24 32 40 48 56 60 #64
		do
		    RTUPLES=$(($MINRELSIZE*$SM))
		    STUPLES=$(($MINRELSIZE*$SM))
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
