#!/bin/bash
#
# Benchmark AVX & Scalar Multi-way Merge with different parameters
# Run from the scripts/ folder.
#
BIN=../src/bench_multiwaymerge
OUTPUTDIR=../measurements/

OUTFILE=$OUTPUTDIR/mergebench.txt
LOGFILE=$OUTPUTDIR/mergebench.log
L1SIZE=$((64*1024))
L3SIZE=$((8*1024*1024))
MINCHUNKSIZE=$((16*1024))
MAXNUMTUPLES=134217728

echo "#Benchmark of AVX Multi-way Merge with different parameters" >> $OUTFILE
echo "#TOTALSIZE NCHUNKS BUFSIZE CHUNKSIZE AVX-MWAY-TIME(usecs) AVX-MWAY-TPUT AVX-MWAY-MB/s SCALAR-TIME(usecs) SCALAR-TPUT SCALAR-MB/s MEMCPY-TIME(usecs) MEMCPY-TPUT MEMCPY-MB/s" >> $OUTFILE

for TOTALSIZE in 4194304 8388608 16777216
do
    for NClog2 in $(seq 2 11) # 4 --> 2048
    do
		NCHUNKS=$((1<<$NClog2))
		CHUNKSIZE=$(($TOTALSIZE/$NCHUNKS))
		
		for BS in $(seq 1 8)
		do
	    	BUFSIZE=$(($L1SIZE*(1<<$BS)))
	
	    	if [[ "$TOTALSIZE" -le "$MAXNUMTUPLES" ]]
	    	then
				EXE="$BIN $CHUNKSIZE $NCHUNKS $BUFSIZE"
				echo "========================================" >> $LOGFILE
				echo $EXE >> $LOGFILE
				#echo -n "$NCHUNKS $BUFSIZE $CHUNKSIZE $TOTALSIZE " >> $OUTFILE
				echo -n "$TOTALSIZE " >> $OUTFILE
				$EXE 2>> $LOGFILE 1>> $OUTFILE
	    	fi
		done
    done
done