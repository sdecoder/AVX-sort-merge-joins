#!/bin/bash
#
# CONFIGURE="./configure --disable-key8B"
# build the binary
OUTPUTDIR=../measurements/
cd ../
# $CONFIGURE
cd ./src/

OUTPUT=$OUTPUTDIR/partitioning-bench-8Btuples.txt

echo "#NUMR RADIXBITS NORMAL-RADIX(usecs,tput) RDX-SWBUFFERS RDX-SWBUFFERS++ HIST+MEMCPY PLAIN-MEMCPY" >> $OUTPUT

# echo $makecmd
# makecmd="make clean; make partitioningbench"
# $makecmd

for n in 8 #16 32 64 128 256 512
do
	NUMR=$((n*1024*1024))

	for RDXBITS in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
	do
		echo -n "$NUMR $RDXBITS " >> $OUTPUT

		for WHAT in 0 1 2 3 4
		do
			cmd="./bench_partitioning $NUMR $WHAT $RDXBITS"
			echo $cmd
			$cmd 2>> $OUTPUT
		done

		echo " " >> $OUTPUT
	done
done
