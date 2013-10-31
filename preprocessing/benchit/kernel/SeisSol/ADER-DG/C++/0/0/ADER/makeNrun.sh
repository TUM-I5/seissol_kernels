#!/bin/sh

VECS="sse avx avx256"
NUMBASES="4 10 20 35 56"

rm *.out

for VEC in $VECS
do
	for NUMBASE in $NUMBASES
	do
		make clean
		make NUMBASE=${NUMBASE} VEC=${VEC}
		./driver.exe 40000 40000 1 | tee ${VEC}_${NUMBASE}.out
	done
done
