#!/bin/sh

VEC=sse3
START=1
END=80

make clean
make

for i in `seq $START 1 $END`
do
        echo "START TEST WITH MATRIX DIMENISON $i"
	rm -rf gen_matmul_dense.hpp
	./generator.exe dense gen_matmul_dense.hpp dense_test $i 9 $i $i $i $i > /dev/null
	icpc -O3 -m$VEC -ansi-alias -DNDEBUG -DORDER_NUMBER=$i dgemm_test.cpp
	./a.out
	echo "END TEST WITH MATRIX DIMENISON $i"
done
