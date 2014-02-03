#!/bin/sh

VEC=avx
START=1
END=84
ORDERS="4 10 20 35 56 84"
#START=56
#END=56

make clean
make

#for i in `seq $START 1 $END`
for i in ${ORDERS}
do
        echo "START TEST WITH MATRIX DIMENISON $i"
	rm -rf gen_matmul_dense.hpp
	./generator.exe dense gen_matmul_dense.hpp dense_test_square $i 9 $i $i $i $i 1 > /dev/null
        ./generator.exe dense gen_matmul_dense.hpp dense_test_rect $i 9 9 $i 9 $i 1 > /dev/null
#	icpc -O3 -xCORE-AVX2 -fma -ansi-alias -DNDEBUG -DORDER_NUMBER=$i dgemm_test.cpp
	icpc -O3 -m$VEC -ansi-alias -DNDEBUG -DORDER_NUMBER=$i dgemm_test.cpp
	./a.out
	echo "END TEST WITH MATRIX DIMENISON $i"
done
