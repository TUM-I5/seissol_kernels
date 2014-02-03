#!/bin/sh

VEC=avx
ORDERS="4 10 20 35 56 84"

make clean
make

for i in ${ORDERS}
do
        echo "START TEST WITH MATRIX DIMENISON $i"
	rm -rf gen_matmul_dense.hpp
	./generator.exe denseCK gen_matmul_dense.hpp dense_test $i 9 $i $i $i $i 1 > /dev/null
#	icpc -O3 -xCORE-AVX2 -fma -ansi-alias -DNDEBUG -DORDER_NUMBER=$i dgemm_test.cpp
	icpc -O3 -m$VEC -ansi-alias -DNDEBUG -DORDER_NUMBER=$i dgemm_test_CK.cpp
	./a.out
	echo "END TEST WITH MATRIX DIMENISON $i"
done
