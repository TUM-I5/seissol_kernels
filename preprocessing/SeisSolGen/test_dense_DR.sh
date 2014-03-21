#!/bin/sh

VEC=avx
ORDERS="56"

make clean
make

for i in ${ORDERS}
do
        echo "START TEST WITH MATRIX DIMENISON $i"
	rm -rf gen_matmul_dense.hpp
	./generator.exe dense gen_matmul_dense.hpp dense_test_square 52 36 $i 52 $i 52 1 > /dev/null
        ./generator.exe dense gen_matmul_dense.hpp dense_test_rect $i 6 9 $i 9 $i 1 > /dev/null
	icpc -O3 -m$VEC -ansi-alias -DNDEBUG -DORDER_NUMBER=$i dgemm_test_DR.cpp
#       g++ -O3 -mavx2 -mfma -funroll-loops -DNDEBUG -DORDER_NUMBER=$i dgemm_test.cpp
#	icpc -O2 -mavx -D__AVX2__ -fma -unroll-loops -DNDEBUG -DORDER_NUMBER=$i dgemm_test.cpp
	./a.out
	echo "END TEST WITH MATRIX DIMENISON $i"
done
