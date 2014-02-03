#!/bin/sh

ORDERS="56"

make clean
make

for i in ${ORDERS}
do
        echo "START TEST WITH MATRIX DIMENISON $i"
	rm -rf gen_matmul_dense.hpp
	./generator.exe denseCK gen_matmul_dense.hpp dense_test $i 9 $i $i $i $i 1 > /dev/null
	icpc -O3 -mmic -fma -ansi-alias -DNDEBUG -DORDER_NUMBER=$i dgemm_test_CK.cpp
        scp a.out mic0:
	ssh mic0 "./a.out"
	echo "END TEST WITH MATRIX DIMENISON $i"
done
