CC=icpc
LD=ifort

#Westmere
CFLAGS=-O2 -msse3 -Wall -DDEBUG -ansi -vec-report3 -no-opt-matmul -ipo
FORTRAN90FLAGS=-O2 -msse3 -DDEBUG -no-opt-matmul -r8 -fpp -ipo
#Sandy Brdige (-r8 is mandatory due to missing precision in the SeisSol routines)
#CFLAGS=-O2 -mavx -Wall -ansi -vec-report3 -no-opt-matmul -ipo
LDFLAGS=-cxxlib -nofor_main

all : generator

gen_objs = GeneratorCSC.o ReaderCSC.o GeneratorDense.o

GeneratorDriver.o : GeneratorDriver.cpp
	$(CC) $(CFLAGS) -c GeneratorDriver.cpp
		
GeneratorCSC.o : GeneratorCSC.cpp
	$(CC) $(CFLAGS) -c GeneratorCSC.cpp

GeneratorDense.o : GeneratorDense.cpp
	$(CC) $(CFLAGS) -c GeneratorDense.cpp
	
ReaderCSC.o : ReaderCSC.cpp
	$(CC) $(CFLAGS) -c ReaderCSC.cpp
	
generator : $(gen_objs) GeneratorDriver.o
	$(CC) $(CFLAGS) $(gen_objs) GeneratorDriver.o -o generator.exe
	
clean :
	rm -rfv *.o *.mod *.exe
