NUMBASE=4
VEC=sse

CC=icpc
CFLAGS:=-ipo -ip -O3 -c -DNUMBEROFVARIABLES=9 -DNUMBEROFBASISFUNCTIONS=$(NUMBASE)
LDFLAGS=-ipo -ip -O3

ifeq ($(VEC),sse)
CFLAGS:=$(CFLAGS) -msse3
endif
ifeq ($(VEC),avx256)
CFLAGS:=$(CFLAGS) -mavx -D__AVX256__
endif
ifeq ($(VEC),avx)
CFLAGS:=$(CFLAGS) -mavx
endif

OBJS=cauchyKovalewskiTimeIntegration.o\
     SeisSolStopwatch.o\
     main.o
     
default : $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o driver.exe
	
cauchyKovalewskiTimeIntegration.o : seissol_src/Solver/cauchyKovalewskiTimeIntegration.cpp
	$(CC) $(CFLAGS) -o cauchyKovalewskiTimeIntegration.o seissol_src/Solver/cauchyKovalewskiTimeIntegration.cpp
	
SeisSolStopwatch.o : seissol_src/Monitoring/SeisSolStopwatch.cpp
	$(CC) $(CFLAGS) -o SeisSolStopwatch.o seissol_src/Monitoring/SeisSolStopwatch.cpp
	
main.o : main.cpp
	$(CC) $(CFLAGS) -o main.o main.cpp
	
clean : 
	rm -fvr *.o driver.exe
