CC=gcc
CXX=g++
MEX=/opt/MATLAB/R2018b/bin/mex
CDEBUG=-Wall
COPTFLAG=-O3
LIBS=-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lgomp -lpthread -lgfortran -lm -ldl
LIBS_STATIC= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
CFLAG = -m64 -std=c99 -I${MKLROOT}/include -I../include -fopenmp -fPIC
CXXFLAG = -std=c++11 -I../include -fPIC

OBJS=rbr.o kmeans.o metric.o spmat.o utils.o contrib/solve_sub_U.o
CONTRIB_OBJS=contrib/kb07ad.o contrib/dsort.o

.PHONY: clean exec contrib

all: exec

exec: ex01 ex02 ex03 ex04

lib: $(OBJS)
	ar -rcs libDCRBR.a $^

mex: mex_rbr.c lib
	$(MEX) CFLAGS='-std=c99 -I../include' CLIBS='-L. -lDCRBR ${LIBS_STATIC}' -O -largeArrayDims $<


ex01: test_synth.o $(OBJS)
	$(CC) -o ex01 $^ $(CDEBUG) $(COPTFLAG) $(LIBS) $(CFLAG)

ex02: test_small_real.o $(OBJS)
	$(CC) -o ex02 $^ $(CDEBUG) $(COPTFLAG) $(LIBS) $(CFLAG)

ex03: test_large_real.o $(OBJS)
	$(CC) -o ex03 $^ $(CDEBUG) $(COPTFLAG) $(LIBS) $(CFLAG)

ex04: test_large_real_p.o $(OBJS)
	$(CC) -o ex04 $^ $(CDEBUG) $(COPTFLAG) $(LIBS) $(CFLAG)

ex05: test_large_real_k.o DCRBR_kmeans.o utils.o kmeans.o spmat.o contrib/solve_sub_U.o
	$(CC) -o ex05 $^ $(CDEBUG) $(COPTFLAG) $(LIBS) $(CFLAG)

kmeans: kmeans.c
	$(CC) -o $@ $^ -g $(LIBS) $(CFLAG)

contrib/solve_sub_U.o: contrib/solve_sub_U.cpp
	$(CXX) -c -o $@ $< $(CDEBUG) $(COPTFLAG) $(CXXFLAG)

%.o: %.c
	$(CC) -c $< $(CDEBUG) $(COPTFLAG) $(CFLAG)

clean: 
	rm -f *.o *.a ex01 ex02 ex03 ex04 contrib/solve_sub_U.o

	
