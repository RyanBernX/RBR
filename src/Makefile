#--------------------------------------------------------------------------
# CONFIGURATION AREA
#--------------------------------------------------------------------------

# name of your C compiler
CC  = gcc

# name of your C++ compiler
CXX = g++

# name of your mex compiler (only needed if you want to compile MATLAB interface) 
MEX = /opt/MATLAB/R2019b/bin/mex

# MKL
MKL = on
MKLROOT ?= 

# CBLAS
CBLAS_HEADER = cblas.h
CBLAS_LIBS = -lblas

#--------------------------------------------------------------------------
# END OF CONFIGURATION AREA
#--------------------------------------------------------------------------

ifeq (${MKL},on)
    LIBS=-Wl,-rpath,${MKLROOT}/lib/intel64 -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lgomp -lpthread -lm -ldl
    LIBS_STATIC= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
    CFLAGS = -Wall -O3 -m64 -std=c99 -DHAVE_MKL -I${MKLROOT}/include -I../include -fopenmp -fPIC
    CXXFLAGS = -Wall -O3 -std=c++11 -DHAVE_MKL -I../include -fPIC
else
    LIBS=${CBLAS_LIBS} -lgomp -lm -ldl
    LIBS_STATIC=${LIBS}
    CFLAGS = -Wall -O3 -m64 -std=c99 -D_CBLAS_HEADER=\"${CBLAS_HEADER}\" -I../include -fopenmp -fPIC
    CXXFLAGS = -Wall -O3 -std=c++11 -D_CBLAS_HEADER=\"${CBLAS_HEADER}\" -I../include -fPIC
endif

OBJS=rbr.o kmeans.o metric.o spmat.o utils.o solve_sub_U.o normal_dist.o

.PHONY: clean exec 

all: lib exec

lib: $(OBJS)
	ar -rcs libDCRBR.a $^

mex: mex_rbr.c lib
	$(MEX) -outdir matlab CFLAGS='-std=c99 -I../include' CLIBS='-L. -lDCRBR ${LIBS_STATIC}' -O -largeArrayDims $<


exec: main.o $(OBJS)
	$(CC) -o rbr $^ $(CFLAGS) $(LIBS)

solve_sub_U.o: solve_sub_U.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

normal_dist.o: normal_dist.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

%.o: %.c
	$(CC) -c $< $(CFLAGS)

clean: 
	rm -f *.o *.a rbr matlab/mex_rbr.mexa64

