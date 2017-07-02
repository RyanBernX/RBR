CC=gcc
CDEBUG=-Wall
COPTFLAG=-O3
LIBS=-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lgomp -lpthread -lgfortran -lm -ldl
CFLAG = -m64 -std=c99 -I${MKLROOT}/include -I../include -fopenmp -fPIC

OBJS=DCRBR_rounding.o spmat.o
CONTRIB_OBJS=contrib/kb07ad.o contrib/dsort.o

.PHONY: clean exec lib contrib

all: exec

exec: ex01 ex02 ex03 ex04 contrib_test

lib: $(OBJS)
	ar -rcs libDCRBR.a $^

ex01: test_synth.o $(OBJS)
	$(CC) -o ex01 $^ $(CONTRIB_OBJS) $(CDEBUG) $(COPTFLAG) $(LIBS) $(CFLAG)

ex02: test_small_real.o $(OBJS)
	$(CC) -o ex02 $^ $(CONTRIB_OBJS) $(CDEBUG) $(COPTFLAG) $(LIBS) $(CFLAG)

ex03: test_large_real.o $(OBJS)
	$(CC) -o ex03 $^ $(CONTRIB_OBJS) $(CDEBUG) $(COPTFLAG) $(LIBS) $(CFLAG)

ex04: test_large_real_p.o $(OBJS)
	$(CC) -o ex04 $^ $(CONTRIB_OBJS) $(CDEBUG) $(COPTFLAG) $(LIBS) $(CFLAG)

contrib:
	cd contrib && make all && cd ../

%.o: %.c
	$(CC) -c $< $(CDEBUG) $(COPTFLAG) $(CFLAG)

clean: 
	rm -f *.o *.a ex01 ex02 ex03 ex04
	cd contrib && make clean && cd ../
	
