CC=gcc
FC=gfortran
FLAGS=-O2 -fPIC

.PHONY: all clean

all: dsort.o kb07ad.o

dsort.o: kb07ad.c
	${CC} ${FLAGS} -c -o $@ $<

kb07ad.o: kb07ad.f
	${FC} ${FLAGS} -c -o $@ $<

clean:
	rm -f dsort.o kb07ad.o
