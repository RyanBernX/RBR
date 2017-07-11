CC=gcc
FC=gfortran
FLAGS=-O2 -fPIC

.PHONY: all clean

all: dsort.o kb07ad.o

dselect.o: dsort.c
	${CC} -O3 -c -o $@ $<

solve_sub_U.o: solve_sub_U.cpp
	${CXX} -Wall -O3 -c -o $@ $< --std=c++11 -I../../include

dsort.o: kb07ad.c
	${CC} ${FLAGS} -c -o $@ $<

kb07ad.o: kb07ad.f
	${FC} ${FLAGS} -c -o $@ $<

ex01: ex01.c dselect.o dsort.o kb07ad.o solve_sub_U.o
	${CC} -Wall -o ex01 -O3 $^ -lgfortran -I../../include

clean:
	rm -f dsort.o kb07ad.o dselect.o
