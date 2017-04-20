CC= gcc
#CFLAGS= -std=c99 -Wall -g -O0
#CFLAGS= -std=c99 -Wall -O3 -g -pg -fprofile-arcs -ftest-coverage
CFLAGS= -std=c99 -Wall -O3

all: index seed

index: index.c divsufsort.o index.h
	$(CC) $(CFLAGS) index.c divsufsort.o -o index

seed: seed.c index.h
	$(CC) $(CFLAGS) seed.c -o seed


do: src2.c divsufsort.o
	$(CC) $(CFLAGS) src2.c divsufsort.o -o do

clean:
	rm divsufsort.o index seed
