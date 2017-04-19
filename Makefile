CC= gcc
CFLAGS= -std=c99 -Wall -g -O0

#all: do
#	./do

index: index.c divsufsort.o index.h
	$(CC) $(CFLAGS) index.c divsufsort.o -o index


do: src2.c divsufsort.o
	$(CC) $(CFLAGS) src2.c divsufsort.o -o do

clean:
	rm do divsufsort.o index
