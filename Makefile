P= index seed

CC= gcc
CFLAGS= -std=c99 -Wall -DASMAIN

all: CFLAGS += -DNDEBUG -O3
all: $(P)

debug: CFLAGS += -DDEBUG -g -O0
debug: $(P)

index: index.c divsufsort.o search.o bwt.h
	$(CC) $(CFLAGS) index.c divsufsort.o search.o -o index

seed: seed.c search.o bwt.h
	$(CC) $(CFLAGS) seed.c search.o -o seed

clean:
	rm -f divsufsort.o search.o $(P)
