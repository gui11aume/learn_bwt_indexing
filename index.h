#define _GNU_SOURCE
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>


#ifndef _BWT_INDEX_H
#define _BWT_INDEX_H

// ------- Diverse definitions ------- //

// Size of the alphabet. Note that the code will break if
// the value is changed. It is indicated here in order
// to highlight where the alphabet size is important.
#define AZ 4

#define ALL64 ((uint64_t) 0xFFFFFFFFFFFFFFFF)

// The 'mmap()' option 'MAP_POPULATE' is available
// only on Linux and from kern version 2.6.23.
#if __linux__
  #include <linux/version.h>
  #if LINUX_VERSION_CODE > KERNEL_VERSION(2,6,22)
    #define _MAP_POPULATE_AVAILABLE
  #endif
#endif

#ifdef _MAP_POPULATE_AVAILABLE
  #define MMAP_FLAGS (MAP_PRIVATE | MAP_POPULATE)
#else
  #define MMAP_FLAGS MAP_PRIVATE
#endif


// ------- Error handling macros ------- //

#define exit_if_null(x) \
   do { if ((x) == NULL) { fprintf(stderr, "memory error %s:%d:%s()\n", \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); }} while(0)

#define exit_cannot_open(x) \
   do { fprintf(stderr, "cannot open file '%s' %s:%d:%s()\n", (x), \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); } while(0)

#define exit_if(x) \
   do { if (x) { fprintf(stderr, "%s %s:%d:%s()\n", #x, \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); }} while(0)


// ------- Type definitions ------- //

// Entries of the Occ table combine an Occ value (the cumulative
// number of occurrences of the character in the BWT up to and
// including the given position), and a bitfield where the value
// is set to 1 if and only if the BWT has the character at this
// position.
//
// Example 'block_t' data:
//
//  |  .smpl (uint32_t)  |          .bits (uint32_t)          |
//  |       146883       | 0b00100011100000000110000001000000 |
//
// The next .smpl value in an array of 'block_t' must be 146890
// (the popcount of the contiguous .bits value).
//
// To avoid confusion, we will use the term "index" when referring
// to 'block_t' structs of an array, and "position" when referring
// to the text or the BWT. Thus, the rank is taken on a position
// and not an index.
//
// There is one 'block_t' per 32 letters of the BWT. Since each
// 'block_t' occupies 64 bits, an array of 'block_t' occupies
// 2 bits per letter of the BWT. Since there is one array per
// symbol of the alphabet, an Occ table occupies '2*AZ' bits per
// letter of the BWT (i.e. one byte when 'AZ' is 4).
//
// Both .smpl and .bits are retrieved in a single memory reference
// (a 'block_t' occupies 8 bytes and cache lines are 64 bytes on
// x86), so both values are available for the price of a single
// cache miss.
struct block_t {
   unsigned int smpl : 32;    // Sampled Occ values.
   unsigned int bits : 32;    // Bitfield Occ values.
};

// Return type for range queries.
struct range_t {
   size_t bot;
   size_t top;
};

// The 'Occ_t' struct contains a size variable 'sz', followed by
// 'AZ' arrays of 'block_t', where 'AZ' is the number of letters in
// the alphabet. Note that 'sz' is not the number of 'block_t' in
// the arrays, but the size of the BWT, including the termination
// character.
#define HSTUB 12
#define MERS 1 << (2*HSTUB)
struct Occ_t {
          size_t   yz;         // 'strlen(txt) + 1'.
          size_t   nb;         // Number of entries.
          size_t   C[AZ+1];    // The 'C' array.
   struct range_t  stub[MERS]; // The stub table.
   struct block_t  rows[0];    // Entries.
};

// The (compressed) suffix array.
struct SA_t {
   size_t   yz;
   size_t   nb;
   size_t   nbits;
   int64_t  bitf[0];
};

// The Burrow-Wheeler transform.
struct BWT_t {
   size_t   yz;        // 'strlen(txt) + 1'.
   size_t   nb;        // Number of slots.
   size_t   zero;      // Position of '$'.
   uint8_t  slots[0];  // Characters.
};


typedef struct block_t block_t;
typedef struct BWT_t   BWT_t;
typedef struct Occ_t   Occ_t;
typedef struct range_t range_t;
typedef struct SA_t    SA_t;


// ------- Text manipulation ------- //

const char ALPHABET[4] = "ACGT";
const char ENCODE[256] = { ['C'] = 1, ['G'] = 2, ['T'] = 3 };

const uint8_t NONALPHABET[256] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
};

const char REVCOMP[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0,'T',0,'G',0, 0, 0,'C',0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0,'A',0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0,'T',0,'G',0, 0, 0,'C',0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0,'A',0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
};

#endif
