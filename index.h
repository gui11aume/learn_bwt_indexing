#ifndef _BWT_INDEX_H
#define _BWT_INDEX_H

// Size of the alphabet.
#define AZ 4

#define exit_if_null(x) \
   do { if ((x) == NULL) { fprintf(stderr, "memory error %s:%d:%s()\n", \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); }} while(0)

#define exit_cannot_open(x) \
   do { fprintf(stderr, "cannot open file '%s' %s:%d:%s()\n", (x), \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); } while(0)

#define exit_if(x) \
   do { if (x) { fprintf(stderr, "%s %s:%d:%s()\n", #x, \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); }} while(0)

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


// The 'Occ_t' struct contains a size variable 'sz', followed by
// 'AZ' arrays of 'block_t', where 'AZ' is the number of letters in
// the alphabet. Note that 'sz' is not the number of 'block_t' in
// the arrays, but the size of the BWT, including the termination
// character.
struct Occ_t {
          size_t   yz;       // Size of the text.
          size_t   nb;       // Number of entries.
          size_t   C[AZ+1];  // The 'C' array.
   struct block_t  rows[0];  // Entries.
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
   size_t   yz;
   size_t   nb;
   size_t   zero;
   uint8_t  slots[0];
};

typedef struct block_t block_t;
typedef struct BWT_t   BWT_t;
typedef struct Occ_t   Occ_t;
typedef struct SA_t    SA_t;


const char ALPHABET[4] = "ACGT";

const uint8_t NONALPHABET[256] = {
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1, 1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,
    1,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1, 1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
};


#define A 'A'
#define C 'C'
#define G 'G'
#define T 'T'

const char REVCOMP[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,T,0,G,0,0,0,C,0,0,0,0,0,0,0,0, 0,0,0,0,A,0,0,0,0,0,0,0,0,0,0,0,
    0,T,0,G,0,0,0,C,0,0,0,0,0,0,0,0, 0,0,0,0,A,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

#undef A
#undef C
#undef G
#undef T

#endif
