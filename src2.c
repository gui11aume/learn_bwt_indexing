#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "divsufsort.h"

// Size of the alphabet.
#define N 4

// Useful macros.
#define ALL64 ((uint64_t) 0xFFFFFFFFFFFFFFFF)

const char ENCODE[256] = { ['C'] = 1, ['G'] = 2, ['T'] = 3 };

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
// symbol of the alphabet, an Occ table occupies '2*N' bits per
// letter of the BWT (i.e. one byte when 'N' is 4).
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
// 'N' arrays of 'block_t', where 'N' is the number of letters in
// the alphabet. Note that 'sz' is not the number of 'block_t' in
// the arrays, but the size of the BWT, including the termination
// character.

struct Occ_t {
          size_t    sz;       // Size of the text.
   struct block_t * rows[N];  // Occ entries.
};


// The 'range_t' struct stores the result of a range query.

struct range_t {
   size_t bot;
   size_t top;
};


// The suffix array.

struct SA_t {
   size_t    sz;
   int64_t * _;
};


// The compressed suffix array.

struct cSA_t {
   size_t      sz;
   uint8_t     nbits;
   uint64_t  * bitf;
};

// The Burrow-Wheeler transform.

struct BWT_t {
   size_t    zero;
   size_t    sz;
   uint8_t * _;
};

typedef struct Occ_t Occ_t;
typedef struct block_t block_t;
typedef struct range_t range_t;
typedef struct cSA_t cSA_t;
typedef struct SA_t SA_t;
typedef struct BWT_t BWT_t;


SA_t
compute_SA
(
   const char * txt
)
{
   SA_t SA = {0};
   SA.sz = strlen(txt) + 1;
   SA._ = malloc(SA.sz * sizeof(int64_t));
   if (SA._ == NULL) exit(1);
   divsufsort((const unsigned char *) txt, SA._, SA.sz);
   return SA;
}


Occ_t
new_Occ
(
   size_t sz
)
// Create a new Occ table, where 'sz' is the size of the text
// (not the number of blocks in the Occ table).
{
   Occ_t new = {0};
   new.sz = sz;
   for (int i = 0 ; i < N ; i++) {
      new.rows[i] = calloc(1+(sz-1)/32, sizeof(block_t));
      if (new.rows[i] == NULL) exit(1);
   }
   return new;
}


int
popcount
(
   uint32_t bits
)
// Compute the popcount of an 'int' of 32 bits. Divide it in 4
// segments of 8 bits and make 4 references to the precomputed
// lookup array POPOCOUNT. The lookup is small and fits in only
// 4 cache lines, so most of the references to POPCOUNT during
// the backward search will be cache hits (i.e > 10 times faster
// than the references to 'block_t', which are mostly misses).
{
  const uint8_t POPCOUNT[256] = {
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4, 1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5, 2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5, 2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6, 3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5, 2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6, 3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6, 3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7, 4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
  };

  // Shift repeatedly to the right compute the
  // popcount of the lowest 8 bits and add.
  return POPCOUNT[bits       & 0b11111111] + 
         POPCOUNT[bits >> 8  & 0b11111111] + 
         POPCOUNT[bits >> 16 & 0b11111111] + 
         POPCOUNT[bits >> 24 & 0b11111111];
}

// The text has 13 characters, but in C it contains an
// extra byte equal to '\0' at the end. This will act
// as the terminator '$'.



BWT_t
compute_BWT
(
   const char   * TXT,
         SA_t     SA
)
// Used to compute the BWT from the SA and the text.
// The characters are encoded as numbers between 0 and 3
// i.e. with 2 bits.
{

   BWT_t BWT = {0};
   BWT.sz = (SA.sz + (4-1)) / 4;
   BWT._ = calloc(BWT.sz, sizeof(uint8_t));
   if (BWT._ == NULL) exit(1);

   for (size_t i = 0 ; i < SA.sz ; i++) {
      if (SA._[i] > 0) { 
         uint8_t c = ENCODE[(uint8_t) TXT[SA._[i]-1]];
         BWT._[i/4] |= c << 2*(i%4);
      }
      else {
         // Record the position of the zero.
         BWT.zero = i;
      }
   }

   return BWT;

}


void
write_blocks
(
   Occ_t      Occ,
   uint32_t * smpl,
   uint32_t * bits,
   size_t     idx    // Index of 'block_t' in array.
)
// Write 'N' smpl/bits blocks to the 'block_t' arrays of 'Occ'
// at position 'pos' (the array index and not the position in
// the BWT).
{
   for (int i = 0 ; i < N ; i++) {
      Occ.rows[i][idx].smpl = smpl[i] - popcount(bits[i]);
      Occ.rows[i][idx].bits = bits[i];
   }
}


void
compute_C_and_Occ
(
   BWT_t   BWT,
   int   * C,
   Occ_t   Occ
)
{

   uint32_t smpl[N] = {0};
   uint32_t bits[N] = {0};

   for (int i = 0 ; i < Occ.sz ; i++) {
      // Extract symbol at position 'i' from BWT.
      uint8_t c = BWT._[i/4] >> 2*(i%4) & 0b11;
      if (i != BWT.zero) {
         smpl[c]++;
         bits[c] |= (1 << (31 - i % 32));
      }
      if (i % 4 == 31) { // Write every 32 entries.
         write_blocks(Occ, smpl, bits, i/32);
         bzero(bits, sizeof(bits));
      }
   }

   write_blocks(Occ, smpl, bits, (Occ.sz-1)/32);

   C[0] = 1;
   for (int i = 1 ; i < N+1 ; i++) {
      C[i] = C[i-1] + smpl[i-1];
   }

}


size_t
get_rank
(
   block_t * row,
   size_t    pos
)
{
   if (pos == -1) return 0;
   uint32_t smpl = row[pos/32].smpl;
   uint32_t bits = row[pos/32].bits;
   return smpl + popcount(bits >> (31 - pos % 32));
}


range_t
backward_search
(
         char   * query,
   const int    * C,
         Occ_t    Occ
)
// Used to search a substring using 'Occ' and 'C'.
// Return (0,0) in case the query is not found.
{
   range_t range = { .bot = 0, .top = C[N]};
   for (int i = strlen(query)-1 ; i >= 0 ; i--) {
      int c = ENCODE[(uint8_t) query[i]];
      range.bot = C[c] + get_rank(Occ.rows[c], range.bot - 1);
      range.top = C[c] + get_rank(Occ.rows[c], range.top) -1;
      if (range.top < range.top) return (range_t) {0};
   }
   return range;
}


cSA_t
downsample_SA
(
   SA_t SA
)
{

   cSA_t cSA = {0};

	// Get the length of the suffix array
	// from the first position.
   size_t initsz = SA.sz;

   // Compute the number of bits required
   // for each entry.
   while (initsz > ((uint64_t) 1 << cSA.nbits)) cSA.nbits++;

   // Set a mask for the 'nbits' lower bits.
   uint64_t mask = ALL64 >> (64-cSA.nbits);
   
   // Clear upper bits of array[0].
   SA._[0] &= mask;
   uint8_t  lastbit = 0;

   // Sample every 16-th value.
   for (size_t i = 0 ; i < initsz; i += 16) {
      // Save the current value.
      int32_t current = SA._[i];
      // Store the compact version.
      SA._[cSA.sz] |= (current & mask) << lastbit;
      // Update bit offset.
      lastbit += cSA.nbits;
      // Word is full.
      if (lastbit >= 64) {
         lastbit = lastbit - 64;
         // Complete with remainder or set to 0 (if lastbit = 0).
         // This will clear the upper bits of array.
         SA._[++cSA.sz] = (current & mask) >> (cSA.nbits - lastbit);
      }
   }

   cSA.sz += (lastbit > 0);
   cSA.bitf = realloc(SA._, cSA.sz * sizeof(uint64_t));

   return cSA;

}


size_t
query_cSA
(
   cSA_t     cSA,
   BWT_t      BWT,
   int      * C,
   Occ_t      Occ,
   size_t     pos
)
{
   if (pos == BWT.zero) return 0;
   if (pos % 16 == 0) {
      // Value is sampled. Extract it.
      size_t  idx = pos / 16;
      size_t  lo = cSA.nbits * idx;
      size_t  hi = cSA.nbits * (idx+1)-1;
      uint64_t mask = ALL64 >> (64-cSA.nbits);
      if (lo/64 == hi/64) {
         // Entry fits in a single 'uint64_t'.
         return cSA.bitf[lo/64] >> lo%64 & mask;
      }
      else {
         // Entry is split between two 'uint64_t'.
         return (cSA.bitf[lo/64] >> lo%64 |
                  cSA.bitf[hi/64] << lo%64) & mask;
      }
   }
   uint8_t c = BWT._[pos/4] >> 2*(pos%4) & 0b11;
   size_t nextpos = C[c] + get_rank(Occ.rows[c], pos) - 1;
   return query_cSA(cSA, BWT, C, Occ, nextpos) + 1;
}


int main(void) {

   const char TXT[] = "GATGCGAGACTCGAGATG";

   SA_t  SA  = compute_SA(TXT);
   BWT_t BWT = compute_BWT(TXT, SA);
   cSA_t cSA = downsample_SA(SA);

   int C[N+1] = {0};
   Occ_t Occ = new_Occ(19);
   compute_C_and_Occ(BWT, C, Occ);

   range_t range = backward_search("GAGA", C, Occ);
   fprintf(stdout, "SA range %zu:%zu\n", range.bot, range.top);
   for (size_t i = range.bot ; i <= range.top ; i++) {
      fprintf(stdout, "Text position: %zu\n",
          query_cSA(cSA, BWT, C, Occ, i));
   }

   free(cSA.bitf);
   free(BWT._);
   for (int i = 0 ; i < N ; i++) free(Occ.rows[i]);

}
