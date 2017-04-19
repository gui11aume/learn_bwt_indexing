#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

// Size of the alphabet.
#define N 4

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


typedef struct Occ_t Occ_t;
typedef struct block_t block_t;
typedef struct range_t range_t;



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

const char TXT[] = "GATGCGAGAGATG";


int compare_suffixes
(
   const void * a,
   const void * b
)
// Used to put the suffix in lexicographical order.
// Do not use this approach on large texts. Use for instance
// https://github.com/y-256/libdivsufsort
{
   const char * suff_a = &TXT[*(int *)a];
   const char * suff_b = &TXT[*(int *)b];
   return strcmp(suff_a, suff_b);
}

void down_sample_SA
(
   const int  SA[14],
         char dSA[2]
)
// Used to down-sample the suffix array by a factor 4.
// To store numbers up to 13 we need 4 bits (56 bits total).
// We downsample 4 times, so we need 14 bits or 2 bytes.
{
   int dSA_index = 0;
   int dSA_shift = 0;
   for (int i = 0 ; i < 14 ; i += 4) {
      dSA[dSA_index] |= SA[i] << dSA_shift; 
      if (dSA_shift == 0) { dSA_shift = 4; }
      else                { dSA_shift = 0; dSA_index++; }
   }
}


void compute_BWT
(
   const int    * SA,
         char   * BWT,
         size_t   sz
)
// Used to compute the BWT from the SA and the text.
// The characters are encoded as numbers between 0 and 3.
// The terminator is encoded as '$' (could be anything else).
{
   for (int i = 0 ; i < sz ; i++) {
      if (SA[i] > 0) BWT[i] = ENCODE[(uint8_t) TXT[SA[i]-1]];
      else           BWT[i] = '$';
   }
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
   const char  * BWT,
         int   * C,
         Occ_t   Occ
)
{

   uint32_t smpl[N] = {0};
   uint32_t bits[N] = {0};

   for (int i = 0 ; i < Occ.sz ; i++) {
      if (BWT[i] < N) { // Skip the '$' symbol.
         smpl[(uint8_t) BWT[i]]++;
         bits[(uint8_t) BWT[i]] |= (1 << (31 - i % 32));
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


size_t
query_dSA
(
   const char   * dSA,
   const char   * BWT,
         int    * C,
         Occ_t    Occ,
         size_t   pos
)
{
   if (pos % 4 == 0) {
      return dSA[pos/8] >> (pos % 8) & 0b1111;
   }
   int c = BWT[pos];
   size_t nextpos = C[c] + get_rank(Occ.rows[c], pos) - 1;
   return query_dSA(dSA, BWT, C, Occ, nextpos) + 1;
}


int main(void) {

   int SA[14] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
   qsort(SA, 14, sizeof(int), compare_suffixes);
   
   char BWT[14] = {0};
   compute_BWT(SA, BWT, 14);

   char dSA[2] = {0};
   down_sample_SA(SA, dSA);

   // From that point we don't need 'SA' anymore.

   int C[N+1] = {0};
   Occ_t Occ = new_Occ(14);
   compute_C_and_Occ(BWT, C, Occ);

   range_t range = backward_search("GAGA", C, Occ);
   fprintf(stdout, "%s\n", TXT);
   fprintf(stdout, "SA range %zu:%zu\n", range.bot, range.top);
   for (size_t i = range.bot ; i <= range.top ; i++) {
      fprintf(stdout, "Text position: %zu\n",
          query_dSA(dSA, BWT, C, Occ, i));
   }

}
