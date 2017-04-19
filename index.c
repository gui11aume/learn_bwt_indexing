#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "divsufsort.h"
#include "index.h"

struct range_t {
   size_t bot;
   size_t top;
};

typedef struct range_t range_t;

// Useful macros.
#define ALL64 ((uint64_t) 0xFFFFFFFFFFFFFFFF)

const char ENCODE[256] = { ['C'] = 1, ['G'] = 2, ['T'] = 3 };

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


SA_t *
create_SA
(
   const char * txt
)
{
   const size_t yz = strlen(txt) + 1;
   const size_t extra = yz * sizeof(uint64_t);

   SA_t *SA = malloc(sizeof(SA_t) + extra);
   exit_if_null(SA);

   SA->yz = yz;
   SA->nb = yz;
   SA->nbits = 64;

   divsufsort((const unsigned char *) txt, SA->bitf, yz);

   return SA;

}



BWT_t *
create_BWT
(
   const char * TXT,
         SA_t * SA
)
{

   // Do not use on compressed suffix arrays.
   if (SA->nbits != 64) return NULL;

   // Allocate new 'BWT_t'.
   const size_t yz = SA->yz;
   const size_t nb = (yz + (4-1)) / 4;
   const size_t extra = nb * sizeof(uint8_t);

   BWT_t *BWT = calloc(1, sizeof(BWT_t) + extra);
   exit_if_null(BWT);

   BWT->yz = yz;
   BWT->nb = nb;

   for (size_t pos = 0 ; pos < yz ; pos++) {
      if (SA->bitf[pos] > 0) { 
         uint8_t c = ENCODE[(uint8_t) TXT[SA->bitf[pos]-1]];
         BWT->slots[pos/4] |= c << 2*(pos % 4);
      }
      else {
         // Record the position of the zero.
         BWT->zero = pos;
      }
   }

   return BWT;

}


void
write_Occ_blocks
(
   Occ_t    * Occ,
   uint32_t * smpl,
   uint32_t * bits,
   size_t     idx    // Index of 'block_t' in array.
)
// Write 'AZ' smpl/bits blocks to the 'block_t' arrays of 'Occ'
// at position 'pos' (the array index and not the position in
// the BWT).
{
   for (int i = 0 ; i < AZ ; i++) {
      Occ->rows[Occ->nb*i + idx].smpl = smpl[i] - popcount(bits[i]);
      Occ->rows[Occ->nb*i + idx].bits = bits[i];
   }
}


Occ_t *
create_Occ
(
   BWT_t * BWT
)
{

   // Allocate new 'Occ_t'.
   const size_t yz = BWT->yz;
   const size_t nb = (yz + (32-1)) / 32;
   const size_t extra = AZ*nb * sizeof(block_t);

   Occ_t * Occ = malloc(sizeof(Occ_t) + extra);
   exit_if_null(Occ);

   Occ->nb = nb;
   Occ->yz = yz;

   uint32_t smpl[AZ] = {0};
   uint32_t bits[AZ] = {0};

   for (size_t pos = 0 ; pos < BWT->yz ; pos++) {
      // Extract symbol at position 'i' from BWT.
      uint8_t c = BWT->slots[pos/4] >> 2*(pos % 4) & 0b11;
      if (pos != BWT->zero) { // Skip the '$' symbol.
         smpl[c]++;
         bits[c] |= (1 << (31 - pos % 32));
      }
      if (pos % 32 == 31) { // Write every 32 entries.
         write_Occ_blocks(Occ, smpl, bits, pos/32);
         bzero(bits, sizeof(bits));
      }
   }

   write_Occ_blocks(Occ, smpl, bits, (BWT->yz-1)/32);

   // Write 'nb' and 'C'.
   Occ->C[0] = 1;
   for (int i = 1 ; i < AZ+1 ; i++) {
      Occ->C[i] = Occ->C[i-1] + smpl[i-1];
   }

   return Occ;

}


size_t
get_rank
(
   Occ_t   * Occ,
   uint8_t   c,
   size_t    pos
)
{
   if (pos == -1) return 1;
   uint32_t smpl = Occ->rows[c*Occ->nb + pos/32].smpl;
   uint32_t bits = Occ->rows[c*Occ->nb + pos/32].bits;
   return Occ->C[c] + smpl + popcount(bits >> (31 - pos % 32));
}


range_t
backward_search
(
         char   * query,
         Occ_t  * Occ
)
// Used to search a substring using 'Occ' and 'C'.
// Return (0,0) in case the query is not found.
{
   range_t range = { .bot = 0, .top = Occ->C[AZ]-1 };
   for (int i = strlen(query)-1 ; i >= 0 ; i--) {
      int c = ENCODE[(uint8_t) query[i]];
      range.bot = get_rank(Occ, c, range.bot - 1);
      range.top = get_rank(Occ, c, range.top) - 1;
      if (range.top < range.top) return (range_t) {0};
   }
   return range;
}


SA_t *
compress_SA
(
   SA_t * SA
)
{

   // Do not run on compressed 'SA_t'.
   if (SA->nbits != 64) return NULL;

   // Compute the number of required bits.
   size_t nbits = 0;
   while (SA->yz > ((uint64_t) 1 << nbits)) nbits++;

   // Set a mask for the 'nbits' lower bits.
   uint64_t mask = ALL64 >> (64-nbits);
   
   uint8_t lastbit = 0;
   size_t  nb = 0;

   // Sample every 16-th value.
   for (size_t pos = 0 ; pos < SA->yz ; pos += 16) {
      // Save the current value.
      int64_t current = SA->bitf[pos];
      // Store the compact version.
      SA->bitf[nb] |= (current & mask) << lastbit;
      // Update bit offset.
      lastbit += nbits;
      // Word is full.
      if (lastbit >= 64) {
         lastbit = lastbit - 64;
         // Complete with remainder or set to 0 (if lastbit = 0).
         // This will clear the upper bits of array.
         SA->bitf[++nb] = (current & mask) >> (nbits - lastbit);
      }
   }

   nb += (lastbit > 0);

   SA->nb = nb;
   SA->nbits = nbits;

   // Reallocate the 'SA_t'.
   size_t extra = (nb * nbits + (64-1)) / 64 * sizeof(uint64_t);
   return realloc(SA, sizeof(SA_t) + extra);

}


size_t
query_SA
(
   SA_t   * SA,
   BWT_t  * BWT,
   Occ_t  * Occ,
   size_t   pos
)
{
   if (pos == BWT->zero) return 0;
   if (pos % 16 == 0) {
      // Value is sampled. Extract it.
      size_t  idx = pos / 16;
      size_t  lo = SA->nbits * idx;
      size_t  hi = SA->nbits * (idx+1)-1;
      uint64_t mask = ALL64 >> (64-SA->nbits);
      if (lo/64 == hi/64) {
         // Entry fits in a single 'uint64_t'.
         return SA->bitf[lo/64] >> lo % 64 & mask;
      }
      else {
         // Entry is split between two 'uint64_t'.
         return (SA->bitf[lo/64] >> lo % 64 |
                 SA->bitf[hi/64] << lo % 64) & mask;
      }
   }
   uint8_t c = BWT->slots[pos/4] >> 2*(pos % 4) & 0b11;
   size_t nextpos = get_rank(Occ, c, pos) - 1;
   return query_SA(SA, BWT, Occ, nextpos) + 1;
}


int main(void) {

   const char TXT[] = "GATGCGAGACTCGAGATG";

   SA_t  * SA  = create_SA(TXT);
   BWT_t * BWT = create_BWT(TXT, SA);
           SA  = compress_SA(SA);
   Occ_t * Occ = create_Occ(BWT);

   range_t range = backward_search("GAGA", Occ);

   fprintf(stdout, "SA range %zu:%zu\n", range.bot, range.top);
   for (size_t i = range.bot ; i <= range.top ; i++) {
      fprintf(stdout, "Text position: %zu\n",
          query_SA(SA, BWT, Occ, i));
   }

   free(SA);
   free(BWT);
   free(Occ);

}
