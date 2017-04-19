#define _GNU_SOURCE
#include <fcntl.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <unistd.h>

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


int main(int argc, char ** argv) {

   // Sanity checks.
   exit_if(argc != 2);
   exit_if(strlen(argv[1]) > 250);

   // Open seq file.
   FILE * fasta = fopen(argv[1], "r");
   if (fasta == NULL) exit_cannot_open(argv[1]);

   range_t range = backward_search("GAGA", Occ);

   fprintf(stdout, "SA range %zu:%zu\n", range.bot, range.top);
   for (size_t i = range.bot ; i <= range.top ; i++) {
      fprintf(stdout, "Text position: %zu\n",
          query_SA(SA, BWT, Occ, i));
   }

   // Clean up.
   free(SA);
   free(BWT);
   free(Occ);

}
