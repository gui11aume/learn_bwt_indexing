#include "bwt.h"


const char ALPHABET[4] = "ACGT";
const char ENCODE[256] = { ['c'] = 1, ['g'] = 2, ['t'] = 3,
   ['C'] = 1, ['G'] = 2, ['T'] = 3 };

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

size_t
get_rank
(
   const Occ_t   * Occ,
         uint8_t   c,
         size_t    pos
)
{
   if (pos == -1) return Occ->C[c];
   uint32_t smpl = Occ->rows[c*Occ->nb + pos/32].smpl;
   uint32_t bits = Occ->rows[c*Occ->nb + pos/32].bits;
   return Occ->C[c] + smpl + __builtin_popcountl(bits >> (31 - pos % 32));

}


range_t
backward_search
(
         char   * query,
         size_t   len,
   const Occ_t  * Occ
)
// Used to search a substring using 'Occ' and 'C'.
// Return (0,0) in case the query is not found.
{

   range_t range = { .bot = 0, .top = Occ->yz-1 };
   size_t pos = len-1;

   if (len >= HSTUB) {
      size_t merid = 0;
      for ( ; pos >= len-HSTUB ; pos--) {
         merid = (merid << 2) + ENCODE[(uint8_t) query[pos]];
      }
      range = Occ->stub[merid];
   }

   for ( ; pos != -1 ; pos--) {
      int c = ENCODE[(uint8_t) query[pos]];
      range.bot = get_rank(Occ, c, range.bot - 1);
      range.top = get_rank(Occ, c, range.top) - 1;
      if (range.top < range.bot) return (range_t) {0};
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
         return (size_t) SA->bitf[lo/64] >> lo % 64 & mask;
      }
      else {
         // Entry is split between two 'uint64_t'.
         size_t lo_bits = (size_t) SA->bitf[lo/64] >> lo % 64;
         size_t hi_bits = (size_t) SA->bitf[hi/64] << (64-lo) % 64;
         return (lo_bits | hi_bits) & mask;
      }
   }
   uint8_t c = BWT->slots[pos/4] >> 2*(pos % 4) & 0b11;
   size_t nextpos = get_rank(Occ, c, pos) - 1;
   return query_SA(SA, BWT, Occ, nextpos) + 1;
}
