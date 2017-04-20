#define _GNU_SOURCE
#include <fcntl.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

#include "index.h"

// Useful macros.
#define ALL64 ((uint64_t) 0xFFFFFFFFFFFFFFFF)

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
   // This option is 10-15% slower than built in popcount.
   //return Occ->C[c] + smpl + popcount(bits >> (31 - pos % 32));
   return Occ->C[c] + smpl + __builtin_popcountl(bits >> (31 - pos % 32));

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

   size_t m = strlen(query);
   int i0 = m-1;

   range_t range = { .bot = 0, .top = Occ->yz-1 };

   if (strlen(query) >= HSTUB) {
      size_t merid = 0;
      for (int j = 0 ; j < HSTUB ; j++) {
         merid = (merid << 2) + ENCODE[(uint8_t) query[m-1-j]];
      }
      range = Occ->stub[merid];
      i0 = m-1 - HSTUB;
   }

   for (int i = i0 ; i >= 0 ; i--) {
      int c = ENCODE[(uint8_t) query[i]];
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


   // Load index files.

   BWT_t  * BWT;
   Occ_t  * Occ;
   SA_t   * SA;

   size_t mmsz;

   int fsar = open("index.sa", O_RDONLY);
   if (fsar < 0) exit_cannot_open("index.sa");

   mmsz = lseek(fsar, 0, SEEK_END);
   SA = (SA_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, fsar, 0);
   exit_if(SA == NULL);
   close(fsar);


   int fbwt = open("index.bwt", O_RDONLY);
   if (fbwt < 0) exit_cannot_open("index.bwt");

   mmsz = lseek(fbwt, 0, SEEK_END);
   BWT = (BWT_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, fbwt, 0);
   exit_if(BWT == NULL);
   close(fbwt);


   int focc = open("index.occ", O_RDONLY);
   if (focc < 0) exit_cannot_open("index.occ");

   mmsz = lseek(focc, 0, SEEK_END);
   Occ = (Occ_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, focc, 0);
   exit_if(Occ == NULL);
   close(focc);


   // Open seq file.

   FILE * fseq = fopen(argv[1], "r");
   if (fseq == NULL) exit_cannot_open(argv[1]);
   
   // Read file line by line.
   ssize_t rlen;
   size_t sz = 64; 
   char * buffer = malloc(64);
   exit_if_null(buffer);

   while ((rlen = getline(&buffer, &sz, fseq)) != -1) {
      buffer[rlen-1] = '\0'; 
      range_t range = backward_search(buffer, Occ);
      fprintf(stdout, "%s %zu:%zu\n", buffer, range.bot, range.top);
//      fprintf(stdout, "%zu\n", query_SA(SA, BWT, Occ, range.bot));
   }

}
