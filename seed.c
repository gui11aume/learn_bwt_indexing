#include "index.h"


size_t
get_rank
(
   Occ_t   * Occ,
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
   Occ_t  * Occ
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


int main(int argc, char ** argv) {

   // Sanity checks.
   exit_if(argc != 3);
   exit_if(strlen(argv[1]) > 250);


   // Load index files.
   BWT_t  * BWT;
   Occ_t  * Occ;
   SA_t   * SA;

   size_t mmsz;
   char buff[256];

   sprintf(buff, "%s.sa", argv[1]);
   int fsar = open(buff, O_RDONLY);
   if (fsar < 0) exit_cannot_open(buff);

   mmsz = lseek(fsar, 0, SEEK_END);
   SA = (SA_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, fsar, 0);
   exit_if(SA == NULL);
   close(fsar);


   sprintf(buff, "%s.bwt", argv[1]);
   int fbwt = open(buff, O_RDONLY);
   if (fbwt < 0) exit_cannot_open(buff);

   mmsz = lseek(fbwt, 0, SEEK_END);
   BWT = (BWT_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, fbwt, 0);
   exit_if(BWT == NULL);
   close(fbwt);


   sprintf(buff, "%s.occ", argv[1]);
   int focc = open(buff, O_RDONLY);
   if (focc < 0) exit_cannot_open(buff);

   mmsz = lseek(focc, 0, SEEK_END);
   Occ = (Occ_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, focc, 0);
   exit_if(Occ == NULL);
   close(focc);

   // Open seq file.
   FILE * fseq = fopen(argv[2], "r");
   if (fseq == NULL) exit_cannot_open(argv[2]);
   
   // Read file line by line.
   ssize_t rlen;
   size_t sz = 64; 
   char * buffer = malloc(64);
   exit_if_null(buffer);

   while ((rlen = getline(&buffer, &sz, fseq)) != -1) {
      size_t truth;
      buffer[rlen-1] = '\0'; 
      if (buffer[0] == '>') {
         fprintf(stdout, "%s\n", buffer);
         truth = atoi(&buffer[1]);
         continue;
      }
      range_t range = backward_search(buffer, 20, Occ);
      fprintf(stdout, "%s %zu:%zu\n", buffer, range.bot, range.top);
      for (size_t pos = range.bot ; pos <= range.top ; pos++) {
         fprintf(stdout, "pos: %zu\n", query_SA(SA, BWT, Occ, pos));
      }
   }

}
