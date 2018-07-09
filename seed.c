#include "bwt.h"

int main(int argc, char ** argv) {

   // Sanity checks.
   exit_if(argc != 2);
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

   // Make all 12-mers.
   //range_t range = backward_search("ATGCTGATGTGATGTGCTGAGA", 12, Occ);
      range_t range = backward_search("AATCAAAAAAAA", 11, Occ);
      fprintf(stdout, "%ld, %ld\n", range.bot, range.top);

}
