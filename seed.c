#include "bwt.h"

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
