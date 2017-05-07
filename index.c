#define _GNU_SOURCE
#include "divsufsort.h"
#include "bwt.h"

/*

NKK_t *
create_NKK
(
   const char  * genome,
   const SA_t  * SA,
   const BWT_t * BWT,
   const Occ_t * Occ
)
{

   size_t yz = SA->yz;
   size_t extra = yz * sizeof(uint8_t);

   NKK_t * NKK = calloc(1, sizeof(NKK_t) + extra);
   exit_if_null(NKK);

   NKK->yz = yz;
   NKK->nb = yz;

   for (size_t pos = 0 ; pos < yz ; pos++) {

      if (pos % 1000 == 0) {
         fprintf(stderr, "%ld\r", pos);
      }

      // Get k-mer from genome.
      char buff[21] = {0};
      memcpy(buff, &genome[pos], 20);
      
      for (int mut = 19 ; mut >= 0 ; mut--) {
         // Scan all the positions, one at a time.
         for (uint8_t c = 0 ; c < 4 ; c++) {
            if (c == genome[pos + mut]) continue;
            buff[mut] = ENCODE[c];
            range_t range = backward_search(buff, 20, Occ);
            if (range.bot && range.top - range.bot < 1) {
               for (size_t idx = range.bot ; idx <= range.top ; idx++) {
                  size_t gpos = SA->bitf[idx];
                  NKK->byte[gpos] = 1;
               }
            }
         }
         // Reset the buffer.
         buff[mut] = genome[pos + mut];
      }
   }

   return NKK;

}
*/


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
   if (SA->nb != SA->yz) return NULL;

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
      Occ->rows[Occ->nb*i + idx].smpl = smpl[i];
      Occ->rows[Occ->nb*i + idx].bits = bits[i];
   }
}


void
fill_stub
(
   Occ_t   * Occ,
   range_t   range,
   size_t    depth,
   size_t    merid
)
{
   if (depth >= HSTUB) {
      Occ->stub[merid] = range;
      return;
   }
   for (uint8_t c = 0 ; c < AZ ; c++) {
      size_t bot = get_rank(Occ, c, range.bot - 1);
      size_t top = get_rank(Occ, c, range.top) - 1;
      fill_stub(Occ, (range_t) { .bot=bot, .top=top },
            depth+1, c + (merid << 2));
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
   uint32_t diff[AZ] = {0};
   uint32_t bits[AZ] = {0};

   for (size_t pos = 0 ; pos < BWT->yz ; pos++) {
      // Extract symbol at position 'i' from BWT.
      uint8_t c = BWT->slots[pos/4] >> 2*(pos % 4) & 0b11;
      if (pos != BWT->zero) { // Skip the '$' symbol.
         diff[c]++;
         bits[c] |= (1 << (31 - pos % 32));
      }
      if (pos % 32 == 31) { // Write every 32 entries.
         write_Occ_blocks(Occ, smpl, bits, pos/32);
         memcpy(smpl, diff, AZ * sizeof(uint32_t));
         bzero(bits, sizeof(bits));
      }
   }

   write_Occ_blocks(Occ, smpl, bits, (BWT->yz-1)/32);

   // Write 'nb' and 'C'.
   Occ->C[0] = 1;
   for (int i = 1 ; i < AZ+1 ; i++) {
      Occ->C[i] = Occ->C[i-1] + diff[i-1];
   }

   range_t range = {.bot = 0, .top = yz-1};
   fill_stub(Occ, range, 0, 0);

   return Occ;

}



SA_t *
compress_SA
(
   SA_t * SA
)
{

   // Do not run on compressed 'SA_t'.
   if (SA->nb != SA->yz) return NULL;

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
   size_t extra = nb * sizeof(int64_t);
   return realloc(SA, sizeof(SA_t) + extra);

}


char *
normalize_genome
(
   FILE   * inputf
)
{

   // Read variables.
   size_t sz = 64;
   ssize_t rlen;
   char * buffer = malloc(64);
   exit_if_null(buffer);

   // Genome storage.
   size_t gsize = 0;
   size_t gbufsize = 64;
   char * genome = malloc(64); 
   exit_if_null(genome);

   while ((rlen = getline(&buffer, &sz, inputf)) != -1) {
      if (buffer[0] == '>') rlen = 1; // Use '>' as separator.
      if (gbufsize < gsize + rlen) {
         while (gbufsize < gsize + rlen) gbufsize *= 2;
         char * rsz = realloc(genome, gbufsize);
         exit_if_null(rsz);
         genome = rsz;
      }
      int one_if_newline = (buffer[rlen-1] == '\n');
      strncpy(genome + gsize, buffer, rlen - one_if_newline);
      gsize += rlen - one_if_newline;
   }

   // Normalize (use only alphabet letters).
   for (size_t pos = 0; pos < gsize ; pos++) {
      int iter = 0;
      if (NONALPHABET[(uint8_t) genome[pos]]) {
        genome[pos] = ALPHABET[iter++ % 4];
      }
   }

   // Realloc buffer.
   char * rsz = realloc(genome, 2*gsize + 1);
   exit_if_null(rsz);
   genome = rsz;

   // Reverse complement.
   size_t div = gsize;
   for (size_t pos = 0 ; pos < div ; pos++)
      genome[div + pos] = REVCOMP[(uint8_t) genome[div-pos-1]];

   gsize = 2*gsize + 1;

   // Add the terminator.
   genome[2*div] = '\0';

   // Clean up.
   free(buffer);

   return genome;

}

int main(int argc, char ** argv) {

   // Sanity checks.
   exit_if(argc != 2);
   exit_if(strlen(argv[1]) > 250);

   // Open fasta file.
   FILE * fasta = fopen(argv[1], "r");
   if (fasta == NULL) exit_cannot_open(argv[1]);

   // Read and normalize genome
   char * genome = normalize_genome(fasta);

   SA_t  * SA  = create_SA(genome);
   BWT_t * BWT = create_BWT(genome, SA);
   Occ_t * Occ = create_Occ(BWT);
//   NKK_t * NKK = create_NKK(genome, SA, BWT, Occ);
           SA  = compress_SA(SA);

   // Write files
   char buff[256];
   char * data;
   ssize_t ws;
   size_t sz;

   // Write the suffix array file.
   sprintf(buff, "%s.sa", argv[1]);
   int fsar = creat(buff, 0644);
   if (fsar < 0) exit_cannot_open(buff);
   
   ws = 0;
   sz = sizeof(SA_t) + SA->nb * sizeof(int64_t);
   data = (char *) SA;
   while (ws < sz) ws += write(fsar, data + ws, sz - ws);
   close(fsar);


   // Write the Burrows-Wheeler transform.
   sprintf(buff, "%s.bwt", argv[1]);
   int fbwt = creat(buff, 0644);
   if (fbwt < 0) exit_cannot_open(buff);

   ws = 0;
   sz = sizeof(BWT_t) + BWT->nb * sizeof(uint8_t);
   data = (char *) BWT;
   while (ws < sz) ws += write(fbwt, data + ws, sz - ws);
   close(fbwt);


   // Write the Burrows-Wheeler transform.
   sprintf(buff, "%s.occ", argv[1]);
   int focc = creat(buff, 0644);
   if (focc < 0) exit_cannot_open(buff);

   ws = 0;
   sz = sizeof(Occ_t) + Occ->nb * AZ * sizeof(block_t);
   data = (char *) Occ;
   while (ws < sz) ws += write(focc, data + ws, sz - ws);
   close(focc);

   // Clean up.
   free(SA);
   free(BWT);
   free(Occ);

}
