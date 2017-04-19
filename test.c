#include <stdio.h>
#include <stdlib.h>

struct tmp_t {
  char _[4];
};

int main(void) {
  struct tmp_t x = {0};
  fprintf(stderr, "%d\n", x._[3]);
}
