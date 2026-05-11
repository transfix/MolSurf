#include <stdio.h>
#include <math.h>

main() {
  int i;

  i = 2147483644;

  while (i < i + 1) {
    i = i + 1;
  }
  printf("i = %d\n", i);
}
