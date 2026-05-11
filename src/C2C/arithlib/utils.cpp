#include <Utility/utility.h>
#include <C2C/arithlib/utils.h>

float MyPower(int exp) {
  float tmp = 1;
  for (int k = 0; k < exp; k++) {
    tmp *= 2;
  }
  return tmp;
}

float maximum(float *array, int n) {
  assert(n >= 1);
  float max = array[0];
  for (int i = 1; i < n; i++) {
    if (array[i] > max) {
      max = array[i];
    }
  }
  return max;
}

float minimum(float *array, int n) {
  assert(n >= 1);
  float min = array[0];
  for (int i = 1; i < n; i++) {
    if (array[i] < min) {
      min = array[i];
    }
  }
  return min;
}

int Comp_Factor(float f) {
  int ret = 0;
  if (f < 0.5) {
    return 0;
  }
  while (f >= 0.5) {
    f = (float)(f / 2.0);
    ret++;
  }
  return ret;
}

int Comp_Bits(int maxv) {
  int tmp = 2, ret = 1;
  assert(maxv >= 0);
  if (maxv < 0) {
    fprintf(stderr, "Comp_Bits(): argument should be positive\n");
    exit(1);
  }
  while (tmp <= maxv) {
    tmp = 2 * tmp;
    ret++;
  }
  return ret;
}

int fquantize(float x, uint pow) {
  if (x >= 1.0) {
    x = 1.0;
  } else if (x <= -1.0) {
    x = -1.0;
  }
  float y = (x + 1) / 2;
  if (floor(y * pow) >= pow - 1) {
    return pow - 1;
  }
  return (int)floor(y * pow) + (((y * pow - floor(y * pow)) > 0.5) ? 1 : 0);
}

float unfquantize(int n, uint pow) {
  return (float)(2 * (n + 0.5) / pow - 1.0);
}
