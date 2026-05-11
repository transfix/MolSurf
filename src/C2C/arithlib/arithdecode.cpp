#include <math.h>
#include <iostream>
#include <C2C/arithlib/encode.h>
#include <C2C/arithlib/utils.h>

static inline int SIGN(float x) {
  if (x > 0) {
    return 1;
  }
  if (x < 0) {
    return 2;
  }
  return 0;
}

static int zero_count = 0;

bool decode_vals(BitBuffer *bbuf, float *vals, int nv, float isoval) {
  int i;
#ifdef ZP_CODEC
  bbuf->zp_decode();
#else
  bbuf->arith_decode();
#endif
  BIT *signs = new BIT[nv];
  for (i = 0; i < nv; i++) {
    signs[i] = 1 + bbuf->get_bits(1);
  }
  // get the vals
  float favg = bbuf->get_a_float();
  float frad = bbuf->get_a_float();
  // get the first float
  vals[0] = bbuf->get_a_float(FLOAT_RES);
  float *delta = new float[nv - 1];
  float delta_pow = MyPower(DELTA_RES);
  for (i = 0; i < nv - 1; i++) {
    int quan = bbuf->get_bits(DELTA_RES);
    delta[i] = unfquantize(quan, (int)delta_pow);
  }
  for (i = 1; i < nv; i++) {
    vals[i] = vals[i - 1] + delta[i - 1];
  }
  float shift = 0.5 / delta_pow;
  for (i = 0; i < nv; i++) {
    vals[i] = favg + 2 * vals[i] * frad;
    if (signs[i] != SIGN(vals[i])) {
      if (signs[i] == 1) {
        vals[i] = shift * 2 * frad;
      } else if (signs[i] == 2) {
        vals[i] = -shift * 2 * frad;
      } else {
        vals[i] = 0.0;
      }
    }
    vals[i] += isoval;
  }
  delete[] signs;
  delete[] delta;
  return true;
}

bool decode_vals(BitBuffer *bbuf, u_short *vals, int nv, float isoval) {
  int i;
#ifdef ZP_CODEC
  bbuf->zp_decode();
#else
  bbuf->arith_decode();
#endif
  int res = USHORT_BIT;
  // get the first u_short
  vals[0] = bbuf->get_bits(res);
  int delta;
  for (i = 1; i < nv; i++) {
    delta = bbuf->get_bits(res);
    if (vals[i - 1] + delta >= USHORT_POW) {
      delta -= USHORT_POW;
    }
    vals[i] = vals[i - 1] + delta;
  }
  return true;
}

bool decode_vals(BitBuffer *bbuf, u_char *vals, int nv, float isoval) {
  int i;
#ifdef ZP_CODEC
  bbuf->zp_decode();
#else
  bbuf->arith_decode();
#endif
  int res = UCHAR_BIT;
  // get the first u_char
  vals[0] = bbuf->get_bits(res);
  int delta;
  for (i = 1; i < nv; i++) {
    delta = bbuf->get_bits(res);
    if (vals[i - 1] + delta >= UCHAR_POW) {
      delta -= UCHAR_POW;
    }
    vals[i] = vals[i - 1] + delta;
  }
  return true;
}
