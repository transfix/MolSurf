#ifndef C2C_UTILS_H
#define C2C_UTILS_H
#include <assert.h>
#include <Utility/utility.h>

#ifdef WIN32
#include <C2C/arithlib/unitypes.h>
#else
#include <unistd.h>
#endif

float maximum(float* array, int);
float minimum(float* array, int);
float MyPower(int exp);
int Comp_Factor(float f);
int Comp_Bits(int maxv);
// Quantize a float number x to a int
// @param x a float number in [-1.0, 1.0]
// @param pow quantized int is in [0, pow)
int fquantize(float x, uint pow);
// Unquantize a int in [0, pow) to a float in [-1.0, 1.0]
float unfquantize(int n, uint pow);

#endif
