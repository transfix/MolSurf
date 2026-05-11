#ifndef BIT_BUFFER_H
#define BIT_BUFFER_H
#include <stdio.h>
#include <math.h>
#include <C2C/arithlib/arith_defines.h>
#include <Utility/utility.h> // If you need to break this dependency later it should be really easy.
#include <C2C/arithlib/unitypes.h>

typedef unsigned char BIT;

/**
   BitBuffer: class to hold a buffer of bits
   BitBuffer supports various add and remove operations,
   especially it supports arithmatic encodeing and decoding
   of the buffered bits.
*/
class BitBuffer {
public:
  /// constructor
  BitBuffer(int _bc = DEFAULT_BITS_CONTEXT);
  /// construct bitbuffer with a bit array
  BitBuffer(char *bits, int nbyte, int nbit, int _bc = DEFAULT_BITS_CONTEXT);
  /// destructor
  ~BitBuffer();

  /// curByte(): accessor of current byte
  int curByte() const { return currbyte; }
  /// curBit(): accessor of current bit in current byte
  int curBit() const { return currbit; }
  /**
   *     getNumBytes(): return number of bytes in the buffer
   */
  int getNumBytes() const {
    int nbits = currbyte + ((currbit > 0) ? 1 : 0);
    return nbits;
  }

  /// writeFile(): dump bit buffer to file
  void writeFile(FILE *ofp);
  /// readFile(): read bit biffer from file
  void readFile(FILE *ifp);

  /// return the raw bit array
  unsigned char *getBits() const { return buf; }
  /** get the raw bit array with current byte and bit num info.
      return the number of bytes in the array
  */
  int getBuffer(unsigned char *ecd);

  /// add a bit to buffer
  void put_a_bit(BIT);

  /// put_bits(n, x): put last n bits of int x to bit buffer
  void put_bits(int n, int x);

  /// put_bits(n, bits): put a n bits array to bit buffer
  void put_bits(int n, BIT *bits);

  /// put_sign_bits(n, x): put a signed bit and last n bits of x to bitbuffer
  void put_sign_bits(int n, int x);

  /// put_a_float(x): put float x as array of char into bit buffer
  void put_a_float(float x);

  /// put a singed float x to buffer with resolution res
  void put_a_float(int res, float x);

  /// add an unsigned int to buffer
  void put_a_int(u_int);

  /// add an unsigned short to buffer
  void put_a_short(u_short);

  /// add an unsigned char to buffer
  void put_a_char(u_char);

  /// add a positive integer to buffer in encoded form to save space
  void Encode_Positive(int);

  /// put_reset(): reset current buffer as empty
  void put_reset();

  /// append(): append a char array to current buffer
  void append(int n, u_char *data);

  /// get a bit from buffer at current position
  BIT get_a_bit();

  /// get n bits from buffer
  u_int get_bits(int n);

  /// get a char from buffer
  int get_a_char();

  /// get a unsigned int from buffer
  int get_a_int();

  /// get a unsigned short from buffer
  int get_a_short();

  /// get a float from buffer
  float get_a_float();

  /// get a float with resolution res
  float get_a_float(int res);

  /// reverse of Encode_Positive(), decode a int from buffer
  int Decode_Positive();

  /// reset current get pointer to begin of buffer
  void get_reset() {
    getbyte = 0;
    getbit = 0;
  }

  /// arith_encode(): arithmetic encode current bit buffer
  void arith_encode();

  /// arith_decode(): arithmetic decode current bit buffer
  void arith_decode();

  /// encode the bit buffer with the ZP_Coder from http://djvu.research.att.com
  void zp_encode();

  /// decode the bitbuffer with the ZP_Coder
  void zp_decode();

protected:
  /// make buffer itself as a copy of the argument buffer
  void copy_buffer(BitBuffer *);

private:
  unsigned char *buf;
  int size;
  int currbyte;
  char currbit;
  int getbyte;
  char getbit;
  int bits_context; // num bits in binary contexts

  const static int INIT_BUF_SIZE;
};

#endif
