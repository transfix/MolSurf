#include <C2C/arithlib/coder.h>
#include <C2C/arithlib/stats.h>
#include <C2C/arithlib/utils.h>
#include <C2C/arithlib/bitbuffer.h>
#include <C2C/arithlib/libdjvu++/ZPCodec.h>
#include <Utility/utility.h>

using namespace CVCUtility;

const int BitBuffer::INIT_BUF_SIZE = 10000;

BitBuffer::BitBuffer(int _bc)
{
	buf = NULL;
	size = 0;
	currbyte = 0;
	currbit = 0;
	getbyte = 0;
	getbit = 0;
	bits_context = _bc;
	assert(bits_context <= MAX_CONTEXT_BITS && bits_context >= MIN_CONTEXT_BITS);
}


BitBuffer::BitBuffer(char* bray, int nbyte, int nbit, int _bc)
{
	currbyte = nbyte;
	currbit = nbit;
	size = currbyte + 100;
	buf = (unsigned char*)malloc(size);
	memcpy(buf, bray, currbyte+((nbit > 0)? 1:0));
	getbyte = getbit = 0;
	bits_context = _bc;
	assert(bits_context <= MAX_CONTEXT_BITS && bits_context >= MIN_CONTEXT_BITS);
}

BitBuffer::~BitBuffer()
{
	if(buf != NULL)
	{
		free(buf);
	}
	buf = NULL;
}

void
BitBuffer::writeFile(FILE* ofp)
{
	fwrite(&currbyte, sizeof(int), 1, ofp);
	fwrite(&currbit, sizeof(char), 1, ofp);
	int nbits = currbyte + ((currbit > 0)? 1:0);
	if(nbits > 0)
	{
		fwrite(buf, 1, nbits, ofp);
	}
}

void
BitBuffer::readFile(FILE* ifp)
{
	put_reset();         // clean current bit buffer
	freadSafely(&currbyte, sizeof(int), 1, ifp);
	freadSafely(&currbit, sizeof(char), 1, ifp);
	size = currbyte + 10;
	buf = (u_char*)malloc(size);
	memset(buf, 0, size);
	int nbits = currbyte + ((currbit > 0)? 1:0);
	if(nbits > 0)
	{
		freadSafely(buf, 1, nbits, ifp);
	}
}

int
BitBuffer::getBuffer(unsigned char* ecd)
{
	int nb = getNumBytes() + sizeof(int) + sizeof(char);
	ecd = new unsigned char[nb];
	memcpy(ecd, &currbyte, sizeof(int));
	memcpy(ecd+sizeof(int),&currbit, sizeof(char));
	memcpy(ecd+sizeof(int)+sizeof(char), buf, nb-sizeof(int)-sizeof(char));
	return nb;
}

void
BitBuffer::put_a_float(int res, float f)
{
	if(res <= 4 || fabs(f) > 1.0)
	{
		printf("res = %d, f = %f\n", res, f);
		assert(res > 4 && fabs(f) <= 1.001);
	}
	BIT*  tmp = (BIT*) malloc(res*sizeof(BIT));
	// The first bit is the sign one:  0 -- Nonnegative; 1 -- Positive
	if(f>=0)
	{
		tmp[0] = 0;
	}
	else
	{
		f = -f;
		tmp[0] = 1;
	}
	float factor = 0.5;
	for(int k=1; k<res; k++)
	{
		if(f >= factor)
		{
			tmp[k] = 1;
			f = (f - factor);
		}
		else
		{
			tmp[k] = 0;
		}
		factor = (float)(factor/2.0);
	}
	put_bits(res, tmp);
	free(tmp);
}

void
BitBuffer::Encode_Positive(int n)
{
	assert(n >= 0);
	int  bits1, bits2;
	bits1 = Comp_Bits(n);
	bits2 = Comp_Bits(bits1-1);
	put_bits(3, bits2-1);
	put_bits(bits2, bits1-1);
	put_bits(bits1, n);
}

void
BitBuffer::put_a_bit(BIT b)
{
	if(currbyte == size)
	{
		int prev_size = size;
		if(buf == NULL)
		{
			buf = (u_char*)malloc(INIT_BUF_SIZE);
			size = INIT_BUF_SIZE;
		}
		else
		{
			buf = (u_char*)realloc(buf, size*2);
			size *= 2;
		}
		// set newly allocated memory to 0
		memset(buf+prev_size, 0, size-prev_size);
	}
	buf[currbyte] |= (b << (7-currbit));
	currbit++;
	if(currbit == 8)
	{
		currbit = 0;
		currbyte ++;
	}
}

/******
// put_bits(int n, int x): put last n bits of int x to bit buffer
// with the same order as in int
******/
void
BitBuffer::put_bits(int n, int x)
{
	BIT b;
#ifdef _DEBUG
	assert(n < 32);
#endif
	/*
	    for (int i = n-1; i >= 0; i--) {
	        b = 0x1 & (x >> i);
	        put_a_bit(b);
	    }
	  */
	if(currbyte + (n >> 3) + 1 >= size)
	{
		int prev_size = size;
		if(buf == NULL)
		{
			buf = (unsigned char*)malloc(INIT_BUF_SIZE);
			size = INIT_BUF_SIZE;
		}
		else
		{
			buf = (unsigned char*)realloc(buf, size*2);
			size *= 2;
		}
		memset(buf+prev_size, 0, size-prev_size);
	}
	for(int i = n-1; i >= 0; i--)
	{
		b = 0x1 & (x >> i);
		buf[currbyte] |= (b << (7-currbit));
		currbit ++;
		if(currbit == 8)
		{
			currbit = 0;
			currbyte++;
		}
	}
}

/******
// put a n bits array to bit buffer
******/
void
BitBuffer::put_bits(int n, BIT* bits)
{
	/*    for (int i = 0; i < n; i++) {
	        put_a_bit(bits[i]);
	    }
		*/
	if(currbyte + (n >> 3) + 1 >= size)
	{
		int prev_size = size;
		if(buf == NULL)
		{
			buf = (unsigned char*)malloc(INIT_BUF_SIZE);
			size = INIT_BUF_SIZE;
		}
		else
		{
			buf = (u_char*)realloc(buf, size*2);
			size *= 2;
		}
		memset(buf+prev_size, 0, size-prev_size);
	}
	for(int i = 0; i < n; i++)
	{
		buf[currbyte] |= (bits[i] << (7-currbit));
		currbit ++;
		if(currbit == 8)
		{
			currbit = 0;
			currbyte++;
		}
	}
}

/******
// put a signed int to bit buffer
******/
void
BitBuffer::put_sign_bits(int n, int x)
{
	BIT sign = (x >= 0)? 0:1;
	put_a_bit(sign);
	put_bits(n, abs(x));
}

/******
// put a float as array of char into bit buffer
******/
void
BitBuffer::put_a_float(float x)
{
	float ax = fabs(x);
	int power;
	float nx;
	if(x == 0)
	{
		power = 0;
		nx = 0;
	}
	else
	{
		double lgx = log10(ax);
		assert(lgx < 99);
		if(lgx < -99)    // treat the small number as 0
		{
			power = 0;
			nx = 0;
		}
		else
		{
			power = (int)ceil(lgx);
			nx = x / pow(10.0, power);
			assert(fabs(nx) <= 1.0);
		}
	}
	// put power into 8 bits
	put_sign_bits(7, power);
	// put mantissa into rest 24 bits
	put_a_float(24, nx);
	/*
	  float xx = x;
	  u_char *bytes = (u_char *)&xx;
	  for(int i = 0; i < 4; i++) {
	    int v = bytes[i];
	    put_bits(8, v);
	  }
	*/
}

void BitBuffer::append(int n, u_char* data)
{
	for(int i = 0; i < n; i++)
	{
		put_a_char(data[i]);
	}
}

void BitBuffer::copy_buffer(BitBuffer* thebuf)
{
	if(this == thebuf)
	{
		return;
	}
	if(buf != NULL)
	{
		free(buf);
	}
	size = thebuf->size;
	currbyte = thebuf->currbyte;
	currbit = thebuf->currbit;
	buf = (unsigned char*)malloc(size);
	memcpy(buf, thebuf->buf, size);
}

void
BitBuffer::arith_encode()
{
	// create a new bit buffer as the global buffer
	bit_buffer = new BitBuffer();
	//int mbytes = 10*DEFAULT_MEM;
	f_bits = DEFAULT_F;
	unsigned int non_null_contexts = 0;     // count of contexts used
	int i, cur_context;
	unsigned int    mask, bit, next_break, buffer;
	int         bits_to_go;
	binary_context**  contexts;
	binary_context*  still_coding;
	// initialise context array
	int ncontext = (1 << bits_context);
	contexts = (binary_context**)malloc(sizeof(binary_context*) * ncontext);
	if(contexts == NULL)
	{
		fprintf(stderr, "bits: not enough memory to allocate context array\n");
		exit(1);
	}
	still_coding = create_binary_context();
	// initialise contexts to NULL
	for(i=0; i < 1 << bits_context; i++)
	{
		contexts[i] = NULL;
	}
	// initalise variables
	cur_context = 0;
	mask = (1 << bits_context) - 1;
	next_break = BREAK_INTERVAL;
	start_encode();
	startoutputtingbits();
	// encode the first characters read for testing MAGIC_NO
	for(i=0; i <= currbyte && currbyte < size; i++)
	{
		buffer = buf[i];
		binary_encode(still_coding, 1);
		for(bits_to_go = 7; bits_to_go >= 0; bits_to_go--)
		{
			if(contexts[cur_context] == NULL)
			{
				contexts[cur_context] = create_binary_context();
				non_null_contexts++;
			}
			bit = (buffer >> bits_to_go) & 1;
			//fprintf(stderr, "%d", bit);
			binary_encode(contexts[cur_context], bit);
			cur_context = ((cur_context << 1) | bit) & mask;
		}
		if(next_break-- == 0)
		{
			finish_encode();
			start_encode();
			next_break = BREAK_INTERVAL;
		}
	}
	// encode end of message flag
	binary_encode(still_coding, 0);
	finish_encode();
	doneoutputtingbits();
	//
	// copy the compressed bit buffer
	//
	copy_buffer(bit_buffer);
	delete bit_buffer;
	for(i = 0; i < ncontext; i++)
	{
		if(contexts[i] != NULL)
		{
			free(contexts[i]);
		}
	}
	free(still_coding);
	free(contexts);
}

void
BitBuffer::arith_decode()
{
	// decode the file
	bit_buffer = new BitBuffer();
	bit_buffer->copy_buffer(this);
	//int nbyte = bit_buffer->currbyte;
	//char nbit = bit_buffer->currbit;
	// reset current buffer as empty
	put_reset();
	int i, cur_context, buffer, bits_to_go;
	int mask, bit, next_break;
	binary_context**  contexts;
	binary_context*  still_coding;
	unsigned int non_null_contexts = 0;     // count of contexts used
	f_bits = DEFAULT_F;
	max_frequency = 1<<f_bits;
	// initialise context array
	int ncontext = (1 << bits_context);
	contexts = (binary_context**)malloc(sizeof(binary_context*) * ncontext);
	if(contexts == NULL)
	{
		fprintf(stderr, "bits: unable to malloc %lu bytes\n",
				sizeof(binary_context*) *(1 << bits_context));
		exit(1);
	}
	still_coding = create_binary_context();
	// initialise contexts to NULL
	for(i=0; i < 1 << bits_context; i++)
	{
		contexts[i] = (binary_context*)NULL;
	}
	// initalise variables
	cur_context = 0;
	mask = (1 << bits_context) - 1;
	next_break = BREAK_INTERVAL;
	start_decode();
	startinputtingbits();
	//for(i = 0; i <= nbyte ; i++){
	while(binary_decode(still_coding) && getbyte <= (currbyte+4))
	{
		buffer = 0;
		for(bits_to_go = 7; bits_to_go >= 0; bits_to_go--)
		{
			if(contexts[cur_context] == (binary_context*)NULL)
			{
				contexts[cur_context] = create_binary_context();
				non_null_contexts++;
			}
			bit = binary_decode(contexts[cur_context]);
			buffer = (buffer << 1) | bit;
			cur_context = ((cur_context << 1) | bit) & mask;
		}
		put_bits(BYTE_SIZE, buffer);
		bytes_output++;
		if(next_break-- == 0)
		{
			finish_decode();
			start_decode();
			next_break = BREAK_INTERVAL;
		}
	}
	finish_decode();
	doneinputtingbits();
	delete bit_buffer;
	for(i = 0; i < ncontext; i++)
	{
		if(contexts[i] != NULL)
		{
			free(contexts[i]);
		}
	}
	free(still_coding);
	free(contexts);
}

/* ZP_Coder help functions */
static void encode_8_bits(ZPCodec& zp, int x, BitContext* ctx)
{
	int n = 1;
	for(int i = 0; i < 8; i++)
	{
		int b = (x&0x80)? 1:0;
		x <<= 1;
		zp.encoder(b, ctx[n-1]);
		n = (n << 1) | b;
	}
}

int decode_8_bits(ZPCodec& zp, BitContext* ctx)
{
	int n = 1;
	for(int i = 0; i < 8; i++)
	{
		n = (n << 1) | zp.decoder(ctx[n-1]);
	}
	return n & 0xff;
}

void BitBuffer::zp_encode()
{
	int nb = getNumBytes();
	MemoryByteStream out_bstream;
	out_bstream.write32(nb);
	BitContext ctx[256];
	memset(ctx, 0, sizeof(ctx));
	ZPCodec* p_zp = new ZPCodec(out_bstream, 1);
	for(int i = 0; i < nb; i++)
	{
		int x = buf[i];
		encode_8_bits(*p_zp, x, ctx);
	}
	// adjust the buffer
	delete p_zp;
	out_bstream.seek(0, SEEK_SET);
	int csize = out_bstream.readall(buf, size);
	currbyte = csize;
	currbit = 0;
	getbyte = 0;
	getbit = 0;
}

void BitBuffer::zp_decode()
{
	int nb = getNumBytes();
	unsigned char* in_buf = new unsigned char[nb];
	memcpy(in_buf, buf, nb);
	MemoryByteStream in_bstream(in_buf, nb);
	int orig_nb = in_bstream.read32();
	BitContext ctx[256];
	memset(ctx, 0, sizeof(ctx));
	ZPCodec zp(in_bstream, 0);
	put_reset();
	for(int i = 0; i < orig_nb; i++)
	{
		unsigned char c = (unsigned char) decode_8_bits(zp, ctx);
		put_a_char(c);
	}
	delete[] in_buf;
}

void BitBuffer::put_reset()
{
	if(buf != NULL)
	{
		free(buf);
		buf = NULL;
		size = 0;
	}
	currbyte = 0;
	currbit = 0;
	getbyte = 0;
	getbit = 0;
}

BIT
BitBuffer::get_a_bit()
{
	BIT bit = 0;
	bit = (buf[getbyte] >> (7-getbit)) & 0x1;
	getbit ++;
	if(getbit == 8)
	{
		getbyte ++;
		getbit = 0;
	}
	return bit;
}

void
BitBuffer::put_a_int(u_int x)
{
	put_bits(31, x);
	put_a_bit(0);
}

void
BitBuffer::put_a_short(u_short x)
{
	put_bits(16, x);
}

void
BitBuffer::put_a_char(u_char x)
{
	put_bits(8, x);
}

int
BitBuffer::get_a_int()
{
	int n = get_bits(31);
	// remove one more bit to match the word boundary
	get_a_bit();
	return n;
}

int
BitBuffer::get_a_short()
{
	int n = get_bits(16);
	return n;
}

int
BitBuffer::get_a_char()
{
	//  assert(getbyte <= currbyte);
	//printf("byte : %d, val = %d, \n", getbyte, buf[getbyte]);
	int n = buf[getbyte];
	getbyte++;
	return n;
}

float
BitBuffer::get_a_float()
{
	// get the sign of power
	int signp = get_a_bit();
	// get the amplitude of power from next 7 bits
	int power = get_bits(7);
	if(signp)
	{
		power = -power;
	}
	// get matissa of float from next 24 bits
	float x = get_a_float(24);
	x *= pow(10.0, power);
	/*
	  u_char bytes[4];

	  for(int i = 0; i < 4; i++) {
	    bytes[i] = get_bits(8);
	    //printf("%d", bytes[i]);
	  }

	  float x = *((float *)bytes);
	*/
	return x;
}

float
BitBuffer::get_a_float(int res)
{
	float x = 0;
	BIT sign;
	if(8*getbyte + getbit + res > 8*size)
	{
		fprintf(stderr, "Error in get_a_float(): Out of boundary of buffer\n");
		exit(1);
	}
	// The first is the sign bit;
	sign = get_a_bit();
	float  factor = 0.5;
	for(int k=1; k<res; k++)
	{
		x = x + factor*get_a_bit();
		factor = factor/2.0;
	}
	if(sign)
	{
		x = -x;
	}
	return x;
}
u_int
BitBuffer::get_bits(int n)
{
	// assert(n < 32 && n > 0);
	if(n <=0 || n >= 32)
	{
		fprintf(stderr, "get_bits: n = %d, should between (0, 32)\n", n);
		exit(1);
	}
	u_int x = 0;
	if(getbyte*8 + getbit + n > size*8)
	{
		fprintf(stderr, "Error in get_bits(): Out of boundary of buffer\n");
		fprintf(stderr, "getbyte = %d, getbit = %d, n = %d, size = %d\n",
				getbyte, getbit, n, size);
		exit(1);
	}
	for(int i = 0; i < n; i++)
	{
		x = 2 * x + get_a_bit();
	}
	return x;
}

int
BitBuffer::Decode_Positive()
{
	int bits1, bits2;
	bits2 = get_bits(3);
	bits2 = (bits2+1);
	bits1 = get_bits(bits2);
	bits1 = (bits1+1);
	return get_bits(bits1);
}

