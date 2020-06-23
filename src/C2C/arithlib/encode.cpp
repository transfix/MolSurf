#ifdef WIN32
#include <windows.h>
#else
#endif

#include <C2C/arithlib/encode.h>
#include <C2C/arithlib/utils.h>
#include <Utility/utility.h>

/******
// encode an array of float numbers
// 1. put the sign respective to isovalue. it's necessary to
//    keep the topology after quantization error
// 2. encode the second order differences
******/
static int RelativeSign(float x, float base)
{
	if(x > base)
	{
		return 1;
	}
	return 0;
}
/*
BitBuffer *encode_vals(float *vals, int nv, float isoval)
{
	int i;
	BitBuffer *bbuf = new BitBuffer;

	float fmax = maximum(vals, nv);
	float fmin = minimum(vals, nv);
	float favg = (fmax + fmin) / 2;
	float frad = (fmax - fmin) / 2;

	//normailze the value array to [-1.0, 1.0]
	for(i = 0; i < nv; i++) {
		vals[i] = (vals[i] - favg) / frad;
	}
	isoval = (isoval - favg) / frad;
	printf("isoval = %f\n", isoval);

	// store the parameters
	bbuf->put_a_float(favg);
	bbuf->put_a_float(frad);

	u_int pow = (u_int)MyPower(DELTA_RES);
	int quan = fquantize(vals[0], pow);
	bbuf->put_bits(DELTA_RES, quan);

	int delta;
	BitBuffer *deltabuf = new BitBuffer;
	for(i = 1; i < nv; i++) {
		int q = fquantize(vals[i], pow);
		float x = unfquantize(q, pow);
		// adjust quantization to assure the quantize value has the same
		// relative sign as the original
		if(RelativeSign(vals[i], isoval) != RelativeSign(x, isoval)) {
			printf("sign flipped: %f -> %f, q = %d\n", vals[i], x, q);
			if(RelativeSign(vals[i], isoval) > 0) q += 1;
			else q -= 1;
		}
		delta = q - quan;
		if(delta < 0) bbuf->put_bits(DELTA_RES, delta+pow);
		bbuf->put_bits(DELTA_RES, delta);
		quan = q;
	}
#ifdef ZP_CODEC
    bbuf->zp_encode();
#else
    bbuf->arith_encode();
#endif
	//bbuf->append(deltabuf->getNumBytes(), deltabuf->getBits());
	delete deltabuf;
	return bbuf;
}
*/

BitBuffer* encode_vals(float* vals, int nv, float isoval)
{
	BitBuffer* bbuf = new BitBuffer;
	int i;
	////////
	// put the sign of float x relative to isovalue
	// 0: x >= isoval
	// 1: x < isoval
	////////
	for(i = 0; i < nv; i++)
	{
		BIT bit[2];
		if(vals[i] < isoval)
		{
			bit[0] = 1;
			bit[1] = 0;
		}
		else if(vals[i] > isoval)
		{
			bit[0] = 0;
			bit[1] = 1;
		}
		else
		{
			bit[0] = 0;
			bit[1] = 0;
		}
		bbuf->put_bits(1, bit);
		vals[i] = vals[i] - isoval;
	}
	// return max and min of vals array
	float fmax = maximum(vals, nv);
	float fmin = minimum(vals, nv);
	float favg = (fmax + fmin)/2;
	float frad = (fmax - fmin)/2;
	// scale the vals array to [-0.5, 0.5]
	for(i = 0; i < nv; i++)
	{
		vals[i] = (vals[i]-favg)/(2*frad);
	}
	///////////
	// encode the first float
	//////////
	bbuf->put_a_float(favg);
	bbuf->put_a_float(frad);
	bbuf->put_a_float(FLOAT_RES, vals[0]);
#ifdef DEBUG_AE
	cout << "max = " << fmax << ", min = " << fmin << endl;
	cout << "avg = " << favg << ", rad = " << frad << endl;
	cout << "average = " << favg << endl;
	cout << "vals[0] = " << vals[0] << endl;
#endif
	////////
	// compute the second order difference of input array
	// quantize them using DELTA_RES and put them into bbuf
	////////
	en_second_diff(bbuf, vals, nv);
#ifdef ZP_CODEC
	bbuf->zp_encode();
#else
	bbuf->arith_encode();
#endif
	return bbuf;
}

void en_second_diff(BitBuffer* buf, float* vals, int nv)
{
	int i;
#ifdef DEBUG_AE
	//all vals should be in [-0.5, 0.5]
	for(i = 0; i < nv; i++)
	{
		assert(fabs(vals[i]) <= 0.5);
	}
#endif
	if(nv <= 1)
	{
		return;
	}
	////////
	// put the deltas into bit buffer
	/////////
	float pow = MyPower(DELTA_RES);
	BitBuffer* tmpbuf = new BitBuffer();
	tmpbuf->put_a_float(FLOAT_RES, vals[0]);
	float xo = tmpbuf->get_a_float(FLOAT_RES);
	delete tmpbuf;
	float delta = vals[1] -xo;
	int quan = fquantize(delta, (u_int)pow);
	int prev = quan;
	buf->put_bits(DELTA_RES, quan);
	for(i = 2; i < nv; i++)
	{
		xo = xo + unfquantize(prev, (u_int)pow);
		delta = vals[i] - xo;
		quan = fquantize(delta, (u_int)pow);
		buf->put_bits(DELTA_RES, quan);
		prev = quan;
	}
}

BitBuffer* encode_vals(u_short* vals, int nv, float isoval)
{
	BitBuffer* bbuf = new BitBuffer;
	int res = USHORT_BIT;
	int delta;
	bbuf->put_bits(res, (int)vals[0]);
	// encode the second order differences
	for(int i = 1; i < nv; i++)
	{
		delta = vals[i] - vals[i-1];
		if(delta < 0)
		{
			bbuf->put_bits(res, delta + USHORT_POW);
		}
		else
		{
			bbuf->put_bits(res, delta);
		}
	}
#ifdef ZP_CODEC
	bbuf->zp_encode();
#else
	bbuf->arith_encode();
#endif
	return bbuf;
}

BitBuffer* encode_vals(u_char* vals, int nv, float isoval)
{
	BitBuffer* bbuf = new BitBuffer;
	int res = UCHAR_BIT;
	int delta;
	bbuf->put_bits(res, (int)vals[0]);
	// encode the second order differences
	for(int i = 1; i < nv; i++)
	{
		delta = vals[i] - vals[i-1];
		if(delta < 0)
		{
			delta += UCHAR_POW;
		}
		bbuf->put_bits(res, delta);
	}
#ifdef ZP_CODEC
	bbuf->zp_encode();
#else
	bbuf->arith_encode();
#endif
	/*for(int i = 0; i < nv; i++) {
	    bbuf->put_bits(res, vals[i]);
	}
	bbuf->arith_encode();
	*/
	return bbuf;
}

