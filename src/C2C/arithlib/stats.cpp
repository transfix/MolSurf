/******************************************************************************
File: 		stats.c

Authors: 	John Carpinelli   (johnfc@ecr.mu.oz.au)
	 	Wayne Salamonsen  (wbs@mundil.cs.mu.oz.au)

Purpose:	Data compression using a word-based model and revised
		arithmetic coding method.

Based on: 	P.M. Fenwick, "A new data structure for cumulative probability
		tables", Software- Practice and Experience, 24:327-336,
		March 1994.

		A. Moffat, R. Neal, I.H. Witten, "Arithmetic Coding Revisited",
		Proc. IEEE Data Compression Conference, Snowbird, Utah,
		March 1995.


Copyright 1995 John Carpinelli and Wayne Salamonsen, All Rights Reserved.

These programs are supplied free of charge for research purposes only,
and may not sold or incorporated into any commercial product.  There is
ABSOLUTELY NO WARRANTY of any sort, nor any undertaking that they are
fit for ANY PURPOSE WHATSOEVER.  Use them at your own risk.  If you do
happen to find a bug, or have modifications to suggest, please report
the same to Alistair Moffat, alistair@cs.mu.oz.au.  The copyright
notice above and this statement of conditions must remain an integral
part of each and every copy made of these files.

******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <C2C/arithlib/stats.h>
#include <C2C/arithlib/coder.h>


/* variables local to stats module */
int f_bits = DEFAULT_F;				/* bits in frequency counts */
unsigned int max_frequency = 1 << DEFAULT_F;	/* max. total frequency */


/*
 *
 * create a new frequency table using a binary index tree
 * table may be STATIC or DYNAMIC depending on the type parameter
 * DYNAMIC tables may grow above their intial length as new symbols
 * are installed
 *
 */
context
*create_context(int length, int type)
{
	context*	pContext;
	int		i;
	int		size = 1;
	/*
	 * increment length to accommodate the fact
	 * that symbol 0 is stored at pos 1 in the array.
	 */
	length++;
	/* round length up to next power of two */
	while(size < length)
	{
		size = size << 1;
	}
	/* malloc context structure and array for frequencies */
	if(((pContext = (context*) malloc(sizeof(context))) == NULL) ||
			((pContext->tree = (unsigned long*) malloc((size+1)*sizeof(long)))
			 == NULL))
	{
		fprintf(stderr, "stats: not enough memory to create context\n");
		exit(1);
	}
	pContext->initial_size = size;	/* save for purging later */
	pContext->length = 0;		/* current no. of symbols */
	pContext->total = 0;		/* total frequency */
	pContext->nSymbols = 0;		/* count of symbols set to zero */
	pContext->type = type;		/* is context DYNAMIC or STATIC */
	pContext->max_length = size;	/* no. symbols before growing */
	/* initialise contents of tree array to zero */
	for(i = 0; i < size; i++)
	{
		pContext->tree[i] = 0;
	}
	pContext->incr = 1 << f_bits;	/* increment is initially 2 ^ f */
	if(type  == DYNAMIC)
	{
		pContext->nSingletons = pContext->incr;
	}
	else
	{
		pContext->nSingletons = 0;
	}
	return pContext;	    		/* return a pointer to the context */
}


/*
 *
 * install a new symbol in a context's frequency table
 * returns 0 if successful, TOO_MANY_SYMBOLS or NO_MEMORY if install fails
 *
 */
int
install_symbol(context* pContext, int symbol)
{
	int i;
	symbol++;	/* increment because first symbol at position one */
	/*
	 * if new symbol is greater than current array length then double length
	 * of array
	 */
	while(symbol >= pContext->max_length)
	{
		pContext->tree = (unsigned long*)
						 realloc(pContext->tree,
								 pContext->max_length * GROWTH_RATE * sizeof(long));
		if(pContext->tree == NULL)
		{
			fprintf(stderr, "stats: not enough memory to expand context\n");
			return NO_MEMORY;
		}
		/* clear new part of table to zero */
		for(i=pContext->max_length; i<GROWTH_RATE*pContext->max_length; i++)
		{
			pContext->tree[i] = 0;
		}
		/*
		 * initialize new part by setting first element of top half
		 * to total of bottom half
		 * this method depends on table length being a power of two
		 */
		pContext->tree[pContext->max_length] = pContext->total;
		pContext->max_length *= GROWTH_RATE;
	}
	/* check that we are not installing too many symbols */
	if(((pContext->nSymbols + 1) << 1) >= (int) max_frequency)
		/*
		 * cannot install another symbol as all frequencies will
		 * halve to one and an infinite loop will result
		 */
	{
		return TOO_MANY_SYMBOLS;
	}
	if(symbol > pContext->length)  	/* update length if necessary */
	{
		pContext->length = symbol;
	}
	pContext->nSymbols++;		/* increment count of symbols */
	i = symbol;	    			/* update elements in tree below */
	do
	{
		pContext->tree[i] += pContext->incr;
		i = FORW(i);
	}
	while(i < pContext->max_length);
	/* update the number of singletons if an context is DYNAMIC */
	if(pContext->type == DYNAMIC)
	{
		pContext->nSingletons += pContext->incr;
	}
	pContext->total += pContext->incr;			/* update total */
	/* halve frequency counts if total greater than max_frequency */
	while(pContext->total+pContext->nSingletons > max_frequency)
	{
		halve_context(pContext);
	}
	return 0;
}



/*
 *
 * encode a symbol given its context
 * the lower and upper bounds are determined using the frequency table,
 * and then passed on to the coder
 * if the symbol has zero frequency, code an escape symbol and
 * return NOT_KNOWN otherwise returns 0
 *
 */
int
encode(context* pContext, int symbol)
{
	int low, high=0;
	symbol++;
	if((symbol > 0) && (symbol < pContext->max_length))
	{
		get_interval(pContext, &low, &high, symbol);
	}
	else
	{
		low = high;
	}
	if(low == high)
	{
		if(pContext->nSingletons == 0)
		{
			fprintf(stderr,
					"stats: cannot code zero-probability novel symbol");
			exit(1);
		}
		/* encode the escape symbol if unknown symbol */
		arithmetic_encode(pContext->total, pContext->total+
						  pContext->nSingletons,
						  pContext->total+pContext->nSingletons);
		return NOT_KNOWN;
	}
	/* call the coder with the low, high and total for this symbol */
	arithmetic_encode(low, high, pContext->total+pContext->nSingletons);
	/* update the singleton count if symbol was previously a singleton */
	if(pContext->type == DYNAMIC)
		if(high-low == (int) pContext->incr)
		{
			pContext->nSingletons -= pContext->incr;
		}
	/* increment the symbol's frequency count */
	while(symbol<pContext->max_length)
	{
		pContext->tree[symbol] += pContext->incr;
		symbol = FORW(symbol);
	}
	pContext->total += pContext->incr;
	while(pContext->total+pContext->nSingletons > max_frequency)
	{
		halve_context(pContext);
	}
	return 0;
}




/*
 *
 * decode function is passed a context, and returns a symbol
 *
 */
int
decode(context* pContext)
{
	int	mid, symbol, i, target;
	int low, high;
	target = arithmetic_decode_target(pContext->total+pContext->nSingletons);
	/* check if the escape symbol has been received */
	if(target >= (int) pContext->total)
	{
		arithmetic_decode(pContext->total,
						  pContext->total+pContext->nSingletons,
						  pContext->total+pContext->nSingletons);
		return NOT_KNOWN;
	}
	symbol = 0;
	mid = pContext->max_length / 2;		/* midpoint is half length */
	/* determine symbol from target value */
	while(mid > 0)
	{
		if((int) pContext->tree[symbol+mid] <= target)
		{
			symbol = symbol+mid;
			target = target-pContext->tree[symbol];
		}
		mid /= 2;
	}
	/*
	 * pass in symbol and symbol+1 instead of symbol-1 and symbol to
	 * account for array starting at 1 not 0
	 */
	i = symbol+1;
	get_interval(pContext, &low, &high, i);
	arithmetic_decode(low, high, pContext->total+pContext->nSingletons);
	/* update the singleton count if symbol was previously a singleton */
	if(pContext->type == DYNAMIC)
		if(high-low == (int) pContext->incr)
		{
			pContext->nSingletons -= pContext->incr;
		}
	/* increment the symbol's frequency count */
	pContext->tree[i] += pContext->incr;
	i = FORW(i);
	while(i<pContext->max_length)
	{
		pContext->tree[i] += pContext->incr;
		i = FORW(i);
	}
	pContext->total += pContext->incr;
	/* halve all frequencies if necessary */
	while(pContext->total+pContext->nSingletons > max_frequency)
	{
		halve_context(pContext);
	}
	return symbol;
}



/*
 *
 * get the low and high limits of the frequency interval
 * occupied by a symbol.
 * this function is faster than calculating the upper bound of the two
 * symbols individually as it exploits the shared parents of s and s-1.
 *
 */
void
get_interval(context* pContext, int* pLow, int* pHigh, int symbol)
{
	int low, high, shared, parent;
	/* calculate first part of high path */
	high = pContext->tree[symbol];
	parent = BACK(symbol);
	/* calculate first part of low path */
	symbol--;
	low = 0;
	while(symbol != parent)
	{
		low += pContext->tree[symbol];
		symbol = BACK(symbol);
	}
	/* sum the shared part of the path back to root */
	shared = 0;
	while(symbol > 0)
	{
		shared += pContext->tree[symbol];
		symbol = BACK(symbol);
	}
	*pLow = shared+low;
	*pHigh = shared+high;
}


/*
 *
 * halve_context is responsible for halving all the frequency counts in a
 * context.
 * halves context in linear time by keeping track of the old and new
 * values of certain parts of the array
 * also recalculates the number of singletons in the new halved context.
 *
 */

void
halve_context(context* pContext)
{
	int	old_values[MAX_F_BITS], new_values[MAX_F_BITS];
	int	i, zero_count, temp, sum_old, sum_new;
	pContext->incr = (pContext->incr + MIN_INCR) >> 1;	/* halve increment */
	pContext->nSingletons = pContext->incr;
	for(i = 1; i < pContext->max_length; i++)
	{
		temp = i;
		/* work out position to store element in old and new values arrays */
		for(zero_count = 0; !(temp&1); temp >>= 1)
		{
			zero_count++;
		}
		/* move through context halving as you go */
		old_values[zero_count] = pContext->tree[i];
		for(temp = zero_count-1, sum_old = 0, sum_new = 0; temp >=0; temp--)
		{
			sum_old += old_values[temp];
			sum_new += new_values[temp];
		}
		pContext->tree[i] -= sum_old;
		pContext->total -= (pContext->tree[i]>>1);
		pContext->tree[i] -= (pContext->tree[i]>>1);
		if(pContext->tree[i] == pContext->incr)
		{
			pContext->nSingletons += pContext->incr;
		}
		pContext->tree[i] += sum_new;
		new_values[zero_count] = pContext->tree[i];
	}
	if(pContext->type == STATIC)
	{
		pContext->nSingletons = 0;
	}
}


/*
 *
 * free memory allocated for a context and initialize empty context
 * of original size
 *
 */
void
purge_context(context* pContext)
{
	int i;
	free(pContext->tree);
	/* malloc new tree of original size */
	if((pContext->tree = (unsigned long*)malloc((pContext->initial_size + 1)
						 * sizeof(long))) == NULL)
	{
		fprintf(stderr, "stats: not enough memory to create context\n");
		exit(1);
	}
	pContext->length = 0;
	pContext->total = 0;
	pContext->nSymbols = 0;
	pContext->max_length = pContext->initial_size;
	for(i = 0; i < pContext->initial_size; i++)
	{
		pContext->tree[i] = 0;
	}
	pContext->incr = 1 << f_bits;   	/* increment is initially 2 ^ f */
	if(pContext->type  == DYNAMIC)
	{
		pContext->nSingletons = pContext->incr;
	}
	else
	{
		pContext->nSingletons = 0;
	}
}

/******************************************************************************
*
* functions for binary contexts
*
******************************************************************************/


/*
 *
 * create a binary_context for binary contexts
 * contexts consists of two counts and an increment which
 * is normalized
 *
 */
binary_context* create_binary_context(void)
{
	binary_context* pContext;
	pContext = (binary_context*) malloc(sizeof(binary_context));
	if(pContext == NULL)
	{
		fprintf(stderr, "stats: not enough memory to create context\n");
		exit(1);
	}
	pContext->incr = 1 << (f_bits - 1);		/* start with incr=2^(f-1) */
	pContext->c0 = pContext->incr;
	pContext->c1 = pContext->incr;
	return pContext;
}



/*
 *
 * encode a binary symbol using special binary arithmetic
 * coding functions
 * returns 0 if successful
 *
 */
int
binary_encode(binary_context* pContext, int bit)
{
	binary_arithmetic_encode(pContext->c0, pContext->c1, bit);
	/* increment symbol count */
	if(bit == 0)
	{
		pContext->c0 += pContext->incr;
	}
	else
	{
		pContext->c1 += pContext->incr;
	}
	/* halve frequencies if necessary */
	if(pContext->c0 + pContext->c1 >= (int) max_frequency)
	{
		pContext->c0 = (pContext->c0 + 1) >> 1;
		pContext->c1 = (pContext->c1 + 1) >> 1;
		pContext->incr = (pContext->incr + MIN_INCR) >> 1;
	}
	return 0;
}



/*
 *
 * decode a binary symbol using specialised binary arithmetic
 * coding functions
 *
 */
int
binary_decode(binary_context* pContext)
{
	int bit;
	bit = binary_arithmetic_decode(pContext->c0, pContext->c1);
	/* increment symbol count */
	if(bit == 0)
	{
		pContext->c0 += pContext->incr;
	}
	else
	{
		pContext->c1 += pContext->incr;
	}
	/* halve frequencies if necessary */
	if(pContext->c0 + pContext->c1 >= (int) max_frequency)
	{
		pContext->c0 = (pContext->c0 + 1) >> 1;
		pContext->c1 = (pContext->c1 + 1) >> 1;
		pContext->incr = (pContext->incr + MIN_INCR) >> 1;
	}
	return bit;
}

