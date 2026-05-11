/*
  Copyright 2000-2003 The University of Texas at Austin

        Authors: Xiaoyu Zhang 2000-2002 <xiaoyu@ices.utexas.edu>
                                         John Wiggins 2003 <prok@cs.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of iotree.

  iotree is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  iotree is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with iotree; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/******************************************************************************
File:		stats.h

Authors: 	John Carpenelli   (johnfc@ecr.mu.oz.au)
                Wayne Salamonsen  (wbs@mundil.cs.mu.oz.au)

Purpose:	Data compression using a word-based model and revised
                arithmetic coding method.

Based on: 	A. Moffat, R. Neal, I.H. Witten, "Arithmetic Coding Revisited",
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
#ifndef STATS_H
#define STATS_H

/*
 * macros to add and remove the end '1' bit of a binary number
 * using two's complement arithmetic
 */
#define BACK(i) ((i) & ((i) - 1))
#define FORW(i) ((i) + ((i) & -(i)))

#define NOT_KNOWN -1        /* attempt to code unknown symbol */
#define TOO_MANY_SYMBOLS -1 /* could not install symbol */
#define NO_MEMORY -2        /* install exceeded memory */

#define STATIC 0      /* context cannot grow- no escape */
#define DYNAMIC 1     /* context may grow- escape needed */
#define GROWTH_RATE 2 /* rate at which context grows */

/* memory used per symbol is a function of the GROWTH_RATE */
#define MEM_PER_SYMBOL 4 * sizeof(int)

#define DEFAULT_F 27  /* default value of f_bits */
#define MIN_F_BITS 10 /* minimum no of bits for f bits */
#define MAX_F_BITS 30 /* maximum no of bits for f bits */
#define MIN_INCR 1    /* minimum increment value */

/* context structure used to store frequencies */
typedef struct {
  int initial_size;       /* original length of context */
  int max_length, length; /* length of tree and current length */
  int nSingletons;        /* no. symbols with frequency=1 */
  int type;               /* context may be STATIC or DYNAMIC */
  int nSymbols;           /* count of installed symbols */
  unsigned long total;    /* total of all frequencies */
  unsigned long *tree;    /* Fenwick's binary index tree */
  unsigned long incr;     /* current increment */
} context;

/* context structure for binary contexts */
typedef struct {
  int c0;   /* number of zeroes */
  int c1;   /* number of ones */
  int incr; /* current increment used */
} binary_context;

/* provide external linkage to stats variables */
extern int f_bits;                 /* number of frequency bits */
extern unsigned int max_frequency; /* maximum total frequency count */

/* function prototypes */
context *create_context(int length, int type);
int install_symbol(context *pTree, int symbol);
int encode(context *pContext, int symbol);
int decode(context *pContext);
void get_interval(context *pContext, int *pLow, int *pHigh, int symbol);
void purge_context(context *pContext);
void halve_context(context *pContext);
binary_context *create_binary_context(void);
int binary_encode(binary_context *pContext, int bit);
int binary_decode(binary_context *pContext);

#endif
