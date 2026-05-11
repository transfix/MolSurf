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
File:		coder.h

Authors: 	John Carpinelli   (johnfc@ecr.mu.oz.au)
                Wayne Salamonsen  (wbs@mundil.cs.mu.oz.au)

Purpose:	Data compression using a word-based model and revised
                arithmetic coding method.

Based on: 	A. Moffat, R. Neal, I.H. Witten, "Arithmetic Coding Revisted",
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
#ifndef CODER_H
#define CODER_H

#include <C2C/arithlib/arith_defines.h>

#ifndef BIT_BUFFER_H
class BitBuffer;
#endif

/* provide external linkage to variables */
extern int f_bits;               /* link to f_bits in stats.c */
extern unsigned int bytes_input; /* make available to other modules */
extern unsigned int bytes_output;

extern BitBuffer *bit_buffer; // output buffer

/* function prototypes */
void arithmetic_encode(unsigned int l, unsigned int h, unsigned int t);
unsigned int arithmetic_decode_target(unsigned int t);
void arithmetic_decode(unsigned int l, unsigned int h, unsigned int t);
void binary_arithmetic_encode(int c0, int c1, int bit);
int binary_arithmetic_decode(int c0, int c1);
void start_encode(void);
void finish_encode(void);
void start_decode(void);
void finish_decode(void);
void startoutputtingbits(void);
void doneoutputtingbits(void);
void startinputtingbits(void);
void doneinputtingbits(void);

#endif
