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
#ifndef ARITH_DEFINES_H
#define ARITH_DEFINES_H

#define MAX_CONTEXT_BITS	20    	/* max. number of bits for context */
#define MIN_CONTEXT_BITS	0	/* min. number of bits for context */
#define DEFAULT_BITS_CONTEXT   16	/* default value for bits_context */
#define ENCODE          	0
#define DECODE          	1
#define MAGICNO         	"123b"  /* Magic Number for files */
#define MAGICNO_LENGTH		4	/* length of magic number */
#define BREAK_INTERVAL		10000000	/* every round off to bit boundary */
#define MEGABYTE		(1<<24)	/* size of one megabyte */
#define NOMEMLEFT		-1	/* flag set when mem runs out */
#define DEFAULT_MEM		1	/* default 1 megabyte limit */
#define MIN_MBYTES        	1	/* minimum allowable memory size */
#define MAX_MBYTES        	255	/* maximum no for 8 bit int */


#define		CODE_BITS		32
#define		BYTE_SIZE		8
#define 	MAX_BITS_OUTSTANDING	256
#define 	HALF			((unsigned) 1 << (CODE_BITS-1))
#define 	QUARTER			(1 << (CODE_BITS-2))

#endif

