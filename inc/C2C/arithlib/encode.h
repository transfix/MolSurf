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
#ifndef C2C_ENCODE_H
#define C2C_ENCODE_H
#ifdef WIN32
#include <C2C/arithlib/unitypes.h>
#else
#include <unistd.h>
#endif
#include <C2C/arithlib/arith.h>
#include <C2C/arithlib/bitbuffer.h>

/** isoval is provided as central of array values. it
    can be set as 0
*/
/// encode a float array of length nv
BitBuffer* encode_vals(float* vals, int nv, float isoval = 0);

/// encode a unsigned short array of length nv
BitBuffer* encode_vals(u_short* vals, int nv, float isoval = 0);

/// encode a unsigned char array of length nv
BitBuffer* encode_vals(u_char* vals, int nv, float isoval = 0);

void en_second_diff(BitBuffer*, float* vals, int nv);

/// decode a float array of length nv
bool  decode_vals(BitBuffer*, float* vals, int nv, float isoval = 0);

/// decode a unsigned short array of length nv
bool  decode_vals(BitBuffer*, u_short* vals, int nv, float isoval = 0);

/// decode a unsigned char array of length nv
bool  decode_vals(BitBuffer*, u_char* vals, int nv, float isoval = 0);


#endif

