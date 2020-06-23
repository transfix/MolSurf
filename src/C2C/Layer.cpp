/*
  Copyright 2011 The University of Texas at Austin

	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of MolSurf.

  MolSurf is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.


  MolSurf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with MolSurf; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <C2C/Layer.h>

Layer::Layer(int d1, int d2, u_char* cells)
{
	dim[0] = d1;
	dim[1] = d2;
	codes = (u_char*)malloc(sizeof(u_char)*(d1-1)*(d2-1));
	if(cells != NULL)
	{
		memcpy(codes, cells, (d1-1)*(d2-1));
	}
	else
	{
		memset(codes, 0, (d1-1)*(d2-1));
	}
	nc = 0;
	for(int i = 0; i < (d1-1)*(d2-1); i++)
	{
		if(cubeedges[codes[i]][0] > 0)
		{
			nc++;
		}
	}
}

Layer::Layer(const Layer& lay)
{
	dim[0] = lay.dim[0];
	dim[1] = lay.dim[1];
	nc = lay.nc;
	codes = (u_char*)malloc(sizeof(u_char)*(dim[0]-1)*(dim[1]-1));
	memcpy(codes, lay.codes, (dim[0]-1)*(dim[1]-1));
}

Layer::~Layer()
{
	if(codes != NULL)
	{
		free(codes);
	}
}

BIT* Layer::diffBits(Layer& lay)
{
	if(dim[0] != lay.dim[0] || dim[1] != lay.dim[1])
	{
		fprintf(stderr, "layer size doesn't match\n");
		return NULL;
	}
	BIT* bits = new BIT[(dim[0]-1)*(dim[1]-1)];
	for(int i = 0; i < dim[1]-1; i++)
		for(int j = 0; j < dim[0]-1; j++)
		{
			int n = j + i*(dim[0]-1);
			if((cubeedges[codes[n]][0] == 0 && cubeedges[lay.codes[n]][0] > 0) ||
					(cubeedges[codes[n]][0] > 0 && cubeedges[lay.codes[n]][0] == 0))
			{
				bits[n] = 1;
			}
			else
			{
				bits[n] = 0;
			}
		}
	return bits;
}

void Layer::writeOut(FILE* fp)
{
	fwrite(&nc, sizeof(int), 1, fp);
	printf("# of intersected cells = %d\n", nc);
	if(nc != 0)
	{
		BIT bit;
		BitBuffer* bbuf = new BitBuffer();
		for(int j = 0; j < dim[1]-1; j++)
		{
			for(int i = 0; i < dim[0]-1; i++)
			{
				int n = i + j*(dim[0]-1);
				bit = (cubeedges[codes[n]][0] > 0)? 1:0;
				bbuf->put_a_bit(bit);
			}
		}
#ifdef ZP_CODEC
		bbuf->zp_encode();
#else
		bbuf->arith_encode();
#endif
		bbuf->writeFile(fp);
		delete bbuf;
	}
}
