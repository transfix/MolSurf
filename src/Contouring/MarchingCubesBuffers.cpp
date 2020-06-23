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
#include <Contouring/MarchingCubesBuffers.h>
#include <Contouring/cubes.h>

using CCVContouring::MarchingCubesBuffers;

MarchingCubesBuffers::MarchingCubesBuffers()
{
	setDefaults();
}

MarchingCubesBuffers::~MarchingCubesBuffers()
{
	destroyEdgeBuffers();
}

void MarchingCubesBuffers::setDefaults()
{
	// set default values for all member variables
	m_EdgeCaches[0] = 0;
	m_EdgeCaches[1] = 0;
	m_EdgeCaches[2] = 0;
	m_EdgeCaches[3] = 0;
	m_EdgeCaches[4] = 0;
	m_EdgeCaches[5] = 0;
	m_VertClassifications[0] = 0;
	m_VertClassifications[1] = 0;
	m_AmountAllocated = 0;
}
bool MarchingCubesBuffers::allocateEdgeBuffers(unsigned int width, unsigned int height)
{
	// only allocate memory if our current buffer is not big enough
	if(width*height>m_AmountAllocated)
	{
		destroyEdgeBuffers();
		return forceAllocateEdgeBuffers(width, height);
	}
	else
	{
		return true;
	}
}

bool MarchingCubesBuffers::forceAllocateEdgeBuffers(unsigned int width, unsigned int height)
{
	// allocate the edge buffer without checking to see if we already have a big enough buffer
	m_EdgeCaches[0] = new unsigned int[width*height];
	m_EdgeCaches[1] = new unsigned int[width*height];
	m_EdgeCaches[2] = new unsigned int[width*height];
	m_EdgeCaches[3] = new unsigned int[width*height];
	m_EdgeCaches[4] = new unsigned int[width*height];
	m_VertClassifications[0] = new unsigned int[width*height];
	m_VertClassifications[1] = new unsigned int[width*height];
	if(m_EdgeCaches[0]&&
			m_EdgeCaches[1]&&
			m_EdgeCaches[2]&&
			m_EdgeCaches[3]&&
			m_EdgeCaches[4]&&
			m_VertClassifications[0]&&
			m_VertClassifications[0])
	{
		m_AmountAllocated = width*height;
		return true;
	}
	else
	{
		destroyEdgeBuffers();
		return false;
	}
}

void MarchingCubesBuffers::destroyEdgeBuffers()
{
	// free the edge buffer memory
	unsigned int c;
	for(c=0; c<5; c++)
	{
		delete [] m_EdgeCaches[c];
		m_EdgeCaches[c] = 0;
	}
	delete [] m_VertClassifications[0];
	m_VertClassifications[0] = 0;
	delete [] m_VertClassifications[1];
	m_VertClassifications[1] = 0;
	m_AmountAllocated = 0;
}

void MarchingCubesBuffers::swapEdgeBuffers()
{
	// swap the edges buffers
	unsigned int* temp;
	temp = m_EdgeCaches[XEdgesBack];
	m_EdgeCaches[XEdgesBack] = m_EdgeCaches[XEdgesFront];
	m_EdgeCaches[XEdgesFront] = temp;
	temp = m_EdgeCaches[YEdgesBack];
	m_EdgeCaches[YEdgesBack] = m_EdgeCaches[YEdgesFront];
	m_EdgeCaches[YEdgesFront] = temp;
	temp = m_VertClassifications[0];
	m_VertClassifications[0] = m_VertClassifications[1];
	m_VertClassifications[1] = temp;
}
