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
// This class reads slices from an memory resident block

#ifndef BLOCK_FAC_H
#define BLOCK_FAC_H

#include <C2C/SliceFactory.h>

template <class T>
class BlockFactory : virtual public SliceFactory<T>
{
	public:
		BlockFactory(T* _data, unsigned char* red, unsigned char* green, unsigned char* blue, int dim[3], float orig[3], float span[3]);
		virtual ~BlockFactory();

		virtual Slice<T>* getNextSlice();

		virtual void reset()
		{
			next = 0;
		}

		virtual int getNumSlices()
		{
			return m_dim[2];
		}

		virtual void getDimension(int dim[3], float orig[3], float span[3])
		{
			for(int i = 0; i < 3; i++)
			{
				dim[i] = m_dim[i];
				orig[i] = m_orig[i];
				span[i] = m_span[i];
			}
		}

		virtual int current()
		{
			return next;
		}

	private:

		int m_dim[3];
		float m_orig[3], m_span[3];
		int next;		   		// the number of the next available slice
		T* m_data;
		unsigned char* m_red, *m_green, *m_blue;
};

template <class T>
BlockFactory<T>::BlockFactory(T* _data, unsigned char* red,
							  unsigned char* green, unsigned char* blue,
							  int dim[3], float orig[3], float span[3])
{
	for(int i=0; i<3; i++)
	{
		m_dim[i] = dim[i];
		m_orig[i] = orig[i];
		m_span[i] = span[i];
	}
	m_data = _data;
	m_red = red;
	m_green = green;
	m_blue = blue;
	next = 0;
}

template <class T>
BlockFactory<T>::~BlockFactory()
{
}

template <class T>
Slice<T>* BlockFactory<T>::getNextSlice()
{
	Slice<T>* slc;
	int off;
	// are we past the end of the block?
	if(next >= m_dim[2])
	{
		return NULL;
	}
	// the new offset
	off = next*m_dim[0]*m_dim[1];
	// build a slice
	if(m_red && m_green && m_blue)
		slc = new Slice<T>(&m_data[off],
						   &m_red[off],&m_green[off],&m_blue[off], m_dim);
	else
	{
		slc = new Slice<T>(&m_data[off],NULL,NULL,NULL,m_dim);
	}
	// push next
	next++;
	return slc;
}

#endif
