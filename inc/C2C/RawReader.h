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
// rawreader.h: this class reads slices of raw file and put them into cache

#ifndef C2C_RAW_READER_H
#define C2C_RAW_READER_H

#include <Utility/utility.h>
#include <C2C/DiskIO.h>
#include <C2C/SliceCache.h>

template <class T>
class RawReader
{

	public:
		// constructor
		RawReader(DiskIO* io, int dim[3]) : m_io(io)
		{
			for(int i = 0; i < 3; i++)
			{
				m_dim[i] = dim[i];
			}
			next = 0;
		}

		// destructor
		~RawReader() {}

		// get the nth slice
		// return NULL if n is greater than the last slice number
		Slice<T>* getSlice(int n);

		// slice size
		int size()
		{
			return m_dim[0]*m_dim[1]*sizeof(T);
		}

	private:
		DiskIO* m_io;
		int m_dim[3];
		int next; // current slice
};

template <class T>
Slice<T>* RawReader<T>::getSlice(int n)
{
	if(n >= m_dim[2])
	{
		return NULL;
	}
	if(next != n)
	{
		m_io->seek((n-next)*size(), DiskIO::FROM_HERE);
		next = n;
	}
	T* data = new T[m_dim[0]*m_dim[1]];
	m_io->get(data, m_dim[0]*m_dim[1]);
	Slice<T>* slc = new Slice<T>(data, NULL,NULL,NULL, m_dim);
	next++;
	delete[] data;
	return slc;
}

#endif
