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
// This class reads slices from a rawv file

#ifndef RAWV_SLICE_FAC_H
#define RAWV_SLICE_FAC_H

#include <Utility/utility.h>
#include <C2C/BufferedIO.h>
#include <C2C/RawVReader.h>
#include <C2C/SliceCache.h>
#include <C2C/SliceFactory.h>

template <class T>
class RawVSliceFactory : virtual public SliceFactory<T>
{
	public:
		RawVSliceFactory(const char* rawfile);

		virtual ~RawVSliceFactory();

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
		void readRawVHeader();
		void loadInitPages();
		int m_dim[3];
		float m_orig[3], m_span[3];
		unsigned char m_types[4];
		int m_nvars, m_timesteps, m_datastart;
		int next; // the number of the next available slice
		SliceCache<T>* cache;
		BufferedIO* m_io;
		RawVReader<T>* reader;
};

// Alas, template functions must stay in the header.

template <class T>
RawVSliceFactory<T>::RawVSliceFactory(const char* rawfile)
{
	m_io = new BufferedIO(rawfile, DiskIO::READ);
	if(!m_io->open())
	{
		printf("unable to open volume data file.\n");
		exit(-1);
	}
	readRawVHeader();
	// make sure the color vars are kosher
	if(m_types[0] != 1 || m_types[1] != 1 || m_types[2] != 1)
	{
		printf("Error: bad rawv file. The color variables should be first and/or");
		printf(" unsigned char in type.\n");
		exit(-1);
	}
	reader = new RawVReader<T>(m_io, m_dim, m_types, m_timesteps, m_datastart);
	cache = new SliceCache<T>(m_dim[0]*m_dim[1]*m_dim[2]*(sizeof(T)+3));
	next = 0;
	loadInitPages();
}

template <class T>
RawVSliceFactory<T>::~RawVSliceFactory()
{
	delete cache;
	delete reader;
	m_io->close();
	delete m_io;
}

template <class T>
Slice<T>* RawVSliceFactory<T>::getNextSlice()
{
	if(cache->contains(next))
	{
		Slice<T>* slc = cache->get(next);
		next++;
		return slc;
	}
	return reader->getSlice(next++);
}

template <class T>
void RawVSliceFactory<T>::readRawVHeader()
{
	int dummy, i;
	// retrieve header info
	m_io->get(&dummy, 1);
	if(dummy != 0xBAADBEEF)
	{
		printf("error: rawv magic value not present. aborting.\n");
		exit(-1);
	}
	m_io->get(m_dim, 3);
	m_io->get(&m_timesteps, 1);
	m_io->get(&m_nvars, 1);
	m_io->get(m_orig, 3); // min
	m_io->seek(4, DiskIO::FROM_HERE); // min timestep
	m_io->get(m_span, 3); // max
	m_io->seek(4, DiskIO::FROM_HERE); // max timestep
	// calculate the span from the min and max extents and dimension
	for(i=0; i < 3; i++)
	{
		m_span[i] = (m_span[i]-m_orig[i]) / (float) m_dim[i];
	}
	for(i=0; i < m_nvars; i++)
	{
		if(i < 4)
		{
			m_io->get(&m_types[i], 1);
		}
		else
		{
			m_io->seek(1, DiskIO::FROM_HERE);
		}
		m_io->seek(64, DiskIO::FROM_HERE);
	}
	// calculate the offset of the data portion
	m_datastart = 56+65*m_nvars;
}

template <class T>
void RawVSliceFactory<T>::loadInitPages()
{
	int sz = m_dim[0]*m_dim[1]*sizeof(T);
	int np = (cache->capacity() < sz)? 1:cache->capacity()/sz;
	for(int i = 0; i < np && i < m_dim[2]; i++)
	{
		Slice<T>* slc = reader->getSlice(i);
		cache->add(i, slc);
	}
}

#endif
