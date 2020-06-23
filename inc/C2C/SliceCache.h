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
// slicecache.h: a cache class where we store slices.
// It's implemented using a FIFO algorithm.
// We use a mutex with slice cache to manage synchronization of accessing cache

#if !defined(C2C_SLICE_CACHE_H)
#define C2C_SLICE_CACHE_H

#include <Utility/utility.h>
#include <C2C/Slice.h>

#ifdef MUL_THREAD
#include <cmutex.h>
#endif

// This class implement a cache for slices
template <class T>
class SliceCache
{
	public:
		// constructor, default size 1M Bytes
		// default cache size 1M bytes
		SliceCache(int nbytes = 1000000):cap(nbytes)
		{
			size = 0;
#ifdef MUL_THREAD
			lock = new Mutex();
#endif
		}
		// destructor
		~SliceCache();

		// check if the cache contains a slice with key
		bool contains(int key);
		// add a new slice to cache with key
		int add(int key, Slice<T>* slice);
		// add slice to cache if there is enough space. return true if the slice is added
		bool tryAdd(int key, Slice<T>* slice);
		// remove slice with key from cache, return size of space released
		int remove(int key);
		//  remove the first page in cache
		int remove_first();
		// get a page from cache. return NULL if page not exist
		Slice<T>* get(int key);

		// return the capacity of cache
		int capacity() const
		{
			return cap;
		}
	private:
		bool enoughSpace(int nbytes);
		map<int, Slice<T>*> pages;
		list<int> keys;
		int cap;                               // capacity of cache in bytes
		int size;                              // size of cache in use
#ifdef MUL_THREAD
		Mutex* lock;
#endif
};

template <class T>
SliceCache<T>::~SliceCache()
{
#ifdef MUL_THREAD
	lock->acquire();
#endif
	typename map<int, Slice<T>*>::iterator i;
	for(i = pages.begin(); i != pages.end(); ++i)
	{
		Slice<T>* slice = (*i).second;
		delete slice;
	}
#ifdef MUL_THREAD
	lock->release();
	delete lock;
#endif
}

// contains(): return true if slice with key is in cache
template <class T>
bool SliceCache<T>::contains(int key)
{
#ifdef MUL_THREAD
	lock->acquire();
#endif
	bool existed = (pages.find(key) != pages.end());
#ifdef MUL_THREAD
	lock->release();
#endif
	return existed;
}

// remove(): remove the slice with key from cache, return size of space released
template <class T>
int SliceCache<T>::remove(int key)
{
	// simply return if no page with key exists
	if(!contains(key))
	{
		return 0;
	}
#ifdef MUL_THREAD
	lock->acquire();
#endif
	typename map<int, Slice<T>*>::iterator i = pages.find(key);
	keys.remove(key);            // remove key from key list
	// remove  pair from pages
	Slice<T>* slice = (*i).second;
	int nbytes = slice->size();
	size -= nbytes;
	delete slice;
	pages.erase(i);
#ifdef MUL_THREAD
	lock->release();
#endif
	return nbytes;
}

// remove_first(): remove the first page in cache
template <class T>
int SliceCache<T>::remove_first()
{
#ifdef MUL_THREAD
	lock->acquire();
#endif
	if(keys.size() == 0)
	{
#ifdef MUL_THREAD
		lock->release();
#endif
		return 0;
	}
	int first = keys.front();
	keys.pop_front();
	typename map<int, Slice<T>*>::iterator i = pages.find(first);
	Slice<T>* slice = (*i).second;
	int nbytes = slice->size();
	delete slice;
	size -= nbytes;
	pages.erase(i);
#ifdef MUL_THREAD
	lock->release();
#endif
	return nbytes;
}

// enoughSpace(): if there is enough space to hold nbytes
template <class T>
bool SliceCache<T>::enoughSpace(int nbytes)
{
#ifdef MUL_THREAD
	lock->acquire();
#endif
	bool enough;
	// note: cache must be able to hold at least one slice
	// so we always return true if no slice in cache
	if(size + nbytes <= cap || keys.size() == 0)
	{
		enough = true;
	}
	else
	{
		enough = false;
	}
#ifdef MUL_THREAD
	lock->release();
#endif
	return enough;
}

// add(): add a new slice with key to cache, return size of space added
template <class T>
int SliceCache<T>::add(int key, Slice<T>* slice)
{
	// return if there is a page with same key in cache
	if(contains(key))
	{
		return 0;
	}
	// remove old pages to release enough space in cache
	while(!enoughSpace(slice->size()))
	{
		remove_first();
	}
#ifdef MUL_THREAD
	lock->acquire();
#endif
	pair<int, Slice<T>*> imgpair(key, slice);
	pages.insert(imgpair);
	keys.push_back(key);
	size += slice->size();
#ifdef MUL_THREAD
	lock->release();
#endif
	return slice->size();
}

// tryAdd(): add a new slice to cache if enough space exists.
// return true if page is added. Otherwise return false;
template <class T>
bool SliceCache<T>::tryAdd(int key, Slice<T>* slice)
{
	if(!enoughSpace(slice->size()))
	{
		return false;
	}
	add(key, slice);
	return true;
}

// get(): get a page from cache. return NULL if page not exist
template <class T>
Slice<T>* SliceCache<T>::get(int key)
{
	if(contains(key))  // page in cache
	{
#ifdef MUL_THREAD
		lock->acquire();
#endif
		typename map<int, Slice<T>*>::iterator i = pages.find(key);
		Slice<T>* slice = (*i).second;
		Slice<T>* slc = new Slice<T>(*slice);
#ifdef MUL_THREAD
		lock->release();
#endif
		return slc;
	}
	else
	{
		return NULL;
	}
}

#endif
