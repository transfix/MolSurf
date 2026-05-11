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
// An interface that encapsulates the details of disk file access
// and hides the small and big endian problem.
// You access the stream by using the overloaded get and put methods.

#ifndef XYZ_DISKIO_H
#define XYZ_DISKIO_H

#include <Utility/utility.h>

#define DBSIZE 4096       // disk block size
#define MAX_DB_NUM 300000 // maximum disk blocks difference in one seek
#define MAX_FN_LEN 512    // maximum filename length
#define REP_FAC 6         // replication factor
// #define _LITTLE_ENDIAN 1
#ifdef WIN32
typedef __int64 int64;
#else
typedef long long int64;
#endif                 // WIN32
typedef int64 Pointer; // Disk block pointer(index)
#ifdef WIN32
#define STRDUP(x) _strdup(x)
#else
#define STRDUP(x) strdup(x)
#endif

//  DiskIO abstract base interface
//  It hides the small and big endian problem. Real disk access classes are
//  subclass of this interface
class DiskIO {
public:
  enum Type { READ = 0, WRITE };
  enum SeekMode { FROM_HERE = 0, FROM_START, FROM_END };
  DiskIO(Type mode = READ) { m_mode = mode; }
  virtual ~DiskIO() {}
  virtual int get(char *data, int size) = 0;
  virtual int get(unsigned char *data, int size) = 0;
  virtual int get(short *data, int size) = 0;
  virtual int get(unsigned short *data, int size) = 0;
  virtual int get(int *data, int size) = 0;
  virtual int get(long *data, int size) = 0;
  virtual int get(int64 *data, int size) = 0;
  virtual int get(float *data, int size) = 0;
  virtual int get(double *data, int size) = 0;
  virtual void put(const char *data, int size) = 0;
  virtual void put(const unsigned char *data, int size) = 0;
  virtual void put(const short *data, int size) = 0;
  virtual void put(const unsigned short *data, int size) = 0;
  virtual void put(const int *data, int size) = 0;
  virtual void put(const long *data, int size) = 0;
  virtual void put(const int64 *data, int size) = 0;
  virtual void put(const float *data, int size) = 0;
  virtual void put(const double *data, int size) = 0;
  Type mode() const { return m_mode; }
  // change the file access type
  virtual bool setMode(DiskIO::Type type) { return false; }
  // open the disk IO for read or write
  virtual bool open() = 0;
  // open another file, but cannot change type
  virtual bool reopen(const char *fname) = 0;
  // reopen the current file
  // write files will be opened in append mode
  virtual bool reopen() { return false; }
  // close the disk IO
  //    @param fill: fill the rest of current disk block
  virtual bool close(bool fill = true) = 0;
  // seek
  virtual bool seek(long offset, SeekMode mode) = 0;
  // write out everything in buffer and possiblly move pointer to next available
  // disk block
  // @param fill fill the rest of current disk block
  virtual void flush(bool fill = true) = 0;

protected:
  Type m_mode;
};

#endif
