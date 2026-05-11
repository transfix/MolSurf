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
#ifndef __EXPANDABLEARRAY_H__
#define __EXPANDABLEARRAY_H__

#include <Utility/utility.h>
#include <PDBParser/IntQueue.h>

const unsigned int InitialArraySize = 16;

template <class T> class ExpandableArray {
public:
  ExpandableArray();
  virtual ~ExpandableArray();
  int add(T *object);
  bool set(T *object, unsigned int index);
  T *remove(unsigned int index);
  T *get(unsigned int index);
  T *getIth(unsigned int i) const;
  unsigned int getNumberOfObjects();
  void clear();
  void nonDestructiveClear();

protected:
  void initArrays();
  void initObjectArray();
  void initIndexArray();
  void deleteArrays();
  void doubleObjectArray();
  void doubleIndexArray();
  unsigned int m_NumberOfObjects;
  unsigned int m_SizeOfObjectsArray;
  T **m_T;
  int *m_ObjToIndex;
  unsigned int m_NextIndexEntry;
  unsigned int m_SizeOfIndexToObjectsArray;
  int *m_IndexToObj;
  Queue m_HoleList;
};

template <class T> ExpandableArray<T>::ExpandableArray() { initArrays(); }

template <class T> ExpandableArray<T>::~ExpandableArray() { deleteArrays(); }

template <class T> void ExpandableArray<T>::initArrays() {
  initObjectArray();
  initIndexArray();
}

template <class T> void ExpandableArray<T>::initObjectArray() {
  m_NumberOfObjects = 0;
  m_SizeOfObjectsArray = InitialArraySize;
  m_T = new T *[m_SizeOfObjectsArray];
  m_ObjToIndex = new int[m_SizeOfObjectsArray];
  unsigned int c;
  for (c = 0; c < m_SizeOfObjectsArray; c++) {
    m_T[c] = 0;
  }
}

template <class T> void ExpandableArray<T>::initIndexArray() {
  m_NextIndexEntry = 0;
  m_SizeOfIndexToObjectsArray = InitialArraySize;
  m_IndexToObj = new int[m_SizeOfIndexToObjectsArray];
}

template <class T> void ExpandableArray<T>::deleteArrays() {
  delete[] m_T;
  delete[] m_IndexToObj;
  delete[] m_ObjToIndex;
}

template <class T> void ExpandableArray<T>::doubleObjectArray() {
  unsigned int c;
  if (m_NumberOfObjects >= m_SizeOfObjectsArray) {
    T **oldT = m_T;
    int *oldObjToIndex = m_ObjToIndex;
    m_T = new T *[m_SizeOfObjectsArray * 2];
    m_ObjToIndex = new int[m_SizeOfObjectsArray * 2];
    for (c = 0; c < m_SizeOfObjectsArray; c++) {
      m_T[c] = oldT[c];
      m_ObjToIndex[c] = oldObjToIndex[c];
    }
    m_SizeOfObjectsArray *= 2;
    delete[] oldT;
    delete[] oldObjToIndex;
  }
}

template <class T> void ExpandableArray<T>::doubleIndexArray() {
  unsigned int c;
  if (m_NextIndexEntry >= m_SizeOfIndexToObjectsArray) {
    int *oldIndexToObjectsArray = m_IndexToObj;
    m_IndexToObj = new int[m_SizeOfIndexToObjectsArray * 2];
    for (c = 0; c < m_NextIndexEntry; c++) {
      m_IndexToObj[c] = oldIndexToObjectsArray[c];
    }
    m_SizeOfIndexToObjectsArray *= 2;
    delete[] oldIndexToObjectsArray;
  }
}

template <class T> int ExpandableArray<T>::add(T *t) {
  doubleObjectArray();
  m_T[m_NumberOfObjects] = t;
  int indexPosition;
  if (m_HoleList.isEmpty()) {
    doubleIndexArray();
    indexPosition = m_NextIndexEntry;
    m_NextIndexEntry++;
  } else {
    indexPosition = m_HoleList.deQueue();
  }
  m_ObjToIndex[m_NumberOfObjects] = indexPosition;
  m_IndexToObj[indexPosition] = m_NumberOfObjects;
  m_NumberOfObjects++;
  return m_NumberOfObjects - 1;
}

template <class T> bool ExpandableArray<T>::set(T *t, unsigned int index) {
  index = m_IndexToObj[index];
  if (index < m_NumberOfObjects && index >= 0) {
    m_T[index] = t;
    return true;
  }
  return false;
}

template <class T> T *ExpandableArray<T>::remove(unsigned int index) {
  int object = m_IndexToObj[index];
  if (object < (int)m_NumberOfObjects && object >= 0) {
    m_HoleList.enQueue(index);
    m_NumberOfObjects--;
    T *temp = m_T[object];
    m_T[object] = m_T[m_NumberOfObjects];
    m_ObjToIndex[object] = m_ObjToIndex[m_NumberOfObjects];
    m_IndexToObj[m_ObjToIndex[object]] = object;
    return temp;
  } else {
    return 0;
  }
}

template <class T> T *ExpandableArray<T>::get(unsigned int index) {
  int object = m_IndexToObj[index];
  if (object < (int)m_NumberOfObjects && object >= 0) {
    return m_T[object];
  } else {
    return 0;
  }
}

template <class T> T *ExpandableArray<T>::getIth(unsigned int i) const {
  if (i < m_NumberOfObjects && i >= 0) {
    return m_T[i];
  } else {
    return 0;
  }
}

template <class T> unsigned int ExpandableArray<T>::getNumberOfObjects() {
  return m_NumberOfObjects;
}

template <class T> void ExpandableArray<T>::clear() {
  unsigned int c;
  for (c = 0; c < m_NumberOfObjects; c++) {
    delete m_T[c];
    m_T[c] = 0;
  }
  m_NumberOfObjects = 0;
  m_NextIndexEntry = 0;
  m_HoleList.clearQueue();
}

template <class T> void ExpandableArray<T>::nonDestructiveClear() {
  unsigned int c;
  for (c = 0; c < m_NumberOfObjects; c++) {
    m_T[c] = 0;
  }
  m_NumberOfObjects = 0;
  m_NextIndexEntry = 0;
  m_HoleList.clearQueue();
}

#endif
