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
#if !defined(AFX_MAYAOBJFILE_H__1A1D22F0_B170_4AD3_8BDA_35CF0F0C54B4__INCLUDED_)
#define AFX_MAYAOBJFILE_H__1A1D22F0_B170_4AD3_8BDA_35CF0F0C54B4__INCLUDED_

#include <GeometryFileTypes/GeometryFileType.h>

class MayaOBJFile : public GeometryFileType {
public:
  virtual ~MayaOBJFile();
  virtual Geometry *loadFile(const string &fileName);
  virtual bool checkType(const string &fileName);
  virtual bool saveFile(const Geometry *geometry, const string &fileName);
  virtual string extension() { return "obj"; };
  virtual string filter() { return "Maya OBJ files (*.obj)"; };
  static MayaOBJFile ms_MayaOBJFileRepresentative;
  static GeometryFileType *getRepresentative();

protected:
  MayaOBJFile();
};

#endif
