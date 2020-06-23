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
#if !defined(AFX_RAWCFILE_H__F47E5524_8A5C_43D2_AFD5_8E6E5ACCB23C__INCLUDED_)
#define AFX_RAWCFILE_H__F47E5524_8A5C_43D2_AFD5_8E6E5ACCB23C__INCLUDED_

#include <GeometryFileTypes/GeometryFileType.h>

class RawcFile : public GeometryFileType
{
	public:
		virtual ~RawcFile();
		virtual Geometry* loadFile(const string& fileName);
		virtual bool checkType(const string& fileName);
		virtual bool saveFile(const Geometry* geometry, const string& fileName);
		virtual string extension()
		{
			return "rawc";
		};
		virtual string filter()
		{
			return "Rawc files (*.rawc)";
		};

		static RawcFile ms_RawcFileRepresentative;
		static GeometryFileType* getRepresentative();
	protected:
		RawcFile();

};

#endif
