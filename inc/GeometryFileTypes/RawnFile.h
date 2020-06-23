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
#if !defined(AFX_RAWNFILE_H__693A5AE3_831C_4120_A062_74EB8F04875E__INCLUDED_)
#define AFX_RAWNFILE_H__693A5AE3_831C_4120_A062_74EB8F04875E__INCLUDED_

#include <GeometryFileTypes/GeometryFileType.h>

class RawnFile : public GeometryFileType
{
	public:
		virtual ~RawnFile();
		virtual Geometry* loadFile(const string& fileName);
		virtual bool checkType(const string& fileName);
		virtual bool saveFile(const Geometry* geometry, const string& fileName);
		virtual string extension()
		{
			return "rawn";
		};
		virtual string filter()
		{
			return "Rawn files (*.rawn)";
		};
		static RawnFile ms_RawnFileRepresentative;
		static GeometryFileType* getRepresentative();
	protected:
		RawnFile();
};

#endif
