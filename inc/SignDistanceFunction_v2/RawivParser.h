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
#if !defined(AFX_RAWIVPASER_H__938E0E03_F02E_4298_B690_AD4DE7EF3AB4__INCLUDED_)
#define AFX_RAWIVPASER_H__938E0E03_F02E_4298_B690_AD4DE7EF3AB4__INCLUDED_

#include <Utility/utility.h>
#include <SignDistanceFunction_v2/Reg3Parser.h>
#include <SignDistanceFunction_v2/reg3data.h>
#include <SignDistanceFunction_v2/bio.h>

class RawivParser : public Reg3Parser
{
	public:
		RawivParser();
		virtual ~RawivParser();
		virtual bool parse(Reg3Data<float>* data, const char* fname);
		virtual bool write(const Reg3Data<float>& data, const char* fname);
	private:
		bool isRawivFile(const char* fname)
		{
			int len = (int)strlen(fname);
			return (len > 6 && strcmp(fname+len-6, ".rawiv") == 0);
		}
};

#endif
