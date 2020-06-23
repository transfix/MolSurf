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
#include <SignDistanceFunction_v2/Geom3DParser.h>

using namespace CVCUtility;

Geom3DParser::Geom3DParser(void)
{
}

Geom3DParser::~Geom3DParser(void)
{
}

void Geom3DParser::ParseRawFile(FaceVertSet3D& fvs, const char* fname)
{

	FILE*	fp;
	int	nverts, ntris;
	float	temp[3];
	fp = fileRead(fname);
	bool isRawn = false;
	bool isRawc = false;

	string tmp = fname;
	if(fname[tmp.size()-1] == 'n') {
	  isRawn = true;
	  printf("Rawn input\n");
	}

	if(fname[tmp.size()-1] == 'c') {
	  isRawc = true;
	  printf("Rawn input\n");

	  if(fname[tmp.size()-2] == 'n') {
	    isRawn = true;
	    printf("Rawn input\n");
	  }
	}

	if(fscanf(fp,"%d %d", &nverts, &ntris) == EOF)
	{
		error("Input file is not valid:" + string(fname));
	}
	for(int i=0; i<nverts; i++)
	{
		if(fscanf(fp,"%f %f %f", &temp[0], &temp[1], &temp[2]) == EOF)
		{
			error("Error parsing vertices!");
		}
		fvs.addVert(temp[0], temp[1], temp[2]);

		if (isRawn) {
		  // junk the vertex normal
		  fscanf(fp,"%f %f %f", &temp[0], &temp[1], &temp[2]);
		}

		if (isRawc) {
		  // junk the vertex color
		  fscanf(fp,"%f %f %f", &temp[0], &temp[1], &temp[2]);
		}

	}
	int v1, v2, v3;
	for(int i=0; i<ntris; i++)
	{
		if(fscanf(fp,"%d %d %d", &v1, &v2, &v3) == EOF)
		{
			error("Error parsing triangles!");
		}
		TriId3i t(v1, v2, v3);
		fvs.AddTri(t);
	}
	fclose(fp);
	fvs.buildBBox();
	fvs.computeTriNormals();

	//TODO: add code to correctly orient all triangles.
	//fvs.flipTriNormals();
}
