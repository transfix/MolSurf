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
#include <Utility/utility.h>
#include <GeometryFileTypes/RawcFile.h>
#include <Geometry/Geometry.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace CVCUtility;

RawcFile RawcFile::ms_RawcFileRepresentative;

RawcFile::RawcFile()
{
}

RawcFile::~RawcFile()
{
}

Geometry* RawcFile::loadFile(const string& fileName)
{
	bool zeroSeen = false, maxSeen = false;
	struct stat fst;
	if(stat(fileName.c_str(), &fst) == -1)
	{
		error("File " + fileName + " does not exist!");
	}
	FILE* fp = fopen(fileName.c_str(), "r");
	if(!fp)
	{
		error("Error opening file " + fileName + "!");
	}
	int numverts, numpolygons;
	if(2!=fscanf(fp, "%d %d", &numverts, &numpolygons))
	{
		error("Error reading in number of verts and polygons\n");
	}
	// make sure the number of verts & polygons are positive
	if(numverts<0 || numpolygons<0)
	{
		error("Negative number of verts or polygons\n");
	}
	vector<double> vertPositions;
	vector<double> vertColors;
	int c;
	// read in the verts
	for(c=0; c<numverts; c++)
	{
		float x, y, z, r, g, b;
		// read in a single vert
		if(6!=fscanf(fp, "%f %f %f %f %f %f\n", &x, &y, &z, &r, &g, &b))
		{
			error("Error reading in vert");
		}
		vertPositions.push_back(x);
		vertPositions.push_back(y);
		vertPositions.push_back(z);
		vertColors.push_back(r);
		vertColors.push_back(g);
		vertColors.push_back(b);
	}
	vector<unsigned int> pointIndices;
	vector<unsigned int> lineIndices;
	vector<unsigned int> triangleIndices;
	// read in the polygons
	for(c=0; c<numpolygons; c++)
	{
		char polyIndices[1024];
		//fGetString(polyIndices, 1023, fp);
		fgetsSafely(polyIndices, 1023, fp);
		vector<unsigned int> indexList;
		{
			char* pch;
			pch = strtok(polyIndices," \t\n");
			while(pch != NULL)
			{
				unsigned int idx = atoi(pch);
				if(idx == 0)
				{
					zeroSeen = true;
				}
				if(idx == numverts)
				{
					maxSeen = true;
				}
				if(zeroSeen && maxSeen)
				{
					error("Indices were both 0 and maxNumVerts!\n");
				}
				if(idx > numverts)
				{
					error("Index out of bounds\n");
				}
				indexList.push_back(idx);
				pch = strtok(NULL, " \t\n");
			}
		}
		if(indexList.size() == 0)
		{
		}
		else if(indexList.size() == 1)
		{
			pointIndices.push_back(indexList[0]);
		}
		else if(indexList.size() == 2)
		{
			lineIndices.push_back(indexList[0]);
			lineIndices.push_back(indexList[1]);
		}
		else
		{
			int n = indexList.size();
			int i;
			unsigned int firstIdx = indexList[0];
			unsigned int prevIdx = indexList[1];
			for(i=2; i<n; i++)
			{
				triangleIndices.push_back(firstIdx);
				triangleIndices.push_back(prevIdx);
				triangleIndices.push_back(indexList[i]);
				prevIdx = indexList[i];
			}
		}
	}
	if(maxSeen)
	{
		int i;
		for(i=0; i<pointIndices.size();    i++)
		{
			pointIndices[i]--;
		}
		for(i=0; i<lineIndices.size();     i++)
		{
			lineIndices[i]--;
		}
		for(i=0; i<triangleIndices.size(); i++)
		{
			triangleIndices[i]--;
		}
	}
	if((vertPositions.size()==0) || ((pointIndices.size() == 0) && (lineIndices.size() == 0) && (triangleIndices.size() == 0)))
	{
		return 0;
	}
	fclose(fp);
	// initialize the geometry
	Geometry* geometry = new Geometry;
	if(pointIndices.size())
	{
		geometry->AllocatePoints(pointIndices.size());
		geometry->AllocatePointColors();
	}
	if(lineIndices.size())
	{
		geometry->AllocateLines(numverts, lineIndices.size() / 2);
		geometry->AllocateLineColors();
	}
	if(triangleIndices.size())
	{
		geometry->AllocateTris(numverts, triangleIndices.size() / 3);
		geometry->AllocateTriVertColors();
	}
	// copy vertices, unfortunately, need to duplicate
	if(pointIndices.size())
	{
		int i;
		for(i=0; i<pointIndices.size(); i++)
		{
			unsigned int pIdx = pointIndices[i];
			geometry->m_Points[3*i+0] = vertPositions[3*pIdx+0];
			geometry->m_Points[3*i+1] = vertPositions[3*pIdx+1];
			geometry->m_Points[3*i+2] = vertPositions[3*pIdx+2];
			geometry->m_PointColorsArray[3*i+0] = vertColors[3*pIdx+0];
			geometry->m_PointColorsArray[3*i+1] = vertColors[3*pIdx+1];
			geometry->m_PointColorsArray[3*i+2] = vertColors[3*pIdx+2];
		}
	}
	if(lineIndices.size())
	{
		int i;
		for(i=0; i<vertPositions.size(); i++)
		{
			geometry->m_LineVerts[i] = vertPositions[i];
		}
		for(i=0; i<vertColors.size(); i++)
		{
			geometry->m_LineColors[i] = vertColors[i];
		}
		for(i=0; i<lineIndices.size(); i++)
		{
			geometry->m_Lines[i] = lineIndices[i];
		}
	}
	if(triangleIndices.size())
	{
		int i;
		for(i=0; i<vertPositions.size(); i++)
		{
			geometry->m_TriVerts[i] = vertPositions[i];
		}
		for(i=0; i<geometry->m_NumTriVerts; i++)
		{
			geometry->m_TriVertColorsTransparent[4*i+0] = vertColors[3*i+0];
			geometry->m_TriVertColorsTransparent[4*i+1] = vertColors[3*i+1];
			geometry->m_TriVertColorsTransparent[4*i+2] = vertColors[3*i+2];
		}
		for(i=0; i<triangleIndices.size(); i++)
		{
			geometry->m_Tris[i] = triangleIndices[i];
		}
	}
	return geometry;
}

bool RawcFile::checkType(const string& fileName)
{
	bool zeroSeen = false, maxSeen = false;
	float f1,f2,f3,f4,f5,f6;
	unsigned int u1,u2,u3;
	struct stat fst;
	if(stat(fileName.c_str(), &fst) == -1)
	{
		printf("File does not exist\n");
		return 0;
	}
	FILE* fp = fopen(fileName.c_str(), "r");
	if(!fp)
	{
		printf("Error opening file\n");
		return false;
	}
	// get the number of verts and triangles
	int numverts, numtris;
	if(2!=fscanf(fp, "%d %d", &numverts, &numtris))
	{
		printf("Error reading in number of verts and tris\n");
		fclose(fp);
		return false;
	}
	// make sure the number of verts & tris are positive
	if(numverts<0 || numtris<0)
	{
		printf("Negative number of verts or tris\n");
		fclose(fp);
		return false;
	}
	int c;
	// read in the verts
	for(c=0; c<numverts; c++)
	{
		// read in a single vert, which includes
		// position and normal
		if(6!=fscanf(fp, "%f %f %f %f %f %f",
					 &(f1),
					 &(f2),
					 &(f3),
					 &(f4),
					 &(f5),
					 &(f6)))
		{
			printf("Error reading in vert # %d\n", c);
			fclose(fp);
			return false;
		}
	}
	// read in the triangles
	for(c=0; c<numtris; c++)
	{
		// read in 3 integers for each triangle
		if(3!=fscanf(fp, "%u %u %u",
					 &(u1),
					 &(u2),
					 &(u3)))
		{
			printf("Error reading in tri # %d\n", c);
			fclose(fp);
		}
		// the file might start indexing verts from 1 or 0
		// check if indexes go up to the num of verts or if they
		// start from 0
		if(u1==0 || u2==0 || u3==0)
		{
			zeroSeen = true;
		}
		if(u1==(unsigned int)numverts || u2==(unsigned int)numverts || u3==(unsigned int)numverts)
		{
			maxSeen = true;
		}
		// cant have both!
		if(maxSeen && zeroSeen)
		{
			printf("Found 0 & max in tri # %d\n", c);
			fclose(fp);
			return false;
		}
		// check the bounds on each vert
		if(u1>(unsigned int)numverts || u2>(unsigned int)numverts || u3>(unsigned int)numverts)
		{
			printf("Bounds error reading in tri # %d\n", c);
			fclose(fp);
			return false;
		}
	}
	fclose(fp);
	return true;
}

bool RawcFile::saveFile(const Geometry* geometry, const string& fileName)
{
	FILE* fp;
	// open the file
	fp = fopen(fileName.c_str(), "w");
	if(!fp)
	{
		printf("Error opening the output file\n");
		return false;
	}
	// write the number of verts & tris
	if(0>=fprintf(fp, "%d %d\n", geometry->m_NumTriVerts, geometry->m_NumTris))
	{
		printf("Error writing the number of verts and tris\n");
		return false;
	}
	unsigned int c;
	// write out the verts
	for(c=0; c<geometry->m_NumTriVerts; c++)
	{
		if(geometry->m_TriVertColorsTransparent)
		{
			if(0>=fprintf(fp, "%f %f %f %f %f %f\n",
						  (geometry->m_TriVerts[c*3+0]),
						  (geometry->m_TriVerts[c*3+1]),
						  (geometry->m_TriVerts[c*3+2]),
						  (geometry->m_TriVertColorsTransparent[c*4+0]),
						  (geometry->m_TriVertColorsTransparent[c*4+1]),
						  (geometry->m_TriVertColorsTransparent[c*4+2])))
			{
				printf("Error writing out vert # %d\n", c);
				fclose(fp);
				return false;
			}
		}
		else
		{
			if(0>=fprintf(fp, "%f %f %f %f %f %f\n",
						  (geometry->m_TriVerts[c*3+0]),
						  (geometry->m_TriVerts[c*3+1]),
						  (geometry->m_TriVerts[c*3+2]),
						  (1.0),
						  (1.0),
						  (1.0)))
			{
				printf("Error writing out vert # %d\n", c);
				fclose(fp);
				return false;
			}
		}
	}
	// write out the tris
	for(c=0; c<geometry->m_NumTris; c++)
	{
		if(0>=fprintf(fp, "%d %d %d\n",
					  (geometry->m_Tris[c*3+0]),
					  (geometry->m_Tris[c*3+1]),
					  (geometry->m_Tris[c*3+2])))
		{
			printf("Error writing out tri # %d\n", c);
			fclose(fp);
			return false;
		}
	}
	fclose(fp);
	return true;
}

GeometryFileType* RawcFile::getRepresentative()
{
	return &ms_RawcFileRepresentative;
}
