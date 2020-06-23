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
#include <Contouring/RGBAExtractor.h>

using CCVContouring::RGBAExtractor;
using CCVContouring::ContourGeometry;

RGBAExtractor::RGBAExtractor()
{
	setDefaults();
}

RGBAExtractor::~RGBAExtractor()
{
	setDefaults();
}

void RGBAExtractor::extractContour(ContourGeometry* contourGeometry, float isovalue, float R, float G, float B) const
{
	contourGeometry->allocateVertexBuffers(m_Width*m_Height);
	contourGeometry->allocateTriangleBuffers(m_Width*m_Height);
	contourGeometry->setUseColors(true);
	contourGeometry->setIsovalue(isovalue);
	unsigned int numCellsX = m_Width-1;
	unsigned int numCellsY = m_Height-1;
	unsigned int numCellsZ = m_Depth-1;
	bool cacheAvailable[12];
	unsigned int offsetTable[8];
	unsigned int withinSliceOffsetTable[8];
	unsigned int cacheOffsetTable[12];
	buildOffsetTable(offsetTable, m_Width, m_Height);
	buildWithinSliceOffsetTable(withinSliceOffsetTable, m_Width);
	buildEdgeCacheOffsetTable(cacheOffsetTable, m_Width);
	unsigned int cubeCase;
	unsigned int i,j,k, voxelOffset, offsetInSlice, newEdge, edgeCounter, triangleCounter, edgeID;
	// the new vertex
	unsigned int v1xuint, v1yuint, v1zuint, v2xuint, v2yuint, v2zuint, vertOffset1, vertOffset2;
	GLfloat v1x, v1y, v1z, v2x, v2y, v2z, den1, den2;
	GLfloat n1x, n1y, n1z, n2x, n2y, n2z;
	GLfloat c1x, c1y, c1z, c2x, c2y, c2z;
	// edges involved in a triangle
	unsigned int e1,e2,e3;
	unsigned int rowStart, inSliceRowStart;
	classifyVertices(0, m_Buffers->m_VertClassifications[0], isovalue);
	for(k=0; k<numCellsZ; k++)    // for every slice
	{
		classifyVertices(k+1, m_Buffers->m_VertClassifications[1], isovalue);
		for(j=0; j<numCellsY; j++)    // for every row
		{
			rowStart = m_Width*m_Height*k+m_Width*j;
			inSliceRowStart = m_Width*j;
			for(i=0; i<numCellsX; i++)    // for every cell
			{
				if(i==0 || i==1 || i==numCellsX-1)
				{
					// isCached only changes at the start and end of a scan line
					// as well as the second voxel of every scan line
					isCached(cacheAvailable, i,j,k,m_Width,m_Height);
				}
				voxelOffset = rowStart+i;
				offsetInSlice = inSliceRowStart+i;
				cubeCase = determineCase(m_Data, withinSliceOffsetTable, offsetInSlice, isovalue);
				if(cubeCase!=0 && cubeCase!=255)
				{
					// for each edge involved
					for(edgeCounter=0; edgeCounter<cubeedges[cubeCase][0]; edgeCounter++)
					{
						edgeID = cubeedges[cubeCase][edgeCounter+1];
						// if the edge isnt cached yet, cache it
						if(!cacheAvailable[edgeID])
						{
							// add the edge and get its index
							v1xuint = i + verts[edges[edgeID][0]][2];
							v1yuint = j + verts[edges[edgeID][0]][1];
							v1zuint = k + verts[edges[edgeID][0]][0];
							v2xuint = i + verts[edges[edgeID][1]][2];
							v2yuint = j + verts[edges[edgeID][1]][1];
							v2zuint = k + verts[edges[edgeID][1]][0];
							vertOffset1 = v1xuint + v1yuint*m_Width + v1zuint*m_Width*m_Height;
							vertOffset2 = v2xuint + v2yuint*m_Width + v2zuint*m_Width*m_Height;
							v1x = (GLfloat)(v1xuint);
							v1y = (GLfloat)(v1yuint);
							v1z = (GLfloat)(v1zuint);
							v2x = (GLfloat)(v2xuint);
							v2y = (GLfloat)(v2yuint);
							v2z = (GLfloat)(v2zuint);
							c1x = (GLfloat)((float)m_Red[vertOffset1]/255.0f);
							c1y = (GLfloat)((float)m_Green[vertOffset1]/255.0f);
							c1z = (GLfloat)((float)m_Blue[vertOffset1]/255.0f);
							c2x = (GLfloat)((float)m_Red[vertOffset2]/255.0f);
							c2y = (GLfloat)((float)m_Green[vertOffset2]/255.0f);
							c2z = (GLfloat)((float)m_Blue[vertOffset2]/255.0f);
							den1 = m_Data[voxelOffset+offsetTable[edges[edgeID][0]]];
							den2 = m_Data[voxelOffset+offsetTable[edges[edgeID][1]]];
							getNormal(m_Data,
									  i + verts[edges[edgeID][0]][2],
									  j + verts[edges[edgeID][0]][1],
									  k + verts[edges[edgeID][0]][0],
									  n1x, n1y, n1z);
							getNormal(m_Data,
									  i + verts[edges[edgeID][1]][2],
									  j + verts[edges[edgeID][1]][1],
									  k + verts[edges[edgeID][1]][0],
									  n2x, n2y, n2z);
							contourGeometry->addEdge(newEdge, v1x, v1y, v1z, n1x, n1y, n1z, c1x, c1y, c1z,
													 v2x, v2y, v2z, n2x, n2y, n2z, c2x, c2y, c2z, den1, den2);
							// The location in the cache is determined by two table lookups:
							// FromEdgeToChacheLookup[edgeID] finds which table the edge is cached in
							// offsetInSlice+cacheOffsetTable[edgeID] determines the location in that table
							m_Buffers->m_EdgeCaches[FromEdgeToChacheLookup[edgeID]][offsetInSlice+cacheOffsetTable[edgeID]] = newEdge;
						}
					}
					// All appropriate edges are now cached
					// Build the triangles, using indexes from the cache
					// for each triangle
					for(triangleCounter=0; triangleCounter<cubes[cubeCase][0]; triangleCounter++)
					{
						e1 = cubes[cubeCase][triangleCounter*3+1];
						e2 = cubes[cubeCase][triangleCounter*3+2];
						e3 = cubes[cubeCase][triangleCounter*3+3];
						contourGeometry->addTriangle(
							m_Buffers->m_EdgeCaches[FromEdgeToChacheLookup[e1]][offsetInSlice+cacheOffsetTable[e1]],
							m_Buffers->m_EdgeCaches[FromEdgeToChacheLookup[e2]][offsetInSlice+cacheOffsetTable[e2]],
							m_Buffers->m_EdgeCaches[FromEdgeToChacheLookup[e3]][offsetInSlice+cacheOffsetTable[e3]]);
					}
				}
			}
		}
		// swap the edge buffers
		m_Buffers->swapEdgeBuffers();
	}
}

static inline double maxOfThree(double n1, double n2, double n3)
{
	double max = (n1>n2?n1:n2);
	return (max>n3?max:n3);
}

void RGBAExtractor::setData(unsigned char* data, unsigned char* red, unsigned char* green, unsigned char* blue,
							unsigned int width, unsigned int height, unsigned int depth,
							double aspectX, double aspectY, double aspectZ,
							double subMinX, double subMinY, double subMinZ,
							double subMaxX, double subMaxY, double subMaxZ)
{
	m_Data = data;
	m_Red = red;
	m_Green = green;
	m_Blue = blue;
	m_Width = width;
	m_Height = height;
	m_Depth = depth;
	//m_Buffers = new MarchingCubesBuffers();
	m_Buffers->allocateEdgeBuffers(width, height);
}

void RGBAExtractor::setDefaults()
{
	// set default values for all member variables
	m_Data = 0;
	m_Width = 0;
	m_Height = 0;
	m_Depth = 0;
	m_Buffers = new MarchingCubesBuffers();
}

void RGBAExtractor::classifyVertices(unsigned int k, unsigned int* cacheMemory, GLfloat isovalue) const
{
	unsigned int i;
	unsigned int sourceOffset = k*m_Width*m_Height;
	unsigned int destOffset = 0;
	for(i=0; i<m_Width*m_Height; i++)
	{
		cacheMemory[destOffset++] = ((GLfloat)m_Data[sourceOffset++] < isovalue);
	}
}

inline unsigned int RGBAExtractor::determineCase(unsigned char* data, unsigned int* offsetTable, unsigned int index, GLfloat isovalue) const
{
	unsigned int cubeCase = 0;
	// determine the marching cube case
	cubeCase |= (m_Buffers->m_VertClassifications[0][index+offsetTable[0]])<<0;
	cubeCase |= (m_Buffers->m_VertClassifications[0][index+offsetTable[1]])<<1;
	cubeCase |= (m_Buffers->m_VertClassifications[1][index+offsetTable[2]])<<2;
	cubeCase |= (m_Buffers->m_VertClassifications[1][index+offsetTable[3]])<<3;
	cubeCase |= (m_Buffers->m_VertClassifications[0][index+offsetTable[4]])<<4;
	cubeCase |= (m_Buffers->m_VertClassifications[0][index+offsetTable[5]])<<5;
	cubeCase |= (m_Buffers->m_VertClassifications[1][index+offsetTable[6]])<<6;
	cubeCase |= (m_Buffers->m_VertClassifications[1][index+offsetTable[7]])<<7;
	return cubeCase;
}

void RGBAExtractor::calcOffsets(unsigned int* offsets, unsigned int width, unsigned int height)
{
	// calculate the table of offsets given the width and height
	offsets[0] = verts[0][0]*width*height + verts[0][1]*width + verts[0][2];
	offsets[1] = verts[1][0]*width*height + verts[1][1]*width + verts[1][2];
	offsets[2] = verts[2][0]*width*height + verts[2][1]*width + verts[2][2];
	offsets[3] = verts[3][0]*width*height + verts[3][1]*width + verts[3][2];
	offsets[4] = verts[4][0]*width*height + verts[4][1]*width + verts[4][2];
	offsets[5] = verts[5][0]*width*height + verts[5][1]*width + verts[5][2];
	offsets[6] = verts[6][0]*width*height + verts[6][1]*width + verts[6][2];
	offsets[7] = verts[7][0]*width*height + verts[7][1]*width + verts[7][2];
	offsets[8] = verts[8][0]*width*height + verts[8][1]*width + verts[8][2];
}
