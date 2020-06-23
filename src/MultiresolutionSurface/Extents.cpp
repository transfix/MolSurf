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
#include <MultiresolutionSurface/Extents.h>

Extents::Extents()
{
	setExtents(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5);
}

Extents::Extents(
	double xMin, double xMax,
	double yMin, double yMax,
	double zMin, double zMax
)
{
	setExtents(
		xMin, xMax,
		yMin, yMax,
		zMin, zMax);
}

Extents::~Extents()
{
}

void Extents::setExtents(
	double xMin, double xMax,
	double yMin, double yMax,
	double zMin, double zMax
)
{
	m_XMin = xMin;
	m_XMax = xMax;
	m_YMin = yMin;
	m_YMax = yMax;
	m_ZMin = zMin;
	m_ZMax = zMax;
}

Vector Extents::getOrigin() const
{
	return Vector(
			   (float)((m_XMin + m_XMax)/2.0),
			   (float)((m_YMin + m_YMax)/2.0),
			   (float)((m_ZMin + m_ZMax)/2.0),
			   1.0f
		   );
}

void Extents::setOrigin(Vector vector, const Extents& boundaryExtents)
{
	double width = (m_XMax-m_XMin)/2.0;
	double height = (m_YMax-m_YMin)/2.0;
	double depth = (m_ZMax-m_ZMin)/2.0;
	vector[3] = 1.0;
	if(vector[0]+width > boundaryExtents.m_XMax)
	{
		vector[0] = (float)(boundaryExtents.m_XMax-width);
	}
	else if(vector[0]-width < boundaryExtents.m_XMin)
	{
		vector[0] = (float)(boundaryExtents.m_XMin+width);
	}
	if(vector[1]+height > boundaryExtents.m_YMax)
	{
		vector[1] = (float)(boundaryExtents.m_YMax-height);
	}
	else if(vector[1]-height < boundaryExtents.m_YMin)
	{
		vector[1] = (float)(boundaryExtents.m_YMin+height);
	}
	if(vector[2]+depth > boundaryExtents.m_ZMax)
	{
		vector[2] = (float)(boundaryExtents.m_ZMax-depth);
	}
	else if(vector[2]-depth < boundaryExtents.m_ZMin)
	{
		vector[2] = (float)(boundaryExtents.m_ZMin+depth);
	}
	m_XMin = vector[0]-width;
	m_YMin = vector[1]-height;
	m_ZMin = vector[2]-depth;
	m_XMax = vector[0]+width;
	m_YMax = vector[1]+height;
	m_ZMax = vector[2]+depth;
}

void Extents::move(const Vector& vector)
{
	m_XMin += vector[0];
	m_YMin += vector[1];
	m_ZMin += vector[2];
	m_XMax += vector[0];
	m_YMax += vector[1];
	m_ZMax += vector[2];
}

bool Extents::withinCube(const Vector& vector) const
{
	return
		vector[0]>=m_XMin && vector[0]<=m_XMax &&
		vector[1]>=m_YMin && vector[1]<=m_YMax &&
		vector[2]>=m_ZMin && vector[2]<=m_ZMax;
}

void Extents::clampTo(const Extents& boundaryExtents)
{
	m_XMax = (m_XMax<boundaryExtents.m_XMax?m_XMax:boundaryExtents.m_XMax);
	m_YMax = (m_YMax<boundaryExtents.m_YMax?m_YMax:boundaryExtents.m_YMax);
	m_ZMax = (m_ZMax<boundaryExtents.m_ZMax?m_ZMax:boundaryExtents.m_ZMax);
	m_XMin = (m_XMin>boundaryExtents.m_XMin?m_XMin:boundaryExtents.m_XMin);
	m_YMin = (m_YMin>boundaryExtents.m_YMin?m_YMin:boundaryExtents.m_YMin);
	m_ZMin = (m_ZMin>boundaryExtents.m_ZMin?m_ZMin:boundaryExtents.m_ZMin);
	// lets check to see if its a well formed box
	if(m_XMin >= m_XMax || m_YMin >= m_YMax || m_ZMin >= m_ZMax)    // bad box
	{
		m_XMin = boundaryExtents.m_XMin * 0.75 + boundaryExtents.m_XMax * 0.25;
		m_YMin = boundaryExtents.m_YMin * 0.75 + boundaryExtents.m_YMax * 0.25;
		m_ZMin = boundaryExtents.m_ZMin * 0.75 + boundaryExtents.m_ZMax * 0.25;
		m_XMax = boundaryExtents.m_XMin * 0.25 + boundaryExtents.m_XMax * 0.75;
		m_YMax = boundaryExtents.m_YMin * 0.25 + boundaryExtents.m_YMax * 0.75;
		m_ZMax = boundaryExtents.m_ZMin * 0.25 + boundaryExtents.m_ZMax * 0.75;
	}
}

double Extents::getXMin() const
{
	return m_XMin;
}

double Extents::getYMin() const
{
	return m_YMin;
}

double Extents::getZMin() const
{
	return m_ZMin;
}

double Extents::getXMax() const
{
	return m_XMax;
}

double Extents::getYMax() const
{
	return m_YMax;
}

double Extents::getZMax() const
{
	return m_ZMax;
}

void Extents::setXMin(double xMin)
{
	m_XMin = xMin;
}

void Extents::setYMin(double yMin)
{
	m_YMin = yMin;
}

void Extents::setZMin(double zMin)
{
	m_ZMin = zMin;
}

void Extents::setXMax(double xMax)
{
	m_XMax = xMax;
}

void Extents::setYMax(double yMax)
{
	m_YMax = yMax;
}

void Extents::setZMax(double zMax)
{
	m_ZMax = zMax;
}
