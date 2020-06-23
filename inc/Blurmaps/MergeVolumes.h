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
#ifndef CCV_TEXMOL_MERGE_VOLUMES_H
#define CCV_TEXMOL_MERGE_VOLUMES_H

#include <Utility/utility.h>
#include <UsefulMath/Matrix.h>
#include <UsefulMath/Vector.h>
#include <UsefulMath/Quaternion.h>
#include <SimpleVolumeData/SimpleVolumeData.h>

class SimpleVolumeData;

class MergeVolumes
{
	public:
		MergeVolumes();
		virtual ~MergeVolumes();

		void setTransformation(float c1x, float c1y, float c1z,
							   float n1x, float n1y, float n1z,
							   float c2x, float c2y, float c2z,
							   float n2x, float n2y, float n2z);
		bool mergeVolumes(SimpleVolumeData* input1, SimpleVolumeData* input2);

	protected:
		void transformPoint(float* point);
		void mergePoint(float* point2, float value2, float* data1, float* origin1, float* span1, unsigned int* dim1);
		float getCombinedSDF(float value1, float value2);
		float getInterpValue(float* data, float* orig, float* span, unsigned int* dim, float* point);
		float getValueAt(int i, int j, int k, float* data, unsigned int* dim, bool* ok);
		CCVOpenGLMath::Matrix m_Transformation;
		bool m_Transform;
};

#endif
