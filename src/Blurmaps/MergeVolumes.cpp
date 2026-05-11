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
#include <Blurmaps/MergeVolumes.h>

MergeVolumes::MergeVolumes() { m_Transform = false; }

MergeVolumes::~MergeVolumes() {}

void MergeVolumes::setTransformation(float c1x, float c1y, float c1z, float n1x,
                                     float n1y, float n1z, float c2x, float c2y,
                                     float c2z, float n2x, float n2y,
                                     float n2z) {
  CCVOpenGLMath::Vector old_center(c1x, c1y, c1z, 1);
  CCVOpenGLMath::Vector old_normal(n1x, n1y, n1z, 0);
  CCVOpenGLMath::Vector new_center(c2x, c2y, c2z, 1);
  CCVOpenGLMath::Vector new_normal(n2x, n2y, n2z, 0);
  // shift center to new one.
  // rotate to align quaternions
  // for each mol, we need a predefined center and normal.
  // translate mol_center to input_center
  // now to rotate:
  //     cross product gives axis of rotation
  //     dot product gives the angle
  old_normal.normalize();
  new_normal.normalize();
  CCVOpenGLMath::Vector axis_of_rotation = old_normal.cross(new_normal);
  axis_of_rotation.normalize();
  double angle_of_rotation = acos(new_normal.dot(old_normal));
  float w = cos(angle_of_rotation / 2.0);
  float x = axis_of_rotation[0] * sin(angle_of_rotation / 2.0);
  float y = axis_of_rotation[1] * sin(angle_of_rotation / 2.0);
  float z = axis_of_rotation[2] * sin(angle_of_rotation / 2.0);
  CCVOpenGLMath::Quaternion quaternion(w, x, y, z);
  CCVOpenGLMath::Matrix rotation = quaternion.buildMatrix();
  CCVOpenGLMath::Matrix translation =
      CCVOpenGLMath::Matrix::translation(new_center - old_center);
  CCVOpenGLMath::Matrix trToOrigin =
      CCVOpenGLMath::Matrix::translation(-c1x, -c1y, -c1z);
  CCVOpenGLMath::Matrix trFromOrigin =
      CCVOpenGLMath::Matrix::translation(c1x, c1y, c1z);
  m_Transformation = trToOrigin;
  m_Transformation = m_Transformation.preMultiplication(rotation);
  m_Transformation = m_Transformation.preMultiplication(trFromOrigin);
  m_Transformation = m_Transformation.preMultiplication(translation);
  // m_Transformation = translation;
  m_Transform = true;
}

void MergeVolumes::transformPoint(float *point) {
  CCVOpenGLMath::Vector p(point[0], point[1], point[2], 1);
  CCVOpenGLMath::Vector ret = m_Transformation * p;
  point[0] = ret[0] / ret[3];
  point[1] = ret[1] / ret[3];
  point[2] = ret[2] / ret[3];
}

float MergeVolumes::getCombinedSDF(float value1, float value2) {
  float newvalue = 0;
  if (value1 >= 0 && value2 >= 0) {
    newvalue = (value1 < value2) ? value1 : value2;
  }
  if (value1 <= 0 && value2 <= 0) {
    newvalue = (value1 > value2) ? value1 : value2;
  }
  if (value1 >= 0 && value2 <= 0) {
    newvalue = value1;
  }
  if (value1 <= 0 && value2 >= 0) {
    newvalue = value2;
  }
  return newvalue;
}

float MergeVolumes::getValueAt(int i, int j, int k, float *data,
                               unsigned int *dim, bool *ok) {
  if (i < 0 || i >= dim[0]) {
    *ok = false;
    return 0;
  }
  if (j < 0 || j >= dim[1]) {
    *ok = false;
    return 0;
  }
  if (k < 0 || k >= dim[2]) {
    *ok = false;
    return 0;
  }
  int index = i * dim[1] * dim[2] + j * dim[2] + k;
  *ok = true;
  return data[index];
}

// wont be correct at boundaries unless boundaries are 0
float MergeVolumes::getInterpValue(float *data, float *orig, float *span,
                                   unsigned int *dim, float *point) {
  float interpValue = 0;
  int numberOfInfluences = 0;
  float v = 0;
  bool ok = false;
  int closestIndex[3];
  closestIndex[0] = (point[0] - orig[0]) / span[0];
  closestIndex[1] = (point[1] - orig[1]) / span[1];
  closestIndex[2] = (point[2] - orig[2]) / span[2];
  v = getValueAt(closestIndex[0], closestIndex[1], closestIndex[2], data, dim,
                 &ok);
  if (ok) {
    interpValue += v;
    numberOfInfluences++;
  }
  v = getValueAt(closestIndex[0] + 1, closestIndex[1], closestIndex[2], data,
                 dim, &ok);
  if (ok) {
    interpValue += v;
    numberOfInfluences++;
  }
  v = getValueAt(closestIndex[0], closestIndex[1] + 1, closestIndex[2], data,
                 dim, &ok);
  if (ok) {
    interpValue += v;
    numberOfInfluences++;
  }
  v = getValueAt(closestIndex[0] + 1, closestIndex[1] + 1, closestIndex[2],
                 data, dim, &ok);
  if (ok) {
    interpValue += v;
    numberOfInfluences++;
  }
  v = getValueAt(closestIndex[0], closestIndex[1], closestIndex[2] + 1, data,
                 dim, &ok);
  if (ok) {
    interpValue += v;
    numberOfInfluences++;
  }
  v = getValueAt(closestIndex[0] + 1, closestIndex[1], closestIndex[2] + 1,
                 data, dim, &ok);
  if (ok) {
    interpValue += v;
    numberOfInfluences++;
  }
  v = getValueAt(closestIndex[0], closestIndex[1] + 1, closestIndex[2] + 1,
                 data, dim, &ok);
  if (ok) {
    interpValue += v;
    numberOfInfluences++;
  }
  v = getValueAt(closestIndex[0] + 1, closestIndex[1] + 1, closestIndex[2] + 1,
                 data, dim, &ok);
  if (ok) {
    interpValue += v;
    numberOfInfluences++;
  }
  if (numberOfInfluences == 0) {
    return 0;
  }
  return interpValue / ((float)numberOfInfluences);
}

void MergeVolumes::mergePoint(float *point2, float value2, float *data1,
                              float *origin1, float *span1,
                              unsigned int *dim1) {
  int closestIndex[3];
  closestIndex[0] = (point2[0] - origin1[0]) / span1[0];
  closestIndex[1] = (point2[1] - origin1[1]) / span1[1];
  closestIndex[2] = (point2[2] - origin1[2]) / span1[2];
  if (closestIndex[0] < 0) {
    return;
  }
  if (closestIndex[1] < 0) {
    return;
  }
  if (closestIndex[2] < 0) {
    return;
  }
  if (closestIndex[0] >= dim1[0]) {
    return;
  }
  if (closestIndex[1] >= dim1[1]) {
    return;
  }
  if (closestIndex[2] >= dim1[2]) {
    return;
  }
  int index = closestIndex[0] * dim1[0] * dim1[1] + closestIndex[1] * dim1[1] +
              closestIndex[2];
  float value1 = data1[index];
  data1[index] = getCombinedSDF(value1, value2);
}

/// merges input2 into input1, assuming both are SDF volumes
bool MergeVolumes::mergeVolumes(SimpleVolumeData *input1,
                                SimpleVolumeData *input2) {
  if (!input1 || !input2) {
    return false;
  }
  unsigned int dim1[3];
  unsigned int dim2[3];
  float span1[3];
  float span2[3];
  float origin1[3];
  float origin2[3];
  float *data1;
  float *data2;
  dim1[0] = input1->getWidth();
  dim1[1] = input1->getHeight();
  dim1[2] = input1->getDepth();
  dim2[0] = input2->getWidth();
  dim2[1] = input2->getHeight();
  dim2[2] = input2->getDepth();
  span1[0] = input1->getSpanX();
  span1[1] = input1->getSpanY();
  span1[2] = input1->getSpanZ();
  span2[0] = input2->getSpanX();
  span2[1] = input2->getSpanY();
  span2[2] = input2->getSpanZ();
  origin1[0] = input1->getMinX();
  origin1[1] = input1->getMinY();
  origin1[2] = input1->getMinZ();
  origin2[0] = input2->getMinX();
  origin2[1] = input2->getMinY();
  origin2[2] = input2->getMinZ();
  data1 = (float *)(input1->getData(0));
  data2 = (float *)(input2->getData(0));
  int i, j, k, c = 0;
  for (i = 0; i < dim2[0]; i++) {
    for (j = 0; j < dim2[1]; j++) {
      for (k = 0; k < dim2[2]; k++) {
        float point2[3];
        point2[0] = origin2[0] + i * span2[0];
        point2[1] = origin2[1] + j * span2[1];
        point2[2] = origin2[2] + k * span2[2];
        float value2 = data2[c];
        c++;
        if (m_Transform) {
          transformPoint(point2);
        }
        mergePoint(point2, value2, data1, origin1, span1, dim1);
      }
    }
  }
  return true;
}
