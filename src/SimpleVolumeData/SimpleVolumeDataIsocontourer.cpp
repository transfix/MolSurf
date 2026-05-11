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
#include <SimpleVolumeData/SimpleVolumeDataIsocontourer.h>

SimpleVolumeDataIsocontourer::SimpleVolumeDataIsocontourer() {}

SimpleVolumeDataIsocontourer::~SimpleVolumeDataIsocontourer() {}

Geometry *SimpleVolumeDataIsocontourer::getIsocontour(SimpleVolumeData *sData,
                                                      float isovalue) {
  // extract an isocontour
  Geometry *geometry = 0;
  {
    CCVContouring::MultiContour *contourExtractor =
        new CCVContouring::MultiContour();
    int width = sData->getWidth();
    int height = sData->getHeight();
    int depth = sData->getDepth();
    // convert sData to unsigned char for isocontouring
    if (sData->getNumberOfVariables() < 4) {
      // RawIV data set
      unsigned char *data1 = new unsigned char[width * height * depth];
      SimpleVolumeData::getNormalizedUnsignedCharData(
          sData->getData(0), sData->getType(0), data1, width, height, depth);
      contourExtractor->setData(
          data1, sData->getWidth(), sData->getHeight(), sData->getDepth(),
          sData->getMaxX() - sData->getMinX(),
          sData->getMaxY() - sData->getMinY(),
          sData->getMaxZ() - sData->getMinZ(), 0, 0, 0, 1, 1, 1,
          sData->getMinX(), sData->getMinY(), sData->getMinZ(),
          sData->getMaxX(), sData->getMaxY(), sData->getMaxZ());
      double minVal = 0.0;
      double maxVal = 1.0;
      sData->getMinMax(0, &minVal, &maxVal);
      contourExtractor->addContour(
          0, (float)((isovalue - minVal) / (maxVal - minVal) * 255), 0.8f, 0.8f,
          0.8f);
      geometry = contourExtractor->getGeometry();
      delete[] data1;
    } else {
      // RawV data set
      unsigned char *data1 = new unsigned char[width * height * depth];
      unsigned char *data2 = new unsigned char[width * height * depth];
      unsigned char *data3 = new unsigned char[width * height * depth];
      unsigned char *data4 = new unsigned char[width * height * depth];
      SimpleVolumeData::getUnsignedCharData(
          sData->getData(0), sData->getType(0), data1, width, height, depth);
      SimpleVolumeData::getUnsignedCharData(
          sData->getData(1), sData->getType(1), data2, width, height, depth);
      SimpleVolumeData::getUnsignedCharData(
          sData->getData(2), sData->getType(2), data3, width, height, depth);
      SimpleVolumeData::getUnsignedCharData(
          sData->getData(3), sData->getType(3), data4, width, height, depth);
      contourExtractor->setData(
          data4, data1, data2, data3, sData->getWidth(), sData->getHeight(),
          sData->getDepth(), sData->getMaxX() - sData->getMinX(),
          sData->getMaxY() - sData->getMinY(),
          sData->getMaxZ() - sData->getMinZ(), 0, 0, 0, 1, 1, 1,
          sData->getMinX(), sData->getMinY(), sData->getMinZ(),
          sData->getMaxX(), sData->getMaxY(), sData->getMaxZ());
      double minVal = 0.0;
      double maxVal = 1.0;
      sData->getMinMax(0, &minVal, &maxVal);
      contourExtractor->addContour(
          0, (float)((isovalue - minVal) / (maxVal - minVal) * 255), 0.8f, 0.8f,
          0.8f);
      geometry = contourExtractor->getGeometry();
      delete[] data1;
      delete[] data2;
      delete[] data3;
      delete[] data4;
    }
    delete contourExtractor;
  }
  if (geometry == 0) {
    // this was returning an object set to 0, which is not really standard, and
    // annoying to check for hopefuly I haven't broken anyone's code
    exit(-1);
  }
  return geometry;
}

bool SimpleVolumeDataIsocontourer::colorIsocontour(Geometry *geometry,
                                                   SimpleVolumeData *sData) {
  if (!geometry || !sData) {
    return false;
  }
  // get extent of the values for RGB
  double minVals[3];
  double maxVals[3];
  minVals[0] = minVals[1] = minVals[2] = 0;
  maxVals[0] = maxVals[1] = maxVals[2] = 1;
  if (sData->getNumberOfVariables() < 4) {
    sData->getMinMax(0, &(minVals[0]), &(maxVals[0]));
  } else {
    sData->getMinMax(0, &(minVals[0]), &(maxVals[0]));
    sData->getMinMax(1, &(minVals[1]), &(maxVals[1]));
    sData->getMinMax(2, &(minVals[2]), &(maxVals[2]));
  }
  unsigned int i;
  for (i = 0; i < geometry->m_NumTriVerts; i++) {
    double x = geometry->m_TriVerts[i * 3 + 0];
    double y = geometry->m_TriVerts[i * 3 + 1];
    double z = geometry->m_TriVerts[i * 3 + 2];
    geometry->m_TriVertColorsTransparent[i * 4 + 0] = 0;
    geometry->m_TriVertColorsTransparent[i * 4 + 1] = 0;
    geometry->m_TriVertColorsTransparent[i * 4 + 2] = 0;
    if (!sData->inVolume(x, y, z)) {
      continue;
    }
    if (sData->getNumberOfVariables() < 4) {
      double funcVal = sData->getValueAt(0, x, y, z);
      // there are 3 cases for the ranges!
      // case 1, no negatives
      if (minVals[0] >= 0) {
        double color = 0;
        double val = (funcVal - minVals[0]) / (maxVals[0] - minVals[0]);
        color = 1.0 - (exp(1) - exp(val)) / exp(1);
        geometry->m_TriVertColorsTransparent[i * 4 + 0] = 1.0;
        geometry->m_TriVertColorsTransparent[i * 4 + 1] = (float)color;
        geometry->m_TriVertColorsTransparent[i * 4 + 2] = (float)color;
      }
      // case 2, no positives
      else if (maxVals[0] <= 0) {
        double color = 0;
        double val = (funcVal - minVals[0]) / (maxVals[0] - minVals[0]);
        color = 1.0 - (exp(1) - exp(val)) / exp(1);
        geometry->m_TriVertColorsTransparent[i * 4 + 0] = 1.0;
        geometry->m_TriVertColorsTransparent[i * 4 + 1] = (float)color;
        geometry->m_TriVertColorsTransparent[i * 4 + 2] = (float)color;
      }
      // case 3: both positives and negatives
      else {
        if (funcVal >= 0) {
          double color = 0;
          double val = (funcVal - 0) / (maxVals[0]);
          color = (exp(1) - exp(val)) / exp(1);
          color = 1.0 - color;
          geometry->m_TriVertColorsTransparent[i * 4 + 0] = 1.0;
          geometry->m_TriVertColorsTransparent[i * 4 + 1] = (float)color;
          geometry->m_TriVertColorsTransparent[i * 4 + 2] = (float)color;
        } else {
          double color = 0;
          double val = (funcVal) / (minVals[0]);
          color = (exp(1) - exp(val)) / exp(1);
          color = 1.0 - color;
          geometry->m_TriVertColorsTransparent[i * 4 + 0] = (float)color;
          geometry->m_TriVertColorsTransparent[i * 4 + 1] = 1.0;
          geometry->m_TriVertColorsTransparent[i * 4 + 2] = (float)color;
        }
      }
    } else {
      double rFuncVal = sData->getValueAt(0, x, y, z);
      double gFuncVal = sData->getValueAt(1, x, y, z);
      double bFuncVal = sData->getValueAt(2, x, y, z);
      geometry->m_TriVertColorsTransparent[i * 4 + 0] =
          (float)((rFuncVal - minVals[0]) / (maxVals[0] - minVals[0]));
      geometry->m_TriVertColorsTransparent[i * 4 + 1] =
          (float)((gFuncVal - minVals[1]) / (maxVals[1] - minVals[1]));
      geometry->m_TriVertColorsTransparent[i * 4 + 2] =
          (float)((bFuncVal - minVals[2]) / (maxVals[2] - minVals[2]));
    }
  }
  return true;
}

bool SimpleVolumeDataIsocontourer::getFunctionValues(SimpleVolumeData *sData,
                                                     Geometry *geometry,
                                                     double *funcVals) {
  if (!sData || !geometry || !funcVals) {
    return false;
  }
  unsigned int i;
  for (i = 0; i < geometry->m_NumTriVerts; i++) {
    double x = geometry->m_TriVerts[i * 3 + 0];
    double y = geometry->m_TriVerts[i * 3 + 1];
    double z = geometry->m_TriVerts[i * 3 + 2];
    if (!sData->inVolume(x, y, z)) {
      continue;
    }
    if (sData->getNumberOfVariables() < 4) {
      funcVals[i] = sData->getValueAt(0, x, y, z);
    } else {
      funcVals[i] = sData->getValueAt(3, x, y, z);
    }
  }
  return true;
}
