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
#include <SignDistanceFunction_v2/RawivParser.h>

#include <SimpleVolumeData/SimpleVolumeData.h>
#include <VolumeFileTypes/VolumeLoader.h>

using namespace CVCUtility;

RawivParser::RawivParser() {}

RawivParser::~RawivParser() {}

bool RawivParser::parse(Reg3Data<float> *data, const char *fname) {
  int nverts, ncells;
  float minext[3], maxext[3];
  if (!isRawivFile(fname)) {
    error(string(fname) + " is not a rawiv file!");
  }
  // determine file data type
  struct stat filestat;
  if (stat(fname, &filestat) < 0) {
    error("cannot find data file " + string(fname));
  }
  int sz = filestat.st_size;
  DiskIO *pio = new BufferedIO(fname);
  if (!pio->open()) {
    error("Data File Open Failed");
  }
  pio->get(minext, 3);
  pio->get(maxext, 3);
  pio->get(&nverts, 1);
  pio->get(&ncells, 1);
  pio->get(data->m_dim, 3);
  pio->get(data->m_orig, 3);
  pio->get(data->m_span, 3);
  nverts = (data->m_dim[0]) * (data->m_dim[1]) * (data->m_dim[2]);
  if (data->p_data) {
    delete[] data->p_data;
  }
  data->p_data = new float[nverts];
  pio->get(data->p_data, nverts);
  data->init();
  pio->close();
  delete pio;
  return true;
}

/*
bool RawivParser::write(const Reg3Data<float>& data, const char* fname)
{
        DiskIO* pio = new BufferedIO(fname, DiskIO::WRITE);
        if(!pio->open())
        {
                error("Cannot Open Data File to write");
        }
        float maxext[3];
        maxext[0] = data.m_orig[0] + (data.m_dim[0]-1)*data.m_span[0];
        maxext[1] = data.m_orig[1] + (data.m_dim[1]-1)*data.m_span[1];
        maxext[2] = data.m_orig[2] + (data.m_dim[2]-1)*data.m_span[2];
        int nverts = data.getNVerts();
        int ncells = data.getNCells();
        pio->put(data.m_orig, 3);
        pio->put(maxext, 3);
        pio->put(&nverts, 1);
        pio->put(&ncells, 1);
        pio->put(data.m_dim, 3);
        pio->put(data.m_orig, 3);
        pio->put(data.m_span, 3);
        pio->put(data.p_data, nverts);
        pio->close(false);
        return true;
}
*/

bool RawivParser::write(const Reg3Data<float> &data, const char *fname) {
  /*	DiskIO* pio = new BufferedIO(fname, DiskIO::WRITE);
          if(!pio->open())
          {
                  error("Cannot Open Data File to write");
          } */
  float maxext[3];
  maxext[0] = data.m_orig[0] + (data.m_dim[0] - 1) * data.m_span[0];
  maxext[1] = data.m_orig[1] + (data.m_dim[1] - 1) * data.m_span[1];
  maxext[2] = data.m_orig[2] + (data.m_dim[2] - 1) * data.m_span[2];
  /*	int nverts = data.getNVerts();
          int ncells = data.getNCells();
          pio->put(data.m_orig, 3);
          pio->put(maxext, 3);
          pio->put(&nverts, 1);
          pio->put(&ncells, 1);
          pio->put(data.m_dim, 3);
          pio->put(data.m_orig, 3);
          pio->put(data.m_span, 3);
          pio->put(data.p_data, nverts);
          pio->close(false);
          return true; */
  SimpleVolumeData *sData = new SimpleVolumeData((unsigned int *)data.m_dim);
  sData->setDimensions((unsigned int *)data.m_dim);
  sData->setNumberOfVariables(1);
  sData->setData(0, data.p_data);
  sData->setType(0, SimpleVolumeData::FLOAT);
  sData->setName(0, "Signed Distance Function V_2");
  sData->setMinExtent((float *)data.m_orig);
  sData->setMaxExtent(maxext);

  VolumeLoader().saveFile(fname, sData);
}
