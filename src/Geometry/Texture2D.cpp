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
#include <Geometry/Texture2D.h>

Texture2D::Texture2D(int width, int height, unsigned char *data,
                     int numBytesPerEntry) {
  m_Width = width;
  m_Height = height;
  m_NumBytesPerEntry = numBytesPerEntry;
  m_Data = new unsigned char[width * height * m_NumBytesPerEntry];
  int i;
  for (i = 0; i < width * height * m_NumBytesPerEntry; i++) {
    m_Data[i] = data[i];
  }
  m_Id = -1;
}

Texture2D::~Texture2D() {
  delete[] m_Data;
  m_Data = 0;
  if (m_Id != -1) {
    glDeleteTextures(1, &m_Id);
  }
  m_Id = -1;
}

Texture2D::Texture2D(const Texture2D &copy) { set(copy); }

Texture2D &Texture2D::operator=(const Texture2D &copy) { return set(copy); }

Texture2D &Texture2D::set(const Texture2D &copy) {
  if (this != &copy) {
    if (m_Data) {
      delete[] m_Data;
      m_Data = 0;
    }
    m_Width = copy.m_Width;
    m_Height = copy.m_Height;
    m_NumBytesPerEntry = copy.m_NumBytesPerEntry;
    if (copy.m_Data) {
      m_Data = new unsigned char[m_Width * m_Height * m_NumBytesPerEntry];
      int i;
      for (i = 0; i < m_Width * m_Height * m_NumBytesPerEntry; i++) {
        m_Data[i] = copy.m_Data[i];
      }
    }
  }
  return *this;
}

bool Texture2D::enable(bool enable) {
  if (enable) {
    if (!m_Data) {
      return false;
    }
    if (m_NumBytesPerEntry < 1) {
      return false;
    }
    if (m_Id == -1) {
      glGenTextures(1, &m_Id);
      glBindTexture(GL_TEXTURE_2D, m_Id);
      glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      switch (m_NumBytesPerEntry) {
      case 1:
        glTexImage2D(GL_TEXTURE_2D, 0, m_NumBytesPerEntry, m_Width, m_Height, 0,
                     GL_LUMINANCE, GL_UNSIGNED_BYTE, m_Data);
        break;
      case 2:
        glTexImage2D(GL_TEXTURE_2D, 0, m_NumBytesPerEntry, m_Width, m_Height, 0,
                     GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, m_Data);
        break;
      case 3:
        glTexImage2D(GL_TEXTURE_2D, 0, m_NumBytesPerEntry, m_Width, m_Height, 0,
                     GL_RGB, GL_UNSIGNED_BYTE, m_Data);
        break;
      case 4:
        glTexImage2D(GL_TEXTURE_2D, 0, m_NumBytesPerEntry, m_Width, m_Height, 0,
                     GL_RGBA, GL_UNSIGNED_BYTE, m_Data);
        break;
      default: // should be safe to just use 1
        glTexImage2D(GL_TEXTURE_2D, 0, m_NumBytesPerEntry, m_Width, m_Height, 0,
                     GL_LUMINANCE, GL_UNSIGNED_BYTE, m_Data);
        break;
      }
    }
  }
  return true;
}
