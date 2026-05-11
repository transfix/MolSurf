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
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <Utility/utility.h>
#include <GeometryFileTypes/GeometryFileType.h>
#include <GeometryFileTypes/GeometryLoader.h>
#include <GeometryFileTypes/MayaOBJFile.h>
#include <GeometryFileTypes/RawcFile.h>
#include <GeometryFileTypes/RawFile.h>
#include <GeometryFileTypes/RawncFile.h>
#include <GeometryFileTypes/RawnFile.h>

GeometryLoader::GeometryLoader() {
  // add each geometry type to the maps
  addGeometryFileType(MayaOBJFile::getRepresentative());
  addGeometryFileType(RawncFile::getRepresentative());
  addGeometryFileType(RawcFile::getRepresentative());
  addGeometryFileType(RawnFile::getRepresentative());
  addGeometryFileType(RawFile::getRepresentative());
}

GeometryLoader::~GeometryLoader() {}

bool GeometryLoader::endsWith(string str, string substr) {
  if (str.length() < substr.length()) {
    return false;
  }
  if (str.substr(str.length() - substr.length(), substr.length()) == substr) {
    return true;
  }
  return false;
}

bool GeometryLoader::saveFile(const string &fileName, Geometry *geometry) {
  if (endsWith(fileName, "rawnc")) {
    return saveFile(fileName, "Rawnc files (*.rawnc)", geometry);
  }
  if (endsWith(fileName, "rawn")) {
    return saveFile(fileName, "Rawn files (*.rawn)", geometry);
  }
  if (endsWith(fileName, "rawc")) {
    return saveFile(fileName, "Rawc files (*.rawc)", geometry);
  }
  if (endsWith(fileName, "raw")) {
    return saveFile(fileName, "Raw files (*.raw)", geometry);
  }
  if (endsWith(fileName, "c2c")) {
    return saveFile(fileName, "c2c files (*.c2c)", geometry);
  }
  if (endsWith(fileName, "obj")) {
    return saveFile(fileName, "Maya OBJ files (*.obj)", geometry);
  }
  return saveFile(fileName, "Rawnc files (*.rawnc)", geometry);
}

bool GeometryLoader::saveFile(const string &fileName,
                              const string &selectedFilter,
                              Geometry *geometry) {
  if (m_FilterMap.count(selectedFilter) != 0) {
    string extension;
    string longName;
    if (fileName.rfind('.') != string::npos) {
      extension = fileName.substr(fileName.rfind('.') + 1);
    }
    GeometryFileType *type = m_FilterMap[selectedFilter];
    // if no extension, add one
    if (extension.empty()) {
      longName = fileName + "." + type->extension();
    } else if (!(extension == type->extension())) {
      longName =
          fileName.substr(0, fileName.rfind('.') + 1) + type->extension();
    }
    { longName = fileName; }
    return type->saveFile(geometry, longName);
  } else {
    return false;
  }
}

Geometry *GeometryLoader::loadFile(const string &fileName) {
  string extension;
  if (fileName.rfind('.') != string::npos) {
    extension = fileName.substr(fileName.rfind('.') + 1);
  }
  if (extension.empty() || m_ExtensionMap.count(extension) == 0) {
    // test every file type to find the correct one
    return tryAll(fileName);
  } else {
    // try to load the file with the correct file type
    GeometryFileType *type = m_ExtensionMap[extension];
    Geometry *geometry = type->loadFile(fileName);

    if (!geometry) // failed, try the other loaders
    {
      return tryAll(fileName);
    } else {
      // success
      return geometry;
    }
  }
}

string GeometryLoader::getLoadFilterString() {
  string str("All Geometry Files ");
  str.append(getAllExtensions());
  // iterate through each loader and combine all the filters
  std::map<string, GeometryFileType *>::iterator it;
  for (it = m_FilterMap.begin(); it != m_FilterMap.end(); ++it) {
    str.append(";;" + (it->first));
  }
  return str;
}

string GeometryLoader::getSaveFilterString() {
  string str;
  bool first = true;
  // iterate through each loader and combine all the filters
  std::map<string, GeometryFileType *>::iterator it;
  for (it = m_FilterMap.begin(); it != m_FilterMap.end(); ++it) {
    if (first) {
      first = false;
      str = it->first;
    } else {
      str.append(";;" + (it->first));
    }
  }
  return str;
}

string GeometryLoader::getAllExtensions() {
  string str("(");
  bool first = true;
  // iterate through each loader and combine all the filters
  std::map<string, GeometryFileType *>::iterator it;
  for (it = m_FilterMap.begin(); it != m_FilterMap.end(); ++it) {
    if (first) {
      first = false;
      str.append("*." + (it->first));
    } else {
      str.append(" *." + (it->first));
    }
  }
  str.append(")");
  return str;
}

Geometry *GeometryLoader::tryAll(const string &fileName) {
  // iterate through each loader and call checkType to determine
  // which loader can load the file
  std::map<string, GeometryFileType *>::iterator it;
  for (it = m_ExtensionMap.begin(); it != m_ExtensionMap.end(); ++it) {
    if ((it->second)->checkType(fileName)) // found it
    {
      return (it->second)->loadFile(fileName);
    }
  }
  // didnt find the right loader
  return 0;
}

void GeometryLoader::addGeometryFileType(GeometryFileType *type) {
  m_ExtensionMap[type->extension()] = type;
  m_FilterMap[type->filter()] = type;
}

bool GeometryLoader::isValidExtension(string extension) {
  std::map<string, GeometryFileType *>::iterator it;
  for (it = m_ExtensionMap.begin(); it != m_ExtensionMap.end(); ++it) {
    if (strcmp((it->first).c_str(), extension.c_str()) == 0) {
      return true;
    }
  }
  return false;
}
