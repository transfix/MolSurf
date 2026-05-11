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
#include <GOAFileTypes/GOAFileType.h>
#include <GOAFileTypes/GOALoader.h>
#include <GOAFileTypes/PDB1File.h>
#include <GOAFileTypes/PDBFile.h>
#include <GOAFileTypes/PQRFile.h>
#include <GOAFileTypes/PTSFile.h>
#include <GOAFileTypes/XYZFile.h>
#include <GOAFileTypes/XYZRFile.h>
#include <PDBParser/GroupOfAtoms.h>

GOALoader::GOALoader() {
  // add each GOA type to the maps
  addGOAFileType(PTSFile::getRepresentative());
  addGOAFileType(XYZFile::getRepresentative());
  addGOAFileType(XYZRFile::getRepresentative());
  addGOAFileType(PDBFile::getRepresentative());
  addGOAFileType(PDB1File::getRepresentative());
  addGOAFileType(PQRFile::getRepresentative());
}

GOALoader::~GOALoader() {}

bool GOALoader::endsWith(string str, string substr) {
  if (str.length() < substr.length()) {
    return false;
  }
  if (str.substr(str.length() - substr.length(), substr.length()) == substr) {
    return true;
  }
  return false;
}

bool GOALoader::saveFile(const string &fileName,
                         PDBParser::GroupOfAtoms *groupOfAtoms,
                         unsigned int level,
                         CCVOpenGLMath::Matrix *transformation) {
  if (endsWith(fileName, "pts")) {
    return saveFile(fileName, "PTS files (*.pts)", groupOfAtoms, level,
                    transformation);
  }
  if (endsWith(fileName, "xyz")) {
    return saveFile(fileName, "XYZ files (*.xyz)", groupOfAtoms, level,
                    transformation);
  }
  if (endsWith(fileName, "xyzr")) {
    return saveFile(fileName, "XYZR files (*.xyzr)", groupOfAtoms, level,
                    transformation);
  }
  if (endsWith(fileName, "pdb")) {
    return saveFile(fileName, "PDB files (*.pdb)", groupOfAtoms, level,
                    transformation);
  }
  if (endsWith(fileName, "pdb1")) {
    return saveFile(fileName, "PDB1 files (*.pdb1)", groupOfAtoms, level,
                    transformation);
  }
  if (endsWith(fileName, "pqr")) {
    return saveFile(fileName, "PQR files (*.pqr)", groupOfAtoms, level,
                    transformation);
  }
  // try simple one
  return saveFile(fileName, "PTS files (*.pts)", groupOfAtoms, level,
                  transformation);
}

bool GOALoader::saveFile(const string &fileName, const string &selectedFilter,
                         PDBParser::GroupOfAtoms *groupOfAtoms,
                         unsigned int level,
                         CCVOpenGLMath::Matrix *transformation) {
  if (m_FilterMap.count(selectedFilter) != 0) {
    string extension;
    string longName;
    if (fileName.rfind('.') != string::npos) {
      extension = fileName.substr(fileName.rfind('.') + 1);
    }
    GOAFileType *type = m_FilterMap[selectedFilter];
    // if no extension, add one
    if (extension.empty()) {
      longName = fileName + "." + type->extension();
    } else {
      longName = fileName;
    }
    return type->saveFile(groupOfAtoms, longName, level, transformation);
  } else {
    return false;
  }
}

PDBParser::GroupOfAtoms *GOALoader::loadFile(const string &fileName,
                                             bool deleteWater) {
  string extension;
  if (fileName.rfind('.') != string::npos) {
    extension = fileName.substr(fileName.rfind('.') + 1);
  }
  if (extension.empty() || m_ExtensionMap.count(extension) == 0) {
    // test every file type to find the correct one
    return tryAll(fileName, deleteWater);
  } else {
    // try to load the file with the correct file type
    GOAFileType *type = m_ExtensionMap[extension];
    PDBParser::GroupOfAtoms *groupOfAtoms =
        type->loadFile(fileName, deleteWater);

    if (!groupOfAtoms) // failed, try the other loaders
    {
      return tryAll(fileName, deleteWater);
    } else {
      // success
      return groupOfAtoms;
    }
  }
}

string GOALoader::getLoadFilterString() {
  string str("All GOA Files ");
  str.append(getAllExtensions());
  // iterate through each loader and combine all the filters
  std::map<string, GOAFileType *>::iterator it;
  for (it = m_FilterMap.begin(); it != m_FilterMap.end(); ++it) {
    str.append(";;" + (it->first));
  }
  return str;
}

string GOALoader::getSaveFilterString() {
  string str;
  bool first = true;
  // iterate through each loader and combine all the filters
  std::map<string, GOAFileType *>::iterator it;
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

string GOALoader::getAllExtensions() {
  string str("(");
  bool first = true;
  // iterate through each loader and combine all the filters
  std::map<string, GOAFileType *>::iterator it;
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

PDBParser::GroupOfAtoms *GOALoader::tryAll(const string &fileName,
                                           bool deleteWater) {
  // iterate through each loader and call checkType to determine
  // which loader can load the file
  std::map<string, GOAFileType *>::iterator it;
  for (it = m_ExtensionMap.begin(); it != m_ExtensionMap.end(); ++it) {
    if ((it->second)->checkType(fileName)) // found it
    {
      return (it->second)->loadFile(fileName, deleteWater);
    }
  }
  // didnt find the right loader
  return 0;
}

void GOALoader::addGOAFileType(GOAFileType *type) {
  m_ExtensionMap[type->extension()] = type;
  m_FilterMap[type->filter()] = type;
}

bool GOALoader::isValidExtension(string extension) {
  std::map<string, GOAFileType *>::iterator it;
  for (it = m_ExtensionMap.begin(); it != m_ExtensionMap.end(); ++it) {
    if (strcmp((it->first).c_str(), extension.c_str()) == 0) {
      return true;
    }
  }
  return false;
}
