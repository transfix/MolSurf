/*
  Copyright 2011 The University of Texas at Austin

	Authors: Alex Rand <arand@ices.utexas.edu>
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
#include <SimpleCuration/SimpleCuration.h>
#include <SimpleCuration/SimpleCurationAPI.h>
#include <GeometryFileTypes/GeometryLoader.h>
#include <Geometry/Geometry.h>

void usageRemoveInteriorPockets()
{
	cout	<< endl << "Usage: MolSurf -removeInteriorPockets <raw input file> <raw output file>"
		<< endl << "    Remove small disconnected components of the mesh."
		<< endl ;
}

bool removeInteriorPockets(int argc, char * argv[]) {

  if(argc != 4) {
    usageRemoveInteriorPockets();
    return false;
  }
  
  // load the geometry
  Geometry* geometry = GeometryLoader().loadFile(argv[2]);
  // run the code
  geometry = RemoveInteriorPockets(geometry);
  // write the output
  GeometryLoader().saveFile(argv[3], geometry);
  
  return true;
}
