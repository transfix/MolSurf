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
#include <Pocket/Pocket.h>

using namespace CVCUtility;

Pocket::Pocket()
{
}

Pocket::~Pocket()
{
}

// ...from a GOA
void Pocket::constructPockets(PDBParser::GroupOfAtoms* goa, string fileName, float distance, Geometry** pocketGeometry, int size)
{
	Geometry* geometry = MolecularSurface::getMolecularSurface(goa, size);
	constructPockets(geometry, fileName, distance, pocketGeometry, size);
}

// ...from a Geometry
void Pocket::constructPockets(Geometry* geometry, string fileName, float distance, Geometry** pocketGeometry, int size)
{
	// Generate our filenames
	string saveFileName	= fileName + ".rawiv";
	string surfFileName	= fileName + "_ms.raw";
	string sdfMsFile	= fileName + "_ms.rawiv";
	string pocketFileName	= fileName + "_pks.raw";
	string pocketFileNameiv	= fileName + "_pks.rawiv";
	string shellFileName	= fileName + "_shell.raw";

	// Save our geometry and work with it from disk.. why?
	GeometryLoader().saveFile(surfFileName, geometry);
	delete geometry;

	// Read the molecular surface and construct SDF volume V1
	FaceVertSet3D fvs;
	Geom3DParser().ParseRawFile(fvs, surfFileName.c_str());
	int idim[3] = {size, size, size};
	DistanceTransform dt1(fvs, idim, distance);
	dt1.transform();
	dt1.writeRawiv(sdfMsFile.c_str());
	SimpleVolumeData* pV1 = VolumeLoader().loadFile(sdfMsFile);

	// Extract and save the shell
	Geometry* pShell = SimpleVolumeDataIsocontourer::getIsocontour(pV1, distance);
	GeometryLoader().saveFile(shellFileName, pShell);

	// Read the shell surface and construct SDF volume V2
	FaceVertSet3D fvs2;
	Geom3DParser().ParseRawFile(fvs2, shellFileName.c_str());
	DistanceTransform dt2(fvs2, dt1.getReg3Data());
	dt2.transform();
	dt2.writeRawiv("v2.rawiv");
	SimpleVolumeData* pV2 = VolumeLoader().loadFile("v2.rawiv");

	// Take the minium of the values in pV1 and pV2
	uint nv = size * size * size;
	for(uint i = 0; i < nv; i++)
	{
		double f1 = pV1->getValueAt(0, i);
		double f2 = pV2->getValueAt(0, i) - distance;
		f1 = minimum(f1, f2);
		pV1->setValueAt(0, i, f1);
	}

	// Extract and save pocket as level set
	(*pocketGeometry) = SimpleVolumeDataIsocontourer::getIsocontour(pV1, 1.5);
	GeometryLoader().saveFile(pocketFileName.c_str(), (*pocketGeometry));

	// Save and exit
	VolumeLoader().saveFile(pocketFileNameiv, pV1);
	delete pShell;
	delete pV2;
	delete pV1;
}
