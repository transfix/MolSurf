/*
  Copyright 2011 The University of Texas at Austin

        Author: Muhibur Rasheed <muhib@ices.utexas.edu>
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
#include<Hausdorff/Hausdorff.h>
#include <iostream>

bool hausdorff(int argc, char ** argv) 
{
	if (argc < 5) 
	{
		cout << "Usage: MolSurf -hausdorff <input raw/rawn file 1> <input raw/rawn file 2> <1=raw / 2 = rawn>" << endl;
		return false;
	}

	bool RAWN;
	int r = atoi(argv[4]);
	if(r==2) RAWN = true;
	else if(r==1) RAWN = false;
	else
	{

		cout << "Usage: MolSurf -hausdorff <input raw/rawn file 1> <input raw/rawn file 2> <1=raw / 2 = rawn>" << endl;
		return false;
	}

	HausdorffDistance H(argv[2],argv[3],RAWN);
  
//	cout<<"Hausdorff distance = "<< H.getHausdorffDistance() <<endl;
//	cout<<"Average distance = "<< H.getAverageDistance() <<endl; 

	cout<<H.getHausdorffDistance()<< " " << H.getAverageDistance() << endl;

	return true; 
}
