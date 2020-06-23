/*
  Copyright 2011 The University of Texas at Austin

	Authors: Muhibur Rasheed <muhibur@ices.utexas.edu>
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
#include <InterfaceStats/CapriCriteria.h>
#include <UsefulMath/Matrix.h>

int main()
{
	Matrix mtx;

	mtx.set(0, 0, 0.975); 	mtx.set(0, 1, 0.147); 	mtx.set(0, 2, 0.168); 	mtx.set(0, 3, -3.947);

	mtx.set(1, 0, -0.109);	mtx.set(1, 1, 0.970);	mtx.set(1, 2, -0.219);	mtx.set(1, 3, -1.498);

	mtx.set(2, 0, -0.195);	mtx.set(2, 1, 0.195);	mtx.set(2, 2, 0.961);	mtx.set(2, 3, 0.590);

	mtx.set(3, 0, 0);	mtx.set(3, 1, 0);	mtx.set(3, 2, 0);	mtx.set(3, 3, 1);

	CapriCriteria *capri = new CapriCriteria( "1SYX_R_U.pqr", "1SYX_L_U.pqr", mtx );

	if(capri->computeCapriScore())
	{
		cout<<"L-RMSD = "<<capri->getLRMSD()<<endl;
		cout<<"I-RMSD = "<<capri->getIRMSD()<<endl;
		cout<<"F_NAT = "<<capri->getFnat()<<endl;
		cout<<"F_NONNAT = "<<capri->getFnonnat()<<endl;
		cout<<"VERDICT = "<<capri->getVerdict()<<endl;
	}

	CapriCriteria *capri2 = new CapriCriteria( "1SYX_R_U.pqr", "1SYX_L_U.pqr", "1SYX_L_U.pqr" );

	if(capri2->computeCapriScore())
	{
		cout<<"L-RMSD = "<<capri2->getLRMSD()<<endl;
		cout<<"I-RMSD = "<<capri2->getIRMSD()<<endl;
		cout<<"F_NAT = "<<capri2->getFnat()<<endl;
		cout<<"F_NONNAT = "<<capri2->getFnonnat()<<endl;
		cout<<"VERDICT = "<<capri2->getVerdict()<<endl;
	}

	
	return 1;
}
