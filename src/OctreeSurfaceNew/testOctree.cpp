
#include <OctreeSurfaceNew/Octree.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>

double getTime( void )
{
#ifdef _WIN32
   time_t ltime;
   _timeb tstruct;
   time( &ltime );
   _ftime( &tstruct );
   return ( double ) ( ltime + 1e-3 * ( tstruct.millitm ) );
#else
   struct timeval t;
   gettimeofday( &t, NULL );
   return ( double )( t.tv_sec + 1e-6 * t.tv_usec );
#endif
}

int main(int argc, char ** argv) 
{
	if (argc < 7) 
	{
		cout << "Usage: Octree <input pqr file> <output raw file> <minDepth> <maxDepth> <isovalue> <blobbyness>" << endl;
	}

	double startT = getTime();

	Octree T (argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atof(argv[5]), atof(argv[6]) );

	double endT = getTime();

	cout<<"Total Time "<<endT-startT<<endl;
}
