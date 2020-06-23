#ifndef _XUDECIMATE_H_
#define _XUDECIMATE_H_


void CallXuDecimate(
	int niter,              // number of iteration
	float controllor,       // ? (*PI radian measure)
	float variation,        // ? (*PI radian measure)
	float innerangle,       // ? (*PI radian measure)
	int showpoly,           // show poly or not (doesn't matter)
	int grouping,           // ?? some flag
	int retriways,          // if (retriways  == 1) MinimalAngleConnection
	// if (retriways  == 2) MinimalDistanceConnection
	// if (retriways  == 3) MinimalPathConnection;
	int savenorm,            // save file flag
	const char* fname
);

void generateNormal(int item_no);

#endif
