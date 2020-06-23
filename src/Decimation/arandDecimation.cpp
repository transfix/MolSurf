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
#include <Decimation/arandDecimation.h>

using namespace CVCUtility;

void ArandDecimation::arandNeighbors(int v, set<int> &neighbors)
{
	for(set<int>::iterator iter = vert2tri[v].begin(); iter != vert2tri[v].end(); iter++)
	{
		int tri = *iter;
		if(triangulation[tri].a != v)
		{
			neighbors.insert(triangulation[tri].a);
		}
		if(triangulation[tri].b != v)
		{
			neighbors.insert(triangulation[tri].b);
		}
		if(triangulation[tri].c != v)
		{
			neighbors.insert(triangulation[tri].c);
		}
	}
}

// ordered list of neighbors
void ArandDecimation::arandNeighbors(int v, vector<int> &neighbors)
{
	set<int>::iterator iter = vert2tri[v].begin();
	int tri = *iter;
	int first = triangulation[tri].a;
	if(first == v)
	{
		first = triangulation[tri].c;
	}
	int last = triangulation[tri].b;
	if(last == v)
	{
		last = triangulation[tri].c;
	}
	int current = first;
	int prev    = last;
	bool success;
	do
	{
		neighbors.push_back(current);
		success = false;
		for(set<int>::iterator myit = vert2tri[v].begin(); myit != vert2tri[v].end(); myit++)
		{
			tri = *myit;
			if     ((triangulation[tri].b == current && triangulation[tri].c == v) ||
				(triangulation[tri].c == current && triangulation[tri].b == v))
			{
				if(prev !=  triangulation[tri].a)
				{
					prev = current;
					current = triangulation[tri].a;
					success = true;
					break;
				}
			}
			else if((triangulation[tri].c == current && triangulation[tri].a == v) ||
				(triangulation[tri].a == current && triangulation[tri].c == v))
			{
				if(prev !=  triangulation[tri].b)
				{
					prev = current;
					current = triangulation[tri].b;
					success = true;
					break;
				}
			}
			else if((triangulation[tri].b == current && triangulation[tri].a == v) ||
				(triangulation[tri].a == current && triangulation[tri].b == v))
			{
				if(prev !=  triangulation[tri].c)
				{
					prev = current;
					current = triangulation[tri].c;
					success = true;
					break;
				}
			}
		}
		if(!success)
		{
			 error("arandNeighbors was not successful");
		}
	}
	while(current != first);
}

void ArandDecimation::arandRetriangulate(int v, vector<arandTri> &newTri)
{
	vector<int> neighbors;
	arandNeighbors(v, neighbors);
	// Create and populate directions
	vector<Vector_3> directions;
	for(int i=0; i<neighbors.size(); i++)
	{
		int iplus = (i+1)%neighbors.size();
		Vector_3 direction = vertices[neighbors[iplus]] - vertices[neighbors[i]];
		direction.normalize();
		directions.push_back(direction);
	}
	// ?
	Vector_3 normal = (vertices[neighbors[0]] - vertices[v])^(vertices[neighbors[1]] - vertices[v]);
	int sometri = *(vert2tri[v].begin());
	Vector_3 orientationCheck =
		(vertices[triangulation[sometri].c] - vertices[triangulation[sometri].a])^
		(vertices[triangulation[sometri].b] - vertices[triangulation[sometri].a]);
	while(neighbors.size() > 3)
	{
		int best = -1;
		double bestangle = 500.0;
		normal.normalize();
		for(int i=0; i<directions.size(); i++)
		{
			int iplus = (i+1)%neighbors.size();
			int iplus2 = (i+2)%neighbors.size();
			Vector_3 myNormal = (vertices[neighbors[iplus]] - vertices[neighbors[i]])^(vertices[neighbors[iplus2]] - vertices[neighbors[i]]);
			if(myNormal*normal > 0)
			{
				double angle = PI - acos(directions[i]*directions[iplus]);
				if(angle < bestangle)
				{
					best      = iplus;
					bestangle = angle;
				}
			}
		}
		int bestminus,bestplus;
		bestminus = (best+neighbors.size()-1)%neighbors.size();
		bestplus = (best+1)%neighbors.size();
		Vector_3 orientationCheck2 =
			(vertices[neighbors[bestplus]] - vertices[neighbors[bestminus]])^
			(vertices[neighbors[best]] - vertices[neighbors[bestminus]]);
		arandTri nt;
		nt.a = neighbors[bestminus];
		nt.b = neighbors[best];
		nt.c = neighbors[bestplus];
		if(orientationCheck*orientationCheck2 < 0)
		{
			nt.b = neighbors[bestplus];
			nt.c = neighbors[best];
		}
		newTri.push_back(nt);
		neighbors.erase(neighbors.begin()+best);
		directions.erase(directions.begin()+best);
		bestminus = (best+neighbors.size()-1)%neighbors.size();
		best = best%neighbors.size();
		directions[bestminus] = vertices[neighbors[best]] - vertices[neighbors[bestminus]];
		directions[bestminus].normalize();
	}

	arandTri nt;
	Vector_3 orientationCheck3 =
		(vertices[neighbors[2]] - vertices[neighbors[0]])^
		(vertices[neighbors[1]] - vertices[neighbors[0]]);
	nt.a = neighbors[0];
	nt.b = neighbors[1];
	nt.c = neighbors[2];
	if(orientationCheck*orientationCheck3 < 0)
	{
		nt.b = neighbors[2];
		nt.c = neighbors[1];
	}
	newTri.push_back(nt);
}


Vector_3 ArandDecimation::getNormal(arandTri& tri)
{
	Vector_3 v1,v2;
	v1 = vertices[tri.b] - vertices[tri.a];
	v2 = vertices[tri.c] - vertices[tri.a];
	v1.normalize();
	v2.normalize();
	return v2^v1;
}

double ArandDecimation::minAngle(arandTri& tri)
{
	double thisAngle,myminAngle;
	Vector_3 v1,v2;
	v1 = vertices[tri.b] - vertices[tri.a];
	v2 = vertices[tri.c] - vertices[tri.a];
	v1.normalize();
	v2.normalize();
	myminAngle = acos(v1*v2)*180.0/PI;
	v1 = vertices[tri.a] - vertices[tri.b];
	v2 = vertices[tri.c] - vertices[tri.b];
	v1.normalize();
	v2.normalize();
	thisAngle = acos(v1*v2)*180.0/PI;
	if(thisAngle < myminAngle)
	{
		myminAngle = thisAngle;
	}
	v1 = vertices[tri.a] - vertices[tri.c];
	v2 = vertices[tri.b] - vertices[tri.c];
	v1.normalize();
	v2.normalize();
	thisAngle = acos(v1*v2)*180.0/PI;
	if(thisAngle < myminAngle)
	{
		myminAngle = thisAngle;
	}
	return myminAngle;
}


// Future additions: prioritize the decimation order
void ArandDecimation::arandDecimate(Geoframe& geo, int nIt, float percentage, float min_angle, float max_variation)
{
	int target = percentage*geo.getNumVerts();
	//cout << "arandDecimate --- " << endl;
	//cout << "  minimum allowable angle: " << min_angle << endl;
	//cout << "  maximum normal variation: " << max_variation << endl;
	//cout << "  maximum number of iterations: " << nIt << endl;
	//cout << "  target mesh size: " << target << endl;
	int nV = geo.getNumVerts();
	int nT = geo.getNumTris();
	vert2tri.resize(nV);
	vertices.resize(nV);
	normals.resize(nV);
	canDecimate.resize(nV);
	for(int i=0; i<nV; i++)
	{
		canDecimate[i] = 0;
		Point_3 p(geo.verts[i][0],geo.verts[i][1],geo.verts[i][2]);
		vertices[i] = p;
		Vector_3 n(0.0,0.0,0.0);
		normals[i] = n;
	}
	for(int i=0; i<geo.getNumTris(); i++)
	{
		arandTri tmp;
		tmp.a = geo.triangles[i][0];
		tmp.b = geo.triangles[i][1];
		tmp.c = geo.triangles[i][2];
		triangulation.push_back(tmp);
		vert2tri[tmp.a].insert(i);
		vert2tri[tmp.b].insert(i);
		vert2tri[tmp.c].insert(i);
		triExists.push_back(true);
		normals[tmp.a] += getNormal(tmp);
		normals[tmp.b] += getNormal(tmp);
		normals[tmp.c] += getNormal(tmp);
	}
	for(int i=0; i<nV; i++)
	{
		normals[i].normalize();
	}
	int curVerts = geo.getNumVerts();
	double max_var_cur, min_angle_cur;

	for(int iIt=0; iIt<nIt && target < curVerts; iIt++)
	{
		max_var_cur = max_variation*(double(iIt+1)/(double)nIt);
		if(min_angle < 30)
		{
			min_angle_cur = 30.0 - (double)(iIt+1)*(30.0 - min_angle)/(double)nIt;
		}
		else
		{
			min_angle_cur = min_angle;
		}
		for(int i=0; i<nV; i++)
		{
			// check if this vertex is valid for decimation...
			if(canDecimate[i] >= 0)
			{
				canDecimate[i] = 0;
			}
			if(canDecimate[i] == 0)
			{
				vector<arandTri> newTri;
				arandRetriangulate(i, newTri);
				canDecimate[i] = 1;
				// check if this retriangulation is acceptable...
				// 0. Check normals of all the neighbors. They should all be similar.
				set<int> myNeighbors;
				arandNeighbors(i, myNeighbors);
				for(set<int>::iterator sIt = myNeighbors.begin(); sIt != myNeighbors.end(); sIt++)
				{
					if(acos(normals[i]*normals[*sIt])*180.0/PI > max_var_cur)
					{
						canDecimate[i] = 0;
					}
				}
				// 1. Check if the new triangles have small angles
				for(int iTri=0; iTri < newTri.size() && canDecimate[i] == 1; iTri++)
				{
					double myminAngle = minAngle(newTri[iTri]);
					if(myminAngle < min_angle_cur)
					{
						canDecimate[i] = 0;
					}
				}
			}
		}
		for(int i=0; i<nV; i++)
		{
			if(canDecimate[i] == 1)
			{
				// decimate this vertex
				vector<arandTri> newTri;
				arandRetriangulate(i, newTri);
				// remove decimated triangles from vert2tri...
				for(set<int>::iterator sIt=vert2tri[i].begin(); sIt != vert2tri[i].end(); sIt++)
				{
					int a,b,c;
					a = triangulation[*sIt].a;
					b = triangulation[*sIt].b;
					c = triangulation[*sIt].c;
					if(a != i)
					{
						vert2tri[a].erase(*sIt);
					}
					if(b != i)
					{
						vert2tri[b].erase(*sIt);
					}
					if(c != i)
					{
						vert2tri[c].erase(*sIt);
					}
				}
				int j=0;
				// add these new triangles...
				for(set<int>::iterator sIt=vert2tri[i].begin(); sIt != vert2tri[i].end(); sIt++)
				{
					if(j < newTri.size())
					{
						triangulation[*sIt].a = newTri[j].a;
						triangulation[*sIt].b = newTri[j].b;
						triangulation[*sIt].c = newTri[j].c;
						vert2tri[newTri[j].a].insert(*sIt);
						vert2tri[newTri[j].b].insert(*sIt);
						vert2tri[newTri[j].c].insert(*sIt);
						canDecimate[newTri[j].a] = 0;
						canDecimate[newTri[j].b] = 0;
						canDecimate[newTri[j].c] = 0;
					}
					else
					{
						// this triangle is now unused...
						triExists[*sIt] = false;
					}
					j++;
				}
				curVerts--;
				canDecimate[i] = -1;
				//break;
			}
		}
	}	
}

void ArandDecimation::printRaw(string fname)
{
	map<int,int> vertMap;
	int numVertices=0;
	for(int i=0; i<vertices.size(); i++)
	{
		if(canDecimate[i] >= 0)
		{
			vertMap[i] = numVertices;
			numVertices++;
		}
	}
	int numTriangles=0;
	for(int i=0; i<triangulation.size(); i++)
	{
		if(triExists[i])
		{
			numTriangles++;
		}
	}
	ofstream fout(fname.c_str());
	//cout << "Printing Triangulation: " << numVertices << " vertices, " << numTriangles << " triangles" << endl;
	fout << numVertices << " " << numTriangles << endl;
	for(int i=0; i<vertices.size(); i++)
	{
		if(canDecimate[i] >= 0)
		{
			fout	<< vertices[i].xyz[0] << " "
				<< vertices[i].xyz[1] << " "
				<< vertices[i].xyz[2];
			if(fname[fname.size() - 1] == 'n')
			{
				fout	<< " " << normals[i].xyz[0]
					<< " " << normals[i].xyz[1]
					<< " " << normals[i].xyz[2];
			}
			fout << endl;
		}
	}
	for(int i=0; i<triangulation.size(); i++)
	{
		if(triExists[i])
		{
			fout	<< vertMap[triangulation[i].a] << " "
				<< vertMap[triangulation[i].c] << " "
				<< vertMap[triangulation[i].b] << endl;
		}
	}
	//cout << vertices.size() << " " << triangulation.size() << endl;
}
