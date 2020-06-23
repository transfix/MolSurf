#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

int gridSize;


using namespace std;

int fnum = 7;
string fnames[] = {"101M", "1GPD", "1ROB", "1SNC", "2ACK", "2NPX", "2POR"};

//int fnum = 1;
//string fnames[] = {"101M"};

double objective(double a_tar, double v_tar, double a, double v)
{
	return a - a_tar  + v - v_tar;
	//return v - v_tar;
	//return a - a_tar;
	//double aerr = a_tar - a;
	//double verr = v - v_tar;
	//return aerr + verr;
}

void getAreaAndVolume(string fn, double& area, double& volume)
{
	// compute area and volume
	string command = "../MolSurf -getArea Temp/"+fn+".raw Temp/area.txt";
	system(command.c_str());
	command = "../MolSurf -getVolume Temp/"+fn+".raw Temp/volume.txt";
	system(command.c_str());
	ifstream fin;
	fin.open("Temp/area.txt");
	fin >> area;
	fin.close();
	fin.open("Temp/volume.txt");
	fin >> volume;
	fin.close();
}


void getAreaAndVolume(string fn, double b, double iso, double& area, double& volume)
{
	ostringstream cmd;
	cmd << "../MolSurf -getSurfaceUsingGaussianBlur "
		<< " PDB/" << fn << ".pdb Temp/" << fn << "_gb.raw " << iso
		<< " " << gridSize << " " << gridSize << " " << gridSize << " "
		<< b << " 0 1 false";
	// compute adaptive grid surface
	string command = cmd.str();
	cout << command << endl;
	system(command.c_str());
	getAreaAndVolume(fn+"_gb", area, volume);
	cout << "Run " << b << " "  << iso << " " << area << " " << volume << endl;
}

void bestParams(string fn, double& b, double& iso)
{
	cout << endl << "Processing file........ " << fn << endl;
	// compute adaptive grid surface
	ostringstream cmd;
	cmd << "../MolSurf -getSurfaceUsingAdaptiveGrid PDB/"+fn
		+".pdb Temp/"+fn+"_ag.raw "
		<< gridSize << " " << gridSize << " " << gridSize;
	// compute adaptive grid surface
	string command = cmd.str();
	//string command = "../MolSurf -getSurfaceUsingAdaptiveGrid PDB/"+fn
	//  +".pdb Temp/"+fn+"_ag.raw 128 128 128";
	system(command.c_str());
	double area, volume,a_tar,v_tar;
	getAreaAndVolume(fn+"_ag", a_tar, v_tar);
	cout << "Target " << a_tar << " " << v_tar << endl;
	double isolarge, isosmall, isom;
	double betalarge, betasmall, beta;
	beta = -1;
	for (int i=0; i<2; i++)
	{
		isolarge = 10;
		isosmall = .1;
		double ofl,ofs, ofm;
		getAreaAndVolume(fn, beta, isolarge, area, volume);
		ofl = objective(a_tar,v_tar,area,volume);
		getAreaAndVolume(fn, beta, isosmall, area, volume);
		ofs = objective(a_tar,v_tar,area,volume);
		if ((ofl > 0.0 && ofs > 0.0) ||
				(ofl < 0.0 && ofs < 0.0))
		{
			cout << "ERROR: Problem initializing parameters" << endl;
			cout << "ERROR: " << ofl << " " << ofs << endl;
		}
		for (int i=0; i<15; i++)
		{
			cout << ofl << " " << ofs << endl;
			isom = (isolarge+isosmall)/2.0;
			getAreaAndVolume(fn, beta, isom, area, volume);
			ofm = objective(a_tar,v_tar,area,volume);
			if (ofm > 0.0)
			{
				if (ofl > 0.0)
				{
					ofl = ofm;
					isolarge = isom;
				}
				else
				{
					ofs = ofm;
					isosmall = isom;
				}
			}
			else
			{
				if (ofl < 0.0)
				{
					ofl = ofm;
					isolarge = isom;
				}
				else
				{
					ofs = ofm;
					isosmall = isom;
				}
			}
		}
		betalarge = -5.0;
		betasmall = -.5;
		getAreaAndVolume(fn, betalarge, isom, area, volume);
		ofl = objective(a_tar,v_tar,area,volume);
		getAreaAndVolume(fn, betasmall, isom, area, volume);
		ofs = objective(a_tar,v_tar,area,volume);
		if ((ofl > 0.0 && ofs > 0.0) ||
				(ofl < 0.0 && ofs < 0.0))
		{
			cout << "ERROR: Problem initializing beta parameters" << endl;
			cout << "ERROR: " << ofl << " " << ofs << endl;
		}
		for (int i=0; i<15; i++)
		{
			cout << ofl << " " << ofs << endl;
			beta = (betalarge+betasmall)/2.0;
			getAreaAndVolume(fn, beta, isom, area, volume);
			ofm = objective(a_tar,v_tar,area,volume);
			if (ofm > 0.0)
			{
				if (ofl > 0.0)
				{
					ofl = ofm;
					betalarge = beta;
				}
				else
				{
					ofs = ofm;
					betasmall = beta;
				}
			}
			else
			{
				if (ofl < 0.0)
				{
					ofl = ofm;
					betalarge = beta;
				}
				else
				{
					ofs = ofm;
					betasmall = beta;
				}
			}
		}
	}
	b = beta;
	iso = isom;
}


int main()
{
	ofstream fout("param_results.txt");
	gridSize = 64;
	for (int i=0; i<fnum; i++)
	{
		double b,iso;
		bestParams(fnames[i],b,iso);
		fout << fnames[i] << " : " << b << " " << iso << endl;
	}
	return 0;
}
