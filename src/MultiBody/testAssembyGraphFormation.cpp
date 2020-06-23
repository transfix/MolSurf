
#include <MultiBody/assembler.h>
#include <iostream>

using namespace std;


int main(int argc, char** argv)
{
	if(argc < 5)
		cout<<"Usage testAssembyGraphFormation <nodeDataFile> <edgeDataFile> <nodeWeightSettingsFile> <edgeWeightSettingsFile>\n";

	Assembler* assemb = new Assembler(true);

	assemb->makeGraph(argv[1], argv[2], argv[3], argv[4]);

	assemb->assemble();
//	assemb->print();
}


