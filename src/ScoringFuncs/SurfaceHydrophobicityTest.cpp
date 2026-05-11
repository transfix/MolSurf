
#include <ScoringFuncs/SurfaceHydrophobicity.h>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
  if (argc < 3) {
    cout << "usage: SurfaceHydrophobicityTest <pqr> <quad>" << endl;
  }

  SurfaceHydrophobicity *sh = new SurfaceHydrophobicity(argv[1], argv[2]);

  cout << sh->getSurfaceHydrophobicity() << endl;

  delete sh;

  return 0;
}
