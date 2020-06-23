#/bin/bash

rm -rf MolSurfExamples
rm MolSurfExamples.tgz

mkdir MolSurfExamples
cd MolSurfExamples

mkdir PDB
mkdir RAW
mkdir RAWN
mkdir TXT

cp ../MolSurfExamples.sh ./
cp ../PDB/2POR.pdb PDB

cd ..

tar czvf MolSurfExamples.tgz MolSurfExamples