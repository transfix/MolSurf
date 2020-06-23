
../bin/MolSurf -getSurfaceUsingAdaptiveGrid PDB/$1.pdb RAW/$1_ag_tri.raw $2 $2 $2
../bin/MolSurf -geoImprove RAW/$1_ag_tri.raw RAW/$1_ag_imp_tri.raw
../bin/MolSurf -decimate RAW/$1_ag_imp_tri.raw RAW/$1_ag_dec_tri.raw 5 20 15
#../bin/MolSurf -geoImprove RAW/$1_ag_dec_tri.raw RAW/$1_ag_imp_tri.raw
#../bin/MolSurf -decimate RAW/$1_ag_imp_tri.raw RAW/$1_ag_dec_tri.raw 10 20 15






#../bin/MolSurf -aspline -raw RAW/$1_ag_dec_tri.raw RAWN/$1_ag_aspline.rawn 3
