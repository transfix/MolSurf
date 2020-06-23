#/bin/bash

# set line below to point to the the MolSurf executable
MOLSURF=/net/electron/workspace/arand/Svn/MolSurf/bin/MolSurf

# below are a number of example calls to the MolSurf program

${MOLSURF} -surfaceUsingAdaptiveGrid PDB/2POR.pdb RAW/2POR_ag.raw 128

${MOLSURF} -surfaceUsingGaussianBlur PDB/2POR.pdb RAW/2POR_gb.raw 1 128 -.5 0 1 false

${MOLSURF} -surfaceUsingHLS PDB/2POR.pdb RAW/2POR_ls.raw 128

${MOLSURF} -normals -average RAW/2POR_ag.raw RAWN/2POR_ag.rawn


${MOLSURF} -area RAW/2POR_ag.raw TXT/2POR_ag.area.txt
${MOLSURF} -volume RAW/2POR_ag.raw TXT/2POR_ag.volume.txt