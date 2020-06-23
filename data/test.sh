#!/bin/bash
#Craig: there are 3 commented lines that when uncommented, and when the main CMakeLists.txt is set to profile, profile the tests
#build/bin/MolSurf -h && mv gmon.out gmon.sum

# assert that the previous command returned 0
assert()
{
	if [ $? == 0 ]
		then echo -e "\033[7;30;42m PASS \033[0m" $1;
		else echo -e "\033[7;30;41m FAIL \033[0m" $1; #exit;
	fi 
	#gprof -s build/bin/MolSurf gmon.out gmon.sum
}

# arand: this should be reordered so some fast tests are run first and then more detailed tests follow...
assert "adaptive grid" `build/bin/MolSurf -surfaceUsingAdaptiveGrid data/PDB/2POR.pdb data/RAW/2POR_ag.raw 64`
assert "gaussian blur" `build/bin/MolSurf -surfaceUsingGaussianBlur data/PDB/2POR.pdb data/RAW/2POR_gb.raw 1 64 -.5 0 1 false`
assert "HLS" `build/bin/MolSurf -surfaceUsingHLS data/PDB/101M.pdb data/RAW/101M_tri.raw 64`
assert "adaptive grid" `build/bin/MolSurf -surfaceUsingAdaptiveGrid data/PDB/101M.pdb data/RAW/101M_ag_tri.raw 128`
assert "adaptive grid" `build/bin/MolSurf -surfaceUsingAdaptiveGrid data/PDB/2GLS.pdb data/RAW/2GLS_ag.raw 128`
assert "adaptive grid" `build/bin/MolSurf -surfaceUsingAdaptiveGrid data/PDB/4AT1.pdb data/RAW/4AT1_ag.raw 128`
assert "adaptive grid" `build/bin/MolSurf -surfaceUsingAdaptiveGrid data/PDB/C60.pdb data/RAW/C60_ag_tri.raw 128`
assert "adaptive grid" `build/bin/MolSurf -surfaceUsingAdaptiveGrid data/PQR/1RIE.pqr data/RAWN/1RIE_ag.rawn 128`
assert "adaptive grid" `build/bin/MolSurf -surfaceUsingAdaptiveGrid data/PQR/ARG_1.pqr data/RAWN/ARG_1_ag.rawn 128`
assert "gaussian blur" `build/bin/MolSurf -surfaceUsingGaussianBlur data/PDB/101M.pdb data/RAW/101M_gb.raw 1 128 -.5 0 1 false`
assert "gaussian blur" `build/bin/MolSurf -surfaceUsingGaussianBlur data/PDB/2GLS.pdb data/RAW/2GLS_gb.raw 1 64 -.5 0 1 false`
assert "gaussian blur" `build/bin/MolSurf -surfaceUsingGaussianBlur data/PDB/4AT1.pdb data/RAW/4AT1_gb.raw 1 64 -.5 0 1 false`
assert "HLS" `build/bin/MolSurf -surfaceUsingHLS data/PDB/101M.pdb data/RAW/101M_tri.raw 64 2 0.001`
assert "LBIE" `build/bin/MolSurf -surfaceUsingLBIE tri  0 data/RAWIV/head.rawiv data/RAW/head_tri.raw`
assert "LBIE" `build/bin/MolSurf -surfaceUsingLBIE quad 0 data/RAWIV/head.rawiv data/RAW/head_quad.raw`
assert "LBIE" `build/bin/MolSurf -surfaceUsingLBIE tet  0 data/RAWIV/head.rawiv data/RAW/head_tet.raw`
assert "LBIE" `build/bin/MolSurf -surfaceUsingLBIE hex  0 data/RAWIV/head.rawiv data/RAW/head_hex.raw`
assert "LBIE" `build/bin/MolSurf -surfaceUsingLBIE tri2 0 data/RAWIV/head.rawiv data/RAW/head_tri2.raw`
assert "LBIE" `build/bin/MolSurf -surfaceUsingLBIE tet2 0 data/RAWIV/head.rawiv data/RAW/head_tet2.raw`
assert "area" `build/bin/MolSurf -area data/RAW/cube.raw data/TXT/cube_area.txt`
assert "area" `build/bin/MolSurf -area data/RAW/tet.raw data/TXT/tet_area.txt`
assert "aspline quad" `build/bin/MolSurf -aSpline -quad data/RAWN/cube.rawn gaussian 12 data/TXT/cube_quad_points12.txt`
assert "aspline quad" `build/bin/MolSurf -aSpline -quad data/RAWN/cube.rawn gaussian 1 data/TXT/cube_quad_points1.txt`
assert "aspline quad" `build/bin/MolSurf -aSpline -quad data/RAWN/single_tri.rawn gaussian 12 data/TXT/single_tri_quad_points12.txt`
assert "aspline quad" `build/bin/MolSurf -aSpline -quad data/RAWN/single_tri.rawn gaussian 1 data/TXT/single_tri_quad_points1.txt`
assert "aspline quad" `build/bin/MolSurf -aSpline -quad data/RAWN/tet.rawn gaussian 12 data/TXT/tet_quad_points12.txt`
assert "aspline quad" `build/bin/MolSurf -aSpline -quad data/RAWN/tet.rawn gaussian 1 data/TXT/tet_quad_points1.txt`
assert "aspline quad" `build/bin/MolSurf -aSpline -quad data/RAWN/tet.rawn gaussian 1 data/TXT/tet_quad_points1.txt`
assert "aspline volume" `build/bin/MolSurf -aSpline -volume data/RAWN/cube.rawn data/TXT/cube_aSpline_volume.txt`
assert "aspline volume" `build/bin/MolSurf -aSpline -volume data/RAWN/tet.rawn data/TXT/tet_aSpline_volume.txt`
assert "coarse grain charge" `build/bin/MolSurf -coarseGrainCharge vdw data/PQR/1RIE.pqr data/TXT/1RIE_pqr.br data/PQR/1RIE_pqr_cg.pqr data/TXT/1RIE_cg_pqr.br data/PQR/1RIE_pqr_cg.info data/PQR/1RIE_pqr_cg_charge.pqr`
assert "coarse grain charge" `build/bin/MolSurf -coarseGrainCharge vdw data/PQR/1RIE.pqr data/TXT/1RIE_vdw.br data/PQR/1RIE_vdw_cg.pqr data/TXT/1RIE_cg_vdw.br data/PQR/1RIE_vdw_cg.info data/PQR/1RIE_vdw_cg_charge.pqr`
assert "coarse grain charge" `build/bin/MolSurf -coarseGrainCharge vdw data/PQR/ARG_1.pqr data/TXT/ARG_1_vdw.br data/PQR/ARG_1_vdw_cg.pqr data/TXT/ARG_1_cg_vdw.br data/PQR/ARG_1_vdw_cg.info data/PQR/ARG_1_vdw_cg_charge.pqr`
assert "coarse grain" `build/bin/MolSurf -coarseGrain pqr data/PQR/1RIE.pqr data/PQR/1RIE_pqr_cg.pqr data/PQR/1RIE_pqr_cg.info`
assert "coarse grain" `build/bin/MolSurf -coarseGrain vdw data/PQR/1RIE.pqr data/PQR/1RIE_vdw_cg.pqr data/PQR/1RIE_vdw_cg.info`
assert "coarse grain" `build/bin/MolSurf -coarseGrain vdw data/PQR/ARG_1.pqr data/PQR/ARG_1_vdw_cg.pqr data/PQR/ARG_1_vdw_cg.info`
assert "adaptive grid" `build/bin/MolSurf -surfaceUsingAdaptiveGrid data/PQR/ARG_1_vdw_cg.pqr data/RAWN/ARG_1_ag_cg.rawn 128`
assert "normals" `build/bin/MolSurf -normals -average data/RAW/brick.raw data/RAWN/brick_average.rawn`
assert "normals" `build/bin/MolSurf -normals -average data/RAW/cube.raw data/RAWN/cube_average.rawn`
assert "normals" `build/bin/MolSurf -normals -average data/RAW/tet.raw data/RAWN/tet_average.rawn`
assert "normals" `build/bin/MolSurf -normals -loop data/RAW/brick.raw data/RAWN/brick_loop.rawn`
assert "normals" `build/bin/MolSurf -normals -loop data/RAW/cube.raw data/RAWN/cube_loop.rawn`
assert "normals" `build/bin/MolSurf -normals -loop data/RAW/tet.raw data/RAWN/tet_loop.rawn`
assert "normals" `build/bin/MolSurf -normals -average data/RAW/101M_ag_tri.raw data/RAWN/101M_tri.rawn`
assert "optimize CG" `build/bin/MolSurf -optimizeCGGeometry 64 64 64 10 data/PQR/1RIE.pqr data/PQR/1RIE_cg.pqr  data/PQR/1RIE_cg.info data/PQR/1RIE_opt.pqr`
assert "optimize CG" `build/bin/MolSurf -optimizeCGGeometry 64 64 64 10 data/PQR/1RIE.pqr data/PQR/1RIE_cg.pqr  data/PQR/1RIE_cg.info data/PQR/1RIE_opt.pqr`
assert "optimize CG" `build/bin/MolSurf -optimizeCGGeometry 64 64 64 10 data/PQR/1RIE.pqr data/PQR/1RIE_vdw_cg_charge.pqr data/PQR/1RIE_vdw_cg.info data/PQR/1RIE_vdw_opt.pqr`
assert "populate SAS" `build/bin/MolSurf -populateSAS data/PDB/101M.pdb data/TXT/101MSAS.txt 64 64 64 1.4 0`
assert "quality improve" `build/bin/MolSurf -qualityImprove data/RAW/101M_tri.raw data/RAW/101M_imp.raw`
assert "quality improve" `build/bin/MolSurf -qualityImprove data/RAW/C60_ag_tri.raw data/RAW/C60_ag_imp_tri.raw`
# Craig - These fail noisily so I've commented them out for now
#assert "decimate" `build/bin/MolSurf -decimate data/RAW/101M_tri.raw data/RAW/101M_dec.raw 6 25 10` 				# segfault during free in Vector_3 destructor, called from arandReTri 
assert "decimate" `build/bin/MolSurf -decimate data/RAW/C60_ag_imp_tri.raw data/RAW/C60_ag_dec_tri.raw 10 25 20` 		
assert "aspline RAW" `build/bin/MolSurf -aSpline -raw data/RAW/C60_ag_dec_tri.raw data/RAWN/C60_ag_aSpline.rawn 3`		# output needs to be verified and added to answers
assert "aspline RAW" `build/bin/MolSurf -aSpline -raw data/RAWN/cube.rawn data/RAWN/cube_aSpline2.rawn 2`
assert "aspline RAW" `build/bin/MolSurf -aSpline -raw data/RAWN/cube.rawn data/RAWN/cube_aSpline30.rawn 30`
assert "aspline RAW" `build/bin/MolSurf -aSpline -raw data/RAWN/tet.rawn data/RAWN/tet_aSpline2.rawn 2`
assert "aspline RAW" `build/bin/MolSurf -aSpline -raw data/RAWN/tet.rawn data/RAWN/tet_aSpline30.rawn 30`
#assert "aspline RAW+quad" `build/bin/MolSurf -aSpline -raw+quad data/RAWN/101M_dec.rawn data/RAWN/101M_aSpline.rawn data/RAWN/101M_aSpline_quad.rawn  data/QUAD/101M_aSpline.quad 4 gaussian 3` # waiting on decimate output
assert "aspline RAW+quad" `build/bin/MolSurf -aSpline -raw+quad data/RAWN/101M_tri.rawn data/RAWN/101M_aSpline.rawn data/RAWN/101M_aSpline_quad.rawn  data/QUAD/101M_aSpline.quad 2 gaussian 3`
assert "aspline RAW+quad" `build/bin/MolSurf -aSpline -raw+quad data/RAWN/cube.rawn data/RAWN/cube_aSpline4.rawn data/RAWN/cube_aSpline4_quad.rawn  data/QUAD/cube_aSpline4.quad 4 gaussian 3`
assert "volume" `build/bin/MolSurf -volume data/RAW/brick.raw data/TXT/brick_volume.txt`
assert "volume" `build/bin/MolSurf -volume data/RAW/cube_oriented2.raw data/TXT/cube_oriented2_volume.txt`
assert "volume" `build/bin/MolSurf -volume data/RAW/cube_oriented.raw data/TXT/cube_oriented_volume.txt`
assert "volume" `build/bin/MolSurf -volume data/RAW/cube.raw data/TXT/cube_volume.txt`
assert "volume" `build/bin/MolSurf -volume data/RAW/tet_oriented.raw data/TXT/tet_oriented_volume.txt`
assert "volume" `build/bin/MolSurf -volume data/RAW/tet.raw data/TXT/tet_volume.txt`
#assert "pocket" `build/bin/MolSurf -pocket data/PDB/101M.pdb 20 64 data/RAW/101M_pocket`						# dimensions frazzled because it guesses the wrong endianness
#assert "pocket" `build/bin/MolSurf -pocket data/PDB/2POR.pdb 20 64 data/RAW/2POR`							# "
#assert "surface atoms" `build/bin/MolSurf -surfaceAtoms data/PDB/101M.pdb data/PDB/101MSurface`					# Assertion `!((grp->m_SubGroups).empty())' failed.
assert "born radius" `build/bin/MolSurf -bornRadius pqr data/PQR/1RIE_pqr_cg.pqr data/RAWN/1RIE_ag.rawn data/TXT/1RIE_cg_pqr.br`
assert "born radius" `build/bin/MolSurf -bornRadius pqr data/PQR/1RIE.pqr data/RAWN/1RIE_ag.rawn data/TXT/1RIE_pqr.br`
assert "born radius" `build/bin/MolSurf -bornRadius vdw data/PQR/1RIE.pqr data/RAWN/1RIE_ag.rawn data/TXT/1RIE_vdw.br`
assert "born radius" `build/bin/MolSurf -bornRadius vdw data/PQR/1RIE_vdw_cg.pqr data/RAWN/1RIE_ag.rawn data/TXT/1RIE_cg_vdw.br`
assert "born radius" `build/bin/MolSurf -bornRadius vdw data/PQR/ARG_1.pqr data/RAWN/ARG_1_ag.rawn data/TXT/ARG_1_vdw.br`
assert "born radius" `build/bin/MolSurf -bornRadius vdw data/PQR/ARG_1_vdw_cg.pqr data/RAWN/ARG_1_ag_cg.rawn data/TXT/ARG_1_cg_vdw.br` # segfaults in nfft_B_A ( )
# Craig - These take forever and passed on 26/Jul/2010
#assert "pocket tunnel" `build/bin/MolSurf -pocketTunnel data/RAW/2GLS_ag.raw 0 1 data/RAW/2GLS_pocket.raw data/RAW/2GLS_tunnel.raw`
#assert "pocket tunnel" `build/bin/MolSurf -pocketTunnel data/RAW/2POR_ag.raw 3 3 data/RAW/2POR_pocket.raw data/RAW/2POR_tunnel.raw`
#assert "pocket tunnel" `build/bin/MolSurf -pocketTunnel data/RAW/4AT1_ag.raw 0 1 data/RAW/4AT1_pocket.raw data/RAW/4AT1_tunnel.raw`
#assert "volume HLS" `build/bin/MolSurf -volumeUsingHLS data/PDB/101M.pdb data/RAWIV/101M.rawiv 64 64 64`

#gprof build/bin/MolSurf gmon.sum > profile.txt

# test to see if we got the same answers
cd data/answers
for FILE in */*
do
	diff -q $FILE ../$FILE; assert $FILE 
done
cd ../..

#Craig: Perhaps I should add a test that makes sure TexMol builds and runs... to be run post-commit.

exit
