--------------------------------------------------------------------------------------
Package for Calculating Leaflet Ka Moduli from All-Atom MD Simulaions

Milka Doktorova, package compiled March 2019

Method based on Doktorova et al. 2019, Biophysical Journal 116, 487-502

If you use the results from the code in a publication, please cite the article above.
--------------------------------------------------------------------------------------

--------------------------------------------------------------------------------------
Calculating leaflet and bilayer Ka
--------------------------------------------------------------------------------------
1. Extract the contents of Ka_LTF.zip

2. cd into Ka_LTF

3. Place the trajectory for analysis (with the bilayer centered in the simulation box)
in the directory and call it traj.dcd

4. Place the structure file (for all-atom simulation, the .psf file) in the directory
and call it struct.psf

5. Create the atomNames.txt file for each leaflet (see file description below and
example files in the example directory)

6. Specify the correct path to VMD's executable in run.sh in heights_top and
heights_bot directories

7. Execute the calculation of interpolated heights by running run.sh, e.g.
./run.sh > log.out &

8. Create two new directories, Ka_top and Ka_bot

9. Copy heightmap.dat for each leaflet in the respective Ka_xxx directory

10. Copy atomNames.txt for each leaflet in the respective Ka_xxx directory, then
open the file, insert empty lines between the groups of atoms (headgroup, chain 1,
chain 2, etc.) and save it as atomGroups.txt (see example files in example/)

11. In each Ka_xxx directory, create apl.txt and temperature.txt files containing
the average area per lipid in the leaflet and the temperature of the simulation in
Kelvin, respectively.

12. cd into each Ka_xxx directory and run calculate_Ka.m script

13. To calculate bilayer Ka from the leaflet moduli, run get_bilayerKa_stats.m

--------------------------------------------------------------------------------------
Files
--------------------------------------------------------------------------------------
heights_top, heights_bot
	directories for the two leaflets in which interpolated heights are calculated
	with VMD

atomNames.txt
	list of atoms whose z-positions are to be analyzed
	every line contains the name of 1 atom or a selection of atoms
	whose z-pos will be interpolated together as a single surface
	(see examples in example/atomNames/)
	
	NOTE: File should not end with an empty line.
	NOTE: Always have at least one headgroup atom such as the phosphate 'P'

atomGroups.txt
	same as atomNames.txt file but with empty lines separating the different sets of atoms (headgroup, chain 1, chain 2, etc.)
	
get_thickdx.tcl
	tcl script for the analysis, uses the customized membraneheights package
	from MEMBPLUGIN (in Ka_LTF/scripts/vmdplugin)
	In addition to all available input/output options provided by MEMBPLUGIN
	(https://sourceforge.net/p/membplugin/wiki/Home/)
	user can specify the following Ka-specific parameters:
	-sel 	 -- atoms to use to distinguish between the two leaflets and define grid boundaries (using x/y min/max values). Default: "name P"
  	-res 	 -- grid resolution (Default: 8.0)
        -o	 -- output file path without extension. A tabulated file is generated. Default: heightmap.dat
  	-leaf    -- leaflet selection, 1 (top) or 0 (bot). Default: 1
  	-rad	 -- radius of interpolation. Default: 150.0 (assuming a small box size, consider all atoms)

run.sh
	executes get_thickdx.tcl by starting VMD in a text-only mode
	
	!!!IMPORTANT!!! Verify correct path to VMD executable

apl.txt
	average area per lipid in the leaflet

temperature.txt
	temperature of the simulation in Kelvin

heightmap.dat
	output file from get_thickdx.tcl script. Has the format:
	{frame #} {grid_x} {grid_y} {interpolated z-pos of all atom selections in atomNames.txt separated by space}
