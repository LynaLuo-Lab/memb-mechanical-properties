# memb-mechanical-properties 


Leaflet area compressibility (Ka) calculation:
------------------------

For multi-component bilayer, if the majority is POPC, POPC will be selected only for the atomNames.txt. 
heightmap.dat is not updated as input file here which is too big (over 100MB) for github.

Package for Calculating Leaflet Ka Moduli from All-Atom MD Simulaions:
Milka Doktorova, package compiled March 2019
Method based on Doktorova et al. 2019, Biophysical Journal 116, 487-502.
If you use the results from the code in a publication, please cite the article above.
Package Link from Milka's paper(doi.org/10.1016/j.bpj.2018.12.016): https://www.dropbox.com/sh/0cag8x7bj0sxebl/AADdrXX7aT0Thinyc2dnWc8ra?dl=0



Leaflet tension calculation:
------------------------
I refer to Milka's paper(doi.org/10.1016/j.bpj.2018.12.016) provided in pressure_profiles_NAMD.zip, which is also uploaded in the bilayer_tension/ folder.
We ran 10000 frames for each system to get a pretty convergent pressure profile. And an example of POPC:POPC bilayer is included.
To avoid uploading large files, run_press.out and run_ewald.out are not included. They can be generated when running NAMD post-analysis for lateral pressure


Bending rigidity modulus (Kc) calculation:
------------------------
The Kc calculation mainly refers to two papers (DOI: 10.1039/c7cp01921a and DOI 10.1186/s12859-016-1003-z), which use real-space fluctuation analysis of molecular dynamics simulations and present an open source implementation of the method as a set of Python modules using the computational framework OpenStructure. The modules are freely available through GitHub at https://github.com/njohner/ost_pymodules/ while OpenStructure can be obtained at http://www.openstructure.org. My calculation of Kc is using docker platform of OpenStructure with Johner's lipid analysis scripts(ost_pymodules).

Here input and output pdb and dcd files are not uploaded here for size issues. Basically, input pdb and dcd files can be prepared according to ost_pymodules example folder. The Kc_lower_pip2_1000frames_modified folder calculates only Kc for lower leaflet of POPC:PIP2 bilayer system using 1000 frames with my modified lipid_analysis.py script, which generates Nij and Ntot for the number of pairs of POPC-POPC, POPC-PIP2 and PIP2-PIP2. Kc_lower_pip2_1000frames_standard folder calculates both tilt (Kt) and bending rigidity modulus (Kc) but without printing the number of pairs (Nij and Ntot); A lower.pdb file of POPC:PIP2 bilayer system is also included for an example of the lipid name selection in pip2_elastic_properties.py


To run the code, python3 pip2_elastic_properties.py 


Shear viscosity calculation:
------------------------
The main process 1) run shear deformation simulation in GROMACS 2016 version, the shear.mdp file is uploaded which includes the deformation = [0, 0, 0, rz, 0, 0 ], which deforms in xy-plane and gives gradient along y-axis. 2) calculation <Pxy> and use Eq.7 from our paper to calculate surface shear viscosity. Details about this method can be found in paper: DOI: 10.1021/acs.jctc.9b00683
