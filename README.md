# memb-mechanical-properties 


leaflet area compressibility Ka calculation:
------------------------

For multi-component bilayer, if the majority is POPC, POPC will be selected only for the atomNames.txt. 
heightmap.dat is not updated as input file here which is too big (over 100MB) for github.

Package for Calculating Leaflet Ka Moduli from All-Atom MD Simulaions:
Milka Doktorova, package compiled March 2019
Method based on Doktorova et al. 2019, Biophysical Journal 116, 487-502.
If you use the results from the code in a publication, please cite the article above.
Package Link from Milka's paper(doi.org/10.1016/j.bpj.2018.12.016): https://www.dropbox.com/sh/0cag8x7bj0sxebl/AADdrXX7aT0Thinyc2dnWc8ra?dl=0



leaflet tension calculation:
------------------------
I refer to Milka's paper(doi.org/10.1016/j.bpj.2018.12.016) provided in pressure_profiles_NAMD.zip, which is also uploaded in the bilayer_tension/ folder.
We ran 10000 frames for each system to get a pretty convergent pressure profile. And an example of POPC:POPC bilayer is included.
To avoid uploading large files, run_press.out and run_ewald.out are not included. They can be generated when running NAMD post-analysis for lateral pressure
