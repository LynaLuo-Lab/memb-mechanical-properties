package require membraneheights

heightmap -structure ../step5_assembly.psf -trajectory ../traj.center.dcd -res 8.0 -sel "name P" -leaf 2
#heightmap -structure ../1.psf -trajectory ../1.dcd -res 8.0 -sel "name P" -leaf 0

exit
