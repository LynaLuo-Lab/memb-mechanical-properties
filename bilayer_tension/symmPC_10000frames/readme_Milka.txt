-------------------------------------------
To run the calculation, need the following:

- trajectory centered so that mean z position of all terminal methyl groups in the bilayer is at z=0 [traj_center.dcd]. Example script for centering the trajectory is included in the scripts/ directory.

- psf file for the trajectory [struct.psf]

- pdb file of the first frame in the traj_center.dcd trajectory [first_frame.pdb]

- stored velocities (optional) [traj_vel.dcd]

- file with simulation box dimensions [box2.txt]. Example script for calculating the box dimensions is included in the scripts/ directory.

- the pressure directory included in this package


-------------------------------------------
In the pressure directory, modify/edit the following:

In press1:

- path to the parameter files
- temperature
- parameter files names
- configuration parameters (DO NOT change nonbondedFreq and fullElectFrequency, leave them set to 1)
- thermostat and barostat parameters
- pressureProfileSlabs (estimate the number of slabs by dividing the mean z position of the simulation box by 0.8 and taking the closest larger odd number)

In get_press:

- path to vmd
- path to namd

In ew/ewald1:

- same as in press1
- set pressureProfileEwaldX, Y and Z to a number that is smaller than half the simulation box length. The larger the number, the more accurate but the more computationally expensive.

In ew/get_press:

- paths to vmd and namd


-------------------------------------------
If no stored velocities exist, modify the following:

- in press1 and ew/ewald1 substitute the line "velocities vel_frame.pdb" with:

temperature $temperature
reinitvels $temperature

- in get_dcd.tcl and ew/get_dcd.tcl, delete the following lines:

animate read dcd ../traj_vel.dcd beg $frame0 end $frame0 waitfor all
set sel [atomselect top "all"]
$sel set {x} [vecscale 20.45482706 [$sel get {x}]]
$sel set {y} [vecscale 20.45482706 [$sel get {y}]]
$sel set {z} [vecscale 20.45482706 [$sel get {z}]]
$sel writepdb vel_frame.pdb


-------------------------------------------
To run the pressure profile calculation, type:

./get_press > log.out & 

If submitting the job to a cluster, modify run.sh as necessary and execute ./run.sh instead.

After the calculation has finished, make a new folder, pressure_analysis, and copy the following two files there:

pressure/run_press.out
pressure/ew/run_ewald.out

cd into pressure_analysis and to process the files:

- run get_press_files
- run write_pressProfile_toFile.m

To get/plot the pressure profile in MATLAB, run the following commands:

>> pp = load('press_profile_lateral');
>> ppmean = get_press(mean(pp));
>> bins = length(pp(1,:));
>> slabs = load('press_slabs.txt');
>> ST = mean(slabs);
>> plot(-ST*(bins-1)/2:ST:ST*(bins-1)/2,ppmean)