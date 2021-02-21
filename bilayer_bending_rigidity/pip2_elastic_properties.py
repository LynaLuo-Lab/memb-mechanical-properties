"""
Example written by Niklaus Johner (niklaus.johner@a3.epfl.ch)
This is an example of how to calculate the membrane elastic properties from a simulation.
We start here with an aligned trajectory as making the alignment (see the align_trajectory.py example)
is slow.
This example should be run from within the "example" directory, otherwise the path
to the python modules and other files will have to be set.
Everything that gets generated will be put in the "tmp" directory
"""
from ost import *
import os,sys
import numpy as npy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

indir="./data"
outdir="./tmp_lower_pip2"
sys.path.append("../")
import lipid_analysis,trajectory_utilities

#First we load the structure and trajectory
p=io.IOProfile(processor=conop.HeuristicProcessor(),dialect='CHARMM',fault_tolerant=True)
#eh=io.LoadPDB(os.path.join(indir,"symm_1stFrame_A.pdb"),profile=p)
#t=io.LoadCHARMMTraj(eh,os.path.join(indir,'symm_last1000f.dcd'),stride=4)
eh=io.LoadPDB(os.path.join(indir,"lower.pdb"),profile=p)
t=io.LoadCHARMMTraj(eh,os.path.join(indir,'5th_lower.dcd'),stride=5)
# eh=io.LoadPDB(os.path.join(indir,"sample.pdb"),profile=p)
# t=io.LoadCHARMMTraj(eh,os.path.join(indir,'sample.dcd'),stride=1)


####################################
# As this is a coarse grained trajectory, atoms are not properly recognised so that masses
# were not set when loading the structure. We need to set the masses and element to make sure
# that the algorithm calculating the density used to determine the membrane interface will work properly.
# for a in eh.atoms:
#  a.SetElement('Ge')#Has the correct atomic weight
#  a.SetRadius(2.35)
#  a.SetMass(72)

####################################
# For this example the trajectory used was already aligned, but it contains only one unit cell
# To properly treat boundary conditions we replicate the unit cell all around the original one.
# And then wrap it around the central unit cell
# See the aign_trajectory.py example for more information.

extension_directions=[[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]]
t=trajectory_utilities.ExtendTrajectoryToNeighboringUnitCells(t,extension_directions,(2,2,2))
eh=t.GetEntity()
centers=mol.alg.AnalyzeCenterOfMassPos(t,eh.Select("cname=A"))
trajectory_utilities.WrapTrajectoryInPeriodicCell(t,centers)
#print centers
#We save the trajectory for visualization
#io.SaveCHARMMTraj(t,os.path.join(outdir,"extended.pdb"),os.path.join(outdir,"extended.dcd"),profile=p)

###################################
# We prepare the dictionaries with information about the selections defining headgroups and tails
# for the different lipids. They will be used in the calculation of tilts and splays.

# 1. The system in the example has two types of lipids: cholesterol (CHO) and DOPE (DOP)
#    We need a list with residue names of all lipids that will be analysed
lipid_names=['OPC','API']

# 2. Define the residue name of the water molecules.
#    This is used together with the lipid names to calculate the density maps
#    Of water and lipids used to extract the membrane interface and orient the normals on the interface.
water_name='WAT,POT,CLA'

# 3. For each lipid type, we define the selection for headgroups and tails and the selection
#    that will be used to calculate distances between lipids. This should be atoms lying at
#    the neutral plane.
head_group_dict={'OPC':'aname=P,C2','API':'aname=P,C2'}
tail_dict={'OPC':'aname=6C21,7C21,8C21,4C31,5C31,6C31','API':'aname=0C22,9C21,8C21,6C31,7C31,8C31'}
distance_sele_dict={'OPC':'aname=C21,C22,C23,C31,C32,C33','API':'aname=C21,C22,C23,C31,C32,C33'}


# 4. We define the different cutoffs used by the algorithms. More information can be found
#    In the documentation
#max distance between two lipids to be considered in the splay calculation
distance_cutoff=10.0 
#max tilt angle with respect to normal for splay calculations
angle_cutoff=0.175 
# radius of area used to calculate the normals to the membrane interface. If too small, normals will be noisy.
within_size_normals=10.0
#Density cutoff used when calculating the interface.
#This can be set to 0, but calculations will be slower.
density_cutoff=0.3
#Stride used when calculating the water and lipid densities.
#When set to 1, all frames are considered, which can be slow.
density_stride=1

#5. Periodic boundaries are best treated by replicating the simulation box around
#   the original unit cell and then calculating tilts and splays only for lipids from
#   the original unit cell, but using the surrounding unit cells to find all neighbors 
#   of a lipid for the splay calculation and to avoid boundary effects on the interfaces
#   and hence on the normal vectors used both for the tilt and splay calculations.
#   We therefore need to tell the function which lipids belong to the central unit cell.
#   This is done by setting a bool property for the corresponding residues.

tilt_bool_prop='do_tilt'
splay_bool_prop='do_splay'
v=eh.Select('cname=A and rname=OPC,API')
for r in v.residues:
  r.SetBoolProp(tilt_bool_prop,True)
  r.SetBoolProp(splay_bool_prop,True)


v=eh.Select('cname!=A and rname=OPC,API')
for r in v.residues:
  r.SetBoolProp(tilt_bool_prop,False)
  r.SetBoolProp(splay_bool_prop,False)


#6. Other parameters
prot_sele=None
sele_dict={}
filename_basis='tilt_splay_'
#7. We calculate the tilts and splays:
(lipid_tilt_dict,lipid_normal_dict,splay_dict,b_eh)=lipid_analysis.AnalyzeLipidTiltAndSplay(t,
lipid_names,head_group_dict,tail_dict,distance_cutoff,within_size_normals,distance_sele_dict,water_name,   
outdir,density_cutoff,prot_sele,density_stride,tilt_bool_prop,splay_bool_prop,filename_basis,sele_dict)


########################################################
# Now we analyze the lipid tilts and splays, fitting the corresponding analytical functions
# to extract the elastic constants. See documentation and papers cited therein for more information.
# This can be done with a single function call to ExtractTiltAndSplayModuli.
# The function will make one fit for each lipid type, and then calculate the overall tilt modulus
# By taking a weighted average of the individual contributions.

# The number of bins used when making the histograms for the tilts
nbins=100

# In this example the area per lipid is 49.5 A^2. For a lipid bilayer the area can be obtained
# Using the lipid_analysis.AnalyzeAreaPerLipid function
# symm area per lipid=64.4
lipid_area=64.3


#Now we extract the constants
k_dict=lipid_analysis.ExtractTiltAndSplayModuli(lipid_tilt_dict,splay_dict,lipid_area,outdir,nbins)
