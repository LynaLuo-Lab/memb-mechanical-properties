#!/bin/bash
#SBATCH -J symmTOP
#SBATCH --get-user-env
#SBATCH --partition=gpus
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=1500:00:00
 

export TCLLIBPATH="/home/wenjuan/work/protein/Piezo/pure_lipid_calc/Ka_LTF/Ka_LTF/scripts/vmdplugin/membplugin-release-1.1 $TCLLIBPATH"

#source /cm/shared/apps/vmd-1.9.3beta1/bin/vmd.sh

/cm/shared/apps/vmd1.9.3 -dispdev text -e ./get_thickdx.tcl
