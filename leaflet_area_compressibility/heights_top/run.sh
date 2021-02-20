#!/bin/bash -l

export TCLLIBPATH="/home/wenjuan/work/protein/Piezo/pure_lipid_calc/Ka_LTF/Ka_LTF/scripts/vmdplugin/membplugin-release-1.1 $TCLLIBPATH"

vmd -dispdev text -e ./get_thickdx.tcl
