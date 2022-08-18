bash correct_bonding.sh > bonding.tcl
vmd -m ../{bonded,non-bonded,dummy}/coor.xyz

#then source bonding.tcl

pbc set [exec cat parm] -molid 0
