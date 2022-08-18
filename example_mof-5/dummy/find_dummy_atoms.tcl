set ZnX_dist 0.6
set fw [open "X_coor.coor" {w}]

# Load the Zn1 and carboxylate atoms
mol new 1x1x1_Zn1+3x3x3_O1orO2.xyz type {xyz} first 0 last -1 step 1 waitfor 1

set ZnSel [atomselect top "name Zn"]
set numZnAtoms [$ZnSel num]
puts "Number of Zn atoms: $numZnAtoms"
set ZnIndexList [$ZnSel get index]

for {set i 0} {$i < $numZnAtoms} {incr i} {
 #select O2 atoms bonded to every Zn atom
set ZnIndex [lindex $ZnIndexList $i]
set ZnSel [atomselect top "index $ZnIndex"]
set ZnCoor [lindex [$ZnSel get {x y z}] 0]
set OSel [atomselect top "name O and (within 3 of index $ZnIndex)"]
set OCoorList [$OSel get {x y z}]
$ZnSel delete
$OSel delete

foreach OCoor $OCoorList {
 set OZnVec [vecsub $OCoor $ZnCoor]
 set OZnUV [vecscale [expr 1./[veclength $OZnVec]] $OZnVec]
 puts $fw "X [vecadd $ZnCoor [vecscale $ZnX_dist $OZnUV]]"
}

}

close $fw

