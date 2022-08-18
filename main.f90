PROGRAM main
 USE GlobalDataModule
 USE CalculateModule
 USE IOModule

 IMPLICIT NONE
 !Read input

 CALL ReadInputFile ()

 CALL WriteInputFileToLog ()

 CALL PopulateAtomicData ()

 CALL WritePDBFile ()

 IF ( WantBonds ) THEN
  CALL GetBondList ()
 ENDIF

 IF ( WantAngles ) THEN
  CALL GetAngleList ()
 ENDIF

 IF ( WantDihedrals ) THEN
  CALL GetDihedralList ()
 ENDIF

 IF ( WantImpropers ) THEN
  CALL GetImproperList ()
 ENDIF
 CALL NumberOfBADIForEachAtom () !!number of bonds, angles, dihedrals, and impropers for each atom
 CALL WriteAtomicDataToLog () 
 CALL WriteLammpsDataFile () !number of bonds, angles, dihedrals, and impropers for each atom must be written
 CALL VerboseData ()

!read atomtypes
!read box parameters
!alterations to the vdW radii

 !identify the bonds
  !calculate the bond lengths
 !identify the angles from bond information
  !calculate the angles
 !identify the dihedrals from angle information
  !calculate the dihedral angles
 !identify the improper from bond information
  !calculate the improper dihedral angles
  !keep the ones which satisfy the criterion
 !CLOSE (61)
ENDPROGRAM main
