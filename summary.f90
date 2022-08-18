!----------------------------------------------------------------------------------
! atomic_data.f90 
!----------------------------------------------------------------------------------
1:MODULE AtomicDataModule
63:!  FUNCTION SymboltoAtomicNumber ()
!----------------------------------------------------------------------------------
! calculate.f90 
!----------------------------------------------------------------------------------
1:MODULE CalculateModule
8:SUBROUTINE PopulateAtomicData ()
25: CALL CalculateBoxQuantities() !calculate HMatrix , IHMatrix and LammpsBoxParameters
61: SUBROUTINE GetBondList()
79:  CALL GetBondType ()
82: FUNCTION IdentifyAngleFromBonds ( i , j ) RESULT (res)
114: SUBROUTINE GetAngleList()
133:  CALL GetAngleType ()
136: FUNCTION IdentifyDihedralFromAngles ( i , j ) RESULT (res)
181: SUBROUTINE GetDihedralList()
200:  CALL GetDihedralType ()
203: FUNCTION IdentifyImproperFromAngles ( i , j ) RESULT (res)
241: SUBROUTINE GetImproperList()
286:  CALL GetImproperType ()
290: SUBROUTINE GetBondType ()
323: SUBROUTINE GetAngleType ()
354: SUBROUTINE GetDihedralType ()
386: SUBROUTINE GetImproperType ()
419: FUNCTION SameArray ( i , j ) RESULT ( res )
433: FUNCTION SameImproper ( i , j ) RESULT ( res )
451: SUBROUTINE NumberOfBADIForEachAtom ()
!----------------------------------------------------------------------------------
! data_structure.f90 
!----------------------------------------------------------------------------------
1:MODULE DataStructureModule
!----------------------------------------------------------------------------------
! global_data.f90 
!----------------------------------------------------------------------------------
1:MODULE GlobalDataModule
2:!Module contains Common data to be shared by the main program and other modules
20: TYPE (FullAtomicData) , ALLOCATABLE :: AtomicData (:)
21: TYPE (FullBondData) , ALLOCATABLE :: BondData (:)
22: TYPE (FullAngleData) , ALLOCATABLE :: AngleData (:)
23: TYPE (FullDihedralData) , ALLOCATABLE :: DihedralData (:)
24: TYPE (FullImproperData) , ALLOCATABLE :: ImproperData (:) , ImproperData2 (:)
26: TYPE (BondType) , ALLOCATABLE :: BondTypeData (:)
27: TYPE (AngleType) , ALLOCATABLE :: AngleTypeData (:)
28: TYPE (DihedralType) , ALLOCATABLE :: DihedralTypeData (:)
29: TYPE (ImproperType) , ALLOCATABLE :: ImproperTypeData (:)
!----------------------------------------------------------------------------------
! io.f90 
!----------------------------------------------------------------------------------
1:MODULE IOModule
5: SUBROUTINE ReadInputFile ()
40: SUBROUTINE WriteInputFileToLog ()
68: SUBROUTINE WriteAtomicDataToLog ()
79: SUBROUTINE WriteLammpsDataFile ()
218: SUBROUTINE VerboseData ()
231:   CALL WriteHeaderDataFile ( NBonds )
254:   CALL WriteHeaderDataFile ( NBonds )
282:   CALL WriteHeaderDataFile ( NBonds )
315:   CALL WriteHeaderDataFile ( NBonds )
344: SUBROUTINE WriteHeaderDataFile (NBonds)
!----------------------------------------------------------------------------------
! main.f90 
!----------------------------------------------------------------------------------
9: CALL ReadInputFile ()
11: CALL WriteInputFileToLog ()
13: CALL PopulateAtomicData ()
17:  CALL GetBondList ()
21:  CALL GetAngleList ()
25:  CALL GetDihedralList ()
29:  CALL GetImproperList ()
32: CALL NumberOfBADIForEachAtom () !!number of bonds, angles, dihedrals, and impropers for each atom
33: CALL WriteAtomicDataToLog () 
34: CALL WriteLammpsDataFile () !number of bonds, angles, dihedrals, and impropers for each atom must be written
35: CALL VerboseData ()
!----------------------------------------------------------------------------------
! pbc.f90 
!----------------------------------------------------------------------------------
1:MODULE PBCModule
2:!Module contains Common data to be shared by the main program and other modules
8: FUNCTION CrossProduct ( a , b) RESULT ( res )
16: SUBROUTINE CalculateBoxQuantities()
70: FUNCTION PBCDistance ( i , j ) RESULT ( res )
79: FUNCTION GetMinImage ( j , i ) RESULT ( res ) 
90: FUNCTION CalculateAngle ( a , b ) RESULT ( res )
102: FUNCTION PBCAngle ( i , j , k ) RESULT (res)
114: FUNCTION PBCDihedralAngle ( i , j , k , l) RESULT (res)
131: FUNCTION PBCImproperAngle ( i , j , k , l) RESULT (res)
