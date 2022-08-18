MODULE GlobalDataModule
!Module contains Common data to be shared by the main program and other modules
 USE DataStructureModule
 IMPLICIT NONE
 REAL( KIND = 8 ), PARAMETER :: degtorad = DATAN( 1.0D0 )/45.0D0 , radtodeg = 1.0D0/degtorad
 CHARACTER ( LEN = 100 ) :: XYZFileName , ParameterFileName , OutputFilesPrefix
 LOGICAL                 :: WantBonds , WantAngles , WantDihedrals , WantImpropers
 REAL ( KIND = 8 )       :: BondLengthTolerance , AngleAngleTolerance , DihedralAngleTolerance , ImproperAngleTolerance , &
                            ImproperConsiderationTolerance
 INTEGER                 :: TotalNumberOfAtomTypes , NumberOfManualvdWRadiusChanges
 CHARACTER ( LEN = 5 ) , ALLOCATABLE &
                         :: TypeAndSymbolList ( : , : )
 REAL ( KIND = 8 ) , ALLOCATABLE &
                         :: vdWRadiusList ( : )

 CHARACTER ( LEN = 100 ) :: VariableFileName , OutputLogFileName, LammpsDataFileName , LammpsParameterFileName

 INTEGER :: TotalNumberOfAtoms

 TYPE (FullAtomicData) , ALLOCATABLE :: AtomicData (:)
 TYPE (FullBondData) , ALLOCATABLE :: BondData (:)
 TYPE (FullAngleData) , ALLOCATABLE :: AngleData (:)
 TYPE (FullDihedralData) , ALLOCATABLE :: DihedralData (:)
 TYPE (FullImproperData) , ALLOCATABLE :: ImproperData (:) , ImproperData2 (:)

 TYPE (BondType) , ALLOCATABLE :: BondTypeData (:)
 TYPE (AngleType) , ALLOCATABLE :: AngleTypeData (:)
 TYPE (DihedralType) , ALLOCATABLE :: DihedralTypeData (:)
 TYPE (ImproperType) , ALLOCATABLE :: ImproperTypeData (:)

 REAL ( KIND = 8 )       :: BoxParameters (6) , LammpsBoxParameters ( 9 )
 REAL ( KIND = 8 ) , DIMENSION ( 3 , 3 ) :: HMatrix, IHMatrix
 INTEGER :: TotalNumberOfBonds , TotalNumberOfAngles , TotalNumberOfDihedrals , TotalNumberOfImpropers, TotalNumberOfImpropers2
 INTEGER :: TotalNumberOfBondTypes , TotalNumberOfAngleTypes , TotalNumberOfDihedralTypes , TotalNumberOfImproperTypes

ENDMODULE GlobalDataModule
