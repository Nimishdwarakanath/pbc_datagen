MODULE DataStructureModule
 IMPLICIT NONE
 TYPE :: FullAtomicData
  INTEGER                :: AtomicTypeIndex
  CHARACTER ( LEN = 10 ) :: AtomicType
  CHARACTER ( LEN = 3 )  :: AtomicSymbol
  INTEGER                :: AtomicNumber
  REAL ( KIND = 8 )      :: AtomicMass
  REAL ( KIND = 8 )      :: AtomicCharge
  REAL ( KIND = 8 )      :: AtomicvdWRadius
  REAL ( KIND = 8 )      :: AtomicCoordinates  (3)
  REAL ( KIND = 8 )      :: AtomicScaledCoordinates (3)
  INTEGER                :: NBonds
  INTEGER                :: NAngles
  INTEGER                :: NDihedrals
  INTEGER                :: NImpropers
 ENDTYPE FullAtomicData

 TYPE :: AtomicType
  INTEGER                :: AtomicTypeIndex
  CHARACTER ( LEN = 10 ) :: AtomicType
  CHARACTER ( LEN = 3 )  :: AtomicSymbol
  REAL ( KIND = 8 )      :: AtomicvdWRadius
 ENDTYPE AtomicType

 TYPE :: FullBondData
  INTEGER                :: BondType
  INTEGER                :: BondAtom(2)
  REAL ( KIND = 8 )      :: BondDistance
 ENDTYPE FullBondData

 TYPE :: BondType
  INTEGER                :: BondAtomType(2)
  REAL ( KIND = 8 )      :: BondDistance
 ENDTYPE BondType

 TYPE :: FullAngleData
  INTEGER                :: AngleType
  INTEGER                :: AngleAtom(3)
  REAL ( KIND = 8 )      :: AngleAngle
 ENDTYPE FullAngleData

 TYPE :: AngleType
  INTEGER                :: AngleAtomType(3)
  REAL ( KIND = 8 )      :: AngleAngle
 ENDTYPE AngleType

 TYPE :: FullDihedralData
  INTEGER                :: DihedralType
  INTEGER                :: DihedralAtom(4)
  REAL ( KIND = 8 )      :: DihedralAngle
 ENDTYPE FullDihedralData

 TYPE :: DihedralType
  INTEGER                :: DihedralAtomType(4)
  REAL ( KIND = 8 )      :: DihedralAngle
 ENDTYPE DihedralType

 TYPE :: FullImproperData
  INTEGER                :: ImproperType
  INTEGER                :: ImproperAtom(4)
  REAL ( KIND = 8 )      :: ImproperAngle
 ENDTYPE FullImproperData

 TYPE :: ImproperType
  INTEGER                :: ImproperAtomType(4)
  REAL ( KIND = 8 )      :: ImproperAngle
 ENDTYPE ImproperType

ENDMODULE DataStructureModule
