MODULE CalculateModule
 USE GlobalDataModule
 USE AtomicDataModule
 USE PBCModule
 IMPLICIT NONE
 CONTAINS

SUBROUTINE PopulateAtomicData ()
 INTEGER                 :: i , j
!Read the XYZ, parameter files and populate AtomicData
 OPEN ( 52 , FILE = TRIM (ParameterFileName) , ACTION = 'READ' )
 READ (52, * ) BoxParameters
 CLOSE (52)

 OPEN ( 51 , FILE = TRIM (XYZFileName) , ACTION = 'READ' )
  READ ( 51, * ) TotalNumberOfAtoms
  READ ( 51, * ) 
  !Read AtomicType and AtomicCoordinates for each atom
  ALLOCATE ( AtomicData ( TotalNumberOfAtoms ) )
  DO i = 1 , TotalNumberOfAtoms
   READ ( 51 , * ) AtomicData(i)%AtomicType, AtomicData(i)%AtomicCoordinates(:)
  ENDDO
 CLOSE ( 51 )

 CALL CalculateBoxQuantities() !calculate HMatrix , IHMatrix and LammpsBoxParameters

 DO i = 1 , TotalNumberOfAtoms
  AtomicData(i)%AtomicScaledCoordinates(:) = MATMUL ( IHMatrix , AtomicData(i)%AtomicCoordinates )
 ENDDO

 DO i = 1 , TotalNumberOfAtoms
  !Find AtomicTypeIndex, AtomicSymbol for each atom
  DO j = 1 , TotalNumberOfAtomTypes
   IF ( TRIM ( AtomicData(i)%AtomicType ) .EQ. TRIM (TypeAndSymbolList ( j , 1 ) ) ) THEN
    AtomicData(i)%AtomicTypeIndex = j
    AtomicData(i)%AtomicSymbol = TypeAndSymbolList ( j , 2 )
    EXIT
   ENDIF
  ENDDO

  !Find AtomicNumber for each atom
  DO j = 1 , NumberOfElements
   IF ( TRIM ( AtomicData(i)%AtomicSymbol ) .EQ. TRIM ( AtomicSymbolArray (j) ) ) AtomicData(i)%AtomicNumber = j
  ENDDO

  !Assign AtomicvdWRadius and AtomicMass for each atom
  IF ( vdWRadiusList (AtomicData(i)%AtomicTypeIndex) .GT. 1.0D+1) THEN
   AtomicData(i)%AtomicvdWRadius = vdWRadiusArray (AtomicData(i)%AtomicNumber)
  ELSE
   AtomicData(i)%AtomicvdWRadius = vdWRadiusList (AtomicData(i)%AtomicTypeIndex)
  ENDIF

    AtomicData(i)%AtomicMass = AtomicMassArray(AtomicData(i)%AtomicNumber)

 ENDDO

 !DEALLOCATE ( vdWRadiusList , TypeAndSymbolList )
ENDSUBROUTINE PopulateAtomicData


 SUBROUTINE GetBondList()
  INTEGER :: i , j
  REAL ( KIND = 8 ) :: distance
  TotalNumberOfBonds = 0
  ALLOCATE (BondData ( 20 * TotalNumberOfAtoms ) )
  BondData%BondType = 0
  !Uses o(N^2) algorithm.
  DO i = 1 , TotalNumberOfAtoms - 1
   DO j = i + 1 , TotalNumberOfAtoms
     distance = PBCDistance ( i , j )
    IF ( distance .LT. 0.6D0 * ( AtomicData(i)%AtomicvdWRadius + AtomicData(j)%AtomicvdWRadius) ) THEN
     TotalNumberOfBonds = TotalNumberOfBonds + 1
     BondData(TotalNumberOfBonds)%BondAtom(1) = i
     BondData(TotalNumberOfBonds)%BondAtom(2) = j
     BondData(TotalNumberOfBonds)%BondDistance = distance
    ENDIF
   ENDDO
  ENDDO
  CALL GetBondType ()
 ENDSUBROUTINE GetBondList

 FUNCTION IdentifyAngleFromBonds ( i , j ) RESULT (res)
  INTEGER , INTENT ( IN ) :: i , j !ith and jth bond data
  INTEGER  :: res ( 3 ) , t
  
  IF (BondData(i)%BondAtom(1) .EQ. BondData(j)%BondAtom(1)) THEN
   res (1) = BondData(i)%BondAtom(2)
   res (2) = BondData(i)%BondAtom(1)
   res (3) = BondData(j)%BondAtom(2)
  ELSEIF (BondData(i)%BondAtom(1) .EQ. BondData(j)%BondAtom(2)) THEN
   res (1) = BondData(i)%BondAtom(2)
   res (2) = BondData(i)%BondAtom(1)
   res (3) = BondData(j)%BondAtom(1)
  ELSEIF (BondData(i)%BondAtom(2) .EQ. BondData(j)%BondAtom(1)) THEN
   res (1) = BondData(i)%BondAtom(1)
   res (2) = BondData(i)%BondAtom(2)
   res (3) = BondData(j)%BondAtom(2)
  ELSEIF (BondData(i)%BondAtom(2) .EQ. BondData(j)%BondAtom(2)) THEN
   res (1) = BondData(i)%BondAtom(1)
   res (2) = BondData(i)%BondAtom(2)
   res (3) = BondData(j)%BondAtom(1)
  ELSE
   res = 0
  ENDIF

  IF ( res (3) .LT. res (1) )THEN
   t = res ( 3 )
   res ( 3 ) = res ( 1 )
   res ( 1 ) = t
  ENDIF
   
 ENDFUNCTION IdentifyAngleFromBonds

 SUBROUTINE GetAngleList()
  INTEGER :: i , j 
  INTEGER :: angle (3)
  ALLOCATE ( AngleData ( TotalNumberOfAtoms * 20 ) )
  TotalNumberOfAngles = 0
  DO i = 1 , TotalNumberOfBonds - 1
   DO j = i + 1 , TotalNumberOfBonds
    angle = IdentifyAngleFromBonds ( i , j)
    IF ( .NOT. ALL ( angle .EQ. 0 )) THEN
    TotalNumberOfAngles = TotalNumberOfAngles + 1
    AngleData(TotalNumberOfAngles)%AngleAtom(1) = angle (1)
    AngleData(TotalNumberOfAngles)%AngleAtom(2) = angle (2)
    AngleData(TotalNumberOfAngles)%AngleAtom(3) = angle (3)
    ENDIF
   ENDDO
  ENDDO
  DO i = 1 , TotalNumberOfAngles
   AngleData(i)%AngleAngle = PBCAngle ( AngleData(i)%AngleAtom(1),AngleData(i)%AngleAtom(2),AngleData(i)%AngleAtom(3) )
  ENDDO
  CALL GetAngleType ()
 ENDSUBROUTINE GetAngleList

 FUNCTION IdentifyDihedralFromAngles ( i , j ) RESULT (res)
  INTEGER , INTENT ( IN ) :: i , j !ith and jth angle data
  INTEGER  :: res ( 4 ) , t
  REAL (KIND=8) :: anglei, anglej

  anglei = degtorad * PBCAngle ( AngleData(i)%AngleAtom(1),AngleData(i)%AngleAtom(2),AngleData(i)%AngleAtom(3) )
  anglej = degtorad * PBCAngle ( AngleData(j)%AngleAtom(1),AngleData(j)%AngleAtom(2),AngleData(j)%AngleAtom(3) )

  IF ( &
      ( ABS ( DCOS (anglei) ) .GT. 0.985D0 )  &
      .OR. &
      ( ABS ( DCOS (anglej) ) .GT. 0.985D0 )  &
     )THEN
   res = 0
   RETURN
  ENDIF
 
  IF ( ( AngleData(i)%AngleAtom(2) == AngleData(j)%AngleAtom(3) ) &
                           .AND. &
       ( AngleData(i)%AngleAtom(3) == AngleData(j)%AngleAtom(2) ) ) THEN
   res ( 1 ) = AngleData(i)%AngleAtom(1)
   res ( 2 ) = AngleData(i)%AngleAtom(2)
   res ( 3 ) = AngleData(i)%AngleAtom(3)
   res ( 4 ) = AngleData(j)%AngleAtom(1)
  ELSEIF ( ( AngleData(i)%AngleAtom(2) == AngleData(j)%AngleAtom(1) ) &
                           .AND. &
       ( AngleData(i)%AngleAtom(3) == AngleData(j)%AngleAtom(2) ) ) THEN
   res ( 1 ) = AngleData(i)%AngleAtom(1)
   res ( 2 ) = AngleData(i)%AngleAtom(2)
   res ( 3 ) = AngleData(i)%AngleAtom(3)
   res ( 4 ) = AngleData(j)%AngleAtom(3)
  ELSEIF ( ( AngleData(i)%AngleAtom(2) == AngleData(j)%AngleAtom(3) ) &
                           .AND. &
       ( AngleData(i)%AngleAtom(1) == AngleData(j)%AngleAtom(2) ) ) THEN
   res ( 1 ) = AngleData(i)%AngleAtom(3)
   res ( 2 ) = AngleData(i)%AngleAtom(2)
   res ( 3 ) = AngleData(i)%AngleAtom(1)
   res ( 4 ) = AngleData(j)%AngleAtom(1)
  ELSEIF ( ( AngleData(i)%AngleAtom(2) == AngleData(j)%AngleAtom(1) ) &
                           .AND. &
       ( AngleData(i)%AngleAtom(1) == AngleData(j)%AngleAtom(2) ) ) THEN
   res ( 1 ) = AngleData(i)%AngleAtom(3)
   res ( 2 ) = AngleData(i)%AngleAtom(2)
   res ( 3 ) = AngleData(i)%AngleAtom(1)
   res ( 4 ) = AngleData(j)%AngleAtom(3)
  ELSE
   res = 0
   RETURN
  ENDIF

!  IF ( res (4) .LT. res (1) )THEN
!   t = res ( 4 )
!   res ( 4 ) = res ( 1 )
!   res ( 1 ) = t
!   t = res ( 2 )
!   res ( 2 ) = res ( 3 )
!   res ( 3 ) = t
!  ENDIF
 ENDFUNCTION IdentifyDihedralFromAngles

 SUBROUTINE GetDihedralList()
  INTEGER :: i , j 
  INTEGER :: dihedral (4)
  ALLOCATE ( DihedralData ( TotalNumberOfAtoms * 20 ) )
  TotalNumberOfDihedrals = 0
  DO i = 1 , TotalNumberOfAngles - 1
   DO j = i + 1 , TotalNumberOfAngles
    dihedral = IdentifyDihedralFromAngles ( i , j)
    IF ( .NOT. ANY ( dihedral .EQ. 0 )) THEN
     TotalNumberOfDihedrals = TotalNumberOfDihedrals + 1
     DihedralData(TotalNumberOfDihedrals)%DihedralAtom = dihedral
    ENDIF
   ENDDO
  ENDDO
  DO i = 1 , TotalNumberOfDihedrals
   DihedralData(i)%DihedralAngle = PBCDihedralAngle ( DihedralData(i)%DihedralAtom(1), &
                                   DihedralData(i)%DihedralAtom(2),DihedralData(i)%DihedralAtom(3), &
                                   DihedralData(i)%DihedralAtom(4) )
  ENDDO
  CALL GetDihedralType ()
 ENDSUBROUTINE GetDihedralList

 FUNCTION IdentifyImproperFromBonds ( i , j , k ) RESULT (res)
  INTEGER , INTENT ( IN ) :: i , j , k !ith, jth and kth bond data
  INTEGER  :: res ( 4 ) , t
  INTEGER , DIMENSION(2) :: bondi , bondj, bondk
  INTEGER , DIMENSION(3) :: plane123 , plane423
  INTEGER :: ii , jj , kk
  REAL (KIND=8) :: angle123 , angle423

  bondi = BondData(i)%BondAtom(:) 
  bondj = BondData(j)%BondAtom(:) 
  bondk = BondData(k)%BondAtom(:) 

  IF ( ALL (bondi .EQ. bondj) .OR. ALL ( bondj .EQ. bondk ) .OR. ALL (bondk .EQ. bondi)) THEN
   res = 0
   RETURN
  ENDIF

  !selection of improper torsion based purely on connectivity
  res = 0
  DO ii = 1 , 2
   DO jj = 1 , 2
    IF ( bondi(ii) .EQ. bondj (jj) ) THEN
     DO kk = 1 , 2
      IF ( bondi(ii) .EQ. bondk(kk) ) THEN
       res=(/ bondi(ii) , bondi(3-ii) , bondj(3-jj) , bondk(3-kk) /)
       plane123 = (/ bondi(ii) , bondi(3-ii) , bondj(3-jj)/)
       plane423 = (/ bondk(3-kk) , bondi(3-ii) , bondj(3-jj)/)
       res (1:3) = plane123 (1:3)
       res (4)   = plane423 (1)
       ! discard improper if atoms in plane 123 or 423 are colinear
       angle123 = degtorad * PBCAngle ( plane123(1) ,plane123(2) ,plane123(3) )
       angle423 = degtorad * PBCAngle ( plane423(1) ,plane423(2) ,plane423(3) )
       IF ( &
            (ABS ( DCOS (angle123) ) .GT. 0.985D0) &
             .OR. &
            (ABS ( DCOS (angle423) ) .GT. 0.985D0) &
          ) THEN
           res = 0
           RETURN
       ENDIF
      ENDIF
     ENDDO
    ENDIF
   ENDDO
  ENDDO



  !calculate the angle between two planes and return the improper if it's within the tolerance, else return 0
  !call this function from the required segment
 ENDFUNCTION IdentifyImproperFromBonds

! FUNCTION IdentifyImproperFromAngles ( i , j ) RESULT (res)
!  INTEGER , INTENT ( IN ) :: i , j !ith and jth angle data
!  INTEGER  :: res ( 4 ) , t
!  REAL (KIND=8) :: anglei, anglej
!
!  anglei = degtorad * PBCAngle ( AngleData(i)%AngleAtom(1),AngleData(i)%AngleAtom(2),AngleData(i)%AngleAtom(3) )
!  anglej = degtorad * PBCAngle ( AngleData(j)%AngleAtom(1),AngleData(j)%AngleAtom(2),AngleData(j)%AngleAtom(3) )
!
!  IF ( &
!      ( ABS ( DCOS (anglei) ) .GT. 0.985D0 )  &
!      .OR. &
!      ( ABS ( DCOS (anglej) ) .GT. 0.985D0 )  &
!     )THEN
!   res = 0
!   RETURN
!  ENDIF
! 
!
!      IF ( ( AngleData(i)%AngleAtom(2) == AngleData(j)%AngleAtom(2) ) &
!                               .AND. &
!           ( AngleData(i)%AngleAtom(3) == AngleData(j)%AngleAtom(3) ) ) THEN
!   res ( 1 ) = AngleData(i)%AngleAtom(2) !i2
!   res ( 2 ) = AngleData(j)%AngleAtom(1) !j1
!   res ( 3 ) = AngleData(i)%AngleAtom(3) !i3
!   res ( 4 ) = AngleData(i)%AngleAtom(1) !i1
!  ELSEIF ( ( AngleData(i)%AngleAtom(2) == AngleData(j)%AngleAtom(2) ) &
!                           .AND. &
!       ( AngleData(i)%AngleAtom(3) == AngleData(j)%AngleAtom(1) ) ) THEN
!   res ( 1 ) = AngleData(i)%AngleAtom(2) !i2
!   res ( 2 ) = AngleData(j)%AngleAtom(3) !j3
!   res ( 3 ) = AngleData(i)%AngleAtom(3) !i3
!   res ( 4 ) = AngleData(i)%AngleAtom(1) !i1
!  ELSEIF ( ( AngleData(i)%AngleAtom(2) == AngleData(j)%AngleAtom(2) ) &
!                           .AND. &
!       ( AngleData(i)%AngleAtom(1) == AngleData(j)%AngleAtom(3) ) ) THEN
!   res ( 1 ) = AngleData(i)%AngleAtom(2) !i2
!   res ( 2 ) = AngleData(j)%AngleAtom(1) !j1
!   res ( 3 ) = AngleData(i)%AngleAtom(1) !i1
!   res ( 4 ) = AngleData(i)%AngleAtom(3) !i3
!  ELSEIF ( ( AngleData(i)%AngleAtom(2) == AngleData(j)%AngleAtom(2) ) &
!                           .AND. &
!       ( AngleData(i)%AngleAtom(1) == AngleData(j)%AngleAtom(1) ) ) THEN
!   res ( 1 ) = AngleData(i)%AngleAtom(2) !i2
!   res ( 2 ) = AngleData(j)%AngleAtom(3) !j3
!   res ( 3 ) = AngleData(i)%AngleAtom(1) !i1
!   res ( 4 ) = AngleData(i)%AngleAtom(3) !i3
!  ELSE
!   res = 0
!  ENDIF
! ENDFUNCTION IdentifyImproperFromAngles

 ! Appears to output repetitions
 SUBROUTINE GetImproperList()
  INTEGER :: i , j , k
  INTEGER :: improper (4) , i4 (4) , j4 (4)
  LOGICAL :: Reassign
  REAL ( KIND = 8 ) :: anglei, anglej
  ALLOCATE ( ImproperData2 ( TotalNumberOfAtoms * 2000 ) )
  ! using impropers from bonds
  !TotalNumberOfImpropers2 = 0
  !DO i = 1 , TotalNumberOfAngles - 1
  ! DO j = i + 1 , TotalNumberOfAngles
  !  improper = IdentifyImproperFromAngles ( i , j)
  !  IF ( .NOT. ANY ( improper .EQ. 0 )) THEN
  !   TotalNumberOfImpropers2 = TotalNumberOfImpropers2 + 1
  !   ImproperData2(TotalNumberOfImpropers2)%ImproperAtom = improper
  !  ENDIF
  ! ENDDO
  !ENDDO
  TotalNumberOfImpropers2 = 0
  DO i = 1 , TotalNumberOfBonds
   DO j = 1 , TotalNumberOfBonds
    DO k = 1 , TotalNumberOfBonds
     improper = IdentifyImproperFromBonds ( i , j , k)
     IF ( .NOT. ANY ( improper .EQ. 0 )) THEN
      TotalNumberOfImpropers2 = TotalNumberOfImpropers2 + 1
      ImproperData2(TotalNumberOfImpropers2)%ImproperAtom = improper
     ENDIF
    ENDDO
   ENDDO
  ENDDO

  DO i = 1 , TotalNumberOfImpropers2
   ImproperData2(i)%ImproperAngle = PBCImproperAngle ( ImproperData2(i)%ImproperAtom(1), &
                                   ImproperData2(i)%ImproperAtom(2),ImproperData2(i)%ImproperAtom(3), &
                                   ImproperData2(i)%ImproperAtom(4) )
  ENDDO
  !remove repetitions
  ALLOCATE ( ImproperData ( SIZE ( ImproperData2 ) ) )
  TotalNumberOfImpropers = 0
  DO i = 1 , TotalNumberOfImpropers2
   i4 = ImproperData2(i)%ImproperAtom
   anglei = ImproperData2(i)%ImproperAngle
   Reassign = .TRUE.
   DO j = 1 , TotalNumberOfImpropers
    j4 = ImproperData(j)%ImproperAtom
    anglej = ImproperData(j)%ImproperAngle
    IF ( SameImproper ( i4 , j4 ) )THEN
     Reassign = .FALSE.
     EXIT
    ENDIF
   ENDDO
   IF ( Reassign .AND. &
   ( ( DABS (anglei - 180.0D0 ) .LT. ImproperConsiderationTolerance  )  &
                            .OR.&
     ( DABS (anglei - 0.0D0 ) .LT. ImproperConsiderationTolerance  ) ) ) THEN
    TotalNumberOfImpropers = TotalNumberOfImpropers + 1
    ImproperData(TotalNumberOfImpropers) = ImproperData2 ( i )
   ENDIF
  ENDDO
  CALL GetImproperType ()
 ENDSUBROUTINE GetImproperList

 SUBROUTINE GetBondType ()
  INTEGER :: i , j , i2(2) , j2(2)
  LOGICAL :: Reassign
  IF ( .NOT. ALLOCATED ( BondTypeData ) ) THEN
   ALLOCATE (BondTypeData ( TotalNumberOfBonds ))
   TotalNumberOfBondTypes = 0
  ENDIF
  DO i = 1 , TotalNumberOfBonds
   i2 ( 1 ) = AtomicData(BondData(i)%BondAtom(1))%AtomicTypeIndex
   i2 ( 2 ) = AtomicData(BondData(i)%BondAtom(2))%AtomicTypeIndex
   Reassign = .TRUE.
   DO j = 1 , TotalNumberOfBondTypes
    j2 = BondTypeData(j)%BondAtomType
    IF ( SameArray ( i2 , j2 ) .AND. &
    ( DABS (BondData(i)%BondDistance - BondTypeData(j)%BondDistance ) &
                             .LT. BondLengthTolerance  ) ) THEN
     BondData(i)%BondType = j
     Reassign = .FALSE.
     EXIT
    ENDIF
   ENDDO

   IF ( Reassign ) THEN
    TotalNumberOfBondTypes = TotalNumberOfBondTypes + 1
    BondTypeData(TotalNumberOfBondTypes)%BondAtomType = i2
    BondTypeData(TotalNumberOfBondTypes)%BondDistance = BondData(i)%BondDistance
    BondData(i)%BondType = TotalNumberOfBondTypes
   ENDIF
  ENDDO
 ENDSUBROUTINE GetBondType

 !get angle, dihedral and improper types 

 SUBROUTINE GetAngleType ()
  INTEGER :: i , j , i3(3) , j3(3)
  LOGICAL :: Reassign
  IF ( .NOT. ALLOCATED ( AngleTypeData ) ) THEN
   ALLOCATE (AngleTypeData ( TotalNumberOfAngles ))
   TotalNumberOfAngleTypes = 0
  ENDIF
  DO i = 1 , TotalNumberOfAngles
   i3 ( 1 ) = AtomicData(AngleData(i)%AngleAtom(1))%AtomicTypeIndex
   i3 ( 2 ) = AtomicData(AngleData(i)%AngleAtom(2))%AtomicTypeIndex
   i3 ( 3 ) = AtomicData(AngleData(i)%AngleAtom(3))%AtomicTypeIndex
   Reassign = .TRUE.
   DO j = 1 , TotalNumberOfAngleTypes
    j3 = AngleTypeData(j)%AngleAtomType
    IF ( SameArray ( i3 , j3 ) .AND. &
    ( DABS (AngleData(i)%AngleAngle - AngleTypeData(j)%AngleAngle ) &
                             .LT. AngleAngleTolerance  ) ) THEN
     AngleData(i)%AngleType = j
     Reassign = .FALSE.
     EXIT
    ENDIF
   ENDDO
   IF ( Reassign ) THEN
    TotalNumberOfAngleTypes = TotalNumberOfAngleTypes + 1
    AngleTypeData(TotalNumberOfAngleTypes)%AngleAtomType = i3
    AngleTypeData(TotalNumberOfAngleTypes)%AngleAngle = AngleData(i)%AngleAngle
    AngleData(i)%AngleType = TotalNumberOfAngleTypes
   ENDIF
  ENDDO
 ENDSUBROUTINE GetAngleType

 SUBROUTINE GetDihedralType ()
  INTEGER :: i , j , i4(4) , j4(4)
  LOGICAL :: Reassign
  IF ( .NOT. ALLOCATED ( DihedralTypeData ) ) THEN
   ALLOCATE (DihedralTypeData ( TotalNumberOfDihedrals ))
   TotalNumberOfDihedralTypes = 0
  ENDIF
  DO i = 1 , TotalNumberOfDihedrals
   i4 ( 1 ) = AtomicData(DihedralData(i)%DihedralAtom(1))%AtomicTypeIndex
   i4 ( 2 ) = AtomicData(DihedralData(i)%DihedralAtom(2))%AtomicTypeIndex
   i4 ( 3 ) = AtomicData(DihedralData(i)%DihedralAtom(3))%AtomicTypeIndex
   i4 ( 4 ) = AtomicData(DihedralData(i)%DihedralAtom(4))%AtomicTypeIndex
   Reassign = .TRUE.
   DO j = 1 , TotalNumberOfDihedralTypes
    j4 = DihedralTypeData(j)%DihedralAtomType
    IF ( SameArray ( i4 , j4 ) .AND. &
    ( DABS (DihedralData(i)%DihedralAngle - DihedralTypeData(j)%DihedralAngle ) &
                             .LT. DihedralAngleTolerance  ) ) THEN
     DihedralData(i)%DihedralType = j
     Reassign = .FALSE.
     EXIT
    ENDIF
   ENDDO
   IF ( Reassign ) THEN
    TotalNumberOfDihedralTypes = TotalNumberOfDihedralTypes + 1
    DihedralTypeData(TotalNumberOfDihedralTypes)%DihedralAtomType = i4
    DihedralTypeData(TotalNumberOfDihedralTypes)%DihedralAngle = DihedralData(i)%DihedralAngle
    DihedralData(i)%DihedralType = TotalNumberOfDihedralTypes
   ENDIF
  ENDDO
 ENDSUBROUTINE GetDihedralType

 SUBROUTINE GetImproperType ()
  INTEGER :: i , j , i4(4) , j4(4)
  LOGICAL :: Reassign
  IF ( .NOT. ALLOCATED ( ImproperTypeData ) ) THEN
   ALLOCATE (ImproperTypeData ( TotalNumberOfImpropers ))
   TotalNumberOfImproperTypes = 0
  ENDIF
  DO i = 1 , TotalNumberOfImpropers
   i4 ( 1 ) = AtomicData(ImproperData(i)%ImproperAtom(1))%AtomicTypeIndex
   i4 ( 2 ) = AtomicData(ImproperData(i)%ImproperAtom(2))%AtomicTypeIndex
   i4 ( 3 ) = AtomicData(ImproperData(i)%ImproperAtom(3))%AtomicTypeIndex
   i4 ( 4 ) = AtomicData(ImproperData(i)%ImproperAtom(4))%AtomicTypeIndex
   Reassign = .TRUE.
   DO j = 1 , TotalNumberOfImproperTypes
    j4 = ImproperTypeData(j)%ImproperAtomType
    IF ( SameImproper ( i4 , j4 ) .AND. &
    ( DABS (ImproperData(i)%ImproperAngle - ImproperTypeData(j)%ImproperAngle ) &
                             .LT. ImproperAngleTolerance  ) ) THEN
     ImproperData(i)%ImproperType = j
     Reassign = .FALSE.
     EXIT
    ENDIF
   ENDDO
   IF ( Reassign ) THEN
    TotalNumberOfImproperTypes = TotalNumberOfImproperTypes + 1
    ImproperTypeData(TotalNumberOfImproperTypes)%ImproperAtomType = i4
    ImproperTypeData(TotalNumberOfImproperTypes)%ImproperAngle = ImproperData(i)%ImproperAngle
    ImproperData(i)%ImproperType = TotalNumberOfImproperTypes
   ENDIF
  ENDDO
 ENDSUBROUTINE GetImproperType


 FUNCTION SameArray ( i , j ) RESULT ( res )
  IMPLICIT NONE
  INTEGER , INTENT (IN) :: i (:) , j (:)
  INTEGER , ALLOCATABLE :: k (:)
  INTEGER               :: m
  LOGICAL :: t1, t2, res
  res = .FALSE.
  ALLOCATE ( k ( SIZE ( i ) ) )
  k = (/ (j (m) , m = SIZE (j) , 1 , -1) /)
  t1 = ALL ( i .EQ. j )
  t2 = ALL ( i .EQ. k )
  res = t1 .OR. t2
 ENDFUNCTION SameArray 

 FUNCTION SameImproper ( i , j ) RESULT ( res )
  IMPLICIT NONE
  INTEGER , INTENT (IN) :: i (4) , j (4)
  INTEGER               :: m , n , o
  LOGICAL :: t1, res
  res = .FALSE.
  o = 0
  t1 = i (1) .EQ. j (1)
  IF ( t1 ) THEN
   DO m = 2 , 4
    DO n = 2 , 4
     IF (i ( m ) == j ( n )) o = o + 1
    ENDDO
   ENDDO
   IF ( o .GE. 3 ) res = .TRUE.
  ENDIF
 ENDFUNCTION SameImproper 

 SUBROUTINE NumberOfBADIForEachAtom ()
  INTEGER :: i , j
  AtomicData(:)%NBonds = 0
  AtomicData(:)%NAngles = 0
  AtomicData(:)%NDihedrals = 0
  AtomicData(:)%NImpropers = 0
  DO j = 1 , TotalNumberOfBonds
   DO i = 1 , TotalNumberOfAtoms
    IF ( ANY (BondData(j)%BondAtom .EQ. i ) ) AtomicData(i)%NBonds = AtomicData(i)%NBonds + 1 
   ENDDO
  ENDDO
  DO j = 1 , TotalNumberOfAngles
   DO i = 1 , TotalNumberOfAtoms
    IF ( ANY (AngleData(j)%AngleAtom .EQ. i ) ) AtomicData(i)%NAngles = AtomicData(i)%NAngles + 1 
   ENDDO
  ENDDO
  DO j = 1 , TotalNumberOfDihedrals
   DO i = 1 , TotalNumberOfAtoms
    IF ( ANY (DihedralData(j)%DihedralAtom .EQ. i ) ) AtomicData(i)%NDihedrals = AtomicData(i)%NDihedrals + 1 
   ENDDO
  ENDDO
  DO j = 1 , TotalNumberOfImpropers
   DO i = 1 , TotalNumberOfAtoms
    IF ( ANY (ImproperData(j)%ImproperAtom .EQ. i ) ) AtomicData(i)%NImpropers = AtomicData(i)%NImpropers + 1 
   ENDDO
  ENDDO
 ENDSUBROUTINE NumberOfBADIForEachAtom
ENDMODULE CalculateModule
