MODULE PBCModule
!Module contains Common data to be shared by the main program and other modules
 USE GlobalDataModule
 IMPLICIT NONE

 CONTAINS 

 FUNCTION CrossProduct ( a , b) RESULT ( res )
  REAL ( KIND = 8 ) , DIMENSION (3) :: a , b
  REAL ( KIND = 8 )                 :: res (3)
  res ( 1 ) = a (2) * b (3) - a (3) * b (2)
  res ( 2 ) = a (3) * b (1) - a (1) * b (3)
  res ( 3 ) = a (1) * b (2) - a (2) * b (1)
 ENDFUNCTION CrossProduct

 SUBROUTINE CalculateBoxQuantities()
 !Calculates hmatrix and inverse of it. Modify to calculate the box extents
  IMPLICIT NONE
  REAL( KIND = 8) :: ralp,rbet, rgam
  REAL( KIND = 8) :: ca,sa,cb,sb,cg,sg,h23,h33,omg
  REAL( KIND = 8) :: ih13 , ih23 , ih33
  REAL( KIND = 8) :: lx , ly , lz , xy , xz , yz , xlo , ylo, zlo, xhi, yhi, zhi
  
  ralp = BoxParameters(4) * degtorad
  rbet = BoxParameters(5) * degtorad
  rgam = BoxParameters(6) * degtorad

  sa = DSIN(ralp)
  ca = DCOS(ralp)
  sb = DSIN(rbet)
  cb = DCOS(rbet)
  sg = DSIN(rgam)
  cg = DCOS(rgam)
  omg = (DSQRT( 1 - ca**2 - cb**2 - cg**2 + 2 * ca * cb * cg))

  h23 = (ca-cb*cg)/sg
  h33 = omg/sg

  HMatrix    = RESHAPE( (/ &
      BoxParameters(1)    , 0.0D0 , 0.0D0 ,&
      BoxParameters(2)*cg , BoxParameters(2)*sg  , 0.0D0 ,&
      BoxParameters(3)*cb , BoxParameters(3)*h23 , BoxParameters(3)*h33  &
      /) , (/3,3/) )   

  ih13 = (ca*cg-cb)/(omg*sg)
  ih23 = (cb*cg-ca)/(omg*sg)
  ih33 = sg/omg

  IHMatrix   = RESHAPE( (/ &
      1.0D0/BoxParameters(1)    ,  0.0D0                       , 0.0D0 ,                &
      -cg/(BoxParameters(1)*sg) , 1.0D0/(BoxParameters(2)*sg)  , 0.0D0 ,                &
      ih13/BoxParameters(1)      , ih23/BoxParameters(2)         ,ih33/BoxParameters(3)    &
      /) , (/3,3/) )   

  xlo = MINVAL ( AtomicData(:)%AtomicCoordinates(1) )
  ylo = MINVAL ( AtomicData(:)%AtomicCoordinates(2) )
  zlo = MINVAL ( AtomicData(:)%AtomicCoordinates(3) )
  xlo = 0.0D0
  ylo = 0.0D0
  zlo = 0.0D0
  lx  = BoxParameters(1)
  xy  = BoxParameters(2)*cg
  xz  = BoxParameters(3)*cb
  ly  = DSQRT( BoxParameters(2)**2 - xy**2)
  yz  = (BoxParameters(2)*BoxParameters(3)*ca - xy * xz)/ly
  lz  = DSQRT(BoxParameters(3)**2 - xz**2 - yz**2)
  xhi = xlo + lx
  yhi = ylo + ly
  zhi = zlo + lz
  LammpsBoxParameters = (/xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz/)
 ENDSUBROUTINE CalculateBoxQuantities

 FUNCTION PBCDistance ( i , j ) RESULT ( res )
  INTEGER , INTENT ( IN ) :: i , j
  REAL ( KIND = 8 ) :: res , sdisp ( 3 ) , disp ( 3 )
  sdisp = AtomicData(j)%AtomicScaledCoordinates - AtomicData(i)%AtomicScaledCoordinates &
          - ANINT ( AtomicData(j)%AtomicScaledCoordinates - AtomicData(i)%AtomicScaledCoordinates )
  disp = MATMUL ( HMatrix , sdisp )
  res = DSQRT ( DOT_PRODUCT ( disp , disp ) )
 ENDFUNCTION PBCDistance

 FUNCTION GetMinImage ( j , i ) RESULT ( res ) 
 ! get the coordinates of the minimum image of jth atom w.r.t ith atom
  INTEGER , INTENT ( IN ) :: i , j
  REAL ( KIND = 8 )       :: t (3) , res (3)
  t = AtomicData(j)%AtomicScaledCoordinates &
     - ANINT (  AtomicData(j)%AtomicScaledCoordinates &
              - AtomicData(i)%AtomicScaledCoordinates )
  res = MATMUL ( HMatrix , t )
 ENDFUNCTION GetMinImage

 !get rid of the two if blocks
 FUNCTION CalculateAngle ( a , b ) RESULT ( res )
 ! Calculate the angle between vectors a and b
  REAL ( KIND = 8 ) , DIMENSION (3) , INTENT (IN) :: a , b
  REAL ( KIND = 8 )                               :: res ,cosres
  cosres  = DOT_PRODUCT ( a , b ) &
         / DSQRT (DOT_PRODUCT ( a , a ) * DOT_PRODUCT ( b , b ) )
  IF ( cosres .GT.  1.0D0 ) cosres =  1.0D0
  IF ( cosres .LT. -1.0D0 ) cosres = -1.0D0
  res  = DACOS ( cosres ) * radtodeg
 ENDFUNCTION CalculateAngle

 ! min image w.r.t j or i?
 FUNCTION PBCAngle ( i , j , k ) RESULT (res)
 !calculate the angle between i , j and kth atoms with pbc
  INTEGER , INTENT ( IN ) :: i , j , k
  REAL ( KIND = 8 )       :: res
  REAL ( KIND = 8 ), DIMENSION (3)  :: a , b , c
  a = AtomicData(i)%AtomicCoordinates
  b = GetMinImage ( j , i )
  c = GetMinImage ( k , i )! min image w.r.t j or i?
  res = CalculateAngle ( a - b , c - b )
 ENDFUNCTION PBCAngle

 ! min image w.r.t j or i?
 FUNCTION PBCDihedralAngle ( i , j , k , l) RESULT (res)
 !calculate the dihedral angle formed by atoms i, j, k and l
  INTEGER , INTENT ( IN ) :: i , j , k , l
  REAL ( KIND = 8 )       :: res
  REAL ( KIND = 8 ), DIMENSION (3)  :: a , b , c , d
  REAL ( KIND = 8 ), DIMENSION (3)  :: ab , cb, bc, dc
  a = AtomicData(i)%AtomicCoordinates
  b = GetMinImage ( j , i )  ! 
  c = GetMinImage ( k , i )  ! min image w.r.t j or i?
  d = GetMinImage ( l , i )  !
  ab = a - b
  cb = c - b
  bc = b - c
  dc = d - c
  res = CalculateAngle ( CrossProduct (ab , cb ) , CrossProduct (bc , dc) )
 ENDFUNCTION PBCDihedralAngle

 FUNCTION PBCImproperAngle ( i , j , k , l) RESULT (res)
 !calculate the improper angle formed by atoms i, j, k and l (i: central atom)
  INTEGER , INTENT ( IN ) :: i , j , k , l
  REAL ( KIND = 8 )       :: res
  REAL ( KIND = 8 ), DIMENSION (3)  :: a , b , c , d
  REAL ( KIND = 8 ), DIMENSION (3)  :: ba , ca , da
  a = AtomicData(i)%AtomicCoordinates
  b = GetMinImage ( j , i )  ! 
  c = GetMinImage ( k , i )  ! min image w.r.t j or i?
  d = GetMinImage ( l , i )  !
  ba = b - a
  ca = c - a
  da = d - a
  res = CalculateAngle ( CrossProduct (ba , ca ) , CrossProduct (da , ca) )
 ENDFUNCTION PBCImproperAngle

ENDMODULE PBCModule
