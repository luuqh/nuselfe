MODULE data_type_mod
IMPLICIT NONE;

INTEGER, PARAMETER :: realKind =8;
INTEGER, PARAMETER :: intKind =4;
INTEGER, PARAMETER :: CharLen = 200;
!INTEGER, PARAMETER :: INTKIND = intKind;
!INTEGER, PARAMETER :: realKind= realKind;

CHARACTER(CharLen) ::cBuff;

INTEGER(intKind),PARAMETER ::sph_no_error =0;
INTEGER(intKind),PARAMETER ::sph_error_warning   =-2;
INTEGER(intKind),PARAMETER ::sph_error_error     =-1;
!REAL(realKind)  ,PARAMETER ::pi      = 3.14159265358979323846d0;
END MODULE