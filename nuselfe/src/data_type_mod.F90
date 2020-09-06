!-----------------------------------------------------------------------------------
!c    COPYRIGHT 2010: prof. Pavel Tkalich(TMSI), Xu Haihua(TMSI), Dao My Ha (TMSI)
!c    1. This subroutine is developed by Xu Haihua(TMSI), Dao My Ha(TMSI) and supervised
!c       by prof. Pavel Tkalich(TMSI)
!     2. TMSI: Tropical Marine Science Institute, Singapore
!c    3. For any enquiry pls contact tmsxh@nus.edu.sg
!-----------------------------------------------------------------------------------

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


!INTEGER(intKind)         :: info2;
END MODULE