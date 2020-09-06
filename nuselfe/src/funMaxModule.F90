!-----------------------------------------------------------------------------------
!c    COPYRIGHT 2010: prof. Pavel Tkalich(TMSI), Xu Haihua(TMSI), Dao My Ha (TMSI)
!c    1. This subroutine is developed by Xu Haihua(TMSI), Dao My Ha(TMSI) and supervised
!c       by prof. Pavel Tkalich(TMSI)
!     2. TMSI: Tropical Marine Science Institute, Singapore
!c    3. For any enquiry pls contact tmsxh@nus.edu.sg
!-----------------------------------------------------------------------------------

MODULE funMaxModule
!use mpi;
IMPLICIT NONE;


INTERFACE max
    MODULE PROCEDURE max2i 
END INTERFACE

CONTAINS

INTEGER FUNCTION max2i(x,y)
INTEGER ::x,y;


max2i= max(x,y);

END FUNCTION




END MODULE