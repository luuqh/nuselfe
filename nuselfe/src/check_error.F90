!-----------------------------------------------------------------------------------
!c    COPYRIGHT 2010: prof. Pavel Tkalich(TMSI), Xu Haihua(TMSI), Dao My Ha (TMSI)
!c    1. This subroutine is developed by Xu Haihua(TMSI), Dao My Ha(TMSI) and supervised
!c       by prof. Pavel Tkalich(TMSI)
!     2. TMSI: Tropical Marine Science Institute, Singapore
!c    3. For any enquiry pls contact tmsxh@nus.edu.sg
!-----------------------------------------------------------------------------------


MODULE check_error_mod
USE data_type_mod
!--------------------------------------------------
!windows
!--------------------------------------------------
#ifdef USE_WINDOWS
USE netcdf90
USE funMaxModule
#else
USE netcdf
#endif
!--------------------------------------------------
!Linux.
!--------------------------------------------------
!USE netcdf

IMPLICIT NONE;

INTEGER, PARAMETER :: NETCDF_ERROR_STOP_CODE  =2;
INTEGER, PARAMETER :: PROGRAM_ERROR_STOP_CODE =2;




CONTAINS

!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   check_error
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : check the error of calling the netcdf function. 
!                       if any error occur print out the error message and stop the program.
!                  
!                 
!                 
!  Input        : status (I) : the value returned by the netcdf function
!                  
!                 
!   
!                 
!
!  Input/output : 
!
!  Output       : 
!			      
!
!  Routines     : 
!                 
!
!  Remarks      :
!
!-----------------------------------------------------------------    
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------------------------------------------------

SUBROUTINE check_error(status)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind), INTENT (IN   ):: status
!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
IF( status /= nf90_noerr ) THEN

  WRITE(*,*) TRIM(nf90_strerror(status))
  STOP NETCDF_ERROR_STOP_CODE ;
ENDIF


END SUBROUTINE




!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :            check_results      
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : get the arry index in TMH grid by the relative poisiton
!                 
!  Input        : min      : min coord of TMH MODULE
!                 max      : max coord of TMH MODULE
!                 gridSize : grid size in x and y direction 
!                 point    : point coordint going to inquery
!                 array_index :arry index in x and y direction 
!                 info     : return value. if 0 mean point in the grid
!                                          if -1 mean point out of the grid.
!                 
!
!  Input/output : 
!
!  Output       : 
!			      
!
!  Routines     : 
!                 
!
!  Remarks      :
!
!-----------------------------------------------------------------    
!  References   :
!
!  Revisions    :
!------------------------------------------------------------------------------------------------------------------
SUBROUTINE check_results(status,caller,err_mesg)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)       , INTENT (IN   ):: status;
CHARACTER(*)      , INTENT (IN   ):: caller,err_mesg;

IF(status == sph_no_error) THEN

ELSEIF(status == sph_error_warning)THEN

    WRITE(*,"(4A)") "<WARNING> in subroutine ", TRIM(caller), "error msg:", TRIM(err_mesg);

ELSEIF(status == sph_error_error) THEN
    WRITE(*,"(4A)") "<ERROR> in subroutine ", TRIM(caller), "error msg:", TRIM(err_mesg);
    STOP PROGRAM_ERROR_STOP_CODE;
ENDIF




END SUBROUTINE

END MODULE
