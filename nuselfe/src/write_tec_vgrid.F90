!------------------------------------------------------------------------------------------------------------------
!write the vertical profile of the simulation data to tecplot format.
!------------------------------------------------------------------------------------------------------------------
MODULE SELFE_TEC_IO_vgrid
USE post_vars

CONTAINS


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :         write_TEC_vGRID 
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : write results to tecplot format the var is scalar or vector. 
!                       IF var is vector and have x and z component 
!                           var_name =' "var_x","var_z" '
!                   
!                  
!                 
!                 
!  Input        :  var_name: should container
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
SUBROUTINE write_TEC_vGRID(fileName,simtime,xver,zver,out_data,var_name,aux_data,append)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------

IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
CHARACTER(*)                   ,INTENT(IN   )::fileName;
REAL(rkind)                    ,INTENT(IN   )::simtime;
REAL(rkind)   ,DIMENSION(:  )  ,INTENT(IN   )::xver   ;
REAL(rkind)   ,DIMENSION(:,:)  ,INTENT(IN   )::zver   ;
REAL(rkind)   ,DIMENSION(:,:,:),INTENT(IN   )::out_data;
CHARACTER(*)                   ,INTENT(IN   )::var_name,aux_data;
INTEGER                        ,INTENT(IN   )::append;
!REAL(realKind),DIMENSION(:)  ,INTENT(IN   )::out_data;
!---------------------------------------------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------------------------------------------------------
INTEGER             ::i,j,UNIT,info2,temp,node_count;

UNIT = 100;
IF(append == 0 )THEN
    OPEN(UNIT=UNIT,FILE=fileName,ACTION='WRITE',IOSTAT=info2);
ELSEIF(append == 1 )THEN
    OPEN(UNIT=UNIT,FILE=fileName,ACTION='WRITE',POSITION='APPEND', IOSTAT=info2);
ENDIF

IF(append == 0 )THEN
!-----------------------------------------------------------------------
!IF append ==0 write the tecfile header.
!-----------------------------------------------------------------------
    WRITE(UNIT,"(A)")'TITLE = "SELFE OUTPUT" ';
    WRITE(UNIT,"(A,A,A)")'VARIABLES = "X", "Y", ',TRIM(var_name),',';
    WRITE(unit,"(A,I10,A,I10,A)")'ZONE T="T_1", DATAPACKING=POINT,I= ',n_vprof_node, ', J= ', nvrt,', ZONETYPE=Ordered,DATAPACKING=POINT,';
    WRITE(UNIT,"(A,F10.5)")'SolutionTime = ',simtime;      !write the output time
    WRITE(unit,"(A)")'DT=(SINGLE SINGLE SINGLE SINGLE)'
    WRITE(UNIT,"(A,A,A)")'AUXDATA Time = " ' ,TRIM(aux_data), ' " ';

    DO j = 1,nvrt
        DO i =1,n_vprof_node
            
            WRITE(UNIT,104)xver(i),zver(j,i),out_data(i,j,:);
        
        ENDDO
    ENDDO        
         

 
ENDIF
!-----------------------------------------------------------------------
!IF append ==1 write the results to the end of the file.
!-----------------------------------------------------------------------
IF(append == 1 )THEN

    WRITE(unit,"(A,I10,A,I10,A)")'ZONE T="T_1", DATAPACKING=POINT, I= ',n_vprof_node, ', J= ', nvrt,',ZONETYPE=Ordered,DATAPACKING=POINT';
    WRITE(UNIT,"(A,F10.5)")'SolutionTime = ',simtime;      !write the output time
    WRITE(unit,"(A)")'DT=(SINGLE SINGLE SINGLE SINGLE)'
    WRITE(UNIT,"(A,A,A)")'AUXDATA Time = " ' ,TRIM(aux_data), ' " ';
    
    DO j = 1,nvrt
        DO i =1,n_vprof_node
            
            WRITE(UNIT,104)xver(i),zver(j,i),out_data(i,j,:);
        
        ENDDO
    ENDDO           


ENDIF



CLOSE(UNIT);

104 FORMAT(2E20.8,20E15.6)
105 FORMAT(3I10);
106 FORMAT(20E15.6)


END SUBROUTINE




END MODULE;