MODULE SELFE_TEC_IO
USE elfe_glbl

CONTAINS


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :         write_TEC_hGRID 
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : write results to tecplot format the var is vector. 
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
SUBROUTINE write_TEC_hGRID(fileName,simtime,out_data,var_name,aux_data,append)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------

IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
CHARACTER(*)                   ,INTENT(IN   )::fileName;
REAL(rkind)                    ,INTENT(IN   )::simtime;
REAL(rkind)   ,DIMENSION(:,:)  ,INTENT(IN   )::out_data;
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

IF(info2/=0) THEN
    WRITE(*,*) "Failed to open file ",fileName;
ENDIF

IF(append == 0 )THEN
!-----------------------------------------------------------------------
!IF append ==0 write the tecfile header.
!-----------------------------------------------------------------------
    WRITE(UNIT,*)'TITLE = "SELFE OUTPUT" ';
    WRITE(UNIT,*)'VARIABLES = "X", "Y", ',TRIM(var_name),',';
    WRITE(UNIT,*)'ZONE T="T_1", DATAPACKING=POINT, NODES=',np,', ELEMENTS=',ne,', ZONETYPE=FETRIANGLE';
    WRITE(UNIT,*)'SolutionTime = ',simtime;      !write the output time
    WRITE(UNIT,*)'DT=(DOUBLE,DOUBLE,DOUBLE)';
    WRITE(UNIT,*)'AUXDATA Time = " ' ,TRIM(aux_data), ' " ';

        !---------write the node.
    DO i=1,np
        WRITE(UNIT,104)x(i),y(i),out_data(i,:);
    ENDDO

        !---------write element.
    DO i=1,ne
        WRITE(UNIT,105)nm(i,:);
    ENDDO
    
ENDIF
!-----------------------------------------------------------------------
!IF append ==1 write the results to the end of the file.
!-----------------------------------------------------------------------
IF(append == 1 )THEN

    WRITE(UNIT,*)'ZONE T="T_1", DATAPACKING=POINT, NODES=',np,', ELEMENTS=',ne,', ZONETYPE=FETRIANGLE';
    WRITE(UNIT,*)'SolutionTime = ',simtime;      !write the output time
    WRITE(UNIT,*)'VARSHARELIST = ([1, 2]=1), CONNECTIVITYSHAREZONE = 1';
    WRITE(UNIT,*)'DT=(DOUBLE,DOUBLE,DOUBLE)';
    WRITE(UNIT,*)'AUXDATA Time = " ' ,TRIM(aux_data), ' " '; 
    
    !---------as the mesh location didn't change, only out put the data-------
    DO i=1,np
        WRITE(UNIT,106)out_data(i,:);
    ENDDO       


ENDIF


CLOSE(UNIT);

104 FORMAT(2E16.8,20E15.6)
105 FORMAT(3I10);
106 FORMAT(20E15.6)
END SUBROUTINE

!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :         write_TEC_hGRID_scalar   
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : write results to tecplot format the var is scalar. 
!                   
!                  
!                 
!                 
!  Input        : 
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
SUBROUTINE write_TEC_hGRID_scalar(fileName,simtime,out_data,var_name,aux_data,append)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------

IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
CHARACTER(*)                   ,INTENT(IN   )::fileName;
REAL(rkind)                    ,INTENT(IN   )::simtime;
REAL(rkind)   ,DIMENSION(:)    ,INTENT(IN   )::out_data;
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
    WRITE(UNIT,*)'TITLE = "SELFE OUTPUT" ';
    WRITE(UNIT,*)'VARIABLES = "X", "Y"  " ',TRIM(var_name),' " ';

    WRITE(UNIT,*)'ZONE T="T_1", DATAPACKING=POINT, NODES=',np,', ELEMENTS=',ne,', ZONETYPE=FETRIANGLE';
    WRITE(UNIT,*)'SolutionTime = ',simtime, ",";      !write the output time
    WRITE(UNIT,*)'DT=(DOUBLE,DOUBLE,DOUBLE)';
    WRITE(UNIT,*)'AUXDATA Time = " ' ,TRIM(aux_data), ' " ';

    !---------write the node.
    DO i=1,np
        WRITE(UNIT,104)x(i),y(i),out_data(i);
    ENDDO

    !---------write element.
    DO i=1,ne
        WRITE(UNIT,105)nm(i,:);
    ENDDO
ENDIF
!-----------------------------------------------------------------------
!IF append ==1 write the results to the end of the file.
!-----------------------------------------------------------------------
IF(append == 1 )THEN
    WRITE(UNIT,*)'ZONE T="T_1", DATAPACKING=POINT, NODES=',np,', ELEMENTS=',ne,', ZONETYPE=FETRIANGLE';
    WRITE(UNIT,*)'SolutionTime = ',simtime;      !write the output time
    WRITE(UNIT,*)'VARSHARELIST = ([1, 2]=1), CONNECTIVITYSHAREZONE = 1';
    WRITE(UNIT,*)'DT=(DOUBLE,DOUBLE,DOUBLE)';
    WRITE(UNIT,*)'AUXDATA Time = " ' ,TRIM(aux_data), ' " '; 
    
    !---------as the mesh location didn't change, only out put the data-------
    DO i=1,np
        WRITE(UNIT,106)out_data(i);
    ENDDO       


ENDIF



CLOSE(UNIT);

104 FORMAT(2E25.10,E15.6)
105 FORMAT(3I10);
106 FORMAT(E15.6)
END SUBROUTINE


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :         write_TEC_hGRID_vector   
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : write results to tecplot format the var is vector. 
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
SUBROUTINE write_TEC_hGRID_vector(fileName,simtime,out_data,var_name,aux_data,append)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------

IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
CHARACTER(*)                   ,INTENT(IN   )::fileName;
REAL(rkind)                    ,INTENT(IN   )::simtime;
REAL(rkind)   ,DIMENSION(:,:)  ,INTENT(IN   )::out_data;
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
    WRITE(UNIT,*)'TITLE = "SELFE OUTPUT" ';
    WRITE(UNIT,*)'VARIABLES = "X", "Y"  " ',TRIM(var_name),' ", " ';

    WRITE(UNIT,*)'ZONE T="T_1", DATAPACKING=POINT, NODES=',np,', ELEMENTS=',ne,', ZONETYPE=FETRIANGLE';
    WRITE(UNIT,*)'SolutionTime = ',simtime;      !write the output time
    WRITE(UNIT,*)'DT=(DOUBLE,DOUBLE,DOUBLE)';
    WRITE(UNIT,*)'AUXDATA Time = " ' ,TRIM(aux_data), ' " ';

    !---------write the node.
    DO i=1,np
        WRITE(UNIT,104)x(i),y(i),out_data(i,1:);
    ENDDO

    !---------write element.
    DO i=1,ne
        WRITE(UNIT,105)nm(i,:);
    ENDDO
ENDIF
!-----------------------------------------------------------------------
!IF append ==1 write the results to the end of the file.
!-----------------------------------------------------------------------
IF(append == 1 )THEN
    WRITE(UNIT,*)'ZONE T="T_1", DATAPACKING=POINT, NODES=',np,', ELEMENTS=',ne,', ZONETYPE=FETRIANGLE';
    WRITE(UNIT,*)'SolutionTime = ',simtime;      !write the output time
    WRITE(UNIT,*)'VARSHARELIST = ([1, 2]=1), CONNECTIVITYSHAREZONE = 1';
    WRITE(UNIT,*)'DT=(DOUBLE,DOUBLE,DOUBLE)';
    WRITE(UNIT,*)'AUXDATA Time = " ' ,TRIM(aux_data), ' " '; 
    
    !---------as the mesh location didn't change, only out put the data-------
    DO i=1,np
        WRITE(UNIT,106)out_data(i,:);
    ENDDO       


ENDIF



CLOSE(UNIT);

104 FORMAT(2E25.15,E15.6)
105 FORMAT(3I10);
106 FORMAT(E15.6)
END SUBROUTINE

!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :         write_TEC_hGRID_bd   
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : write open and land boundary  to tecplot format 
!                   
!                  
!                 
!                 
!  Input        : 
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
SUBROUTINE write_TEC_hGRID_bd(fileName)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------

IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
CHARACTER(*)                   ,INTENT(IN   )::fileName;
!REAL(rkind)                    ,INTENT(IN   )::simtime;
!REAL(rkind)   ,DIMENSION(:)    ,INTENT(IN   )::out_data;
!CHARACTER(*)                   ,INTENT(IN   )::var_name,aux_data;
!INTEGER                        ,INTENT(IN   )::append;
!REAL(realKind),DIMENSION(:)  ,INTENT(IN   )::out_data;
!---------------------------------------------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------------------------------------------------------
INTEGER                     ::i,j,UNIT,info2;
CHARACTER(LEN=50)           ::zoneName;
UNIT = 100;
!IF(append == 0 )THEN
!    OPEN(UNIT=UNIT,FILE=fileName,ACTION='WRITE',IOSTAT=info2);
!ELSEIF(append == 1 )THEN
!    OPEN(UNIT=UNIT,FILE=fileName,ACTION='WRITE',POSITION='APPEND', IOSTAT=info2);
!ENDIF
OPEN(UNIT=UNIT,FILE=fileName,ACTION='WRITE',IOSTAT=info2);
!-----------------------------------------------------------------------
!write the tecfile header.
!-----------------------------------------------------------------------
WRITE(UNIT,*)'TITLE = "SELFE OUTPUT" ';
WRITE(UNIT,*)'VARIABLES = "X", "Y"  "dd " ';

!-----------------------------------------------------------------------
!Write the open boundary.....
!-----------------------------------------------------------------------
zoneName ="openBd_";
DO i = 1,nope
    WRITE(zoneName(8:11),100)i;
    WRITE(UNIT,101)'ZONE T= " ', TRIM(zoneName), ' ",I= ',nond(i),',F=POINT';
    DO j = 1,nond(i)
        WRITE(UNIT,102)x(iond(i,j)),y(iond(i,j)),1.0d0;
    ENDDO
ENDDO

!-----------------------------------------------------------------------
!Write the land boundary.....
!-----------------------------------------------------------------------
zoneName ="landBd_";
DO i = 1,nope
    WRITE(zoneName(8:11),100)i;
    WRITE(UNIT,101)' ZONE T= " ', TRIM(zoneName), ' ",I= ',nlnd(i),',F=POINT';
    DO j = 1,nlnd(i)
        WRITE(UNIT,102)x(ilnd(i,j)),y(ilnd(i,j)),1.0d0;
    ENDDO
ENDDO

CLOSE(UNIT);

100  FORMAT(I3);
101  FORMAT(A,A,A,I6,A);
102  FORMAT(2E25.15,E15.6);

END SUBROUTINE


END MODULE;