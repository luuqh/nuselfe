!---------------------------------------------------
!write the results to sms data format.
!---------------------------------------------------
MODULE write_sms_grid_mod
USE post_vars
IMPLICIT NONE

INTEGER,PRIVATE             :: iVersion,ID100,iObjectType,ID110, iSFLT,ID120,iSFLG,IDDim,iISTAT;
INTEGER,PRIVATE             :: ID150,iVECTYPE,ID160,iOBJID,ID170,iNUMDATA,ID180,iNUMCELLS,ID190;
INTEGER,PRIVATE             :: ID200,ID210;
character(LEN=40),PRIVATE   :: iNAME;
INTEGER,PRIVATE,PARAMETER   :: sms_byte =4;

CONTAINS


!--------------------------------------------------
!open the file and return the file unit and record location.
!--------------------------------------------------
SUBROUTINE write_sms_head(fileName,unit,ilrec)
IMPLICIT NONE
!---------------------------------------------------
!ARGUMENT....................
!---------------------------------------------------
CHARACTER(LEN=*),INTENT(INOUT):: fileName;
INTEGER         ,INTENT(INOUT):: unit;
INTEGER         ,INTENT(INOUT):: ilrec;
!---------------------------------------------------
!local varialbe
!---------------------------------------------------
INTEGER             :: m,info2;

!nbyte     = sms_byte;
!--------------------------------------------------
!init the variables.
!--------------------------------------------------

iVersion   =3000
iObjectType=3               !3-mesh, 5-scatter points
iSFLT      =sms_byte
iSFLG      =sms_byte
iVECTYPE   =0
iOBJID     =0
iNUMDATA   =np
iNUMCELLS  =ne
iNAME(1:)  =file63(1:lfile63-3);
iISTAT     =0
ID100      =100
ID110      =110
ID120      =120
IDDim      =120+ivs*10      !-- 130-scalar, 140-vector --
ID150      =150
ID160      =160
ID170      =170
ID180      =180
ID190      =190
ID200      =200
ID210      =210

!--------------------------
OPEN(unit,file=TRIM(fileName),status='replace',access='direct',recl=onbyte,iostat=info2);
IF(info2/=0) WRITE(*,*) "Failed to open file :", TRIM(fileName);

ilrec  = 0;

write(unit,rec=ilrec+1)  iVersion
write(unit,rec=ilrec+2)  ID100
write(unit,rec=ilrec+3)  iObjectType
write(unit,rec=ilrec+4)  ID110
write(unit,rec=ilrec+5)  iSFLT
write(unit,rec=ilrec+6)  ID120
write(unit,rec=ilrec+7)  iSFLG
write(unit,rec=ilrec+8)  IDDim
write(unit,rec=ilrec+9)  ID150
write(unit,rec=ilrec+10) iVECTYPE
write(unit,rec=ilrec+11) ID160
write(unit,rec=ilrec+12) iOBJID
write(unit,rec=ilrec+13) ID170
write(unit,rec=ilrec+14) iNUMDATA
write(unit,rec=ilrec+15) ID180
write(unit,rec=ilrec+16) iNUMCELLS
write(unit,rec=ilrec+17) ID190
ilrec=ilrec+17
do m=1,10
 write(unit,rec=ilrec+m) iNAME(onbyte*(m-1)+1:onbyte*m)
enddo
ilrec=ilrec+10

END SUBROUTINE


!--------------------------------------------------
!open the file and return the file unit and record location.
!--------------------------------------------------
SUBROUTINE write_sms_time_data(fileName,unit,ilrec,out_data)
IMPLICIT NONE
!---------------------------------------------------
!ARGUMENT....................
!---------------------------------------------------
CHARACTER(LEN=*),INTENT(INOUT) :: fileName;
INTEGER         ,INTENT(INOUT) :: unit;
INTEGER         ,INTENT(INOUT) :: ilrec;
REAL(rkind)     ,DIMENSION(:,:):: out_data; 

!---------------------------------------------------
!Local variable.
!---------------------------------------------------
INTEGER                       :: i,m,ils(np);

WRITE(unit,rec=ilrec+1) ID200
WRITE(unit,rec=ilrec+2) iISTAT
IF(onbyte == 1) THEN
    WRITE(unit,rec=ilrec+3) REAL(time,sms_byte);
ENDIF

ilrec=ilrec+3

DO i=1,np
    ils(i)=(i-1)*ivs
ENDDO	    

DO i=1, np
    DO m=1,ivs 
     
        IF(onbyte == 1) THEN
            WRITE(unit,rec=ilrec+ils(i)+m) REAL(out_data(i,m),sms_byte);
        ENDIF          
!        WRITE(unit,rec=ilrec+ils(i)+m) tec_out(i,m)
    ENDDO !m
ENDDO

!update the record length.

ilrec=ilrec+np*ivs

END SUBROUTINE


SUBROUTINE close_sms_file(fileName,unit,ilrec)

IMPLICIT NONE
!---------------------------------------------------
!ARGUMENT....................
!---------------------------------------------------
CHARACTER(LEN=*),INTENT(INOUT):: fileName;
INTEGER         ,INTENT(INOUT):: unit;
INTEGER         ,INTENT(INOUT):: ilrec;


!--------write the data finish flag.
WRITE(unit,rec=ilrec+1) ID210

!--------close the file--------
CLOSE(unit);



END SUBROUTINE

END MODULE



