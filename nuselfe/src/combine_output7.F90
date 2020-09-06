MODULE combine_output7
USE post_vars
IMPLICIT NONE

INTEGER,DIMENSION(:),ALLOCATABLE        :: ihot_len;

REAL(orkind),DIMENSION(:,:),ALLOCATABLE :: temp_outb(:,:,:);
REAL(orkind),DIMENSION(:  ),ALLOCATABLE :: temp_eta2;
REAL(orkind)                            :: temp_time;
INTEGER                                 :: temp_it;




CONTAINS


SUBROUTINE combine_results()

IMPLICIT NONE;

CHARACTER(LEN=200)  :: fileName,cday,crank;
CHARACTER(LEN=200)  :: outFileName;
INTEGER             :: iday,irank,ispool;
INTEGER             :: com_unit = 201;


CALL get_lg_map()

!----------------------------------------------------------
!READ data for each day step.
!----------------------------------------------------------
DO iday = ibgn,iend

    !compose the file Name
    WRITE(cday,'(i12)')iday;
    cday  = TRIM(ADJUSTL(cday))//'_';

    fileName = TRIM(inFolder)//TRIM(cday);
    
    outFileName = TRIM(outFolder)//TRIM(cday)//TRIM(file63);
    CALL out_combine_head(outFileName,com_unit);
    
    DO ispool = 1,nrec
        CALL get_time_data(fileName,ispool);
        CALL out_combine_time_data(com_unit)
        WRITE(*,*) "time = ", time,"it =",it;
        
    ENDDO
  
    CLOSE(com_unit);  
ENDDO




END SUBROUTINE


!--------------------------------------------------------
!get the input parameters 
!--------------------------------------------------------
SUBROUTINE  get_pos_parm()

!--------------------------------------------------------
!Local variables.
!--------------------------------------------------------

!--------------------------------------------------------
!--------------------------------------------------------


open(10,file='combine_output.in',status='old')
read(10,'(a30)') file63   !file name to be combined.
read(10,*) ibgn,iend      !begin and end day
read(10,*) inFolder;      !input data folder.
read(10,*) outFolder;     !output data foler.
close(10)


!inFolder = TRIM(inFolder)//"/";
!outFolder =TRIM(outFolder)//"/";
WRITE(*,*)"inFolder:",TRIM(inFolder);
WRITE(*,*)"outFolder:",TRIM(outFolder);

!-------------------------------------------------
!get the combine varable type, scalar to vector.
!-------------------------------------------------
lfile63=len_trim(file63)
if(file63((lfile63-1):lfile63).eq.'61'.or.file63((lfile63-1):lfile63).eq.'63') then
    ivs=1
else if(file63((lfile63-1):lfile63).eq.'62'.or.file63((lfile63-1):lfile63).eq.'64') then
    ivs=2
else
    print*, 'Unknown file type:',file63
endif

if(file63((lfile63-1):lfile63).eq.'61'.or.file63((lfile63-1):lfile63).eq.'62') then
    i23d=2
else
    i23d=3
endif




END SUBROUTINE




!--------------------------------------------------------
!get local to global mapping 
!--------------------------------------------------------
SUBROUTINE  get_lg_map()
!--------------------------------------------------------
!Local variables.
!--------------------------------------------------------
INTEGER            :: irank,istat,itmp,i,j,k,ntmp;
INTEGER            :: m,mm,ipgb,iegb;
CHARACTER(LEN=200) :: fileName;
INTEGER            :: ne_local_max,np_local_max;

! Read local_to_global_0000 for global info
fileName = TRIM(inFolder)//'local_to_global_0000'
open(10,file=TRIM(fileName),status='old')
    read(10,*)ne_global,np_global,nvrt,nproc 
close(10)


!-----------------------------------------------------------
! Read in local-global mappings from all ranks
!-----------------------------------------------------------

ne_local_max = 0;
np_local_max = 0;

fdb='local_to_global_0000'
lfdb=len_trim(fdb)

DO  irank=0,nproc-1
  
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=TRIM(inFolder)//fdb,status='old')

    read(10,*) !global info

    read(10,*) !info
    
    read(10,*)itmp  !ne_local;
    
    ne_local_max = max(ne_local_max,itmp);
    
    do i=1,itmp
      read(10,*)
    enddo !i
    
    read(10,*)itmp  !np_local;
    
    np_local_max = max(np_local_max,itmp);

    CLOSE(10);
    
ENDDO

!--------------------------------------------------------
!allocate the memory
!--------------------------------------------------------
IF (ALLOCATED(ne_local))    DEALLOCATE(ne_local)
IF (ALLOCATED(ielg_local))  DEALLOCATE(ielg_local)
ALLOCATE(ne_local(0:nproc-1),ielg_local(0:nproc-1,ne_local_max),stat=istat)
if(istat/=0) stop 'Allocation error: ne,ielg';


IF (ALLOCATED(np_local))    DEALLOCATE(np_local)
IF (ALLOCATED(iplg_local))  DEALLOCATE(iplg_local)
ALLOCATE(np_local(0:nproc-1),iplg_local(0:nproc-1,np_local_max),stat=istat)
if(istat/=0) stop 'Allocation error: np,iplg';


IF (ALLOCATED(kbp00_local))    DEALLOCATE(kbp00_local)
ALLOCATE(kbp00_local(0:nproc-1,np_local_max),stat=istat)
if(istat/=0) stop 'Allocation error: kbp00_local';

IF (ALLOCATED(ztot))     DEALLOCATE(ztot)
IF (ALLOCATED(sigma))    DEALLOCATE(sigma)
ALLOCATE(ztot(nvrt),sigma(nvrt),stat=istat)
if(istat/=0) stop 'Allocation error: ztot,sigma';

IF (ALLOCATED(nm_local))     DEALLOCATE(nm_local)
ALLOCATE(nm_local(ne_local_max,3),stat=istat)
if(istat/=0) stop 'Allocation error: nm_local';
!--------------------------------------------------------
!allocate the memory
!--------------------------------------------------------
IF (ALLOCATED(x))     DEALLOCATE(x)
IF (ALLOCATED(y))     DEALLOCATE(y)
IF (ALLOCATED(dp))    DEALLOCATE(dp)
IF (ALLOCATED(kbp)) DEALLOCATE(kbp)
ALLOCATE(x(np_global),y(np_global),dp(np_global),kbp(np_global),stat=istat)
if(istat/=0) stop 'Allocation error: x,y,dp,kbp';

IF (ALLOCATED(nm)) DEALLOCATE(nm)
ALLOCATE(nm(ne_global,3),stat=istat)
if(istat/=0) stop 'Allocation error: nm';
!-----------------------------------------------------------
! Read in local-global mappings from all ranks
!-----------------------------------------------------------

  fdb='local_to_global_0000'
  lfdb=len_trim(fdb)

  do irank=0,nproc-1
  
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=TRIM(inFolder)//fdb,status='old')

    read(10,*) !global info

    read(10,*) !info
    
    read(10,*)ne_local(irank)
    do i=1,ne_local(irank)
      read(10,*)j,ielg_local(irank,i)
    enddo !i
    
    read(10,*)np_local(irank)
    do i=1,np_local(irank)
      read(10,*)j,iplg_local(irank,i)
    enddo
    
    read(10,*)itmp !sides
    do i=1,itmp
      read(10,*)
    enddo
!-----------------------------------------------------------
!    print*, 'Mapping:',irank,ielg(irank,ne(irank))
!-----------------------------------------------------------
    read(10,*) !'Header:'
    read(10,'(a)')data_format,version,start_time
    read(10,*)nrec,dtout,nspool,nvrt,kz,h0,h_s,h_c,theta_b,theta_f
    read(10,*)(ztot(k),k=1,kz-1),(sigma(k),k=1,nvrt-kz+1)
    read(10,*)np_local(irank),ne_local(irank),  &
              (x(iplg_local(irank,m)),y(iplg_local(irank,m)),dp(iplg_local(irank,m)),kbp00_local(irank,m),m=1,np_local(irank)), &
    &         (ntmp,(nm_local(m,mm),mm=1,3),m=1,ne_local(irank))

    close(10)

!-----------------------------------------------------------
!   Compute kbp (to avoid mismatch of indices)
!-----------------------------------------------------------
    do m=1,np_local(irank)
      ipgb=iplg_local(irank,m)
      kbp(ipgb)=kbp00_local(irank,m)
    enddo !m
 
!-----------------------------------------------------------
!   Reconstruct connectivity table-------------------
!-----------------------------------------------------------
    do m=1,ne_local(irank)
    
      iegb=ielg_local(irank,m)
      if(iegb>ne_global) stop 'Overflow!'
    
      do mm=1,3
    
        itmp=nm_local(m,mm)
    
!        if(itmp>np(irank).or.itmp<=0) then
!          write(*,*)'Overflow:',m,mm,itmp
!          stop
!        endif
        
        nm(iegb,mm)=iplg_local(irank,itmp)
      enddo !mm
    enddo !m
enddo !irank=0,nproc-1



!--------------------------------------------------------
!--------------------------------------------------------
IF (ALLOCATED(ihot_len)) DEALLOCATE(ihot_len)
ALLOCATE(ihot_len(0:nproc-1),stat=istat)
if(istat/=0) stop 'Allocation error: ihot_len';

! Compute record length for each rank-specific binary output per time step
do irank=0,nproc-1

!    ihot_len(irank)=nbyte*(2+np_local(irank))
    
!    if(i23d==2) then
!      ihot_len(irank)=ihot_len(irank)+nbyte*ivs*np_local(irank)
!      
!    else
!      do i=1,np_local(irank)
!            do k=max0(1,kbp00_local(irank,i)),nvrt
!              do m=1,ivs
!                ihot_len(irank)=ihot_len(irank)+nbyte;
!              enddo !m
!            enddo !k
!      enddo !i
!    endif

    ihot_len(irank)=onbyte*(2+np_local(irank))
    
    if(i23d==2) then
      ihot_len(irank)=ihot_len(irank)+onbyte*ivs*np_local(irank)
      
    else
      do i=1,np_local(irank)
            do k=max0(1,kbp00_local(irank,i)),nvrt
              do m=1,ivs
                ihot_len(irank)=ihot_len(irank)+onbyte;
              enddo !m
            enddo !k
      enddo !i
    endif

enddo !irank


!--------------------------------------------------------
!allocate temp memory.
!--------------------------------------------------------
!IF (ALLOCATED(temp_outb)) DEALLOCATE(temp_outb)
!ALLOCATE(temp_outb(np_local_max,nvrt,ivs),stat=istat)
!if(istat/=0) stop 'Allocation error: temp_outb';

IF (ALLOCATED(eta2)) DEALLOCATE(eta2)
IF (ALLOCATED(outb)) DEALLOCATE(outb)
ALLOCATE(eta2(np_global),outb(np_global,nvrt,ivs),stat=istat)
if(istat/=0) stop 'Allocation error: eta2 , outb';

IF (ALLOCATED(temp_eta2)) DEALLOCATE(temp_eta2)
IF (ALLOCATED(temp_outb)) DEALLOCATE(temp_outb);
ALLOCATE(temp_eta2(np_global),temp_outb(np_global,nvrt,ivs),stat=istat)
if(istat/=0) stop 'Allocation error: eta2 , outb';

END SUBROUTINE


!--------------------------------------------------------
!read SELFE output data (SELFE ORIGINAL OUTPUT) at particular day and time
!
!--------------------------------------------------------
SUBROUTINE get_time_data(fileName_prex,ispool)
IMPLICIT NONE
!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
CHARACTER(LEN=*)           :: fileName_prex;
INTEGER                    :: ispool;
!--------------------------------------------------------
!Local variables.
!--------------------------------------------------------
INTEGER            :: irank,lfgb,i,m,k;
CHARACTER(LEN=200)  :: fgb2;
!--------------------------------------------------------
!--------------------------------------------------------

do irank=0,nproc-1
    !Open input file
    fgb2 = fileName_prex
    lfgb = len_trim(fgb2);
    write(fgb2(lfgb+1:),'(i4.4)') irank

    open(63,file=TRIM(fgb2)//'_'//file63,access='direct',recl=ihot_len(irank),status='old')

    if(i23d==2) then
!      read(63,rec=ispool)time,it,(eta2(iplg(irank,i)),i=1,np_local(irank)),&
!                &        ((temp_outb(iplg(irank,i),1,m),m=1,ivs),i=1,np_local(irank));
      read(63,rec=ispool)temp_time,temp_it,(temp_eta2(iplg_local(irank,i)),i=1,np_local(irank)),&
                &        ((temp_outb(iplg_local(irank,i),1,m),m=1,ivs),i=1,np_local(irank));

    else !3D
!      read(63,rec=ispool)time,it,(eta2(iplg(irank,i)),i=1,np(irank)), &
!               &    (((outb(iplg(irank,i),k,m),m=1,ivs),k=max0(1,kbp01(irank,i)),nvrt),i=1,np(irank))
      read(63,rec=ispool)temp_time,temp_it,(temp_eta2(iplg_local(irank,i)),i=1,np_local(irank)), &
               &    (((temp_outb(iplg_local(irank,i),k,m),m=1,ivs),k=max0(1,kbp00_local(irank,i)),nvrt),i=1,np_local(irank))


    endif

    ! Close input file
    close(63)

enddo !irank

time  = REAL(temp_time,rKind);;
it    = REAL(temp_it,rKind);
eta2  = REAL(temp_eta2,rKind);
outb  = REAL(temp_outb,rKind);



END SUBROUTINE 



!--------------------------------------------------------
!subroutine output combine head
!--------------------------------------------------------
SUBROUTINE out_combine_head(fileName,UNIT)
!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
CHARACTER(LEN=*)                 :: fileName;
INTEGER           ,INTENT(INOUT) :: UNIT
!--------------------------------------------------------
!Local variables. m,mm;
!--------------------------------------------------------
INTEGER                    :: k,kin,m,mm,stat;
!--------------------------------------------------------
!out put the head of the combined file
!--------------------------------------------------------
open(UNIT,file=TRIM(fileName),status='replace',FORM='BINARY',iostat=stat);
IF(stat/=0) WRITE(*,*)"Fail to open file: ", TRIM(fileName);

data_format ='DataFormat v5.0'
variable_nm =file63    !not important
variable_dim=file63

write(Unit) data_format  ! data format.
write(Unit) version      !
write(Unit) start_time
write(Unit) variable_nm
write(Unit) variable_dim

!------------------------------------------ 
!write the record ........
!------------------------------------------
write(Unit) nrec
write(Unit) dtout
write(Unit) nspool
write(Unit) ivs
write(Unit) i23d

!------------------------------------------ 
!Vertical grid
!------------------------------------------ 
write(Unit) nvrt
write(Unit) kz
write(Unit) h0
write(Unit) h_s
write(Unit) h_c
write(Unit) theta_b
write(Unit) theta_f

!---------------Z layer----------------
do k=1,kz-1
  write(Unit) ztot(k)
enddo

!---------------S layer------------
do k=kz,nvrt
  kin=k-kz+1
  write(Unit) sigma(kin)
enddo !k 

!------------------------------------------ 
!    !Horizontal grid
!------------------------------------------ 
write(Unit) np_global
write(Unit) ne_global

do m=1,np_global
  write(Unit) x(m)
  write(Unit) y(m)
  write(Unit) dp(m)
  write(Unit) kbp(m)
enddo !m=1,np
    
do m=1,ne_global
  write(Unit) 3
  do mm=1,3
    write(Unit) nm(m,mm)
  enddo !mm
enddo !m

!------------------------------------------ 
!finsihed close the file
!------------------------------------------ 
!CLOSE(65)
END SUBROUTINE


!--------------------------------------------------------
!subroutine output combine head
!--------------------------------------------------------
SUBROUTINE out_combine_time_data(UNIT)
!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
!CHARACTER(LEN=*)           :: fileName;
INTEGER          ,INTENT(IN):: UNIT;
!--------------------------------------------------------
!Local variables. m,mm;
!--------------------------------------------------------
INTEGER                    :: i,k,kin,m,mm;
!--------------------------------------------------------
!out put the head of the combined file
!--------------------------------------------------------
!open(65,file=TRIM(fileName),ACTION="WRITE",POSITION ='APPEND',FORM='BINARY');


write(Unit) time ;
write(Unit) it   ;

do i=1,np_global
    write(Unit) eta2(i)
enddo    !i

do i=1,np_global
    if(i23d==2) then
      do m=1,ivs
        write(Unit)   outb(i,1,m)
      enddo !m
    else !i23d=3 
      do k=max0(1,kbp(i)),nvrt
        do m=1,ivs
          write(Unit) outb(i,k,m)
        enddo !m
      enddo !k
    endif !i23d
enddo !i

!------------------------------------------ 
!finsihed close the file
!------------------------------------------ 
!CLOSE(65)

END SUBROUTINE


!--------------------------------------------------------
!Read the head file of the combined results
!    OPEN THE FILE AND READ THE HEAD DATA, the Unit will kept.
!--------------------------------------------------------
SUBROUTINE in_combine_head(fileName,Unit)
!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
CHARACTER(LEN=*)              :: fileName;
INTEGER        ,INTENT(INOUT) :: Unit;
!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
!CHARACTER(LEN=*)           :: fileName;
!--------------------------------------------------------
!Local variables. m,mm;
!--------------------------------------------------------
INTEGER                    :: k,kin,m,mm,temp_int,istat;
!--------------------------------------------------------
!out put the head of the combined file
!--------------------------------------------------------

OPEN(Unit,file=TRIM(fileName),ACTION="READ",status='OLD',FORM='BINARY',IOSTAT=istat );
IF(istat/=0)WRITE(*,*)"Failed to open file: ", TRIM(fileName);

READ(Unit) data_format  ! data format.
READ(Unit) version      !
READ(Unit) start_time
READ(Unit) variable_nm
READ(Unit) variable_dim

!------------------------------------------ 
!write the record ........
!------------------------------------------
READ(Unit) nrec
READ(Unit) dtout
READ(Unit) nspool
READ(Unit) ivs
READ(Unit) i23d

!------------------------------------------ 
!Vertical grid
!------------------------------------------ 
READ(Unit) nvrt
READ(Unit) kz
READ(Unit) h0
READ(Unit) h_s
READ(Unit) h_c
READ(Unit) theta_b
READ(Unit) theta_f

!-------------------------------------------------
!ALLOCATE memory
!-------------------------------------------------
IF(ALLOCATED(ztot))  DEALLOCATE(ztot);
IF(ALLOCATED(sigma)) DEALLOCATE(sigma);
ALLOCATE(ztot(nvrt),sigma(nvrt),stat=istat);
if(istat/=0) stop 'Allocation error: ztot,sigma';

do k=1,kz-1
  READ(Unit) ztot(k)
enddo

do k=kz,nvrt
  kin=k-kz+1
  READ(Unit) sigma(kin)
enddo !k 

!------------------------------------------ 
!    !Horizontal grid
!------------------------------------------ 
READ(Unit) np_global
READ(Unit) ne_global

np = np_global;
ne = ne_global;

!-------------------------------------------------
!ALLOCATE memory
!-------------------------------------------------
IF(ALLOCATED(x))     DEALLOCATE(x);
IF(ALLOCATED(y))     DEALLOCATE(y);
IF(ALLOCATED(dp))    DEALLOCATE(dp);
IF(ALLOCATED(kbp))   DEALLOCATE(kbp);

ALLOCATE(x(np_global),y(np_global),dp(np_global),kbp(np_global),stat=istat);
if(istat/=0) stop 'Allocation error: x,y,dp,kbp';


IF(ALLOCATED(nm)) DEALLOCATE(nm); 
ALLOCATE(nm(ne_global,3),stat=istat);
if(istat/=0) stop 'Allocation error: mn';


do m=1,np_global
  READ(Unit) x(m)
  READ(Unit) y(m)
  READ(Unit) dp(m)
  READ(Unit) kbp(m)
enddo !m=1,np
    
do m=1,ne_global
  READ(Unit) temp_int
  do mm=1,3
    READ(Unit) nm(m,mm)
  enddo !mm
enddo !m

!------------------------------------------ 
!finsihed close the file
!------------------------------------------ 
!CLOSE(65)
!---------------------------------------------------------
!ALLOCATE MEMORY eta2
!---------------------------------------------------------
IF(ALLOCATED(eta2)) DEALLOCATE(eta2);
IF(ALLOCATED(outb)) DEALLOCATE(outb);

IF(i23d == 2)THEN
    ALLOCATE(eta2(np_global),outb(np_global,1,ivs),stat=istat);
ELSE
    ALLOCATE(eta2(np_global),outb(np_global,nvrt,ivs),stat=istat);
ENDIF

if(istat/=0) stop 'Allocation error:eta2,outb';



END SUBROUTINE



!--------------------------------------------------------
!Read the head file of the combined results
!    OPEN THE FILE AND READ THE HEAD DATA, the Unit will kept.
!--------------------------------------------------------
SUBROUTINE in_combine_head2(fileName,Unit)
!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
CHARACTER(LEN=*)              :: fileName;
INTEGER        ,INTENT(INOUT) :: Unit;
!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
!CHARACTER(LEN=*)           :: fileName;
!--------------------------------------------------------
!Local variables. m,mm;
!--------------------------------------------------------
INTEGER                    :: k,kin,m,mm,temp_int,istat;
!--------------------------------------------------------
!out put the head of the combined file
!--------------------------------------------------------

OPEN(Unit,file=TRIM(fileName),ACTION="READ",status='OLD',FORM='BINARY',IOSTAT=istat );
IF(istat/=0)WRITE(*,*)"Failed to open file: ", TRIM(fileName);

READ(Unit) data_format  ! data format.
READ(Unit) version      !
READ(Unit) start_time
READ(Unit) variable_nm
READ(Unit) variable_dim

!------------------------------------------ 
!write the record ........
!------------------------------------------
READ(Unit) nrec
READ(Unit) dtout
READ(Unit) nspool
READ(Unit) ivs
READ(Unit) i23d

!------------------------------------------ 
!Vertical grid
!------------------------------------------ 
READ(Unit) nvrt
READ(Unit) kz
READ(Unit) h0
READ(Unit) h_s
READ(Unit) h_c
READ(Unit) theta_b
READ(Unit) theta_f

!-------------------------------------------------
!ALLOCATE memory
!-------------------------------------------------
!IF(ALLOCATED(ztot))  DEALLOCATE(ztot);
!IF(ALLOCATED(sigma)) DEALLOCATE(sigma);
!ALLOCATE(ztot(nvrt),sigma(nvrt),stat=istat);
!if(istat/=0) stop 'Allocation error: ztot,sigma';

do k=1,kz-1
  READ(Unit) ztot(k)
enddo

do k=kz,nvrt
  kin=k-kz+1
  READ(Unit) sigma(kin)
enddo !k 

!------------------------------------------ 
!    !Horizontal grid
!------------------------------------------ 
READ(Unit) np_global
READ(Unit) ne_global

np = np_global;
ne = ne_global;

!!-------------------------------------------------
!!ALLOCATE memory
!!-------------------------------------------------
!IF(ALLOCATED(x))     DEALLOCATE(x);
!IF(ALLOCATED(y))     DEALLOCATE(y);
!IF(ALLOCATED(dp))    DEALLOCATE(dp);
!IF(ALLOCATED(kbp))   DEALLOCATE(kbp);
!
!ALLOCATE(x(np_global),y(np_global),dp(np_global),kbp(np_global),stat=istat);
!if(istat/=0) stop 'Allocation error: x,y,dp,kbp';
!
!
!IF(ALLOCATED(nm)) DEALLOCATE(nm); 
!ALLOCATE(nm(ne_global,3),stat=istat);
!if(istat/=0) stop 'Allocation error: mn';


do m=1,np_global
!  READ(Unit) x(m)
!  READ(Unit) y(m)
!  READ(Unit) dp(m)
!  READ(Unit) kbp(m)
  READ(Unit) x(m),y(m),dp(m),kbp(m);
  
enddo !m=1,np
    
do m=1,ne_global
!  READ(Unit) temp_int
!  do mm=1,3
!    READ(Unit) nm(m,mm)
!  enddo !mm
  READ(Unit) temp_int, nm(m,1:3);

enddo !m

!------------------------------------------ 
!finsihed close the file
!------------------------------------------ 
!CLOSE(65)
!---------------------------------------------------------
!ALLOCATE MEMORY eta2
!---------------------------------------------------------
!IF(ALLOCATED(eta2)) DEALLOCATE(eta2);
!IF(ALLOCATED(outb)) DEALLOCATE(outb);
!
!IF(i23d == 2)THEN
!    ALLOCATE(eta2(np_global),outb(np_global,1,ivs),stat=istat);
!ELSE
!    ALLOCATE(eta2(np_global),outb(np_global,nvrt,ivs),stat=istat);
!ENDIF
!
!if(istat/=0) stop 'Allocation error:eta2,outb';



END SUBROUTINE



!--------------------------------------------------------
!Read the head file of the combined results
!    OPEN THE FILE AND READ THE HEAD DATA, the Unit will kept.
!--------------------------------------------------------
SUBROUTINE acquire_geometry_comb(fileName)
!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
CHARACTER(LEN=*)              :: fileName;

!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
!CHARACTER(LEN=*)           :: fileName;
!--------------------------------------------------------
!Local variables. m,mm;
!--------------------------------------------------------
INTEGER                    :: k,kin,m,mm,temp_int,istat;
INTEGER                    :: Unit =200;
!--------------------------------------------------------
!out put the head of the combined file
!--------------------------------------------------------

OPEN(Unit,file=TRIM(fileName),ACTION="READ",status='OLD',FORM='BINARY',IOSTAT=istat );
IF(istat/=0)WRITE(*,*)"Failed to open file: ", TRIM(fileName);

READ(Unit) data_format  ! data format.
READ(Unit) version      !
READ(Unit) start_time
READ(Unit) variable_nm
READ(Unit) variable_dim

!------------------------------------------ 
!write the record ........
!------------------------------------------
READ(Unit) nrec
READ(Unit) dtout
READ(Unit) nspool
READ(Unit) ivs
READ(Unit) i23d

!------------------------------------------ 
!Vertical grid
!------------------------------------------ 
READ(Unit) nvrt
READ(Unit) kz
READ(Unit) h0
READ(Unit) h_s
READ(Unit) h_c
READ(Unit) theta_b
READ(Unit) theta_f

!-------------------------------------------------
!ALLOCATE memory
!-------------------------------------------------
IF(ALLOCATED(ztot))  DEALLOCATE(ztot);
IF(ALLOCATED(sigma)) DEALLOCATE(sigma);
ALLOCATE(ztot(nvrt),sigma(nvrt),stat=istat);
if(istat/=0) stop 'Allocation error: ztot,sigma';

do k=1,kz-1
  READ(Unit) ztot(k)
enddo

do k=kz,nvrt
  kin=k-kz+1
  READ(Unit) sigma(kin)
enddo !k 

!------------------------------------------ 
!    !Horizontal grid
!------------------------------------------ 
READ(Unit) np_global;
READ(Unit) ne_global;
CLOSE(Unit);

np_global = np;
ne_global = ne;

!
!np = np_global;
!ne = ne_global;
!
!!-------------------------------------------------
!!ALLOCATE memory
!!-------------------------------------------------
!IF(ALLOCATED(x))     DEALLOCATE(x);
!IF(ALLOCATED(y))     DEALLOCATE(y);
!IF(ALLOCATED(dp))    DEALLOCATE(dp);
!IF(ALLOCATED(kbp))   DEALLOCATE(kbp);
!
!ALLOCATE(x(np_global),y(np_global),dp(np_global),kbp(np_global),stat=istat);
!if(istat/=0) stop 'Allocation error: x,y,dp,kbp';
!
!
!IF(ALLOCATED(nm)) DEALLOCATE(nm); 
!ALLOCATE(nm(ne_global,3),stat=istat);
!if(istat/=0) stop 'Allocation error: mn';
!
!
!do m=1,np_global
!  READ(Unit) x(m)
!  READ(Unit) y(m)
!  READ(Unit) dp(m)
!  READ(Unit) kbp(m)
!enddo !m=1,np
!    
!do m=1,ne_global
!  READ(Unit) temp_int
!  do mm=1,3
!    READ(Unit) nm(m,mm)
!  enddo !mm
!enddo !m
!
!!------------------------------------------ 
!!finsihed close the file
!!------------------------------------------ 
!!CLOSE(65)
!!---------------------------------------------------------
!!ALLOCATE MEMORY eta2
!!---------------------------------------------------------
!IF(ALLOCATED(eta2)) DEALLOCATE(eta2);
!IF(ALLOCATED(outb)) DEALLOCATE(outb);
!
!IF(i23d == 2)THEN
!    ALLOCATE(eta2(np_global),outb(np_global,1,ivs),stat=istat);
!ELSE
!    ALLOCATE(eta2(np_global),outb(np_global,nvrt,ivs),stat=istat);
!ENDIF
!
!if(istat/=0) stop 'Allocation error:eta2,outb';



END SUBROUTINE

!--------------------------------------------------------
!subroutine read combine data.
! For sequential read, unit is improtant.
!--------------------------------------------------------
SUBROUTINE in_combine_time_data(Unit)
!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
INTEGER     ,INTENT(IN  )  :: Unit;
!--------------------------------------------------------
!Local variables. m,mm;
!--------------------------------------------------------
INTEGER                    :: i,k,kin,m,mm;
!--------------------------------------------------------
!out put the head of the combined file
!--------------------------------------------------------
!open(65,file=TRIM(fileName),ACCESS='APPEND',FORM='BINARY');
!open(65,file=TRIM(fileName),ACTION="READ",FORM='BINARY');

READ(Unit) time ;
READ(Unit) it   ;

DO i=1,np_global
    READ(Unit) eta2(i)
ENDDO    !i

do i=1,np_global
    if(i23d==2) then
      do m=1,ivs
          READ(Unit) outb(i,1,m)
      enddo !m
    else !i23d=3 
      do k=max0(1,kbp(i)),nvrt
        do m=1,ivs
          READ(Unit) outb(i,k,m)
        enddo !m
      enddo !k
    endif !i23d
enddo !i

!------------------------------------------ 
!finsihed close the file
!------------------------------------------ 
END SUBROUTINE


END MODULE