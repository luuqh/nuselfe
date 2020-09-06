!-------------------------------------------------------------
!1. get horizontal distribution of the scalar and vertical distribution at free surface
!        bottom and averaged values. 
!2. get vertical profile distribution of the scalar and vertical distribution at free surface
!        bottom and averaged values.           
!3. get time series of the scalar and vertical distribution at free surface
!        bottom and averaged values.     
!-------------------------------------------------------------
MODULE post_proc_mod
USE post_vars
USE SELFE_TEC_IO
USE SELFE_TEC_IO_vgrid
USE write_sms_grid_mod
USE combine_output7
USE nsf_self_vars_mod
IMPLICIT NONE  

INTEGER              :: post_proc_out_unit =201;

CONTAINS

!!--------------------------------------------------------------
!!extrac the data and output the data to the required format.
!!--------------------------------------------------------------
!SUBROUTINE data_extraction()
!
!CHARACTER(LEN=200)  :: in_fileName,out_fileName,cday,crank,tec_file,tec_file_prx;
!CHARACTER(LEN=200)  :: outFileName;
!INTEGER             :: iday,irank,ispool,tec_appended,info2;
!INTEGER             :: in_unit = 200;
!INTEGER             :: out_unit= 201;
!CHARACTER(LEN=200)  :: varName,cit;
!INTEGER             :: out_it,istat;
!INTEGER             :: tec_form  = 2 ! =1 the resutls are in separate files (one time step per file).
!                                     ! =2 the results are in one file 
!INTEGER             :: ilrec ;
!
!!---------------------------------------------------------
!!get the vertical output nodelist.
!!---------------------------------------------------------
!IF(i_mode == i_mode_vprof_comb   .OR. &
!   i_mode == i_mode_pTime_comb   .OR. &
!   i_mode == i_mode_vprof_netcdf .OR. &
!   i_mode == i_mode_pTime_netcdf  ) THEN
!   
!    CALL read_node_list();  
!    
!ENDIF
!             
!!---------------------------------------------------------
!!require the number of nodes, elements and vrt layers.
!!get the geometry and allocate the memory.
!!---------------------------------------------------------
!WRITE(cday,'(i12)')ibgn;
!cday  = TRIM(ADJUSTL(cday))//'_';
!in_fileName = TRIM(inFolder)//TRIM(cday)//TRIM(file63);
!CALL acquire_geometry_comb(in_fileName);
!
!!---------------------------------------------------------
!!ALLOCATE THE MEMORY...............
!!---------------------------------------------------------
!CALL alloc_memory();
!
!
!!---------------------------------------------------------
!!start to read the data and process the data...
!!---------------------------------------------------------
!!----------------------------------------------------------
!!READ data for each day step.
!!----------------------------------------------------------
!
!DO iday = ibgn,iend
!
!!-----------------------------------------------------------
!!Read the combined head file.
!!----------------------------------------------------------- 
!!compose the file Name
!    WRITE(cday,'(i12)')iday;
!    cday  = TRIM(ADJUSTL(cday))//'_';
!    in_fileName = TRIM(inFolder)//TRIM(cday)//TRIM(file63);
!    
!    CALL in_combine_head2(in_fileName,in_unit);
!    
!!-----------------------------------------------------------
!! read data for each time step.
!!-----------------------------------------------------------
!
!    DO ispool = 1,nrec
!    
!!-----------------------------------------------------------    
!! READ the time dependent data.
!!-----------------------------------------------------------    
!        CALL in_combine_time_data(in_unit)
!        time = time / (24*60*60);  !-----------change the time to day
!
!!-----------------------------------------------------------    
!!update the z coordinate.    
!!-----------------------------------------------------------    
!
!        CALL update_z();
!        
!!----------------------------------------------------------- 
!!extract the data we need.   
!!-----------------------------------------------------------    
!!-----------horizontal distribution------------------------
!        IF(i_mode  == i_mode_horz_comb   .OR. & 
!           i_mode  == i_mode_horz_netcdf ) THEN
!          
!           CALL get_tec_hgrid_horz(data_out(:,1,:)) ;
!           
!           xx_out       = x;
!           yy_out(:,1)  = y;
!           
!!----------------vertical profile------------------------
!        ELSEIF(i_mode  == i_mode_vprof_comb   .OR. & 
!           i_mode  == i_mode_vprof_netcdf ) THEN           
!                    
!
!
!            
!!----------------------------------------------------------- 
!!output the data......
!!----------------------------------------------------------- 
!!----------------------------------------------------------- 
!!get the tecplot out put var name  Not including x and y
!!----------------------------------------------------------- 
!        IF(ivs ==1) THEN
!            varName = ' " '//TRIM(file63(1:lfile63-3))//' " ';
!        ELSEIF(ivs ==2 ) THEN
!            varName = ' " '//TRIM(file63(1:lfile63-3))//'_x " ,'//' " '//TRIM(file63(1:lfile63-3))//'_y " ';
!        ENDIF        
!
!
!!-----------------------------------------------------------    
!! Get the output tecplot file Name;
!!-----------------------------------------------------------    
!        IF(tec_form  == 1) THEN
!            
!            out_it = out_it + 1;
!            WRITE(cit,"(I12)") out_it;
!            cit = TRIM(ADJUSTL(cit)); 
!            
!            tec_file = TRIM(outFolder)//TRIM(file63(1:lfile63-3))//"_"//TRIM(cit)//"prof.dat"; 
!            tec_appended = 0              
!               
!         ELSE
!            
!            tec_file = TRIM(outFolder)//TRIM(file63(1:lfile63-3))//"prof.dat";
!            
!            IF(iday == ibgn .AND. ispool ==1) THEN
!                tec_appended = 0   
!            ELSE
!                tec_appended = 1   
!            ENDIF               
!          ENDIF
!
!
!!--------extract and proces the data----------------------------------
!!        CALL get_tec_hgrid_horz(out_proc)       
!!        xver(1:n_vprof_node)        = x(iout_node_list(1:n_vprof_node));
!!        zver(1:nvrt,1:n_vprof_node) = z(1:nvrt,iout_node_list(1:n_vprof_node));
!        CALL write_TEC_vGRID(tec_file,time,x(iout_node_list(1:n_vprof_node)) , &
!                z(1:nvrt,iout_node_list(1:n_vprof_node)),               &
!                outb(iout_node_list(1:n_vprof_node),1:nvrt,1:ivs) ,     &
!                TRIM(varName),"xx",tec_appended) ;
!
!        WRITE(*,*) "it =",it,"time = ", time; 
!    
!        ENDDO
!                         
!    ENDDO
!    
!    CLOSE(in_unit);  ! close the combined file.....
!
! 
!ENDDO
!
!!------------------------------------------------------
!!CLOSE THE FILE.....
!!------------------------------------------------------
!!IF(i_out_form == i_out_form_SMS) THEN
!!    CALL close_sms_file (out_fileName,out_unit,ilrec);   
!!ELSE
!!    CLOSE(out_unit);
!!ENDIF  
!!np = np_global;
!!ne = ne_global;
!
!102 FORMAT(100E20.5);
!
!
!
!
!
!
!
!END SUBROUTINE



!--------------------------------------------------------------
!extract_vprof_data: extract the vertical profile data 
!  and output to the tecplot file. 
!   for scalar data no projection needed, for vector data 
!   projection is needed.
!NOTE: extract the velocity profile need use subroutine
!             extract_vprof_vel
!-------------------------------------------------------------

SUBROUTINE extract_vprof_data()
IMPLICIT NONE        ;

CHARACTER(LEN=200)  :: in_fileName,out_fileName,cday,crank,tec_file,tec_file_prx;
CHARACTER(LEN=200)  :: outFileName;
INTEGER             :: iday,irank,ispool,tec_appended,info2;
INTEGER             :: in_unit = 200;
INTEGER             :: out_unit= 201;
CHARACTER(LEN=200)  :: varName,cit;
INTEGER             :: out_it,istat;
INTEGER             :: tec_form  = 2 ! =1 the resutls are in separate files (one time step per file).
                                     ! =2 the results are in one file 
INTEGER             :: ilrec ;

tec_form  = 2 
!---------------------------------------------------------
!read the parameters from file combine_output.in
!---------------------------------------------------------
!CALL get_horz_parm();

!---------------------------------------------------------
!get the vertical output nodelist.
!---------------------------------------------------------
CALL read_node_list();

n_vprof_node = nout_node_list;
!---------------------------------------------------------
!get the geometry and allocate the memory.
!---------------------------------------------------------
WRITE(cday,'(i12)')ibgn;
cday  = TRIM(ADJUSTL(cday))//'_';
in_fileName = TRIM(inFolder)//TRIM(cday)//TRIM(file63);
CALL in_combine_head(in_fileName,in_unit);
CLOSE(in_unit);


!---------------------------------------------------------
!Allocate the memory..
!---------------------------------------------------------


IF(ALLOCATED(xvert))  DEALLOCATE(xvert);
IF(ALLOCATED(zvert))  DEALLOCATE(zvert);
ALLOCATE(xvert(n_vprof_node),zvert(nvrt,n_vprof_node),stat=istat)
IF(istat/=0) STOP 'Allocation error: xvert,zvert ';



!---------------------------------------------------------
!------set the tecplot file name-------------
!---------------------------------------------------------
out_fileName  = TRIM(outFolder)//"out";
out_it        = 0;

!----------------------------------------------------------
!READ data for each day step.
!----------------------------------------------------------

DO iday = ibgn,iend

!compose the file Name
    WRITE(cday,'(i12)')iday;
    cday  = TRIM(ADJUSTL(cday))//'_';

    in_fileName = TRIM(inFolder)//TRIM(cday)//TRIM(file63);
    
!-----------------------------------------------------------
!Read the combined head file.
!----------------------------------------------------------- 
    CALL in_combine_head(in_fileName,in_unit);
    

  
    DO ispool = 1,nrec
    

!-----------------------------------------------------------    
! READ the time dependent data.
!-----------------------------------------------------------    
        CALL in_combine_time_data(in_unit)
        time = time / (24*60*60);  !-----------change the time to day
        
!------------------------------------------------------------ 
!       allocate the memory       
!------------------------------------------------------------        
        IF(iday == ibgn .AND. ispool == 1) THEN
            IF(ALLOCATED(hmod))  DEALLOCATE(hmod);
            IF(ALLOCATED(cs))    DEALLOCATE(cs);
            IF(ALLOCATED(kfp))   DEALLOCATE(kfp);
            IF(ALLOCATED(htot))  DEALLOCATE(htot);
            IF(ALLOCATED(z))     DEALLOCATE(z);
            IF(ALLOCATED(delz))     DEALLOCATE(delz);           
            IF(ALLOCATED(out_proc)) DEALLOCATE(out_proc);
            
            ALLOCATE(hmod(np_global),cs(nvrt),kfp(np_global),htot(np_global),stat=istat)
            if(istat/=0) stop 'Allocation error: hmod,cs,kfp,htot';
            
            ALLOCATE(z(nvrt,np_global),delz(nvrt,np_global),stat=istat)
            if(istat/=0) stop 'Allocation error: z,delz';
                       
            ALLOCATE(out_proc(np_global,ivs),stat=istat)
            if(istat/=0) stop 'Allocation error: outb'; 
        ENDIF
!-----------------------------------------------------------    
!process the data................
!-----------------------------------------------------------    
!--------update the z coordinate.    
        CALL update_z;
        
            
!----------------------------------------------------------- 
!output the data......
!----------------------------------------------------------- 
!----------------------------------------------------------- 
!get the tecplot out put var name  Not including x and y
!----------------------------------------------------------- 
        IF(ivs ==1) THEN
            varName = ' " '//TRIM(file63(1:lfile63-3))//' " ';
        ELSEIF(ivs ==2 ) THEN
            varName = ' " '//TRIM(file63(1:lfile63-3))//'_x " ,'//' " '//TRIM(file63(1:lfile63-3))//'_y " ';
        ENDIF        


!-----------------------------------------------------------    
! Get the output tecplot file Name;
!-----------------------------------------------------------    
        IF(tec_form  == 1) THEN
            
            out_it = out_it + 1;
            WRITE(cit,"(I12)") out_it;
            cit = TRIM(ADJUSTL(cit)); 
            
            tec_file = TRIM(outFolder)//TRIM(file63(1:lfile63-3))//"_"//TRIM(cit)//"prof.dat"; 
            tec_appended = 0              
               
         ELSE
            
            tec_file = TRIM(outFolder)//TRIM(file63(1:lfile63-3))//"prof.dat";
            
            IF(iday == ibgn .AND. ispool ==1) THEN
                tec_appended = 0   
            ELSE
                tec_appended = 1   
            ENDIF               
          ENDIF


!--------extract and proces the data----------------------------------
!        CALL get_tec_hgrid_horz(out_proc)       
!        xver(1:n_vprof_node)        = x(iout_node_list(1:n_vprof_node));
!        zver(1:nvrt,1:n_vprof_node) = z(1:nvrt,iout_node_list(1:n_vprof_node));
        CALL write_TEC_vGRID(tec_file,time,x(iout_node_list(1:n_vprof_node)) , &
                z(1:nvrt,iout_node_list(1:n_vprof_node)),               &
                outb(iout_node_list(1:n_vprof_node),1:nvrt,1:ivs) ,     &
                TRIM(varName),"xx",tec_appended) ;

        WRITE(*,*) "it =",it,"time = ", time;          
    ENDDO
    
    CLOSE(in_unit);  ! close the combined file.....

 
ENDDO

!------------------------------------------------------
!CLOSE THE FILE.....
!------------------------------------------------------
!IF(i_out_form == i_out_form_SMS) THEN
!    CALL close_sms_file (out_fileName,out_unit,ilrec);   
!ELSE
!    CLOSE(out_unit);
!ENDIF  
!np = np_global;
!ne = ne_global;

102 FORMAT(100E20.5);

END SUBROUTINE

!--------------------------------------------------------------
!extract_vprof_vel: extract the velocity at the vertical profile
!      1.  the x component of the velocity is the magnitude or
!          the x component of the horizontal velocity, the z 
!          component is  vertical velcoity.
!      2.  the x location can be the x location of the node or 
!          the distance of the profile line from this point.
!-------------------------------------------------------------

SUBROUTINE extract_vprof_vel()
IMPLICIT NONE        ;

CHARACTER(LEN=200)  :: in_fileName,out_fileName,cday,crank,tec_file,tec_file_prx;
CHARACTER(LEN=200)  :: outFileName;
INTEGER             :: iday,irank,ispool,tec_appended,info2;
INTEGER             :: in_unit = 200;
INTEGER             :: out_unit= 201;
CHARACTER(LEN=200)  :: varName,cit;
INTEGER             :: out_it,istat;
INTEGER             :: tec_form  = 1 ! =1 the resutls are in separate files (one time step per file).
                                     ! =2 the results are in one file 
INTEGER             :: ilrec ;
CHARACTER(LEN=200)  :: hvel, vert;
INTEGER             :: hvel_unit,vert_unit;

hvel = 'hvel.64'; !horizontal velocity file name.
vert = 'vert.63'; !vertical velocity file name.
hvel_unit = 201;
vert_unit = 202;

tec_form  = 2
!---------------------------------------------------------
!read the parameters from file combine_output.in
!---------------------------------------------------------
!CALL get_horz_parm();

!---------------------------------------------------------
!get the vertical output nodelist.
!---------------------------------------------------------
CALL read_node_list();

n_vprof_node = nout_node_list;

!---------------------------------------------------------
!get the geometry and allocate the memory.
!---------------------------------------------------------
file63 = hvel;
WRITE(cday,'(i12)')ibgn;
cday  = TRIM(ADJUSTL(cday))//'_';
in_fileName = TRIM(inFolder)//TRIM(cday)//TRIM(file63);
CALL in_combine_head(in_fileName,hvel_unit);
CLOSE(hvel_unit);


!---------------------------------------------------------
!Allocate the memory..
!---------------------------------------------------------


IF(ALLOCATED(xvert))  DEALLOCATE(xvert);
IF(ALLOCATED(zvert))  DEALLOCATE(zvert);
IF(ALLOCATED(disvert_vel))    DEALLOCATE(disvert_vel) ;
IF(ALLOCATED(disvert))         DEALLOCATE(disvert) ;

ALLOCATE(xvert(n_vprof_node),zvert(nvrt,n_vprof_node),stat=istat);
ALLOCATE(disvert_vel(n_vprof_node,nvrt,2),&
         disvert(0:n_vprof_node+1),stat=istat);
IF(istat/=0) STOP 'Allocation error: xvert,zvert ';



!---------------------------------------------------------
!------set the tecplot file name-------------
!---------------------------------------------------------
out_fileName  = TRIM(outFolder)//"out";
out_it        = 0;

!----------------------------------------------------------
!READ data for each day step.
!----------------------------------------------------------

DO iday = ibgn,iend

!-----------------------------------------------------------
!deal with the vertical velocity . 
!-----------------------------------------------------------
!compose the file Name
    file63  = vert;
    WRITE(cday,'(i12)')iday;
    cday  = TRIM(ADJUSTL(cday))//'_';
    in_fileName = TRIM(inFolder)//TRIM(cday)//TRIM(file63);
!-----------------------------------------------------------
!Read the combined head file.
!----------------------------------------------------------- 
    CALL in_combine_head(in_fileName,vert_unit);


!-----------------------------------------------------------
!deal with the horizontal velocity . 
!-----------------------------------------------------------
!compose the file Name
    file63  = hvel;
    WRITE(cday,'(i12)')iday;
    cday  = TRIM(ADJUSTL(cday))//'_';
    in_fileName = TRIM(inFolder)//TRIM(cday)//TRIM(file63);
!-----------------------------------------------------------
!Read the combined head file.
!----------------------------------------------------------- 
    CALL in_combine_head(in_fileName,hvel_unit);
    


  
    DO ispool = 1,nrec
    

!-----------------------------------------------------------    
! READ the time dependent horizontal velocity.
!-----------------------------------------------------------
        ivs = 2;    
        CALL in_combine_time_data(hvel_unit)
        time = time / (24*60*60);  !-----------change the time to day

!------------------------------------------------------------ 
!get the velocity magnitude.
!------------------------------------------------------------ 
        disvert_vel(:,:,1) = SQRT(outb(iout_node_list(1:n_vprof_node),:,1)**2 + &
                                     outb(iout_node_list(1:n_vprof_node),:,2)**2 );                   

!-----------------------------------------------------------    
! READ the time dependent vertical velocity.
!-----------------------------------------------------------
        ivs = 1;    
        CALL in_combine_time_data(vert_unit)
        time = time / (24*60*60);  !-----------change the time to day

!------------------------------------------------------------ 
!get the velocity magnitude.
!------------------------------------------------------------ 
        disvert_vel(:,:,2) = outb(iout_node_list(1:n_vprof_node),:,1)

!------------------------------------------------------------ 
!       allocate the memory       
!------------------------------------------------------------        
        IF(iday == ibgn .AND. ispool == 1) THEN
            IF(ALLOCATED(hmod))  DEALLOCATE(hmod);
            IF(ALLOCATED(cs))    DEALLOCATE(cs);
            IF(ALLOCATED(kfp))   DEALLOCATE(kfp);
            IF(ALLOCATED(htot))  DEALLOCATE(htot);
            IF(ALLOCATED(z))     DEALLOCATE(z);
            IF(ALLOCATED(delz))     DEALLOCATE(delz);           
            IF(ALLOCATED(out_proc)) DEALLOCATE(out_proc);
            
            ALLOCATE(hmod(np_global),cs(nvrt),kfp(np_global),htot(np_global),stat=istat)
            if(istat/=0) stop 'Allocation error: hmod,cs,kfp,htot';
            
            ALLOCATE(z(nvrt,np_global),delz(nvrt,np_global),stat=istat)
            if(istat/=0) stop 'Allocation error: z,delz';
                       
            ALLOCATE(out_proc(np_global,ivs),stat=istat)
            if(istat/=0) stop 'Allocation error: outb'; 
        ENDIF
       
!-----------------------------------------------------------    
!process the data................
!-----------------------------------------------------------    
!--------update the z coordinate.    
        CALL update_z;
        
!-----------------------------------------------------------    
!get the x and z coordinate of the node....
! for project case or the profile length cases.
!-----------------------------------------------------------    
        xvert(1:n_vprof_node) = x(iout_node_list(1:n_vprof_node));
        zvert(1:nvrt,1:n_vprof_node) = z(1:nvrt,iout_node_list(1:n_vprof_node));
           
!----------------------------------------------------------- 
!output the data......
!----------------------------------------------------------- 
!----------------------------------------------------------- 
!get the tecplot out put var name  Not including x and y
!----------------------------------------------------------- 
!        IF(ivs ==1) THEN
!            varName = ' " '//TRIM(file63(1:lfile63-3))//' " ';
!        ELSEIF(ivs ==2 ) THEN
!            varName = ' " '//TRIM(file63(1:lfile63-3))//'_x " ,'//' " '//TRIM(file63(1:lfile63-3))//'_y " ';
!        ENDIF        
        varName = ' " '//TRIM(hvel)//' " ,'//' " '//TRIM(vert)//'" ';
!-----------------------------------------------------------    
! Get the output tecplot file Name;
!-----------------------------------------------------------    
        IF(tec_form  == 1) THEN
            
            out_it = out_it + 1;
            WRITE(cit,"(I12)") out_it;
            cit = TRIM(ADJUSTL(cit)); 
            
            tec_file = TRIM(outFolder)//TRIM(file63(1:lfile63-3))//"_"//TRIM(cit)//".dat"; 
            tec_appended = 0              
               
         ELSE
            
            tec_file = TRIM(outFolder)//TRIM(file63(1:lfile63-3))//".dat";
            
            IF(iday == ibgn .AND. ispool ==1) THEN
                tec_appended = 0   
            ELSE
                tec_appended = 1   
            ENDIF               
          ENDIF


!--------extract and proces the data----------------------------------
!        CALL get_tec_hgrid_horz(out_proc)       
!        xver(1:n_vprof_node)        = x(iout_node_list(1:n_vprof_node));
!        zver(1:nvrt,1:n_vprof_node) = z(1:nvrt,iout_node_list(1:n_vprof_node));
        ivs = 2;
        CALL write_TEC_vGRID(tec_file,time,xvert,zvert,disvert_vel,     &
                TRIM(varName),"xx",tec_appended) ;

        WRITE(*,*) "it =",it,"time = ", time;          
    ENDDO
    
    CLOSE(in_unit);  ! close the combined file.....

 
ENDDO

!------------------------------------------------------
!CLOSE THE FILE.....
!------------------------------------------------------
!IF(i_out_form == i_out_form_SMS) THEN
!    CALL close_sms_file (out_fileName,out_unit,ilrec);   
!ELSE
!    CLOSE(out_unit);
!ENDIF  
!np = np_global;
!ne = ne_global;

102 FORMAT(100E20.5);

END SUBROUTINE


!--------------------------------------------------------------
!extract the horizontal data from the combined results and 
! out put them to tecplot format.
!-------------------------------------------------------------
SUBROUTINE extract_horz_data()
IMPLICIT NONE        ;

CHARACTER(LEN=200)  :: in_fileName,out_fileName,out_fileName_base,cday,crank,tec_file,tec_file_prx;
CHARACTER(LEN=200)  :: outFileName;
INTEGER             :: iday,irank,ispool,tec_appended,info2;
INTEGER             :: in_unit = 200;
INTEGER             :: out_unit= 201;
CHARACTER(LEN=200)  :: varName,cit;
INTEGER             :: out_it,istat;
!INTEGER             :: tec_form  = 1 ! =1 the resutls are in separate files (one time step per file).
!                                     ! =2 the results are in one file 
INTEGER             :: ilrec ;
INTEGER             :: data_type

!---------------------------------------------------------
!read the parameters from file combine_output.in
!---------------------------------------------------------
!CALL get_horz_parm();

!---------------------------------------------------------
!get the specified layer.
!---------------------------------------------------------
!------set the tecplot file name-------------


!varName   = file63(1:lfile63-3);
!----------------------------------------------------------
!READ data for each day step.
!----------------------------------------------------------
out_it             = 0;    
out_fileName_base  = TRIM(outFolder)//TRIM(file63(1:lfile63-3));

DO iday = ibgn,iend

!compose the file Name
    WRITE(cday,'(i12)')iday;
    cday  = TRIM(ADJUSTL(cday))//'_';

    in_fileName = TRIM(inFolder)//TRIM(cday)//TRIM(file63);
    
!-----------------------------------------------------------
!Read the combined head file.
!----------------------------------------------------------- 
    CALL in_combine_head(in_fileName,in_unit);
    

    IF(iday == ibgn ) THEN
    
        IF(i_out_form == i_out_form_SMS)THEN
        
            out_fileName  = TRIM(out_fileName)//'_horz_sms.dat';
            CALL write_sms_head(out_fileName,out_unit,ilrec);
        ENDIF
        
    ENDIF        

   
    DO ispool = 1,nrec
    

!-----------------------------------------------------------    
! READ the time dependent data.
!-----------------------------------------------------------    
        CALL in_combine_time_data(in_unit)
        time = time / (24*60*60);  !-----------change the time to day
        
!------------------------------------------------------------ 
!       allocate the memory       
!------------------------------------------------------------        
        IF(iday == ibgn .AND. ispool == 1) THEN
            IF(ALLOCATED(hmod))  DEALLOCATE(hmod);
            IF(ALLOCATED(cs))    DEALLOCATE(cs);
            IF(ALLOCATED(kfp))   DEALLOCATE(kfp);
            IF(ALLOCATED(htot))  DEALLOCATE(htot);
            IF(ALLOCATED(z))     DEALLOCATE(z);
            IF(ALLOCATED(delz))     DEALLOCATE(delz);           
            IF(ALLOCATED(out_proc)) DEALLOCATE(out_proc);
            
            ALLOCATE(hmod(np_global),cs(nvrt),kfp(np_global),htot(np_global),stat=istat)
            if(istat/=0) stop 'Allocation error: hmod,cs,kfp,htot';
            
            ALLOCATE(z(nvrt,np_global),delz(nvrt,np_global),stat=istat)
            if(istat/=0) stop 'Allocation error: z,delz';
                       
            ALLOCATE(out_proc(np_global,ivs),stat=istat)
            if(istat/=0) stop 'Allocation error: outb'; 
        ENDIF
!-----------------------------------------------------------    
!process the data................
!-----------------------------------------------------------    
!--------update the z coordinate.    
        CALL update_z;
        
!--------extract and proces the data----------------------------------
        CALL get_tec_hgrid_horz(out_proc)       


!------------------------------------------------------------------
!output the data.
!------------------------------------------------------------------    
        IF(i_out_form == i_out_form_SMS)THEN ! Output to sms format.
        
            CALL write_sms_time_data(out_fileName,out_unit,ilrec,out_proc);
        
        ELSE
!-----------------------------------------------------------    
! Get the output tecplot file Name;
!-----------------------------------------------------------    
            IF(tec_form  == 1) THEN
            
                out_it = out_it + 1;
                WRITE(cit,"(I12)") out_it;
                cit = TRIM(ADJUSTL(cit)); 
                
                tec_file = TRIM(out_fileName_base)//TRIM(cit)//"_horz_tec.dat"; 
                tec_appended = 0              
               
            ELSE
            
                tec_file = TRIM(out_fileName_base)//"_horz_tec.dat";
                
                IF(iday == ibgn .AND. ispool ==1) THEN
                    tec_appended = 0   
                ELSE
                    tec_appended = 1   
                ENDIF               
            
            ENDIF

        
!----------------------------------------------------------- 
!get the tecplot out put var name
!  Not including x and y
!----------------------------------------------------------- 
        IF(ivs ==1) THEN
            varName = ' " '//TRIM(file63(1:lfile63-3))//' " ';
        ELSEIF(ivs ==2 ) THEN
            varName = ' " '//TRIM(file63(1:lfile63-3))//'_x " ,'//' " '//TRIM(file63(1:lfile63-3))//'_y " ';
        ENDIF

!----------------------------------------------------------- 
!set the out put data.......................
!----------------------------------------------------------- 
!        tec_out  => outb(:,1,1:ivs);
        IF(nout_node_list ==0 ) THEN !out put all the data..........
            CALL write_TEC_hGRID(tec_file,time,out_proc(:,1:ivs),TRIM(varName),"xx",tec_appended);
        
        ELSE ! out the the time series to text file.
            OPEN(UNIT=101,FILE=TRIM(outFolder)//TRIM(file63)//"_point.txt",ACTION='WRITE', &
                 & POSITION='APPEND', IOSTAT=info2);   
            WRITE(101,102)time, out_proc(iout_node_list(1:nout_node_list),1:ivs); 
            CLOSE(101);                   
        ENDIF
        
!        WRITE(*,*) "it =",it,"time = ", time;
        
     ENDIF
     
     
    WRITE(*,*) "it =",it,"time = ", time;        
    
    ENDDO     !DO ispool = 1,nrec



    
    CLOSE(in_unit);  ! close the combined file.....

 
ENDDO

!------------------------------------------------------
!CLOSE THE FILE.....
!------------------------------------------------------
IF(i_out_form == i_out_form_SMS) THEN
    CALL close_sms_file (out_fileName,out_unit,ilrec);   
ELSE
!    CLOSE(out_unit);
ENDIF  
!np = np_global;
!ne = ne_global;

102 FORMAT(100E20.5);

END SUBROUTINE



!--------------------------------------------------------------
!extract the horizontal data from the combined results and 
! out put them to tecplot format.
!-------------------------------------------------------------
SUBROUTINE extract_horz_netcdf()
IMPLICIT NONE        ;

CHARACTER(LEN=200)  :: in_fileName,out_fileName,out_fileName_base,cday,crank,tec_file,tec_file_prx;
CHARACTER(LEN=200)  :: outFileName;
INTEGER             :: iday,irank,ispool,tec_appended,info2;
INTEGER             :: in_unit = 200;
INTEGER             :: out_unit= 201;
CHARACTER(LEN=200)  :: varName,cit;
INTEGER             :: out_it,istat;
!INTEGER             :: tec_form  = 1 ! =1 the resutls are in separate files (one time step per file).
!                                     ! =2 the results are in one file 
INTEGER             :: ilrec ;
INTEGER             :: data_type
INTEGER             :: ndata(1:3);
!---------------------------------------------------------
!read the parameters from file combine_output.in
!---------------------------------------------------------
!CALL get_horz_parm();

!---------------------------------------------------------
!get the specified layer.
!---------------------------------------------------------
!------set the tecplot file name-------------


!varName   = file63(1:lfile63-3);
!----------------------------------------------------------
!READ data for each day step.
!----------------------------------------------------------
out_it             = 0;    
out_fileName_base  = TRIM(outFolder)//TRIM(file63(1:lfile63-3));

DO iday = ibgn,iend

!compose the file Name
    WRITE(cday,'(i12)')iday;
    cday  = TRIM(ADJUSTL(cday));

    in_fileName = TRIM(inFolder)//TRIM(cday)//'.nc';

!!-----------------------------------------------------------
!!Read the combined head file.
!!----------------------------------------------------------- 
!    CALL in_combine_head(in_fileName,in_unit);
    CALL nfs_open_file(nsf_self,TRIM(in_fileName),NF90_NOWRITE) ;
    
    CALL nsf_selfe_acquire_dims();

!----------------------------------------------------------- 
!Allocate the memory..........................................
!----------------------------------------------------------- 
    CALL alloc_memory();

!----------------------------------------------------------- 
!get the geometry information...
!----------------------------------------------------------- 
    CALL nsf_selfe_in_grid();
       

    IF(iday == ibgn ) THEN
    
        IF(i_out_form == i_out_form_SMS)THEN
        
            out_fileName  = TRIM(out_fileName)//'_horz_sms.dat';
            CALL write_sms_head(out_fileName,out_unit,ilrec);
        ENDIF
        
    ENDIF        

   
    DO ispool = 1,nrec


!-----------------------------------------------------------    
!read simulation time 
!-----------------------------------------------------------    
        CALL nsf_self_netcdf_time_time('time',time,ispool);   
!-----------------------------------------------------------    
!read the elevation from the file.
!-----------------------------------------------------------
        ndata(1) = np;    
        CALL nsf_self_netcdf_time_data_in('elev',eta2,ispool,ndata);

!--------------------------------------------------------------------
!READ the time dependent data from the file . 
!-------------------------------------------------------------------- 
        IF(i23d ==2) THEN  !2D data;
             ndata(1) = np;            
             CALL nsf_self_netcdf_time_data_in(TRIM(file63(1:lfile63-3)),outb(:,1,1),ispool,ndata);
        ELSEIF(i23d == 3) THEN  !3D data;
        
            IF(ivs == 1 )THEN   !3D scalar.
               
               ndata(1:2) = (/np,nvrt/);  
               CALL nsf_self_netcdf_time_data_in(TRIM(file63(1:lfile63-3)),outb(:,:,1),ispool,ndata);
             
            
            ELSE
                ndata(1:3) = (/nvrt,np,ivs/);  
                CALL nsf_self_netcdf_time_data_in(TRIM(file63(1:lfile63-3)),outb(:,:,:),ispool,ndata);
           
            ENDIF
        ENDIF

        
!-----------------------------------------------------------    
!process the data................
!-----------------------------------------------------------    
!--------update the z coordinate.    
        kbp = 1;
        CALL update_z;
        
!--------extract and proces the data----------------------------------
        CALL get_tec_hgrid_horz(out_proc)       


!------------------------------------------------------------------
!output the data.
!------------------------------------------------------------------    
        IF(i_out_form == i_out_form_SMS)THEN ! Output to sms format.
        
            CALL write_sms_time_data(out_fileName,out_unit,ilrec,out_proc);
        
        ELSE
!-----------------------------------------------------------    
! Get the output tecplot file Name;
!-----------------------------------------------------------    
            IF(tec_form  == 1) THEN
            
                out_it = out_it + 1;
                WRITE(cit,"(I12)") out_it;
                cit = TRIM(ADJUSTL(cit)); 
                
                tec_file = TRIM(out_fileName_base)//TRIM(cit)//"_horz_tec.dat"; 
                tec_appended = 0              
               
            ELSE
            
                tec_file = TRIM(out_fileName_base)//"_horz_tec.dat";
                
                IF(iday == ibgn .AND. ispool ==1) THEN
                    tec_appended = 0   
                ELSE
                    tec_appended = 1   
                ENDIF               
            
            ENDIF

        
!----------------------------------------------------------- 
!get the tecplot out put var name
!  Not including x and y
!----------------------------------------------------------- 
        IF(ivs ==1) THEN
            varName = ' " '//TRIM(file63(1:lfile63-3))//' " ';
        ELSEIF(ivs ==2 ) THEN
            varName = ' " '//TRIM(file63(1:lfile63-3))//'_x " ,'//' " '//TRIM(file63(1:lfile63-3))//'_y " ';
        ENDIF

!----------------------------------------------------------- 
!set the out put data.......................
!----------------------------------------------------------- 
!        tec_out  => outb(:,1,1:ivs);
        IF(nout_node_list ==0 ) THEN !out put all the data..........
            CALL write_TEC_hGRID(tec_file,time,out_proc(:,1:ivs),TRIM(varName),"xx",tec_appended);
        
        ELSE ! out the the time series to text file.
            OPEN(UNIT=101,FILE=TRIM(outFolder)//TRIM(file63)//"_point.txt",ACTION='WRITE', &
                 & POSITION='APPEND', IOSTAT=info2);   
            WRITE(101,102)time, out_proc(iout_node_list(1:nout_node_list),1:ivs); 
            CLOSE(101);                   
        ENDIF
        
!        WRITE(*,*) "it =",it,"time = ", time;
        
     ENDIF
     
     
    WRITE(*,*) "it =",it,"time = ", time;        
    
    ENDDO     !DO ispool = 1,nrec



    
!    CLOSE(in_unit);  ! close the combined file.....
    
    CALL nfs_close_file(nsf_self);
 
ENDDO

!------------------------------------------------------
!CLOSE THE FILE.....
!------------------------------------------------------
IF(i_out_form == i_out_form_SMS) THEN
    CALL close_sms_file (out_fileName,out_unit,ilrec);   
ELSE
!    CLOSE(out_unit);
ENDIF  
!np = np_global;
!ne = ne_global;

102 FORMAT(100E20.5);

END SUBROUTINE



!--------------------------------------------------------------
!extract the horizontal data from the combined results and 
! out put them to tecplot format.
!-------------------------------------------------------------
SUBROUTINE extract_vprof_netcdf_data()
IMPLICIT NONE        ;

CHARACTER(LEN=200)  :: in_fileName,out_fileName,out_fileName_base,cday,crank,tec_file,tec_file_prx;
CHARACTER(LEN=200)  :: outFileName;
INTEGER             :: iday,irank,ispool,tec_appended,info2;
INTEGER             :: in_unit = 200;
INTEGER             :: out_unit= 201;
CHARACTER(LEN=200)  :: varName,cit;
INTEGER             :: out_it,istat;
!INTEGER             :: tec_form  = 1 ! =1 the resutls are in separate files (one time step per file).
!                                     ! =2 the results are in one file 
INTEGER             :: ilrec ;
INTEGER             :: data_type
INTEGER             :: ndata(1:3);
!---------------------------------------------------------
!read the parameters from file combine_output.in
!---------------------------------------------------------
!CALL get_horz_parm();

!---------------------------------------------------------
!get the specified layer.
!---------------------------------------------------------
!------set the tecplot file name-------------


!varName   = file63(1:lfile63-3);
!----------------------------------------------------------
!READ data for each day step.
!----------------------------------------------------------
out_it             = 0;    
out_fileName_base  = TRIM(outFolder)//TRIM(file63(1:lfile63-3));

n_vprof_node  = nout_node_list;

DO iday = ibgn,iend

!compose the file Name
    WRITE(cday,'(i12)')iday;
    cday  = TRIM(ADJUSTL(cday));

    in_fileName = TRIM(inFolder)//TRIM(cday)//'.nc';

!!-----------------------------------------------------------
!!Read the combined head file.
!!----------------------------------------------------------- 
!    CALL in_combine_head(in_fileName,in_unit);
    CALL nfs_open_file(nsf_self,TRIM(in_fileName),NF90_NOWRITE) ;
    
    CALL nsf_selfe_acquire_dims();

!----------------------------------------------------------- 
!Allocate the memory..........................................
!----------------------------------------------------------- 
    CALL alloc_memory();

!----------------------------------------------------------- 
!get the geometry information...
!----------------------------------------------------------- 
    CALL nsf_selfe_in_grid();
       

    IF(iday == ibgn ) THEN
    
        IF(i_out_form == i_out_form_SMS)THEN
        
            out_fileName  = TRIM(out_fileName)//'_horz_sms.dat';
            CALL write_sms_head(out_fileName,out_unit,ilrec);
        ENDIF
        
    ENDIF        

   
    DO ispool = 1,nrec


!-----------------------------------------------------------    
!read simulation time 
!-----------------------------------------------------------    
        CALL nsf_self_netcdf_time_time('time',time,ispool);   
!-----------------------------------------------------------    
!read the elevation from the file.
!-----------------------------------------------------------
        ndata(1) = np;    
        CALL nsf_self_netcdf_time_data_in('elev',eta2,ispool,ndata);

!--------------------------------------------------------------------
!READ the time dependent data from the file . 
!-------------------------------------------------------------------- 
        IF(i23d ==2) THEN  !2D data;
             ndata(1) = np;            
             CALL nsf_self_netcdf_time_data_in(TRIM(file63(1:lfile63-3)),outb(:,1,1),ispool,ndata);
        ELSEIF(i23d == 3) THEN  !3D data;
        
            IF(ivs == 1 )THEN   !3D scalar.
               
               ndata(1:2) = (/np,nvrt/);  
               CALL nsf_self_netcdf_time_data_in(TRIM(file63(1:lfile63-3)),outb(:,:,1),ispool,ndata);
             
            
            ELSE
                ndata(1:3) = (/nvrt,np,ivs/);  
                CALL nsf_self_netcdf_time_data_in(TRIM(file63(1:lfile63-3)),outb(:,:,:),ispool,ndata);
           
            ENDIF
        ENDIF

        
!-----------------------------------------------------------    
!process the data................
!-----------------------------------------------------------    
!--------update the z coordinate.    
!        kbp = 1;
        CALL update_z;

!----------------------------------------------------------- 
!get the tecplot out put var name  Not including x and y
!----------------------------------------------------------- 
        IF(ivs ==1) THEN
            varName = ' " '//TRIM(file63(1:lfile63-3))//' " ';
        ELSEIF(ivs ==2 ) THEN
            varName = ' " '//TRIM(file63(1:lfile63-3))//'_x " ,'//' " '//TRIM(file63(1:lfile63-3))//'_y " ';
        ENDIF 
                
!-----------------------------------------------------------    
! Get the output tecplot file Name;
!-----------------------------------------------------------    
        IF(tec_form  == 1) THEN
            
            out_it = out_it + 1;
            WRITE(cit,"(I12)") out_it;
            cit = TRIM(ADJUSTL(cit)); 
            
            tec_file = TRIM(outFolder)//TRIM(file63(1:lfile63-3))//"_"//TRIM(cit)//"prof.dat"; 
            tec_appended = 0              
               
         ELSE
            
            tec_file = TRIM(outFolder)//TRIM(file63(1:lfile63-3))//"_prof.dat";
            
            IF(iday == ibgn .AND. ispool ==1) THEN
                tec_appended = 0   
            ELSE
                tec_appended = 1   
            ENDIF               
          ENDIF        
!--------extract and proces the data----------------------------------
!        CALL get_tec_hgrid_horz(out_proc)       
!        xver(1:n_vprof_node)        = x(iout_node_list(1:n_vprof_node));
!        zver(1:nvrt,1:n_vprof_node) = z(1:nvrt,iout_node_list(1:n_vprof_node));
        CALL write_TEC_vGRID(tec_file,time,x(iout_node_list(1:n_vprof_node)) , &
                z(1:nvrt,iout_node_list(1:n_vprof_node)),               &
                outb(iout_node_list(1:n_vprof_node),1:nvrt,1:ivs) ,     &
                TRIM(varName),"xx",tec_appended) ;


!!------------------------------------------------------------------
!!output the data.
!!------------------------------------------------------------------    
!        IF(i_out_form == i_out_form_SMS)THEN ! Output to sms format.
!        
!            CALL write_sms_time_data(out_fileName,out_unit,ilrec,out_proc);
!        
!        ELSE
!!-----------------------------------------------------------    
!! Get the output tecplot file Name;
!!-----------------------------------------------------------    
!            IF(tec_form  == 1) THEN
!            
!                out_it = out_it + 1;
!                WRITE(cit,"(I12)") out_it;
!                cit = TRIM(ADJUSTL(cit)); 
!                
!                tec_file = TRIM(out_fileName_base)//TRIM(cit)//"_horz_tec.dat"; 
!                tec_appended = 0              
!               
!            ELSE
!            
!                tec_file = TRIM(out_fileName_base)//"_horz_tec.dat";
!                
!                IF(iday == ibgn .AND. ispool ==1) THEN
!                    tec_appended = 0   
!                ELSE
!                    tec_appended = 1   
!                ENDIF               
!            
!            ENDIF
!
!        
!!----------------------------------------------------------- 
!!get the tecplot out put var name
!!  Not including x and y
!!----------------------------------------------------------- 
!        IF(ivs ==1) THEN
!            varName = ' " '//TRIM(file63(1:lfile63-3))//' " ';
!        ELSEIF(ivs ==2 ) THEN
!            varName = ' " '//TRIM(file63(1:lfile63-3))//'_x " ,'//' " '//TRIM(file63(1:lfile63-3))//'_y " ';
!        ENDIF
!
!!----------------------------------------------------------- 
!!set the out put data.......................
!!----------------------------------------------------------- 
!!        tec_out  => outb(:,1,1:ivs);
!        IF(nout_node_list ==0 ) THEN !out put all the data..........
!            CALL write_TEC_hGRID(tec_file,time,out_proc(:,1:ivs),TRIM(varName),"xx",tec_appended);
!        
!        ELSE ! out the the time series to text file.
!            OPEN(UNIT=101,FILE=TRIM(outFolder)//TRIM(file63)//"_point.txt",ACTION='WRITE', &
!                 & POSITION='APPEND', IOSTAT=info2);   
!            WRITE(101,102)time, out_proc(iout_node_list(1:nout_node_list),1:ivs); 
!            CLOSE(101);                   
!        ENDIF
!        
!!        WRITE(*,*) "it =",it,"time = ", time;
        
!    ENDIF
     
     
    WRITE(*,*) "it =",it,"time = ", time;        
    
    ENDDO     !DO ispool = 1,nrec



    
!    CLOSE(in_unit);  ! close the combined file.....
    
    CALL nfs_close_file(nsf_self);
 
ENDDO

!------------------------------------------------------
!CLOSE THE FILE.....
!------------------------------------------------------
IF(i_out_form == i_out_form_SMS) THEN
    CALL close_sms_file (out_fileName,out_unit,ilrec);   
ELSE
!    CLOSE(out_unit);
ENDIF  
!np = np_global;
!ne = ne_global;

102 FORMAT(100E20.5);

END SUBROUTINE



!--------------------------------------------------------------
!extract the horizontal data from the combined results and 
! out put them to tecplot format.
!-------------------------------------------------------------
SUBROUTINE extract_vprof_netcdf_vel()
IMPLICIT NONE        ;

CHARACTER(LEN=200)  :: in_fileName,out_fileName,out_fileName_base,cday,crank,tec_file,tec_file_prx;
CHARACTER(LEN=200)  :: outFileName;
INTEGER             :: iday,irank,ispool,tec_appended,info2;
INTEGER             :: in_unit = 200;
INTEGER             :: out_unit= 201;
CHARACTER(LEN=200)  :: varName,cit;
INTEGER             :: out_it,istat;
!INTEGER             :: tec_form  = 1 ! =1 the resutls are in separate files (one time step per file).
!                                     ! =2 the results are in one file 
INTEGER             :: ilrec,i,j ;
INTEGER             :: data_type
INTEGER             :: ndata(1:3);
!---------------------------------------------------------
!read the parameters from file combine_output.in
!---------------------------------------------------------
!CALL get_horz_parm();

!---------------------------------------------------------
!get the specified layer.
!---------------------------------------------------------
!------set the tecplot file name-------------


!varName   = file63(1:lfile63-3);
!----------------------------------------------------------
!READ data for each day step.
!----------------------------------------------------------
out_it             = 0;    
out_fileName_base  = TRIM(outFolder)//TRIM(file63(1:lfile63-3));

n_vprof_node  = nout_node_list;

DO iday = ibgn,iend

!compose the file Name
    WRITE(cday,'(i12)')iday;
    cday  = TRIM(ADJUSTL(cday));

    in_fileName = TRIM(inFolder)//TRIM(cday)//'.nc';

!!-----------------------------------------------------------
!!Read the combined head file.
!!----------------------------------------------------------- 
!    CALL in_combine_head(in_fileName,in_unit);
    CALL nfs_open_file(nsf_self,TRIM(in_fileName),NF90_NOWRITE) ;
    
    CALL nsf_selfe_acquire_dims();

!----------------------------------------------------------- 
!Allocate the memory..........................................
!----------------------------------------------------------- 
    ivs = 2;
    CALL alloc_memory();

!----------------------------------------------------------- 
!get the geometry information...
!----------------------------------------------------------- 
    CALL nsf_selfe_in_grid();
       

    IF(iday == ibgn ) THEN
    
        IF(i_out_form == i_out_form_SMS)THEN
        
            out_fileName  = TRIM(out_fileName)//'_horz_sms.dat';
            CALL write_sms_head(out_fileName,out_unit,ilrec);
        ENDIF
        
    ENDIF        

   
    DO ispool = 1,nrec


!-----------------------------------------------------------    
!read simulation time 
!-----------------------------------------------------------    
        CALL nsf_self_netcdf_time_time('time',time,ispool);  
         
!-----------------------------------------------------------    
!read the elevation from the file.
!-----------------------------------------------------------
        ndata(1) = np;    
        CALL nsf_self_netcdf_time_data_in('elev',eta2,ispool,ndata);


!-----------------------------------------------------------
!read the horizontal velocity 
!-----------------------------------------------------------
       ndata(1:2) = (/np,nvrt/);  
       CALL nsf_self_netcdf_time_data_in('u',outb(:,:,1),ispool,ndata);
       CALL nsf_self_netcdf_time_data_in('v',outb(:,:,2),ispool,ndata);       
!-----------------------------------------------------------
!compute the magnitude of the velocity
!-----------------------------------------------------------

          DO j = 1,nvrt
                 DO i =1,np
                    outb(i,j,1) = SQRT(outb(i,j,1)*outb(i,j,1) + outb(i,j,2)*outb(i,j,2));
                ENDDO
          ENDDO 
!       outb(:,:,1) = SQRT(outb(:,:,1)*outb(:,:,1) + outb(:,:,2)*outb(:,:,1));
!       outb(:,:,1) =  outb(:,:,2)
!-----------------------------------------------------------
!read the vertical velocity 
!-----------------------------------------------------------
       CALL nsf_self_netcdf_time_data_in('w',outb(:,:,2),ispool,ndata);
         
!-----------------------------------------------------------    
!process the data................
!-----------------------------------------------------------    
!--------update the z coordinate.    
        kbp = 1;
        CALL update_z;

!----------------------------------------------------------- 
!get the tecplot out put var name  Not including x and y
!----------------------------------------------------------- 
        IF(ivs ==1) THEN
            varName = ' " '//TRIM(file63(1:lfile63-3))//' " ';
        ELSEIF(ivs ==2 ) THEN
            varName = ' "horz_mag", "w" ';
        ENDIF 
                
!-----------------------------------------------------------    
! Get the output tecplot file Name;
!-----------------------------------------------------------    
        IF(tec_form  == 1) THEN
            
            out_it = out_it + 1;
            WRITE(cit,"(I12)") out_it;
            cit = TRIM(ADJUSTL(cit)); 
            
            tec_file = TRIM(outFolder)//TRIM(file63(1:lfile63-3))//"_"//TRIM(cit)//"prof.dat"; 
            tec_appended = 0              
               
         ELSE
            
            tec_file = TRIM(outFolder)//TRIM(file63(1:lfile63-3))//"_prof.dat";
            
            IF(iday == ibgn .AND. ispool ==1) THEN
                tec_appended = 0   
            ELSE
                tec_appended = 1   
            ENDIF               
          ENDIF        
!--------extract and proces the data----------------------------------
!        CALL get_tec_hgrid_horz(out_proc)       
!        xver(1:n_vprof_node)        = x(iout_node_list(1:n_vprof_node));
!        zver(1:nvrt,1:n_vprof_node) = z(1:nvrt,iout_node_list(1:n_vprof_node));
        CALL write_TEC_vGRID(tec_file,time,x(iout_node_list(1:n_vprof_node)) , &
                z(1:nvrt,iout_node_list(1:n_vprof_node)),               &
                outb(iout_node_list(1:n_vprof_node),1:nvrt,1:ivs) ,     &
                TRIM(varName),"xx",tec_appended) ;


     
     
    WRITE(*,*) "it =",it,"time = ", time;        
    
    ENDDO     !DO ispool = 1,nrec



    
!    CLOSE(in_unit);  ! close the combined file.....
    
    CALL nfs_close_file(nsf_self);
 
ENDDO

!------------------------------------------------------
!CLOSE THE FILE.....
!------------------------------------------------------
IF(i_out_form == i_out_form_SMS) THEN
    CALL close_sms_file (out_fileName,out_unit,ilrec);   
ELSE
!    CLOSE(out_unit);
ENDIF  
!np = np_global;
!ne = ne_global;

102 FORMAT(100E20.5);

END SUBROUTINE
!--------------------------------------------------------
!allocate the memory for the variables.
!--------------------------------------------------------
SUBROUTINE alloc_memory()

INTEGER          :: istat;
INTEGER          :: xx_1,yy_1,yy_2,data_1,data_2,data_3;



!-------------------------------------------------
!ALLOCATE memory
!-------------------------------------------------
IF(ALLOCATED(ztot))  DEALLOCATE(ztot);
IF(ALLOCATED(sigma)) DEALLOCATE(sigma);
ALLOCATE(ztot(nvrt),sigma(nvrt),stat=istat);
if(istat/=0) stop 'Allocation error: ztot,sigma';


!-------------------------------------------------
!ALLOCATE memory
!-------------------------------------------------
IF(ALLOCATED(x))     DEALLOCATE(x);
IF(ALLOCATED(y))     DEALLOCATE(y);
IF(ALLOCATED(dp))    DEALLOCATE(dp);
IF(ALLOCATED(kbp))   DEALLOCATE(kbp);
!IF(ALLOCATED(kbp))   DEALLOCATE(kbp);
ALLOCATE(x(np_global),y(np_global),dp(np_global),kbp(np_global),stat=istat);
if(istat/=0) stop 'Allocation error: x,y,dp,kbp';
!
!ALLOCATE(kbp00(np_global),stat=istat);
!if(istat/=0) stop 'Allocation error: kbp00';

IF(ALLOCATED(nm)) DEALLOCATE(nm); 
ALLOCATE(nm(ne_global,3),stat=istat);
if(istat/=0) stop 'Allocation error: mn';

    

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














IF(ALLOCATED(hmod))      DEALLOCATE(hmod);
IF(ALLOCATED(cs))        DEALLOCATE(cs);
IF(ALLOCATED(kfp))       DEALLOCATE(kfp);
IF(ALLOCATED(htot))      DEALLOCATE(htot);
IF(ALLOCATED(z))         DEALLOCATE(z);
IF(ALLOCATED(delz))      DEALLOCATE(delz);           
IF(ALLOCATED(out_proc))  DEALLOCATE(out_proc);

ALLOCATE(hmod(np_global),cs(nvrt),kfp(np_global),htot(np_global),stat=istat)
if(istat/=0) stop 'Allocation error: hmod,cs,kfp,htot';

ALLOCATE(z(nvrt,np_global),delz(nvrt,np_global),stat=istat)
if(istat/=0) stop 'Allocation error: z,delz';
           
ALLOCATE(out_proc(np_global,ivs),stat=istat)
if(istat/=0) stop 'Allocation error: outb'; 


!--------------------------------------------------------
!allocate the output memory
!--------------------------------------------------------
IF(ALLOCATED(xx_out))DEALLOCATE(xx_out);
IF(ALLOCATED(yy_out))DEALLOCATE(yy_out);
IF(ALLOCATED(data_out))DEALLOCATE(data_out);

IF(i_mode == i_mode_horz_comb .OR. &
   i_mode == i_mode_horz_netcdf ) THEN
   
   xx_1 = np_global;
   yy_1 = np_global;
   yy_1 = 1;
   
   data_1 = np_global;
   data_2 = 1;
   data_3 = ivs; 
   
ELSEIF(i_mode == i_mode_vprof_comb .OR. &
       i_mode == i_mode_vprof_netcdf ) THEN
       
   xx_1 = nout_node_list; 
        
   yy_1 = nout_node_list;
   yy_2 = nvrt;
   
   data_1 = nout_node_list;
   data_2 = nvrt;
   data_3 = ivs;

ELSEIF(i_mode == i_mode_pTime_comb .OR. &
       i_mode == i_mode_pTime_netcdf ) THEN
       

   xx_1   = nout_node_list; 
        
   yy_1   = nout_node_list;
   yy_2   = 1;
   
   data_1 = nout_node_list;
   data_2 = 1;
   data_3 = ivs;
          
ENDIF


ALLOCATE(xx_out(xx_1),yy_out(yy_1,yy_2),data_out(data_1,data_2,data_3),stat=istat);



END SUBROUTINE


!--------------------------------------------------------
!get the input parameters 
!--------------------------------------------------------
SUBROUTINE  get_in_parm()

!--------------------------------------------------------
!Local variables.
!--------------------------------------------------------

!--------------------------------------------------------
!--------------------------------------------------------
open(10,file='combine_output.in',status='old')
READ(10,*) i_mode      ;   !i_mode: =0 combine the original selfe output, =1: output combined horz ret, =2: output vertical prof rst, =3: output point time series.
read(10,*) file63;   !file name to be combined.
read(10,*) ibgn,iend   ;   !begin and end day
read(10,*) i_out_form  ;   !output file form : SMS, TECPLOT, TXT;
read(10,*) tec_form    ;   !tec_form   =1:in separate TEC file, =2: in single TEC file; 
read(10,*) i_vprof_mode;   !i_vprof_mode =1: exttract the scalar, =2: extract the velocity (horz and vertical)
read(10,*) inFolder    ;   !input data folder.
read(10,*) outFolder   ;   !output data foler.
read(10,*) ilev        ;   !-2: vertical average; -1: BOTTOM_LAYER; 0: surface layer; 1~nvrt: particluar lever;
read(10,*) mark_dry,dry;   !mark_dry =1: mark the cell, =0 don't mark the dry cell;  dry=: value;
close(10)

IF(i_mode == i_mode_pTime_comb .OR. &
   i_mode == i_mode_pTime_netcdf )THEN
    
    i_out_form = i_out_form_TXT  ;

ENDIF


WRITE(*,*)"inFolder:" ,TRIM(inFolder);
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
!get the input parameters 
!--------------------------------------------------------
SUBROUTINE  get_horz_parm()

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
read(10,*) ilev     ;     !-2: vertical average; -1: BOTTOM_LAYER; 0: surface layer; 1~nvrt: particluar lever;
read(10,*) mark_dry,dry;  !mark_dry =1: mark the cell, =0 don't mark the dry cell;  dry=: value;
READ(10,*) i_out_form;    !i_out_form =1:TECPLOT ; =2: SMS;
close(10)


!inFolder = TRIM(inFolder)//"/";
!outFolder =TRIM(outFolder)//"/";
WRITE(*,*)"inFolder:" ,TRIM(inFolder);
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



!------------------------------------------------------------
!comput the z coordinate and delz for each node at each level.
!  Z(nvrt,np), delz(nvrt-1,np), according the eta2;
! variables needed: hmod,cs,sigma,kfp,htot,kbp
! variables output: z(nvrt,np), delz(nrt-1,np)
!------------------------------------------------------------
SUBROUTINE update_z
IMPLICIT NONE
INTEGER             :: i,k,iformat,kin;

!------------------------------------------------------
!compute cs;
!------------------------------------------------------
do k=1,nvrt-kz+1
    cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
    &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
enddo !k=1,nvrt

!-------------------------------------------------------
! compute cs for each node, iformat == 5
!-------------------------------------------------------
do i=1,np_global

    kfp(i) =nvrt
    htot(i)=dp(i)+eta2(i)
    if(htot(i).le.h0) kfp(i)=0

    hmod(i)=min(dp(i),h_s)  ;  !compute hmod 

    if(dp(i)+eta2(i)<=h0) then !dry
!      idry_new(i)=1 
      if(dp(i)>=h_s) then
!       write(errmsg,*)'Deep depth dry:',i
!       call parallel_abort(errmsg)
      endif
      kbp(i)=0       ! for dry node, kbp(i) ==0;
    else                        !wet
!      idry_new(i)=0

!-----------------------S-levels
      do k=kz,nvrt
        kin=k-kz+1

        if(hmod(i)<=h_c) then
!          iback(i)=1
          z(k,i)=sigma(kin)*(hmod(i)+eta2(i))+eta2(i)
        else if(eta2(i)<=-h_c-(hmod(i)-h_c)*theta_f/s_con1) then !hmod(i)>h_c>=0
!          write(errmsg,*)'Pls choose a larger h_c (1):',eta2(i),h_c
!          call parallel_abort(errmsg)
        else
          z(k,i)=eta2(i)*(1+sigma(kin))+h_c*sigma(kin)+(hmod(i)-h_c)*cs(kin)
        endif
      enddo !k=kz,nvrt

!------------------------z-levels
      if(dp(i)<=h_s) then
        kbp(i)=kz
      else !bottom index 
!        if(imm==1.or.it==iths) then
!          kbp(i)=0 !flag
!          do k=1,kz-1
!            if(-dp(i)>=ztot(k).and.-dp(i)<ztot(k+1)) then
!              kbp(i)=k
!              exit
!            endif
!          enddo !k
!          if(kbp(i)==0) then
!            write(errmsg,*)'Cannot find a bottom level for node (3):',i
!
!            call parallel_abort(errmsg)
!          endif
!        endif !imm

        if(kbp(i)>=kz.or.kbp(i)<1) then
!          write(errmsg,*)'Impossible 92:',kbp(i),kz,i
!          call parallel_abort(errmsg)
        endif
        
        z(kbp(i),i)=-dp(i)
        
        do k=kbp(i)+1,kz-1
          z(k,i)=ztot(k)
        enddo !k
        
      endif

      do k=kbp(i)+1,nvrt     !----------------check------------------------
        if(z(k,i)-z(k-1,i)<=0) then
!          write(errmsg,*)'Inverted z-levels at:',i,k,z(k,i)-z(k-1,i),eta2(i),hmod(i)
!          call parallel_abort(errmsg)
        endif
      enddo !k
      

       DO k=max(1,kbp(i)),nvrt-1

           delz(k,i)=z(k+1,i)-z(k,i)

       ENDDO !wet ot dry   
          
    endif !wet ot dry
enddo !i=1,npa


END SUBROUTINE



!-----------------------------------------------------
!get the tecplot horiztonal data.
!ilev    : specify the level of the data,
!mark_dry: specify mark the dry point ;
!output data : tec_out;
!-----------------------------------------------------
SUBROUTINE get_tec_hgrid_horz(out_data)
!-----------------------------------------------------
!arguments.
!-----------------------------------------------------
REAL(rkind)   ,DIMENSION(:,:),INTENT(OUT)::out_data;
!-----------------------------------------------------
!Local variables.
!----------------------------------------------------- 
INTEGER               :: i,k,m;
REAL(rkind)           :: summ(ivs);

IF(i23d == 2) THEN
   ilev = SELFE_BOTTOM_LAYER ; ! non 3D data equal to bottom layer.
ENDIF

DO i=1, np_global !number of node
            	  	  
  IF(kfp(i) == 0 ) THEN  !dry node;
  
    IF( mark_dry == 1) THEN
           out_data(i,:)    = dry;
    ELSE
        DO m=1,ivs
           out_data(i,m) = outb(i,1,m);  !only assign the bottom data.
        ENDDO
    ENDIF
  
  ELSE  ! wet node.
    
   IF(ilev == SELFE_VERTICAL_AVERAGE ) THEN     !vertical averaging
   
        DO m=1,ivs
           summ(m)=0.0
        ENDDO
       
        DO k=max(1,kbp(i)),kfp(i)-1   !if not dry kfp == nvrt;
           
             DO m=1,ivs
               summ(m)=summ(m)+(outb(i,k,m)+outb(i,k+1,m))/2.0d0*delz(k,i)
             ENDDO !m
           
        ENDDO !k
           
        DO m=1,ivs
 !             tec_out(i,m)=summ(m)/(z(kfp(i),i)-z(kbp(i),i))
              out_data(i,m)=summ(m)/(z(kfp(i),i)-z(max(1,kbp(i)),i))
        ENDDO
         
 
       
!   elseif (ilev.eq.-1) then !bottom
    ELSEIF( ilev == SELFE_BOTTOM_LAYER) THEN ! bottom layer.
    
        DO m=1,ivs
           out_data(i,m) = outb(i,max(1,kbp(i)),m); 
        ENDDO             


        
    ELSEIF( ilev == SELFE_SURFACE_LAYER) THEN ! SURFACE layer.

        DO m=1,ivs
           out_data(i,m) = outb(i,max(1,kfp(i)),m); 
        ENDDO
        
    ELSE  ! sepcified layer   :ilev
    
        IF(ilev > kfp(i) ) THEN   ! for air.
        
            DO m=1,ivs
               out_data(i,m) = air;
            ENDDO
            
        ELSEIF(ilev < kbp(i) ) THEN  ! for soil.
        
            DO m=1,ivs
               out_data(i,m) = soil;
            ENDDO   
               
        ELSE                      ! at particular level. 
         
            DO m=1,ivs  
                out_data(i,m) = outb(i,ilev,m);
            ENDDO
        
        ENDIF 
         
    ENDIF
    
ENDIF
                
          	  

ENDDO !DO i=1, np_global !number of node

END SUBROUTINE



SUBROUTINE read_node_list()

INTEGER ::  i,info2,istat;
INTEGER ::  unit =200;

OPEN (UNIT = unit,FILE="node_list.dat",ACTION="READ",IOSTAT=info2);   

READ(unit,*) nout_node_list;

IF(nout_node_list <=0) THEN 
    RETURN;
ENDIF 

IF(ALLOCATED(iout_node_list)) DEALLOCATE(iout_node_list);
ALLOCATE(iout_node_list(nout_node_list),stat=istat);
if(istat/=0) stop 'Allocation error: iout_node_list';

DO i =1,nout_node_list
    READ(unit,*) iout_node_list(i);
ENDDO
 
CLOSE(unit);

END SUBROUTINE

END MODULE