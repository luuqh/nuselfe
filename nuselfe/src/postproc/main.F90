PROGRAM main
USE combine_output7
USE SELFE_TEC_IO
USE post_proc_mod
IMPLICIT NONE;

CHARACTER(LEN=200)  :: fileName,cday,crank,tec_file;
CHARACTER(LEN=200)  :: outFileName;
INTEGER             :: iday,irank,ispool,tec_appended;
INTEGER             :: unit = 200;
CHARACTER(LEN=200)  :: varName;

!-------------------------------------------------------
!NOTE: For netcdf format theta_b, theta_f h_s mus specified.
!-----------------------------------------------------

!theta_b = 1.;
!theta_f = 1.e-4;
!h_s     = 150;
CALL get_in_parm();

!CALL extract_horz_netcdf();
IF( i_mode == i_mode_combine) THEN
    
    CALL combine_results();
    
ELSEIF(i_mode == i_mode_horz_comb) THEN

    CALL extract_horz_data();

ELSEIF(i_mode == i_mode_vprof_comb) THEN

    IF(i_vprof_mode == i_vprof_scalar ) THEN
        CALL extract_vprof_data();
    ELSEIF(i_vprof_mode ==i_vprof_velocity ) THEN

        CALL extract_vprof_vel();
    ENDIF

ELSEIF(i_mode == i_mode_pTime_comb) THEN

!---------------------------------------------------------
!get the vertical output nodelist.
!---------------------------------------------------------
    CALL read_node_list();

    CALL extract_horz_data();
ELSEIF(i_mode == i_mode_horz_netcdf) THEN
    
    CALL extract_horz_netcdf();
    

ELSEIF(i_mode == i_mode_vprof_netcdf) THEN

    CALL read_node_list();

    IF(i_vprof_mode == i_vprof_scalar ) THEN
        CALL extract_vprof_netcdf_data();
    ELSEIF(i_vprof_mode ==i_vprof_velocity ) THEN

        CALL extract_vprof_netcdf_vel();        
    ENDIF
    
ENDIF
!-----------------------------------------------------------
!
!-----------------------------------------------------------

!CALL extract_vprof_vel();
!CALL combine_results();
!CALL rst_hgrid_tec_out();
!CALL get_pos_parm();
!CALL get_lg_map()

!!----------------------------------------------------------
!!READ data for each day step.
!!----------------------------------------------------------
!tec_file  = TRIM(outFolder)//"tec.dat";
!varName   = file63(1:lfile63-3);
!
!
!DO iday = ibgn,iend
!
!    !compose the file Name
!    WRITE(cday,'(i12)')iday;
!    cday  = TRIM(ADJUSTL(cday))//'_';
!
!    fileName = TRIM(inFolder)//TRIM(cday)//TRIM(file63);
!    
!!-----------------------------------------------------------
!!Read the combined head file.
!!----------------------------------------------------------- 
!   
!    CALL in_combine_head(fileName,Unit)   ;
!    
!    DO ispool = 1,nrec
!    
!!-----------------------------------------------------------    
!! READ the time dependent data.
!!-----------------------------------------------------------    
!        CALL in_combine_time_data(Unit)
!
!!-----------------------------------------------------------    
!! Write the data to the results.
!!-----------------------------------------------------------    
!        IF(iday == ibgn .AND. ispool ==1) THEN
!            tec_appended = 0   
!        ELSE
!            tec_appended = 1   
!        ENDIF
!!SUBROUTINE write_TEC_hGRID_scalar(fileName,simtime,outdata,var_name,aux_data,append)        
!        CALL write_TEC_hGRID_scalar(tec_file,time,outb(:,1,1),TRIM(varName),TRIM(varName),tec_appended);
!
!
!        WRITE(*,*) "time = ", time,"it =",it;
!        
!    ENDDO
!    
!    CLOSE(Unit); 
!ENDDO
!
!
!np = np_global;
!ne = ne_global;
!CALL write_TEC_hGRID_scalar("test.dat",0.0d0,dp,"depth","depth",0);

END PROGRAM








!!---------------------------------------------------------------------
!!output the results to tecplot to visualize the data.
!!---------------------------------------------------------------------
!SUBROUTINE rst_hgrid_tec_out()
!USE combine_output7
!USE SELFE_TEC_IO
!IMPLICIT NONE;
!
!CHARACTER(LEN=200)  :: fileName,cday,crank,tec_file,tec_file_prx;
!CHARACTER(LEN=200)  :: outFileName;
!INTEGER             :: iday,irank,ispool,tec_appended;
!INTEGER             :: unit = 200;
!CHARACTER(LEN=200)  :: varName,cit;
!INTEGER             :: out_it;
!INTEGER             :: tec_form  = 1 ! =1 the resutls are in separate files (one time step per file).
!                                     ! =2 the results are in one file 
!
!!CALL combine_results();
!
!CALL get_pos_parm();
!!CALL get_lg_map()
!
!!----------------------------------------------------------
!!READ data for each day step.
!!----------------------------------------------------------
!tec_file_prx  = TRIM(outFolder)//"tec";
!!varName   = file63(1:lfile63-3);
!
!out_it    = 0;
!
!DO iday = ibgn,iend
!
!    !compose the file Name
!    WRITE(cday,'(i12)')iday;
!    cday  = TRIM(ADJUSTL(cday))//'_';
!
!    fileName = TRIM(inFolder)//TRIM(cday)//TRIM(file63);
!    
!!-----------------------------------------------------------
!!Read the combined head file.
!!----------------------------------------------------------- 
!    CALL in_combine_head(fileName,Unit);
!    
!   
!    DO ispool = 1,nrec
!    
!
!!-----------------------------------------------------------    
!! READ the time dependent data.
!!-----------------------------------------------------------    
!        CALL in_combine_time_data(Unit)
!        time = time / (24*60*60);  !-----------change the time to day
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
!            tec_file = TRIM(tec_file_prx)//"_"//TRIM(cit); 
!            tec_appended = 0              
!           
!        ELSE
!        
!            tec_file = tec_file_prx;
!            IF(iday == ibgn .AND. ispool ==1) THEN
!                tec_appended = 0   
!            ELSE
!                tec_appended = 1   
!            ENDIF               
!        
!        ENDIF
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
!        
!        CALL write_TEC_hGRID(tec_file,time,tec_out(:,1:ivs),TRIM(varName),"xx",tec_appended);
!
!        WRITE(*,*) "it =",it,"time = ", time;
!        
!    ENDDO     !DO ispool = 1,nrec
!    
!    CLOSE(Unit);  ! close the combined file.....
!ENDDO
!
!
!!np = np_global;
!!ne = ne_global;
!
!
!END SUBROUTINE 