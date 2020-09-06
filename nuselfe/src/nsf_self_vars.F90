!-----------------------------------------------------------------------------------
!c    COPYRIGHT 2010: prof. Pavel Tkalich(TMSI), Xu Haihua(TMSI), Dao My Ha (TMSI)
!c    1. This subroutine is developed by Xu Haihua(TMSI), Dao My Ha(TMSI) and supervised
!c       by prof. Pavel Tkalich(TMSI)
!     2. TMSI: Tropical Marine Science Institute, Singapore
!c    3. For any enquiry pls contact tmsxh@nus.edu.sg
!-----------------------------------------------------------------------------------

MODULE nsf_self_vars_mod
USE nsf_container_mod
USE elfe_glbl
USE elfe_msgp

IMPLICIT NONE;

integer                                 :: nhotout_recd = 0,info2;

CHARACTER(LEN=200)                      :: nc_file_name;  

INTEGER(intKind)                        :: nsf_selfe_dims_count       =10;
INTEGER(intKind)                        :: nsf_selfe_global_att_count =10;
INTEGER(intKind)                        :: nsf_selfe_output_count     =26; !only time dependent variables.
INTEGER(intKind)                        :: nsf_selfe_total_vars          ;
INTEGER(intKind)                        :: global_id;
TYPE(Tnsf_container)                    :: nsf_self;
REAL(realKind),DIMENSION(:),ALLOCATABLE :: temp_out(:),temp_out2(:);
REAL(realKind)                          :: nsf_selfe_time_start,nsf_selfe_time_end;

!--------------------------------------------------------------
!nsf_np        : used to store number of particles in each process.
!nsf_np_index  : the start index of nsf_ipgl for each proc;
!nsf_buffer    : buffer to store the output data.
!nsf_data      : output data after assemble.
!--------------------------------------------------------------
INTEGER(intKind),DIMENSION(:),ALLOCATABLE           :: nsf_np,nsf_np_index;
INTEGER(intKind),DIMENSION(:),ALLOCATABLE           :: nsf_ipgl;              
REAL(realKind)  ,DIMENSION(:,:),ALLOCATABLE         :: nsf_data3d; 
REAL(realKind)  ,DIMENSION(:)  ,ALLOCATABLE         :: nsf_data2d,nsf_buffer;
!-----------------------------------------------
!np_array  : store the number of point in each procs.
!index_base: the base index for each proc 
!   index_base=sum(np_array(1:myrank));
!ip_index      : the ip in the global array index. 
!nsf_out_ntime : output time steps.   
!-----------------------------------------------
INTEGER(intKind),DIMENSION(:),ALLOCATABLE :: np_array     !------------------------
INTEGER(intKind)                          :: index_base ,ip_index  !
INTEGER(intKind)                          :: nsf_out_ntime = 0;
INTEGER(intKind)                          :: temp_var_id,temp_var_id2;
!REAL(MK)                                 :: temp_var_val,temp_var_val2;

INTEGER(intKind),DIMENSION(50)            :: iof_netcdf;

!!!2. ID
INTEGER(intKind):: nsf_dim_node_id    ;
INTEGER(intKind):: nsf_dim_nele_id    ;
INTEGER(intKind):: nsf_dim_nbnd_id    ;
INTEGER(intKind):: nsf_dim_nface_id   ;
INTEGER(intKind):: nsf_dim_nbi_id     ;
INTEGER(intKind):: nsf_dim_nzlayer_id ;
INTEGER(intKind):: nsf_dim_nslayer_id ;
INTEGER(intKind):: nsf_dim_nlayer_id  ;
INTEGER(intKind):: nsf_dim_time_id    ;
INTEGER(intKind):: nsf_dim_two_id     ;

INTEGER(intKind),POINTER:: nsf_dim_node_val   ;
INTEGER(intKind),POINTER:: nsf_dim_nele_val   ;
INTEGER(intKind),POINTER:: nsf_dim_nbnd_val   ;
INTEGER(intKind),POINTER:: nsf_dim_nface_val  ;
INTEGER(intKind),POINTER:: nsf_dim_nbi_val    ;
INTEGER(intKind),POINTER:: nsf_dim_nzlayer_val;
INTEGER(intKind),POINTER:: nsf_dim_nslayer_val;
INTEGER(intKind),POINTER:: nsf_dim_nlayer_val ;
INTEGER(intKind),POINTER:: nsf_dim_time_val   ;

INTEGER(intKind)        :: nsf_self_ncid    ;
!-----------------------------------------------
!SELFE grid varialbe id.
!-----------------------------------------------
INTEGER(intKind),POINTER:: nsf_var_x_id;
INTEGER(intKind),POINTER:: nsf_var_y_id;
INTEGER(intKind),POINTER:: nsf_var_depth_id;

INTEGER(intKind):: nsf_var_time_id 
INTEGER(intKind):: nsf_var_elev_id 
INTEGER(intKind):: nsf_var_pres_id 
INTEGER(intKind):: nsf_var_airt_id 
INTEGER(intKind):: nsf_var_shum_id 
INTEGER(intKind):: nsf_var_srad_id 
INTEGER(intKind):: nsf_var_flsu_id 
INTEGER(intKind):: nsf_var_fllu_id  
INTEGER(intKind):: nsf_var_radu_id   
INTEGER(intKind):: nsf_var_radd_id   
INTEGER(intKind):: nsf_var_flux_id   
INTEGER(intKind):: nsf_var_evap_id  
INTEGER(intKind):: nsf_var_prcp_id 
INTEGER(intKind):: nsf_var_windx_id
INTEGER(intKind):: nsf_var_windy_id                                 
INTEGER(intKind):: nsf_var_wistx_id
INTEGER(intKind):: nsf_var_wisty_id
INTEGER(intKind):: nsf_var_dahvx_id
INTEGER(intKind):: nsf_var_dahvy_id             
INTEGER(intKind):: nsf_var_vert_id                  
INTEGER(intKind):: nsf_var_temp_id 
INTEGER(intKind):: nsf_var_salt_id 
INTEGER(intKind):: nsf_var_conc_id 
INTEGER(intKind):: nsf_var_tdff_id 
INTEGER(intKind):: nsf_var_vdff_id 
INTEGER(intKind):: nsf_var_kine_id 
INTEGER(intKind):: nsf_var_mixl_id 
INTEGER(intKind):: nsf_var_zcor_id 
INTEGER(intKind):: nsf_var_qnon_id 
INTEGER(intKind):: nsf_var_hvelx_id
INTEGER(intKind):: nsf_var_hvely_id     
INTEGER(intKind):: nsf_var_u_id;
INTEGER(intKind):: nsf_var_v_id;         

INTEGER(intKind)  ,DIMENSION(100):: nsf_var_output_ids;
CHARACTER(CharLen),DIMENSION(100):: nsf_var_output_names;

CHARACTER(CharLen),DIMENSION(100):: nsf_selfe_dim_names;
INTEGER(intKind)  ,DIMENSION(100):: nsf_selfe_dim_vals;


INTERFACE nsf_self_netcdf_time_data_in
    MODULE PROCEDURE nsf_self_netcdf_time_data_in_1D
    MODULE PROCEDURE nsf_self_netcdf_time_data_in_2D
    MODULE PROCEDURE nsf_self_netcdf_time_data_in_3D 
END INTERFACE            
  
CONTAINS



 
 
!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nfs_self_create_file
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : create the output netcdf file for SELF.
!                  
!                 
!                 
!  Input        : 
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
SUBROUTINE nfs_self_create_file(container,path,mode,noutput,iof)
!---------------------------------------------------------
!define the interface of subroutine ppInteraction
!---------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_container)                       :: container;
CHARACTER(*)                  ,INTENT(IN  ):: path;
INTEGER(intKind)              ,INTENT(IN  ):: mode;
INTEGER(intKind)              ,INTENT(IN  ):: noutput;
INTEGER(intKind) ,DIMENSION(:),INTENT(IN  ):: iof;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                 :: i,j,status;
INTEGER(intKind)                 :: varid;
INTEGER(intKind)                 :: n_vars,old_mode;
INTEGER(intKind)                 :: SELFE_OUT = 26;

nsf_selfe_output_count = noutput;

!---------------------------------------------------------------------------------------------------------------------
!1.estimate the total number of ouput varialbes in the self
!---------------------------------------------------------------------------------------------------------------------
!CALL nfs_self_est_nvars(iof,nsf_selfe_total_vars);

!---------------------------------------------------------------------------------------------------------------------
!2.init the container first when n_att,n_dim and n_var are known.
!---------------------------------------------------------------------------------------------------------------------
CALL init_nsf_container(container);
!CALL init_nsf_container(container,nsf_selfe_global_att_count, &
!                        nsf_selfe_dims_count,nsf_selfe_total_vars);

!---------------------------------------------------------------------------------------------------------------------
!init nsf_self
!---------------------------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------------------------
!3.init the global attribute, dimension and output varialbes.
!---------------------------------------------------------------------------------------------------------------------
!CALL nsf_selfe_init(noutput);
CALL nfs_self_set_global_atts(container);
CALL nfs_self_set_dims(container);
CALL nfs_self_set_vars(container,iof);

!---------------------------------------------------------------------------------------------------------------------
!4.all are done except create the file.
!---------------------------------------------------------------------------------------------------------------------
CALL nfs_create_file(container,path,mode);

!---------------------------------------------------------------------------------------------------------------------
!5.get the time dependent variable ids.
!---------------------------------------------------------------------------------------------------------------------
!DO i =1 ,nsf_selfe_output_count
!    IF( i == 26) THEN
!        CALL nsf_get_var_id(container,"u",nsf_var_u_id);
!        CALL nsf_get_var_id(container,"v",nsf_var_v_id);
!    ELSE
!        CALL nsf_get_var_id(container,TRIM(nsf_var_output_names(i)),nsf_var_output_ids(i));
!    ENDIF
!ENDDO

END SUBROUTINE


!!------------------------------------------------------------------------------------------------------------------
!!  Subroutine   :                   nfs_self_write_var2file_2d_scalar
!!------------------------------------------------------------------------------------------------------------------
!!
!!  Purpose      : write the 2 dimensional SELFE scalar varialbe to netcdf file, the variable is scalar.
!!                  The subroutine collect the data then write to the netcdf file by the rank 0 processor.
!!                  
!!                 
!!                 
!!  Input        : ncid    : the id of the netcdf file.
!!                 varname : the variable name  
!!                 var     : array contain the variable.
!!                 iterm   : the record number in the netcdf file.
!!                  
!!                 
!!   
!!                 
!!
!!  Input/output : 
!!
!!  Output       : 
!!			      
!!
!!  Routines     : 
!!                 
!!
!!  Remarks      :
!!
!!-----------------------------------------------------------------    
!!  References   :
!!
!!  Revisions    :
!!------------------------------------------------------------------------------------------------------------------
!SUBROUTINE nfs_self_write_var2file_2d_scalar(ncid,varname,var,iterm)
!!---------------------------------------------------------
!!define the interface of subroutine ppInteraction
!!---------------------------------------------------------
!IMPLICIT NONE
!!---------------------------------------------------------------------------------------------------------------------
!!  Modules 
!!---------------------------------------------------------------------------------------------------------------------
!!---------------------------------------------------------------------------------------------------------------------
!!  Arguments     
!!---------------------------------------------------------------------------------------------------------------------
!INTEGER(intKind)               ,INTENT(IN  )::ncid;
!CHARACTER(LEN=*)               ,INTENT(IN  )::varname;
!REAL(kind=rkind),DIMENSION(:)  ,INTENT(IN  )::var;
!INTEGER(intKind)               ,INTENT(IN  )::var_id,ierr,status;
!
!
!
!REAL(kind=rkind)(:,:)               :: data2d;
!
!!Allocate the memory before collect the data
!ALLOCATE(data2d(np_global),stat =ierr );
!IF(ierr /=0) STOP "Failed to allocate memory data2d in subroutine nfs_self_write_var2file_2d_scalar";
!
!!Collect the data......................
!CALL nfs_collect_rdata_2d(var,data2d);
!
!
!!---------rank 0 output the data to the netcdf file.
!IF( myrank ==0) THEN
!!--------------get the var_id; 
!    status = nf90_inq_varid(ncid,TRIM(varName), var_id);
!    IF(status /=0) THEN
!        WRITE(*,*) "Failed to get the variable Id in netcdf file, variable name = ",TRIM(varName);
!    ENDIF
!                                                                   
!    IF(nsf_float == NF90_FLOAT) THEN     
!      status = nf90_put_var(ncid,var_id,REAL(data2d,4), &
!                       start=(/1,iterm/),count=(/np_global,1/));
!    ELSEIF(nsf_float == NF90_DOUBLE) THEN 
!      status = nf90_put_var(ncid,var_id,REAL(data2d,8), &
!                       start=(/1,iterm/),count=(/np_global,1/));
!
!    ENDIF
!    
!    CALL  check_error(status)
!ENDIF !IF( myrank ==0) THEN  
!
!!----Deallocate the memory after the subroutine...............
!DEALLOCATE(data2d);
!
!END SUBROUTINE
!
!
!!------------------------------------------------------------------------------------------------------------------
!!  Subroutine   :                   nfs_self_write_var2file_3d_scalar
!!------------------------------------------------------------------------------------------------------------------
!!
!!  Purpose      : write the 3 dimensional SELFE scalar varialbe to netcdf file, the variable is scalar.
!!                  The subroutine collect the data then write to the netcdf file by the rank 0 processor.
!!                  
!!                 
!!                 
!!  Input        : ncid    : the id of the netcdf file.
!!                 varname : the variable name  
!!                 var     : array contain the variable. the variable need in the order of (nvrt,np);
!!                 iterm   : the record number in the netcdf file.
!!                  
!!                 
!!   
!!                 
!!
!!  Input/output : 
!!
!!  Output       : 
!!			      
!!
!!  Routines     : 
!!                 
!!
!!  Remarks      :
!!
!!-----------------------------------------------------------------    
!!  References   :
!!
!!  Revisions    :
!!------------------------------------------------------------------------------------------------------------------
!SUBROUTINE nfs_self_write_var2file_3d_scalar(ncid,varname,var,iterm)
!!---------------------------------------------------------
!!define the interface of subroutine ppInteraction
!!---------------------------------------------------------
!IMPLICIT NONE
!!---------------------------------------------------------------------------------------------------------------------
!!  Modules 
!!---------------------------------------------------------------------------------------------------------------------
!!---------------------------------------------------------------------------------------------------------------------
!!  Arguments     
!!---------------------------------------------------------------------------------------------------------------------
!INTEGER(intKind)               ,INTENT(IN  )::ncid;
!CHARACTER(LEN=*)               ,INTENT(IN  )::varname;
!REAL(kind=rkind),DIMENSION(:,:),INTENT(IN  )::var;
!INTEGER(intKind)               ,INTENT(IN  )::var_id,ierr,status;
!
!
!
!REAL(kind=rkind)(:,:)               :: data3d;
!
!
!!Allocate the memory before collect the data
!ALLOCATE(data3d(nvrt,np_global),stat =ierr );
!IF(ierr /=0) STOP "Failed to allocate memory data3d in subroutine nfs_self_write_var2file_3d_scalar";
!
!!Collect the data......................
!CALL nfs_collect_rdata_3d(var,2,data3d);
!
!
!!---------rank 0 output the data to the netcdf file.
!IF( myrank ==0) THEN
!!--------------get the var_id; 
!    status = nf90_inq_varid(ncid,TRIM(varName), var_id);
!    IF(status /=0) THEN
!        WRITE(*,*) "Failed to get the variable Id in netcdf file, variable name = ",TRIM(varName);
!    ENDIF
!                                                                   
!    IF(nsf_float == NF90_FLOAT) THEN     
!      status = nf90_put_var(ncid,var_id,REAL(data3d,4), &
!                       start=(/1,iterm/),count=(/np_global,1/));
!    ELSEIF(nsf_float == NF90_DOUBLE) THEN 
!      status = nf90_put_var(ncid,var_id,REAL(data3d,8), &
!                       start=(/1,iterm/),count=(/np_global,1/));
!
!    ENDIF
!    
!    CALL  check_error(status)
!ENDIF !IF( myrank ==0) THEN  
!
!!----Deallocate the memory after the subroutine...............
!DEALLOCATE(data3d);
!
!END SUBROUTINE



!---------------------------------------------------------
!get time dependent variables ids. 
!---------------------------------------------------------
SUBROUTINE get_nfs_self_vars_id(container)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_container)                       :: container;
!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                 :: i;


!---------------------------------------------------------------------------------------------------------------------
!5.get the time dependent variable ids.
!---------------------------------------------------------------------------------------------------------------------
DO i =1 ,nsf_selfe_output_count

    IF(iof(i) ==0) CYCLE;
    
    IF( i == 26) THEN
        CALL nsf_get_var_id(container,"u",nsf_var_u_id);
        CALL nsf_get_var_id(container,"v",nsf_var_v_id);
    ELSE
        CALL nsf_get_var_id(container,TRIM(nsf_var_output_names(i)),nsf_var_output_ids(i));
    ENDIF
!   WRITE(*,*) i,TRIM(nsf_var_output_names(i));
ENDDO

CALL nsf_get_var_id(container,"time",nsf_var_time_id);

END SUBROUTINE


!---------------------------------------------------------------------------------------------------------------------
!init the varialbes name id for netcdf ile
!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE nsf_selfe_nectcdf_out_init(noutput)

IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                ::i,noutput,ierr;






nsf_out_ntime              =0;



            
IF(ALLOCATED(temp_out))  DEALLOCATE (temp_out);
IF(ALLOCATED(temp_out2)) DEALLOCATE (temp_out2);

ALLOCATE(temp_out(nvrt),temp_out2(nvrt),stat =ierr );
IF(ierr /=0) STOP "Failed to allocate memory temp_out & temp_out2 "

IF(ALLOCATED(nsf_np))DEALLOCATE(nsf_np);
IF(ALLOCATED(nsf_np_index))DEALLOCATE(nsf_np_index);
ALLOCATE(nsf_np(nproc),nsf_np_index(nproc+1),stat =ierr );
IF(ierr /=0) STOP "Failed to allocate memory nsf_np & nsf_np_index "

!---------------------------------------------------------------------------------------------------------------------
!get the number of nodes on each proc.
!---------------------------------------------------------------------------------------------------------------------
!int MPI_Allgather ( void *sendbuf, int sendcount, MPI_Datatype sendtype,
!                    void *recvbuf, int recvcount, MPI_Datatype recvtype, 
!                   MPI_Comm comm )
call mpi_allgather(np,1,itype, &
                   nsf_np,1,itype,comm,ierr);
                   
!write(myrank+10000,*)"nsf_np:", nsf_np; 

                   
nsf_np_index(1) = 1;                   
DO i =1,nproc  !get the start index 
    nsf_np_index(i+1) = nsf_np_index(i) + nsf_np(i)
ENDDO

!for test 
!write(myrank+10000,*)"nsf_np_index:", nsf_np_index;

!---------------------------------------------------------------------------------------------------------------------
!ALLOCATE THE MEMORY FOR nsf_ipgl, buff transfer data out put
!---------------------------------------------------------------------------------------------------------------------
IF(myrank ==0) THEN

    IF(ALLOCATED(nsf_ipgl))    DEALLOCATE(nsf_ipgl);
    IF(ALLOCATED(nsf_data2d))  DEALLOCATE(nsf_data2d);
    IF(ALLOCATED(nsf_data3d))  DEALLOCATE(nsf_data3d);    
    IF(ALLOCATED(nsf_buffer))  DEALLOCATE(nsf_buffer);

    ALLOCATE(nsf_ipgl(nsf_np_index(nproc+1)-1),stat =ierr );
    IF(ierr /=0) STOP "Failed to allocate memory nsf_ipgl "
    
    ALLOCATE(nsf_data2d(np_global),stat =ierr );
    IF(ierr /=0) STOP "Failed to allocate memory nsf_data2d ";

    ALLOCATE(nsf_data3d(nvrt,np_global),stat =ierr );
    IF(ierr /=0) STOP "Failed to allocate memory nsf_data2d ";
    
    ALLOCATE(nsf_buffer(nvrt*(nsf_np_index(nproc+1)-1)),stat =ierr );                
    IF(ierr /=0) STOP "Failed to allocate memory nsf_buffer ";

ENDIF

!---------------------------------------------------------------------------------------------------------------------
!collect the local to global index to each proc.
!---------------------------------------------------------------------------------------------------------------------
!int MPI_Gatherv ( void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
!                  void *recvbuf, int *recvcnts, int *displs, 
!                 MPI_Datatype recvtype, 
!                  int root, MPI_Comm comm )
!nsf_ipgl = 1;
call MPI_Gatherv(iplg, np ,itype, &
                 nsf_ipgl,nsf_np,nsf_np_index-1, &
                 itype,0,comm,ierr)


!IF(myrank ==0) THEN
!
!ENDIF
!IF(myrank ==0) THEN                 
!    WRITE(100000,*)   nsf_ipgl;
!ENDIF              
!IF(myrank ==0) THEN
!    write(1001,*)nproc,nsf_np_index(nproc+1);
!    DO i =1,nsf_np_index(nproc+1)-1
!        WRITE(1001,*)nsf_ipgl(i);
!    ENDDO
!ENDIF
!                 
!stop;



!IF(myrank ==0)THEN
!    WRITE(*,*)nsf_nnp;
!ENDIF


END SUBROUTINE



!---------------------------------------------------------------------------------------------------------------------
!free the memory allocated for netcdf output
!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE nsf_selfe_nectcdf_out_free(noutput)

IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                ::i,noutput,ierr;



IF(ALLOCATED(temp_out))  DEALLOCATE (temp_out);
IF(ALLOCATED(temp_out2)) DEALLOCATE (temp_out2);

IF(ALLOCATED(nsf_np))      DEALLOCATE(nsf_np);
IF(ALLOCATED(nsf_np_index))DEALLOCATE(nsf_np_index);



!---------------------------------------------------------------------------------------------------------------------
!ALLOCATE THE MEMORY FOR nsf_ipgl, buff transfer data out put
!---------------------------------------------------------------------------------------------------------------------
IF(myrank ==0) THEN

    IF(ALLOCATED(nsf_ipgl))    DEALLOCATE(nsf_ipgl);
    IF(ALLOCATED(nsf_data2d))  DEALLOCATE(nsf_data2d);
    IF(ALLOCATED(nsf_data3d))  DEALLOCATE(nsf_data3d);    
    IF(ALLOCATED(nsf_buffer))  DEALLOCATE(nsf_buffer);


ENDIF


END SUBROUTINE


!---------------------------------------------------------------------------------------------------------------------
! Collect the 2d data,
!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE  nfs_collect_rdata_2d(indata,outdata)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
REAL(rkind),DIMENSION(:),INTENT(IN   )  :: indata ;
REAL(rkind),DIMENSION(:),INTENT(  OUT)  :: outdata;

!---------------------------------------------------------------------------------------------------------------------
!Local variable.
!---------------------------------------------------------------------------------------------------------------------
INTEGER                            ::i,j;
!---------------------------------------------------------------------------------------------------------------------
!1. gather the data.
!---------------------------------------------------------------------------------------------------------------------
call MPI_Gatherv(indata    , nsf_np(myrank+1),rtype, &
                 nsf_buffer, nsf_np,nsf_np_index-1, &
                 rtype,0,comm,ierr)

!---------------------------------------------------------------------------------------------------------------------
!for rank 0: assemble the proces according nsf_iplg
!---------------------------------------------------------------------------------------------------------------------
IF(myrank ==0) THEN

    DO i = 1,nproc
        outdata(nsf_ipgl(nsf_np_index(i):nsf_np_index(i+1)-1)) = &
            nsf_buffer(nsf_np_index(i):nsf_np_index(i+1)-1);
    ENDDO

ENDIF


END SUBROUTINE

!---------------------------------------------------------------------------------------------------------------------
! Collect the 3d data,
!  order : =1 the data is stored in (np,vrt) order
!          =2 the data is stored in (vrt,np) order
!  the out put data is stored in (vrt,np) order consistent with netcdf format
!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE  nfs_collect_rdata_3d(indata,order,outdata)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
REAL(rkind),DIMENSION(:,:),INTENT(IN   )  :: indata ;
INTEGER                                   :: order
REAL(rkind),DIMENSION(:,:),INTENT(  OUT)  :: outdata;

!---------------------------------------------------------------------------------------------------------------------
!Local variable.
!---------------------------------------------------------------------------------------------------------------------
INTEGER                            ::iProc,j,ip,id;
!---------------------------------------------------------------------------------------------------------------------
!1. gather the data.
!---------------------------------------------------------------------------------------------------------------------
call MPI_Gatherv(indata    , nsf_np(myrank+1)*nvrt,rtype, &
                 nsf_buffer, nsf_np*nvrt,(nsf_np_index-1)*nvrt, &
                 rtype,0,comm,ierr)

!---------------------------------------------------------------------------------------------------------------------
!for rank 0: assemble the proces according nsf_iplg
!---------------------------------------------------------------------------------------------------------------------
IF( myrank ==0 .AND. order ==1) THEN

        DO iProc = 1,nproc
        
            DO ip = nsf_np_index(iProc),nsf_np_index(iProc+1)-1 
                
               id   =nsf_ipgl(ip);
    !assemble according to (vrt,np) order              
               outdata(:,id) = nsf_buffer(ip:(ip+nsf_np(iProc)*nvrt):nsf_np(iProc)) 
            
            ENDDO !DO iProc = 1,nproc

        ENDDO
        
ELSEIF( myrank ==0 .AND. order ==2 ) THEN
!---------------------------------------------------------------------------------------------------------------------
!for rank 0: assemble the proces according nsf_iplg
!---------------------------------------------------------------------------------------------------------------------

            DO iProc = 1,nproc
            
                DO ip = nsf_np_index(iProc),nsf_np_index(iProc+1)-1 
                    
                   id   =nsf_ipgl(ip);
        !assemble according to (vrt,np) order              
                   outdata(1:nvrt,id) = nsf_buffer((ip-1)*nvrt+1:(ip-1)*nvrt+1+nvrt) 
                
                ENDDO

            ENDDO

ENDIF !IF( order ==1) THEN

END SUBROUTINE



SUBROUTINE  nfs_self_set_global_atts(nsf_self)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_container),INTENT(INOUT) :: nsf_self;
!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
CHARACTER(CharLen)       :: version,ctime;
INTEGER(intKind)         :: status ,date_time(8);


version(1:80)  = nf90_inq_libvers();
!CALL check_error(status);
version = "NETCDF "//TRIM(version);


CALL DATE_AND_TIME(values=date_time);
WRITE(ctime,101)date_time(1:3),date_time(5:7);


nsf_selfe_global_att_count =10

!---------------------------------------------------------------------------------------------------------------------
!only need to set the dim name and values.
!---------------------------------------------------------------------------------------------------------------------
!!1.Name
nsf_self%atts(1)%name      ="conventions";
nsf_self%atts(2)%name      ="grid_type";
nsf_self%atts(3)%name      ="model";
nsf_self%atts(4)%name      ="title";
nsf_self%atts(5)%name      ="comment";
nsf_self%atts(6)%name      ="source"
nsf_self%atts(7)%name      ="institution";
nsf_self%atts(8)%name      ="history";
nsf_self%atts(9)%name      ="references" ;
nsf_self%atts(10)%name     ="creation_date";

CALL nsf_set_attribute(nsf_self%atts(1),version          )
CALL nsf_set_attribute(nsf_self%atts(2),"Triangular"     )
CALL nsf_set_attribute(nsf_self%atts(3),"SELFE"          )
CALL nsf_set_attribute(nsf_self%atts(4),"CER"            )
CALL nsf_set_attribute(nsf_self%atts(5),"testing"        )
CALL nsf_set_attribute(nsf_self%atts(6),"Fortran script" )
CALL nsf_set_attribute(nsf_self%atts(7),"PORL/TMSI/NUS"  )
CALL nsf_set_attribute(nsf_self%atts(8),"original"       )
CALL nsf_set_attribute(nsf_self%atts(9),"XXXXXX"         )
CALL nsf_set_attribute(nsf_self%atts(10),ctime           )


nsf_self%atts_count = nsf_selfe_global_att_count;


101 FORMAT(I4,"-",I2.2,"-",I2.2, " ",I2.2, ":", I2.2, ":",I2.2)
END SUBROUTINE


SUBROUTINE  nfs_self_set_dims(nsf_self)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_container),INTENT(INOUT) :: nsf_self;
!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                :: i;


nsf_selfe_dims_count         =10 ;
nsf_selfe_dim_names(1:nsf_selfe_dims_count) =&
["node"      ,&
"nele"      ,&
"nbnd"      ,&
"nface"     ,&
"nbi"       ,&
"nzlayer"   ,&
"nslayer"   ,&
"nlayer"    ,&
"time"      ,&
"two"        ]  
      
nsf_selfe_dim_vals(1:nsf_selfe_dims_count) =[np_global,ne_global, max(mnond_global,mnlnd_global) ,& 
                                       3,nope_global + nland_global,kz ,nsig,nvrt,nf90_unlimited,2];            
nsf_dim_node_id     =1;   
nsf_dim_nele_id     =2;
nsf_dim_nbnd_id     =3;
nsf_dim_nface_id    =4;
nsf_dim_nbi_id      =5;
nsf_dim_nzlayer_id  =6;
nsf_dim_nslayer_id  =7;
nsf_dim_nlayer_id   =8;
nsf_dim_time_id     =9;
nsf_dim_two_id      =10;           

nsf_self%dims_count = nsf_selfe_dims_count;

DO i =1,nsf_self%dims_count
   nsf_self%dims(i)%name =nsf_selfe_dim_names(i);
   nsf_self%dims(i)%val  =nsf_selfe_dim_vals(i);     
ENDDO

END SUBROUTINE


SUBROUTINE  nfs_self_set_vars(nsf_self,iof)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_container)         ,INTENT(INOUT) :: nsf_self;
INTEGER(intKind),DIMENSION(:),INTENT(IN   ) :: iof;
!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                :: i,dimCount,attCount;
INTEGER(intKind)                :: iof_index;
TYPE(Tnsf_attribute)            :: long_name_att,std_name_att,unit_att,fill_val_int_att,fill_val_float_att;
TYPE(Tnsf_attribute)            :: positive_att,base_date_att,location_att;

!---------------------------------------------------------------------------------------------------------------------
!Initialization for the attribute.
!-----------------------------------------------------------------------------------------------------------!----------
nsf_selfe_output_count     =noutput
nsf_var_output_names(1:nsf_selfe_output_count) = &
["elev",&
"pres" ,&
"airt" ,&
"shum" ,&
"srad" ,&
"flsu" ,&
"fllu" ,&
"radu" ,&
"radd" ,&
"flux" ,&
"evap" ,&
"prcp" ,&
"wind" ,&
"wist" ,&
"dahv" ,&
"w"    ,&
"temp" ,&
"salt" ,&
"conc" ,&
"tdff" ,&
"vdff" ,&
"kine" ,&
"mixl" ,&
"zcor" ,&
"qnon" ,&
"hvel" ];


long_name_att%name      ="long_name";
std_name_att%name       ="standard_name";
base_date_att%name      ="base_date";
unit_att%name           ="unit";
fill_val_int_att%name   ="_FillValue";
fill_val_float_att%name ="_FillValue";
positive_att%name       ="positive";
location_att%name       ="location";

CALL nsf_set_attribute(fill_val_int_att,INT(0,4));
IF(nsf_float == nf90_float )THEN
    CALL nsf_set_attribute(fill_val_float_att,REAL(-9999,4));
ELSEIF(nsf_float == nf90_DOUBLE )THEN
    CALL nsf_set_attribute(fill_val_float_att,REAL(-9999,8));
ENDIF

CALL nsf_set_attribute(positive_att,"down");
CALL nsf_set_attribute(location_att,"node");

!---------------------------------------------------------------------------------------------------------------------
!Add non-time dependented varialbes.
!---------------------------------------------------------------------------------------------------------------------
i=0
!-------------------------------------------
!1.ele;
i = i+1 ;
nsf_self%vars(i)%name       = "ele";
nsf_self%vars(i)%data_type  = nsf_int;
nsf_self%vars(i)%dims_count = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nface_id,nsf_dim_nele_id/);

nsf_self%vars(i)%atts_count  = 3;

CALL nsf_set_attribute(long_name_att,"Horizontal Triangular Element Incidence List");
CALL nsf_set_attribute(unit_att,"index_start_1");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,fill_val_int_att];

                                                        
!-------------------------------------------
!2.0.bndi;
i = i+1 ;
nsf_self%vars(i)%name        = "bndi";
nsf_self%vars(i)%data_type   = nsf_int;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nbi_id/);

nsf_self%vars(i)%atts_count  = 4;

CALL nsf_set_attribute(long_name_att,"Boundary Segment Type List");
CALL nsf_set_attribute(unit_att,"index_start_1");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,fill_val_float_att,fill_val_int_att];

!-------------------------------------------
!2.1.bnd;
i = i+1 ;
nsf_self%vars(i)%name        = "bnd";
nsf_self%vars(i)%data_type   = nsf_int;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nbi_id,nsf_dim_nbnd_id/);

nsf_self%vars(i)%atts_count  = 4;

CALL nsf_set_attribute(long_name_att,"Boundary Segment Node List");
CALL nsf_set_attribute(unit_att,"index_start_1");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,fill_val_float_att,fill_val_int_att];

!-------------------------------------------
!3.lon;
i = i+1 ;
nsf_self%vars(i)%name        = "lon";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id/);

nsf_self%vars(i)%atts_count  = 4;

CALL nsf_set_attribute(long_name_att,"Nodal Longitude");
CALL nsf_set_attribute(unit_att,"degrees_east");
CALL nsf_set_attribute(std_name_att,"longitude");
nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,std_name_att,fill_val_float_att];

!-------------------------------------------
!4.lat;
i = i+1 ;
nsf_self%vars(i)%name        = "lat";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id/);

nsf_self%vars(i)%atts_count  = 4;

CALL nsf_set_attribute(long_name_att,"Nodal Latitude");
CALL nsf_set_attribute(unit_att,"degrees_north");
CALL nsf_set_attribute(std_name_att,"latitude");
nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,std_name_att,fill_val_float_att];

!-------------------------------------------
!5.x;
i = i+1 ;
nsf_self%vars(i)%name        = "x";
nsf_self%vars(i)%data_type   =nsf_float;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id/);

nsf_self%vars(i)%atts_count  = 4;

CALL nsf_set_attribute(long_name_att,"Nodal x-coordinate");
CALL nsf_set_attribute(unit_att,"m");
CALL nsf_set_attribute(std_name_att,"x_coordinate");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,std_name_att,fill_val_float_att];


!-------------------------------------------
!6.y;
i = i+1 ;
nsf_self%vars(i)%name        = "y";
nsf_self%vars(i)%data_type   =nsf_float;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id/);

nsf_self%vars(i)%atts_count  = 4;

CALL nsf_set_attribute(long_name_att,"Nodal y-coordinate");
CALL nsf_set_attribute(unit_att,"m");
CALL nsf_set_attribute(std_name_att,"y_coordinate");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,std_name_att,fill_val_float_att];


!7.depth
i = i+1 ;
nsf_self%vars(i)%name        = "depth";
nsf_self%vars(i)%data_type   =nsf_float;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id/);

nsf_self%vars(i)%atts_count  = 4;

CALL nsf_set_attribute(long_name_att,"Bathymetry");
CALL nsf_set_attribute(unit_att,"m");
CALL nsf_set_attribute(std_name_att,"still_water_depth");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,positive_att,std_name_att,fill_val_float_att];


!8.sigma
i = i+1 ;
nsf_self%vars(i)%name        = "sigma";
nsf_self%vars(i)%data_type   =nsf_float;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nslayer_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Sigma Stretched Vertical Coordinate at Nodes");
CALL nsf_set_attribute(unit_att,"sigma_layer");
CALL nsf_set_attribute(std_name_att,"ocean_sigma_coordinate");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,positive_att,std_name_att,fill_val_float_att];

!nsf_self%var_natt_Array(i) = 4;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Sigma Stretched Vertical Coordinate at Nodes" , &
!                                                         & "ocean_sigma_coordinate" ,"sigma_layer","down"/);

!9.zeta
i = i+1 ;
nsf_self%vars(i)%name        = "zeta";
nsf_self%vars(i)%data_type   =nsf_float;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nzlayer_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Zeta Vertical Coordinate at Nodes");
CALL nsf_set_attribute(unit_att,"m");
CALL nsf_set_attribute(std_name_att,"zeta_coordinate");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,positive_att,std_name_att,fill_val_float_att];


!nsf_self%var_natt_Array(i) = 4;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Zeta Vertical Coordinate at Nodes" , &
!                                                         & "zeta_coordinate" ,"m","down"/);
                                                        
!10.time
i = i+1 ;
nsf_self%vars(i)%name        = "time";
nsf_self%vars(i)%data_type   =nsf_float;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Time");
CALL nsf_set_attribute(unit_att,"seconds since "//TRIM(start_time));
CALL nsf_set_attribute(std_name_att,"time");
CALL nsf_set_attribute(base_date_att,TRIM(start_time));

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,base_date_att,std_name_att,fill_val_float_att];

!-------------------------------------------
!add time depend varialbes.
!-------------------------------------------


!11.elev 
iof_index  =  0;
IF( iof(1) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "elev";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Water Surface Elevation");
CALL nsf_set_attribute(unit_att,"m");
CALL nsf_set_attribute(std_name_att,"sea_surface_elevation");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,positive_att,std_name_att,fill_val_float_att];
ENDIF

!12.pres
IF( iof(2) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "pres";
nsf_self%vars(i)%data_type   =nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Atmospheric pressure at sea surface");
CALL nsf_set_attribute(unit_att,"Pa");
CALL nsf_set_attribute(std_name_att,"air_pressure_at_sea_level");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];

ENDIF

!airt 
IF( iof(3) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "airt";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Air temperature 2m above sea surface");
CALL nsf_set_attribute(unit_att,"Celcius");
CALL nsf_set_attribute(std_name_att,"air_temperature_at_2m_above_sea_surface");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];
ENDIF 
      
!shum
IF( iof(4) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "shum";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Specific humidity 2m above sea surface");
CALL nsf_set_attribute(unit_att,"kg/kg");
CALL nsf_set_attribute(std_name_att,"specific_humidity_at_2m_above_sea_surface");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];
                          
ENDIF 

!srad
IF( iof(5) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "srad";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Solar radiation");
CALL nsf_set_attribute(unit_att,"W m-2");
CALL nsf_set_attribute(std_name_att,"solar_radiation");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];
ENDIF 

!flsu
IF( iof(6) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "flsu";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Turbulent upwelling flux of sensible heat");
CALL nsf_set_attribute(unit_att,"W m-2");
CALL nsf_set_attribute(std_name_att,"turbulent_upwelling_flux_of_sensible_heat");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];

ENDIF 

!fllu  
IF( iof(7) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "fllu";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Turbulent upwelling flux of latent heat");
CALL nsf_set_attribute(unit_att,"W m-2");
CALL nsf_set_attribute(std_name_att,"turbulent_upwelling_flux_of_latent_heat");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];
ENDIF    

IF( iof(8) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "radu";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Upwelling longwave flux at surface");
CALL nsf_set_attribute(unit_att,"W m-2");
CALL nsf_set_attribute(std_name_att,"upwelling_longwave_flux_at_surface");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];
ENDIF                                                     

!radd
IF( iof(9) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "radd";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Downwelling longwave flux at surface");
CALL nsf_set_attribute(unit_att,"W m-2");
CALL nsf_set_attribute(std_name_att,"downwelling_longwave_flux_at_surface");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];

ENDIF   

IF( iof(10) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "flux";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Net heat flux");
CALL nsf_set_attribute(unit_att,"W m-2");
CALL nsf_set_attribute(std_name_att,"net_heat_flux");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];

ENDIF   

!evap
IF( iof(11) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "evap";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Evaporation rate");
CALL nsf_set_attribute(unit_att,"kg m-2 s-1");
CALL nsf_set_attribute(std_name_att,"evaporation_rate");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];
ENDIF   

!prcp
IF( iof(12) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "prcp";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Precipitation rate");
CALL nsf_set_attribute(unit_att,"kg m-2 s-1");
CALL nsf_set_attribute(std_name_att,"precipitation_rate");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];
!                                                        & "m" ,"down","node","-9999"/);
ENDIF   

!wind
IF( iof(13) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "windx";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_two_id,nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Wind speed 10m above sea surface in x direction");
CALL nsf_set_attribute(unit_att,"m/s");
CALL nsf_set_attribute(std_name_att,"wind_speed_10m_above_sea_surface_x ");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];

i = i+1 ;
nsf_self%vars(i)%name        = "windy";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_two_id,nsf_dim_node_id,nsf_dim_time_id/);


nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Wind speed 10m above sea surface in y direction");
CALL nsf_set_attribute(unit_att,"m/s");
CALL nsf_set_attribute(std_name_att,"wind_speed_10m_above_sea_surface_y ");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];

ENDIF 
!wist
IF( iof(14) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "wist";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_two_id,nsf_dim_node_id,nsf_dim_time_id/);
nsf_self%vars(i)%atts_count  = 0;

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Wind stress in x direction");
CALL nsf_set_attribute(unit_att,"m/s");
CALL nsf_set_attribute(std_name_att,"wind_stress_x ");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];
!
i = i+1 ;
nsf_self%vars(i)%name        = "wisty";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);


nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Wind stress in y direction");
CALL nsf_set_attribute(unit_att,"m/s");
CALL nsf_set_attribute(std_name_att,"wind_stress_y ");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];
!                                                         & "m" ,"down","node","-9999"/);

ENDIF 

!dahv
IF( iof(15) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "dahvx";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_two_id,nsf_dim_node_id,nsf_dim_time_id/);
nsf_self%vars(i)%atts_count  = 0;

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Depth averaged horizontal water velocity in x direction");
CALL nsf_set_attribute(unit_att,"m/s");
CALL nsf_set_attribute(std_name_att,"depth_averaged_horizontal_water_velocity_x");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];

i = i+1 ;
nsf_self%vars(i)%name        = "dahvy";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_two_id,nsf_dim_node_id,nsf_dim_time_id/);
nsf_self%vars(i)%atts_count  = 0;

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Depth averaged horizontal water velocity in y direction");
CALL nsf_set_attribute(unit_att,"m/s");
CALL nsf_set_attribute(std_name_att,"depth_averaged_horizontal_water_velocity_y");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,std_name_att,fill_val_float_att];


ENDIF 
 
IF( iof(16) == 1 ) THEN   
i = i+1 ;
!nsf_self%vars(i)%name        = "vert";
nsf_self%vars(i)%name        = "w";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;
CALL nsf_set_attribute(long_name_att,"Upward Water Velocity");
CALL nsf_set_attribute(unit_att,"m/s");
CALL nsf_set_attribute(std_name_att,"upward_vertical_sea_water_velocity");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];


!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
ENDIF   

IF( iof(17) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "temp";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);
nsf_self%vars(i)%atts_count  = 0;

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Temperature");
CALL nsf_set_attribute(unit_att,"Celsius");
CALL nsf_set_attribute(std_name_att,"sea_water_temperature");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];

ENDIF      

IF( iof(18) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "salt";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Salinity");
CALL nsf_set_attribute(unit_att,"ppt");
CALL nsf_set_attribute(std_name_att,"sea_water_salinity");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];
ENDIF   

IF( iof(19) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "conc";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Density");
CALL nsf_set_attribute(unit_att,"ppt");
CALL nsf_set_attribute(std_name_att,"sea_water_salinity");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];

ENDIF    

IF( iof(20) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "tdff";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);


nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Horizontal Eddy Diffusivity");
CALL nsf_set_attribute(unit_att,"m^2/s");
CALL nsf_set_attribute(std_name_att,"horizontal_eddy_diffusivity");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];

ENDIF                                          
          
IF( iof(21) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "vdff";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count        = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Vertical Eddy Diffusivity");
CALL nsf_set_attribute(unit_att,"m^2/s");
CALL nsf_set_attribute(std_name_att,"vertical_eddy_diffusivity");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];


ENDIF     

IF( iof(22) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "kine";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Turbulent Kinetic Energy");
CALL nsf_set_attribute(unit_att,"m^2/s^2");
CALL nsf_set_attribute(std_name_att,"turbulent_kinetic_energy");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];



ENDIF       

IF( iof(23) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "mixl";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);
nsf_self%vars(i)%atts_count  = 0;

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Macroscale Mixing Length");
CALL nsf_set_attribute(unit_att,"m");
CALL nsf_set_attribute(std_name_att,"macroscale_mixing_length");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];

ENDIF   
 
IF( iof(24) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "zcor";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count        = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Height of Vertical Layers");
CALL nsf_set_attribute(unit_att,"m");
CALL nsf_set_attribute(std_name_att,"layer_height");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];
!                                                         & "m" ,"down","node","-9999"/);
ENDIF      

IF( iof(25) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "qnon";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count        = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Nonhydrostatic Pressure");
CALL nsf_set_attribute(unit_att,"N/m^2");
CALL nsf_set_attribute(std_name_att,"nonhydrostatic_pressure");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];


ENDIF    

IF( iof(26) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "u";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"Eastward Water Velocity");
CALL nsf_set_attribute(unit_att,"m/s");
CALL nsf_set_attribute(std_name_att,"eastward_sea_water_velocity");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];


i = i+1 ;
nsf_self%vars(i)%name        = "v";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nlayer_id,nsf_dim_node_id,nsf_dim_time_id/);
nsf_self%vars(i)%atts_count  = 0;

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"northward_sea_water_velocity");
CALL nsf_set_attribute(unit_att,"m/s");
CALL nsf_set_attribute(std_name_att,"northward_sea_water_velocity");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];

ENDIF     

!kbp00
i = i+1 ;
nsf_self%vars(i)%name        = "kbp00";  ! start (bottom) index of the vertical layer
nsf_self%vars(i)%data_type   = nsf_int;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id/);
nsf_self%vars(i)%atts_count  = 0;

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"bottom vertical index");
CALL nsf_set_attribute(unit_att,"NON");
CALL nsf_set_attribute(std_name_att,"bottom vertical index");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_int_att,std_name_att];


!------------------------------------------------------------
!write theta_b, theta_f and h_s
!------------------------------------------------------------
!theta_b
i = i+1 ;
nsf_self%vars(i)%name        = "theta_b";  ! theta_f used to compute Slayer
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_two_id/);
nsf_self%vars(i)%atts_count  = 0;

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"theta_b");
CALL nsf_set_attribute(unit_att,"NON");
CALL nsf_set_attribute(std_name_att,"theta_b");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];

!theta_f
i = i+1 ;
nsf_self%vars(i)%name        = "theta_f";  ! theta_f used to compute Slayer
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_two_id/);
nsf_self%vars(i)%atts_count  = 0;

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"theta_f");
CALL nsf_set_attribute(unit_att,"NON");
CALL nsf_set_attribute(std_name_att,"theta_f");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];

!h_s
i = i+1 ;
nsf_self%vars(i)%name        = "h_s";  ! h_s used to compute Slayer
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 1;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_two_id/);
nsf_self%vars(i)%atts_count  = 0;

nsf_self%vars(i)%atts_count  = 5;

CALL nsf_set_attribute(long_name_att,"h_s");
CALL nsf_set_attribute(unit_att,"NON");
CALL nsf_set_attribute(std_name_att,"h_s");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,location_att,fill_val_float_att,std_name_att];



nsf_self%vars_count = i;
                   
END SUBROUTINE

!estimate number of varialbes.
SUBROUTINE  nfs_self_est_nvars(iof,nvars)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind),DIMENSION(:),INTENT(IN   ) :: iof;
INTEGER(intKind)             ,INTENT(  OUT) :: nvars;
!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)           :: i;

nvars   = 0;
i=0
!---------------------------------------------------------------------------------------------------------------------
!Add non-time dependented varialbes.
! 11 non-time dependent varialbes.
!---------------------------------------------------------------------------------------------------------------------
!nvars = nvars + 11 + 3 theate_f, theta_b, h_s
nvars = nvars + 15;
!-------------------------------------------
!add time depend varialbes.
!-------------------------------------------
nvars = nvars + sum(iof(1:nsf_selfe_output_count));
!---26 have u and v;
IF(iof(26) == 1 ) THEN
nvars = nvars + 1;
ENDIF
                                                       
END SUBROUTINE



!------------------------------------------------------------------------------------------------------------------
!output the grid information ordered by the global index......
!This procedur is complex and can be replace by only output by rank 0 and using the data read form grd3
! and vgrid.
!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE nsf_selfe_out_grid()
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                       :: temp_var_id,index_base,status,index_base_ne,iproc;
INTEGER(intKind)                       :: i,istat;
INTEGER(intKind),DIMENSION(nproc)      :: np_array,ne_array;
INTEGER(intKind)                       :: global_index;
!---------------------------------------------------------------------------------------------------------------------
!out put node information.
!---------------------------------------------------------------------------------------------------------------------


!CALL MPI_ALLGATHER(np,1,itype,np_array,1,itype,comm,istat);
!CALL MPI_ALLGATHER(ne,1,itype,ne_array,1,itype,comm,istat);
!IF(myrank ==0 )THEN
!      index_base    = 1;
!      index_base_ne = 1;
!ELSE
!      index_base    = SUM(np_array(1:myrank)) + 1 ;
!      index_base_ne = SUM(ne_array(1:myrank)) + 1 ;
!ENDIF
!-------------STAR TO OUT PUT TO NETCDF FILE
! form proc 0 to nproc-1;
!-------------------open the file-------------------   

DO iproc = 1,nproc

    IF(myrank+1 == iproc) THEN !CURRENT proc to output
 
!    WRITE(*,*) "rank ",myrank, " start to open file";  
!    WRITE(*,*) "rank ",myrank, "file Name ", TRIM(nc_file_name); 
!------open the file for myrank------
     CALL nfs_open_file(nsf_self,TRIM(nc_file_name), NF90_WRITE);
     nsf_self_ncid = nsf_self%ncid;    
    
!    WRITE(*,*) "myrank :" ,nsf_self_ncid;
    
!    out put the element.
       CALL nsf_get_var_id(nsf_self,"ele",temp_var_id) 
       DO i = 1,ne     
        global_index = ielg(i)
        istat = nf90_put_var(nsf_self_ncid,temp_var_id,nmgb(global_index,1:3), &
                           (/1,global_index/),(/3,1/));                    
        CALL  check_error(istat)!    out put the node x 
       ENDDO     

!    out put the x 
      CALL nsf_get_var_id(nsf_self,"x",temp_var_id)
      DO i = 1,np
          global_index =iplg(i);
!          istat = nf90_put_var(nsf_self_ncid,temp_var_id,REAL(xnd(i:i),4), & !version 3.1h
!                            (/global_index/),(/1/));  
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,REAL(x(i:i),4), &    !version 3.0b
                            (/global_index/),(/1/));                  
          CALL  check_error(istat)      
      ENDDO
      
!    out put the y 
      CALL nsf_get_var_id(nsf_self,"y",temp_var_id)
      DO i = 1,np
          global_index =iplg(i);
!          istat = nf90_put_var(nsf_self_ncid,temp_var_id,REAL(ynd(i:i),4), & !version 3.1h
!                           (/global_index/),(/1/));  
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,REAL(y(i:i),4), &    !version 3.0b
                           (/global_index/),(/1/));                             
          CALL  check_error(istat)      
      ENDDO

!    out put the depth 
      CALL nsf_get_var_id(nsf_self,"depth",temp_var_id)
      DO i = 1,np
          global_index =iplg(i);
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,REAL(dp00(i:i),4), &
                           (/global_index/),(/1/));                    
          CALL  check_error(istat)      
      ENDDO        


!    out put the node xlon 
      CALL nsf_get_var_id(nsf_self,"lon",temp_var_id)
      DO i = 1,np
          global_index =iplg(i);
!change the xlon from radius to degress          
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,REAL(xlon(i:i)*180/pi,4), &
                           (/global_index/),(/1/));                    
          CALL  check_error(istat) 
      ENDDO

!    out put the node ylat 
      CALL nsf_get_var_id(nsf_self,"lat",temp_var_id)
      DO i = 1,np
          global_index =iplg(i);
!change the ylon from radius to degress                
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,REAL(ylat(i:i)*180/pi,4), &
                           (/global_index/),(/1/));                    
          CALL  check_error(istat) 
      ENDDO
      
!--------------------------------------------------------------------
!bottom layer index..kbp00
!--------------------------------------------------------------------
      CALL nsf_get_var_id(nsf_self,"kbp00",temp_var_id)
      DO i = 1,np
          global_index =iplg(i);
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,REAL(kbp00(i:i),4), &
                           (/global_index/),(/1/));                    
          CALL  check_error(istat) 
      ENDDO

!-------------------------------------------------------------------------
!sigma layer and zlayer only need to be outputed by rank 0
!-------------------------------------------------------------------------
      IF(myrank ==0) THEN
      
!-------------------------------------------------------------------------
!  out put the open boundary ......
!-------------------------------------------------------------------------
          DO i=1, nope_global
!----------------bounday type -----          
              CALL nsf_get_var_id(nsf_self,"bndi",temp_var_id);
              istat = nf90_put_var(nsf_self_ncid,temp_var_id,(/0/), &
                               (/i/),(/1/));                    
              CALL  check_error(istat)  
!----------------bounday node index.          
              CALL nsf_get_var_id(nsf_self,"bnd",temp_var_id);
              istat = nf90_put_var(nsf_self_ncid,temp_var_id,iond_global(i,:), &
                               (/i,1/),(/1,nond_global(i)/));                    
              CALL  check_error(istat)                
                    
          ENDDO
          
!-------------------------------------------------------------------------
!  out put the open boundary ......
!-------------------------------------------------------------------------
          DO i=1, nland_global
!----------------bounday type -----          
              CALL nsf_get_var_id(nsf_self,"bndi",temp_var_id);
              istat = nf90_put_var(nsf_self_ncid,temp_var_id,(/1/), &
                               (/i+nope_global/),(/1/));                    
              CALL  check_error(istat)  
!----------------bounday node index.          
              CALL nsf_get_var_id(nsf_self,"bnd",temp_var_id);
              istat = nf90_put_var(nsf_self_ncid,temp_var_id,ilnd_global(i,:), &
                               (/i+nope_global,1/),(/1,nlnd_global(i)/));                    
              CALL  check_error(istat)                
                    
          ENDDO          
               
      
!   out put sigma layer.      
          CALL nsf_get_var_id(nsf_self,"sigma",temp_var_id);
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,REAL(sigma(kz:nvrt),4), &
                           (/1/),(/nsig/));                    
          CALL  check_error(istat)  
!   output zlayer........
          CALL nsf_get_var_id(nsf_self,"zeta",temp_var_id)
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,REAL(ztot(1:kz),4), &
                           (/1/),(/kz/));                    
          CALL  check_error(istat)                

!--------------------------------------------------------------------
!bottom layer index..theate_f
!--------------------------------------------------------------------
         CALL nsf_get_var_id(nsf_self,"theta_f",temp_var_id)
              istat = nf90_put_var(nsf_self_ncid,temp_var_id,(/REAL(theta_f,4)/), &
                               (/1/),(/1/)); 
         CALL  check_error(istat) 
                                              
      
         CALL nsf_get_var_id(nsf_self,"theta_b",temp_var_id)
              istat = nf90_put_var(nsf_self_ncid,temp_var_id,(/REAL(theta_b,4)/), &
                                (/1/),(/1/));          
         CALL  check_error(istat) 
     
         CALL nsf_get_var_id(nsf_self,"h_s",temp_var_id)
              istat = nf90_put_var(nsf_self_ncid,temp_var_id,(/REAL(h_s,4)/), &
                               (/1/),(/1/));      
         CALL  check_error(istat)
 
 

                    
       ENDIF  !IF(myrank ==0) THEN
        
!--------------------------------------------------------------------
!close the file for myrank then it can be opened by another rank.  
     
       istat = nf90_sync(nsf_self_ncid);
       CALL  check_error(istat) 
       CALL  nfs_close_file(nsf_self);     
     ENDIF !IF(myrank+1 == iproc) THEN !CURRENT proc to output

!barrier all the process here...
CALL  mpi_barrier(comm,ierr)   
!WRITE(*,*) "proc rank =: ",myrank,"output grid finished."

ENDDO !DO iproc = 1,nproc



CALL  mpi_barrier(comm,ierr)
IF(myrank ==0)THEN
    WRITE(*,*) "output grid file finished."
ENDIF
END SUBROUTINE


!---------------------------------------------------------------------------------------------------------------------
!read the geometry information from the netcdf file
!NOTE: only for single proc to run 
!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE nsf_selfe_in_grid()
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                       :: temp_var_id,index_base,status,index_base_ne,iproc;
INTEGER(intKind)                       :: i,istat;
INTEGER(intKind),DIMENSION(nproc)      :: np_array,ne_array;
INTEGER(intKind)                       :: global_index;

REAL(rkind)                            :: rtemp(1);
!INTEGER(intKind)                       :: temp_nm(3,ne);
!---------------------------------------------------------------------------------------------------------------------
!get node information.
!---------------------------------------------------------------------------------------------------------------------
nsf_self_ncid = nsf_self%ncid;

CALL nsf_get_var_id(nsf_self,"ele",temp_var_id)

!istat = nf90_get_var(nsf_self_ncid,temp_var_id,nm,(/1,1/),(/3,ne/)); 
!CALL  check_error(istat)  
!
!nm = TRANSPOSE(nm);
DO i = 1,ne
    istat = nf90_get_var(nsf_self_ncid,temp_var_id,nm(i,:),(/1,i/),(/3,1/)); 
ENDDO
CALL  check_error(istat)  
!nm = TRANSPOSE(temp_nm);

!    get the x 
CALL nsf_get_var_id(nsf_self,"x",temp_var_id)
!CALL nsf_get_var_id(nsf_self,"lon",temp_var_id)
!istat =  nf90_get_var(nsf_self_ncid,temp_var_id,xnd,(/1/),(/np/));  !version 3.1h
istat =  nf90_get_var(nsf_self_ncid,temp_var_id,x,(/1/),(/np/));   !version 3.0
CALL check_error(istat)  
    
!    get the y 
CALL nsf_get_var_id(nsf_self,"y",temp_var_id)
!CALL nsf_get_var_id(nsf_self,"lat",temp_var_id)
!istat =  nf90_get_var(nsf_self_ncid,temp_var_id,ynd,(/1/),(/np/)); !version 3.1h
istat =  nf90_get_var(nsf_self_ncid,temp_var_id,y,(/1/),(/np/));    !version 3.0
CALL check_error(istat)  

!    get the depth 
CALL nsf_get_var_id(nsf_self,"depth",temp_var_id)
istat =  nf90_get_var(nsf_self_ncid,temp_var_id,dp,(/1/),(/np/)); 
CALL check_error(istat)       

!   get sigma layer.      
CALL nsf_get_var_id(nsf_self,"sigma",temp_var_id);
istat = nf90_get_var(nsf_self_ncid,temp_var_id,sigma(kz:nvrt), &
               (/1/),(/nsig/));   
 
!   get zlayer layer.
CALL nsf_get_var_id(nsf_self,"zeta",temp_var_id)
istat = nf90_get_var(nsf_self_ncid,temp_var_id,ztot, &
               (/1/),(/kz/));                    
CALL  check_error(istat)   
                                          

!--------------------------------------------------------------------
!get bottom layer index..kbp00
!--------------------------------------------------------------------
CALL nsf_get_var_id(nsf_self,"kbp00",temp_var_id)
istat = nf90_get_var(nsf_self_ncid,temp_var_id,kbp, &
               (/1/),(/np/)); 
CALL  check_error(istat)  
         
         
         
!--------------------------------------------------------------------
!get bottom layer index theta_f;
!--------------------------------------------------------------------
CALL nsf_get_var_id(nsf_self,"theta_f",temp_var_id)
istat = nf90_get_var(nsf_self_ncid,temp_var_id,rtemp, &
               (/1/),(/1/)); 
CALL  check_error(istat)  
theta_f = rtemp(1);

!--------------------------------------------------------------------
!get bottom layer index theta_b;
!--------------------------------------------------------------------
CALL nsf_get_var_id(nsf_self,"theta_b",temp_var_id)
istat = nf90_get_var(nsf_self_ncid,temp_var_id,rtemp, &
               (/1/),(/1/)); 
CALL  check_error(istat)  
theta_b = rtemp(1);

!--------------------------------------------------------------------
!get bottom layer index theta_b;
!--------------------------------------------------------------------
CALL nsf_get_var_id(nsf_self,"h_s",temp_var_id)
istat = nf90_get_var(nsf_self_ncid,temp_var_id,rtemp, &
               (/1/),(/1/)); 
CALL  check_error(istat)  
h_s = rtemp(1);
                      
END SUBROUTINE



!--------------------------------------------------------
! read the dimesnion of the array in netcdf and 
!   
!--------------------------------------------------------
SUBROUTINE nsf_selfe_acquire_dims()
IMPLICIT NONE;
!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
!CHARACTER(LEN=*)              :: fileName;

!--------------------------------------------------------
!Arguments.
!--------------------------------------------------------
!CHARACTER(LEN=*)           :: fileName;
!--------------------------------------------------------
!Local variables. m,mm;
!--------------------------------------------------------
INTEGER                    :: temp_dim_id,temp_len,istat;

!---------------------------------------------------------------------------------------------------------------------
!out put node information.
!---------------------------------------------------------------------------------------------------------------------
nsf_self_ncid = nsf_self%ncid;

istat = nf90_inq_dimid(nsf_self_ncid,"node", temp_dim_id)
istat =nf90_Inquire_Dimension(nsf_self_ncid, temp_dim_id, len=temp_len)
np        = temp_len;
np_global = temp_len;


istat = nf90_inq_dimid(nsf_self_ncid,"nele", temp_dim_id)
istat =nf90_Inquire_Dimension(nsf_self_ncid, temp_dim_id, len=temp_len)
ne        = temp_len;
ne_global = temp_len;


istat = nf90_inq_dimid(nsf_self_ncid,"nzlayer", temp_dim_id)
istat =nf90_Inquire_Dimension(nsf_self_ncid, temp_dim_id, len=temp_len)
kz        = temp_len;

istat = nf90_inq_dimid(nsf_self_ncid,"nslayer", temp_dim_id)
istat =nf90_Inquire_Dimension(nsf_self_ncid, temp_dim_id, len=temp_len)
nsig        = temp_len;

istat = nf90_inq_dimid(nsf_self_ncid,"nlayer", temp_dim_id)
istat =nf90_Inquire_Dimension(nsf_self_ncid, temp_dim_id, len=temp_len)
nvrt      = temp_len;


istat = nf90_inq_dimid(nsf_self_ncid,"time", temp_dim_id)
istat =nf90_Inquire_Dimension(nsf_self_ncid, temp_dim_id, len=temp_len)
nrec      = temp_len;

END SUBROUTINE



!--------------------------------------------------------
! nsf_self_netcdf_time_data_in
!--------------------------------------------------------
SUBROUTINE nsf_self_netcdf_time_data_in_1D(varName,outb,itime,ndata)
IMPLICIT NONE;

CHARACTER(LEN=*)                          :: varName;
REAL(rKind) ,  DIMENSION(:),INTENT(  OUT) :: outb;
INTEGER                    ,INTENT(IN   ) :: itime;
INTEGER     ,  DIMENSION(:),INTENT(IN   ) :: ndata;

!--------------------------------------------------------
!Local varialbes.
!--------------------------------------------------------
INTEGER          :: temp_var_id,istat;
REAL(4)          :: temp_data(ndata(1));

temp_var_id = 0;

CALL nsf_get_var_id(nsf_self,TRIM(varName),temp_var_id)

istat =  nf90_get_var(nsf_self_ncid,temp_var_id,outb,&
                      (/1,itime/),(/ndata(1),1/)); 
CALL check_error(istat);

!IF( nsf_float == NF90_float .AND.  rKind     == 8 )THEN 
!
!    istat =  nf90_get_var(nsf_self_ncid,temp_var_id,temp_data,&
!                          (/1,itime/),(/ndata(1),1/)); 
!    outb  =temp_data;
!ELSEIF(nsf_float == NF90_double .AND. rKind ==8 ) THEN
!    istat =  nf90_get_var(nsf_self_ncid,temp_var_id,outb,&
!                         (/1,itime/),(/ndata(1),1/)); 
!ELSE
!    WRITE(*,*) "nsf_float and rkind are not condiserned right now"
!ENDIF
!

                          
CALL check_error(istat)  


END SUBROUTINE 



!--------------------------------------------------------
! get simulation time at itime step
!--------------------------------------------------------
SUBROUTINE nsf_self_netcdf_time_time(varName,out_time,itime)
IMPLICIT NONE;

CHARACTER(LEN=*)                          :: varName;
REAL(rKind)                ,INTENT(  OUT) :: out_time;
INTEGER                    ,INTENT(IN   ) :: itime;
!INTEGER     ,  DIMENSION(:),INTENT(IN   ) :: ndata;

!--------------------------------------------------------
!Local varialbes.
!--------------------------------------------------------
INTEGER          :: temp_var_id,istat;
REAL(4)          :: temp_data(1);

temp_var_id = 0;

CALL nsf_get_var_id(nsf_self,TRIM(varName),temp_var_id)

istat =  nf90_get_var(nsf_self_ncid,temp_var_id,temp_data,&
                      (/itime/),(/1/)); 
CALL check_error(istat);

out_time = temp_data(1);

!IF( nsf_float == NF90_float .AND.  rKind     == 8 )THEN 
!
!    istat =  nf90_get_var(nsf_self_ncid,temp_var_id,temp_data,&
!                          (/1,itime/),(/ndata(1),1/)); 
!    outb  =temp_data;
!ELSEIF(nsf_float == NF90_double .AND. rKind ==8 ) THEN
!    istat =  nf90_get_var(nsf_self_ncid,temp_var_id,outb,&
!                         (/1,itime/),(/ndata(1),1/)); 
!ELSE
!    WRITE(*,*) "nsf_float and rkind are not condiserned right now"
!ENDIF
!

                          
CALL check_error(istat)  


END SUBROUTINE 

!--------------------------------------------------------
! nsf_self_netcdf_time_data_in
!--------------------------------------------------------
SUBROUTINE nsf_self_netcdf_time_data_in_2D(varName,outb,itime,ndata)
IMPLICIT NONE;

CHARACTER(LEN=*)                            :: varName;
REAL(rKind) ,  DIMENSION(:,:),INTENT(  OUT) :: outb;
INTEGER                      ,INTENT(IN   ) :: itime;
INTEGER     ,  DIMENSION(1:2),INTENT(IN   ) :: ndata;

!--------------------------------------------------------
!Local varialbes.
!--------------------------------------------------------
INTEGER                    :: temp_var_id,istat,i;
REAL(4)                    :: temp_data(ndata(2));

temp_var_id = 0;
CALL nsf_get_var_id(nsf_self,TRIM(varName),temp_var_id)


DO i =1,ndata(1) ! np_global

    istat =  nf90_get_var(nsf_self_ncid,temp_var_id,outb(i,:), &
                           start=(/1,i,itime/),count=(/ndata(2),1,1/));   
    CALL check_error(istat)
   
ENDDO

!IF( nsf_float == NF90_float .AND.  rKind     == 8 )THEN 
!    
!
!!    istat =  nf90_get_var(nsf_self_ncid,temp_var_id,temp_data,&
!!                           (itime-1)*ndata(1:2)+1,ndata(1:2)); 
!    DO i =1,ndata(1) ! np_global
!        istat =  nf90_get_var(nsf_self_ncid,temp_var_id,temp_data, &
!                               start=(/1,i,itime/),count=(/ndata(2),1,1/));   
!        CALL check_error(istat)
!        outb(i,:)= temp_data;     
!    ENDDO
!!        istat =  nf90_get_var(nsf_self_ncid,temp_var_id,temp_data, &
!!                               start=(/1,1,itime/),count=(/ndata(2),ndata(1),1/));   
!!        
!
!ELSEIF(nsf_float == NF90_double .AND. rKind ==8 ) THEN
!!    istat =  nf90_get_var(nsf_self_ncid,temp_var_id,outb,&
!!                          (itime-1)*ndata(1:2)+1,ndata(1:2)); 
!ELSE
!    WRITE(*,*) "nsf_float and rkind are not condiserned right now"
!ENDIF


END SUBROUTINE 


!--------------------------------------------------------
! nsf_self_netcdf_time_data_in
!--------------------------------------------------------
SUBROUTINE nsf_self_netcdf_time_data_in_3D(varName,outb,itime,ndata)
IMPLICIT NONE;

CHARACTER(LEN=*)                              :: varName;
REAL(rKind) ,  DIMENSION(:,:,:),INTENT(  OUT) :: outb;
INTEGER                        ,INTENT(IN   ) :: itime;
INTEGER     ,  DIMENSION(:)    ,INTENT(IN   ) :: ndata;

!--------------------------------------------------------
!Local varialbes.
!--------------------------------------------------------
INTEGER                    :: temp_var_id,i,j,istat;
REAL(4)                    :: temp_data(ndata(3));

temp_var_id = 0;

CALL nsf_get_var_id(nsf_self,TRIM(varName),temp_var_id)

DO j = 1,ndata(1)    ! np_global
    DO i =1,ndata(2) ! nvrt
        istat =  nf90_get_var(nsf_self_ncid,temp_var_id,outb(j,i,:), &
                               start=(/1,i,j,itime/),count=(/ndata(3),1,1,1/));   
        CALL check_error(istat)
    !        outb(i,:)= temp_data;   
    ENDDO  
ENDDO

!IF( nsf_float == NF90_float .AND.  rKind     == 4 )THEN 
!    istat =  nf90_get_var(nsf_self_ncid,temp_var_id,temp_data,&
!                           (itime-1)*ndata(1:3)+1,ndata(1:3)); 
!    outb  =temp_data;
!ELSEIF(nsf_float == NF90_double .AND. rKind ==8 ) THEN
!    istat =  nf90_get_var(nsf_self_ncid,temp_var_id,outb,&
!                           (itime-1)*ndata(1:3)+1,ndata(1:3)); 
!ELSE
!    WRITE(*,*) "nsf_float and rkind are not condiserned right now"
!ENDIF


END SUBROUTINE 




END MODULE 

