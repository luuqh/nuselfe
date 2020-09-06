MODULE nsf_self_vars_mod
USE nsf_container_mod
USE elfe_glbl
USE elfe_msgp
USE funMaxModule
IMPLICIT NONE;

CHARACTER(LEN=200)                   :: nc_file_name;  


INTEGER(intKind)                     :: nsf_selfe_dims_count       =10;
INTEGER(intKind)                     :: nsf_selfe_global_att_count =10;
INTEGER(intKind)                     :: nsf_selfe_output_count     =26; !only time dependent variables.
INTEGER(intKind)                     :: nsf_selfe_total_vars          ;

INTEGER(intKind)                     :: global_id;

TYPE(Tnsf_container)                 :: nsf_self;
REAL(realKind)                       :: temp_out(1),temp_out2(1);
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
!REAL(MK)                             :: temp_var_val,temp_var_val2;

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
INTEGER(intKind):: nsf_dim_two_id    ;

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
CALL nfs_self_est_nvars(iof,nsf_selfe_total_vars);

!---------------------------------------------------------------------------------------------------------------------
!2.init the container first when n_att,n_dim and n_var are known.
!---------------------------------------------------------------------------------------------------------------------
CALL init_nsf_container(container,nsf_selfe_global_att_count, &
                        nsf_selfe_dims_count,nsf_selfe_total_vars);

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
DO i =1 ,nsf_selfe_output_count
    IF( i == 26) THEN
        CALL nsf_get_var_id(container,"u",nsf_var_u_id);
        CALL nsf_get_var_id(container,"v",nsf_var_v_id);
    ELSE
        CALL nsf_get_var_id(container,TRIM(nsf_var_output_names(i)),nsf_var_output_ids(i));
    ENDIF
ENDDO

END SUBROUTINE

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
!
!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE nsf_selfe_nectcdf_out_init(noutput)
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
INTEGER(intKind)                ::noutput;


nsf_selfe_dims_count       =10
nsf_selfe_global_att_count =10
nsf_selfe_output_count     =noutput

nsf_out_ntime          =0;


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
      
nsf_selfe_dim_vals(1:nsf_selfe_dims_count) =[np_global,ne_global,max(neta_global,1) ,& 
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
INTEGER(intKind)                :: i;
INTEGER(intKind)                :: iof_index;
TYPE(Tnsf_attribute)            :: long_name_att,std_name_att,unit_att,fill_val_int_att,fill_val_float_att;
TYPE(Tnsf_attribute)            :: positive_att,base_date_att,location_att;

!---------------------------------------------------------------------------------------------------------------------
!Initialization for the attribute.
!---------------------------------------------------------------------------------------------------------------------
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
!2.bnd;
i = i+1 ;
nsf_self%vars(i)%name        = "bnd";
nsf_self%vars(i)%data_type   = nsf_int;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_nbi_id,nsf_dim_nbi_id/);

nsf_self%vars(i)%atts_count  = 4;

CALL nsf_set_attribute(long_name_att,"Boundary Segment Type List");
CALL nsf_set_attribute(unit_att,"index_start_1");

nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
                          [long_name_att,unit_att,fill_val_float_att,fill_val_int_att];

!nsf_self%var_natt_Array(i) = 4;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  = (/"long_name","standard_name","units","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Boundary Segment Node List" , &
!                                                       & "mesh element" ,"index_start_1","-9999"/);

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

!nsf_self%var_natt_Array(i)   = 3;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Nodal Longitude" , "longitude" ,"degrees_east" /);


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

!nsf_self%var_natt_Array(i) = 3;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Nodal Latitude" , "latitude" ,"degrees_north" /);


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


!nsf_self%var_natt_Array(i) = 3;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Nodal x-coordinate" , "x_coordinate" ,"m" /);

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


!nsf_self%var_natt_Array(i) = 3;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Nodal y-coordinate" , "y_coordinate" ,"m" /);

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
CALL nsf_set_attribute(unit_att,"days since 2008-08-22 12:00:00 00:00");
CALL nsf_set_attribute(std_name_att,"time");
CALL nsf_set_attribute(base_date_att,"2008, 8, 22, 12");

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

CALL nsf_set_attribute(long_name_att,"Time");
CALL nsf_set_attribute(unit_att,"days since 2008-08-22 12:00:00 00:00");
CALL nsf_set_attribute(std_name_att,"time");

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

nsf_self%vars(i)%atts_count  = 0;

!CALL nsf_set_attribute(long_name_att,"Time");
!CALL nsf_set_attribute(unit_att,"days since 2008-08-22 12:00:00 00:00");
!CALL nsf_set_attribute(std_name_att,"time");
!
!nsf_self%vars(i)%atts(1:nsf_self%vars(i)%atts_count) = &
!                          [long_name_att,unit_att,location_att,positive_att,std_name_att];
!                                                         & "m" ,"down","node","-9999"/);
ENDIF
 
IF( iof(3) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "airt";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
ENDIF 
      

IF( iof(4) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "shum";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
ENDIF 

IF( iof(5) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "srad";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
ENDIF 

IF( iof(6) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "flsu";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);
nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
ENDIF 
  
IF( iof(7) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "fllu";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
ENDIF    

IF( iof(8) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "radu";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
ENDIF                                                     

IF( iof(9) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "radd";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 0;

!!nsf_self%var_natt_Array(i) = 6;
!!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!!                                                         & "m" ,"down","node","-9999"/);
ENDIF   

IF( iof(10) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "flux";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
ENDIF   


IF( iof(11) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "evap";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
ENDIF   


IF( iof(12) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "prcp";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);

nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
ENDIF   

IF( iof(13) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "wind";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 3;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_two_id,nsf_dim_node_id,nsf_dim_time_id/);
nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);

!i = i+1 ;
!nsf_self%vars(i)%name        = "windy";
!nsf_self%vars(i)%data_type   = nsf_float;
!nsf_self%vars(i)%dims_count        = 2;
!nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);
!nsf_self%vars(i)%atts_count  = 0;
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);

ENDIF 

IF( iof(14) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "wist";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_two_id,nsf_dim_node_id,nsf_dim_time_id/);
nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!
!i = i+1 ;
!nsf_self%vars(i)%name        = "wisty";
!nsf_self%vars(i)%data_type   = nsf_float;
!nsf_self%vars(i)%dims_count        = 2;
!nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);
!nsf_self%vars(i)%atts_count  = 0;
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);

ENDIF 

IF( iof(15) == 1 ) THEN   
i = i+1 ;
nsf_self%vars(i)%name        = "dahv";
nsf_self%vars(i)%data_type   = nsf_float;
nsf_self%vars(i)%dims_count  = 2;
nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_two_id,nsf_dim_node_id,nsf_dim_time_id/);
nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);

!i = i+1 ;
!nsf_self%vars(i)%name        = "dahvy";
!nsf_self%vars(i)%data_type   = nsf_float;
!nsf_self%vars(i)%dims_count        = 2;
!nsf_self%vars(i)%dims_ids(1:nsf_self%vars(i)%dims_count)     = (/nsf_dim_node_id,nsf_dim_time_id/);
!nsf_self%vars(i)%atts_count  = 0;
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);

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
nsf_self%vars(i)%dims_count        = 3;
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

nsf_self%vars(i)%atts_count  = 0;

!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
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
! 10 non-time dependent varialbes.
!---------------------------------------------------------------------------------------------------------------------
nvars = nvars + 10

!-------------------------------------------
!add time depend varialbes.
!-------------------------------------------
nvars = nvars + sum(iof(1:nsf_selfe_output_count));
!---26 have u and v;
IF(iof(26) == 1 ) THEN
nvars = nvars + 1;
ENDIF
                                                       
END SUBROUTINE



!---------------------------------------------------------------------------------------------------------------------
!output the grid information ordered by the global index......
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
nsf_self_ncid = nsf_self%ncid;

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
DO iproc = 1,nproc

    IF(myrank+1 == iproc) THEN !CURRENT proc to output
    
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
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,x(i:i), &
                            (/global_index/),(/1/));                    
          CALL  check_error(istat)      
      ENDDO
      
!    out put the y 
      CALL nsf_get_var_id(nsf_self,"y",temp_var_id)
      DO i = 1,np
          global_index =iplg(i);
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,y(i:i), &
                           (/global_index/),(/1/));                    
          CALL  check_error(istat)      
      ENDDO

!    out put the depth 
      CALL nsf_get_var_id(nsf_self,"depth",temp_var_id)
      DO i = 1,np
          global_index =iplg(i);
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,dp00(i:i), &
                           (/global_index/),(/1/));                    
          CALL  check_error(istat)      
      ENDDO        


!    out put the node xlon 
      CALL nsf_get_var_id(nsf_self,"lon",temp_var_id)
      DO i = 1,np
          global_index =iplg(i);
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,xlon(i:i), &
                           (/global_index/),(/1/));                    
          CALL  check_error(istat) 
      ENDDO

!    out put the node ylat 
      CALL nsf_get_var_id(nsf_self,"lat",temp_var_id)
      DO i = 1,np
          global_index =iplg(i);
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,ylat(i:i), &
                           (/global_index/),(/1/));                    
          CALL  check_error(istat) 
      ENDDO
      

!-------------------------------------------------------------------------
!sigma layer and zlayer only need to be outputed by rank 0
!-------------------------------------------------------------------------
      IF(myrank ==0) THEN
!   out put sigma layer.      
          CALL nsf_get_var_id(nsf_self,"sigma",temp_var_id);
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,sigma(kz:nvrt), &
                           (/1/),(/nsig/));                    
          CALL  check_error(istat)  
!   output zlayer........
          CALL nsf_get_var_id(nsf_self,"zeta",temp_var_id)
          istat = nf90_put_var(nsf_self_ncid,temp_var_id,ztot(1:kz), &
                           (/1/),(/kz/));                    
          CALL  check_error(istat)                
      ENDIF

            
     ENDIF !IF(myrank+1 == iproc) THEN !CURRENT proc to output

!sync before barrier.......................    
    istat = nf90_sync(nsf_self_ncid);
    CALL  check_error(istat)      
    !--------------------------------------------------------------------
    !For parallel code marrier the procs.       
#ifdef USE_MPIMODULE     
    ! barrier the proc
     
          CALL  mpi_barrier(comm,ierr)   
#endif
ENDDO


END SUBROUTINE



!SUBROUTINE  nfs_self_set_atts(nsf_self)
!!---------------------------------------------------------------------------------------------------------------------
!!  Modules 
!!---------------------------------------------------------------------------------------------------------------------
!IMPLICIT NONE
!!---------------------------------------------------------------------------------------------------------------------
!!  Arguments     
!!---------------------------------------------------------------------------------------------------------------------
!TYPE(Tnsf_container),INTENT(INOUT) :: nsf_self;
!!---------------------------------------------------------------------------------------------------------------------
!!  Local variables 
!!---------------------------------------------------------------------------------------------------------------------
!CHARACTER(CharLen)       :: version,ctime;
!INTEGER(intKind)              :: status ,date_time(8);
!
!
!version(1:80)  = nf90_inq_libvers();
!!CALL check_error(status);
!version = "NETCDF "//TRIM(version);
!
!
!CALL DATE_AND_TIME(values=date_time);
!WRITE(ctime,101)date_time(1:3),date_time(5:7);
!!---------------------------------------------------------------------------------------------------------------------
!!only need to set the dim name and values.
!!---------------------------------------------------------------------------------------------------------------------
!!!1.Name
!nsf_self%att_name_Array(1) ="conventions";
!nsf_self%att_name_Array(2) ="grid_type";
!nsf_self%att_name_Array(3) ="model";
!nsf_self%att_name_Array(4) ="title";
!nsf_self%att_name_Array(5) ="comment";
!nsf_self%att_name_Array(6) ="source"
!nsf_self%att_name_Array(7) ="institution";
!nsf_self%att_name_Array(8) ="history";
!nsf_self%att_name_Array(9) ="references" ;
!nsf_self%att_name_Array(10) ="creation_date";
!
!!!2.Value
!nsf_self%att_val_Array(1)  =version;
!nsf_self%att_val_Array(2)  ="Triangular";
!nsf_self%att_val_Array(3)  ="SELFE";
!nsf_self%att_val_Array(4)  ="CER";
!nsf_self%att_val_Array(5)  ="testing";
!nsf_self%att_val_Array(6)  ="Fortran script" ;
!nsf_self%att_val_Array(7)  ="PORL/TMSI/NUS";
!nsf_self%att_val_Array(8)  ="original" ;
!nsf_self%att_val_Array(9)  ="XXXXXX";
!nsf_self%att_val_Array(10) =ctime;
!
!
!101 FORMAT(I4,"-",I2.2,"-",I2.2, " ",I2.2, ":", I2.2, ":",I2.2)
!END SUBROUTINE








!SUBROUTINE  nfs_self_set_vars(nsf_self,iof)
!!---------------------------------------------------------------------------------------------------------------------
!!  Modules 
!!---------------------------------------------------------------------------------------------------------------------
!IMPLICIT NONE
!!---------------------------------------------------------------------------------------------------------------------
!!  Arguments     
!!---------------------------------------------------------------------------------------------------------------------
!TYPE(Tnsf_container)         ,INTENT(INOUT) :: nsf_self;
!INTEGER(intKind),DIMENSION(:),INTENT(IN   ) :: iof;
!!---------------------------------------------------------------------------------------------------------------------
!!  Local variables 
!!---------------------------------------------------------------------------------------------------------------------
!INTEGER(intKind)                :: i;
!INTEGER(intKind)                :: iof_index;
!i=0
!!---------------------------------------------------------------------------------------------------------------------
!!Add non-time dependented varialbes.
!!---------------------------------------------------------------------------------------------------------------------
!
!!-------------------------------------------
!!1.ele;
!i = i+1 ;
!nsf_self%var_name_Array(i)       = "ele";
!nsf_self%var_data_type_array(i)  = nsf_int;
!nsf_self%var_nDim_array(i)       = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nface_id,nsf_dim_nele_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nface_id,nsf_dim_nele_id/);
!
!nsf_self%var_natt_Array(i) = 3;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  = (/"long_name","standard_name","units"/)
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) = (/"Horizontal Triangular Element Incidence List", &
!                                                         & "mesh element" ,"index_start_1"/);
!                                                        
!!-------------------------------------------
!!2.bnd;
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "bnd";
!nsf_self%var_data_type_array(i)   =nsf_int;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nbi_val,nsf_dim_nbnd_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nbi_id,nsf_dim_nbnd_id/);
!
!nsf_self%var_natt_Array(i) = 4;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  = (/"long_name","standard_name","units","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Boundary Segment Node List" , &
!                                                       & "mesh element" ,"index_start_1","-9999"/);
!
!!-------------------------------------------
!!3.lon;
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "lon";
!nsf_self%var_data_type_array(i)   =nsf_float;
!nsf_self%var_nDim_array(i)        = 1;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id/);
!
!nsf_self%var_natt_Array(i) = 3;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Nodal Longitude" , "longitude" ,"degrees_east" /);
!
!
!!-------------------------------------------
!!4.lat;
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "lat";
!nsf_self%var_data_type_array(i)   =nsf_float;
!nsf_self%var_nDim_array(i)        = 1;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id/);
!
!nsf_self%var_natt_Array(i) = 3;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Nodal Latitude" , "latitude" ,"degrees_north" /);
!
!
!!-------------------------------------------
!!5.x;
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "x";
!nsf_self%var_data_type_array(i)   =nsf_float;
!nsf_self%var_nDim_array(i)        = 1;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id/);
!
!nsf_self%var_natt_Array(i) = 3;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Nodal x-coordinate" , "x_coordinate" ,"m" /);
!
!!-------------------------------------------
!!6.y;
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "y";
!nsf_self%var_data_type_array(i)   =nsf_float;
!nsf_self%var_nDim_array(i)        = 1;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id/);
!
!nsf_self%var_natt_Array(i) = 3;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Nodal y-coordinate" , "y_coordinate" ,"m" /);
!
!!7.depth
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "depth";
!nsf_self%var_data_type_array(i)   =nsf_float;
!nsf_self%var_nDim_array(i)        = 1;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id/);
!
!nsf_self%var_natt_Array(i) = 4;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Bathymetry" , "still_water_depth" ,"m","down" /);
!
!!8.sigma
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "sigma";
!nsf_self%var_data_type_array(i)   =nsf_float;
!nsf_self%var_nDim_array(i)        = 1;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nslayer_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nslayer_id/);
!
!nsf_self%var_natt_Array(i) = 4;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Sigma Stretched Vertical Coordinate at Nodes" , &
!                                                         & "ocean_sigma_coordinate" ,"sigma_layer","down"/);
!
!!9.zeta
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "zeta";
!nsf_self%var_data_type_array(i)   =nsf_float;
!nsf_self%var_nDim_array(i)        = 1;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nzlayer_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nzlayer_id/);
!
!nsf_self%var_natt_Array(i) = 4;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Zeta Vertical Coordinate at Nodes" , &
!                                                         & "zeta_coordinate" ,"m","down"/);
!                                                        
!!10.time
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "time";
!nsf_self%var_data_type_array(i)   =nsf_float;
!nsf_self%var_nDim_array(i)        = 1;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 4;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","base_date "/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"time" , &
!                                                         & "time" ,"days since 2008-08-22 12:00:00 00:00","2008, 8, 22, 12 "/);
!
!
!!-------------------------------------------
!!add time depend varialbes.
!!-------------------------------------------
!
!
!!11.elev 
!iof_index  =  0;
!IF( iof(1) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "elev";
!nsf_self%var_data_type_array(i)   =nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF
!
!!12.pres
!IF( iof(2) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "pres";
!nsf_self%var_data_type_array(i)   =nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF
! 
!IF( iof(3) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "airt";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF 
!      
!
!IF( iof(4) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "shum";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF 
!
!IF( iof(5) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "srad";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF 
!
!IF( iof(6) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "flsu";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF 
!  
!IF( iof(7) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "fllu";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF    
!
!IF( iof(8) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "radu";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF                                                     
!
!IF( iof(9) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "radd";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF   
!
!IF( iof(10) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "flux";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF   
!
!
!IF( iof(11) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "evap";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF   
!
!
!IF( iof(12) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "prcp";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF   
!
!IF( iof(13) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "windx";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "windy";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!
!ENDIF 
!
!IF( iof(14) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "wistx";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "wisty";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!
!ENDIF 
!
!IF( iof(15) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "dahvx";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "dahvy";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 2;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!
!ENDIF 
! 
!IF( iof(16) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "vert";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF   
!
!IF( iof(17) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "temp";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF      
!
!IF( iof(18) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "salt";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF   
!
!IF( iof(19) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "conc";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF    
!
!IF( iof(20) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "tdff";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF                                          
!          
!IF( iof(21) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "vdff";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF     
!
!IF( iof(22) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "kine";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF       
!
!IF( iof(23) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "mixl";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF   
! 
!IF( iof(24) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "zcor";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF      
!
!IF( iof(25) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "qnon";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF    
!
!IF( iof(26) == 1 ) THEN   
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "hvelx";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!
!i = i+1 ;
!nsf_self%var_name_Array(i)        = "hvely";
!nsf_self%var_data_type_array(i)   = nsf_float;
!nsf_self%var_nDim_array(i)        = 3;
!nsf_self%var_dims_Array(1:nsf_self%var_nDim_array(i),i)     = (/nsf_dim_nlayer_val,nsf_dim_node_val,nsf_dim_time_val/);
!nsf_self%var_dims_id_Array(1:nsf_self%var_nDim_array(i),i)  = (/nsf_dim_nlayer_id,nsf_dim_node_id ,nsf_dim_time_id/);
!
!nsf_self%var_natt_Array(i) = 6;
!nsf_self%var_att_name_Array(1:nsf_self%var_natt_Array(i),i)  =(/"long_name","standard_name","units","positive","location","fill_value"/);
!nsf_self%var_att_val_Array(1:nsf_self%var_natt_Array(i),i) =(/"Water Surface Elevation" , "sea_surface_elevation", &
!                                                         & "m" ,"down","node","-9999"/);
!ENDIF                        
!END SUBROUTINE


END MODULE 

