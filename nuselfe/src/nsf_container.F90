!-----------------------------------------------------------------------------------
!c    COPYRIGHT 2010: prof. Pavel Tkalich(TMSI), Xu Haihua(TMSI), Dao My Ha (TMSI)
!c    1. This subroutine is developed by Xu Haihua(TMSI), Dao My Ha(TMSI) and supervised
!c       by prof. Pavel Tkalich(TMSI)
!     2. TMSI: Tropical Marine Science Institute, Singapore
!c    3. For any enquiry pls contact tmsxh@nus.edu.sg
!-----------------------------------------------------------------------------------

MODULE nsf_container_mod
USE data_type_mod
USE check_error_mod
IMPLICIT NONE;

!----------all the variable will start as nsf (mean ntcdf self form)
INTEGER, PARAMETER  :: NSF_MAXIMUM_ATTS = 20;
INTEGER, PARAMETER  :: NSF_MAXIMUM_DIMS = 20;

INTEGER , PARAMETER ::  nsf_int   = NF90_INT;
INTEGER , PARAMETER ::  nsf_float = NF90_float;

INTEGER(intKind) , PARAMETER ::  nsf_int_default_fill_value  = 0;
REAL(realKind)   , PARAMETER ::  nsf_real_default_fill_value = -9999;

INTEGER , PARAMETER ::  nsf_maximum_att_count =100;
INTEGER , PARAMETER ::  nsf_maximum_dim_count =100;
INTEGER , PARAMETER ::  nsf_maximum_var_count =100;
INTEGER , PARAMETER ::  nsf_maximum_var_dim_count =5;


!----------------------------------------------------------
!NETCDF dimension type
!name : the name of the dimension.
!id   : the Id of the dimension. (if created first then id =1, second id=2..)
!val  : dimension value of the id.
!----------------------------------------------------------

TYPE Tnsf_dimension
    CHARACTER(CharLen)   :: name                 ="";
    INTEGER(intKind)     :: id                   =0;  
    INTEGER(intKind)     :: val                  =0;       
END TYPE


!----------------------------------------------------------
!NETCDF attribute type
!name      : name of the attribute.
!data_type : data type of the attribute.
!val_*     : value of the attribute.
!----------------------------------------------------------
TYPE Tnsf_attribute 

    CHARACTER(CharLen)   :: name            ="";
    INTEGER(intKind)     :: data_type       =0 ;
    CHARACTER(CharLen)   :: val_char        ="";
    INTEGER(2)           :: val_int_short   =0 ;    
    INTEGER(4)           :: val_int_long    =0 ;
    REAL(4)              :: val_real_float  =0.0;     
    REAL(8)              :: val_real_double =0.0; 
           
END TYPE



!----------------------------------------------------------
!NETCDF variable type
!name      : name of the variable.
!data_type : data type of the variable.
!id        : id of the variable.
!dims_count: number of dimensions of the varialbe.
!dims_ids  : dimension id of each dimension.
!dims      : lenght of each dimension.
!att_count : number of attribute.
!att       : attribute.
!----------------------------------------------------------
TYPE Tnsf_variable
    CHARACTER(CharLen)   :: name                       ="";
    INTEGER(intKind)     :: data_type                  =0; 
    INTEGER(intKind)     :: id                         =0;   
    INTEGER(intKind)     :: dims_count                 =0;  
    INTEGER(intKind)     :: dims_ids(NSF_MAXIMUM_DIMS) =0;    
    INTEGER(intKind)     :: dims(NSF_MAXIMUM_DIMS)     =0;        
    INTEGER(intKind)     :: atts_count                  =0;  
    TYPE(Tnsf_attribute) :: atts(NSF_MAXIMUM_ATTS)      ;
END TYPE


TYPE Tnsf_container
!----------------------------------------------------------
! ncid     : ID of the netcdf file.
! file_path: the file path of the netcdf file.
! file_mode: the file mode of the netcdf file. (NF90_NOCLOBBER, NF90_SHARE, and NF90_64BIT_OFFSET)
!----------------------------------------------------------
INTEGER(intKind)                           ::  ncid      =0;
CHARACTER(CharLen)                         ::  file_path ="";
INTEGER(intKind)                           ::  file_mode =-1;

!----------------------------------------------------------
!SELF NETCDF FILE GLOBAL ATTRIBUT DATA,
! att_conventions_val, att_creation_date_val : need to be
!              get at the run time.
!NOTE: use the function: nf90_put_att (NF90_GLobal)  to set the global attribute
!            ;
!----------------------------------------------------------
INTEGER(intKind)                            ::  atts_count    =0;
TYPE(Tnsf_attribute)   ,DIMENSION(:),ALLOCATABLE::  atts        ;             

!----------------------------------
!NETCDF dimension part.
INTEGER(intKind)                            ::  dims_count    =0;
TYPE(Tnsf_dimension)   ,DIMENSION(:),ALLOCATABLE::  dims 

!---------------------
!NETCDF variables part.
INTEGER(intKind)                            :: vars_count     =0;
TYPE(Tnsf_variable)    ,DIMENSION(:),ALLOCATABLE:: vars 


INTEGER(intKind)                            :: unlimitedDimId =0;

END TYPE 


INTERFACE nsf_put_attribute
    MODULE PROCEDURE nsf_set_attribute_netcdf;
END INTERFACE

INTERFACE nsf_set_attribute
    MODULE PROCEDURE nsf_set_attribute_netcdf;
    MODULE PROCEDURE nsf_set_attribute_char
    MODULE PROCEDURE nsf_set_attribute_int2
    MODULE PROCEDURE nsf_set_attribute_int4
    MODULE PROCEDURE nsf_set_attribute_real4
    MODULE PROCEDURE nsf_set_attribute_real8
END INTERFACE                


INTERFACE nsf_get_attribute
    MODULE PROCEDURE nsf_get_attribute_netcdf;
    MODULE PROCEDURE nsf_get_attribute_char
    MODULE PROCEDURE nsf_get_attribute_int2
    MODULE PROCEDURE nsf_get_attribute_int4
    MODULE PROCEDURE nsf_get_attribute_real4
    MODULE PROCEDURE nsf_get_attribute_real8
END INTERFACE   





CONTAINS



!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   init_nsf_container
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : allocate the memory for nsf_vars
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
SUBROUTINE init_nsf_container(container,n_global_atts,n_dims,n_vars)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_container)         ,INTENT(INOUT) :: container;
INTEGER(intKind),OPTIONAL    ,INTENT(IN   ) :: n_global_atts,n_dims,n_vars;
!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                   :: status;
CHARACTER(CharLen)                 :: cSub ="init_nsf_container";

!-------------init the arrary-------
IF(PRESENT(n_global_atts)) THEN
    container%atts_count  = n_global_atts;
ELSE
    container%atts_count  = nsf_maximum_att_count;
ENDIF

IF(PRESENT(n_dims)) THEN
    container%dims_count  = n_dims;
ELSE
    container%dims_count  = nsf_maximum_dim_count;
ENDIF

IF(PRESENT(n_vars)) THEN
    container%vars_count  = n_vars;
ELSE
    container%vars_count  = nsf_maximum_var_count;
ENDIF

ALLOCATE(container%atts(container%atts_count) , &
         container%dims(container%dims_count) , &
         container%vars(container%vars_count) , &
         STAT = status);
CALL check_results(status,cSub,"Failed to allocate array");         
         
END SUBROUTINE


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nfs_create_file
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
SUBROUTINE nfs_create_file(container,path,mode)
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
TYPE(Tnsf_container)                :: container;
CHARACTER(*)        , INTENT (IN   ):: path;
INTEGER(intKind)    , INTENT (IN   ):: mode;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                 :: i,j,status;
INTEGER(intKind)                 :: varid,old_mode;
CHARACTER(CharLen)               :: temp_char;
INTEGER(2)                       :: temp_int2;
INTEGER(4)                       :: temp_int4;
REAL(4)                          :: temp_real4;
REAL(8)                          :: temp_real8;
!---------------------------------------------------------------------------------------------------------------------
!1.create the file.
!---------------------------------------------------------------------------------------------------------------------
container%file_path = path ;
container%file_mode = mode ;

status = nf90_create(path = container%file_path, cmode = container%file_mode, ncid = container%ncid);
CALL check_error(status);

!---------------------------------------------------------------------------------------------------------------------
!2.assign the global attribute.
!---------------------------------------------------------------------------------------------------------------------
! 1.conventions

DO i = 1,container%atts_count
    CALL nsf_put_attribute(container%ncid,NF90_GLOBAL,container%atts(i));
ENDDO

!---------------------------------------------------------------------------------------------------------------------
!3. create dimesnions in netcdf
!---------------------------------------------------------------------------------------------------------------------
DO i = 1,container%dims_count

    status = nf90_def_dim(container%ncid,TRIM(container%dims(i)%name),container%dims(i)%val, container%dims(i)%id);
    CALL check_error(status);

ENDDO


!---------------------------------------------------------------------------------------------------------------------
!4.create varialbes and attributes.
!---------------------------------------------------------------------------------------------------------------------
DO i = 1,container%vars_count

    CALL nfs_create_variable(container%ncid,container%vars(i))

ENDDO
!-------------------------------------------------------
!set to be fill model.
!-------------------------------------------------------
!status = nf90_set_fill(container%ncid, NF90_FILL, old_mode)
!
!------------------------------------------------------
!5.end define and ready to put the data.
!--------------------------------------------------
status = nf90_enddef(container%ncid);
CALL check_error(status);

END SUBROUTINE


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nfs_close_file
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
SUBROUTINE nfs_close_file(container)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_container)        :: container;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                 :: i,status;

!------de_allocate the memory........
IF(ALLOCATED(container%atts)) DEALLOCATE(container%atts);
IF(ALLOCATED(container%dims)) DEALLOCATE(container%dims);
IF(ALLOCATED(container%vars)) DEALLOCATE(container%vars);

!------close the file........

status = nf90_close(container%ncid);
CALL check_error(status);

END SUBROUTINE


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nfs_open_file
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
SUBROUTINE nfs_open_file(container,path,mode)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_container)                 :: container;
CHARACTER(*)         , INTENT (IN   ):: path;
INTEGER(intKind)     , INTENT (IN   ):: mode;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                 :: i,j,status;
INTEGER(intKind)                 :: varid;

INTEGER(intKind)                 :: xtype, len, attnum;

CHARACTER(CharLen)               :: att_char;
INTEGER(2)                       :: att_int2
INTEGER(4)                       :: att_int4;
REAL(4)                          :: att_float;
REAL(8)                          :: att_double;
!---------------------------------------------------------------------------------------------------------------------
!1.open the file.
!---------------------------------------------------------------------------------------------------------------------
container%file_path = path ;
container%file_mode = mode ;

status = nf90_open(path = container%file_path, mode = container%file_mode, ncid = container%ncid);
CALL check_error(status);

!WRITE(*,*) "Open success.."
!---------------------------------------------------------------------------------------------------------------------
!2. inquiry the file..
!---------------------------------------------------------------------------------------------------------------------
status = nf90_Inquire(ncid = container%ncid, nDimensions = container%dims_count, &
                 &    nVariables = container%vars_count, nAttributes =container%atts_count,unlimitedDimId = container%unlimitedDimId);

!WRITE(*,*) "Onf90_Inquire success.."
!---------------------------------------------------------------------------------------------------------------------
!3. allocate the variables.
!---------------------------------------------------------------------------------------------------------------------
!CALL alloc_nsf_vars(container%att_count,container%dim_count,container%var_count);
CALL init_nsf_container(container,container%atts_count,container%dims_count,container%vars_count);

!WRITE(*,*) "init_nsf_container success.."
!---------------------------------------------------------------------------------------------------------------------
!4. get the global attribute Name and value.
!---------------------------------------------------------------------------------------------------------------------
DO i = 1, container%atts_count
    status = nf90_inq_attname(container%ncid, NF90_GLOBAL, i, container%atts(i)%name);
    CALL check_error(status);
    CALL nsf_get_attribute(container%ncid,NF90_GLOBAL,container%atts(i));
ENDDO

!WRITE(*,*) "nf90_inq_attname success.."

!---------------------------------------------------------------------------------------------------------------------
!5. get the dimension.
!---------------------------------------------------------------------------------------------------------------------
DO i = 1, container%dims_count

    container%dims(i)%id  = i;
    status = nf90_Inquire_Dimension(container%ncid, container%dims(i)%id ,&
                   & container%dims(i)%name,container%dims(i)%val) ;
    
ENDDO

!WRITE(*,*) "nf90_Inquire_Dimension success.."

!---------------------------------------------------------------------------------------------------------------------
!6. get the variables name id, data types and so on 
!---------------------------------------------------------------------------------------------------------------------
DO  i = 1, container%vars_count

    container%vars(i)%id  = i;

    status = nf90_Inquire_Variable(container%ncid,  container%vars(i)%id ,&
                    container%vars(i)%name,container%vars(i)%data_type, &
                    container%vars(i)%dims_count,container%vars(i)%dims_ids, &
                    container%vars(i)%atts_count);
    CALL check_error(status);                   

!---------get the length of each dimension;
    DO j = 1, container%vars(i)%dims_count
      container%vars(i)%dims(j) =  container%dims(container%vars(i)%dims_ids(j))%val;
    ENDDO 
    
!---------get the attributes of each dimension;
    DO j = 1,container%vars(i)%atts_count; 
    
        status = nf90_inq_attname(container%ncid,container%vars(i)%id,j, container%vars(i)%atts(j)%name);
        CALL check_error(status);
 
        CALL nsf_get_attribute(container%ncid,container%vars(i)%id,container%vars(i)%atts(j));
    END DO

                   
ENDDO

!WRITE(*,*) "nf90_Inquire_Variable success.."

CALL check_error(status);
END SUBROUTINE

!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nfs_create_variable
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
SUBROUTINE nfs_create_variable(ncid,var)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                 , INTENT(IN   )::ncid;
TYPE(Tnsf_variable)              , INTENT(INOUT)::var;
!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                 :: i,status;
CHARACTER(CharLen)               :: temp_char;
INTEGER(2)                       :: temp_int2;
INTEGER(4)                       :: temp_int4;
REAL(4)                          :: temp_real4;
REAL(8)                          :: temp_real8;

!---------------------------------------------------------------------------------------------------------------------
!define the array
!---------------------------------------------------------------------------------------------------------------------
status = nf90_def_var(ncid,TRIM(var%name), var%data_type, &
                      var%dims_ids(1:var%dims_count),var%id);
CALL check_error(status);

!---------------------------------------------------------------------------------------------------------------------
!define the attribute array
!---------------------------------------------------------------------------------------------------------------------
DO i = 1, var%atts_count
    CALL nsf_put_attribute(ncid,var%id,var%atts(i));
ENDDO 


END SUBROUTINE



!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_set_attribute_netcdf
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the attribute by netcdf function.
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
SUBROUTINE nsf_set_attribute_netcdf(ncid,varid,att)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)    ,INTENT(IN   ) :: ncid,varid;
TYPE(Tnsf_attribute),INTENT(INOUT) :: att;


!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                 :: status
CHARACTER(CharLen)               :: temp_char;
INTEGER(2)                       :: temp_int2;
INTEGER(4)                       :: temp_int4;
REAL(4)                          :: temp_real4;
REAL(8)                          :: temp_real8;

temp_char ="";
IF(att%data_type == NF90_CHAR) THEN
    CALL nsf_get_attribute(att,temp_char);
    status = nf90_put_att(ncid,varid,TRIM(att%name),TRIM(temp_char));
    CALL check_error(status);
ELSEIF(att%data_type == NF90_SHORT) THEN
    CALL nsf_get_attribute(att,temp_int2);
    status = nf90_put_att(ncid,varid,TRIM(att%name), temp_int2);
    CALL check_error(status);    
ELSEIF(att%data_type == NF90_INT) THEN
    CALL nsf_get_attribute(att,temp_int4);
    status = nf90_put_att(ncid,varid,TRIM(att%name), temp_int4);
    CALL check_error(status);    
ELSEIF(att%data_type == NF90_FLOAT) THEN
    CALL nsf_get_attribute(att,temp_real4);
    status = nf90_put_att(ncid,varid,TRIM(att%name), temp_real4);
    CALL check_error(status);    
ELSEIF(att%data_type == NF90_DOUBLE) THEN
    CALL nsf_get_attribute(att,temp_real8);
    status = nf90_put_att(ncid,varid,TRIM(att%name), temp_real8);
    CALL check_error(status);    
ELSE
    WRITE(*,*) "Attribute type not defined.",att%data_type;
    status = sph_error_error;
ENDIF

    CALL check_error(status);    

END SUBROUTINE


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_set_attribute_char
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the character value of the attribute.
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
SUBROUTINE nsf_set_attribute_char(att,val)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_attribute),INTENT(INOUT) :: att;
CHARACTER(*)        ,INTENT(IN   ) :: val;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
att%data_type = NF90_CHAR;
att%val_char  = TRIM(val);

END SUBROUTINE


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_set_attribute_int2
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the character value of the attribute.
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
SUBROUTINE nsf_set_attribute_int2(att,val)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_attribute),INTENT(INOUT) :: att;
INTEGER(2)          ,INTENT(IN   ) :: val;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
att%data_type = NF90_SHORT;
att%val_int_short  = val;

END SUBROUTINE


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_set_attribute_int4
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the character value of the attribute.
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
SUBROUTINE nsf_set_attribute_int4(att,val)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_attribute),INTENT(INOUT) :: att;
INTEGER(4)          ,INTENT(IN   ) :: val;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
att%data_type     = NF90_INT;
att%val_int_long  = val;

END SUBROUTINE


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_set_attribute_real4
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the character value of the attribute.
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
SUBROUTINE nsf_set_attribute_real4(att,val)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_attribute),INTENT(INOUT) :: att;
REAL(4)             ,INTENT(IN   ) :: val;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
att%data_type = NF90_FLOAT;
att%val_real_float  = val;

END SUBROUTINE


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_set_attribute_real8
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the character value of the attribute.
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
SUBROUTINE nsf_set_attribute_real8(att,val)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_attribute),INTENT(INOUT) :: att;
REAL(8)             ,INTENT(IN   ) :: val;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
att%data_type       = NF90_DOUBLE;
att%val_real_double = val;


END SUBROUTINE

!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_get_attribute_netcdf
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the attribute by netcdf function.
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
SUBROUTINE nsf_get_attribute_netcdf(ncid,varid,att)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)    ,INTENT(IN   ) :: ncid,varid;
TYPE(Tnsf_attribute),INTENT(INOUT) :: att;


!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)                 :: status
CHARACTER(CharLen)               :: temp_char;
INTEGER(2)                       :: temp_int2;
INTEGER(4)                       :: temp_int4;
REAL(4)                          :: temp_real4;
REAL(8)                          :: temp_real8;

temp_char ="";

status=nf90_inquire_attribute(ncid,varid, TRIM(att%name), xtype = att%data_type);

CALL check_error(status);    
 
IF(att%data_type == NF90_CHAR) THEN
   status = nf90_get_att(ncid ,varid, TRIM(att%name),temp_char); 
   CALL check_error(status); 
   CALL nsf_set_attribute(att, TRIM(temp_char)); 

ELSEIF(att%data_type == NF90_SHORT) THEN
   status = nf90_get_att(ncid ,varid, TRIM(att%name),temp_int2); 
   CALL check_error(status); 
   CALL nsf_set_attribute(att, temp_int2); 
    
ELSEIF(att%data_type == NF90_INT) THEN
   status = nf90_get_att(ncid ,varid, TRIM(att%name),temp_int4); 
   CALL check_error(status); 
   CALL nsf_set_attribute(att, temp_int4); 
   
ELSEIF(att%data_type == NF90_FLOAT) THEN
   status = nf90_get_att(ncid ,varid, TRIM(att%name),temp_real4); 
   CALL check_error(status); 
   CALL nsf_set_attribute(att, temp_real4); 
   
ELSEIF(att%data_type == NF90_DOUBLE) THEN
   status = nf90_get_att(ncid ,varid, TRIM(att%name),temp_real8); 
   CALL check_error(status); 
   CALL nsf_set_attribute(att, temp_real8); 
      
ELSE
    WRITE(*,*) "Attribute type not defined.",att%data_type;
    status = sph_error_error;
ENDIF

    CALL check_error(status);    

END SUBROUTINE




!---------------------------------------------------------------------------------------------------------------------
!get the var id of selfe data.
!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE nsf_get_var_id(container,varName,varId)
!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_container),INTENT(IN ) :: container
CHARACTER(*)      ,INTENT(IN   ) :: varName;
INTEGER(intKind)  ,INTENT(  OUT) :: varId  ;
!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(intKind)      :: i,status;

!varid = -1;
!DO i = 1,container%vars_count
!    IF(TRIM(container%vars(i)%name)==TRIM(varName) )THEN
!        varid = i;
!    ENDIF
!ENDDO

status = nf90_inq_varid(container%ncid,TRIM(varName), varid)
IF(status /=0) THEN
    varid = -1;
    WRITE(*,*) "Failed to inquire the variable Id for var :",TRIM(varName);
ENDIF



END SUBROUTINE

!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_get_attribute_char
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the character value of the attribute.
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
SUBROUTINE nsf_get_attribute_char(att,val)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_attribute),INTENT(IN   ) :: att;
CHARACTER(*)        ,INTENT(  OUT) :: val;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
val = att%val_char  ;

END SUBROUTINE


!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_get_attribute_char
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the character value of the attribute.
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
SUBROUTINE nsf_get_attribute_int2(att,val)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_attribute),INTENT(IN   ) :: att;
INTEGER(2)          ,INTENT(  OUT) :: val;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
val = att%val_int_short ;

END SUBROUTINE



!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_get_attribute_char
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the character value of the attribute.
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
SUBROUTINE nsf_get_attribute_int4(att,val)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_attribute),INTENT(INOUT) :: att;
INTEGER(4)          ,INTENT(  OUT) :: val;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
val = att%val_int_long ;

END SUBROUTINE



!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_get_attribute_char
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the character value of the attribute.
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
SUBROUTINE nsf_get_attribute_real4(att,val)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_attribute),INTENT(IN   ) :: att;
REAL(4)             ,INTENT(  OUT) :: val;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
val= att%val_real_float ;

END SUBROUTINE



!------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   nsf_get_attribute_char
!------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : set the character value of the attribute.
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
SUBROUTINE nsf_get_attribute_real8(att,val)

!---------------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
TYPE(Tnsf_attribute),INTENT(INOUT) :: att;
REAL(8)             ,INTENT(  OUT) :: val;

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
val = att%val_real_double;


END  SUBROUTINE


END MODULE