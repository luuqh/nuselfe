MODULE post_vars
USE elfe_glbl
IMPLICIT NONE

!------------------------------------------------------
!In windows onbyte == 1;
!In Linux  onbyte == 4;
!------------------------------------------------------
INTEGER, PARAMETER        :: onbyte  = 1;

!------------------------------------------------------
!orkind: SELFE out put real type kind default real(4)
!------------------------------------------------------
INTEGER,PARAMETER         :: orkind = 4;

!------------------------------------------------------
!file63  :file name need to process. which contain the name and ext
!         name specify the varialbe and ext specify the var type (scalar or vector, 2D or 3D data)
!lfile63 :length of the file name. 
!------------------------------------------------------
CHARACTER(LEN=200)        :: file63,inFolder,outFolder;
INTEGER                   :: lfile63;
!------------------------------------------------------
!ibgn,iend: begin day and end day of the simulation.
!ivs      : data type (1: scalar, 2: vector)
!i23d     : 2d or 3d data (2d: mean no vertical distribution: elevation, 3d: mean 
!------------------------------------------------------
INTEGER                   :: ibgn,iend;
INTEGER                   :: ivs  ! variable type =1 scalar, =2 vector, 
INTEGER                   :: i23d ! data dimesnion 2d or 3d;

!------------------------------------------------------
!number of proc used for the SELFE simulation.
!------------------------------------------------------
INTEGER                   :: nproc  ! number of proc on the file.


!------------------------------------------------------
!iplg, np,ne,nm,ielg on each proc
!------------------------------------------------------
INTEGER,ALLOCATABLE, DIMENSION(:,:):: iplg_local,nm_local,ielg_local;
INTEGER,ALLOCATABLE, DIMENSION(:,:):: kbp00_local;
INTEGER,ALLOCATABLE, DIMENSION(:  ):: np_local,ne_local

!------------------------------------------------------
!htot: store the total water detph;
!delz: store the layer thickness.
!------------------------------------------------------
REAL(rkind),ALLOCATABLE,DIMENSION(:  ):: htot;
REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: delz;


!------------------------------------------------------
!outb    : the data read from the file.
!dtout   : 
!out_proc: the data after processing, data extract from
!             the bottom layer,top layer or the averaged.
!------------------------------------------------------
REAL(rkind),ALLOCATABLE,DIMENSION(:,:,:),TARGET :: outb;
REAL(rkind)                                     :: dtout
REAL(rkind),ALLOCATABLE,DIMENSION(:,:)          :: out_proc;
REAL(rkind),ALLOCATABLE,DIMENSION(:,:)          :: tec_out

REAL(rkind)                                     :: time;
INTEGER                                         :: it;

!------------------------------------------------------
!out put data to tecplot SMS or TXT other type
!------------------------------------------------------
REAL(rkind),ALLOCATABLE,DIMENSION(:)            :: xx_out;
REAL(rkind),ALLOCATABLE,DIMENSION(:,:)          :: yy_out;
REAL(rkind),ALLOCATABLE,DIMENSION(:,:,:)        :: data_out;


INTEGER,PARAMETER                     :: SELFE_VERTICAL_AVERAGE  = -2;
INTEGER,PARAMETER                     :: SELFE_BOTTOM_LAYER      = -1;
INTEGER,PARAMETER                     :: SELFE_SURFACE_LAYER     =  0;


INTEGER                               :: mark_dry  =1;
INTEGER                               :: ilev      =0;
REAL(rkind)                           :: dry       =-9999.0;
REAL(rkind)                           :: air       =-9998.0;
REAL(rkind)                           :: soil      =-9997.0;


INTEGER                               :: nout_node_list;
INTEGER,ALLOCATABLE,DIMENSION(:)      :: iout_node_list;



!--------------------------------------------------------------------
!Nodes in vertical profile..
!--------------------------------------------------------------------
INTEGER                                 ::n_vprof_node  ;
REAL(rkind),DIMENSION(:  ),ALLOCATABLE  ::xvert;
REAL(rkind),DIMENSION(:,:),ALLOCATABLE  ::zvert;

!--------------------------------------------------------------------------
!disvert: distance of the node from the first node of the vertical profile.
!--------------------------------------------------------------------------
REAL(rkind),DIMENSION(:  ),ALLOCATABLE  ::disvert;

!--------------------------------------------------------------------------
!disvert_hvel_mag : magnitude of horizontal velocity.
!disvert_vvel     : magnitude of vertical velocity.
!--------------------------------------------------------------------------
REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE ::disvert_vel;

!i_mode_combine   :: combine the results.
!i_mode_horz_comb :: output the horziontal distribution from the combined results.
!i_mode_vprof_comb:: output the vertical profile distribution from the combined results.
!i_mode_pTime_comb:: output the point time series from the combined results.

!i_mode_horz_netcdf:: output the horziontal distribution from the netcdf results.
!i_mode_vprof_netcdf:: output the vertical profile distribution from the netcdf results.
!i_mode_pTime_netcdf:: output the point time series from the netcdf results.

INTEGER,PARAMETER                        :: i_mode_combine      = 0;
INTEGER,PARAMETER                        :: i_mode_horz_comb    = 1;
INTEGER,PARAMETER                        :: i_mode_vprof_comb   = 2;
INTEGER,PARAMETER                        :: i_mode_pTime_comb   = 3;
INTEGER,PARAMETER                        :: i_mode_horz_netcdf  = 4;
INTEGER,PARAMETER                        :: i_mode_vprof_netcdf = 5;
INTEGER,PARAMETER                        :: i_mode_pTime_netcdf = 6;

INTEGER                                  :: i_mode = i_mode_combine ;


INTEGER,PARAMETER                        :: i_out_form_TECPLOT = 1;
INTEGER,PARAMETER                        :: i_out_form_SMS     = 2;
INTEGER,PARAMETER                        :: i_out_form_TXT     = 3;
INTEGER                                  :: i_out_form         = i_out_form_SMS;;


INTEGER,PARAMETER                        :: tec_form_sep   = 1;
INTEGER,PARAMETER                        :: tec_form_comb  = 2;
INTEGER                                  :: tec_form = tec_form_sep;


INTEGER,PARAMETER                        :: i_vprof_scalar   = 1;
INTEGER,PARAMETER                        :: i_vprof_velocity = 2;
INTEGER                                  :: i_vprof_mode = tec_form_sep;
END MODULE