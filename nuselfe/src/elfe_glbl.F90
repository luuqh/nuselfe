!===============================================================================
!===============================================================================
! SELFE GLOBAL DATA MODULE
!===============================================================================
!===============================================================================
module elfe_glbl
  implicit none
  public  !Default scope is public

#ifdef USE_SINGLE
  integer,parameter :: rkind = 4  
#else
  integer,parameter :: rkind = 8      ! Default real datatype
#endif

  ! ADT for global-to-local linked-lists
  type :: llist_type
    integer                      :: rank      ! Processor rank assignment
    integer                      :: id=0      ! Local index on processor "rank"
    type(llist_type),pointer :: next=>null()  ! Next entry in linked-list
  end type llist_type

  ! ADT for inter-subdomain backtracking
  type :: bt_type
    integer :: rank          ! Originating processor rank
    integer :: l0            ! Originating node or side or element index (1<=l0<=3)
    integer :: i0gb          ! Originating node or side or element _global_ index 
    integer :: isbndy        ! Flag for originating node or side on the boundary (for Kriging)
    integer :: j0            ! Originating vertical level
    integer :: adv            ! Original advection flag (0-2)
    integer :: ndt           ! Number of backtracking sub-steps fro Euler tracking
    integer :: iegb          ! Global index of current encompassing element
    integer :: jvrt          ! Current vertical level
    real(rkind) :: dtbm      ! max. step used in R-K tracking
    real(rkind) :: vis      ! vis_coe
    real(rkind) :: rt        ! time remaining
    real(rkind) :: ut,vt,wt  ! Current backtracking sub-step velocity
    real(rkind) :: xt,yt,zt  ! Current backtracking sub-step point
    real(rkind) :: tt,st     ! Backtracked temperature & salinity
  end type bt_type
  integer,save :: bt_mpitype ! MPI datatype for inter-subdomain backtracking
  real(rkind),save :: s1_mxnbt      ! scale used in mxnbt fro btlist()
  real(rkind),save :: s2_mxnbt     ! scale used in the routine inter_btrack()
  integer,save :: mxnbt      ! Dimension of btlist()

  ! Vertical layer data
  integer,save :: ivcor                   ! SZ coordinates; use for sflux routines
  integer,save :: nvrt                    ! Number of vertical layers
  integer,save :: kz                      ! Number of Z levels
  integer,save :: nsig                     ! Number of S levels
  real(rkind),save :: h_s,h_c,theta_b,theta_f,s_con1 !constants used in vgrid.in
  real(rkind),save,allocatable :: ztot(:) ! Z coord. of Z levels
  real(rkind),save,allocatable :: sigma(:) ! sigma coordinates
  real(rkind),save,allocatable :: cs(:) ! function in S-coordinate
  real(rkind),save,allocatable :: dcs(:) ! derivative of cs()

  ! Element geometry data
  integer,save :: ne_global    ! Global number of elements
  integer,save :: ne           ! Local number of resident elements
  integer,save :: neg          ! Local number of ghost elements
  integer,save :: nea          ! Local number of elements in augmented subdomain (ne+neg)
  integer,save,allocatable :: ielg(:)      ! Local-to-global element index table (augmented)
  type(llist_type),save,pointer :: iegl(:) ! Global-to-local element index table (augmented)
  integer,save,allocatable :: iegrpv(:)    ! Global element to resident processor vector
  integer,save :: nx(3,2)                       ! Cyclic index of nodes in an element
  integer,save,allocatable :: nm(:,:)      ! Element-node tables,i34
  integer,save,allocatable :: nmgb(:,:)  ! Global element-node tables;
  integer,save,allocatable :: iself(:,:)          ! Index of node in element-node table
  integer,save,allocatable :: ic3(:,:)            ! Element-side-element tables
  integer,save,allocatable :: ic3gb(:,:)          ! Global element-side-element table
  integer,save,allocatable :: js(:,:)             ! Element-side tables
  real(rkind),save,allocatable :: ssign(:,:)      ! Sign associated with each side of an element
  real(rkind),save,allocatable :: area(:)        ! Element areas
  real(rkind),save,allocatable :: radiel(:)       ! Element equivalent radii
  real(rkind),save,allocatable :: xctr(:),yctr(:) ! Cartesian coordinates of element centers
  real(rkind),save,allocatable :: xlon_el(:),ylat_el(:) ! Element center lat/lon coordinates in degrees
  real(rkind),save,allocatable :: dpe(:)          ! Depth at element centers
  integer,save,allocatable :: kbe(:)       ! Element bottom vertical indices
  integer,save,allocatable :: idry_e(:),idry_e0(:)       ! wet/dry flag
  integer,save,allocatable :: interpol(:)       ! interpolation mode
  integer,save,allocatable :: lqk(:)       ! interpolation for S,T in btrack
  integer,save,allocatable :: ie_kr(:)       ! used in Kriging
  integer,save,allocatable :: krvel(:)       ! used in Kriging
  real(rkind),save,allocatable :: ze(:,:)         ! z-coord
  real(rkind),save,allocatable :: dl(:,:,:)         ! derivatives of shape function
!  real(rkind),save,allocatable :: dc(:,:)         ! Layer thickness at each element
!  real(rkind),save,allocatable :: dchalf(:,:)     ! Half-layer thickness at each element

  ! Node geometry data
  integer,save :: np_global    ! Global number of nodes
  integer,save :: np           ! Local number of resident nodes
  integer,save :: npg          ! Local number of ghost nodes
  integer,save :: npa          ! Local number of nodes in augmented subdomain (np+npg)
  integer,save,allocatable :: iplg(:)      ! Local-to-global node index table (augmented)
  type(llist_type),save,pointer :: ipgl(:) ! Global-to-local node index table (augmented)
  real(rkind),save,allocatable :: x(:),y(:)       ! Node cartesian coordinates
  real(rkind),save,allocatable :: xlon(:),ylat(:) ! Node lat/lon coordinates in radians
  real(rkind),save,allocatable :: dp(:),dp00(:)           ! Node depths
!  integer,save,allocatable :: ibad(:)             ! Reliable bndry elevation flag
  integer,save,allocatable :: nnegb(:),inegb(:,:) ! Global node-element tables
  integer,save,allocatable :: nne(:),ine(:,:)     ! Node-element tables
  integer,save,allocatable :: nnp(:),inp(:,:)     ! Node-node tables
  integer,save,allocatable :: isbnd(:,:)          ! local node to _global_ open bndry segment flags
  integer,save,allocatable :: isbnd_global(:) ! Node to open bndry segment flags (global)
  integer,save,allocatable :: kfp(:),kbp(:),kbp00(:),kbp_e(:) ! Node surface & bottom vertical indices; kfp used only for sflux routines
  integer,save,allocatable :: idry(:)        ! wet/dry flag
  integer,save,allocatable :: iback(:)        ! back-up flag for abnormal cases in S-coord.
  real(rkind),save,allocatable :: hmod(:)        ! constrained depth
  real(rkind),save,allocatable :: z(:,:)        ! z-coord
!  real(rkind),save,allocatable :: dzp(:,:)        ! Layer thickness at each node
!  real(rkind),save,allocatable :: dzphalf(:,:)    ! Half-layer thickness at each node

  ! Side geometry data
  integer,save :: ns_global    ! Global number of sides
  integer,save :: ns           ! Local number of resident sides
  integer,save :: nsg          ! Local number of ghost sides
  integer,save :: nsa          ! Local number of sides in augmented subdomain (ns+nsg)
  integer,save,allocatable :: islg(:)      ! Local-to-global side index table (augmented)
  type(llist_type),save,pointer :: isgl(:) ! Global-to-local side index table (augmented)
  integer,save,allocatable :: is(:,:)             ! Side-element tables
  integer,save,allocatable :: isidenode(:,:)      ! Side-node tables
  real(rkind),save,allocatable :: xcj(:),ycj(:)   ! Cartesian coordinates of side centers
  real(rkind),save,allocatable :: dps(:)          ! Depth at side centers
  real(rkind),save,allocatable :: distj(:)        ! Side lengths
  ! Distance between adjacent elements of an internal side; used only in horizontal diffusion
  real(rkind),save,allocatable :: delj(:)        
  real(rkind),save,allocatable :: snx(:)          ! Cosine of local orientation angle for sides
  real(rkind),save,allocatable :: sny(:)          ! Sine of local orientation angle for sides
  integer,save,allocatable :: isbs(:)           ! local side to _global_ open bndry segment mapping
  integer,save,allocatable :: isbs_global(:)    ! Side to open bndry segment mapping (global)
  integer,save,allocatable :: kbs(:)       ! Side bottom vertical indices
  integer,save,allocatable :: idry_s(:)        ! wet/dry flag
  integer,save,allocatable :: isidenei(:,:),isidenei2(:,:)        !side neighborhood 
  real(rkind),save,allocatable :: zs(:,:)         ! z-coord.
  real(rkind),save,allocatable :: side_ac(:,:,:)         ! used in horizontal viscosity
  real(rkind),save,allocatable :: side_x(:,:)         ! used in horizontal viscosity

  ! Open boundary segment data
  integer,save :: nope_global                  ! Global number of local open bndry segments
  integer,save :: neta_global                  ! Global number of local open bndry nodes
  integer,save :: nope                         ! Local number of local open bndry segments
  integer,save :: neta                         ! Local number of local open bndry nodes
  integer,save :: mnond                        ! Max # nodes per open bndry segment
  integer,save :: mnond_global                 ! Max # nodes per open bndry segment (global)
  integer,save,allocatable :: iopelg(:)        ! Local-to-global open bndry segment table
  integer,save,allocatable :: iopegl(:,:)      ! Global-to-Local open bndry segment table
  integer,save,allocatable :: nond(:)          ! Number of nodes in each open bndry segment
  integer,save,allocatable :: iond(:,:)        ! Node list for each open bndry segment
  integer,save,allocatable :: nond_global(:)   ! Number of nodes in each open global bndry segment
  integer,save,allocatable :: iond_global(:,:) ! Node list for each open bndry segment (global)
  real(rkind),save,allocatable :: cwidth(:)    ! length of each global open bnd segment for imposing discharge

  ! Land boundary segment data
  integer,save :: nland_global                 ! Global number of local land bndry segments
  integer,save :: nvel_global                  ! Global number of local land bndry nodes
  integer,save :: nland                        ! Local number of local land bndry segments
  integer,save :: nvel                         ! Local number of local land bndry nodes
  integer,save :: mnlnd                        ! Max # nodes per land bndry segment
  integer,save :: mnlnd_global                 ! Max # nodes per land bndry segment (global)
  integer,save,allocatable :: nlnd_global(:)   ! Number of nodes in each land bndry segment (global)
  integer,save,allocatable :: ilnd_global(:,:) ! Node list for each land bndry segment (global)
  integer,save,allocatable :: nlnd(:)          ! Number of nodes in each land bndry segment
  integer,save,allocatable :: ilnd(:,:)        ! Node list for each land bndry segment

  ! Dynamical quantities
  real(rkind),save,allocatable :: tsel(:,:,:) ! S,T at elements and half levels for upwind & TVD scheme
  real(rkind),save,allocatable :: tr_el(:,:,:) !tracer converntration @ prism center; used as temp. storage 
  real(rkind),save,allocatable :: bdy_frc(:,:,:) !body force at prism center Q_{i,k}
  real(rkind),save,allocatable :: flx_sf(:,:) !surface b.c. \kappa*dC/dz = flx_sf (at element center)
  real(rkind),save,allocatable :: flx_bt(:,:) !bottom b.c.
  real(rkind),save,allocatable :: hdif(:,:) !horizontal diffusivity
  real(rkind),save,allocatable :: tem0(:,:) ! Initial temperature at nodes
  real(rkind),save,allocatable :: sal0(:,:) ! Initial salinity at nodes
  real(rkind),save,allocatable :: eta1(:)   ! Elevation at nodes at previous timestep
  real(rkind),save,allocatable :: eta2(:)   ! Elevation at nodes at current timestep
  real(rkind),save,allocatable :: we(:,:)   ! Vertical velocity at element centers & whole levels
  !Vertical velocity at element centers & whole levels, calculated using F.V.M. For hydrostatic 
  !model, this is the same as we(); for non-hydrostatic model, this is only used in upwind transport
  real(rkind),save,allocatable :: we_fv(:,:) 
  real(rkind),save,allocatable :: su2(:,:)  ! x-velocity at side centers & whole levels
  real(rkind),save,allocatable :: sv2(:,:)  ! y-velocity at side centers & whole levels
  real(rkind),save,allocatable :: ufg(:,:,:)  ! x-velocity at nodes in an element
  real(rkind),save,allocatable :: vfg(:,:,:)  ! y-velocity at nodes in an element
  real(rkind),save,allocatable :: tsd(:,:)  ! Temperature at side centers & whole levels
  real(rkind),save,allocatable :: ssd(:,:)  ! Salinity at side centers & whole levels
  real(rkind),save,allocatable :: tnd(:,:)  ! Temperature at nodes & whole levels
  real(rkind),save,allocatable :: snd(:,:)  ! Salinity at nodes & whole levels
  real(rkind),save,allocatable :: prho(:,:) ! Density at whole levels and either nodes or elements
!  real(rkind),save,allocatable :: sig_t(:,:) ! density anomaly
  real(rkind),save,allocatable :: q2(:,:)   ! Turbulent kinetic energy at sides & half levels
  real(rkind),save,allocatable :: xl(:,:)   ! Turbulent mixing length at sides & half levels
  real(rkind),save,allocatable :: xlmin2(:) ! ??? Turbulent mixing length
!  real(rkind),save,allocatable :: uu1(:,:)  ! X-Velocity at nodes & whole levels at previous timestep
!  real(rkind),save,allocatable :: vv1(:,:)  ! Y-Velocity at nodes & whole levels at previous timestep
!  real(rkind),save,allocatable :: ww1(:,:)  ! Z-Velocity at nodes & whole levels at previous timestep
  real(rkind),save,allocatable :: uu2(:,:)  ! X-Velocity at nodes & whole levels at current timestep
  real(rkind),save,allocatable :: vv2(:,:)  ! Y-Velocity at nodes & whole levels at current timestep
  real(rkind),save,allocatable :: ww2(:,:)  ! Z-Velocity at nodes & whole levels at current timestep
  real(rkind),save,allocatable :: bdef(:)   !bottom deformation
  real(rkind),save,allocatable :: bdef1(:)   !bottom deformation
  real(rkind),save,allocatable :: bdef2(:)   !bottom deformation
  real(rkind),save,allocatable :: dfh(:,:) !diffusivity
  integer,save,allocatable :: itier_nd(:,:) !multi-tier neighborhood; used in Kriging
  real(rkind),save,allocatable :: akrmat_nd(:,:,:)         ! Kriging matrix
  real(rkind),save,allocatable :: albedo(:)         ! albedo
  real(rkind),save,allocatable :: z_r(:)         ! z-cor. used in ts.ic
  real(rkind),save,allocatable :: tem1(:)         ! T profile in ts.ic
  real(rkind),save,allocatable :: sal1(:)         ! S profile in ts.ic
  real(rkind),save,allocatable :: cspline_ypp(:,:)         ! 2nd derivative in cubic spline for tem1,sal1
  real(rkind),save,allocatable :: rho_mean(:,:)         ! 2nd derivative in cubic spline for tem1,sal1
  real(rkind),save,allocatable :: Cdp(:)         ! drag at node

  ! Jet model. a - ambient, j - jet, 0 - original
  integer,save :: nofe,Ntj
  real(rkind),save :: dtj
  integer,save,allocatable :: ofe(:) !list of element near the outfall
  integer,save,allocatable :: ej(:,:),iej(:,:) !element and layer containing jet element k
  real(rkind),save,allocatable :: xMT(:,:),xMS(:,:),xMC(:,:),xMV(:,:)
!  real(rkind),save,allocatable :: xj(:),yj(:),zj(:),bj(:),hj(:),Tj(:),Sj(:),Cj(:),Mj(:),dMj(:),Mjg(:),Qjg(:) ! jet element parameters
  real(rkind),save,allocatable :: Volj(:),Taj(:),Saj(:),Caj(:),rhoaj(:) !,Maj(:),Qaj(:) ! ambient parameters
 
  ! Variables for global output files
  integer, parameter :: nbyte=4          !# bytes for output record size
  integer, parameter :: mnout=100        !max. # of output files
  integer, parameter :: mirec=1109000000 !max. record # to prevent output ~> 4GB
  integer,save :: iwrite
  character(len=48),save :: start_time,version,data_format='DataFormat v5.0'
  character(len=12),save :: ifile_char
  character(len=48),save,dimension(mnout) :: outfile,variable_nm,variable_dim
  integer,save :: ihot,ihfskip,nrec,nspool,igmp,noutgm,ifile,ifile_len,noutput,ifort12(100)
  integer,save,dimension(mnout) :: ichan,irec,iof !,mrec
!  real(rkind),save,dimension(mnout) :: vpos
  character(len=48),save :: a_48
  character(len=16),save :: a_16
  character(len= 8),save :: a_8
  character(len= 4),save :: a_4
        
  ! Error message string
  character(len=1000),save :: errmsg

  ! Constants used in UB closure
  character(len=2),save :: mid,stab
  real(rkind),save :: ubd0,ubd1,ubd2,ubd3,ubd4,ubd5, &
                      ubs0,ubs1,ubs2,ubs4,ubs5,ubs6, &
                      a2_cm03,schk,schpsi

  ! Miscellaneous variables
  integer,save :: nonhydro,iupwind_t,nz_r,ieqstate,imm,kr_co,indvel,ihconsv,isconsv,ihdif, &
                 &ntracers,mntr,ntracers2 !ntracers2=ntracers+2
  ! Max number of neighboring elements surrounding a node; max. number of srrounding nodes is mnei+1
  integer,save :: mnei 
  integer,save :: mnei_kr   ! Max # of Kriging nodes
  real(rkind),save :: h0,q2min,dt,tempmin,tempmax,saltmin,saltmax,vis_coe1,vis_coe2,h_bcc1, &
                     &velmin_btrack

  ! Some constants
  real(rkind),parameter :: small1=1.e-6 !small non-negative number
  real(rkind),parameter :: small2=small1*100 !slightly larger number
  real(rkind),parameter :: pi=3.1415926
  real(rkind),parameter :: rho0=1000. !1025. !ref. density for S=33 and T=10C

  ! For timing
  integer,parameter :: mxtimer=20          ! Max number of wallclock timers
  real(rkind),save :: wtimer(0:mxtimer,2)  ! Array of wallclock timers
                                           ! (:,1)=execution time
                                           ! (:,2)=communication time

  ! For debugging
  character(72) :: fdb  ! Name of debugging file
  integer :: lfdb       ! Length of debugging file name

  ! Tracers
  character(len=48) :: inputfile
  integer :: flag_model,flag_ic

contains


  subroutine release_gl(n,gl_ll)
  ! Free memory associated with global-to-local linked-list
    implicit none
    integer,intent(in) :: n
    type(llist_type),pointer :: gl_ll(:)
    integer i

    do i=1,n
      call release_llist(gl_ll(i)%next)
    enddo
    deallocate(gl_ll)
    nullify(gl_ll)

  end subroutine release_gl

  recursive subroutine release_llist(llentry)
  ! Free memory associated with linked-list
    implicit none
    type(llist_type),pointer :: llentry

    if(associated(llentry)) then
      call release_llist(llentry%next)
      deallocate(llentry)
      nullify(llentry)
    endif

  end subroutine release_llist


end module elfe_glbl
!===============================================================================
!===============================================================================
! END SELFE GLOBAL DATA MODULE
!===============================================================================
!===============================================================================
