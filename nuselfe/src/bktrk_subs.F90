!===============================================================================
!===============================================================================
! ELCIRC BACKTRACKING SUBROUTINES
!
! subroutine init_inter_btrack
! subroutine inter_btrack
! subroutine btrack
! subroutine quicksearch
! subroutine intersect2
! subroutine area_coord
!
!===============================================================================
!===============================================================================

subroutine init_inter_btrack
!-------------------------------------------------------------------------------
! Initialize data-types for inter-subdomain backtracking.
!-------------------------------------------------------------------------------
#ifdef USE_MPIMODULE
  use mpi
#endif
  use elfe_glbl
  use elfe_msgp
  implicit none
#ifndef USE_MPIMODULE
  include 'mpif.h'
#endif
  integer :: blockl(2),types(2),nmm
#if MPIVERSION==1
  integer :: displ(2),base
#elif MPIVERSION==2
  integer(kind=MPI_ADDRESS_KIND) :: displ(2),base
#endif
  type(bt_type) :: bttmp
!-------------------------------------------------------------------------------

  ! Dimension of inter-subdomain btrack arrays for sending and receiving
  call mpi_allreduce(nsa,nmm,1,itype,MPI_MAX,comm,ierr)
  mxnbt=s1_mxnbt*nmm*nvrt

  ! First part of bt_type is block of 9 integers
  ! (starting at bttmp%rank)
  blockl(1)=9
#if MPIVERSION==1
! displ(1) is the address
  call mpi_address(bttmp%rank,displ(1),ierr)
#elif MPIVERSION==2
  call mpi_get_address(bttmp%rank,displ(1),ierr)
#endif
  if(ierr/=MPI_SUCCESS) call parallel_abort('INIT_INTER_BTRACK: mpi_get_address',ierr)
  types(1)=itype

  ! Second part of bt_type is block of 11 doubles
  ! (starting at bttmp%dtbm)
  blockl(2)=11
#if MPIVERSION==1
  call mpi_address(bttmp%dtbm,displ(2),ierr)
#elif MPIVERSION==2
  call mpi_get_address(bttmp%dtbm,displ(2),ierr)
#endif
  if(ierr/=MPI_SUCCESS) call parallel_abort('INIT_INTER_BTRACK: mpi_get_address',ierr)
  types(2)=rtype

  ! Shift displ to compute actual displacements
  base=displ(1)
  displ(1)=displ(1)-base
  displ(2)=displ(2)-base

  ! MPI datatype for bt_type
#if MPIVERSION==1
  call mpi_type_struct(2,blockl,displ,types,bt_mpitype,ierr)
#elif MPIVERSION==2
  call mpi_type_create_struct(2,blockl,displ,types,bt_mpitype,ierr)
#endif
  if(ierr/=MPI_SUCCESS) call parallel_abort('INIT_INTER_BTRACK: type_create',ierr)
  call mpi_type_commit(bt_mpitype,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INIT_INTER_BTRACK: type_commit',ierr)

end subroutine init_inter_btrack

!===============================================================================
!===============================================================================

subroutine inter_btrack(itime,nbt,btlist)
!-------------------------------------------------------------------------------
! Routine for completing inter-subdomain backtracking.
!
! Input:
!   itime: global time stepping # (info only);
!   nbt: number of inter-subdomain trajectories
!   btlist: list of inter-subdomain trajectories
!
! Output:
!   btlist: list of completed inter-subdomain trajectories; not in original order
!-------------------------------------------------------------------------------
#ifdef USE_MPIMODULE
  use mpi
#endif
  use elfe_glbl
  use elfe_msgp
  implicit none
#ifndef USE_MPIMODULE
  include 'mpif.h'
#endif

  integer,intent(in) :: itime
  integer,intent(in) :: nbt
  type(bt_type),intent(inout) :: btlist(mxnbt)

  integer :: stat,i,ii,j,ie,irank,nnbrq,inbr,nbts,nbtd
  integer :: mxbtsend,mxbtrecv,mnbt
  real(rkind) :: xt,yt,zt,uuint,vvint,wwint,ttint,ssint
  logical :: lexit,bt_donel,bt_done
  integer :: icw
  real(rkind) :: cwtmp

  integer :: ncmplt,icmplt(nproc)
  integer,allocatable :: nbtsend(:),ibtsend(:,:),nbtrecv(:),ibtrecv(:,:)
#if MPIVERSION==1
  integer,allocatable :: bbtsend(:),bbtrecv(:)
#endif
  integer,allocatable :: btsend_type(:),btsend_rqst(:),btsend_stat(:,:)
  integer,allocatable :: btrecv_type(:),btrecv_rqst(:),btrecv_stat(:,:)
  type(bt_type),allocatable :: btsendq(:),btrecvq(:),bttmp(:),btdone(:)
#ifdef DEBUG
    integer,save :: ncalls=0
#endif
!-------------------------------------------------------------------------------

#ifdef DEBUG
  ncalls=ncalls+1
  fdb='interbtrack_0000'
  lfdb=len_trim(fdb)
  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
  if(ncalls==1) then
    open(30,file='outputs/'//fdb,status='replace')
  else
    open(30,file='outputs/'//fdb,status='old',position='append')
  endif
  write(30,'(a,3i6)') 'INTER_BTRACK START: ',itime,nbt
#endif

  ! Index of wall-timer
!  if(imode==1) then
    icw=4
!  else
!    icw=9
!  endif

  ! Compute max nbt for dimension parameter
#ifdef INCLUDE_TIMING
  cwtmp=mpi_wtime()
#endif
  call mpi_allreduce(nbt,mnbt,1,itype,MPI_MAX,comm,ierr) !mnbt>=1
#ifdef INCLUDE_TIMING
  wtimer(4,2)=wtimer(4,2)+mpi_wtime()-cwtmp
#endif
  mnbt=mnbt*s2_mxnbt !add a scale
 
  ! Allocate type bt_type arrays for sending and receiving
  allocate(btsendq(mnbt),btrecvq(mnbt*nnbr),bttmp(mnbt),btdone(mnbt),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: type bt_type allocation failure')

  ! Allocate communication data structures
  allocate(nbtsend(nnbr),ibtsend(mnbt,nnbr), &
  nbtrecv(nnbr),ibtrecv(mnbt,nnbr), &
  btsend_type(nnbr),btsend_rqst(nnbr),btsend_stat(MPI_STATUS_SIZE,nnbr), &
  btrecv_type(nnbr),btrecv_rqst(nnbr),btrecv_stat(MPI_STATUS_SIZE,nnbr),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: comm data allocation failure')
#if MPIVERSION==1
  allocate(bbtsend(mnbt),bbtrecv(mnbt),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bbtsend/recv allocation failure')
  bbtsend=1; bbtrecv=1; !blocksize is always 1
#endif

!  write(30,*)'Mark 1',nnbr,mxnbt,iegrpv(btlist(1:nbt)%iegb)

  ! Initialize send queue
  nbts=nbt
  btsendq(1:nbts)=btlist(1:nbts)

  ! Initialize completed bt count
  nbtd=0

  !-----------------------------------------------------------------------------
  ! Outer loop:
  ! > All ranks participate until all inter-subdomain backtracked trajectories
  !   are completed.
  ! > Completed trajectories are placed in btdone list.
  !-----------------------------------------------------------------------------
  outer_loop: do

#ifdef INCLUDE_TIMING
  ! Init communication timer
  cwtmp=mpi_wtime()
#endif

!  write(30,*)'Mark 2'

  ! Count and index sends
  nbtsend=0
  do i=1,nbts
    irank=iegrpv(btsendq(i)%iegb)
    inbr=ranknbr(irank)
    if(inbr==0) then
      write(errmsg,*) 'INTER_BTRACK: bt to non-neighbor!',irank
      call parallel_abort(errmsg)
    endif
    nbtsend(inbr)=nbtsend(inbr)+1
    if(nbtsend(inbr)>mnbt) call parallel_abort('bktrk_subs: overflow (1)')
    ibtsend(nbtsend(inbr),inbr)=i-1 !displacement
  enddo

  ! Set MPI bt send datatypes
  do inbr=1,nnbr
    if(nbtsend(inbr)/=0) then
#if MPIVERSION==1
      call mpi_type_indexed(nbtsend(inbr),bbtsend,ibtsend(1,inbr),bt_mpitype, &
      btsend_type(inbr),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nbtsend(inbr),1,ibtsend(1,inbr),bt_mpitype, &
      btsend_type(inbr),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: create btsend_type',ierr)
      call mpi_type_commit(btsend_type(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: commit btsend_type',ierr)
    endif
  enddo !inbr

  ! Post recvs for bt counts
  do inbr=1,nnbr
    call mpi_irecv(nbtrecv(inbr),1,itype,nbrrank(inbr),700,comm,btrecv_rqst(inbr),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: irecv 700',ierr)
  enddo

  ! Post sends for bt counts
  do inbr=1,nnbr
    call mpi_isend(nbtsend(inbr),1,itype,nbrrank(inbr),700,comm,btsend_rqst(inbr),ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: isend 700',ierr)
  enddo

!  write(30,*)'Mark 4'

  ! Wait for recvs to complete
  call mpi_waitall(nnbr,btrecv_rqst,btrecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitall recv 700',ierr)
  ! Wait for sends to complete
  call mpi_waitall(nnbr,btsend_rqst,btsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitall send 700',ierr)

!  write(30,*)'Mark 4.5'

  ! Set MPI bt recv datatypes
  i=0 !total # of recv from all neighbors 
  nnbrq=0; !# of "active" neighbors
  do inbr=1,nnbr
    if(nbtrecv(inbr)/=0) then
      nnbrq=nnbrq+1
      if(nbtrecv(inbr)>mnbt) call parallel_abort('bktrk_subs: overflow (3)')
      do j=1,nbtrecv(inbr); ibtrecv(j,inbr)=i+j-1; enddo; !displacement
      i=i+nbtrecv(inbr)
#if MPIVERSION==1
      call mpi_type_indexed(nbtrecv(inbr),bbtrecv,ibtrecv(1,inbr),bt_mpitype, &
      btrecv_type(inbr),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nbtrecv(inbr),1,ibtrecv(1,inbr),bt_mpitype, &
      btrecv_type(inbr),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: create btrecv_type',ierr)
      call mpi_type_commit(btrecv_type(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: commit btrecv_type',ierr)
!      write(30,*)'recev from',nbrrank(inbr),nbtrecv(inbr)
    endif
  enddo !inbr
 
  ! Check bound for btrecvq
  if(i>mnbt*nnbr) call parallel_abort('bktrk_subs: overflow (2)')

!  write(30,*)'Mark 5'

  ! Post sends for bt data
  do inbr=1,nnbr
    if(nbtsend(inbr)/=0) then
      ! btsendq(1)%rank is the starting address
      call mpi_isend(btsendq(1)%rank,1,btsend_type(inbr),nbrrank(inbr),701, &
      comm,btsend_rqst(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: isend 701',ierr)
    else
      btsend_rqst(inbr)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post recvs for bt data
  do inbr=1,nnbr
    if(nbtrecv(inbr)/=0) then
      call mpi_irecv(btrecvq(1)%rank,1,btrecv_type(inbr),nbrrank(inbr),701, &
      comm,btrecv_rqst(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: irecv 701',ierr)
    else
      btrecv_rqst(inbr)=MPI_REQUEST_NULL
    endif
  enddo

!  write(30,*)'Mark 6'

#ifdef INCLUDE_TIMING
  ! Add to communication timer
  wtimer(icw,2)=wtimer(icw,2)+mpi_wtime()-cwtmp
#endif

  !-----------------------------------------------------------------------------
  ! Inner loop: (for efficiency)
  ! > Process inter-subdomain backtracking receive queue until empty
  ! > Completed trajectories placed in btdone list
  ! > Exited trajectories placed in "temporary" send queue bttmp()
  !-----------------------------------------------------------------------------
  nbts=0 !current # of requests from myrank to all neighbors for inter-domain tracking
  inner_loop: do
  if(nnbrq==0) exit inner_loop

#ifdef DEBUG
  write(30,'(a)') 'INNER LOOP'
#endif

  ! Wait for some bt recvs to complete
  ! Parameters of mpi_waitsome:
  ! Inputs: nnbr - # of requests (dimension of btrecv_rqst); btrecv_rqst - array of requests;
  ! Outputs: ncmplt - # of completed requests; icmplt - array of indices of completed operations (integer);
  !          btrecv_stat - array of status objects for completed operations.

#ifdef INCLUDE_TIMING
  cwtmp=mpi_wtime()
#endif

  call mpi_waitsome(nnbr,btrecv_rqst,ncmplt,icmplt,btrecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitsome recv 701',ierr)
#ifdef INCLUDE_TIMING
  wtimer(icw,2)=wtimer(icw,2)+mpi_wtime()-cwtmp
#endif

!  write(30,*)'Mark 7'

  ! Perform local backtracking for received trajectories
  ! Process several neighbors at a time as soon as the info is received from them
  do ii=1,ncmplt
    inbr=icmplt(ii)
    do j=1,nbtrecv(inbr)
      i=ibtrecv(j,inbr)+1
      ie=iegl(btrecvq(i)%iegb)%id

!      write(30,*)'btrack #',ii,ie,btrecvq(i)%jvrt,btrecvq(i)%rank

      call btrack(btrecvq(i)%l0,btrecvq(i)%i0gb,btrecvq(i)%isbndy,btrecvq(i)%j0,btrecvq(i)%adv,btrecvq(i)%ndt, &
      btrecvq(i)%dtbm,btrecvq(i)%vis,btrecvq(i)%rt,btrecvq(i)%ut,btrecvq(i)%vt,btrecvq(i)%wt, &
      ie,btrecvq(i)%jvrt,btrecvq(i)%xt,btrecvq(i)%yt,btrecvq(i)%zt, &
      btrecvq(i)%tt,btrecvq(i)%st,lexit)

!      write(30,*)'btrack #',ii,ie,btrecvq(i)%jvrt,btrecvq(i)%rank

      if(lexit) then !backtracking exits augmented subdomain
        !Move point to "temporary" send queue
        nbts=nbts+1
        btrecvq(i)%iegb=ielg(ie)
        if(nbts>mnbt) call parallel_abort('bktrk_subs: overflow (5)')
        bttmp(nbts)=btrecvq(i)
      else !backtracking completed within augmented subdomain
        !Move point to done list
        nbtd=nbtd+1
        btrecvq(i)%iegb=ielg(ie)
        if(nbtd>mnbt) call parallel_abort('bktrk_subs: overflow (4)')
        btdone(nbtd)=btrecvq(i)
      endif
    enddo !j
  enddo !ii

  ! Decrement nnbrq according to number of recvs processed
  nnbrq=nnbrq-ncmplt

  !-----------------------------------------------------------------------------
  ! End inner loop
  !-----------------------------------------------------------------------------
  enddo inner_loop
 
#ifdef DEBUG
  write(30,'(a)') 'DONE INNER LOOP'
#endif

  ! All recv's are complete

#ifdef INCLUDE_TIMING
  ! Init communication timer
  cwtmp=mpi_wtime()
#endif

  ! Wait for sends to complete
  call mpi_waitall(nnbr,btsend_rqst,btsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitall send 701',ierr)

  ! Free MPI bt send datatypes
  do inbr=1,nnbr
    if(nbtsend(inbr)/=0) then
      call mpi_type_free(btsend_type(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: free btsend_type',ierr)
    endif
  enddo !inbr

  ! Free MPI bt recv datatypes
  do inbr=1,nnbr
    if(nbtrecv(inbr)/=0) then
      call mpi_type_free(btrecv_type(inbr),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: free btrecv_type',ierr)
    endif
  enddo !inbr

#ifdef DEBUG
  write(30,'(a,4i6)') 'CYCLE OUTER: ',itime,nbts,nbtd
#endif

  ! Exit outer loop if all backtracks are completed
  bt_donel=(nbts==0)
  call mpi_allreduce(bt_donel,bt_done,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
#ifdef INCLUDE_TIMING
  ! Add to communication timer
  wtimer(icw,2)=wtimer(icw,2)+mpi_wtime()-cwtmp
#endif
  if(bt_done) exit outer_loop

  ! Initialize send queue for next cycle of outer loop
  btsendq(1:nbts)=bttmp(1:nbts)


  !-----------------------------------------------------------------------------
  ! End outer loop
  !-----------------------------------------------------------------------------
  enddo outer_loop

  ! Deallocate some type bt_type arrays
  deallocate(btsendq,btrecvq,bttmp,stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bt type deallocation failure (1)')

  !-----------------------------------------------------------------------------
  ! All inter-subdomain backtracked trajectories completed. Communicate
  ! completed trajectories back to originating subdomain. Note that originating
  ! subdomain is not necessarily a neighbor to the subdomain where the
  ! trajectory was completed.
  ! Avoid communicating to itself as trajectory can come back!
  !-----------------------------------------------------------------------------

#ifdef INCLUDE_TIMING
  ! Init communication timer
  cwtmp=mpi_wtime()
#endif

#ifdef DEBUG
  write(30,'(a)') 'Start all-rank communication'
#endif

  ! Deallocate communication data structures
  deallocate(nbtsend,ibtsend,nbtrecv,ibtrecv, &
  btsend_type,btsend_rqst,btsend_stat, &
  btrecv_type,btrecv_rqst,btrecv_stat)
#if MPIVERSION==1
  deallocate(bbtsend,bbtrecv)
#endif

  ! nbtsend: # of sends from myrank to each proc excluding myrank
  allocate(nbtsend(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: nbtsend allocation failure')
  nbtsend=0
  do i=1,nbtd
    irank=btdone(i)%rank !destination rank
    if(irank/=myrank) nbtsend(irank)=nbtsend(irank)+1
  enddo !i
  mxbtsend=0; do irank=0,nproc-1; mxbtsend=max(mxbtsend,nbtsend(irank)); enddo;

  ! ibtsend: displacement or position in btdone(1:nbtd) list
  allocate(ibtsend(mxbtsend,0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: ibtsend allocation failure')
#if MPIVERSION==1
  allocate(bbtsend(mxbtsend),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bbtsend allocation failure')
  bbtsend=1 !blocksize is always 1
#endif
  nbtsend=0
  do i=1,nbtd
    irank=btdone(i)%rank
    if(irank/=myrank) then
      nbtsend(irank)=nbtsend(irank)+1
      ibtsend(nbtsend(irank),irank)=i-1 !displacement
    endif
  enddo !i

#ifdef DEBUG
  write(30,'(a,66i6)') 'INTER_BTRACK -- NBTSEND: ', &
  itime,(nbtsend(irank),irank=0,nproc-1)
#endif

  ! Allocate recv count array
  ! nbtrecv: # of recv's to myrank from each proc excluding myrank
  allocate(nbtrecv(0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: nbtrecv allocation failure')

  ! All to all scatter/gather of send counts
  call mpi_alltoall(nbtsend,1,itype,nbtrecv,1,itype,comm,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: alltoall nbtsend',ierr)

#ifdef DEBUG
  write(30,'(a,66i6)') 'INTER_BTRACK -- NBTRECV: ', &
  itime,(nbtrecv(irank),irank=0,nproc-1)
#endif

  ! Count max number of completed recv trajectories per rank and allocate ibtrecv
  mxbtrecv=0; do irank=0,nproc-1; mxbtrecv=max(mxbtrecv,nbtrecv(irank)); enddo;
  allocate(ibtrecv(mxbtrecv,0:nproc-1),stat=stat) !position in blist
  if(stat/=0) call parallel_abort('INTER_BTRACK: ibtrecv allocation failure')
#if MPIVERSION==1
  allocate(bbtrecv(mxbtrecv),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bbtrecv allocation failure')
  bbtrecv=1 !blocksize is always 1
#endif

  ! Allocate remaining communication data structures
  allocate(btsend_type(0:nproc-1),btsend_rqst(0:nproc-1), &
  btsend_stat(MPI_STATUS_SIZE,0:nproc-1), &
  btrecv_type(0:nproc-1),btrecv_rqst(0:nproc-1), &
  btrecv_stat(MPI_STATUS_SIZE,0:nproc-1),stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bt type/rqst/stat allocation failure')

  ! Set MPI bt send datatypes
  do irank=0,nproc-1
    if(nbtsend(irank)/=0) then
      if(irank==myrank) call parallel_abort('INTER_BTRACK: self communication (1)')
#if MPIVERSION==1
      call mpi_type_indexed(nbtsend(irank),bbtsend,ibtsend(1,irank),bt_mpitype, &
      btsend_type(irank),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nbtsend(irank),1,ibtsend(1,irank),bt_mpitype, &
      btsend_type(irank),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: create btsend_type',ierr)
      call mpi_type_commit(btsend_type(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: commit btsend_type',ierr)
    endif
  enddo !irank

  ! Set MPI bt recv datatypes
  i=0 !total # of recv's
  do irank=0,nproc-1
    if(nbtrecv(irank)/=0) then
      if(irank==myrank) call parallel_abort('INTER_BTRACK: self communication (2)')
      do j=1,nbtrecv(irank); ibtrecv(j,irank)=i+j-1; enddo; !displacement
      i=i+nbtrecv(irank)
#if MPIVERSION==1
      call mpi_type_indexed(nbtrecv(irank),bbtrecv,ibtrecv(1,irank),bt_mpitype, &
      btrecv_type(irank),ierr)
#elif MPIVERSION==2
      call mpi_type_create_indexed_block(nbtrecv(irank),1,ibtrecv(1,irank),bt_mpitype, &
      btrecv_type(irank),ierr)
#endif
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: create btrecv_type',ierr)
      call mpi_type_commit(btrecv_type(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: commit btrecv_type',ierr)
    endif
  enddo !irank

  ! Check bound for btlist
  if(i>mxnbt) call parallel_abort('INTER_BTRACK: overflow (6)')  

  ! Post sends for bt data
  do irank=0,nproc-1
    if(nbtsend(irank)/=0) then
      ! irank is checked to not be myrank
      call mpi_isend(btdone(1)%rank,1,btsend_type(irank),irank,711, &
      comm,btsend_rqst(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: isend 711',ierr)
    else
      btsend_rqst(irank)=MPI_REQUEST_NULL
    endif
  enddo

  ! Post recvs for bt data
  do irank=0,nproc-1
    if(nbtrecv(irank)/=0) then
      ! irank is checked to not be myrank
      call mpi_irecv(btlist(1)%rank,1,btrecv_type(irank),irank,711, &
      comm,btrecv_rqst(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: irecv 711',ierr)
    else
      btrecv_rqst(irank)=MPI_REQUEST_NULL
    endif
  enddo

  ! Wait for all bt recvs to complete
  call mpi_waitall(nproc,btrecv_rqst,btrecv_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitall recv 711',ierr)

  ! Wait for all bt sends to complete
  call mpi_waitall(nproc,btsend_rqst,btsend_stat,ierr)
  if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: waitall send 711',ierr)

  ! Append btlist with left-over from btdone list
  do ii=1,nbtd
    irank=btdone(ii)%rank
    if(irank==myrank) then
      i=i+1
      if(i>mxnbt) call parallel_abort('INTER_BTRACK: overflow (7)')  
      btlist(i)=btdone(ii)
#ifdef DEBUG
      write(30,*)'Back to myself!'
#endif
    endif
  enddo !ii

  ! Check if the total # fo recv's is nbt
  if(i/=nbt) call parallel_abort('bktrk_subs: mismatch (1)')

  ! Free MPI bt send datatypes
  do irank=0,nproc-1
    if(nbtsend(irank)/=0) then
      call mpi_type_free(btsend_type(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: free btsend_type',ierr)
    endif
  enddo !irank

  ! Free MPI bt recv datatypes
  do irank=0,nproc-1
    if(nbtrecv(irank)/=0) then
      call mpi_type_free(btrecv_type(irank),ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort('INTER_BTRACK: free btrecv_type',ierr)
    endif
  enddo !irank

  ! Deallocate btdone
  deallocate(btdone,stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: btdone deallocation failure')

  ! Deallocate communication data structures
  deallocate(nbtsend,ibtsend,nbtrecv,ibtrecv, &
  btsend_type,btsend_rqst,btsend_stat, &
  btrecv_type,btrecv_rqst,btrecv_stat,stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: rqst/stat deallocation failure')
#if MPIVERSION==1
  deallocate(bbtsend,bbtrecv,stat=stat)
  if(stat/=0) call parallel_abort('INTER_BTRACK: bbtsend/recv deallocation failure')
#endif

#ifdef INCLUDE_TIMING
  ! Add to communication timer
  wtimer(icw,2)=wtimer(icw,2)+mpi_wtime()-cwtmp
#endif

#ifdef DEBUG
  close(99)
#endif

end subroutine inter_btrack

!===============================================================================
!===============================================================================

      subroutine btrack(l_ns,ipsgb,ifl_bnd,j0,iadvf,ndelt,dtb_max,vis_coe,time_rm,uuint,vvint,wwint, &
                        nnel,jlev,xt,yt,zt,ttint,ssint,iexit)
!
!***************************************************************************
!									   
!     			Routine for backtracking. 			   
!       Input:
!             l_ns: node if l_ns=1; side if l_ns=2; element if l_ns=3;
!             ipsgb: global originating node or side or element index (info only);
!             ifl_bnd: Flag for originating node or side on the boundary (for Kriging);
!             j0: Originating vertical level (info only);
!             iadvf: advection flag (0 or 1: use Euler tracking; 2: R-K tracking) associated with originating position;
!             ndelt: # of subdivisions in Euler tracking (not needed for R-K);
!             dtb_max: max. tracking step for R-K (not needed for Euler);
!             vis_coe: weighting factor between continuous and discontinuous vel.
!             
!       Input/output:
!             time_rm: remaining time (<=dt);
!             uuint,vvint,wwint: starting and interpolated vel.;
!             xt,yt,zt: starting and current coordinates;
!             nnel,jlev: initial and final element and level;
!
!       Output:
!             ttint,ssint: interpolated T,S (only if iexit=.false.);
!             iexit: logical flag indicating backtracking exits augmented subdomain. If
!                    iexit=.true., nnel is inside the aug. domain and should also be inside
!                    one of the neighboring processes. (xt,yt) is inside nnel.
!									   
!***************************************************************************
!
      use elfe_glbl
      use elfe_msgp, only : parallel_abort,myrank
      implicit real(rkind)(a-h,o-z),integer(i-n)
      real(rkind), parameter :: per1=1.e-3 !used to check error in tracking
      real(rkind), parameter :: safety=0.8 !safyty factor in estimating time step

      integer, intent(in) :: l_ns,ipsgb,ifl_bnd,j0,iadvf,ndelt
      real(rkind), intent(in) :: dtb_max,vis_coe
      real(rkind), intent(inout) :: time_rm,uuint,vvint,wwint,xt,yt,zt
      integer, intent(inout) :: nnel,jlev
      real(rkind), intent(out) :: ttint,ssint
      logical, intent(out) :: iexit

      real(rkind) :: vxl(3,2),vyl(3,2),vzl(3,2),vxn(3),vyn(3),vzn(3)
      real(rkind) :: arco(3),t_xi(6),s_xi(6),sig(3),subrat(4),ztmp(nvrt), &
     &swild(10),swild2(nvrt,10),swild3(nvrt) 
      real(rkind) :: al_beta(mnei_kr+3,4),uvdata(mnei_kr,2) !,tsdata(mnei_kr,2)
      real(rkind) :: a(7),b(7,6),dc(6),rk(6,3)
      logical :: lrk

!     Constants used in 5th order R-K
      a(2)=0.2; a(3)=0.3; a(4)=0.6; a(5)=1; a(6)=0.875; a(7)=1
      b(2,1)=0.2; b(3,1)=0.075; b(3,2)=0.225; b(4,1)=0.3; b(4,2)=-0.9; b(4,3)=1.2
      b(5,1)=-11./54; b(5,2)=2.5; b(5,3)=-70./27; b(5,4)=35./27
      b(6,1)=1631./55296; b(6,2)=175./512; b(6,3)=575./13824; b(6,4)=44275./110592; b(6,5)=253./4096
!     b(7,*) are c(*)
      b(7,1)=37./378; b(7,2)=0; b(7,3)=250./621; b(7,4)=125./594; b(7,5)=0; b(7,6)=512./1771
!     dc() are c_i-c_i^*
      dc(1)=b(7,1)-2825./27648; dc(2)=0; dc(3)=b(7,3)-18575./48384
      dc(4)=b(7,4)-13525./55296; dc(5)=b(7,5)-277./14336; dc(6)=b(7,6)-0.25

      !Initial exit is false
      iexit=.false.

!      wwint00=wwint !save for quadratic interpolation

!...  Euler tracking
      if(iadvf==1.or.iadvf==0) then
!------------------------------------------------------------------------------------------------
      x0=xt
      y0=yt
      z0=zt
      idt=0 !iteration #
      do 
        idt=idt+1
        dtb=dt/ndelt
!        if(idt>1) dtb=dtb+trm
        dtb=min(dtb,time_rm)
        if(dtb<=0) call parallel_abort('BTRACK: dtb<=0')
        xt=x0-dtb*uuint
        yt=y0-dtb*vvint
        zt=z0-dtb*wwint

        call quicksearch(dtb,x0,y0,z0,nnel,jlev,xt,yt,zt,trm,iflqs1,kbpl,arco,zrat,ztmp,ss)

!       Check aug. exit
!       Note: for iflqs1=3 exit, continue to interpolate vel. and hopefully during next iteration iflqs1=2
!             will become true. In other words, for Euler tracking, the only exit mode is iflqs1=2, in 
!             order to get correct vel. before exiting.
        if(iflqs1>=2) then
!         Pt (xt,yt,zt) reset, which is inside nnel. jlev, uuint,vvint,wwint are unchanged
          if(iflqs1==3) time_rm=time_rm-(dtb-trm)
          iexit=.true.; return
        endif

!	nnel wet
        if(indvel==-1) then !interpolated hvel using P_1^NC
!         Interpolate in vertical 
          do j=1,3 !sides and nodes
            nd=nm(nnel,j)
            isd=js(nnel,j)
            if(interpol(nnel)==1) then
              kbb=kbs(isd)
              swild3(kbb:nvrt)=zs(kbb:nvrt,isd)
              swild2(kbb:nvrt,1)=su2(kbb:nvrt,isd)
              swild2(kbb:nvrt,2)=sv2(kbb:nvrt,isd)
              call vinter(nvrt,10,2,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
              vxn(j)=swild(1)
              vyn(j)=swild(2)

              kbb=kbp(nd)
              swild3(kbb:nvrt)=z(kbb:nvrt,nd)
              swild2(kbb:nvrt,1)=ww2(kbb:nvrt,nd)
              call vinter(nvrt,10,1,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
              vzn(j)=swild(1)
            else !pure S region
              vxn(j)=su2(jlev,isd)*(1-zrat)+su2(jlev-1,isd)*zrat
              vyn(j)=sv2(jlev,isd)*(1-zrat)+sv2(jlev-1,isd)*zrat !side
              vzn(j)=ww2(jlev,nd)*(1-zrat)+ww2(jlev-1,nd)*zrat !node
            endif
          enddo !j

!         Interpolate in horizontal
          uuint=vxn(1)*(1-2*arco(1))+vxn(2)*(1-2*arco(2))+vxn(3)*(1-2*arco(3))
          vvint=vyn(1)*(1-2*arco(1))+vyn(2)*(1-2*arco(2))+vyn(3)*(1-2*arco(3))
          wwint=vzn(1)*arco(1)+vzn(2)*arco(2)+vzn(3)*arco(3)

        else !indvel>=0; interpolated hvel using P_1
!	  No interpolate in time
          do j=1,3
            nd=nm(nnel,j)
            do l=1,2
              lev=jlev+l-2
              vxl(j,l)=(1-vis_coe)*uu2(lev,nd)+vis_coe*ufg(lev,nnel,j)
              vyl(j,l)=(1-vis_coe)*vv2(lev,nd)+vis_coe*vfg(lev,nnel,j)
              vzl(j,l)=ww2(lev,nd)
            enddo !l
          enddo !j

!	  Interpolate in vertical 
          do j=1,3
            if(interpol(nnel)==1) then
              nd=nm(nnel,j)
              kbb=kbp(nd)
              swild3(kbb:nvrt)=z(kbb:nvrt,nd)
              swild2(kbb:nvrt,1)=(1-vis_coe)*uu2(kbb:nvrt,nd)+vis_coe*ufg(kbb:nvrt,nnel,j)
              swild2(kbb:nvrt,2)=(1-vis_coe)*vv2(kbb:nvrt,nd)+vis_coe*vfg(kbb:nvrt,nnel,j)
              swild2(kbb:nvrt,3)=ww2(kbb:nvrt,nd)

              call vinter(nvrt,10,3,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
              vxn(j)=swild(1)
              vyn(j)=swild(2)
              vzn(j)=swild(3)

            else !pure S region
              vxn(j)=vxl(j,2)*(1-zrat)+vxl(j,1)*zrat
              vyn(j)=vyl(j,2)*(1-zrat)+vyl(j,1)*zrat
              vzn(j)=vzl(j,2)*(1-zrat)+vzl(j,1)*zrat
            endif
          enddo !j

!	  Interpolate in horizontal
          uuint=vxn(1)*arco(1)+vxn(2)*arco(2)+vxn(3)*arco(3)
          vvint=vyn(1)*arco(1)+vyn(2)*arco(2)+vyn(3)*arco(3)
          wwint=vzn(1)*arco(1)+vzn(2)*arco(2)+vzn(3)*arco(3)
        endif !indvel

        if(iflqs1==1) exit 
!       Update time_rm
        time_rm=time_rm-(dtb-trm)
        if(time_rm<=1.e-6*dt) exit

        x0=xt
        y0=yt
        z0=zt
      end do !idt
!------------------------------------------------------------------------------------------------
      endif !Euler

!...  5th-order R-K tracking
!     Results may vary with # of processors due to crossing of aug. domain
      if(iadvf==2) then
!------------------------------------------------------------------------------------------------
!      dtb0=min(dtb_max,dt)
      dtb0=min(dtb_max,time_rm)
!      time_rm=dt !time remaining
!     rk(i,1:3) is k_i in the book; 1:3 correspond to u,v,w
      rk(1,1)=-dtb0*uuint
      rk(1,2)=-dtb0*vvint
      rk(1,3)=-dtb0*wwint
      nnel0=nnel
      jlev0=jlev
      x0=xt
      y0=yt
      z0=zt
      uuint0=uuint
      vvint0=vvint
      wwint0=wwint

      icount=0 !adpative iteration # for death trap
      idt=0
!     If one sub-step hits death trap (iflqs1=1), then exit R-K loop successfully and do interpolation
      lrk=.false.
      loop5: do 
        idt=idt+1

!       k2-6 and k1 for the next step
        do k=2,7 !k=7 --> k1
          xt=x0
          yt=y0
          zt=z0
          do l=1,k-1
!           For k=7, (xt,yt,zt) is y_{n+1} in the book
            xt=xt+rk(l,1)*b(k,l)
            yt=yt+rk(l,2)*b(k,l)
            zt=zt+rk(l,3)*b(k,l)
          enddo !l

          nnel=nnel0
          jlev=jlev0
          call quicksearch(dtb0*a(k),x0,y0,z0,nnel,jlev,xt,yt,zt,trm,iflqs1,kbpl,arco,zrat,ztmp,ss)

          if(iflqs1==1) lrk=.true.

!	  nnel wet
!	  Interpolate in time; not used
!          if(k==7) then
!            trat=time_rm/dt !time_rm updated
!          else
!            trat=(time_rm-a(k)*dtb0)/dt !dtb0 not updated for k=2-6
!          endif
          if(indvel==-1) then !interpolated hvel using P_1^NC
!           Interpolate in vertical 
            do j=1,3 !sides and nodes
              nd=nm(nnel,j)
              isd=js(nnel,j)
              if(interpol(nnel)==1) then
                kbb=kbs(isd)
                swild3(kbb:nvrt)=zs(kbb:nvrt,isd)
                swild2(kbb:nvrt,1)=su2(kbb:nvrt,isd)
                swild2(kbb:nvrt,2)=sv2(kbb:nvrt,isd)
                call vinter(nvrt,10,2,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
                vxn(j)=swild(1)
                vyn(j)=swild(2)

                kbb=kbp(nd)
                swild3(kbb:nvrt)=z(kbb:nvrt,nd)
                swild2(kbb:nvrt,1)=ww2(kbb:nvrt,nd)
                call vinter(nvrt,10,1,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
                vzn(j)=swild(1)
              else !pure S region
                vxn(j)=su2(jlev,isd)*(1-zrat)+su2(jlev-1,isd)*zrat
                vyn(j)=sv2(jlev,isd)*(1-zrat)+sv2(jlev-1,isd)*zrat !side
                vzn(j)=ww2(jlev,nd)*(1-zrat)+ww2(jlev-1,nd)*zrat !node
              endif
            enddo !j

!           Interpolate in horizontal
            uuint=vxn(1)*(1-2*arco(1))+vxn(2)*(1-2*arco(2))+vxn(3)*(1-2*arco(3))
            vvint=vyn(1)*(1-2*arco(1))+vyn(2)*(1-2*arco(2))+vyn(3)*(1-2*arco(3))
            wwint=vzn(1)*arco(1)+vzn(2)*arco(2)+vzn(3)*arco(3)

          else !indvel>=0; interpolated hvel using P_1
            do j=1,3
              nd=nm(nnel,j)
              do l=1,2
                lev=jlev+l-2
                vxl(j,l)=(1-vis_coe)*uu2(lev,nd)+vis_coe*ufg(lev,nnel,j)
                vyl(j,l)=(1-vis_coe)*vv2(lev,nd)+vis_coe*vfg(lev,nnel,j)
                vzl(j,l)=ww2(lev,nd)
              enddo !l
            enddo !j

!	    Interpolate in vertical 
            do j=1,3
              if(interpol(nnel)==1) then
                nd=nm(nnel,j)
                kbb=kbp(nd)
                swild3(kbb:nvrt)=z(kbb:nvrt,nd)
                swild2(kbb:nvrt,1)=(1-vis_coe)*uu2(kbb:nvrt,nd)+vis_coe*ufg(kbb:nvrt,nnel,j)
                swild2(kbb:nvrt,2)=(1-vis_coe)*vv2(kbb:nvrt,nd)+vis_coe*vfg(kbb:nvrt,nnel,j)
                swild2(kbb:nvrt,3)=ww2(kbb:nvrt,nd)
                call vinter(nvrt,10,3,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
                vxn(j)=swild(1)
                vyn(j)=swild(2)
                vzn(j)=swild(3)
              else !pure S region
                vxn(j)=vxl(j,2)*(1-zrat)+vxl(j,1)*zrat
                vyn(j)=vyl(j,2)*(1-zrat)+vyl(j,1)*zrat
                vzn(j)=vzl(j,2)*(1-zrat)+vzl(j,1)*zrat
              endif
            enddo !j

!	    Interpolate in horizontal
            uuint=vxn(1)*arco(1)+vxn(2)*arco(2)+vxn(3)*arco(3)
            vvint=vyn(1)*arco(1)+vyn(2)*arco(2)+vyn(3)*arco(3)
            wwint=vzn(1)*arco(1)+vzn(2)*arco(2)+vzn(3)*arco(3)
          endif !indvel

!         Check aug. exit
          if(iflqs1==2) then
!           Pt (xt,yt,zt) reset, which is inside nnel. time_rm, jlev, uuint,vvint,wwint are updated
            iexit=.true.; return
          endif

!         Check aug. exit
          if(iflqs1==3) then
!           If the time step > 5s, redo; otherwise prepare inter-subdomain btrack
!           Each time the bnd is crossed R-K loses some accuracy
            if(dtb0>5) then
              dtb0=5
              rk(1,1)=-dtb0*uuint0
              rk(1,2)=-dtb0*vvint0
              rk(1,3)=-dtb0*wwint0
              cycle loop5
            else
!             (xt,yt,zt), nnel, jlev, uuint,vvint,wwint updated
              time_rm=time_rm-(dtb0*a(k)-trm)
              if(time_rm<=0) call parallel_abort('BTRACK: no time left')
              iexit=.true.; return
            endif
          endif

          if(k==7) then
            dx=0
            dy=0
            dz=0
            do l=1,k-1
              dx=dx+rk(l,1)*dc(l)
              dy=dy+rk(l,2)*dc(l)
              dz=dz+rk(l,3)*dc(l)
            enddo !l
            del_xy=sqrt(dx*dx+dy*dy)
            del_z=abs(dz)
            n1=nm(nnel,1) !wet node
            del0_xy=per1*radiel(nnel)
            jmin=max0(jlev,kbp(n1)+1)
            del0_z=per1*(z(jmin,n1)-z(jmin-1,n1))
            if(del0_z<=0) then
              write(errmsg,*)'BTRACK: Negative layer:',del0_z
              call parallel_abort(errmsg)
            endif
            if(del_xy==0) then
              dtb_xy=dtb0
            else 
              dtb_xy=safety*dtb0*(del0_xy/del_xy)**0.2
            endif
            if(del_z==0) then
              dtb_z=dtb0
            else 
              dtb_z=safety*dtb0*(del0_z/del_z)**0.2
            endif
!           Proposed time step for next iteration
            dtb=min(dtb_xy,dtb_z,time_rm,dtb_max) !dtb_xy,dtb_z,time_rm >0
            if(del_xy>del0_xy.or.del_z>del0_z) then
!             Debug
!             write(90,*)ielem,icount,time_rm,dtb0,dtb,idt

              if(icount<=5) then !max. 5 chances
!               Update dtb0 and k1, and redo current step
                dtb0=dtb
                rk(1,1)=-dtb0*uuint0
                rk(1,2)=-dtb0*vvint0
                rk(1,3)=-dtb0*wwint0
                icount=icount+1
                cycle loop5
              else !trap reached
                if(ifort12(2)==0) then
                  ifort12(2)=1
                  write(12,*)'Adaptivity trap reached after ',icount
                endif
              endif
            endif

!           Successful
            time_rm=time_rm-dtb0
            if(time_rm<0) then
              write(errmsg,*)'BTRACK: Negative time level:',time_rm
              call parallel_abort(errmsg)
            endif
!           dtb0 for next step
            dtb0=min(dtb,time_rm)
            if(dtb0==0.and.time_rm/=0) then
              write(errmsg,*)'BTRACK: Weird btrack:',dtb0,time_rm
              call parallel_abort(errmsg)
            endif
            icount=0 !reset

!           Debug
!            close(90)
!            open(90,file='fort.90')
!            rewind(90)
          endif !k==7

          if(k==7) then
            in1=1
          else
            in1=k
          endif
          rk(in1,1)=-dtb0*uuint !dtb0 updated for k=7
          rk(in1,2)=-dtb0*vvint
          rk(in1,3)=-dtb0*wwint
        enddo !k=2,7

        nnel0=nnel
        jlev0=jlev
        x0=xt
        y0=yt
        z0=zt
        uuint0=uuint
        vvint0=vvint
        wwint0=wwint

        if(lrk) exit loop5
        if(time_rm<=0) exit loop5
      end do loop5
!------------------------------------------------------------------------------------------------
      endif !R-K

!     Return if for element
!     Error: Kriging for wvel as well?
      if(l_ns==3) return

!     Kriging for vel. (excluding bnd nodes/sides)
      if(ifl_bnd/=1.and.krvel(nnel)==1) then
!       Do more inter-domain btrack if necessary to make sure the ending element is resident
        if(nnel>ne) then 
          !Nudge the final point a little; this may create variation using different # of processors
          time_rm=1.e-8*dt
          iexit=.true.; return
        else !nnel resident
!         Prepare data
          ie=ie_kr(nnel) !local index for Kriging
          if(ie==0) then
            write(errmsg,*)'Out of Kriging zone:',ielg(nnel)
            call parallel_abort(errmsg)   
          endif
          npp=itier_nd(ie,0)
          do i=1,npp
            nd=itier_nd(ie,i)
            if(idry(nd)==1) then !i.c.
              uvdata(i,1)=0
              uvdata(i,2)=0
            else !wet
              if(interpol(nnel)==1) then
                kbb=kbp(nd)
                swild3(kbb:nvrt)=z(kbb:nvrt,nd)
                swild2(kbb:nvrt,1)=uu2(kbb:nvrt,nd)
                swild2(kbb:nvrt,2)=vv2(kbb:nvrt,nd)
                call vinter(nvrt,10,2,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
                uvdata(i,1:2)=swild(1:2)
              else !along S
                uvdata(i,1)=uu2(jlev,nd)*(1-zrat)+uu2(jlev-1,nd)*zrat
                uvdata(i,2)=vv2(jlev,nd)*(1-zrat)+vv2(jlev-1,nd)*zrat
              endif
            endif
          enddo !all ball nodes

          do i=1,npp+3
            al_beta(i,1:2)=0
            do j=1,npp
              al_beta(i,1:2)=al_beta(i,1:2)+akrmat_nd(ie,i,j)*uvdata(j,1:2)
            enddo !j
          enddo !i

          uuint=al_beta(npp+1,1)+al_beta(npp+2,1)*xt+al_beta(npp+3,1)*yt
          vvint=al_beta(npp+1,2)+al_beta(npp+2,2)*xt+al_beta(npp+3,2)*yt
          do i=1,npp
            nd=itier_nd(ie,i)
            rr=sqrt((x(nd)-xt)**2+(y(nd)-yt)**2)
            covar2=covar(kr_co,rr)
            uuint=uuint+al_beta(i,1)*covar2
            vvint=vvint+al_beta(i,2)*covar2
          enddo !i
        endif !resident element
      endif !Kriging

!...  Interpolation at the foot for S,T
!...  Not calculated for upwind/TVD
      ttint=0; ssint=0 !initialize
!     nnel wet
      if(zrat<0.or.zrat>1) then
        write(errmsg,*)'BTRACK: zrat wrong:',jlev,zrat
        call parallel_abort(errmsg)
      endif
 
      if(iupwind_t/=0) return

!     Split-linear, quadratic or Kriging
      if(lqk(nnel)==1) then
!-----------------------------------------------------------------------
!     Split-linear
      ifl=0 !flag
      do i=1,4
        if(i<=3) then
          n1=nm(nnel,i)
          n2=js(nnel,nx(i,2))
          n3=js(nnel,nx(i,1))
          aa1=signa(xt,xcj(n2),xcj(n3),yt,ycj(n2),ycj(n3))
          aa2=signa(x(n1),xt,xcj(n3),y(n1),yt,ycj(n3))
          aa3=signa(x(n1),xcj(n2),xt,y(n1),ycj(n2),yt)
          subrat(i)=min(aa1,aa2,aa3)/area(nnel)
          !aa=abs(aa1)+abs(aa2)+abs(aa3)
          !subrat(i)=abs(aa-area(nnel)/4)*4/area(nnel)
          !if(subrat(i)<small2) then
          if(subrat(i)>=-small2) then
            ifl=1
            sig(1)=aa1*4/area(nnel)
            sig(2)=aa2*4/area(nnel)
            sig(1)=max(0._rkind,min(1._rkind,sig(1)))
            sig(2)=max(0._rkind,min(1._rkind,sig(2)))
            if(sig(1)+sig(2)>1) then
              sig(3)=0
              sig(2)=1-sig(1)
            else
              sig(3)=1-sig(1)-sig(2)
            endif

!           S,T extended
!            t_xi(1)=tnd(jlev,n1)*(1-zrat)+tnd(jlev-1,n1)*zrat
!            t_xi(2)=tsd(jlev,n2)*(1-zrat)+tsd(jlev-1,n2)*zrat

!            smax=-99; tmin=100; ibb=0 !flag
!           Prepare for cubic spline (not used)
!            do jj=1,3 !node or side
!              if(jj==1) then
!                vxl(jj,1)=z(kbp(n1),n1)
!                vxl(jj,2)=z(nvrt,n1)
!              else if(jj==2) then
!                vxl(jj,1)=zs(kbs(n2),n2)
!                vxl(jj,2)=zs(nvrt,n2)
!              else
!                vxl(jj,1)=zs(kbs(n3),n3)
!                vxl(jj,2)=zs(nvrt,n3)
!              endif  
!            enddo !jj
!            zmin=maxval(vxl(1:3,1)) !not really used
!            zmax=minval(vxl(1:3,2))

            do jj=1,3 !node or side
              if(jj==1) then
                kbb=kbp(n1)
                swild3(kbb:nvrt)=z(kbb:nvrt,n1)
                swild2(kbb:nvrt,1)=tnd(kbb:nvrt,n1)
                swild2(kbb:nvrt,2)=snd(kbb:nvrt,n1)
              else if(jj==2) then
                kbb=kbs(n2)
                swild3(kbb:nvrt)=zs(kbb:nvrt,n2)
                swild2(kbb:nvrt,1)=tsd(kbb:nvrt,n2)
                swild2(kbb:nvrt,2)=ssd(kbb:nvrt,n2)
              else !=3
                kbb=kbs(n3)
                swild3(kbb:nvrt)=zs(kbb:nvrt,n3)
                swild2(kbb:nvrt,1)=tsd(kbb:nvrt,n3)
                swild2(kbb:nvrt,2)=ssd(kbb:nvrt,n3)
              endif

              call vinter(nvrt,10,2,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
              t_xi(jj)=swild(1); s_xi(jj)=swild(2)

!             Below bottom extrapolation
!              call do_cubic_spline(nvrt-kbb+1,swild3(kbb:nvrt),swild2(kbb:nvrt,1), &
!     &0._rkind,0._rkind,1,zt,1,zmin,zmax,tmp)
!!             Impose bounds
!              tmax=maxval(swild2(kbb:nvrt,1))
!              tmin=minval(swild2(kbb:nvrt,1))
!              t_xi(jj)=min(tmax,max(tmin,tmp))
!              call do_cubic_spline(nvrt-kbb+1,swild3(kbb:nvrt),swild2(kbb:nvrt,2), &
!     &0._rkind,0._rkind,1,zt,1,zmin,zmax,tmp)
!              smax=maxval(swild2(kbb:nvrt,2))
!              smin=minval(swild2(kbb:nvrt,2))
!              s_xi(jj)=min(smax,max(smin,tmp))
            enddo !jj=1,3

            ttint=t_xi(1)*sig(1)+t_xi(2)*sig(2)+t_xi(3)*sig(3)
            ssint=s_xi(1)*sig(1)+s_xi(2)*sig(2)+s_xi(3)*sig(3)
            exit
          endif !subrat(i)<small2
        else !i=4
          n1=js(nnel,1)
          n2=js(nnel,2)
          n3=js(nnel,3)
          aa1=signa(xt,xcj(n2),xcj(n3),yt,ycj(n2),ycj(n3))
          aa2=signa(xcj(n1),xt,xcj(n3),ycj(n1),yt,ycj(n3))
          aa3=signa(xcj(n1),xcj(n2),xt,ycj(n1),ycj(n2),yt)
          subrat(i)=min(aa1,aa2,aa3)/area(nnel)
          !aa=abs(aa1)+abs(aa2)+abs(aa3)
          !subrat(i)=abs(aa-area(nnel)/4)*4/area(nnel)
          !if(subrat(i)<small2) then
          if(subrat(i)>=-small2) then
            ifl=1
            sig(1)=aa1*4/area(nnel)
            sig(2)=aa2*4/area(nnel)
            sig(1)=max(0._rkind,min(1._rkind,sig(1)))
            sig(2)=max(0._rkind,min(1._rkind,sig(2)))
            if(sig(1)+sig(2)>1) then
              sig(3)=0
              sig(2)=1-sig(1)
            else
              sig(3)=1-sig(1)-sig(2)
            endif

!            t_xi(1)=tsd(jlev,n1)*(1-zrat)+tsd(jlev-1,n1)*zrat
!            t_xi(2)=tsd(jlev,n2)*(1-zrat)+tsd(jlev-1,n2)*zrat

!            smax=-99; tmin=100; ibb=0 !flag
!            do jj=1,3 !side
!              isd=js(nnel,jj)
!              vxl(jj,1)=zs(kbs(isd),isd)
!              vxl(jj,2)=zs(nvrt,isd)
!            enddo !jj
!            zmin=maxval(vxl(1:3,1)) !not really used
!            zmax=minval(vxl(1:3,2))

            do jj=1,3 !side
              isd=js(nnel,jj)
              kbb=kbs(isd)
              swild3(kbb:nvrt)=zs(kbb:nvrt,isd)
              swild2(kbb:nvrt,1)=tsd(kbb:nvrt,isd)
              swild2(kbb:nvrt,2)=ssd(kbb:nvrt,isd)
              call vinter(nvrt,10,2,zt,kbb,nvrt,jlev,swild3,swild2,swild,ibelow)
              t_xi(jj)=swild(1); s_xi(jj)=swild(2)

!              call do_cubic_spline(nvrt-kbb+1,swild3(kbb:nvrt),swild2(kbb:nvrt,1), &
!     &0._rkind,0._rkind,1,zt,1,zmin,zmax,tmp)
!!             Impose bounds
!              tmax=maxval(swild2(kbb:nvrt,1))
!              tmin=minval(swild2(kbb:nvrt,1))
!              t_xi(jj)=min(tmax,max(tmin,tmp))
!              call do_cubic_spline(nvrt-kbb+1,swild3(kbb:nvrt),swild2(kbb:nvrt,2), &
!     &0._rkind,0._rkind,1,zt,1,zmin,zmax,tmp)
!              smax=maxval(swild2(kbb:nvrt,2))
!              smin=minval(swild2(kbb:nvrt,2))
!              s_xi(jj)=min(smax,max(smin,tmp))
            enddo !jj=1,3

!            if(ibb==1) then
!              if(smax<-98) then
!                write(11,*)'Max. S failed to exist (2):',zt,nnel
!                stop
!              endif
!              t_xi(1:3)=tmin; s_xi(1:3)=smax
!            endif

            ttint=t_xi(1)*sig(1)+t_xi(2)*sig(2)+t_xi(3)*sig(3)
            ssint=s_xi(1)*sig(1)+s_xi(2)*sig(2)+s_xi(3)*sig(3)
            exit
          endif !subrat(i)
        endif !i<=3
      enddo !i=1,4

      if(ifl==0) then
        write(errmsg,*)'BTRACK: Not in any sub-element',ielg(nnel),(subrat(i),i=1,4),xt,yt
        call parallel_abort(errmsg)
      endif

!-----------------------------------------------------------------------
      else if(lqk(nnel)==2) then
!-----------------------------------------------------------------------
!     Quadratic interplation
!     For pure S only
      if(ss<-1.or.ss>0.or.kbpl/=kz) then
        write(errmsg,*)'BTRACK: ss out of bound:',ss,kbpl
        call parallel_abort(errmsg)
      endif

      !if(wwint00>=0) then
      if(wwint>=0) then
        lin=-1 !lower interval
        if(jlev==kbpl+1) lin=-99
      else
        lin=1
        if(jlev==nvrt) lin=-98
      endif
      outq1=0; outq2=0 
      t_min=100; t_max=-100; s_min=100; s_max=-100

      do i=1,3 !nodes and sides
        nd=nm(nnel,i)
        isd=js(nnel,i)
        in1=nx(i,1)
        in2=nx(i,2)
!       check (range extended)
        if(tnd(jlev,nd)<-98.or.snd(jlev,nd)<-98.or.tsd(jlev,isd)<-98.or.ssd(jlev,isd)<-98) then
          write(errmsg,*)'BTRACK: Wrong S,T:',i,iplg(nd),islg(isd),tnd(jlev,nd),snd(jlev,nd),tsd(jlev,isd),ssd(jlev,isd)
          call parallel_abort(errmsg)
        endif

        if(abs(ss+1)<1.e-4.or.abs(ss)<1.e-4) then !two surfaces
          if(abs(ss+1)<1.e-4) then
            lev=kz
          else
            lev=nvrt
          endif
          t_n=tnd(lev,nd)
          s_n=snd(lev,nd)
          t_s=tsd(lev,isd)
          s_s=ssd(lev,isd)
          temp_min=min(tnd(lev,nd),tsd(lev,isd))
          temp_max=max(tnd(lev,nd),tsd(lev,isd))
          salt_min=min(snd(lev,nd),ssd(lev,isd))
          salt_max=max(snd(lev,nd),ssd(lev,isd))
        else if(lin<=-98) then !constrained bottom or surface
          if(lin==-99) then
!            zrat3=((zt-ztmp(kbpl))/(ztmp(kbpl+1)-ztmp(kbpl)))**2
            srat=((ss+1)/(sigma(2)+1))**2
          else 
!            zrat3=((zt-ztmp(nvrt))/(ztmp(nvrt-1)-ztmp(nvrt)))**2
!            zrat3=1-zrat3 !to put in same form
            srat=(ss/sigma(nsig-1))**2
            srat=1-srat !to put in same form
          endif
          if(srat<0.or.srat>1) then
            write(errmsg,*)'BTRACK: Out of bound (9):',srat
            call parallel_abort(errmsg)
          endif
          t_n=tnd(jlev,nd)*srat+tnd(jlev-1,nd)*(1-srat)
          s_n=snd(jlev,nd)*srat+snd(jlev-1,nd)*(1-srat)
          t_s=tsd(jlev,isd)*srat+tsd(jlev-1,isd)*(1-srat)
          s_s=ssd(jlev,isd)*srat+ssd(jlev-1,isd)*(1-srat)
          temp_min=min(tnd(jlev,nd),tnd(jlev-1,nd),tsd(jlev,isd),tsd(jlev-1,isd))
          temp_max=max(tnd(jlev,nd),tnd(jlev-1,nd),tsd(jlev,isd),tsd(jlev-1,isd))
          salt_min=min(snd(jlev,nd),snd(jlev-1,nd),ssd(jlev,isd),ssd(jlev-1,isd))
          salt_max=max(snd(jlev,nd),snd(jlev-1,nd),ssd(jlev,isd),ssd(jlev-1,isd))
        else !normal
          if(i==1) then !the following is indepdendent of loop i
            if(lin==1) then
              k1=jlev-1
            else !=-1
              k1=jlev-2
            endif           
            k2=k1+1; k3=k2+1
            if(k1<kbpl.or.k3>nvrt) then
              write(errmsg,*)'BTRACK: Weird level:',k1,k2,k3
              call parallel_abort(errmsg)
            endif
            k1s=k1-kz+1; k2s=k2-kz+1; k3s=k3-kz+1 !change to sigma indices
           
            denom=sigma(k3s)-2*sigma(k2s)+sigma(k1s)
            if(abs(denom)<1.e-5) then !degenerate
              xi=2*(ss-sigma(k2s))/(sigma(k3s)-sigma(k1s))
            else
              del=(sigma(k3s)-sigma(k1s))**2-8*(sigma(k2s)-ss)*denom
              if(del<0) then
                write(errmsg,*)'BTRACK: No inverse quadratic mapping:',del
                call parallel_abort(errmsg)
              endif
              icount=0
              vzn(1)=(sigma(k1s)-sigma(k3s)+sqrt(del))/2/denom !!temporary storage
              vzn(2)=(sigma(k1s)-sigma(k3s)-sqrt(del))/2/denom
              xi_m=vzn(1) !for no root scenario
              if(abs(vzn(2))<abs(vzn(1))) xi_m=vzn(2)
              do l=1,2
                if(abs(vzn(l))<=1+1.e-4) then
                  icount=icount+1
                  vxn(icount)=max(-1._rkind,min(1._rkind,vzn(l)))
                endif
              enddo !l
              if(icount==0) then
                write(errmsg,*)'BTRACK: No roots in inverse quadratic mapping:',(vzn(l),l=1,2), &
                                ss,sigma(k1s),sigma(k2s),sigma(k3s)
                call parallel_abort(errmsg)
              else if(icount==2) then
                if(ifort12(14)==0) then
                  ifort12(14)=1
                  write(12,*)'Warning: 2 roots:',(vxn(l),l=1,2)
                  write(12,*)ss,sigma(k1s),sigma(k2s),sigma(k3s),k1,k2,k3,kbpl
                endif
                xi=vxn(1)
              else !=1
                xi=vxn(1)
              endif
            endif !degenerate
           
            phi1=xi*(xi-1)/2; phi2=1-xi*xi; phi3=xi*(xi+1)/2
          endif !i==1

          t_n=tnd(k1,nd)*phi1+tnd(k2,nd)*phi2+tnd(k3,nd)*phi3
          s_n=snd(k1,nd)*phi1+snd(k2,nd)*phi2+snd(k3,nd)*phi3
          t_s=tsd(k1,isd)*phi1+tsd(k2,isd)*phi2+tsd(k3,isd)*phi3
          s_s=ssd(k1,isd)*phi1+ssd(k2,isd)*phi2+ssd(k3,isd)*phi3
          temp_min=min(tnd(k1,nd),tnd(k2,nd),tnd(k3,nd),tsd(k1,isd),tsd(k2,isd),tsd(k3,isd))
          temp_max=max(tnd(k1,nd),tnd(k2,nd),tnd(k3,nd),tsd(k1,isd),tsd(k2,isd),tsd(k3,isd))
          salt_min=min(snd(k1,nd),snd(k2,nd),snd(k3,nd),ssd(k1,isd),ssd(k2,isd),ssd(k3,isd))
          salt_max=max(snd(k1,nd),snd(k2,nd),snd(k3,nd),ssd(k1,isd),ssd(k2,isd),ssd(k3,isd))
        endif !normal

        outq1=outq1+t_n*(2*arco(i)*arco(i)-arco(i))+t_s*4*arco(in1)*arco(in2)
        outq2=outq2+s_n*(2*arco(i)*arco(i)-arco(i))+s_s*4*arco(in1)*arco(in2)
        if(temp_min<t_min) t_min=temp_min
        if(temp_max>t_max) t_max=temp_max
        if(salt_min<s_min) s_min=salt_min
        if(salt_max>s_max) s_max=salt_max
      enddo !i=1,3

      if(t_min>t_max) then
        write(errmsg,*)'BTRACK: Illegal min/max for temp:',t_min,t_max,ielg(nnel)
        call parallel_abort(errmsg)
      endif
      if(s_min>s_max) then
        write(errmsg,*)'BTRACK: Illegal min/max for salt:',s_min,s_max,ielg(nnel)
        call parallel_abort(errmsg)
      endif

      ttint=max(t_min,min(t_max,outq1))
      ssint=max(s_min,min(s_max,outq2))

!-----------------------------------------------------------------------
      endif !linear or quadratic

      end subroutine btrack

!===============================================================================
!===============================================================================

      subroutine quicksearch(time,x0,y0,z0,nnel,jlev,xt,yt,zt,trm,nfl,kbpl,arco,zrat,ztmp,ss)
!
!********************************************************************************
!										*
!     Straightline search algorithm. 
!
!     Inputs: 
!       time: time step from (x0,y0,z0) to (xt,yt,zt);
!       x0,y0,z0:  starting pt. (x0,y0) must be inside nnel;
!       nnel,jlev: starting element and level. nnel must be inside aug. domain.
!       xt,yt,zt: projected end pt;
!       In addition, su2,sv2,ww2 are also used.
! 
!     Outputs:
!      nnel, jlev: end element and level. nnel must be inside aug. domain;
!      (xt,yt,zt):  the updated end pt (if so); 
!      trm: time remaining (trm<=time). trm=0 unless the path crosses aug. bnd (nfl=2 or3);
!      nfl: a flag. nfl=1 if a bnd or dry element is hit and vel. there is small,              
!           or death trap is reached. In this case, all outputs are valid, 
!           and the time stepping in btrack is exited successfully. 
!           If nfl=2 (hit aug. bnd upon entry) or 3 (hit aug. bnd during iteration),
!           nnel, jlev, (xt,yt,zt) and trm are updated, and nnel should be inside 
!           one of the neighbor process. If nfl=2, (xt,yt,zt) are moved to (x0,y0,z0).
!      kbpl: the local bottom level at (xt,yt);
!      arco(3): area coordinates of (xt,yt);
!      zrat: vertical ratio of zt;
!      ztmp(nvrt):  z-coords. at (xt,yt);
!      ss: the sigma-coord of zt if nnel is inside pure S region. Otherwise -99.
!********************************************************************************
!
      use elfe_glbl
      use elfe_msgp, only : myrank,parallel_abort
      implicit real(rkind)(a-h,o-z),integer(i-n)

      real(rkind), intent(in) :: time,x0,y0,z0
      integer, intent(inout) :: nnel,jlev
      real(rkind), intent(inout) :: xt,yt,zt
      integer, intent(out) :: nfl,kbpl
      real(rkind), intent(out) :: trm,arco(3),zrat,ztmp(nvrt),ss
      real(rkind) :: wild(10,2),wild2(10,2)

!     Debug
!      fdb='qs_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(98,file='outputs/'//fdb,status='unknown')

      if(idry_e(nnel)==1) then
        write(errmsg,*)'QUICKSEARCH: Starting element is dry'
        call parallel_abort(errmsg)
      endif

      nfl=0
      trm=time !time remaining

!     Starting element nel
      nel=nnel
      !aa=0
      !aa1=0
      do i=1,3
        n1=nm(nel,i)
        n2=nm(nel,nx(i,1))
        wild(i,1)=signa(x(n1),x(n2),x0,y(n1),y(n2),y0)
        wild(i,2)=signa(x(n1),x(n2),xt,y(n1),y(n2),yt)
        !aa=aa+abs(signa(x(n1),x(n2),x0,y(n1),y(n2),y0))
        !aa1=aa1+abs(signa(x(n1),x(n2),xt,y(n1),y(n2),yt))
      enddo !i
      ar_min1=minval(wild(1:3,1))/area(nel)
      ar_min2=minval(wild(1:3,2))/area(nel)
!      ae=abs(aa-area(nel))/area(nel)
!      if(ae>small1) then
      if(ar_min1<-small1) then
        write(errmsg,*)'QUICKSEARCH: (x0,y0) not in nnel initially',ae,ielg(nnel),x0,y0,ar_min1
        call parallel_abort(errmsg)
      endif

!      ae=abs(aa1-area(nel))/area(nel)
!      if(ae<small1) then
      if(ar_min2>-small1) then
        nnel=nel
        trm=0
        go to 400
      endif

!     Search surrounding elements as well for (xt,yt)
!      do l=1,3 !nodes
!        nd=nm(nel,l)
!        do j=1,nne(nd)
!          ie=ine(nd,j)
!          if(ie>0.and.idry_e(ie)==0.and.ie/=nel) then
!            aa1=0
!            do i=1,3
!              n1=nm(ie,i)
!              n2=nm(ie,nx(i,1))
!              aa1=aa1+abs(signa(x(n1),x(n2),xt,y(n1),y(n2),yt))
!            enddo !i
!            ae=abs(aa1-area(ie))/area(ie)
!            if(ae<small1) then
!              nnel=ie
!              trm=0
!              go to 400
!            endif
!          endif !ie>0.and.ie/=nel
!        enddo !j
!      enddo !l

!     (xt,yt) not in nel, and thus (x0,y0) and (xt,yt) are distinctive
!      xcg=(1-small2)*x0+small2*xctr(nel)
!      ycg=(1-small2)*y0+small2*yctr(nel)
!     Initialize (moving) starting pt
      xcg=x0; ycg=y0

      pathl=sqrt((xt-xcg)**2+(yt-ycg)**2)
      if(pathl==0) then
        write(errmsg,*)'QUICKSEARCH: Zero path',x0,y0,xt,yt,xcg,ycg
        call parallel_abort(errmsg)
      endif

!     Starting edge nel_j
      wild=0; wild2=0 !initialize for debugging output
      nel_j=0
      do j=1,3
        jd1=nm(nel,nx(j,1))
        jd2=nm(nel,nx(j,2))
        ar1=signa(xcg,x(jd1),xt,ycg,y(jd1),yt)    
        ar2=signa(xcg,xt,x(jd2),ycg,yt,y(jd2))    
        wild2(j,1)=ar1; wild2(j,2)=ar2
        if(ar1>=0.and.ar2>=0) then
          call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
          wild(j,1)=tt1; wild(j,2)=tt2; wild(3+j,1)=xin; wild(3+j,2)=yin
          if(iflag/=1) then
            write(errmsg,*)'QUICKSEARCH: Found no intersecting edges (1):', &
     &ielg(nel),xcg,ycg,xt,yt,ar_min2,wild(1:3,1:2),wild(4:6,1:2),ar1,ar2
            call parallel_abort(errmsg)
          else
            nel_j=j; exit
          endif
        endif !ar1>=0.and.ar2>=0
      enddo !j=1,3
      if(nel_j==0) then
        write(errmsg,*)'QUICKSEARCH: no intersecting edge: ',ielg(nel),xcg,ycg,xt,yt,wild2(1:3,1:2)
        call parallel_abort(errmsg)
      endif

!     Check aug. exit
!     nnel, jlev, and trm are unchanged; (xt,yt,zt) moved to (x0,y0,z0)
!     to be ready for inter-subdomain tracking
      if(ic3(nel,nel_j)<0) then
        xt=x0; yt=y0; zt=z0; nnel=nel
        nfl=2; go to 400
      endif

      zin=z0 !intialize
      it=0
      loop4: do
!----------------------------------------------------------------------------------------
      it=it+1

!     Exit loop if death trap is reached
      if(it>1000) then
        if(ifort12(3)==0) then
          ifort12(3)=1
          write(12,*)'QUICKSEARCH: Death trap reached'
        endif
        nfl=1
        xt=xin
        yt=yin
        zt=zin
        nnel=nel
        trm=0 !min(trm,time)
        exit loop4
      endif
      md1=nm(nel,nx(nel_j,1))
      md2=nm(nel,nx(nel_j,2))
      
!     Compute z position 
      dist=sqrt((xin-xt)**2+(yin-yt)**2)

!     Debug
!      write(98,*)dist/pathl,ielg(nel)

      !if(dist/pathl>1+1.0e-4) then
!      if(dist/pathl>1+small2) then
!        write(errmsg,*)'QUICKSEARCH: Path overshot',dist/pathl,xin,yin,xt,yt,ielg(nel)
!        call parallel_abort(errmsg)
!      endif
      tmp=min(1._rkind,dist/pathl)
      zin=zt-tmp*(zt-zin)
      trm=trm*tmp !time remaining
      
      pathl=sqrt((xin-xt)**2+(yin-yt)**2)
      if(pathl==0.or.trm==0) then
        write(errmsg,*)'QUICKSEARCH: Target reached'
        call parallel_abort(errmsg)
      endif

!     Check for aug. exit
      if(ic3(nel,nel_j)<0) then
!       nnel is the last element inside aug. domain
!       xt,yt,zt, jlev, and trm are updated.
        nfl=3
        xt=xin
        yt=yin
        zt=zin
        nnel=nel
        trm=min(trm,time) !>0
        nnel=nel
        exit loop4
      endif

!     Next element is inside aug. domain
      lit=0 !flag
!     For horizontal exit and dry elements, compute tangential vel.,
!     update target (xt,yt,zt) and continue.
      if(ic3(nel,nel_j)==0.or.idry_e(ic3(nel,nel_j))==1) then
        lit=1
        isd=js(nel,nel_j)
        if(isidenode(isd,1)+isidenode(isd,2)/=md1+md2) then
          write(errmsg,*)'QUICKSEARCH: Wrong side'
          call parallel_abort(errmsg)
        endif

!       Nudge intersect (xin,yin), and update starting pt
        eps=100*small2
        xin=(1-eps)*xin+eps*xctr(nel)
        yin=(1-eps)*yin+eps*yctr(nel)
        xcg=xin
        ycg=yin

        vtan=-su2(jlev,isd)*sny(isd)+sv2(jlev,isd)*snx(isd)
        xvel=-vtan*sny(isd)
        yvel=vtan*snx(isd)
        zvel=(ww2(jlev,md1)+ww2(jlev,md2))/2
        xt=xin-xvel*trm
        yt=yin-yvel*trm
        zt=zin-zvel*trm
        hvel=sqrt(xvel**2+yvel**2)
        !if(hvel<1.e-4) then
        if(hvel<=velmin_btrack) then
          nfl=1
          xt=xin
          yt=yin
          zt=zin
          nnel=nel
          trm=0
          exit loop4
        endif
        pathl=hvel*trm
      endif !abnormal cases

!     Search for nel's neighbor with edge nel_j, or in abnormal cases, the same element
      if(lit==0) nel=ic3(nel,nel_j) !next front element

      do i=1,3
        k1=nm(nel,i)
        k2=nm(nel,nx(i,1))
        wild(i,1)=signa(x(k1),x(k2),xt,y(k1),y(k2),yt)
      enddo !i
      ar_min1=minval(wild(1:3,1))/area(nel)
      if(ar_min1>-small1) then
        nnel=nel
        trm=0
        exit loop4
      endif

!     Next intersecting edge
      wild=0; wild2=0 !initialize for output
      nel_j=0
      do j=1,3
        jd1=nm(nel,nx(j,1))
        jd2=nm(nel,nx(j,2))
!       For abnormal case, same side (border side) cannot be hit again
        if(jd1==md1.and.jd2==md2.or.jd2==md1.and.jd1==md2) cycle
        ar1=signa(xcg,x(jd1),xt,ycg,y(jd1),yt)
        ar2=signa(xcg,xt,x(jd2),ycg,yt,y(jd2))
        wild2(j,1)=ar1; wild2(j,2)=ar2
        if(ar1>=0.and.ar2>=0) then
          call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
          wild(j,1)=tt1; wild(j,2)=tt2; wild(3+j,1)=xin; wild(3+j,2)=yin
          if(iflag/=1) then
            write(errmsg,*)'QUICKSEARCH: Failed to find next edge (2):',lit,xcg,ycg,xt,yt,ielg(nel), &
     &iplg(md1),iplg(md2),ar_min1,wild(1:3,1:2),wild(4:6,1:2),ar1,ar2
            call parallel_abort(errmsg)
          else
            nel_j=j; !next front edge
            cycle loop4
          endif
        endif !ar1>=0.and.ar2>=0
      enddo !j
      if(nel_j==0) then
        write(errmsg,*)'QUICKSEARCH: no intersecting edge (2): ',ielg(nel),xcg,ycg,xt,yt,wild2(1:3,1:2)
        call parallel_abort(errmsg)
      endif

!----------------------------------------------------------------------------------------
      end do loop4 

400   continue
!     No vertical exit from domain
      if(idry_e(nnel)==1) then
        write(errmsg,*)'QUICKSEARCH: Ending element is dry:',ielg(nnel)
        call parallel_abort(errmsg)
      endif

!     Compute area & sigma coord.
      call area_coord(nnel,xt,yt,arco)
      n1=nm(nnel,1)
      n2=nm(nnel,2)
      n3=nm(nnel,3)
      etal=eta2(n1)*arco(1)+eta2(n2)*arco(2)+eta2(n3)*arco(3)
      dep=dp(n1)*arco(1)+dp(n2)*arco(2)+dp(n3)*arco(3)
      if(etal+dep<h0) then
        write(errmsg,*)'QUICKSEARCH: Weird wet element in quicksearch:',ielg(nnel),eta2(n1),eta2(n2),eta2(n3)
        call parallel_abort(errmsg)
      endif

!     Compute z-coordinates
      do k=kz,nvrt
        kin=k-kz+1
        hmod2=min(dep,h_s)
        if(hmod2<=h_c) then
          ztmp(k)=sigma(kin)*(hmod2+etal)+etal
        else if(etal<=-h_c-(dep-h_c)*theta_f/s_con1) then
          write(errmsg,*)'QUICKSEARCH: Pls choose a larger h_c (2):',etal,h_c
          call parallel_abort(errmsg)
        else
          ztmp(k)=etal*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
        endif

!       Following to prevent underflow
        if(k==kz) ztmp(k)=-hmod2
        if(k==nvrt) ztmp(k)=etal
      enddo !k

      if(dep<=h_s) then
        kbpl=kz
      else !z levels
!       Find bottom index
        kbpl=0
        do k=1,kz-1
          if(-dep>=ztot(k).and.-dep<ztot(k+1)) then
            kbpl=k
            exit
          endif
        enddo !k
        if(kbpl==0) then
          write(errmsg,*)'QUICKSEARCH: Cannot find a bottom level at foot:',dep
          call parallel_abort(errmsg)
        endif
        ztmp(kbpl)=-dep
        do k=kbpl+1,kz-1
          ztmp(k)=ztot(k)
        enddo !k
      endif

      do k=kbpl+1,nvrt
        if(ztmp(k)-ztmp(k-1)<=0) then
          write(errmsg,*)'QUICKSEARCH: Inverted z-level in quicksearch:',ielg(nnel),etal,dep,ztmp(k)-ztmp(k-1)
          call parallel_abort(errmsg)
        endif
      enddo !k

      if(zt<=ztmp(kbpl)) then
        zt=ztmp(kbpl)
        zrat=1
        jlev=kbpl+1
      else if(zt>=ztmp(nvrt)) then
        zt=ztmp(nvrt)
        zrat=0
        jlev=nvrt
      else
        jlev=0
        do k=kbpl,nvrt-1
          if(zt>=ztmp(k).and.zt<=ztmp(k+1)) then 
            jlev=k+1
            exit
          endif
        enddo !k
        if(jlev==0) then
          write(errmsg,*)'QUICKSEARCH: Cannot find a vert. level:',zt,etal,dep,(ztmp(k),k=kbpl,nvrt)
          call parallel_abort(errmsg)
        endif
        zrat=(ztmp(jlev)-zt)/(ztmp(jlev)-ztmp(jlev-1))
      endif

      if(zrat<0.or.zrat>1) then
        write(errmsg,*)'QUICKSEARCH: Sigma coord. wrong (4):',jlev,zrat
        call parallel_abort(errmsg)
      endif

      if(kbpl==kz) then !in pure S region
        ss=(1-zrat)*sigma(jlev-kz+1)+zrat*sigma(jlev-kz)
      else
        ss=-99
      endif


!      if(ss<sigma(jlev-1).or.ss>sigma(jlev)) then
!        write(11,*)'Sigma coord. wrong (5):',jlev,ss,sigma(jlev-1),sigma(jlev)
!        stop
!      endif

      end subroutine quicksearch

!===============================================================================
!===============================================================================

      subroutine intersect2(x1,x2,x3,x4,y1,y2,y3,y4,iflag,xin,yin,tt1,tt2)
!
!********************************************************************************
!										*
!     Program to detect if two segments (1,2) and (3,4) have common pts   	*
!     Assumption: the 4 pts are distinctive.					*
!     The eqs. for the 2 lines are: X=X1+(X2-X1)*tt1 and X=X3+(X4-X3)*tt2.	*
!     Output: iflag: 0: no intersection or colinear; 1: exactly 1 intersection.	*
!     If iflag=1, (xin,yin) is the intersection.				*
!     Modified to only check tt2, assuming 0<=tt1<=1 is already assured.
!										*
!********************************************************************************
!
      use elfe_glbl, only: rkind !,small2
      implicit real(rkind)(a-h,o-z), integer(i-n)
      real(rkind), parameter :: small=0.0 !small positive number or 0

      real(rkind), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4
      integer, intent(out) :: iflag
      real(rkind), intent(out) :: xin,yin,tt1,tt2

      tt1=-1000
      tt2=-1000
      xin=-1.e25; yin=xin
      iflag=0
      delta=(x2-x1)*(y3-y4)-(y2-y1)*(x3-x4)
      delta1=(x3-x1)*(y3-y4)-(y3-y1)*(x3-x4)
      delta2=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

      if(delta/=0) then
        tt1=delta1/delta
        tt2=delta2/delta
        !if(tt1>=-small.and.tt1<=1+small.and.tt2>=-small.and.tt2<=1+small) then
        if(tt2>=-small.and.tt2<=1+small) then
          iflag=1
          xin=x3+(x4-x3)*tt2
          yin=y3+(y4-y3)*tt2
        endif
      endif

      end subroutine intersect2

!===============================================================================
!===============================================================================

      subroutine area_coord(nnel,xt,yt,arco)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                       !
!       Compute area coordinates of pt (xt,yt), which must be inside element nnel.      !
!       Impose bounds for area coordinates.                                             !
!                                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use elfe_glbl
      implicit real(rkind)(a-h,o-z),integer(i-n)
      integer, intent(in) :: nnel
      real(rkind), intent(in) :: xt,yt
      real(rkind), intent(out) :: arco(3)

      n1=nm(nnel,1)
      n2=nm(nnel,2)
      n3=nm(nnel,3)
      arco(1)=signa(xt,x(n2),x(n3),yt,y(n2),y(n3))/area(nnel)
      arco(2)=signa(x(n1),xt,x(n3),y(n1),yt,y(n3))/area(nnel)
      arco(1)=max(0._rkind,min(1._rkind,arco(1)))
      arco(2)=max(0._rkind,min(1._rkind,arco(2)))
      if(arco(1)+arco(2)>1) then
        arco(3)=0
        arco(1)=min(1._rkind,max(0._rkind,arco(1)))
        arco(2)=1-arco(1)
      else
        arco(3)=1-arco(1)-arco(2)
      endif

      end subroutine area_coord


!===============================================================================
!===============================================================================
! END ELCIRC BACKTRACKING SUBROUTINES
!===============================================================================
!===============================================================================
