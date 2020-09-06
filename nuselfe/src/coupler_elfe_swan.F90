!=================================================================================
!  Coupling SELFE and SWAN
!  To test alone:
!    mpif90  -DMPIVERSION=2  -O2 -o coupler_elfe_swan coupler_elfe_swan.F90
!=================================================================================
      program coupler_elfe_swan
#ifdef USE_MPIMODULE
      use mpi
#endif      
      implicit real(8)(a-h,o-z),integer(i-n)
#ifndef USE_MPIMODULE
      include 'mpif.h'
#endif
      integer :: comm1,comm2,comm3,nproc0,myrank0,ier,myproc_all
      character*120 input_line
      allocatable :: eta_gb(:),dav_gb(:,:)

!     Initialize MPI environment
!     Two models share same processors
      call mpi_init(ier)
      call mpi_comm_dup(MPI_COMM_WORLD,comm1,ier) !for SELFE
      call mpi_comm_dup(MPI_COMM_WORLD,comm2,ier) !for SWAN
      call mpi_comm_dup(MPI_COMM_WORLD,comm3,ier) !for this program
      call mpi_comm_size(comm3,nproc0,ier)
      call mpi_comm_rank(comm3,myrank0,ier)

!     Test
!      call test1(comm1,myproc_all)

!     Read input file coupler_elfe_swan.in
      open(10,file='coupler_elfe_swan.in',status='old') 
      read(10,*)rnday,step_elfe,step_swan
      close(10)
      if(step_elfe>=step_swan) then
        nstep_elfe=1
        nstep_swan=step_elfe/step_swan
        nsteps=rnday*86400/step_elfe+0.5
        dt1=step_elfe; dt2=nstep_swan*step_swan; diff=abs(dt1-dt2)
      else
        nstep_swan=1
        nstep_elfe=step_swan/step_elfe
        nsteps=rnday*86400/step_swan+0.5
        dt1=step_swan; dt2=nstep_elfe*step_elfe; diff=abs(dt1-dt2)
      endif

      if(diff>1.e-5.or.step_elfe<=0.or.step_swan<=0) then
         if(myrank0==0) print*, 'COUPLER: wrong coupling step:',step_elfe,step_swan
         call mpi_abort(comm3,0,ier)
      endif

!     Read in grid info for SELFE
      open(32,file='hgrid.gr3',status='old')
      read(32,*)
      read(32,*)ne_global,np_global
      close(32)

!     Allocate exchange arrays
      allocate(eta_gb(np_global),dav_gb(2,np_global),stat=istat)

!     Run SELFE & SWAN
      istep_elfe=0 !counter for SELFE
      do iter=1,nsteps
!       SELFE (including interpolate from SELFE to SWAN grid)
        istart_elfe=(iter-1)*nstep_elfe+1
        iend_elfe=iter*nstep_elfe
        if(myrank0==0) print*, 'start ELFE',istart_elfe
        call elfe(comm1,istart_elfe,iend_elfe,np_global,eta_gb,dav_gb)
!       Print out for testing
        if(myrank0==0) print*, 'done ELFE'

!       SWAN
!       Prepare inputs
!       INPUT
        if(myrank0==0) then
          open(10,file='INPUT',status='old')
          lines=0
          do
            read(10,'(a)',end=101)input_line
            input_line=adjustl(input_line)
            len_line=len_trim(input_line)
            lines=lines+1
!            if(lines==30.and.iter/=1) write(10,'(a)')"INIT HOTS MULT 'swan_hot2'"
!"
          enddo
101       print*, 'Lines in INPUT=',lines
          close(10)
        endif !myrank0==0
        call mpi_barrier(comm3,ierr)

!       Run SWAN
        call SWAN(comm2)
        call mpi_barrier(comm3,ierr)
      enddo !iter

      call mpi_finalize(ier)

      contains
      subroutine test1(comm1,myproc_all)
      implicit none
      integer :: comm1,comm,mysize,ier,myproc,myproc_all

      comm=comm1
      call mpi_comm_size(comm1,mysize,ier)
      call mpi_comm_rank(comm1,myproc,ier)
      call mpi_allreduce(myproc,myproc_all,1,MPI_INTEGER,MPI_SUM,comm,ier)

      end subroutine test1

      end program coupler_elfe_swan


