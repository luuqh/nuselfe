!===============================================================================
!===============================================================================
! ELFE MISCELLANEOUS SUBROUTINES
!
! subroutine levels0
! subroutine nodalvel
! subroutine vinter
! subroutine eqstate
! subroutine asm
! function rint_lag
! function lindex_s 
! function covar
! subroutine cubic_spline
! subroutine eval_cubic_spline
! subroutine do_cubic_spline
! subroutine mean_density
! function kronecker
! subroutine hgrad_nodes


!===============================================================================
!===============================================================================

      subroutine levels0(iths,it)
!-------------------------------------------------------------------------------
! Routine to update level indices and wetting and drying.
! Use levels1() for better inundation if resolution is fine enough.
!-------------------------------------------------------------------------------
#ifdef USE_MPIMODULE
      use mpi
#endif
      use elfe_glbl
      use elfe_msgp
      implicit real(rkind)(a-h,o-z),integer(i-n)
#ifndef USE_MPIMODULE
      include 'mpif.h'
#endif
      integer, intent(in) :: iths,it
      integer :: idry_new(npa),idry_s_new(nsa)
      real(rkind),allocatable :: swild(:,:,:)
      logical :: srwt_xchng,prwt_xchng
      logical :: srwt_xchng_gb,prwt_xchng_gb
      logical :: cwtime
!-------------------------------------------------------------------------------

! Flag for comm timing
      cwtime=it/=iths

!...  z-coor. for nodes
!...  
      iback=0
      do i=1,npa
        if(dp(i)+eta2(i)<=h0) then !dry
          idry_new(i)=1 
          if(dp(i)>=h_s) then
            write(errmsg,*)'Deep depth dry:',i
            call parallel_abort(errmsg)
          endif
          kbp(i)=0
        else !wet
          idry_new(i)=0
!         S-levels
          do k=kz,nvrt
            kin=k-kz+1

            if(hmod(i)<=h_c) then
!              if(ifort12(12)==0) then
!                ifort12(12)=1
!                write(12,*)'Initial depth too shallow for S:',iplg(i),hmod(i),h_c
!              endif
              iback(i)=1
              z(k,i)=sigma(kin)*(hmod(i)+eta2(i))+eta2(i)
            else if(eta2(i)<=-h_c-(hmod(i)-h_c)*theta_f/s_con1) then !hmod(i)>h_c>=0
              write(errmsg,*)'Pls choose a larger h_c (1):',eta2(i),h_c
              call parallel_abort(errmsg)
            else
              z(k,i)=eta2(i)*(1+sigma(kin))+h_c*sigma(kin)+(hmod(i)-h_c)*cs(kin)
            endif
          enddo !k=kz,nvrt

!         z-levels
          if(dp(i)<=h_s) then
            kbp(i)=kz
          else !bottom index 
            if(imm==1.or.it==iths) then
              kbp(i)=0 !flag
              do k=1,kz-1
                if(-dp(i)>=ztot(k).and.-dp(i)<ztot(k+1)) then
                  kbp(i)=k
                  exit
                endif
              enddo !k
              if(kbp(i)==0) then
                write(errmsg,*)'Cannot find a bottom level for node (3):',i
!'
                call parallel_abort(errmsg)
              endif
            endif !imm

            if(kbp(i)>=kz.or.kbp(i)<1) then
              write(errmsg,*)'Impossible 92:',kbp(i),kz,i
              call parallel_abort(errmsg)
            endif
            z(kbp(i),i)=-dp(i)
            do k=kbp(i)+1,kz-1
              z(k,i)=ztot(k)
            enddo !k
          endif

          do k=kbp(i)+1,nvrt
            if(z(k,i)-z(k-1,i)<=0) then
              write(errmsg,*)'Inverted z-levels at:',i,k,z(k,i)-z(k-1,i),eta2(i),hmod(i)
              call parallel_abort(errmsg)
            endif
          enddo !k
        endif !wet ot dry
      enddo !i=1,npa

!     Debug
!      fdb='dry_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(10,file='outputs/'//fdb,status='unknown')
!      rewind(10)
!      write(10,*)'Time step=',it
!      write(10,*)'Node'
!      do i=1,npa
!        write(10,*)i,iplg(i),dp(i),eta2(i),idry_new(i)
!      enddo !i

!...  Set wet/dry flags for element; element is "dry" if one of nodes is dry; conversely, 
!...  an element is wet if all nodes are wet (and all sides are wet as well)
!...  Weed out fake wet nodes; a node is wet if and only if at least one surrounding element is wet
!...
!      if(it/=iths) idry_e0=idry_e !save only for upwindtrack()

      do i=1,nea
        idry_e(i)=max0(idry_new(nm(i,1)),idry_new(nm(i,2)),idry_new(nm(i,3)))
      enddo !i

!      write(10,*)'Element'
!      do i=1,nea
!        write(10,*)i,ielg(i),idry_e(i)
!      enddo !i

      idry_s_new(:)=idry(:) !temporary save
      idry=1 !dry unless wet
      do i=1,np
        do j=1,nne(i)
          ie=ine(i,j)
          if(idry_e(ie)==0) then
            idry(i)=0; exit
          endif
        enddo !j
      enddo !i

#ifdef INCLUDE_TIMING
      if(cwtime) cwtmp=mpi_wtime()
#endif
      call exchange_p2di(idry) !update ghost values
#ifdef INCLUDE_TIMING
      if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif

!      write(10,*)'nodes'
!      do i=1,npa
!        write(10,*)i,iplg(i),idry(i),np
!      enddo !i

!     Consistency check
#ifdef DEBUG
      do i=1,npa
        if(idry(i)==1) cycle
 
        if(eta2(i)+dp(i)<=h0) then
          write(errmsg,*)'levels0: weird wet node:',iplg(i),eta2(i),dp(i),idry_new(i)
          call parallel_abort(errmsg)
        endif

        if(i>np) cycle !do rest for residents only
        ifl=0
        do j=1,nne(i)
          ie=ine(i,j)
          if(idry_e(ie)==0) then
            ifl=1; exit
          endif 
        enddo !j
        if(ifl==0) then
          write(errmsg,*)'Node-element inconsistency:',iplg(i),idry(i),(idry_e(ine(i,j)),j=1,nne(i))
          call parallel_abort(errmsg)
        endif
      enddo !i=1,npa
#endif

!     Compute element bottom index
      kbe=0
      do i=1,nea
        if(idry_e(i)/=0) cycle

!       Wet
        n1=nm(i,1); n2=nm(i,2); n3=nm(i,3)
        if(idry(n1)/=0.or.idry(n2)/=0.or.idry(n3)/=0) then
          write(errmsg,*)'level0: Element-node inconsistency (0):',ielg(i),idry_e(i), &
                    iplg(nm(i,1:3)),idry(nm(i,1:3)),idry_new(nm(i,1:3))
          call parallel_abort(errmsg)
        endif
        kbe(i)=max0(kbp(n1),kbp(n2),kbp(n3))
        do k=kbe(i),nvrt
          ze(k,i)=(z(k,n1)+z(k,n2)+z(k,n3))/3
          if(k>=kbe(i)+1.and.ze(k,i)-ze(k-1,i)<=0) then
            write(errmsg,*)'Weird element:',k,i,ze(k,i),ze(k-1,i)
            call parallel_abort(errmsg)
          endif
        enddo !k
      enddo !i

!     Compute vel., S,T for re-wetted nodes (q2 and xl are fine)
      prwt_xchng=.false.
      if(it/=iths) then
        do i=1,npa !ghosts not updated
          if(idry_s_new(i)==1.and.idry(i)==0) then
            if(.not.prwt_xchng.and.i>np) prwt_xchng=.true. !ghost rewetted; need exchange
            if(i>np) cycle !do rest for residents

            do k=1,nvrt
              uu2(k,i)=0
              vv2(k,i)=0
              tnd(k,i)=0
              snd(k,i)=0
              icount=0
              do j=1,nnp(i)
                nd=inp(i,j) !must be inside the aug. domain
!               Wet nbrs not affected by this part and so each sub-domain should use same values
                if(idry_s_new(nd)==0) then !all indices extended
                  icount=icount+1
                  uu2(k,i)=uu2(k,i)+uu2(k,nd)
                  vv2(k,i)=vv2(k,i)+vv2(k,nd)
                  tnd(k,i)=tnd(k,i)+tnd(k,nd)
                  snd(k,i)=snd(k,i)+snd(k,nd)
                endif
              enddo !j
              if(icount==0) then
                if(ifort12(7)==0) then
                  ifort12(7)=1
                  write(12,*)'Isolated rewetted node:',iplg(i)
                endif
                tnd(k,i)=tem0(k,i)
                snd(k,i)=sal0(k,i)
              else
                uu2(k,i)=uu2(k,i)/icount
                vv2(k,i)=vv2(k,i)/icount
                tnd(k,i)=tnd(k,i)/icount
                snd(k,i)=snd(k,i)/icount
              endif
            enddo !k=1,nvrt
          endif !rewetted
        enddo !i=1,npa
      endif !it/=iths

!...  z-coor. for sides
!...  A side is wet if and only if at least one of its elements is wet
      idry_s_new=1 !reinitialize to wipe out previous temp. storage
      do i=1,ns
        do j=1,2 !elements
          ie=is(i,j)
          if(ie/=0.and.idry_e(ie)==0) idry_s_new(i)=0
        enddo !j
      enddo !i

#ifdef INCLUDE_TIMING
      if(cwtime) cwtmp=mpi_wtime()
#endif
      call exchange_s2di(idry_s_new) !update ghost values
#ifdef INCLUDE_TIMING
      if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif

!      write(10,*)'Side'
!      do i=1,nsa
!        write(10,*)i,islg(i),idry_s_new(i),ns
!      enddo !i

!     Consistency checks
#ifdef DEBUG
      do i=1,nea
        if(idry_e(i)/=0) cycle
!       Wet
        do j=1,3
          isd=js(i,j)
          if(idry_s_new(isd)/=0) then
            write(errmsg,*)'Element-side inconsistency:',ielg(i),islg(isd),idry_s_new(isd)
            call parallel_abort(errmsg)
          endif
        enddo !j
      enddo !i

      do i=1,ns
        if(idry_s_new(i)==1) cycle

        ifl=0
        do j=1,2
          ie=is(i,j)
          if(ie/=0.and.idry_e(ie)==0) then
            ifl=1; exit
          endif
        enddo !j
        if(ifl==0) then
          write(errmsg,*)'Side-element inconsistency:',islg(i),idry_s_new(i), &
                         (is(i,j),idry_e(is(i,j)),j=1,2)
          call parallel_abort(errmsg)
        endif
      enddo !i
#endif

!     Compute side bottom index
      do i=1,nsa
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        kbs(i)=0 !dry
        if(idry_s_new(i)==0) then !wet side with 2 wet nodes
          if(idry(n1)/=0.or.idry(n2)/=0) then
            write(errmsg,*)'Side-node inconsistency:',it,islg(i),'node:',iplg(n1),iplg(n2), &
              eta2(n1),eta2(n2),idry(n1),idry(n2),';element:', &
              (is(i,j),ielg(is(i,j)),idry_e(is(i,j)),j=1,2)
            call parallel_abort(errmsg)
          endif
          if(dps(i)+(eta2(n1)+eta2(n2))/2<=h0) then
            write(errmsg,*)'Weird side:',islg(i),iplg(n1),iplg(n2),eta2(n1),eta2(n2)
            call parallel_abort(errmsg)
          endif
          kbs(i)=max0(kbp(n1),kbp(n2))
          do k=kbs(i),nvrt
            zs(k,i)=(z(k,n1)+z(k,n2))/2
            if(k>=kbs(i)+1.and.zs(k,i)-zs(k-1,i)<=0) then
              write(errmsg,*)'Weird side:',k,iplg(n1),iplg(n2),z(k,n1),z(k,n2),z(k-1,n1),z(k-1,n2)
              call parallel_abort(errmsg)
            endif
          enddo !k
        endif !wet side
      enddo !i=1,nsa

!     Compute vel., S,T for re-wetted sides 
      srwt_xchng=.false.
      if(it/=iths) then
        do i=1,nsa
          if(idry_s(i)==1.and.idry_s_new(i)==0) then
            if(.not.srwt_xchng.and.i>ns) srwt_xchng=.true. !rewetted ghost side; needs exchange
            if(i>ns) cycle !do the rest only for residents

            n1=isidenode(i,1)
            n2=isidenode(i,2)
            do k=1,nvrt
              su2(k,i)=0
              sv2(k,i)=0
              tsd(k,i)=0
              ssd(k,i)=0
              icount=0
              do j=1,2
                ie=is(i,j)
                if(ie/=0) then
                  if(ie<0) call parallel_abort('levels0: ghost element')
                  do jj=1,3 !side; in the aug. domain
                    isd=js(ie,jj)
                    if(idry_s(isd)==0) then
                      icount=icount+1
                      su2(k,i)=su2(k,i)+su2(k,isd)
                      sv2(k,i)=sv2(k,i)+sv2(k,isd)
                      tsd(k,i)=tsd(k,i)+tsd(k,isd)
                      ssd(k,i)=ssd(k,i)+ssd(k,isd)
                    endif
                  enddo !jj
                endif !ie/=0
              enddo !j
              if(icount==0) then
                if(ifort12(10)==0) then
                  ifort12(10)=1
                  write(12,*)'Isolated rewetted side:',i,iplg(n1),iplg(n2)
                endif
                tsd(k,i)=(tem0(k,n1)+tem0(k,n2))/2
                ssd(k,i)=(sal0(k,n1)+sal0(k,n2))/2
              else
                su2(k,i)=su2(k,i)/icount
                sv2(k,i)=sv2(k,i)/icount
                tsd(k,i)=tsd(k,i)/icount
                ssd(k,i)=ssd(k,i)/icount
              endif
            enddo !k
          endif !rewetted
        enddo !i=1,ns
      endif !it/=iths

      idry_s=idry_s_new
     
      if(nproc>1) then
#ifdef INCLUDE_TIMING
        if(cwtime) cwtmp=mpi_wtime()
#endif
!       See if the node/side exchange is needed
        call mpi_allreduce(prwt_xchng,prwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('levels0: allreduce prwt_xchng_gb',ierr)
        call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('levels0: allreduce srwt_xchng_gb',ierr)
!'
#ifdef INCLUDE_TIMING
        if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif

!       Allocate temporary array
        if(prwt_xchng_gb.or.srwt_xchng_gb) then
          allocate(swild(4,nvrt,nsa),stat=istat)
          if(istat/=0) call parallel_abort('Levels0: fail to allocate swild')
!'
        endif

!       update ghost nodes
        if(prwt_xchng_gb) then
          swild(1,:,1:npa)=uu2(:,:)
          swild(2,:,1:npa)=vv2(:,:)
          swild(3,:,1:npa)=tnd(:,:)
          swild(4,:,1:npa)=snd(:,:)
#ifdef INCLUDE_TIMING
          if(cwtime) cwtmp=mpi_wtime()
#endif
          call exchange_p3d_4(swild)
!          call exchange_p3dw(uu2) 
!          call exchange_p3dw(vv2) 
!          call exchange_p3dw(tnd) 
!          call exchange_p3dw(snd) 
#ifdef INCLUDE_TIMING
          if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
          uu2(:,:)=swild(1,:,1:npa)
          vv2(:,:)=swild(2,:,1:npa)
          tnd(:,:)=swild(3,:,1:npa)
          snd(:,:)=swild(4,:,1:npa)
        endif

!       update ghost sides
        if(srwt_xchng_gb) then
          swild(1,:,:)=su2(:,:)
          swild(2,:,:)=sv2(:,:)
          swild(3,:,:)=tsd(:,:)
          swild(4,:,:)=ssd(:,:)
#ifdef INCLUDE_TIMING
          if(cwtime) cwtmp=mpi_wtime()
#endif
          call exchange_s3d_4(swild)
!          call exchange_s3dw(su2)
!          call exchange_s3dw(sv2)
!          call exchange_s3dw(tsd)
!          call exchange_s3dw(ssd)
#ifdef INCLUDE_TIMING
          if(cwtime) wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
          su2(:,:)=swild(1,:,:)
          sv2(:,:)=swild(2,:,:)
          tsd(:,:)=swild(3,:,:)
          ssd(:,:)=swild(4,:,:)
        endif

        if(prwt_xchng_gb.or.srwt_xchng_gb) deallocate(swild)
      endif !nproc>1

!      close(10)

      end subroutine levels0

!===============================================================================
!===============================================================================

      subroutine nodalvel(ifltype)
!-------------------------------------------------------------------------------
! Convert normal vel. to 3D nodal vel. at WHOLE levels.
!-------------------------------------------------------------------------------
#ifdef USE_MPIMODULE
      use mpi
#endif
      use elfe_glbl
      use elfe_msgp
      implicit real(rkind)(a-h,o-z),integer(i-n)
#ifndef USE_MPIMODULE
      include 'mpif.h'
#endif

      integer, intent(in) :: ifltype(max(1,nope_global))

      logical :: ltmp,ltmp2
      dimension swild(2),swild2(nvrt,2),swild3(nvrt)
      real(rkind), allocatable :: swild4(:,:,:) !swild4 used for exchange

!     Compute discontinuous hvel first (used in btrack)
      ufg=0; vfg=0
      do i=1,nea
        do k=1,nvrt
          do j=1,3
            isd0=js(i,j)
            isd1=js(i,nx(j,1))
            isd2=js(i,nx(j,2))
            ufg(k,i,j)=su2(k,isd1)+su2(k,isd2)-su2(k,isd0)
            vfg(k,i,j)=sv2(k,isd1)+sv2(k,isd2)-sv2(k,isd0)
!           Error: impose bounds for ufg, vfg?
          enddo !j
        enddo !k
      enddo !i=1,nea

!      swild=-99; swild2=-99; swild3=-99 !initialize for calling vinter
      if(indvel<=0) then
!-------------------------------------------------------------------------------
      uu2=0; vv2=0; ww2=0 !initialize and for dry nodes etc.
      do i=1,np !resident only
        if(idry(i)==1) cycle

!       Wet node
        do k=kbp(i),nvrt
          weit_w=0
          icount=0
          do j=1,nne(i)
            ie=ine(i,j)
            id=iself(i,j)
            if(idry_e(ie)==0) then
              icount=icount+1
              uu2(k,i)=uu2(k,i)+ufg(k,ie,id)
              vv2(k,i)=vv2(k,i)+vfg(k,ie,id)
            endif

            if(interpol(ie)==1) then !along Z
              if(idry_e(ie)==1) then
                swild(1)=0
              else !wet element; node i is also wet
                kbb=kbe(ie)
                swild3(kbb:nvrt)=ze(kbb:nvrt,ie) 
                swild2(kbb:nvrt,1)=we(kbb:nvrt,ie)
                call vinter(nvrt,2,1,z(k,i),kbb,nvrt,k,swild3,swild2,swild,ibelow)
              endif
            else !along S
              swild(1)=we(k,ie)
            endif !Z or S

            ww2(k,i)=ww2(k,i)+swild(1)*area(ie)
            weit_w=weit_w+area(ie)
          enddo !j
          if(icount==0) then
            write(errmsg,*)'Isolated wet node (8):',i
            call parallel_abort(errmsg)
          else
            uu2(k,i)=uu2(k,i)/icount
            vv2(k,i)=vv2(k,i)/icount
          endif
          ww2(k,i)=ww2(k,i)/weit_w
        enddo !k=kbp(i),nvrt

!       Extend
        do k=1,kbp(i)-1
          uu2(k,i)=0 !uu2(kbp(i),i) 
          vv2(k,i)=0 !vv2(kbp(i),i) 
          ww2(k,i)=0 !ww2(kbp(i),i) 
        enddo !k
      enddo !i=1,np

!-------------------------------------------------------------------------------
      else !indvel=1: averaging vel.
!-------------------------------------------------------------------------------
      uu2=0; vv2=0; ww2=0 !initialize and for dry nodes etc.
      do i=1,np !resident only
        if(idry(i)==1) cycle

!       Wet node
        icase=2
        do j=1,nne(i)
          ie=ine(i,j)
          if(interpol(ie)==1) icase=1
        enddo !j

        do k=kbp(i),nvrt
          weit=0
          weit_w=0

          do j=1,nne(i)
            ie=ine(i,j)
            id=iself(i,j)

            if(isbnd(1,i)/=0) then !bnd ball
              limit=1
            else !internal ball
              limit=2
            endif
            do l=2,limit,-1
              isd=js(ie,nx(id,l))
!             If i is on an open bnd where vel. is imposed, only the sides with imposed 
!             vel. b.c. are used in the calculation and contributions from other side are 0.
              ltmp=isbnd(1,i)>0.and.ifltype(isbnd(1,i))/=0.or. &
                   isbnd(2,i)>0.and.ifltype(isbnd(2,i))/=0
              if(ltmp) then
                nfac=0
                ltmp2=isbnd(1,i)>0.and.ifltype(isbnd(1,i))/=0.and.isbs(isd)==isbnd(1,i).or. &
                      isbnd(2,i)>0.and.ifltype(isbnd(2,i))/=0.and.isbs(isd)==isbnd(2,i)
                if(ltmp2) nfac=1
              else
                nfac=1
              endif

              if(icase==1) then !along Z
                if(idry_s(isd)==1) then
                  swild(1:2)=0
                else !wet side; node i is also wet
                  kbb=kbs(isd)
                  swild2(kbb:nvrt,1)=su2(kbb:nvrt,isd)
                  swild2(kbb:nvrt,2)=sv2(kbb:nvrt,isd)
                  swild3(kbb:nvrt)=zs(kbb:nvrt,isd)
                  call vinter(nvrt,2,2,z(k,i),kbb,nvrt,k,swild3,swild2,swild,ibelow)
                endif
              else !along S
                swild(1)=su2(k,isd)
                swild(2)=sv2(k,isd)
              endif !Z or S

              uu2(k,i)=uu2(k,i)+swild(1)/distj(isd)*nfac
              vv2(k,i)=vv2(k,i)+swild(2)/distj(isd)*nfac
              weit=weit+1/distj(isd)*nfac
            enddo !l

            if(interpol(ie)==1) then !along Z
              if(idry_e(ie)==1) then
                swild(1)=0
              else !wet element; node i is also wet
                kbb=kbe(ie)
                swild3(kbb:nvrt)=ze(kbb:nvrt,ie) 
                swild2(kbb:nvrt,1)=we(kbb:nvrt,ie)
                call vinter(nvrt,2,1,z(k,i),kbb,nvrt,k,swild3,swild2,swild,ibelow)
              endif
            else !along S
              swild(1)=we(k,ie)
            endif !Z or S

            ww2(k,i)=ww2(k,i)+swild(1)*area(ie)
            weit_w=weit_w+area(ie)
          enddo !j=1,nne(i)

          if(weit==0) then
            write(errmsg,*)'nodalvel: Isolated open bnd node:',iplg(i),isbnd(1:2,i)
            call parallel_abort(errmsg)
          endif
          uu2(k,i)=uu2(k,i)/weit
          vv2(k,i)=vv2(k,i)/weit
          ww2(k,i)=ww2(k,i)/weit_w
        enddo !k=kbp(i),nvrt

!       Extend
        do k=1,kbp(i)-1
          uu2(k,i)=0 !uu2(kbp(i),i) 
          vv2(k,i)=0 !vv2(kbp(i),i) 
          ww2(k,i)=0 !ww2(kbp(i),i) 
        enddo !k
      enddo !i=1,np
!-------------------------------------------------------------------------------
      endif !discontinous or averaging vel.

!     Exchange ghosts
      allocate(swild4(3,nvrt,npa),stat=istat)
      if(istat/=0) call parallel_abort('nodalvel: fail to allocate')
      swild4(1,:,:)=uu2(:,:)
      swild4(2,:,:)=vv2(:,:)
      swild4(3,:,:)=ww2(:,:)
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_p3d_3(swild4)
!      call exchange_p3dw(uu2)
!      call exchange_p3dw(vv2)
!      call exchange_p3dw(ww2)
#ifdef INCLUDE_TIMING
      wtimer(10,2)=wtimer(10,2)+mpi_wtime()-cwtmp
#endif
      uu2(:,:)=swild4(1,:,:)
      vv2(:,:)=swild4(2,:,:)
      ww2(:,:)=swild4(3,:,:)
      deallocate(swild4)

!...  Compute discrepancy between avergaed and elemental vel. vectors 
!      do i=1,np
!	do k=1,nvrt
!	  testa(i,k)=0
!          do j=1,nne(i)
!	    iel=ine(i,j)
!	    index=0
!	    do l=1,3
!	      if(nm(iel,l).eq.i) index=l
!	    enddo !l
!	    if(index.eq.0) then
!	      write(*,*)'Wrong element ball'
!	      stop
!	    endif
!	    testa(i,k)=testa(i,k)+sqrt((uuf(iel,index,k)-uu2(k,i))**2+
!     +(vvf(iel,index,k)-vv2(k,i))**2)/nne(i)
!	  enddo !j
!	enddo !k
!      enddo !i

      end subroutine nodalvel

!===============================================================================
!===============================================================================

      subroutine vinter(nmax1,nmax2,nc,zt,k1,k2,k3,za,sint,sout,ibelow)
!     Routine to do vertical linear interpolation in z
!     Inputs:
!       (nmax1,nmax2) : dimension of sint() in the calling routine
!       nc: actual # of variables
!       k1,k2: lower and upper limits for za, sint
!       k3: initial guess for level index (to speed up)
!       zt: desired interpolation level
!       za(k1:k2): z-cor for sint (must be in ascending order)
!       sint(k1:k2,1:nc): values to be interpolated from; dimensions must match driving program
!     Outputs:
!       sout(1:nc): interpolated value @ z=zt (bottom value if ibelow=1). Constant extrapolation
!                   is used below bottom or above surface.
!       ibelow: flag indicating if zt is below za(k1)
!
!  TODO: change index order for sint()
      use elfe_glbl, only : rkind,errmsg
      use elfe_msgp, only : parallel_abort
      implicit real(rkind)(a-h,o-z), integer(i-n)
      integer, intent(in) :: nmax1,nmax2,nc,k1,k2,k3
      real(rkind), intent(in) :: zt,za(nmax1),sint(nmax1,nmax2)
      real(rkind), dimension(:), intent(out) :: sout(nmax2)
      integer, intent(out) :: ibelow
      logical :: first_call

      first_call=.true.

      if(k1>k2) then !.or.nc>10) then
        write(errmsg,*)'k1>k2 in vinter()'
        call parallel_abort(errmsg)
      endif

      if(zt<za(k1)) then
        ibelow=1
        sout(1:nc)=sint(k1,1:nc)
      else !normal
        ibelow=0
        if(zt==za(k1)) then
          sout(1:nc)=sint(k1,1:nc)
        else if(zt>=za(k2)) then
          sout(1:nc)=sint(k2,1:nc)
        else
          kout=0 !flag
          if(k3<k1.or.k3>k2) then
            l1=k1; l2=k2-1
          else
            if(zt<za(k3)) then
              l1=k1; l2=k3-1
            else
              l1=k3; l2=k2-1
            endif
          endif
          do k=l1,l2
            if(zt>=za(k).and.zt<=za(k+1)) then
              kout=k
              exit
            endif
          enddo !k
          if(kout==0.or.za(kout+1)-za(kout)==0) then
            write(errmsg,*)'Failed to find a level in vinter():',kout,zt,(za(k),k=k1,k2)
            call parallel_abort(errmsg)
          endif
          zrat=(zt-za(kout))/(za(kout+1)-za(kout))
          sout(1:nc)=sint(kout,1:nc)*(1-zrat)+sint(kout+1,1:nc)*zrat
        endif
      endif

      first_call=.false.
      end subroutine vinter

!===============================================================================
!===============================================================================
!
!***************************************************************************
!									   *
!     Solve for the density
!     From Pond and Pickard's book.					   *
!     validity region: T: [0,40], S: [0:42]				   *
!     Inputs: 
!            tem,sal: T,S (assumed to be at wet spots).
!     Output: density.
!									   *
!***************************************************************************
!   
      function eqstate(tem2,sal2)
      use elfe_glbl, only: rkind,tempmin,tempmax,saltmin,saltmax,errmsg,ifort12
      use elfe_msgp, only : parallel_abort
      implicit real(rkind)(a-h,o-z),integer(i-n)
      real(rkind), intent(in) :: tem2,sal2

      den(t,s)=1000-0.157406+6.793952E-2*t-9.095290E-3*t**2 &
     &+1.001685E-4*t**3-1.120083E-6*t**4+6.536332E-9*t**5+ &
     &s*(0.824493-4.0899E-3*t+&
     &7.6438E-5*t**2-8.2467E-7*t**3+5.3875E-9*t**4)+&
     &sqrt(s)**3*(-5.72466E-3+1.0227E-4*t-1.6546E-6*t**2)+&
     &4.8314E-4*s**2

      tem=tem2; sal=sal2
      if(tem<-98.or.sal<-98) then
        write(errmsg,*)'EQSTATE: Impossible dry (7):',tem,sal
        call parallel_abort(errmsg)
      endif
      if(tem<tempmin.or.tem>tempmax.or.sal<saltmin.or.sal>saltmax) then
        if(ifort12(6)==0) then
          ifort12(6)=1
          write(12,*)'Invalid temp. or salinity for density:',tem,sal
        endif
        tem=max(tempmin,min(tem,tempmax))
        sal=max(saltmin,min(sal,saltmax))
      endif

!     Density at one standard atmosphere
      eqstate=den(tem,sal)
      if(eqstate<980) then
        write(errmsg,*)'Weird density:',eqstate,tem,sal
        call parallel_abort(errmsg)
      endif

      end function eqstate

!===============================================================================
!===============================================================================
      subroutine asm(g,i,j,vd,td,qd1,qd2)
!     Algebraic Stress Models
      use elfe_glbl
      use elfe_msgp, only : parallel_abort
      implicit real(rkind)(a-h,o-z), integer(i-n)
      integer, intent(in) :: i,j
      real(rkind), intent(in) :: g
      real(rkind), intent(out) :: vd,td,qd1,qd2

      if(j<kbp(i).or.j>nvrt) then
        write(errmsg,*)'Wrong input level:',j
        call parallel_abort(errmsg)
      endif

!     Wet node i with rho defined; kbp(i)<=j<=nvrt
      if(j==kbp(i).or.j==nvrt) then
        drho_dz=0
      else
        drho_dz=(prho(j+1,i)-prho(j-1,i))/(z(j+1,i)-z(j-1,i))
      endif
      bvf=g/rho0*drho_dz
      Gh=xl(i,j)**2/2/q2(i,j)*bvf
      Gh=min(max(Gh,-0.28_rkind),0.0233_rkind)

      if(stab.eq.'GA') then
        sh=0.49393/(1-34.676*Gh)
        sm=(0.39327-3.0858*Gh)/(1-34.676*Gh)/(1-6.1272*Gh)
        cmiu=sqrt(2.d0)*sm
        cmiup=sqrt(2.d0)*sh
        cmiu1=sqrt(2.d0)*0.2 !for k-eq
        cmiu2=sqrt(2.d0)*0.2 !for psi-eq.
      else if(stab.eq.'KC') then !Kantha and Clayson
!       Error: Warner's paper wrong
!        Ghp=(Gh-(Gh-0.02)**2)/(Gh+0.0233-0.04) !smoothing
        Ghp=Gh
        sh=0.4939/(1-30.19*Ghp)
        sm=(0.392+17.07*sh*Ghp)/(1-6.127*Ghp)
        cmiu=sqrt(2.d0)*sm
        cmiup=sqrt(2.d0)*sh
        cmiu1=cmiu/schk
        cmiu2=cmiu/schpsi
      else
        write(errmsg,*)'Unknown ASM:',mid
        call parallel_abort(errmsg)
      endif

      vd=cmiu*xl(i,j)*sqrt(q2(i,j))
      td=cmiup*xl(i,j)*sqrt(q2(i,j))
      qd1=cmiu1*xl(i,j)*sqrt(q2(i,j))
      qd2=cmiu2*xl(i,j)*sqrt(q2(i,j))

      end subroutine asm

!===============================================================================
!===============================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!											!
!    Generic routine to compute \int_{\sigma_k}^{\sigma_{k+1}} \psi(\sigma)d\sigma,	!
!    where Nmin<=k<=Nmax-1, \sigma & \psi(Nmin:Nmax), using Lagrangian  		!
!    interpolation of order 2*m (i.e., from k-m to k+m).				!
!    mnv: dimensioning parameter from driving routine (input);				!
!    Nmin, Nmax: limits of vertical levels (input);					!
!    m: order of Lagrangian polynormial (input);					!
!    k: input for limits;								!
!    sigma,sigmap,sigma_prod,psi: input (sigmap&sigma_prod are the pre-computed 	!
!                                  powers and products of sigma for speed)		!
!    gam, coef: working arrays (output).						!
!    WARNING: Nmax must =nsig, and 1<=Nmin<=nsig-1 for sigma_prod!!			!
!											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      function rint_lag(mnv,Nmin,Nmax,m,k,sigma,sigmap,sigma_prod,psi,gam,coef)
      use elfe_glbl, only : rkind,errmsg
      use elfe_msgp, only : parallel_abort
      implicit real(rkind)(a-h,o-z), integer(i-n)
      integer, intent(in) :: mnv,Nmin,Nmax,m,k
      real(rkind), intent(in) :: sigma(mnv),sigmap(mnv,10),sigma_prod(mnv,mnv,-4:4),psi(mnv)
      real(rkind), intent(out) :: gam(mnv),coef(0:mnv)

!     Sanity check
      if(Nmin>=Nmax.or.Nmax>mnv.or.Nmin<1) then
        write(errmsg,*)'Check inputs in rint_lag:',Nmin,Nmax
        call parallel_abort(errmsg)
      endif
      if(k>Nmax-1.or.k<Nmin) then
        write(errmsg,*)'Wrong k:',k
        call parallel_abort(errmsg)
      endif
      if(m<1) then
        write(errmsg,*)'m<1',m
        call parallel_abort(errmsg)
      endif
      if(m>3) then
        write(errmsg,*)'m>3 not covered presently' 
        call parallel_abort(errmsg)
      endif
      if(2*m+1>10) then
        write(errmsg,*)'Re-dimension sigmap'
        call parallel_abort(errmsg)
      endif

!     Compute J1,2
      j1=max0(Nmin,k-m)
      j2=min0(Nmax,k+m)
      if(j1>=j2) then
         write(errmsg,*)'Weird indices:',j1,j2
         call parallel_abort(errmsg)
      endif

!     Compute sum
      rint_lag=0
      do i=j1,j2
!       Denominator & assemble working array gam
!        prod=1
        id=0
        do j=j1,j2
          if(j/=i) then
            id=id+1
            gam(id)=-sigma(j)
          endif
        enddo !j
        if(id/=j2-j1.or.id>2*m) then
          write(errmsg,*)'Miscount:',id,j2-j1,m
          call parallel_abort(errmsg)
        endif

!       Inner sum
        if(id==1) then
          coef(0)=gam(1); coef(1)=1
        else if(id==2) then
          coef(0)=gam(1)*gam(2)
          coef(1)=gam(1)+gam(2)
          coef(2)=1
        else if(id==3) then
          coef(0)=gam(1)*gam(2)*gam(3)
          coef(1)=gam(1)*(gam(2)+gam(3))+gam(2)*gam(3)
          coef(2)=gam(1)+gam(2)+gam(3)
          coef(3)=1
        else if(id==4) then
          coef(0)=gam(1)*gam(2)*gam(3)*gam(4)
          coef(1)=gam(1)*gam(2)*(gam(3)+gam(4))+(gam(1)+gam(2))*gam(3)*gam(4)
          coef(2)=gam(1)*(gam(2)+gam(3))+(gam(1)+gam(3))*gam(4)+gam(2)*(gam(3)+gam(4))
!          coef(2)=gam(1)*gam(2)+gam(1)*gam(3)+gam(1)*gam(4)+gam(2)*gam(3)+gam(2)*gam(4)+gam(3)*gam(4)
          coef(3)=gam(1)+gam(2)+gam(3)+gam(4)
          coef(4)=1
        else if(id==5) then
          coef(0)=gam(1)*gam(2)*gam(3)*gam(4)*gam(5)
          coef(1)=gam(1)*gam(2)*gam(3)*gam(4)+gam(1)*gam(2)*gam(3)*gam(5)+gam(1)*gam(2)*gam(4)*gam(5)+ &
     &gam(1)*gam(3)*gam(4)*gam(5)+gam(2)*gam(3)*gam(4)*gam(5)
          coef(2)=gam(1)*gam(2)*gam(3)+gam(1)*gam(2)*gam(4)+gam(1)*gam(2)*gam(5)+gam(1)*gam(3)*gam(4)+ &
     &gam(1)*gam(3)*gam(5)+gam(1)*gam(4)*gam(5)+gam(2)*gam(3)*gam(4)+gam(2)*gam(3)*gam(5)+ &
     &gam(2)*gam(4)*gam(5)+gam(3)*gam(4)*gam(5)
          coef(3)=gam(1)*gam(2)+gam(1)*gam(3)+gam(1)*gam(4)+gam(1)*gam(5)+gam(2)*gam(3)+ &
     &gam(2)*gam(4)+gam(2)*gam(5)+gam(3)*gam(4)+gam(3)*gam(5)+gam(4)*gam(5)
          coef(4)=gam(1)+gam(2)+gam(3)+gam(4)+gam(5)
          coef(5)=1
        else if(id==6) then
          coef(0)=gam(1)*gam(2)*gam(3)*gam(4)*gam(5)*gam(6)
          coef(1)=gam(1)*gam(2)*gam(3)*gam(4)*gam(5)+gam(1)*gam(2)*gam(3)*gam(4)*gam(6)+&
     &gam(1)*gam(2)*gam(3)*gam(5)*gam(6)+gam(1)*gam(2)*gam(4)*gam(5)*gam(6)+ &
     &gam(1)*gam(3)*gam(4)*gam(5)*gam(6)+gam(2)*gam(3)*gam(4)*gam(5)*gam(6)
          coef(2)=gam(1)*gam(2)*gam(3)*gam(4)+gam(1)*gam(2)*gam(3)*gam(5)+gam(1)*gam(2)*gam(3)*gam(6)+ &
     &gam(1)*gam(2)*gam(4)*gam(5)+gam(1)*gam(2)*gam(4)*gam(6)+gam(1)*gam(2)*gam(5)*gam(6)+ &
     &gam(1)*gam(3)*gam(4)*gam(5)+gam(1)*gam(3)*gam(4)*gam(6)+gam(1)*gam(3)*gam(5)*gam(6)+ &
     &gam(1)*gam(4)*gam(5)*gam(6)+gam(2)*gam(3)*gam(4)*gam(5)+gam(2)*gam(3)*gam(4)*gam(6)+ &
     &gam(2)*gam(3)*gam(5)*gam(6)+gam(2)*gam(4)*gam(5)*gam(6)+gam(3)*gam(4)*gam(5)*gam(6)
           coef(3)=gam(1)*gam(2)*gam(3)+gam(1)*gam(2)*gam(4)+gam(1)*gam(2)*gam(5)+ &
     &gam(1)*gam(2)*gam(6)+gam(1)*gam(3)*gam(4)+gam(1)*gam(3)*gam(5)+gam(1)*gam(3)*gam(6)+ &
     &gam(1)*gam(4)*gam(5)+gam(1)*gam(4)*gam(6)+gam(1)*gam(5)*gam(6)+gam(2)*gam(3)*gam(4)+ &
     &gam(2)*gam(3)*gam(5)+gam(2)*gam(3)*gam(6)+gam(2)*gam(4)*gam(5)+gam(2)*gam(4)*gam(6)+ &
     &gam(2)*gam(5)*gam(6)+gam(3)*gam(4)*gam(5)+gam(3)*gam(4)*gam(6)+gam(3)*gam(5)*gam(6)+ &
     &gam(4)*gam(5)*gam(6)
           coef(4)=gam(1)*gam(2)+gam(1)*gam(3)+gam(1)*gam(4)+gam(1)*gam(5)+gam(1)*gam(6)+ &
     &gam(2)*gam(3)+gam(2)*gam(4)+gam(2)*gam(5)+gam(2)*gam(6)+gam(3)*gam(4)+gam(3)*gam(5)+ &
     &gam(3)*gam(6)+gam(4)*gam(5)+gam(4)*gam(6)+gam(5)*gam(6)
           coef(5)=gam(1)+gam(2)+gam(3)+gam(4)+gam(5)+gam(6)
           coef(6)=1
        else
          write(errmsg,*)'Not covered:',id
          call parallel_abort(errmsg)
        endif

        sum1=0
        do l=0,id
          sum1=sum1+coef(l)/(l+1)*(sigmap(k+1,l+1)-sigmap(k,l+1))
        enddo !l

        if(abs(i-k)>4) then
          write(errmsg,*)'sigma_prod index out of bound (2)'
          call parallel_abort(errmsg)
        endif

        rint_lag=rint_lag+psi(i)/sigma_prod(Nmin,k,i-k)*sum1
      enddo !i=j1,j2

      end function rint_lag

!     Compute local index of a side (0 if inside aug. domain)
      function lindex_s(i,ie)
      use elfe_glbl, only : rkind,js
      implicit real(rkind)(a-h,o-z), integer(i-n)
      integer, intent(in) :: i,ie

      l0=0 !local index
      do l=1,3
        if(js(ie,l)==i) then
          l0=l
          exit
        endif
      enddo !l
      lindex_s=l0

      end function lindex_s

      function covar(kr_co,hh)
      use elfe_glbl, only : rkind,errmsg
      use elfe_msgp, only : parallel_abort
      implicit real(rkind)(a-h,o-z), integer(i-n)

      if(hh<0) then
        write(errmsg,*)'Negative hh in covar:',hh
        call parallel_abort(errmsg) 
      endif

      if(kr_co==1) then
        covar=-hh
      else if(kr_co==2) then
        if(hh==0) then
          covar=0
        else
          covar=hh*hh*log(hh)
        endif
      else if(kr_co==3) then !cubic
        covar=hh*hh*hh
      else if(kr_co==4) then !5th
        h2=hh*hh
        covar=-h2*h2*hh
      else
        write(errmsg,*)'Unknown covariance function option:',kr_co
        call parallel_abort(errmsg)
      endif

      end function covar

!===============================================================================
!     Do interpolation with cubic spline
!     Needs coefficients from routine cubic_spline()
!     Inputs: 
!            npts: dimesion for xcor etc; # of data pts;
!            xcor(npts),yy(npts): x and y coordinates of the original function 
!                                 (same as in cubic_spline()); xcor in ascending order;
!            ypp(npts): 2nd deriavtives (output from cubic_spline);
!            npts2: # of output pts;
!            xout(npts2): x coordinates of the output pts (no ordering required);
!            xmax: if xout>xmax, it is reset to xmax;
!            ixmin (0 or 1): bottom option;
!            xmin: if xout<xcor(1), it is either reset to xmin (ixmin=0), or 
!                  to xcor(1) (ixmin=1), i.e. yyout takes the value of yy(1), and
!                  xmin is not used except for debugging messages.
!     Output: 
!            yyout(npts2): output y values; if xmin>xmax, yyout=yy(1).
!===============================================================================
      subroutine eval_cubic_spline(npts,xcor,yy,ypp,npts2,xout,ixmin,xmin,xmax,yyout)
      use elfe_glbl, only : rkind,errmsg
      use elfe_msgp, only : parallel_abort
      implicit real(rkind)(a-h,o-z), integer(i-n)
      integer, intent(in) :: npts,npts2,ixmin
      real(rkind), intent(in) :: xcor(npts),yy(npts),ypp(npts),xout(npts2),xmin,xmax
      real(rkind), intent(out) :: yyout(npts2)

      if(xmin>xmax) then
!        write(errmsg,*)'EVAL_CUBIC: xmin>xmax:',xmin,xmax
!        call parallel_abort(errmsg)
        yyout=yy(1); return
      endif

      do i=1,npts2
        ifl=0 !flag
        xtmp=min(xout(i),xmax)
        if(ixmin==0) then
          xtmp=max(xtmp,xmin)
        else
          if(xout(i)<xcor(1)) then
            yyout(i)=yy(1); cycle
          endif
        endif

        do j=1,npts-1
          if(xtmp>=xcor(j).and.xtmp<=xcor(j+1)) then
            ifl=1
            aa=(xcor(j+1)-xtmp)/(xcor(j+1)-xcor(j))
            bb=1-aa
            cc=(aa*aa*aa-aa)*(xcor(j+1)-xcor(j))/6
            dd=(bb*bb*bb-bb)*(xcor(j+1)-xcor(j))/6
            yyout(i)=aa*yy(j)+bb*yy(j+1)+cc*ypp(j)+dd*ypp(j+1)
            exit
          endif
        enddo !j
        if(ifl==0) then
          write(errmsg,*)'EVAL_CUBIC: Falied to find:',i,xtmp,xmin,xmax
          call parallel_abort(errmsg)
        endif
      enddo !i=1,npts2

      end subroutine eval_cubic_spline

!===============================================================================
!     Generate coefficients (2nd derivatives) for cubic spline for interpolation later
!     Inputs: 
!            npts: dimesion for xcor etc; # of data pts;
!            xcor(npts): x coordinates; must be in ascending order (and distinctive);
!            yy(npts): functional values; 
!            yp1 and yp2: 1st derivatives at xcor(1) and xcor(npts);
!     Output: 
!            ypp(npts): 2nd deriavtives used in interpolation.
!===============================================================================
      subroutine cubic_spline(npts,xcor,yy,yp1,yp2,ypp)
      use elfe_glbl, only : rkind,errmsg
      use elfe_msgp, only : parallel_abort
      implicit real(rkind)(a-h,o-z), integer(i-n)
      integer, intent(in) :: npts
      real(rkind), intent(in) :: xcor(npts),yy(npts),yp1,yp2
      real(rkind), intent(out) :: ypp(npts)
  
      real(rkind) :: alow(npts),bdia(npts),cupp(npts),rrhs(npts,1),gam(npts)

      do k=1,npts
        if(k==1) then
          bdia(k)=(xcor(k+1)-xcor(k))/3
          if(bdia(k)==0) then
            write(errmsg,*)'CUBIC_SP: bottom problem:',xcor(k+1),xcor(k)
            call parallel_abort(errmsg)
          endif
          cupp(k)=bdia(k)/2
          rrhs(k,1)=(yy(k+1)-yy(k))/(xcor(k+1)-xcor(k))-yp1
        else if(k==npts) then
          bdia(k)=(xcor(k)-xcor(k-1))/3
          if(bdia(k)==0) then
            write(errmsg,*)'CUBIC_SP: surface problem:',xcor(k),xcor(k-1)
            call parallel_abort(errmsg)
          endif
          alow(k)=bdia(k)/2
          rrhs(k,1)=-(yy(k)-yy(k-1))/(xcor(k)-xcor(k-1))+yp2
        else
          bdia(k)=(xcor(k+1)-xcor(k-1))/3
          alow(k)=(xcor(k)-xcor(k-1))/6
          cupp(k)=(xcor(k+1)-xcor(k))/6
          if(alow(k)==0.or.cupp(k)==0) then
            write(errmsg,*)'CUBIC_SP: middle problem:',xcor(k),xcor(k-1),xcor(k+1)
            call parallel_abort(errmsg)
          endif
          rrhs(k,1)=(yy(k+1)-yy(k))/(xcor(k+1)-xcor(k))-(yy(k)-yy(k-1))/(xcor(k)-xcor(k-1))
        endif
      enddo !k

      call tridag(npts,1,npts,1,alow,bdia,cupp,rrhs,ypp,gam)
!      ypp(:)=soln(:,1)

      end subroutine cubic_spline

!===============================================================================
!     Do cubic spline with 1 step, i.e., combining cubic_spline and eval_cubic_spline.
!     Inputs: 
!            npts: dimesion for xcor etc; # of data pts;
!            xcor(npts): x coordinates; must be in ascending order (and distinctive);
!            yy(npts): functional values; 
!            yp1 and yp2: 1st derivatives at xcor(1) and xcor(npts);
!            npts2: # of output pts;
!            xout(npts2): x coordinates of the output pts (no ordering required);
!            xmax: if xout>xmax, it is reset to xmax;
!            ixmin (0 or 1): bottom option;
!            xmin: if xout<xcor(1), it is either reset to xmin (ixmin=0), or 
!                  to xcor(1) (ixmin=1), i.e. yyout takes the value of yy(1), and
!                  xmin is not used except for debugging messages.
!     Output: 
!            yyout(npts2): output y values
!===============================================================================
      subroutine do_cubic_spline(npts,xcor,yy,yp1,yp2,npts2,xout,ixmin,xmin,xmax,yyout)
      use elfe_glbl, only : rkind,errmsg
      use elfe_msgp, only : parallel_abort
      implicit real(rkind)(a-h,o-z), integer(i-n)
      integer, intent(in) :: npts,npts2,ixmin
      real(rkind), intent(in) :: xcor(npts),yy(npts),yp1,yp2,xout(npts2),xmin,xmax
      real(rkind), intent(out) :: yyout(npts2)
 
      real(rkind) :: ypp(npts)

      call cubic_spline(npts,xcor,yy,yp1,yp2,ypp)
      call eval_cubic_spline(npts,xcor,yy,ypp,npts2,xout,ixmin,xmin,xmax,yyout)

      end subroutine do_cubic_spline

!===============================================================================
!     Compute mean density (rho_mean) at nodes (whole levels) or elements (half levels) 
!===============================================================================
      subroutine mean_density
      use elfe_glbl
      use elfe_msgp, only : parallel_abort
      implicit real(rkind)(a-h,o-z), integer(i-n)
!      integer, intent(in) :: iupwind_t
      real(rkind) :: swild(nvrt) !,swild2(nvrt,nea,2)
      real(rkind), allocatable :: swild2(:,:,:)

      allocate(swild2(nvrt,nea,2),stat=istat)

      rho_mean=-99
      if(iupwind_t==0) then !ELM
!       T,S @ nodes
        do i=1,npa
          if(idry(i)==1) cycle

!         Wet nodes
!         Extrapolation used above surface
          if(z(kbp(i),i)<z_r(1)) then !.or.z(nvrt,i)>z_r(nz_r)) then
            call parallel_abort('MISC: 2.node depth too big for ts.ic')
          endif 
!         Use swild2 for temporarily saving T,S
          call eval_cubic_spline(nz_r,z_r,tem1,cspline_ypp(1:nz_r,1),nvrt-kbp(i)+1,z(kbp(i):nvrt,i), &
     &0,z_r(1),z_r(nz_r),swild2(kbp(i):nvrt,i,1))
          call eval_cubic_spline(nz_r,z_r,sal1,cspline_ypp(1:nz_r,2),nvrt-kbp(i)+1,z(kbp(i):nvrt,i), &
     &0,z_r(1),z_r(nz_r),swild2(kbp(i):nvrt,i,2))

!         Impose no slip b.c. to be consistent with ELM transport
          if(Cdp(i)/=0) then
            swild2(kbp(i),i,1:2)=swild2(kbp(i)+1,i,1:2)
          endif

!         Extend
          do k=1,kbp(i)-1
            swild2(k,i,1:2)=swild2(kbp(i),i,1:2)
          enddo !k

!         Whole levels
          do k=1,nvrt
            rho_mean(k,i)=eqstate(swild2(k,i,1),swild2(k,i,2))
          enddo !k
        enddo !i=1,npa

      else !upwind
!       T,S @ elements
        do i=1,nea
          if(idry_e(i)==1) cycle

!         Wet element
          if(ze(kbe(i),i)<z_r(1)) then !.or.ze(nvrt,i)>z_r(nz_r)) then
            call parallel_abort('MISC: 2.ele. depth too big for ts.ic')
          endif 

          do k=kbe(i)+1,nvrt
            swild(k)=(ze(k,i)+ze(k-1,i))/2
          enddo !k
          call eval_cubic_spline(nz_r,z_r,tem1,cspline_ypp(1:nz_r,1),nvrt-kbe(i),swild(kbe(i)+1:nvrt), &
     &0,z_r(1),z_r(nz_r),swild2(kbe(i)+1:nvrt,i,1))
          call eval_cubic_spline(nz_r,z_r,sal1,cspline_ypp(1:nz_r,2),nvrt-kbe(i),swild(kbe(i)+1:nvrt), &
     &0,z_r(1),z_r(nz_r),swild2(kbe(i)+1:nvrt,i,2))

!         Extend
          do k=1,kbe(i)
            swild2(k,i,1:2)=swild2(kbe(i)+1,i,1:2)
          enddo !k

!         Half levels
          do k=1,nvrt
            rho_mean(k,i)=eqstate(swild2(k,i,1),swild2(k,i,2))
          enddo !k
        enddo !i=1,nea
      endif !ELM or upwind

      deallocate(swild2)

      end subroutine mean_density

!     Kronecker delta
      function kronecker(i,j)
      implicit integer(i-n)
      integer, intent(in) :: i,j

      if(i==j) then
        kronecker=1
      else
        kronecker=0
      endif

      end function kronecker

!===============================================================================
!     Calculate horizontal gradient at (resident) sides and whole level for variable
!     defined at nodes, using cubic spline
!===============================================================================
      subroutine hgrad_nodes(ihbnd,nvrt1,npa1,nsa1,var_nd,dvar_dxy)
      use elfe_glbl
      use elfe_msgp, only : parallel_abort
      implicit real(rkind)(a-h,o-z), integer(i-n)
      integer, intent(in) :: ihbnd !flag (0: no flux b.c. for horizontal bnd side; 1: use shape function)
      integer, intent(in) :: nvrt1,npa1,nsa1 !dimension parameters (=nvrt,npa,nsa)
      real(rkind), intent(in) :: var_nd(nvrt1,npa1) !variable defined at nodes and whole levels
      real(rkind), intent(out) :: dvar_dxy(2,nvrt1,nsa1) !only resident sides are defined (1: x-derivative; 2: y-derivative)

      !real(rkind) :: hp_int(nvrt1,npa1),swild(nvrt1),swild2(nvrt1,4),nwild(3)
      real(rkind) :: hp_int(nvrt1,npa1),swild(nvrt1),swild2(nvrt1,4)!HUY modified for integer variable
      integer :: nwild(3)
      
      dvar_dxy=0
      hp_int=0 !temporary save of 2nd deriavtives
      do i=1,npa
        if(idry(i)==1) cycle

        call cubic_spline(nvrt-kbp(i)+1,z(kbp(i):nvrt,i),var_nd(kbp(i):nvrt,i),0._rkind,0._rkind,swild)
        hp_int(kbp(i):nvrt,i)=swild(1:(nvrt-kbp(i)+1))
      enddo !i=1,npa

      do i=1,ns
        if(idry_s(i)==1) cycle

!       Wet side; pts 1&2
        node1=isidenode(i,1); node2=isidenode(i,2)
        eta_min=min(z(nvrt,node1),z(nvrt,node2))
        zmax=max(z(kbp(node1),node1),z(kbp(node2),node2)) !for bottom option
        if(-zmax>h_bcc1) then !deep sea
          ibot_fl=0
        else !shallow
          ibot_fl=1
        endif
!       Currently bounds not enforced
        call eval_cubic_spline(nvrt-kbp(node1)+1,z(kbp(node1):nvrt,node1),var_nd(kbp(node1):nvrt,node1), &
     &hp_int(kbp(node1):nvrt,node1),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
        swild2(kbs(i):nvrt,1)=swild(1:(nvrt-kbs(i)+1))
        call eval_cubic_spline(nvrt-kbp(node2)+1,z(kbp(node2):nvrt,node2),var_nd(kbp(node2):nvrt,node2), &
     &hp_int(kbp(node2):nvrt,node2),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
        swild2(kbs(i):nvrt,2)=swild(1:(nvrt-kbs(i)+1))

        !pts 3&4
        if(is(i,2)==0.and.ihbnd==0) then !no flux b.c.
          swild2(kbs(i):nvrt,3:4)=0
          x43=y(node2)-y(node1)
          y43=x(node1)-x(node2)
        else if(is(i,2)==0.and.ihbnd/=0) then !use shape function
          ie=is(i,1)
          node3=sum(nm(ie,1:3))-node1-node2
          if(idry(node3)==1) then
            write(errmsg,*)'hgrad_nodes: node3 dry',iplg(node3),ielg(ie)
            call parallel_abort(errmsg)
          endif
          !Find local indices
          nwild=0
          do j=1,3
            if(j<=2) then
              nd=isidenode(i,j)
            else
              nd=node3
            endif
            do jj=1,3
              if(nm(ie,jj)==nd) then
                nwild(j)=jj; exit
              endif
            enddo !jj
            if(nwild(j)==0) then
              write(errmsg,*)'hgrad_nodes: no index found:',iplg(nd),ielg(ie)
              call parallel_abort(errmsg)
            endif
          enddo !j
          eta_min=z(nvrt,node3)
          zmax=z(kbp(node3),node3)
          if(-zmax>h_bcc1) then !deep sea
            ibot_fl=0
          else !shallow
            ibot_fl=1
          endif
          call eval_cubic_spline(nvrt-kbp(node3)+1,z(kbp(node3):nvrt,node3),var_nd(kbp(node3):nvrt,node3), &
     &hp_int(kbp(node3):nvrt,node3),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
          swild2(kbs(i):nvrt,3)=swild(1:(nvrt-kbs(i)+1))
          do k=kbs(i),nvrt
            do j=1,3
              dvar_dxy(1:2,k,i)=dvar_dxy(1:2,k,i)+swild2(k,j)*dl(ie,nwild(j),1:2)
            enddo !j
          enddo !k          

        else !internal side
          node3=sum(nm(is(i,1),1:3))-node1-node2
          node4=sum(nm(is(i,2),1:3))-node1-node2
          x43=x(node4)-x(node3)
          y43=y(node4)-y(node3)
          if(idry(node3)==1.or.idry(node4)==1) then
            swild2(kbs(i):nvrt,3:4)=0
          else !both wet
            eta_min=min(z(nvrt,node3),z(nvrt,node4))
            zmax=max(z(kbp(node3),node3),z(kbp(node4),node4)) !for bottom option
            if(-zmax>h_bcc1) then !deep sea
              ibot_fl=0
            else !shallow
              ibot_fl=1
            endif

            call eval_cubic_spline(nvrt-kbp(node3)+1,z(kbp(node3):nvrt,node3),var_nd(kbp(node3):nvrt,node3), &
     &hp_int(kbp(node3):nvrt,node3),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
            swild2(kbs(i):nvrt,3)=swild(1:(nvrt-kbs(i)+1))
            call eval_cubic_spline(nvrt-kbp(node4)+1,z(kbp(node4):nvrt,node4),var_nd(kbp(node4):nvrt,node4), &
     &hp_int(kbp(node4):nvrt,node4),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ibot_fl,zmax,eta_min,swild)
            swild2(kbs(i):nvrt,4)=swild(1:(nvrt-kbs(i)+1))
          endif
        endif !bnd side or not

        if(ihbnd==0.or.is(i,2)/=0) then
          delta1=(x(node2)-x(node1))*y43-x43*(y(node2)-y(node1))
          if(delta1==0) then
            write(errmsg,*)'hgrad_nodes failure:',iplg(node1),iplg(node2)
            call parallel_abort(errmsg)
          endif
          do k=kbs(i),nvrt
            dvar_dxy(1,k,i)=(y43*(swild2(k,2)-swild2(k,1))-(y(node2)-y(node1))*(swild2(k,4)-swild2(k,3)))/delta1
            dvar_dxy(2,k,i)=((x(node2)-x(node1))*(swild2(k,4)-swild2(k,3))-x43*(swild2(k,2)-swild2(k,1)))/delta1
          enddo !k
        endif !ihbnd==0.or.is(i,2)/=0
      enddo !i=1,ns

      end subroutine hgrad_nodes
