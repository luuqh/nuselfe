      subroutine jet_model(Qj0,Dj0,Tj0,Sj0,Cj0,xj0,yj0,zj0,phij0,thetaj0,grav,myrank,it)
      use elfe_glbl, only : rkind,Ntj,dtj,xMT,xMS,xMC,xMV,ofe,tsel,iupwind_t,nm,js, &
                            tnd,snd,tsd,ssd,ej,iej,Volj,Taj,Saj,Caj,rhoaj
                            !xj,yj,zj,bj,hj,Tj,Sj,Cj,Mj,dMj,Mjg,Qjg,Maj,Qaj
      implicit real(rkind)(a-h,m,o-z),integer(i,j,k,n)
      real(rkind), intent(in) :: Qj0,Dj0,Tj0,Sj0,Cj0,xj0,yj0,zj0,phij0,thetaj0,grav
      integer, intent(in) :: myrank,it
      real(rkind) :: Tak,Sak,Cak,ua,va,wa
      
      xMT = 0.
      xMS = 0.
      xMC = 0.
      xMV = 0.
      ej = 0
      iej = 0
        
!!    Numerical parameters
      alpha = 0.5;
      niter = 3;
      g = -grav;
      pi=4.*atan(1.)
      k = 1;
      
!!    Computed Jet parameters
      Aj0 = pi*Dj0*Dj0/4.;
      Vj0 = Qj0/Aj0;
      bj0 = Dj0/2.;
      Vtk = Vj0;
      cospk = cosd(phij0);
      sinpk = sind(phij0);
      costk = cosd(thetaj0);
      sintk = sind(thetaj0);
      uk = Vtk*cospk*costk;
      vk = Vtk*cospk*sintk;
      wk = Vtk*sinpk;
      alpha2 = alpha*alpha;
            
!!    Ambient environment parameters
      xk = xj0+0.5*uk*dtj;
      yk = yj0+0.5*vk*dtj;
      zk = zj0+0.5*wk*dtj;
      iej = 0
      call ambient(k,xk,yk,zk,Tak,Sak,Cak,ua,va,wa)
      rhoak = eqstate(Tak,Sak); ! ambient density
      Taj(k) = Tak
      Saj(k) = Sak
      Caj(k) = Cak
      rhoaj(k) = rhoak
      um = sqrt(ua**2+va**2+wa**2);

!!    Computed Jet parameters at discharge point
      rhoj0 = eqstate(Tj0,Sj0);
      drhoj0 = rhoak - rhoj0;
      gj0_prime = g*(rhoak - rhoj0)/rhoak;
      Mj0 = Qj0*rhoj0*dtj;
      drhok = drhoj0;
      rhok = rhoj0;
      Mk = Mj0;
      Tk = Tj0;
      Sk = Sj0;
      Ck = Cj0;
      bk = bj0;
      hk = Vtk*dtj;
      gk_prime = g*drhok/rhoak;
      cospcostk = cospk*costk;
      dcospcostk = 0.;
      dbk = 0.;
      dsk = 0.;
      
!      xj(k) = xk
!      yj(k) = yk
!      zj(k) = zk
!      bj(k) = bk
!      hj(k) = hk
!      Tj(k) = Tk
!      Sj(k) = Sk
!      Cj(k) = Ck
!      Mj(k) = Mk
!      dMj(k) = 0.
!      Mjg(k) = Mk
!      Qjg(k) = Qj0*dtj
!      Maj(k) = 0.
!      Qaj(k) = 0.
      
      i = ej(k,1)
      if (i.gt.0) then
         j = ej(k,2)
         xMT(j,i) = Mk*Tk
         xMS(j,i) = Mk*Sk
         xMC(j,i) = Mk*Ck
         xMV(j,i) = Mk
      endif
      !write(myrank+2000,101) it,i,j,Tak,Tk,Sak,Sk,xk,yk,zk
      
      do k = 2,Ntj
         
         if (max(abs(Tak-Tk),abs(Sak-Sk)).le.1.e-4) goto 997
         
         dbkp1 = dbk;
         dskp1 = dsk;
         dcospcostkp1 = dcospcostk;
       
         F12 = alpha2*(Vtk - um*cospcostk)**2/(gk_prime*bk);
         E = 1.4142*(0.057+0.55*sinpk/F12)/(1.+5.*um*cospcostk/abs(Vtk - um*cospcostk));
         dMsk = rhoak*2.*pi*bk*hk*E - um*cospcostk;
                  
         do i=1,niter-1
            dMfk = rhoak*um*(2.*bk*sqrt(1. - cospcostk**2)*dskp1 + &
                   & pi*bk*dbkp1*cospcostk + 0.5*pi*bk**2*dcospcostkp1)*dtj;
            dMk = dMsk + dMfk;
            Mkp1 = Mk + dMk;
            
            Tkp1 = (Mk*Tk + dMk*Tak)/Mkp1;
            Skp1 = (Mk*Sk + dMk*Sak)/Mkp1;
            rhokp1 = eqstate(Tkp1,Skp1);
            Ckp1 = (Mk*Ck + dMk*Cak)/Mkp1;
            gkp1_prime = (rhokp1-rhoak)/rhoak*g;
            ukp1 = (Mk*uk + dMk*ua)/Mkp1;
            vkp1 = (Mk*vk + dMk*va)/Mkp1;
            wkp1 = (Mk*wk + dMk*wa + Mkp1*gkp1_prime*dtj)/Mkp1;
            velhkp1 = sqrt(ukp1**2 + vkp1**2);
            Vtkp1 = sqrt(velhkp1**2 + wkp1**2);
            hkp1 = Vtkp1*hk/Vtk;
            bkp1 = sqrt(Mkp1/(rhokp1*pi*hkp1));
            cospkp1 = velhkp1/Vtkp1;
            costkp1 = ukp1/velhkp1;
            dskp1 = Vtkp1*dtj;
            dbkp1 = bkp1 - bk;
            dcospcostkp1 = cospkp1*costkp1 - cospcostk;
         enddo

         dMfk = rhoak*um*(2.*bk*sqrt(1. - cospcostk**2)*dskp1 + &
                & pi*bk*dbkp1*cospcostk + 0.5*pi*bk**2*dcospcostkp1)*dtj;
         dMk = dMsk + dMfk;
         Mkp1 = Mk + dMk;

         Tk = (Mk*Tk + dMk*Tak)/Mkp1;
         Sk = (Mk*Sk + dMk*Sak)/Mkp1;
         rhok = eqstate(Tk,Sk);
         Ck = (Mk*Ck + dMk*Cak)/Mkp1;
         gk_prime = (rhok-rhoak)/rhoak*g;
         uk = (Mk*uk + dMk*ua)/Mkp1;
         vk = (Mk*vk + dMk*va)/Mkp1;
         wk = (Mk*wk + dMk*wa + Mkp1*gk_prime*dtj)/Mkp1;
         velhk = sqrt(uk**2 + vk**2);
         Vtkp1 = sqrt(velhk**2 + wk**2);
         hk = Vtkp1*hk/Vtk;
         bkp1 = sqrt(Mkp1/(rhok*pi*hk));
         sinpk = wk/Vtkp1;
         cospk = velhk/Vtkp1;
         costk = uk/velhk;
         dsk = Vtkp1*dtj;
         dbk = bkp1 - bk;
         dcospcostk = cospk*costk - cospcostk;
         cospcostk = cospk*costk;
         Mk = Mkp1;
         Vtk = Vtkp1;
         bk = bkp1;
         
    997  xk = xk + uk*dtj;
         yk = yk + vk*dtj;
         zk = zk + wk*dtj;
                  
         call ambient(k,xk,yk,zk,Tak,Sak,Cak,ua,va,wa)
         rhoak = eqstate(Tak,Sak)
         Taj(k) = Tak
         Saj(k) = Sak
         Caj(k) = Cak
         rhoaj(k) = rhoak
          
!         xj(k) = xk
!         yj(k) = yk
!         zj(k) = zk
!         bj(k) = bk
!         hj(k) = hk
!         Tj(k) = Tk
!         Sj(k) = Sk
!         Cj(k) = Ck
!         Mj(k) = Mk
!         dMj(k) = dMk
!         Mjg(k) = Mk !Mjg(k-1) + dMk
!         Qjg(k) = Mk/rhoak !Qjg(k-1) + dMk/rhoak 
!         Maj(k) = (Ntj+1-k)*dMk
!         Qaj(k) = (Ntj+1-k)*dMk/rhoak
         
         i = ej(k,1)
         if (i.gt.0) then
            j = ej(k,2)
            xMT(j,i) = xMT(j,i) + Mk*Tk !- dMk*Tak
            xMS(j,i) = xMS(j,i) + Mk*Sk !- dMk*Sak
            xMC(j,i) = xMC(j,i) + Mk*Ck !- dMk*Cak
            xMV(j,i) = xMV(j,i) + Mk !- dMk
         endif
         !write(myrank+2000,101) it,i,j,Tak,Tk,Sak,Sk,xk,yk,zk
         
!         if (max(abs(Tak-Tk),abs(Sak-Sk)).le.1.e-4) then
!            Ntk = k
!            goto 998
!         endif   

      enddo !k = 2,Ntj
!  998 continue
      !write(myrank+2000,101) it,i,j,Tak,Tk,Sak,Sk,xk,yk,zk  !last
      
      do k = 1,Ntj
         i = ej(k,1) !one of 1:nofe
         if (i.gt.0) then
            j = ej(k,2)
            if (iej(j,i).eq.1) then
               iej(j,i) = 0  !(j,i) repeated, need to do once only
               xMek = Volj(k)*rhoaj(k)
               Tak = Taj(k)
               Sak = Saj(k)
               Cak = Caj(k) 
               Mnew = xMV(j,i) + xMek
               Tnew = (xMT(j,i) + xMek*Tak)/Mnew
               Snew = max(0.,(xMS(j,i) + xMek*Sak)/Mnew)
               Cnew = max(0.,(xMC(j,i) + xMek*Cak)/Mnew)
               ie=ofe(i)
               if (iupwind_t/=0) then
                  tsel(1,j,ie) = Tnew
                  tsel(2,j,ie) = Snew 
                  tsel(1,j-1,ie) = Tnew
                  tsel(2,j-1,ie) = Snew
               else
                  n1=nm(ie,1)
                  n2=nm(ie,2)
                  n3=nm(ie,3)
                  tnd(j,n1) = Tnew 
                  tnd(j,n2) = Tnew 
                  tnd(j,n3) = Tnew 
                  snd(j,n1) = Snew 
                  snd(j,n2) = Snew 
                  snd(j,n3) = Snew
                  tnd(j-1,n1) = Tnew 
                  tnd(j-1,n2) = Tnew 
                  tnd(j-1,n3) = Tnew 
                  snd(j-1,n1) = Snew 
                  snd(j-1,n2) = Snew 
                  snd(j-1,n3) = Snew 
                  n1=js(ie,1)
                  n2=js(ie,2)
                  n3=js(ie,3)
                  tsd(j,n1) = Tnew 
                  tsd(j,n2) = Tnew 
                  tsd(j,n3) = Tnew 
                  ssd(j,n1) = Snew 
                  ssd(j,n2) = Snew 
                  ssd(j,n3) = Snew 
                  tsd(j-1,n1) = Tnew 
                  tsd(j-1,n2) = Tnew 
                  tsd(j-1,n3) = Tnew 
                  ssd(j-1,n1) = Snew 
                  ssd(j-1,n2) = Snew 
                  ssd(j-1,n3) = Snew 
               endif
               write(myrank+1000,101) it,i,j,Tak,Tnew,Sak,Snew
            endif
         endif
      enddo !k = 1,Ntk
      
  101 format(3i6,4f10.4,3f16.4)
      
      return
      end subroutine

!! ------------------------------------      
      subroutine ambient(k,xk,yk,zk,Tak,Sak,Cak,ua,va,wa)
      use elfe_glbl, only : rkind,nofe,ofe,ej,iej,Volj,nvrt,tsel,iupwind_t, &
                     x,y,nm,area,dpe,ze,tnd,snd,uu2,vv2,ww2
      implicit real(rkind)(a-h,m,o-z),integer(i,j,n)
      integer, intent(in) :: k
      real(rkind), intent(in)  :: xk,yk,zk
      real(rkind), intent(out) :: Tak,Sak,Cak,ua,va,wa
      
      Volj(k) = 0.
      
      Tak = 0.
      Sak = 0.
      Cak = 0.
      ua = 0.
      va = 0.
      wa = 0.
      do jh=1,nofe
         ie=ofe(jh)
         n1=nm(ie,1)
         n2=nm(ie,2)
         n3=nm(ie,3)
         area0=2.*area(ie) 
         x1 = x(n1)
         y1 = y(n1)
         x2 = x(n2)
         y2 = y(n2)
         x3 = x(n3)
         y3 = y(n3)
         area1=abs(x1*y2+x2*yk+xk*y1-x1*yk-xk*y2-x2*y1)
         area2=abs(x1*yk+xk*y3+x3*y1-x1*y3-x3*yk-xk*y1)
         area3=abs(xk*y2+x2*y3+x3*yk-xk*y3-x3*y2-x2*yk)
         area3=area3+area2+area1 
         if (abs(area3-area0).le.1.e-4*area0) then
            ej(k,1) = jh
            dpei = dpe(ie)
            do jv=2,nvrt
               z1=dpei+ze(jv-1,ie)
               z2=dpei+ze(jv,ie)
               if ((zk.ge.z1.and.zk.le.z2).or.(jv.eq.nvrt)) then
                  ej(k,2) = jv
                  iej(jv,jh) = 1
                  Volj(k) = 0.5*area0*abs(z2-z1)
                  if (iupwind_t/=0) then
                     Tak = 0.5*(tsel(1,jv-1,ie)+tsel(1,jv,ie))
                     Sak = 0.5*(tsel(2,jv-1,ie)+tsel(2,jv,ie))
                     Cak = 0.
                  else
                     Tak = (tnd(jv-1,n1)+tnd(jv-1,n2)+tnd(jv-1,n3)+tnd(jv,n1)+tnd(jv,n2)+tnd(jv,n3))/6.
                     Sak = (snd(jv-1,n1)+snd(jv-1,n2)+snd(jv-1,n3)+snd(jv,n1)+snd(jv,n2)+snd(jv,n3))/6.
                     Cak = 0.
                  endif
                  ua  = (uu2(jv-1,n1)+uu2(jv-1,n2)+uu2(jv-1,n3)+uu2(jv,n1)+uu2(jv,n2)+uu2(jv,n3))/6.
                  va  = (vv2(jv-1,n1)+vv2(jv-1,n2)+vv2(jv-1,n3)+vv2(jv,n1)+vv2(jv,n2)+vv2(jv,n3))/6.
                  wa  = (ww2(jv-1,n1)+ww2(jv-1,n2)+ww2(jv-1,n3)+ww2(jv,n1)+ww2(jv,n2)+ww2(jv,n3))/6.
                  goto 999
               endif
            enddo   
         endif
      enddo
 999  continue
       
      !write(1,'(2i6,4f12.6)')jh,jv,dpei,zk,z1,z2
      return
      end subroutine
