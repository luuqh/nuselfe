
!------------------------------------------------------------
!comput the z coordinate and delz for each node at each level.
!  Z(nvrt,np), delz(nvrt-1,np);
!------------------------------------------------------------
SUBROUTINE update_z
USE post_vars
INTEGER             :: i;

!-------------------------------------------------------
! compute cs for each node, iformat == 5
!-------------------------------------------------------
!if(ilev.eq.-2) then

!IF(iformat /=5) THEN
!    WRITE(*,*) "ifromat must be 5 ..."
!ENDIF
iformat  = 6   ! self data format.
!-------------------------------------------------------
! compute cs for each node, iformat == 5
!-------------------------------------------------------
if(iformat.lt.5) then !ELCIRC
  do i=1,np_global
    z(0,i)=zmsl-dp(i)
	z(1,i)=ztot(1)
    delz(1,i)=z(1,i)-z(0,i)
    do k=2,nvrt
	  z(k,i)=ztot(k)
      delz(k,i)=z(k,i)-z(k-1,i)
    enddo
  enddo
else !iformat=5 SELFE

  do i=1,np_global
    hmod(i)=min(dp(i),h_s)
  enddo

  do k=1,nsig
    cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
   &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
  enddo !k=1,nvrt

endif !iformat


!-------------------------------------------------------
! compute z delz,cordinate by eta2, sigma and z levels.
!-------------------------------------------------------
do i=1,np_global

   if(kfp(i).eq.0) then !dry
     cycle
   else !wet
!         S-levels
     do k=kz,nvrt
       kin=k-kz+1
       if(hmod(i)<=h_c) then
         z(k,i)=sigma(kin)*(hmod(i)+eta2(i))+eta2(i)
       else
         z(k,i)=eta2(i)*(1+sigma(kin))+h_c*sigma(kin)+(hmod(i)-h_c)*cs(kin)
       endif
     enddo !k=kz,nvrt
!         Z-levels
     if(dp(i)<=h_s) then
       go to 88 !kbp(i)=kz
     else !bottom index
       z(kbp(i),i)=-dp(i)
      do k=kbp(i)+1,kz-1
        z(k,i)=ztot(k)
      enddo !k
     endif
!         level thicknesses
88	     do k=kbp(i)+1,nvrt
       delz(k,i)=z(k,i)-z(k-1,i)
     enddo
   endif !wet ot dry
enddo !i=1,np


END SUBROUTINE