    
    SUBROUTINE read_eco_input

        
!!======================================================================
!! June, 2007                                                          ! 
!!======================================================Marta Rodrigues=
!!                                                                     !
!! This subroutine reads ecological models inputs                      !
!!                                                                     ! 
!!======================================================================	 

!
      USE bio_param
!      USE global   
      USE biology 

!
      IMPLICIT NONE
      SAVE  

! Local variables

      CHARACTER(len=100) :: var1, var2
      INTEGER :: i, j, ierror      

      
      OPEN(5, file='ecosim.in', status='old')
!      IF(ierror/=0) THEN
!        WRITE(*,*) 'Error opening ecological model input parameters (ecosim.in)!'
!        STOP
!      END IF
      
! Reads input

      DO j=1,18
        READ(5,*) var1
      END DO	
      
     READ(5,*) var1, var2, BioIter
      
      DO j=1,3
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, RtUVR_flag
      
      DO j=1,3
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, NFIX_flag
      
      DO j=1,3
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, Regen_flag
      
      DO j=1,8
        READ(5,*) var1
      END DO
     
!-----------------------------------------------------------------------------
! Phytoplankton group parameters.
!-----------------------------------------------------------------------------
!

     READ(5,*) var1, var2, (HsNO3(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (HsNH4(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (HsSiO(i), i=1,Nphy)  
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (HsPO4(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (HsFe(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (GtALG_max(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (PhyTbase(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (PhyTfac(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (BET_(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (maxC2nALG(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (minC2nALG(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (C2nALGminABS(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (maxC2SiALG(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (minC2SiALG(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (C2SiALGminABS(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (maxC2pALG(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (minC2pALG(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (C2pALGminABS(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (maxC2FeALG(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (minC2FeALG(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (C2FeALGminABS(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (qu_yld(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (E0_comp(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (E0_inhib(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (inhib_fac(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (C2CHL_max(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (mxC2Cl(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (b_C2Cl(i), i=1,Nphy)
      
      DO j=1,6
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (mxC2Cn(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (b_C2Cn(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (mxPacEff(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (b_PacEff(i), i=1,Nphy)
      
      DO j=1,6
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (mxChlB(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (b_ChlB(i), i=1,Nphy)
      
      DO j=1,6
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (mxChlC(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (b_ChlC(i), i=1,Nphy)
      
      DO j=1,6
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (mxPSC(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (b_PSC(i), i=1,Nphy)
      
      DO j=1,6
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (mxPPC(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (b_PPC(i), i=1,Nphy)
       
      DO j=1,6
        READ(5,*) var1
      END DO
            
     READ(5,*) var1, var2, (mxLPUb(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (b_LPUb(i), i=1,Nphy)
      
      DO j=1,6
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (mxHPUb(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (b_HPUb(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (FecDOC(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, ((FecPEL(i,j), i=1,Nphy), j=1,Nfec)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (FecCYC(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (ExALG(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (WS(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (HsGRZ(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (MinRefuge(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (RefugeDep(i), i=1,Nphy)
       
      DO j=1,4
        READ(5,*) var1
      END DO
           
     READ(5,*) var1, var2, (Norm_Vol(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (Norm_Surf(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (HsDOP(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (C2pALKPHOS(i), i=1,Nphy)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (HsDON(i), i=1,Nphy)
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (C2nNupDON(i), i=1,Nphy)
      
      DO j=1,8
        READ(5,*) var1
      END DO
      
!-----------------------------------------------------------------------------
! Bacteria group parameters.
!-----------------------------------------------------------------------------
!
     
      READ(5,*) var1, var2, (HsDOC_ba(i), i=1,Nbac)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (GtBAC_max(i), i=1,Nbac)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (BacTbase(i), i=1,Nbac)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (BacTfac(i), i=1,Nbac) 
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, C2nBAC
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, C2pBAC
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, C2FeBAC
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, BacDOC
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, BacPEL
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, BacCYC
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, ExBAC_c
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, ExBacC2N
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, Bac_Ceff
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, RtNIT
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, HsNIT
      
      DO j=1,8
        READ(5,*) var1
      END DO
      
!-----------------------------------------------------------------------------
! DOM group parameters.
!-----------------------------------------------------------------------------
!
  
     READ(5,*) var1, var2, (cDOCfrac_c(i), i=1,Ndom)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, RtUVR_DIC
      
      DO j=1,5
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, RtUVR_DOC
      
      DO j=1,8
        READ(5,*) var1
      END DO
      
!-----------------------------------------------------------------------------
! Fecal and detritus group parameters.
!-----------------------------------------------------------------------------
!
     READ(5,*) var1, var2, (WF(i), i=1,Nfec)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (RegTbase(i), i=1,Nfec)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (RegTfac(i), i=1,Nfec)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (RegCR(i), i=1,Nfec)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (RegNR(i), i=1,Nfec)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (RegSR(i), i=1,Nfec)
        
      DO j=1,4
        READ(5,*) var1
      END DO
             
     READ(5,*) var1, var2, (RegPR(i), i=1,Nfec)
      
      DO j=1,4
        READ(5,*) var1
      END DO
      
     READ(5,*) var1, var2, (RegFR(i), i=1,Nfec)
      
      DO j=1,10
        READ(5,*) var1
      END DO
         
!-----------------------------------------------------------------------------
! Zooplankton group parameters 
!-----------------------------------------------------------------------------
!          
     READ(5,*) var1, var2, (zoo_sp(i), i=1,Nzoo)
       
      DO j=1,4
        READ(5,*) var1
      END DO
          
     READ(5,*) var1, var2, (ZooDOC(i), i=1,Nzoo)
     
      DO j=1,4
        READ(5,*) var1
      END DO
           
     READ(5,*) var1, var2, ((ZooPEL(i,j), i=1,Nzoo), j=1,Nfec)
      
      DO j=1,4
        READ(5,*) var1
      END DO
           
     READ(5,*) var1, var2, (ZooCYC(i), i=1,Nzoo)
      
      DO j=1,4
        READ(5,*) var1
      END DO
           
     READ(5,*) var1, var2, ((DeltaZoo(i,j), j=1,Nphy), i=1,Nzoo)
      
      DO j=1,4
        READ(5,*) var1
      END DO
           
     READ(5,*) var1, var2, ((EfcCap(i,j), j=1,Nphy), i=1,Nzoo)
      
      DO j=1,4
        READ(5,*) var1
      END DO
           
     READ(5,*) var1, var2, (HsZoo(i), i=1,Nzoo)
      
      DO j=1,4
        READ(5,*) var1
      END DO

     READ(5,*) var1, var2, (EfcPrd(i), i=1,Nzoo)
      
      DO j=1,3
        READ(5,*) var1
      END DO            

     READ(5,*) var1, var2, (ExZoo(i), i=1,Nzoo)
      
      DO j=1,4
        READ(5,*) var1
      END DO             
     
     READ(5,*) var1, var2, (GZ(i), i=1,Nzoo)
           
!-----------------------------------------------------------------------------
! Physical and output parameters 
!-----------------------------------------------------------------------------

!     READ(5,*) var1, var2, (TNU2(i), i=1,NBIT)

!     READ(5,*) var1, var2, (TNU4(i), i=1,NBIT)

!     READ(5,*) var1, var2, (AKT_BAK(i), i=1,NBIT)

!     READ(5,*) var1, var2, (TNUDG(i), i=1,NBIT)

!     READ(5,*) var1, var2, (Hout(idTvar)(i), i=1,NBIT)
         
      CLOSE(5)     
 
      RETURN
      END SUBROUTINE read_eco_input


      
      
      
