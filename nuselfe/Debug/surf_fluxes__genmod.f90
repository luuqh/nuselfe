        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:09 2011
        MODULE SURF_FLUXES__genmod
          INTERFACE 
            SUBROUTINE SURF_FLUXES(TIME,U_AIR,V_AIR,P_AIR,T_AIR,Q_AIR,  &
     &SHORTWAVE_D,SEN_FLUX,LAT_FLUX,LONGWAVE_U,LONGWAVE_D,TAU_XZ,TAU_YZ,&
     &NWS,FLUXSU00,SRAD00)
              USE ELFE_GLBL, ONLY :                                     &
     &          RKIND,                                                  &
     &          NPA,                                                    &
     &          UU2,                                                    &
     &          VV2,                                                    &
     &          TND,                                                    &
     &          SND,                                                    &
     &          KFP,                                                    &
     &          IDRY,                                                   &
     &          NVRT,                                                   &
     &          IVCOR,                                                  &
     &          IPGL,                                                   &
     &          FDB,                                                    &
     &          LFDB
              REAL(KIND=8), INTENT(IN) :: TIME
              REAL(KIND=8), INTENT(IN) :: U_AIR(NPA)
              REAL(KIND=8), INTENT(IN) :: V_AIR(NPA)
              REAL(KIND=8), INTENT(IN) :: P_AIR(NPA)
              REAL(KIND=8), INTENT(IN) :: T_AIR(NPA)
              REAL(KIND=8), INTENT(IN) :: Q_AIR(NPA)
              REAL(KIND=8), INTENT(OUT) :: SHORTWAVE_D(NPA)
              REAL(KIND=8), INTENT(OUT) :: SEN_FLUX(NPA)
              REAL(KIND=8), INTENT(OUT) :: LAT_FLUX(NPA)
              REAL(KIND=8), INTENT(OUT) :: LONGWAVE_U(NPA)
              REAL(KIND=8), INTENT(OUT) :: LONGWAVE_D(NPA)
              REAL(KIND=8), INTENT(OUT) :: TAU_XZ(NPA)
              REAL(KIND=8), INTENT(OUT) :: TAU_YZ(NPA)
              INTEGER(KIND=4), INTENT(IN) :: NWS
              REAL(KIND=8), INTENT(IN) :: FLUXSU00
              REAL(KIND=8), INTENT(IN) :: SRAD00
            END SUBROUTINE SURF_FLUXES
          END INTERFACE 
        END MODULE SURF_FLUXES__genmod
