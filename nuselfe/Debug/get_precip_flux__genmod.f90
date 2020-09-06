        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:09 2011
        MODULE GET_PRECIP_FLUX__genmod
          INTERFACE 
            SUBROUTINE GET_PRECIP_FLUX(TIME,PRECIP_FLUX)
              USE ELFE_GLBL, ONLY :                                     &
     &          RKIND,                                                  &
     &          NPA,                                                    &
     &          FDB,                                                    &
     &          LFDB
              REAL(KIND=8), INTENT(IN) :: TIME
              REAL(KIND=8), INTENT(OUT) :: PRECIP_FLUX(NPA)
            END SUBROUTINE GET_PRECIP_FLUX
          END INTERFACE 
        END MODULE GET_PRECIP_FLUX__genmod
