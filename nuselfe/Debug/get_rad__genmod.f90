        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:09 2011
        MODULE GET_RAD__genmod
          INTERFACE 
            SUBROUTINE GET_RAD(TIME,SHORTWAVE_D,LONGWAVE_D)
              USE ELFE_GLBL, ONLY :                                     &
     &          RKIND,                                                  &
     &          NPA,                                                    &
     &          FDB,                                                    &
     &          LFDB,                                                   &
     &          ALBEDO
              REAL(KIND=8), INTENT(IN) :: TIME
              REAL(KIND=8), INTENT(OUT) :: SHORTWAVE_D(NPA)
              REAL(KIND=8), INTENT(OUT) :: LONGWAVE_D(NPA)
            END SUBROUTINE GET_RAD
          END INTERFACE 
        END MODULE GET_RAD__genmod
