        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:20 2011
        MODULE AREA_COORD__genmod
          INTERFACE 
            SUBROUTINE AREA_COORD(NNEL,XT,YT,ARCO)
              INTEGER(KIND=4), INTENT(IN) :: NNEL
              REAL(KIND=8), INTENT(IN) :: XT
              REAL(KIND=8), INTENT(IN) :: YT
              REAL(KIND=8), INTENT(OUT) :: ARCO(3)
            END SUBROUTINE AREA_COORD
          END INTERFACE 
        END MODULE AREA_COORD__genmod
