        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:10 2011
        MODULE VINTER__genmod
          INTERFACE 
            SUBROUTINE VINTER(NMAX1,NMAX2,NC,ZT,K1,K2,K3,ZA,SINT,SOUT,  &
     &IBELOW)
              INTEGER(KIND=4), INTENT(IN) :: NMAX2
              INTEGER(KIND=4), INTENT(IN) :: NMAX1
              INTEGER(KIND=4), INTENT(IN) :: NC
              REAL(KIND=8), INTENT(IN) :: ZT
              INTEGER(KIND=4), INTENT(IN) :: K1
              INTEGER(KIND=4), INTENT(IN) :: K2
              INTEGER(KIND=4), INTENT(IN) :: K3
              REAL(KIND=8), INTENT(IN) :: ZA(NMAX1)
              REAL(KIND=8), INTENT(IN) :: SINT(NMAX1,NMAX2)
              REAL(KIND=8), INTENT(OUT) :: SOUT(NMAX2)
              INTEGER(KIND=4), INTENT(OUT) :: IBELOW
            END SUBROUTINE VINTER
          END INTERFACE 
        END MODULE VINTER__genmod
