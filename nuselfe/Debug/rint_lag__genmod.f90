        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:11 2011
        MODULE RINT_LAG__genmod
          INTERFACE 
            FUNCTION RINT_LAG(MNV,NMIN,NMAX,M,K,SIGMA,SIGMAP,SIGMA_PROD,&
     &PSI,GAM,COEF) RESULT(RINT_LAG_0)
              INTEGER(KIND=4), INTENT(IN) :: MNV
              INTEGER(KIND=4), INTENT(IN) :: NMIN
              INTEGER(KIND=4), INTENT(IN) :: NMAX
              INTEGER(KIND=4), INTENT(IN) :: M
              INTEGER(KIND=4), INTENT(IN) :: K
              REAL(KIND=8), INTENT(IN) :: SIGMA(MNV)
              REAL(KIND=8), INTENT(IN) :: SIGMAP(MNV,10)
              REAL(KIND=8), INTENT(IN) :: SIGMA_PROD(MNV,MNV,-4:4)
              REAL(KIND=8), INTENT(IN) :: PSI(MNV)
              REAL(KIND=8), INTENT(OUT) :: GAM(MNV)
              REAL(KIND=8), INTENT(OUT) :: COEF(0:MNV)
              REAL(KIND=8) :: RINT_LAG_0
            END FUNCTION RINT_LAG
          END INTERFACE 
        END MODULE RINT_LAG__genmod
