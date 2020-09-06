        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:20 2011
        MODULE QUICKSEARCH__genmod
          INTERFACE 
            SUBROUTINE QUICKSEARCH(TIME,X0,Y0,Z0,NNEL,JLEV,XT,YT,ZT,TRM,&
     &NFL,KBPL,ARCO,ZRAT,ZTMP,SS)
              USE ELFE_GLBL
              REAL(KIND=8), INTENT(IN) :: TIME
              REAL(KIND=8), INTENT(IN) :: X0
              REAL(KIND=8), INTENT(IN) :: Y0
              REAL(KIND=8), INTENT(IN) :: Z0
              INTEGER(KIND=4), INTENT(INOUT) :: NNEL
              INTEGER(KIND=4), INTENT(INOUT) :: JLEV
              REAL(KIND=8), INTENT(INOUT) :: XT
              REAL(KIND=8), INTENT(INOUT) :: YT
              REAL(KIND=8), INTENT(INOUT) :: ZT
              REAL(KIND=8), INTENT(OUT) :: TRM
              INTEGER(KIND=4), INTENT(OUT) :: NFL
              INTEGER(KIND=4), INTENT(OUT) :: KBPL
              REAL(KIND=8), INTENT(OUT) :: ARCO(3)
              REAL(KIND=8), INTENT(OUT) :: ZRAT
              REAL(KIND=8), INTENT(OUT) :: ZTMP(NVRT)
              REAL(KIND=8), INTENT(OUT) :: SS
            END SUBROUTINE QUICKSEARCH
          END INTERFACE 
        END MODULE QUICKSEARCH__genmod
