        !COMPILER-GENERATED INTERFACE MODULE: Sun Mar 20 15:46:17 2011
        MODULE DSPR__genmod
          INTERFACE 
            SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: AP(*)
            END SUBROUTINE DSPR
          END INTERFACE 
        END MODULE DSPR__genmod
