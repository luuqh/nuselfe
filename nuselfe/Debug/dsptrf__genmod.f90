        !COMPILER-GENERATED INTERFACE MODULE: Sun Mar 20 15:46:17 2011
        MODULE DSPTRF__genmod
          INTERFACE 
            SUBROUTINE DSPTRF(UPLO,N,AP,IPIV,INFO)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: AP(*)
              INTEGER(KIND=4) :: IPIV(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DSPTRF
          END INTERFACE 
        END MODULE DSPTRF__genmod
