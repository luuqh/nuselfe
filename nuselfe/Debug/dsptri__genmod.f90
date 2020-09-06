        !COMPILER-GENERATED INTERFACE MODULE: Sun Mar 20 15:46:17 2011
        MODULE DSPTRI__genmod
          INTERFACE 
            SUBROUTINE DSPTRI(UPLO,N,AP,IPIV,WORK,INFO)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: AP(*)
              INTEGER(KIND=4) :: IPIV(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DSPTRI
          END INTERFACE 
        END MODULE DSPTRI__genmod
