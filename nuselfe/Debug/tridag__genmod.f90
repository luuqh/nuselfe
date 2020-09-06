        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:21 2011
        MODULE TRIDAG__genmod
          INTERFACE 
            SUBROUTINE TRIDAG(NMAX,NVEC,N,NC,A,B,C,R,U,GAM)
              INTEGER(KIND=4), INTENT(IN) :: NVEC
              INTEGER(KIND=4), INTENT(IN) :: NMAX
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: NC
              REAL(KIND=8), INTENT(IN) :: A(NMAX)
              REAL(KIND=8), INTENT(IN) :: B(NMAX)
              REAL(KIND=8), INTENT(IN) :: C(NMAX)
              REAL(KIND=8), INTENT(IN) :: R(NMAX,NVEC)
              REAL(KIND=8), INTENT(OUT) :: U(NMAX,NVEC)
              REAL(KIND=8), INTENT(OUT) :: GAM(NMAX)
            END SUBROUTINE TRIDAG
          END INTERFACE 
        END MODULE TRIDAG__genmod
