        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:21 2011
        MODULE SOLVE_JCG__genmod
          INTERFACE 
            SUBROUTINE SOLVE_JCG(ITIME,MOITN,MXITN,RTOL,S,X,B,BC,LBC)
              USE ELFE_GLBL, ONLY :                                     &
     &          RKIND,                                                  &
     &          NP,                                                     &
     &          NPA,                                                    &
     &          WTIMER,                                                 &
     &          IPLG,                                                   &
     &          IPGL,                                                   &
     &          MNEI,                                                   &
     &          NNP,                                                    &
     &          INP,                                                    &
     &          ERRMSG
              INTEGER(KIND=4), INTENT(IN) :: ITIME
              INTEGER(KIND=4), INTENT(IN) :: MOITN
              INTEGER(KIND=4), INTENT(IN) :: MXITN
              REAL(KIND=8), INTENT(IN) :: RTOL
              REAL(KIND=8), INTENT(IN) :: S(NP,0:(MNEI+1))
              REAL(KIND=8), INTENT(INOUT) :: X(NPA)
              REAL(KIND=8), INTENT(IN) :: B(NP)
              REAL(KIND=8), INTENT(IN) :: BC(NPA)
              LOGICAL(KIND=4), INTENT(IN) :: LBC(NPA)
            END SUBROUTINE SOLVE_JCG
          END INTERFACE 
        END MODULE SOLVE_JCG__genmod
