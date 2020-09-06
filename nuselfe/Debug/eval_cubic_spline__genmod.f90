        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:11 2011
        MODULE EVAL_CUBIC_SPLINE__genmod
          INTERFACE 
            SUBROUTINE EVAL_CUBIC_SPLINE(NPTS,XCOR,YY,YPP,NPTS2,XOUT,   &
     &IXMIN,XMIN,XMAX,YYOUT)
              INTEGER(KIND=4), INTENT(IN) :: NPTS2
              INTEGER(KIND=4), INTENT(IN) :: NPTS
              REAL(KIND=8), INTENT(IN) :: XCOR(NPTS)
              REAL(KIND=8), INTENT(IN) :: YY(NPTS)
              REAL(KIND=8), INTENT(IN) :: YPP(NPTS)
              REAL(KIND=8), INTENT(IN) :: XOUT(NPTS2)
              INTEGER(KIND=4), INTENT(IN) :: IXMIN
              REAL(KIND=8), INTENT(IN) :: XMIN
              REAL(KIND=8), INTENT(IN) :: XMAX
              REAL(KIND=8), INTENT(OUT) :: YYOUT(NPTS2)
            END SUBROUTINE EVAL_CUBIC_SPLINE
          END INTERFACE 
        END MODULE EVAL_CUBIC_SPLINE__genmod