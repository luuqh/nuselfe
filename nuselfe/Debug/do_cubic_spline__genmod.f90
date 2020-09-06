        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:11 2011
        MODULE DO_CUBIC_SPLINE__genmod
          INTERFACE 
            SUBROUTINE DO_CUBIC_SPLINE(NPTS,XCOR,YY,YP1,YP2,NPTS2,XOUT, &
     &IXMIN,XMIN,XMAX,YYOUT)
              INTEGER(KIND=4), INTENT(IN) :: NPTS2
              INTEGER(KIND=4), INTENT(IN) :: NPTS
              REAL(KIND=8), INTENT(IN) :: XCOR(NPTS)
              REAL(KIND=8), INTENT(IN) :: YY(NPTS)
              REAL(KIND=8), INTENT(IN) :: YP1
              REAL(KIND=8), INTENT(IN) :: YP2
              REAL(KIND=8), INTENT(IN) :: XOUT(NPTS2)
              INTEGER(KIND=4), INTENT(IN) :: IXMIN
              REAL(KIND=8), INTENT(IN) :: XMIN
              REAL(KIND=8), INTENT(IN) :: XMAX
              REAL(KIND=8), INTENT(OUT) :: YYOUT(NPTS2)
            END SUBROUTINE DO_CUBIC_SPLINE
          END INTERFACE 
        END MODULE DO_CUBIC_SPLINE__genmod
