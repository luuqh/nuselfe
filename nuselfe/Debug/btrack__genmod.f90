        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:20 2011
        MODULE BTRACK__genmod
          INTERFACE 
            SUBROUTINE BTRACK(L_NS,IPSGB,IFL_BND,J0,IADVF,NDELT,DTB_MAX,&
     &VIS_COE,TIME_RM,UUINT,VVINT,WWINT,NNEL,JLEV,XT,YT,ZT,TTINT,SSINT, &
     &IEXIT)
              USE ELFE_GLBL
              INTEGER(KIND=4), INTENT(IN) :: L_NS
              INTEGER(KIND=4), INTENT(IN) :: IPSGB
              INTEGER(KIND=4), INTENT(IN) :: IFL_BND
              INTEGER(KIND=4), INTENT(IN) :: J0
              INTEGER(KIND=4), INTENT(IN) :: IADVF
              INTEGER(KIND=4), INTENT(IN) :: NDELT
              REAL(KIND=8), INTENT(IN) :: DTB_MAX
              REAL(KIND=8), INTENT(IN) :: VIS_COE
              REAL(KIND=8), INTENT(INOUT) :: TIME_RM
              REAL(KIND=8), INTENT(INOUT) :: UUINT
              REAL(KIND=8), INTENT(INOUT) :: VVINT
              REAL(KIND=8), INTENT(INOUT) :: WWINT
              INTEGER(KIND=4), INTENT(INOUT) :: NNEL
              INTEGER(KIND=4), INTENT(INOUT) :: JLEV
              REAL(KIND=8), INTENT(INOUT) :: XT
              REAL(KIND=8), INTENT(INOUT) :: YT
              REAL(KIND=8), INTENT(INOUT) :: ZT
              REAL(KIND=8), INTENT(OUT) :: TTINT
              REAL(KIND=8), INTENT(OUT) :: SSINT
              LOGICAL(KIND=4), INTENT(OUT) :: IEXIT
            END SUBROUTINE BTRACK
          END INTERFACE 
        END MODULE BTRACK__genmod
