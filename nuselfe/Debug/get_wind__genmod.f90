        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:09 2011
        MODULE GET_WIND__genmod
          INTERFACE 
            SUBROUTINE GET_WIND(TIME,U_AIR_NODE,V_AIR_NODE,P_AIR_NODE,  &
     &T_AIR_NODE,Q_AIR_NODE)
              USE ELFE_GLBL, ONLY :                                     &
     &          RKIND,                                                  &
     &          NPA,                                                    &
     &          FDB,                                                    &
     &          LFDB
              REAL(KIND=8), INTENT(IN) :: TIME
              REAL(KIND=8), INTENT(OUT) :: U_AIR_NODE(NPA)
              REAL(KIND=8), INTENT(OUT) :: V_AIR_NODE(NPA)
              REAL(KIND=8), INTENT(OUT) :: P_AIR_NODE(NPA)
              REAL(KIND=8), INTENT(OUT) :: T_AIR_NODE(NPA)
              REAL(KIND=8), INTENT(OUT) :: Q_AIR_NODE(NPA)
            END SUBROUTINE GET_WIND
          END INTERFACE 
        END MODULE GET_WIND__genmod
