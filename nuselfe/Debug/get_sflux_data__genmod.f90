        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:08 2011
        MODULE GET_SFLUX_DATA__genmod
          INTERFACE 
            SUBROUTINE GET_SFLUX_DATA(TIME_NOW,INFO,DATA_NAME,DATA_OUT, &
     &GOT_SUITABLE_BRACKET,NUM_NODES_OUT)
              USE NETCDF_IO
              INTEGER(KIND=4), INTENT(IN) :: NUM_NODES_OUT
              TYPE (DATASET_INFO), INTENT(IN) :: INFO
              REAL(KIND=8), INTENT(IN) :: TIME_NOW
              CHARACTER(*), INTENT(IN) :: DATA_NAME
              REAL(KIND=8), INTENT(OUT) :: DATA_OUT(NUM_NODES_OUT)
              LOGICAL(KIND=4), INTENT(OUT) :: GOT_SUITABLE_BRACKET
            END SUBROUTINE GET_SFLUX_DATA
          END INTERFACE 
        END MODULE GET_SFLUX_DATA__genmod