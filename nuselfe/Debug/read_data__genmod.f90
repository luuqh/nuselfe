        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:08 2011
        MODULE READ_DATA__genmod
          INTERFACE 
            SUBROUTINE READ_DATA(FILE_NAME,DATA_NAME,DATA,NX,NY,TIME_NUM&
     &)
              INTEGER(KIND=4), INTENT(IN) :: NY
              INTEGER(KIND=4), INTENT(IN) :: NX
              CHARACTER(*), INTENT(IN) :: FILE_NAME
              CHARACTER(*), INTENT(IN) :: DATA_NAME
              REAL(KIND=8), INTENT(OUT) :: DATA(NX,NY)
              INTEGER(KIND=4), INTENT(IN) :: TIME_NUM
            END SUBROUTINE READ_DATA
          END INTERFACE 
        END MODULE READ_DATA__genmod