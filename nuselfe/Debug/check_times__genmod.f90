        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:08 2011
        MODULE CHECK_TIMES__genmod
          INTERFACE 
            SUBROUTINE CHECK_TIMES(TEST_TIME,TIMES,NUM_TIMES,REPEAT_NUM,&
     &REPEAT,AT_END)
              INTEGER(KIND=4), INTENT(IN) :: NUM_TIMES
              REAL(KIND=8), INTENT(IN) :: TEST_TIME
              REAL(KIND=8), INTENT(IN) :: TIMES(NUM_TIMES)
              INTEGER(KIND=4), INTENT(OUT) :: REPEAT_NUM
              LOGICAL(KIND=4), INTENT(OUT) :: REPEAT
              LOGICAL(KIND=4), INTENT(OUT) :: AT_END
            END SUBROUTINE CHECK_TIMES
          END INTERFACE 
        END MODULE CHECK_TIMES__genmod
