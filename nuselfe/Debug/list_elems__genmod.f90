        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:08 2011
        MODULE LIST_ELEMS__genmod
          INTERFACE 
            SUBROUTINE LIST_ELEMS(ELEM_NODES,NODE_NUM,NX,NY,NUM_ELEMS)
              INTEGER(KIND=4), INTENT(IN) :: NUM_ELEMS
              INTEGER(KIND=4), INTENT(IN) :: NY
              INTEGER(KIND=4), INTENT(IN) :: NX
              INTEGER(KIND=4), INTENT(OUT) :: ELEM_NODES(NUM_ELEMS,3)
              INTEGER(KIND=4), INTENT(IN) :: NODE_NUM(NX,NY)
            END SUBROUTINE LIST_ELEMS
          END INTERFACE 
        END MODULE LIST_ELEMS__genmod
