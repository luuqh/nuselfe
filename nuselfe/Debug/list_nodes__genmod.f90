        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 27 11:48:08 2011
        MODULE LIST_NODES__genmod
          INTERFACE 
            SUBROUTINE LIST_NODES(NODE_I,NODE_J,NODE_NUM,NUM_NODES,NX,NY&
     &)
              INTEGER(KIND=4), INTENT(IN) :: NY
              INTEGER(KIND=4), INTENT(IN) :: NX
              INTEGER(KIND=4), INTENT(IN) :: NUM_NODES
              INTEGER(KIND=4), INTENT(OUT) :: NODE_I(NUM_NODES)
              INTEGER(KIND=4), INTENT(OUT) :: NODE_J(NUM_NODES)
              INTEGER(KIND=4), INTENT(OUT) :: NODE_NUM(NX,NY)
            END SUBROUTINE LIST_NODES
          END INTERFACE 
        END MODULE LIST_NODES__genmod