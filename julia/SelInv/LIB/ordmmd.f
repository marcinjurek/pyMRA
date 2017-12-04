C***********************************************************************
C***********************************************************************
C
C   Version:        0.3
C   Last modified:  December 27, 1994
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C****     ORDMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE CALLS LIU'S MULTIPLE MINIMUM DEGREE
C               ROUTINE.
C
C     INPUT PARAMETERS -
C        OUTUNT - OUTPUT UNIT.
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE.
C        IWSIZ  - SIZE OF INTEGER WORKING STORAGE.
C
C     OUTPUT PARAMETERS -
C        PERM   - THE MINIMUM DEGREE ORDERING.
C        INVP   - THE INVERSE OF PERM.
C        NNZL   - NUMBER OF NONZEROS IN LOWER TRIAGULAR FACTOR.
C        NOFSUB - NUMBER OF SUBSCRIPTS FOR THE COMPRESSED STORAGE
C                 SCHEME.
C        NOFSUB - NUMBER OF SUBSCRIPTS FOR THE COMPRESSED STORAGE
C                 SCHEME.
C        COLCNT - NUMBER OF NONZEROS IN EACH FACTOR COLUMN, INCLUDING
C                 THE DIAGONAL ENTRY (DEGREE+1).
C        NSUPER - NUMBER OF SUPERNODES.
C        XSUPER - FIRST COLUMN OF EACH SUPERNODE.
C        SNODE  - SUPERNODE MEMBERSHIP OF EACH COLUMN.
C        SFIFLG - SFIFLG=.F. MEANS SKIP SYMBOLIC FACTORIZATION 
C                 INITIALIZATION (SFINIT), SFIFLG=.T. MEANS EXECUTE
C                 SFINIT.
C        IFLAG  - ERROR FLAG.
C                   0: SUCCESSFUL ORDERING
C                  -1: INSUFFICIENT WORKING STORAGE
C                      [IWORK(*)].
C
C     WORKING PARAMETERS -
C        IWORK  - INTEGER WORKSPACE OF LENGTH 4*NEQNS.
C
C***********************************************************************
C
      SUBROUTINE ORDMMD  (  OUTUNT, NEQNS , XADJ  , ADJNCY, INVP  , 
     1                      PERM  , IWSIZ , IWORK , NNZL  , NOFSUB, 
     1                      COLCNT, NSUPER, XSUPER, SNODE , SFIFLG, 
     1                      IFLAG                                   )
C
C***********************************************************************
C
         INTEGER    ADJNCY(1), COLCNT(1), INVP(1), IWORK(1), 
     &              PERM(1), SNODE(1), XSUPER(1)
         INTEGER    XADJ(1)
         INTEGER    DELTA , IFLAG , IWSIZ , MAXINT, NEQNS , 
     &              NNZL  , NOFSUB, NSUPER, OUTUNT
         LOGICAL    SFIFLG
C
C*********************************************************************
C
        IFLAG = 0
        IF  ( IWSIZ .LT. 4*NEQNS )  THEN
            IFLAG = -1
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) '*** INTEGER WORK SPACE = ',
     &                       IWSIZ
            WRITE (OUTUNT,*) '*** IS SMALLER THAN REQUIRED = ',
     &                       4*NEQNS
            RETURN
        ENDIF
C
C       DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
C       MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER
C                (ANY SMALLER ESTIMATE WILL DO) FOR MARKING
C                NODES.
C
        DELTA  = 0
        MAXINT = 32767
        CALL GENMMD  (  NEQNS , XADJ  , ADJNCY, INVP  , PERM  ,
     1                  DELTA , 
     1                  IWORK(1)              ,
     1                  IWORK(NEQNS+1)        ,
     1                  IWORK(2*NEQNS+1)      ,
     1                  IWORK(3*NEQNS+1)      ,
     1                  MAXINT, NNZL  , NOFSUB, COLCNT, NSUPER,
     1                  XSUPER, SNODE                           )
         SFIFLG = .false.
         RETURN
C
      END
