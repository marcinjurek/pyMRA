C***********************************************************************
C***********************************************************************
C
C   Version:        0.3
C   Last modified:  December 27, 1994
C   Authors:        Joseph W.H. Liu
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = GENMMD
C  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
C***********************************************************************
C***********************************************************************
C****     GENMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE IMPLEMENTS THE MINIMUM DEGREE
C        ALGORITHM.  IT MAKES USE OF THE IMPLICIT REPRESENTATION
C        OF ELIMINATION GRAPHS BY QUOTIENT GRAPHS, AND THE
C        NOTION OF INDISTINGUISHABLE NODES.  IT ALSO IMPLEMENTS
C        THE MODIFICATIONS BY MULTIPLE ELIMINATION AND MINIMUM
C        EXTERNAL DEGREE.
C        ---------------------------------------------
C        CAUTION - THE ADJACENCY VECTOR ADJNCY WILL BE
C        DESTROYED.
C        ---------------------------------------------
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE.
C        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
C        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER
C                 (ANY SMALLER ESTIMATE WILL DO) FOR MARKING
C                 NODES.
C
C     OUTPUT PARAMETERS -
C        PERM   - THE MINIMUM DEGREE ORDERING.
C        INVP   - THE INVERSE OF PERM.
C        NNZL   - NUMBER OF NONZEROS IN LOWER TRIAGULAR FACTOR.
C        NOFSUB - NUMBER OF SUBSCRIPTS FOR THE COMPRESSED STORAGE 
C                 SCHEME.
C        COLCNT - NUMBER OF NONZEROS IN EACH FACTOR COLUMN, INCLUDING
C                 THE DIAGONAL ENTRY (DEGREE+1).
C        NSUPER - NUMBER OF SUPERNODES.
C        XSUPER - FIRST COLUMN OF EACH SUPERNODE.
C        SNODE  - SUPERNODE MEMBERSHIP OF EACH COLUMN.
C
C     WORKING PARAMETERS -
C        DHEAD  - VECTOR FOR HEAD OF DEGREE LISTS.
C        INVP   - USED TEMPORARILY FOR DEGREE FORWARD LINK.
C        PERM   - USED TEMPORARILY FOR DEGREE BACKWARD LINK.
C        QSIZE  - VECTOR FOR SIZE OF SUPERNODES.
C        LLIST  - VECTOR FOR TEMPORARY LINKED LISTS.
C        MARKER - A TEMPORARY MARKER VECTOR.
C
C     PROGRAM SUBROUTINES -
C        MMDELM, MMDINT, MMDNUM, MMDUPD.
C
C***********************************************************************
C
      SUBROUTINE  GENMMD ( NEQNS, XADJ, ADJNCY, INVP, PERM,
     1                     DELTA, DHEAD, QSIZE, LLIST, MARKER,
     1                     MAXINT, NNZL, NOFSUB, COLCNT, NSUPER, 
     1                     XSUPER, SNODE )
C
C***********************************************************************
C
         INTEGER    ADJNCY(1), DHEAD(1) , INVP(1)  , LLIST(1) ,
     1              MARKER(1), PERM(1)  , QSIZE(1) , COLCNT(1),
     1              XSUPER(1), SNODE(1)
         INTEGER    XADJ(1)
         INTEGER    CC    , DELTA , EHEAD , FSTCOL, I     , 
     1              KCOL  , KSUPER, LSTCOL, MAXINT, MDEG  ,
     1              MDLMT , MDNODE, NEQNS , NEXTMD, NNZL  ,
     1              NOFSUB, NSUPER, NUM   , TAG
C
C***********************************************************************
C
         IF  ( NEQNS .LE. 0 )  RETURN
C
C        ------------------------------------------------
C        INITIALIZATION FOR THE MINIMUM DEGREE ALGORITHM.
C        ------------------------------------------------
         NSUPER = 0
         CALL  MMDINT ( NEQNS, XADJ, ADJNCY, DHEAD, INVP, PERM,
     1                  QSIZE, LLIST, MARKER )
C
C        ----------------------------------------------
C        NUM COUNTS THE NUMBER OF ORDERED NODES PLUS 1.
C        ----------------------------------------------
         NUM = 1
C
C        -----------------------------
C        ELIMINATE ALL ISOLATED NODES.
C        -----------------------------
         NEXTMD = DHEAD(1)
  100    CONTINUE
             IF  ( NEXTMD .LE. 0 )  GO TO 200
                 MDNODE = NEXTMD
                 NEXTMD = INVP(MDNODE)
                 MARKER(MDNODE) = MAXINT
                 INVP(MDNODE) = - NUM
                 NSUPER = NSUPER + 1
                 XSUPER(NSUPER) = NUM
                 COLCNT(NUM) = 1
                 NUM = NUM + 1
                 GO TO 100
C
  200    CONTINUE
C        ----------------------------------------
C        SEARCH FOR NODE OF THE MINIMUM DEGREE.
C        MDEG IS THE CURRENT MINIMUM DEGREE;
C        TAG IS USED TO FACILITATE MARKING NODES.
C        ----------------------------------------
         IF  ( NUM .GT. NEQNS )  GO TO 1000
         TAG = 1
         DHEAD(1) = 0
         MDEG = 2
  300    CONTINUE
             IF  ( DHEAD(MDEG) .GT. 0 )  GO TO 400
                 MDEG = MDEG + 1
                 GO TO 300
  400        CONTINUE
C            -------------------------------------------------
C            USE VALUE OF DELTA TO SET UP MDLMT, WHICH GOVERNS
C            WHEN A DEGREE UPDATE IS TO BE PERFORMED.
C            -------------------------------------------------
             MDLMT = MDEG + DELTA
             EHEAD = 0
C
  500        CONTINUE
                 MDNODE = DHEAD(MDEG)
                 IF  ( MDNODE .GT. 0 )  GO TO 600
                     MDEG = MDEG + 1
                     IF  ( MDEG .GT. MDLMT )  GO TO 900
                         GO TO 500
  600            CONTINUE
C                ----------------------------------------
C                REMOVE MDNODE FROM THE DEGREE STRUCTURE.
C                ----------------------------------------
                 NEXTMD = INVP(MDNODE)
                 DHEAD(MDEG) = NEXTMD
                 IF  ( NEXTMD .GT. 0 )  PERM(NEXTMD) = - MDEG
                 NSUPER = NSUPER + 1
                 XSUPER(NSUPER) = NUM
                 COLCNT(NUM) = MDEG + QSIZE(MDNODE) - 1
                 INVP(MDNODE) = - NUM
                 IF  ( NUM+QSIZE(MDNODE) .GT. NEQNS )  GO TO 1000
C                ----------------------------------------------
C                ELIMINATE MDNODE AND PERFORM QUOTIENT GRAPH
C                TRANSFORMATION.  RESET TAG VALUE IF NECESSARY.
C                ----------------------------------------------
                 TAG = TAG + 1
                 IF  ( TAG .LT. MAXINT )  GO TO 800
                     TAG = 1
                     DO  700  I = 1, NEQNS
                         IF  ( MARKER(I) .LT. MAXINT )  MARKER(I) = 0
  700                CONTINUE
  800            CONTINUE
                 CALL  MMDELM ( MDNODE, XADJ, ADJNCY, DHEAD, INVP,
     1                          PERM, QSIZE, LLIST, MARKER, MAXINT,
     1                          TAG )
                 NUM = NUM + QSIZE(MDNODE)
                 LLIST(MDNODE) = EHEAD
                 EHEAD = MDNODE
                 IF  ( DELTA .GE. 0 )  GO TO 500
  900        CONTINUE
C            -------------------------------------------
C            UPDATE DEGREES OF THE NODES INVOLVED IN THE
C            MINIMUM DEGREE NODES ELIMINATION.
C            -------------------------------------------
             IF  ( NUM .GT. NEQNS )  GO TO 1000
             CALL  MMDUPD ( EHEAD, NEQNS, XADJ, ADJNCY, DELTA, MDEG,
     1                      DHEAD, INVP, PERM, QSIZE, LLIST, MARKER,
     1                      MAXINT, TAG )
             GO TO 300
C
 1000    CONTINUE
         XSUPER(NSUPER+1) = NEQNS + 1
         NNZL = 0
         NOFSUB = 0
         DO  1200  KSUPER = 1, NSUPER
             FSTCOL = XSUPER(KSUPER)
             LSTCOL = XSUPER(KSUPER+1) - 1
             CC = COLCNT(FSTCOL)
             NOFSUB = NOFSUB + CC
             DO  1100  KCOL = FSTCOL, LSTCOL
                 SNODE(KCOL) = KSUPER
                 COLCNT(KCOL) = CC
                 NNZL = NNZL + CC
                 CC = CC - 1
 1100        CONTINUE
 1200    CONTINUE
         CALL  MMDNUM ( NEQNS, PERM, INVP, QSIZE )
         RETURN
C
      END
