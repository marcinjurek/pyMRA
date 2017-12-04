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
C******     PCHOL .... DENSE PARTIAL CHOLESKY             **************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE PERFORMS CHOLESKY
C               FACTORIZATION ON THE COLUMNS OF A SUPERNODE
C               THAT HAVE RECEIVED ALL UPDATES FROM COLUMNS
C               EXTERNAL TO THE SUPERNODE.
C
C     INPUT PARAMETERS -
C        M      - NUMBER OF ROWS (LENGTH OF THE FIRST COLUMN).
C        N      - NUMBER OF COLUMNS IN THE SUPERNODE.
C        XPNT   - XPNT(J+1) POINTS ONE LOCATION BEYOND THE END
C                 OF THE J-TH COLUMN OF THE SUPERNODE.
C        X(*)   - CONTAINS THE COLUMNS OF OF THE SUPERNODE TO
C                 BE FACTORED.
C        SMXPY  - EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY.
C
C     OUTPUT PARAMETERS -
C        X(*)   - ON OUTPUT, CONTAINS THE FACTORED COLUMNS OF
C                 THE SUPERNODE.
C        IFLAG  - UNCHANGED IF THERE IS NO ERROR.
C                 =1 IF NONPOSITIVE DIAGONAL ENTRY IS ENCOUNTERED.
C
C***********************************************************************
C
      SUBROUTINE  PCHOL  ( M, N, XPNT, X, IFLAG, SMXPY )
C
C***********************************************************************
C
C     -----------
C     PARAMETERS.
C     -----------
C
      EXTERNAL            SMXPY
C
      INTEGER             M, N, IFLAG
C
      INTEGER             XPNT(*)
C
      DOUBLE PRECISION    X(*)
C
C     ----------------
C     LOCAL VARIABLES.
C     ----------------
C
      INTEGER             JPNT  , JCOL  , MM
C
      DOUBLE PRECISION    DIAG
C
C***********************************************************************
C
C       ------------------------------------------
C       FOR EVERY COLUMN JCOL IN THE SUPERNODE ...
C       ------------------------------------------
        MM     = M
        JPNT = XPNT(1)
        DO  100  JCOL = 1, N
C
C           ----------------------------------
C           UPDATE JCOL WITH PREVIOUS COLUMNS.
C           ----------------------------------
            IF  ( JCOL .GT. 1 )  THEN
                CALL SMXPY ( MM, JCOL-1, X(JPNT), XPNT, X )
            ENDIF
C
C           ---------------------------
C           COMPUTE THE DIAGONAL ENTRY.
C           ---------------------------
            DIAG = X(JPNT)
            IF  ( DIAG  .EQ. 0.0D+00 )  THEN
                IFLAG = 1
                RETURN
            ENDIF
            DIAG = 1.0D+00 / DIAG
C
C           ----------------------------------------------------
C           SCALE COLUMN JCOL WITH RECIPROCAL OF DIAGONAL ENTRY.
C           ----------------------------------------------------
            MM = MM - 1
            JPNT = JPNT + 1
            CALL DSCAL ( MM, DIAG, X(JPNT) )
            JPNT = JPNT + MM
C
  100   CONTINUE
C
      RETURN
      END
