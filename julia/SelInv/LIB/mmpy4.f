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
C*************     MMPY4  .... MATRIX-MATRIX MULTIPLY     **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE -
C       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA,
C       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY
C       CODES.
C
C       LOOP UNROLLING: LEVEL 4
C
C   INPUT PARAMETERS -
C       M               -   NUMBER OF ROWS IN X AND IN Y.
C       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS
C                           IN A.
C       Q               -   NUMBER OF COLUMNS IN A AND Y.
C       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE
C                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO
C                           USED TO ACCESS THE ROWS OF A.
C       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A.
C       LDY             -   LENGTH OF FIRST COLUMN OF Y.
C
C   UPDATED PARAMETERS -
C       Y(*)            -   ON OUTPUT, Y = Y + AX.
C
C***********************************************************************
C
      SUBROUTINE  MMPY4  (  M     , N     , Q     , XPNT  , X     ,
     &                      Y     , LDY                             )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
C
        INTEGER             LDY   , M     , N     , Q
        INTEGER             XPNT(*)
        DOUBLE PRECISION    X(*)          , Y(*)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
C
        INTEGER             I1    , I2    , I3    , I4
        INTEGER             J1    , J2    , J3    , J4
        INTEGER             IY    , IYLAST, IYSTRT, IYSTOP, LENY  ,
     &                      MM    , REMAIN, XCOL  , YCOL
        DOUBLE PRECISION    A1    , A2    , A3    , A4
C
        INTEGER             LEVEL
        PARAMETER       (   LEVEL = 4 )
C
C***********************************************************************
C
C       -----------------------------------------------------------
C       INITIAL OFFSETS, COLUMN LENGTHS, AND INDEX RANGE VARIABLES.
C       -----------------------------------------------------------
        REMAIN = MOD ( N, LEVEL ) + 1
        MM = M
        IYLAST = 0
        LENY = LDY
C
C       ------------------------------------
C       TO COMPUTE EACH COLUMN YCOL OF Y ...
C       ------------------------------------
C
        DO  700  YCOL = 1, Q
C
            IYSTRT = IYLAST + 1
            IYSTOP = IYSTRT + MM - 1
            IYLAST = IYLAST + LENY
C
C           --------------------------------------------------
C           ... PERFORM THE APPROPRATE MATRIX VECTOR MULTIPLY:
C               X * A(*,YCOL) WITH LEVEL 4 LOOP-UNROLLING.
C           --------------------------------------------------
C
            GO TO ( 400, 100, 200, 300 ), REMAIN
C
  100           CONTINUE
                J1 = XPNT(1)
                I1 = XPNT(1+1) - MM
                A1 = - X(J1)*X(I1)
                DO  150  IY = IYSTRT, IYSTOP
                    Y(IY) = Y(IY) + A1*X(I1)
                    I1 = I1 + 1
  150           CONTINUE
                GO TO 400
C
  200           CONTINUE
                J1 = XPNT(1)
                J2 = XPNT(1+1)
                I1 = J2        - MM
                I2 = XPNT(1+2) - MM
                A1 = - X(J1)*X(I1)
                A2 = - X(J2)*X(I2)
                DO  250  IY = IYSTRT, IYSTOP
                    Y(IY) = ( (Y(IY))
     &                      + A1*X(I1)) + A2*X(I2)
                    I1 = I1 + 1
                    I2 = I2 + 1
  250           CONTINUE
                GO TO 400
C
  300           CONTINUE
                J1 = XPNT(1)
                J2 = XPNT(1+1)
                J3 = XPNT(1+2)
                I1 = J2        - MM
                I2 = J3        - MM
                I3 = XPNT(1+3) - MM
                A1 = - X(J1)*X(I1)
                A2 = - X(J2)*X(I2)
                A3 = - X(J3)*X(I3)
                DO  350  IY = IYSTRT, IYSTOP
                    Y(IY) = (( (Y(IY))
     &                      + A1*X(I1)) + A2*X(I2))
     &                      + A3*X(I3)
                    I1 = I1 + 1
                    I2 = I2 + 1
                    I3 = I3 + 1
  350           CONTINUE
                GO TO 400
C
  400           CONTINUE
                DO  600  XCOL = REMAIN, N, LEVEL
                    J1 = XPNT(XCOL)
                    J2 = XPNT(XCOL+1)
                    J3 = XPNT(XCOL+2)
                    J4 = XPNT(XCOL+3)
                    I1 = J2           - MM
                    I2 = J3           - MM
                    I3 = J4           - MM
                    I4 = XPNT(XCOL+4) - MM
                    A1 = - X(J1)*X(I1)
                    A2 = - X(J2)*X(I2)
                    A3 = - X(J3)*X(I3)
                    A4 = - X(J4)*X(I4)
                    DO  500  IY = IYSTRT, IYSTOP
                        Y(IY) = ((( (Y(IY))
     &                          + A1*X(I1)) + A2*X(I2))
     &                          + A3*X(I3)) + A4*X(I4)
                        I1 = I1 + 1
                        I2 = I2 + 1
                        I3 = I3 + 1
                        I4 = I4 + 1
  500               CONTINUE
  600           CONTINUE
C
            MM = MM - 1
            LENY = LENY - 1
C
  700   CONTINUE
C
        RETURN
        END
