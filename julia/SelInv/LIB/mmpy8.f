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
C*************     MMPY8  .... MATRIX-MATRIX MULTIPLY     **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE -
C       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA,
C       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY
C       CODES.
C
C       LOOP UNROLLING: LEVEL 8
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
      SUBROUTINE  MMPY8  (  M     , N     , Q     , XPNT  , X     ,
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
        INTEGER             I1    , I2    , I3    , I4    , I5    ,
     &                      I6    , I7    , I8
        INTEGER             J1    , J2    , J3    , J4    , J5    ,
     &                      J6    , J7    , J8
        INTEGER             IY    , IYLAST, IYSTRT, IYSTOP, LENY  ,
     &                      MM    , REMAIN, XCOL  , YCOL
        DOUBLE PRECISION    A1    , A2    , A3    , A4    , A5    ,
     &                      A6    , A7    , A8
C
        INTEGER             LEVEL
        PARAMETER       (   LEVEL = 8 )
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
        DO  1100  YCOL = 1, Q
C
            IYSTRT = IYLAST + 1
            IYSTOP = IYSTRT + MM - 1
            IYLAST = IYLAST + LENY
C
C           --------------------------------------------------
C           ... PERFORM THE APPROPRATE MATRIX VECTOR MULTIPLY:
C               X * A(*,YCOL) WITH LEVEL 8 LOOP-UNROLLING.
C           --------------------------------------------------
C
            GO TO ( 800, 100, 200, 300, 400, 500, 600, 700 ), REMAIN
C
  100           CONTINUE
                J1 = XPNT(1)
                I1 = XPNT(1+1) - MM
                A1  = - X(J1)*X(I1)
                DO  150  IY = IYSTRT, IYSTOP
                    Y(IY) = Y(IY) + A1*X(I1)
                    I1 = I1 + 1
  150           CONTINUE
                GO TO 800
C
  200           CONTINUE
                J1 = XPNT(1)
                J2 = XPNT(1+1)
                I1 = J2        - MM
                I2 = XPNT(1+2) - MM
                A1  = - X(J1)*X(I1)
                A2  = - X(J2)*X(I2)
                DO  250  IY = IYSTRT, IYSTOP
                    Y(IY) = ( (Y(IY))
     &                      + A1*X(I1)) + A2*X(I2)
                    I1 = I1 + 1
                    I2 = I2 + 1
  250           CONTINUE
                GO TO 800
C
  300           CONTINUE
                J1 = XPNT(1)
                J2 = XPNT(1+1)
                J3 = XPNT(1+2)
                I1 = J2        - MM
                I2 = J3        - MM
                I3 = XPNT(1+3) - MM
                A1  = - X(J1)*X(I1)
                A2  = - X(J2)*X(I2)
                A3  = - X(J3)*X(I3)
                DO  350  IY = IYSTRT, IYSTOP
                    Y(IY) = (( (Y(IY))
     &                      + A1*X(I1)) + A2*X(I2))
     &                      + A3*X(I3)
                    I1 = I1 + 1
                    I2 = I2 + 1
                    I3 = I3 + 1
  350           CONTINUE
                GO TO 800
C
  400           CONTINUE
                J1 = XPNT(1)
                J2 = XPNT(1+1)
                J3 = XPNT(1+2)
                J4 = XPNT(1+3)
                I1 = J2        - MM
                I2 = J3        - MM
                I3 = J4        - MM
                I4 = XPNT(1+4) - MM
                A1  = - X(J1)*X(I1)
                A2  = - X(J2)*X(I2)
                A3  = - X(J3)*X(I3)
                A4  = - X(J4)*X(I4)
                DO  450  IY = IYSTRT, IYSTOP
                    Y(IY) = ((( (Y(IY))
     &                      + A1*X(I1)) + A2*X(I2))
     &                      + A3*X(I3)) + A4*X(I4)
                    I1 = I1 + 1
                    I2 = I2 + 1
                    I3 = I3 + 1
                    I4 = I4 + 1
  450           CONTINUE
                GO TO 800
C
  500           CONTINUE
                J1 = XPNT(1)
                J2 = XPNT(1+1)
                J3 = XPNT(1+2)
                J4 = XPNT(1+3)
                J5 = XPNT(1+4)
                I1 = J2        - MM
                I2 = J3        - MM
                I3 = J4        - MM
                I4 = J5        - MM
                I5 = XPNT(1+5) - MM
                A1  = - X(J1)*X(I1)
                A2  = - X(J2)*X(I2)
                A3  = - X(J3)*X(I3)
                A4  = - X(J4)*X(I4)
                A5  = - X(J5)*X(I5)
                DO  550  IY = IYSTRT, IYSTOP
                    Y(IY) = (((( (Y(IY))
     &                      + A1*X(I1)) + A2*X(I2))
     &                      + A3*X(I3)) + A4*X(I4))
     &                      + A5*X(I5)
                    I1 = I1 + 1
                    I2 = I2 + 1
                    I3 = I3 + 1
                    I4 = I4 + 1
                    I5 = I5 + 1
  550           CONTINUE
                GO TO 800
C
  600           CONTINUE
                J1 = XPNT(1)
                J2 = XPNT(1+1)
                J3 = XPNT(1+2)
                J4 = XPNT(1+3)
                J5 = XPNT(1+4)
                J6 = XPNT(1+5)
                I1 = J2        - MM
                I2 = J3        - MM
                I3 = J4        - MM
                I4 = J5        - MM
                I5 = J6        - MM
                I6 = XPNT(1+6) - MM
                A1  = - X(J1)*X(I1)
                A2  = - X(J2)*X(I2)
                A3  = - X(J3)*X(I3)
                A4  = - X(J4)*X(I4)
                A5  = - X(J5)*X(I5)
                A6  = - X(J6)*X(I6)
                DO  650  IY = IYSTRT, IYSTOP
                    Y(IY) = ((((( (Y(IY))
     &                      + A1*X(I1)) + A2*X(I2))
     &                      + A3*X(I3)) + A4*X(I4))
     &                      + A5*X(I5)) + A6*X(I6)
                    I1 = I1 + 1
                    I2 = I2 + 1
                    I3 = I3 + 1
                    I4 = I4 + 1
                    I5 = I5 + 1
                    I6 = I6 + 1
  650           CONTINUE
                GO TO 800
C
  700           CONTINUE
                J1 = XPNT(1)
                J2 = XPNT(1+1)
                J3 = XPNT(1+2)
                J4 = XPNT(1+3)
                J5 = XPNT(1+4)
                J6 = XPNT(1+5)
                J7 = XPNT(1+6)
                I1 = J2        - MM
                I2 = J3        - MM
                I3 = J4        - MM
                I4 = J5        - MM
                I5 = J6        - MM
                I6 = J7        - MM
                I7 = XPNT(1+7) - MM
                A1  = - X(J1)*X(I1)
                A2  = - X(J2)*X(I2)
                A3  = - X(J3)*X(I3)
                A4  = - X(J4)*X(I4)
                A5  = - X(J5)*X(I5)
                A6  = - X(J6)*X(I6)
                A7  = - X(J7)*X(I7)
                DO  750  IY = IYSTRT, IYSTOP
                    Y(IY) = (((((( (Y(IY))
     &                      + A1*X(I1)) + A2*X(I2))
     &                      + A3*X(I3)) + A4*X(I4))
     &                      + A5*X(I5)) + A6*X(I6))
     &                      + A7*X(I7)
                    I1 = I1 + 1
                    I2 = I2 + 1
                    I3 = I3 + 1
                    I4 = I4 + 1
                    I5 = I5 + 1
                    I6 = I6 + 1
                    I7 = I7 + 1
  750           CONTINUE
                GO TO 800
C
  800           CONTINUE
                DO  1000  XCOL = REMAIN, N, LEVEL
                    J1 = XPNT(XCOL)
                    J2 = XPNT(XCOL+1)
                    J3 = XPNT(XCOL+2)
                    J4 = XPNT(XCOL+3)
                    J5 = XPNT(XCOL+4)
                    J6 = XPNT(XCOL+5)
                    J7 = XPNT(XCOL+6)
                    J8 = XPNT(XCOL+7)
                    I1 = J2           - MM
                    I2 = J3           - MM
                    I3 = J4           - MM
                    I4 = J5           - MM
                    I5 = J6           - MM
                    I6 = J7           - MM
                    I7 = J8           - MM
                    I8 = XPNT(XCOL+8) - MM
                    A1  = - X(J1)*X(I1)
                    A2  = - X(J2)*X(I2)
                    A3  = - X(J3)*X(I3)
                    A4  = - X(J4)*X(I4)
                    A5  = - X(J5)*X(I5)
                    A6  = - X(J6)*X(I6)
                    A7  = - X(J7)*X(I7)
                    A8  = - X(J8)*X(I8)
                    DO  900  IY = IYSTRT, IYSTOP
                        Y(IY) = ((((((( (Y(IY))
     &                          + A1*X(I1)) + A2*X(I2))
     &                          + A3*X(I3)) + A4*X(I4))
     &                          + A5*X(I5)) + A6*X(I6))
     &                          + A7*X(I7)) + A8*X(I8)
                        I1 = I1 + 1
                        I2 = I2 + 1
                        I3 = I3 + 1
                        I4 = I4 + 1
                        I5 = I5 + 1
                        I6 = I6 + 1
                        I7 = I7 + 1
                        I8 = I8 + 1
  900               CONTINUE
 1000           CONTINUE
C
            MM = MM - 1
            LENY = LENY - 1
C
 1100   CONTINUE
C
        RETURN
        END
