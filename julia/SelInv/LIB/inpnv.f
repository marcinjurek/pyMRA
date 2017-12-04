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
C
C     ------------------------------------------------------
C     INPUT NUMERICAL VALUES INTO SPARSE DATA STRUCTURES ...
C     ------------------------------------------------------
C
      SUBROUTINE  INPNV  (  OUTUNT, NEQNS, XADJF, ADJF, ANZF, PERM, 
     &                      INVP, NSUPER, XSUPER, XLINDX, LINDX,
     &                      XLNZ, LNZ, IWSIZ, OFFSET, IFLAG )
C
        INTEGER             XADJF(*), ADJF(*)
        DOUBLE PRECISION    ANZF(*)
        INTEGER             PERM(*), INVP(*)
        INTEGER             NEQNS, NSUPER
        INTEGER             XSUPER(*), XLINDX(*), LINDX(*)
        INTEGER             XLNZ(*)
        DOUBLE PRECISION    LNZ(*)
        INTEGER             OFFSET(*)
C
        INTEGER             IFLAG, IWSIZ, OUTUNT
        INTEGER             I, II, J, JLEN, JSUPER, LAST, OLDJ
C
        IFLAG = 0
        IF  ( IWSIZ .LT. NEQNS )  THEN
            IFLAG = -1
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) '*** INTEGER WORK SPACE = ',
     &                       IWSIZ
            WRITE (OUTUNT,*) '*** IS SMALLER THAN REQUIRED = ',
     &                       NEQNS
            RETURN
        ENDIF
C
        DO  500  JSUPER = 1, NSUPER
C
C           ----------------------------------------
C           FOR EACH SUPERNODE, DO THE FOLLOWING ...
C           ----------------------------------------
C
C           -----------------------------------------------
C           FIRST GET OFFSET TO FACILITATE NUMERICAL INPUT.
C           -----------------------------------------------
            JLEN = XLINDX(JSUPER+1) - XLINDX(JSUPER)
            DO  100  II = XLINDX(JSUPER), XLINDX(JSUPER+1)-1
                I = LINDX(II)
                JLEN = JLEN - 1
                OFFSET(I) = JLEN
  100       CONTINUE
C
            DO  400  J = XSUPER(JSUPER), XSUPER(JSUPER+1)-1
C               -----------------------------------------
C               FOR EACH COLUMN IN THE CURRENT SUPERNODE,
C               FIRST INITIALIZE THE DATA STRUCTURE.
C               -----------------------------------------
                DO  200  II = XLNZ(J), XLNZ(J+1)-1
                    LNZ(II) = 0.0D+00
  200           CONTINUE
C
C               -----------------------------------
C               NEXT INPUT THE INDIVIDUAL NONZEROS.
C               -----------------------------------
                OLDJ = PERM(J)
                LAST = XLNZ(J+1) - 1
                DO  300  II = XADJF(OLDJ), XADJF(OLDJ+1)-1
                    I = INVP(ADJF(II))
                    IF  ( I .GE. J )  THEN
                        LNZ(LAST-OFFSET(I)) = ANZF(II)
                    ENDIF
  300           CONTINUE
  400       CONTINUE
C
  500   CONTINUE
        RETURN
      END
