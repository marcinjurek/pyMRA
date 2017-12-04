C***********************************************************************
C***********************************************************************
C
C   Version:        0.3
C   Last modified:  March 6, 1995
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C*********     BLKFCT .....  BLOCK GENERAL SPARSE LDL'         *********
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS SUBROUTINE CALLS THE BLOCK GENERAL SPARSE LDL' ROUTINE,
C       BLKFC2.
C
C   INPUT PARAMETERS:
C       OUTUNT          -   OUTPUT UNIT.
C       NSUPER          -   NUMBER OF SUPERNODES.
C       NUNROL          -   LOOP UNROLLING LEVEL.
C       XSUPER          -   SUPERNODE PARTITION.
C       SNODE           -   MAPS EACH COLUMN TO THE SUPERNODE CONTAINING
C                           IT.
C       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT
C                           INTO CACHE.
C       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE (INCLUDING
C                           THE DIAGONAL ELEMENTS).
C       (XLNZ,LNZ)      -   ON INPUT, CONTAINS MATRIX TO BE FACTORED.
C       IWSIZ           -   SIZE OF INTEGER WORKING STORAGE
C       TMPSIZ          -   SIZE OF FLOATING POINT WORKING STORAGE.
C       MMPYN           -   EXTERNAL ROUTINE: MATRIX-MATRIX MULTIPLY.
C       SMXPY           -   EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY.
C
C   OUTPUT PARAMETERS:
C       LNZ             -   ON OUTPUT, CONTAINS TRIANGULAR FACTOR
C                           AND DIAGONAL MATRIX OF LDL' DECOMPOSITION.
C       DIAG            -   DIAGONAL MATRIX OF LDL' DECOMPOSITION.
C       IFLAG           -   ERROR FLAG.
C                               0: SUCCESSFUL FACTORIZATION.
C                              -1: ZERO DIAGONAL ENCOUNTERED,
C                                  MATRIX IS SINGULAR.
C                              -2: INSUFFICIENT WORKING STORAGE 
C                                  [TEMP(*)].
C                              -3: INSUFFICIENT WORKING STORAGE 
C                                  [IWORK(*)].
C
C   WORKING PARAMETERS:
C       IWORK           -   INTEGER WORKING STORAGE OF LENGTH 
C                           2*NEQNS + 2*NSUPER.
C       TMPVEC          -   DOUBLE PRECISION WORKING STORAGE OF LENGTH
C                           NEQNS.
C       
C***********************************************************************
C
      SUBROUTINE  BLKFCT (  OUTUNT, NEQNS , NSUPER, NUNROL, XSUPER, 
     &                      SNODE , SPLIT , XLINDX, LINDX , XLNZ, 
     &                      LNZ   , DIAG  , IWSIZ , IWORK , TMPSIZ, 
     &                      TMPVEC, IFLAG)
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
C
        EXTERNAL            MMPY1 , SMXPY1,
     &                      MMPY2 , SMXPY2,
     &                      MMPY4 , SMXPY4,
     &                      MMPY8 , SMXPY8
        INTEGER             XLINDX(*)     , XLNZ(*)
        INTEGER             IWORK(*)      , LINDX(*)      , 
     &                      SNODE(*)      , SPLIT(*)      , 
     &                      XSUPER(*)
        INTEGER             IFLAG , IWSIZ , NEQNS , NSUPER, OUTUNT,
     &                      TMPSIZ, NUNROL
        DOUBLE PRECISION    DIAG(*)       , LNZ(*)        , 
     &                      TMPVEC(*)
C
C       ---------------
C       LOCAL VARIABLE.
C       ---------------
        INTEGER             JCOL
C
C*********************************************************************
C
        IFLAG = 0
        IF  ( IWSIZ .LT. 2*NEQNS+2*NSUPER )  THEN
            IFLAG = -3
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*) '*** INTEGER WORK SPACE = ', 
     &                       IWSIZ 
            WRITE (OUTUNT,*) '*** IS SMALLER THAN REQUIRED = ',
     &                       2*NEQNS+2*NSUPER
            RETURN
        ENDIF
        IFLAG = 0
        IF  ( IWSIZ .LT. 2*NEQNS+2*NSUPER )  THEN
            IFLAG = -3
            RETURN
        ENDIF

chao_begin
        if (nunrol .eq. 8) then
           CALL  BLKFC2 (  NSUPER, XSUPER, SNODE , SPLIT , 
     &                     XLINDX, LINDX , XLNZ  , LNZ   , 
     &                     IWORK(1)                      ,
     &                     IWORK(NSUPER+1)               ,
     &                     IWORK(2*NSUPER+1)             ,
     &                     IWORK(2*NSUPER+NEQNS+1)       ,
     &                     TMPSIZ, TMPVEC, IFLAG , MMPY8 ,
     &                     SMXPY8)
        elseif (nunrol .eq. 4) then
           CALL  BLKFC2 (  NSUPER, XSUPER, SNODE , SPLIT , 
     &                     XLINDX, LINDX , XLNZ  , LNZ   , 
     &                     IWORK(1)                      ,
     &                     IWORK(NSUPER+1)               ,
     &                     IWORK(2*NSUPER+1)             ,
     &                     IWORK(2*NSUPER+NEQNS+1)       ,
     &                     TMPSIZ, TMPVEC, IFLAG , MMPY4 ,
     &                     SMXPY4)
        elseif (nunrol .eq. 2) then
           CALL  BLKFC2 (  NSUPER, XSUPER, SNODE , SPLIT , 
     &                     XLINDX, LINDX , XLNZ  , LNZ   , 
     &                     IWORK(1)                      ,
     &                     IWORK(NSUPER+1)               ,
     &                     IWORK(2*NSUPER+1)             ,
     &                     IWORK(2*NSUPER+NEQNS+1)       ,
     &                     TMPSIZ, TMPVEC, IFLAG , MMPY2 ,
     &                     SMXPY2)
        else
           CALL  BLKFC2 (  NSUPER, XSUPER, SNODE , SPLIT , 
     &                     XLINDX, LINDX , XLNZ  , LNZ   , 
     &                     IWORK(1)                      ,
     &                     IWORK(NSUPER+1)               ,
     &                     IWORK(2*NSUPER+1)             ,
     &                     IWORK(2*NSUPER+NEQNS+1)       ,
     &                     TMPSIZ, TMPVEC, IFLAG , MMPY1 ,
     &                     SMXPY1)
        end if
chao_end
        IF  ( IFLAG .EQ. -1 )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*)
     &          '*** MATRIX IS SINGULAR ***'
            WRITE (OUTUNT,*)
     &          '*** ZERO DIAGONAL ENTRY ENCOUNTERED ***'
            RETURN
        ELSEIF ( IFLAG .EQ. -2 )  THEN
            WRITE (OUTUNT,*) ' '
            WRITE (OUTUNT,*)
     &          '*** INSUFFICIENT WORK STORAGE [TMPVEC(*)] ***'
            RETURN
        ENDIF
        IF  ( IFLAG .EQ. 0 )  THEN
            DO  100  JCOL = 1, NEQNS
                DIAG(JCOL) = LNZ(XLNZ(JCOL))
  100       CONTINUE
        ENDIF
        RETURN
      END
