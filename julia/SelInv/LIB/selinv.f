C***********************************************************************
C***********************************************************************
C
C   Version:        0.1
C   Last modified:  October 27, 2009
C   Authors:        
C     Chao Yang,
C     Computer Research Division, Lawrence Berkeley National Lab
C     Lin Lin,
C     Program in Applied and Computational Mathematics, Princeton
C     University
C
C***********************************************************************
C***********************************************************************
C
C       ****************************************************************
C       ****************************************************************
C
C       SELINV.F :  SUBROUTINE THAT COMPUTES SELECTED COMPONENTS OF A
C       SPARSE, SYMMETRIC, NONSINGULAR LINEAR SYSTEM BASED ON ITS SPARSE
C       LDL' FACTORIZATION. THIS PROCEDURE IS ALSO CALLED SELECTED
C       INVERSION. THE FACTORIZATION SUBROUTINE IS PROVIDED BY THE
C       SUPERNODAL LEFT-LOOKING LDL' FACTORIZATION CODE WRITTEN BY
C       ESMOND NG AND BARRY PEYTON.
C
C
C       ****************************************************************
C       ****************************************************************
C
C
C       INPUTS PARAMETERS: 
C
C       NSUPER          -   NUMBER OF SUPERNODES.
C       XSUPER          -   THE STARTING COLUMN OF EACH SUPERNODE.
C       XLINDX          -   THE STARTING INDEX IN LINDX FOR EACH SUPERNODE.
C       LINDX           -   ROW INDICES FOR EACH SUPERNODE.
C       XLNZ            -   THE STARTING INDEX OF EACH COLUMN IN THE
C                           CHOLESKY FACTOR L.
C       LNZ             -   NON-ZERO ELEMENTS OF CHOLESKY FACTOR L.
C       SNODES          -   SUPERNODE INDEX FOR EACH ROW
C       DIAG            -   D MATRIX OF LDL' FACTORIZATION.
C       PERM            -   PERMUTATION MATRIX GENERATED IN ANALYSIS
C                           PHASE.
C       NEQNS           -   THE ROW(COLUMN) NUMBER OF A.
C       DUMPL           -   TEST FLAG. IF DUMPL == 1 THEN L IS OUTPUTED.
C
C       OUTPUT PARAMETERS:
C       DIAG            -   DIAGONAL ELEMENTS OF THE INVERSE OF A.
C       LNZ             -   SELECTED INVERSION OF A.
C
C
      SUBROUTINE  SELINV ( NSUPER, XSUPER, XLINDX, LINDX, XLNZ  ,
     &                        LNZ   , SNODES, DIAG  , PERM , NEQNS ,
     &                        DUMPL)
      IMPLICIT NONE
C
C***********************************************************************
C
        INTEGER             NSUPER, NEQNS
        INTEGER             LINDX(*)      , XSUPER(*), PERM(*)
        INTEGER             XLINDX(*)     , XLNZ(*), SNODES(*)
        DOUBLE PRECISION    LNZ(*)        , DIAG(*)
        INTEGER             DUMPL
C
C***********************************************************************
C
        INTEGER             FJCOL , I     , IPNT  , IX    , IXSTOP,
     &                      IXSTRT, JCOL  , JPNT  , JSUP  , LJCOL
        DOUBLE PRECISION    T, DINV
        INTEGER IROW, SNSIZE, IDXREL, JXSTRT, JXSTOP, JX, ISUP, IDEST
        INTEGER IMAT, IVEC, IMATCH, VECROW, MATROW, IDEL, NZVEC
        REAL                GTIMER, T0, T1
        REAL*8,  ALLOCATABLE, DIMENSION(:)   :: DWORK, YWORK
        INTEGER, ALLOCATABLE, DIMENSION(:)   :: NEWXLNZ, IPIV
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDMAP
        INTEGER  SUPSIZE, NNZLPLUS, COLNNZ, NROWS, NCOLS, IX0, ISUP0
        INTEGER  IERR, MAXSUP, LDSUPI, LDSUPJ, NVROWS, MAXWORK,
     &           MAXLDSUP, MATSIZE, NBLKS, IROW0, IVBLK, IB, LDSUPK, KX,
     &           KXSTRT, KXSTOP, IY, IY0, LDY
        REAL*8   ONE, ZERO
        PARAMETER (ONE = 1.0D0, ZERO = 0.0D0 )
        DOUBLE PRECISION DDOT
        REAL     TDIAGINV0, TDIAGINV1
 
C
C***********************************************************************
C
        T0 = GTIMER()
        IF  ( NSUPER .LE. 0 )  RETURN
        IF (DUMPL .EQ. 1) THEN 
           CALL DUMPLMAT(NSUPER, XSUPER, XLINDX, LINDX, XLNZ, LNZ,
     &                   'DEBUGL.M')  
        ENDIF
        ! 
        ! FIND OUT HOW MUCH EXTRA STORAGE NEEDED
        !
        NNZLPLUS = XLNZ(NEQNS+1)-1
C        WRITE(6,*) 'NUMBER OF NONZEROS = ', NNZLPLUS
        MAXSUP = 0
        DO JSUP = 1, NSUPER
           SUPSIZE = XSUPER(JSUP+1)-XSUPER(JSUP)
           IF (SUPSIZE .GT. MAXSUP) MAXSUP = SUPSIZE
           NNZLPLUS = NNZLPLUS + SUPSIZE*(SUPSIZE-1)/2
        END DO
        ALLOCATE(NEWXLNZ(NEQNS+1))
        ALLOCATE(IPIV(MAXSUP))
        ALLOCATE(DWORK(MAXSUP))
        DO I = 1, MAXSUP
           IPIV(I) = I 
        END DO
        
        !
        ! COPY L IN PLACE AND SETUP THE POINTER
        !  
        TDIAGINV0 = GTIMER()
        NEWXLNZ(NEQNS+1) = NNZLPLUS+1
        MAXWORK = 0
        MAXLDSUP = 0
        DO JSUP = NSUPER,1,-1
           FJCOL = XSUPER(JSUP)
           LJCOL = XSUPER(JSUP+1)-1
           COLNNZ = XLNZ(FJCOL+1) - XLNZ(FJCOL)
           DO JCOL = LJCOL,FJCOL,-1
              IXSTRT = XLNZ(JCOL) 
              IXSTOP = XLNZ(JCOL+1)-1
              DO IX = IXSTOP, IXSTRT, -1
                  !
                  ! DESTINATION MUST BE BELOW THE DIAGONAL
                  !
                  IDEST = NEWXLNZ(JCOL+1)-1
                  LNZ(IDEST+(IX-IXSTOP)) = LNZ(IX)
              END DO 
              NEWXLNZ(JCOL) = NEWXLNZ(JCOL+1) - COLNNZ 
           END DO 
           IF (JSUP .LT. NSUPER) THEN
              MATSIZE = COLNNZ*(LJCOL-FJCOL+1)
              IF (COLNNZ .GT. MAXLDSUP) MAXLDSUP = COLNNZ
              IF (MATSIZE .GT. MAXWORK) MAXWORK = MATSIZE
           ENDIF           
        END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
c$$$        open(unit=10,file="LNZ2.txt") 
c$$$        DO IX = 1, NNZLPLUS
c$$$           write(10,*) LNZ(IX)
c$$$        END DO
c$$$        close(10)
c$$$
c$$$        open(unit=10,file="NEWXLNZ2.txt") 
c$$$        DO IX = 1, NEQNS+1
c$$$           write(10,*) NEWXLNZ(IX)
c$$$        END DO
c$$$        close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


C        WRITE(6,*) 'NEWNNZL = ', NEWXLNZ(NEQNS+1)-1
C        WRITE(6,*) 'MAXWORK = ', MAXWORK
        !
        !  MODIFY ALL SUPERNODES
        !
        DO JSUP = 1, NSUPER
           IERR = 0
           FJCOL = XSUPER(JSUP)
           LJCOL = XSUPER(JSUP+1)-1
           SUPSIZE = LJCOL-FJCOL+1
           COLNNZ = NEWXLNZ(FJCOL+1) - NEWXLNZ(FJCOL)
           NROWS = COLNNZ-SUPSIZE
           NCOLS = SUPSIZE
           IXSTRT = NEWXLNZ(FJCOL) 
           IF (SUPSIZE .GT. 1) THEN 
              ! 
              !  PERFORM A TRIANGULAR SOLVE BELOW THE DIAGONAL BLOCK
              ! 
              CALL DTRSM('RIGHT','LOWER','NOTRANSPOSE','UNIT', 
     &                   NROWS, NCOLS, ONE, LNZ(IXSTRT), COLNNZ, 
     &                   LNZ(IXSTRT+SUPSIZE), COLNNZ)
              !
              !  INVERT THE DIAGONAL BLOCK
              !  
              CALL DSYTRI('LOWER',SUPSIZE,LNZ(IXSTRT),COLNNZ,
     &                    IPIV, DWORK, IERR)
              IF (IERR .NE. 0) THEN
                 WRITE(6,*) 'TRIANGULAR SOLVE FAILED, IERR = ', IERR
                 WRITE(6,*) 'JSUP = ', JSUP, 'SUPSIZE =', SUPSIZE
              ENDIF 
              !
              ! NEED TO STORE THE UPPER TRIANGULAR PART OF THE DIAGONAL
              ! BLOCK FOR FOR DGEMM
              !
              DO IROW = 1,SUPSIZE-1
                 DO JCOL = IROW+1,SUPSIZE
                    LNZ(IXSTRT+(JCOL-1)*COLNNZ+IROW-1) 
     &              = LNZ(IXSTRT+(IROW-1)*COLNNZ + JCOL-1)
                 END DO
              END DO
            ELSE
              LNZ(IXSTRT) = 1/LNZ(IXSTRT)
           ENDIF
        END DO
        TDIAGINV1 = GTIMER()
C        WRITE(6,123) TDIAGINV1 - TDIAGINV0
C 123    FORMAT(1X,' TIME FOR INVERSE OF DIAGONAL BLOCKS= ', 1PE11.3)
    
        ALLOCATE(YWORK(MAXWORK))
C       LL: INITIALIZE THE VARIABLE, VERY IMPORTANT! 
        DO I = 1, MAXWORK
          YWORK(I) = 0
        ENDDO
        ALLOCATE(INDMAP(MAXLDSUP,3))
        DO JSUP = NSUPER-1, 1, -1
           FJCOL = XSUPER(JSUP)
           LJCOL = XSUPER(JSUP+1)-1
           SUPSIZE = LJCOL-FJCOL+1
           JPNT = XLINDX(JSUP)+SUPSIZE ! POINT TO THE BEGINNING OF
                                       ! THE OFF-DIAGONAL BLOCKS  
                                       ! COLUMN
           IXSTRT = NEWXLNZ(FJCOL)+SUPSIZE 
           IXSTOP = NEWXLNZ(FJCOL+1)-1
           LDSUPJ = NEWXLNZ(FJCOL+1)-NEWXLNZ(FJCOL)
           LDY    = IXSTOP - IXSTRT + 1
C          LL: DO NOT DO ANYTHING IF THERE IS NO OFF DIAGONAL           
           IF (LDY .GT. 0) THEN
           !
           ! CONSTRUCT A BLOCKMAP LIST FOR THE JSUP-TH SUPERNODE
           ! A BLOCK ROW IS DEFINED TO BE A SET OF ROWS WITHIN 
           ! THE SAME SUPERNODE THAT HAVE CONSEQUTIVE ROW INDICES
           !
           ! INDMAP(IB,1) --- ROW INDEX OF THE FIRST ROW IN THE BLOCK
           ! INDMAP(IB,2) --- THE SUPERNODE THE BLOCK BELONGS TO
           ! INDMAP(IB,3) --- NUMBER OF ROWS IN THE BLOCK
           !
           NBLKS = 0
           IROW0 = -1
           ISUP0 = -1
           NROWS = 0
           IPNT  = JPNT
           DO IX = IXSTRT, IXSTOP
              IROW = LINDX(IPNT)
              ISUP = SNODES(IROW)
              IF (IROW .EQ. IROW0 + 1 .AND. ISUP .EQ. ISUP0 ) THEN
                 !
                 ! WITHIN THE SAME BLOCK
                 !
                 NROWS = NROWS + 1
              ELSE
                 !
                 ! START A NEW BLOCK
                 !
                 ISUP0 = ISUP
                 NBLKS = NBLKS + 1
                 INDMAP(NBLKS,1) = IROW
                 INDMAP(NBLKS,2) = ISUP
                 IF (NBLKS .GT. 1) THEN
                    INDMAP(NBLKS-1,3) = NROWS
                 ENDIF 
                 NROWS = 1
              ENDIF
              IROW0 = IROW
              IPNT = IPNT + 1               
           ENDDO
           IF (NBLKS .GE. 1) THEN
              INDMAP(NBLKS,3) = NROWS 
           ENDIF 
           !  DO Y = -S*X, 
           !
           !  IX   POINTER TO THE NONZERO VALUE ARRAY
           !  KX   POINTER TO THE NONZERO VALUE ARRAY
           !  IY   POINTER TO THE NONZERO VALUE WORK ARRAY Y
           !  IMAT POINTER TO ROW INDEX ARRAY
           !
           !  STEP THRU BLOCKS OF X, EACH BLOCK BELONGS TO A SUPERNODE
           !
           IF (NBLKS .GT. 0) THEN
              IROW0 = INDMAP(1,1)
              IY0 = 1 
              IX0 = IXSTRT
              DO IB = 1, NBLKS
                 IROW  = INDMAP(IB,1)
                 ISUP  = INDMAP(IB,2)
                 NCOLS = INDMAP(IB,3) ! IT IS THE NUMER OF ROWS IN A NONZERO
                                      ! BLOCK OF X. IT CORRESPONDS TO THE NUMBER
                                      ! OF COLUMNS IN THE TRAILING SCHUR 
                                      ! COMPLEMENT THAT WILL BE USED
                 IVBLK = IB
                 !
                 KXSTRT = NEWXLNZ(IROW) 
                 KXSTOP = NEWXLNZ(IROW+1)-1
                 IMAT   = XLINDX(ISUP)  ! POINTER TO THE ROW INDEX ARRAY THAT 
                                        ! CORRESPONDS TO THE FIRST NONZERO
                                        ! ENTRY OF THE ISUP-TH SUPERNODE
                 LDSUPK = KXSTOP-KXSTRT+1
                 IMATCH = 0
                 IY = IY0
                 IX = IX0
                 KX = KXSTRT
                 !
                 ! WALL DOWN THE IROW-TH COLUMN (WHICH BELONGS TO ISUP-TH
                 ! SUPERNODE) OF THE TRAILING SCHUR COMPLEMENT. PERFORM 
                 ! MATRIX OPERATIONS WHEN THE ROW INDEX OF THE SCHUR COMPLEMENT
                 ! MATCHES THE FIRST ROW OF A NON-ZERO BLOCK OF X
                 !
                 DO WHILE (KX .LE. KXSTOP)
                    MATROW = LINDX(IMAT)
                    VECROW = INDMAP(IVBLK,1)
                    NVROWS = INDMAP(IVBLK,3) 
                    IF (MATROW .LT. VECROW) THEN
                       ! 
                       ! SKIP THIS ROW
                       !
                       IMAT = IMAT + 1
                       KX   = KX + 1 
                    ELSE 
                       ! MATROW MUST MATCH VECROW
                       !
                       ! THE LOWER TRIANGULAR CONTRIBUTION (AXPY)
                       !
                       IMATCH = IMATCH + 1
                       !
                       ! SEPARATE THE FOLLOWING TWO CASES FOR PERFORMANCE
                       ! REASON
                       !
                       IF (SUPSIZE .EQ. 1 ) THEN
                          IF ( NCOLS .EQ. 1) THEN
                             CALL DAXPY(NVROWS,-LNZ(IX), LNZ(KX),
     &                                  1, YWORK(IY), 1)
                          ELSE
                             CALL DGEMV('N',NVROWS,NCOLS,-1.0D0,
     &                                  LNZ(KX),LDSUPK, 
     &                                  LNZ(IX), 1, 1.0D0,
     &                                  YWORK(IY),1)
                          ENDIF 
                       ELSE
                          CALL DGEMM('N','N',NVROWS,SUPSIZE,NCOLS,
     &                               -1.0D0,LNZ(KX),LDSUPK,
     &                               LNZ(IX),LDSUPJ,1.0D0, 
     &                               YWORK(IY), LDY)
                       ENDIF
                       !
                       ! THE UPPER TRIANGULAR CONTRIBUTION (DOT)
                       !
                       IF (IMATCH .GT. 1) THEN
                          !
                          ! SEPARATE THE FOLLOWING TWO CASES
                          ! FOR PERFORMANCE REASON
                          ! 
                          IF (SUPSIZE .EQ. 1 ) THEN
                             IF (NCOLS .EQ. 1) THEN
                                YWORK(IY0) = YWORK(IY0)
     &                                  - DDOT(NVROWS,LNZ(KX),1,
     &                                         LNZ(IX+IY-IY0),1)
                             ELSE
                                CALL DGEMV('T',NVROWS,NCOLS,-1.0D0,
     &                                     LNZ(KX),LDSUPK,
     &                                     LNZ(IX+IY-IY0),1,
     &                                     1.0D0, YWORK(IY0),1)
                             ENDIF 
                          ELSE
                             CALL DGEMM('T','N',NCOLS,SUPSIZE,NVROWS,
     &                                  -1.0D0,LNZ(KX),LDSUPK,
     &                                  LNZ(IX+IY-IY0),LDSUPJ,
     &                                  1.0D0, YWORK(IY0),LDY)
                          ENDIF
                       ENDIF
                       !
                       ! MOVE ON TO THE NEXT BLOCK IN THE SCHUR COMPLEMENT
                       !
                       IVBLK = IVBLK + 1
                       IMAT = IMAT + NVROWS
                       KX   = KX + NVROWS
                       IY   = IY + NVROWS
                       IF (IMATCH .GE. NBLKS-IB+1) GOTO 20
                    ENDIF
                 ENDDO !WHILE KX
 20              CONTINUE
                 !
                 ! REMEMBER THAT NCOLS IS ACTUALLY THE NUMBER OF ROWS
                 ! IN A NONZERO BLOCK OF X THAT HAVE JUST BEEN PROCESSED.
                 !
                 IY0 = IY0 + NCOLS
                 IX0 = IX0 + NCOLS
              END DO ! IB
           ENDIF ! NBLKS > 0
           !
           ! DO G = G + X'*SX
           !
           IX = IXSTRT
           JX = NEWXLNZ(FJCOL)
           IY = 1
           IF (SUPSIZE .EQ. 1) THEN
              LNZ(JX) = LNZ(JX)
     &                   - DDOT(IXSTOP-IXSTRT+1,LNZ(IXSTRT),1,
     &                          YWORK(IY),1)  
           ELSE
              CALL DGEMM('T','N',SUPSIZE,SUPSIZE, IXSTOP-IXSTRT+1,
     &                   -1.0D0, LNZ(IXSTRT), LDSUPJ, YWORK(1),
     &                   LDY, 1.0D0, LNZ(JX), LDSUPJ) 
              ! 
              ! SYMMETRIZE THE DIAGONAL BLOCK IS IMPORTANT  
              ! FOR STABILITY (LL: 12/18/2012)
              !
              DO IROW = 1,SUPSIZE-1
                 DO JCOL = IROW+1,SUPSIZE
                    LNZ(JX+(JCOL-1)*LDSUPJ+IROW-1) = 
     &                0.5d0 * ( LNZ(JX+(IROW-1)*LDSUPJ+JCOL-1) + 
     &                LNZ(JX+(JCOL-1)*LDSUPJ+IROW-1) )
                    LNZ(JX+(IROW-1)*LDSUPJ+JCOL-1) =
     &                LNZ(JX+(JCOL-1)*LDSUPJ+IROW-1)
                 END DO
              END DO
           ENDIF
           CALL DLACPY('ALL',(IXSTOP-IXSTRT+1),SUPSIZE,YWORK,LDY,
     &                 LNZ(IXSTRT),LDSUPJ)
           ! 
           ! MUST ZERO OUT THE WORKSPACE
           !
           DO I = 1, LDY*SUPSIZE
              YWORK(I) = 0.0D0
           ENDDO
           ENDIF  ! IF(LDY .GT. 0)
        ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!        open(unit=10,file="LNZ3.txt") 
!        DO IX = 1, NNZLPLUS
!           write(10,*) LNZ(IX)
!        END DO
!        close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        !
        ! PERMUTE THE DIAGONAL BACK TO ITS ORIGINAL ORDER
        !
        DO JSUP = 1, NSUPER
           FJCOL = XSUPER(JSUP)
           LJCOL = XSUPER(JSUP+1)-1
           DO JCOL = FJCOL,LJCOL
              JX = NEWXLNZ(JCOL)+JCOL-FJCOL
!             DIAG(PERM(JCOL)) = LNZ(JX)
              DIAG(JCOL) = JX
           ENDDO
        ENDDO


!!!!!!!!!!!!!!!!!!
        T1 = GTIMER()
        WRITE(6,333) T1-T0

  333   FORMAT(1X,'TIME FOR SELECTED INVERSION = ',1PE11.3)

        DEALLOCATE(NEWXLNZ)
        DEALLOCATE(IPIV)
        DEALLOCATE(DWORK)
        DEALLOCATE(YWORK)
        DEALLOCATE(INDMAP)
        RETURN
      END
