C***********************************************************************
C***********************************************************************
C*********     BLKDINV ... Compute the diagonal of the inverse *********
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       GIVEN THE CHOLESKY FACTORIZATION OF A SPARSE SYMMETRIC
C       POSITIVE DEFINITE MATRIX, THIS SUBROUTINE computes the
C       diagonal of the inverse of the given matrix.  
C       IT USES OUTPUT FROM BLKFCT.
C
C   INPUT PARAMETERS:
C       NSUPER          -   NUMBER OF SUPERNODES.
C       XSUPER          -   SUPERNODE PARTITION.
C       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE.
C       (XLNZ,LNZ)      -   CHOLESKY FACTOR.
C
C   UPDATED PARAMETERS:
C       DINV            -   OUTPUT, CONTAINS THE diagonal of the inverse
C
C   for j = 1:n
C      solve L*y = e_j
C      w = D\y
C      dinv(j) = y^T*w
C   end for
C
C***********************************************************************
C
      SUBROUTINE  BLKDINV (  NSUPER, XSUPER, XLINDX, LINDX , XLNZ  ,
     &                       LNZ   , DINVP , Y,      mark, nzidx, perm )
      implicit none 
C
C***********************************************************************
C
        INTEGER             NSUPER
        INTEGER             LINDX(*)      , XSUPER(*), perm(*)
        INTEGER             XLINDX(*)     , XLNZ(*)
        DOUBLE PRECISION    LNZ(*)        , DINVP(*)
        REAL*8, allocatable, dimension(:) :: DINV
        DOUBLE PRECISION    Y(*)
C
C***********************************************************************
C
        INTEGER             FJCOL , I     , IPNT  , IX    , IXSTOP,
     &                      IXSTRT, JCOL  , JPNT  , JSUP  , LJCOL
        DOUBLE PRECISION    T
        integer             neqns, jsup2, jcol2
        integer             mark(*), nzidx(*)
        integer             nnz, j, jx

C
C***********************************************************************
C
        IF  ( NSUPER .LE. 0 )  RETURN
C
        neqns = XSUPER(NSUPER+1)-1
        allocate(dinv(neqns))

        DO J = 1, NEQNS
           DINV(J) = 0.0d0 
           DINVP(J) = 0.0d0 
           mark(j) = 0
c           nzidx(j) = -1;
           Y(J) = 0.0d0
        ENDDO

        DO jsup2 = 1,nsuper
           DO jcol2 = XSUPER(jsup2),XSUPER(jsup2+1)-1
C
              Y(jcol2) = 1.0d0
              mark(jcol2) = 1
              nnz = 1
              nzidx(nnz) = jcol2
C
C             ------------------------
C             FORWARD SUBSTITUTION ...
C             ------------------------
C
              FJCOL = jcol2
              DO JSUP = jsup2, NSUPER
                 LJCOL  = XSUPER(JSUP+1) - 1
                 IXSTRT = XLNZ(FJCOL)
                 JPNT   = XLINDX(JSUP)+ FJCOL-XSUPER(JSUP)
                 DO JCOL = FJCOL, LJCOL
                    IXSTOP    = XLNZ(JCOL+1) - 1
                    T         = Y(JCOL)
                    IF (T .NE. 0.0) then
                       IPNT      = JPNT +1
                       DO IX = IXSTRT+1, IXSTOP
                          I    = LINDX(IPNT)
                          Y(I) = Y(I) - T*LNZ(IX)
                          IPNT = IPNT + 1
                          if (mark(i) .eq. 0) then
                             mark(i) = 1
                             nnz = nnz + 1
                             nzidx(nnz) = i
                          endif
                       ENDDO
                    endif
                    IXSTRT = IXSTOP + 1
                    JPNT   = JPNT + 1
                 ENDDO
                 FJCOL = LJCOL + 1
              ENDDO
C
C             ------------------
C             DIAGONAL SOLVE ...
C             ------------------
c              DO  J = JCOL2, neqns
c                  W(J) = Y(J) / LNZ(XLNZ(J))
c              ENDDO
c
c               dinv(jcol2) = 0.0
c               do j = jcol2,neqns
c                  dinv(jcol2) = dinv(jcol2) + y(j)*y(j)/lnz(xlnz(j))
c                  nzidx(j) = -1
cc                  mark(j) = 0
c                  y(j) = 0.0d0
c               enddo 
               do j = 1, nnz
                  i = nzidx(j)
                  dinv(jcol2) = dinv(jcol2) + y(i)*y(i)/lnz(xlnz(i))
c                  nzidx(j) = -1
                  mark(i) = 0
                  Y(i)    = 0.0d0
               end do 
           ENDDO
        ENDDO
        DO JSUP = 1, NSUPER
           FJCOL = XSUPER(JSUP)
           LJCOL = XSUPER(JSUP+1)-1
           DO JCOL = FJCOL,LJCOL
              dinvp(PERM(JCOL)) = dinv(JCOL)
           ENDDO
        ENDDO
        deallocate(dinv)
C
        RETURN
      END
