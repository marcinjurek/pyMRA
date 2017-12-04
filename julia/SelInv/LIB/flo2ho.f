      subroutine flo2ho(n,      colptr, rowind, nzvals, 
     &                  newptr, newind, newvals,ip)
      implicit none
c
c     Purpose:
c
c        This routine converts the lower triangluar part of symmetric
c        sparse matrix to its full representation (i.e. both
c        upper and lower triangular parts are stored.  The diagonal
c        of the matrix is INCLUDED in newvals in this version.
c
c     Arguments:
c
c     N        INTEGER (INPUT)
c              Dimension of the matrix.
c
c     COLPTR   INTEGER array of size N+1. (INPUT)
c              The column pointers of the lower triangular input matrix.
c
c     ROWIND   INTEGER array of size COLPTR(N+1)-1. (INPUT)
c              The row indices of the lower triangular input matrix.
c
c     NZVALS   DOUBLE PRECISION array of size COLPTR(N+1)-1. (INPUT)
c              The nonzero entries in the lower triangular part of the
c              matrix columnwise compressed.
c
c     NEWPTR   INTEGER array of size N+1. (OUTPUT)
c              The column pointers of the converted matrix.
c
c     NEWIND   INTEGER array of size H_NNZ. (OUTPUT)
c              The row indices of the converted matrix.
c
c     NEWVALS  DOUBLE PRECISION array of size H_NNZ (OUTPUT)
c              The nonzero entries of the upper and lower triangular
c              part of the matrix columnwize compressed.  The diagonal
c              of the matrix is not included.
c
c     IP       INTEGER array of size N. (WORK)
c              Work array for storing pointers.
c
      integer           n,           colptr(n+1), rowind(*), 
     &                  newptr(n+1), newind(*),   ip(n)
      double precision  nzvals(*),   newvals(*)
c
      integer           i,           j,           irow, 
     &                  nnz,         ibegin,      iend,
     &                  inew
c
      nnz = newptr(n+1)-1
c
c     === work array to accumulate the non-zero 
c         counts for each row ===
c
*vdir noaltcode
      do i = 1, n
         ip(i) = 0
      end do
c
c     === nonzero count for each row of the 
c         lower triangular part 
c         (including the diagonal) ===
c
c
      do j = 1, n
         ibegin = colptr(j)
         iend   = colptr(j+1)-1
*vdir nodep
         do i = ibegin, iend
            irow = rowind(i)
            ip(irow) = ip(irow) + 1
         end do
      end do
c
c     === nonzero count for each column 
c         (excluding the diagonal)      ===
c
      do j = 1, n
         ibegin = colptr(j)
         iend   = colptr(j+1)-1
         if (iend .gt. ibegin ) then
            ip(j) = ip(j) + iend - ibegin
         endif
      end do
c
c     === compute pointers to the beginning of each column ===
c
      newptr(1) = 1
      do i = 1, n
         newptr(i+1) = newptr(i)+ip(i)
      end do
c
*vdir noaltcode
      do i = 1, n
         ip(i) = newptr(i) 
      end do
c 
c     === copy the upper triangular part ===
c            (excluding the diagonal)
c
      do j = 1, n
         ibegin  = colptr(j)
         iend    = colptr(j+1)-1
         if (ibegin .lt. iend) then
*vdir nodep
            do i = ibegin+1, iend
               irow = rowind(i)
               newind(ip(irow))  = j
               newvals(ip(irow)) = nzvals(i)
               ip(irow) = ip(irow) + 1 
            end do
         end if 
      end do
c
c     === copy the lower triangular part ===
c            (including the diagonal)
c
      do j = 1, n
         ibegin = colptr(j)
         iend   = colptr(j+1)-1
         inew   = ip(j)
         if (ibegin .le. iend) then
            do i = ibegin, iend
               newind(inew)  = rowind(i)
               newvals(inew) = nzvals(i) 
               inew = inew  + 1
            end do
         endif
      end do
      return 
      end
