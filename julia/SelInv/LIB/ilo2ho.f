      subroutine ilo2ho(n, Hnnz, colptr, rowind, newptr, newind, ip)
      implicit none
c
c     Purpose:
c
c        This routine converts the lower triangluar representation 
c        of symmetric sparse matrix pattern to its full representation 
c        (i.e. both upper and lower triangular indices are stored.
c        Diagonal entries are EXCLUDED from NEWPTR and NEWIND.
c
c     Arguments:
c
c     N        INTEGER (INPUT)
c              Dimension of the matrix.
c
c     HNNZ     INTEGER (OUTPUT)
c              The number of nonzeros in the entire (both upper and lower
c              triangular part) of the matrix.
c 
c     COLPTR   INTEGER array of size N+1. (INPUT)
c              The column pointers of the lower triangular input matrix.
c
c     ROWIND   INTEGER array of size COLPTR(N+1)-1. (INPUT)
c              The row indices of the lower triangular input matrix.
c
c     NEWPTR   INTEGER array of size N+1. (OUTPUT)
c              The column pointers of the converted matrix.
c
c     NEWIND   INTEGER array of size H_NNZ. (OUTPUT)
c              The row indices of the converted matrix.
c
c     IP       INTEGER array of size N. (WORK)
c              Work array for pointers.
c
      integer           n,           Hnnz
      integer           colptr(n+1), rowind(*), newptr(n+1),
     &                  newind(*),   ip(n)
c
      integer           i,           j,         irow, 
     &                  nnz,         ibegin,    iend, 
     &                  inew
c
      nnz = colptr(n+1)-1
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
c         lower triangular part ===
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
c     === nonzero count for each column ===
c
      do j = 1, n
         ibegin = colptr(j)
         iend   = colptr(j+1)-1
         if (iend .ge. ibegin ) then
            ip(j) = ip(j) + iend - ibegin - 1 
         endif
      end do
c
c     === compute pointers to the beginning of each column ===
c
      newptr(1) = 1
      do i = 1, n
         newptr(i+1) = newptr(i)+ip(i)
      end do
      Hnnz = newptr(n+1)-1
c
*vdir noaltcode
      do i = 1, n
         ip(i) = newptr(i) 
      end do
c 
c     === copy the upper triangular part ===
c
      do j = 1, n
         ibegin = colptr(j)
         iend   = colptr(j+1)-1
         if (ibegin .lt. iend) then
*vdir nodep
            do i = ibegin+1, iend
               irow = rowind(i)
               newind(ip(irow)) = j
               ip(irow) = ip(irow) + 1 
            end do
         end if 
      end do
c
c     === copy the lower triangular part ===
c
      do j = 1, n
         ibegin = colptr(j)
         iend   = colptr(j+1)-1
         inew   = ip(j)
         if (ibegin .lt. iend) then
            do i = ibegin+1, iend
               newind(inew)  = rowind(i)
               inew = inew  + 1
            end do
         endif
      end do
      return 
      end
