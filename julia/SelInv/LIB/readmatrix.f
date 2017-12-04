      subroutine readmatrix(filename, colptr, rowind, nzvals)
      character*120 filename
      integer n, nnz
      integer colptr(*), rowind(*)
      real*8 nzvals(*)
      integer funit

      funit = 22
      open(unit=funit,file=filename,form='unformatted')
      read(funit) n, nnz
      read(funit) (colptr(i), i = 1, n+1)
      read(funit) (rowind(i), i = 1, nnz)
      read(funit) (nzvals(i), i = 1, nnz)
      close(unit=funit)
      end subroutine 
