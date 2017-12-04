      subroutine readmatrixheader(filename, n, nnz)
      character*120 filename
      integer n, nnz
      integer funit

      funit = 22     
      open(unit=funit,file=filename,form='unformatted')
      read(funit)  n, nnz
      close(unit=funit)
      end subroutine 
