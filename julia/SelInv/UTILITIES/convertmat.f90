!***********************************************************************
!***********************************************************************
!
!   Version:        0.1
!   Last modified:  October 27, 2009
!   Authors:        
!     Chao Yang,
!     Computer Research Division, Lawrence Berkeley National Lab
!     Lin Lin,
!     Program in Applied and Computational Mathematics, Princeton
!     University
!
!***********************************************************************
!***********************************************************************
!
!       ****************************************************************
!       ****************************************************************
!
!       CONVERTMAT.F converts the Matrix Market (MM) format to the 
!       ccf format that selinv.x uses. This subroutine utilizes mmio.f 
!       provided by Matrix Market.
!
!       usage: convertmat.x -infile=<input file name> -outfile=<output file name>
!
!       For example, the MM format of Harwell-Boeing matrix bcsstk14.mtx 
!       is included in the currrent folder. In order to convert .mtx format 
!       to .ccf format, run
!
!       convertmat.x -infile=bcsstk14.mtx -outfile=bcsstk14.ccf

program convertmat
implicit none
integer:: iunit, ounit, rows, cols, entries
character*10 :: rep
character*7 ::field
character*19::  symm
integer:: nnzmax = 40000
integer, allocatable, dimension(:):: indx, jndx
double precision, allocatable, dimension(:) :: rval
complex, allocatable, dimension(:):: cval
integer, allocatable, dimension(:):: ival
integer:: i, j, k

integer, allocatable, dimension(:) :: colptr
character*120:: infilename, outfilename, str
integer:: nargc
! integer, external:: iargc
! NOTE: On some machines iargc should NOT be declared as external, but as
! intrinsic. So if there is error here, try
 integer, intrinsic:: iargc

nargc = iargc()
infilename = ''
outfilename = ''
if( nargc == 0 ) then
  write(0, *) 'invalid argument'
  write(0, *) 'usage convertmat.x -infile=<input file name> -outfile=<output file name>'
  write(0, *) 
  write(0, *) 'For example, the MM format of Harwell-Boeing matrix bcsstk14.mtx'
  write(0, *) 'is included in the currrent folder. In order to convert .mtx format'
  write(0, *) 'to .ccf format, run'
  write(0, *)
  write(0, *) 'convertmat.x -infile=bcsstk14.mtx -outfile=bcsstk14.ccf'
  stop
end if
do i = 1, nargc
  call getarg(i, str)
  if( str(1:8) .eq. '-infile=' ) then
    infilename = str(9:)
  elseif( str(1:9) .eq. '-outfile=' ) then
    outfilename = str(10:)
  else
    write(0, *) 'invalid argument'
    write(0, *) 'usage convertmat.x -infile=<input file name> -outfile=<output file name>'
    write(0, *) 
    write(0, *) 'For example, the MM format of Harwell-Boeing matrix bcsstk14.mtx'
    write(0, *) 'is included in the currrent folder. In order to convert .mtx format'
    write(0, *) 'to .ccf format, run'
    write(0, *)
    write(0, *) 'convertmat.x -infile=bcsstk14.mtx -outfile=bcsstk14.ccf'
    stop
  end if
end do
! ===================================================================
! Read input data in Matrix Market (MM) format
! ===================================================================
iunit = 20
open(unit=iunit, file=infilename,action='read')
call mminfo(iunit, rep, field, symm, rows, cols, entries)
write(6, *) rep
write(6, *)field
write(6, *) symm
write(6, *)rows, cols, entries

nnzmax = entries
allocate(indx(nnzmax), jndx(nnzmax), rval(nnzmax))

call mmread(iunit, rep, field, symm, rows, cols, entries, nnzmax, &
  indx, jndx, ival, rval, cval)
close(iunit)

! ===================================================================
! Generate compressed column format
! ===================================================================

allocate( colptr(cols+1) )
colptr(1) = 1
j = 1
do k = 1, entries
  if( jndx(k) == j+1 ) then
    j = j+1
    colptr(j) = k
  elseif( jndx(k) == j) then
    ! do nothing
  else
    write(6, *) "the matrix is not correctly ordered or contain empty column"
    stop
  end if
end do
colptr(cols+1) = entries + 1

! ===================================================================
! Output in compressed column format
! ===================================================================

ounit = 21
open(unit=ounit, file=outfilename,action='write', form='unformatted')
write(ounit) cols, entries
write(ounit) (colptr(k), k = 1, cols+1)
write(ounit) (indx(k), k = 1, entries)
write(ounit) (rval(k), k = 1, entries)
close(ounit)

deallocate(indx,jndx,rval, colptr)
stop
end program

