      subroutine dumpLmat(nsuper, xsuper, xlindx, lindx, xlnz, 
     &                    lnz, filename)
      implicit none
      integer             nsuper
      integer             lindx(*) , xsuper(*)
      integer             xlindx(*), xlnz(*) 
      double precision    lnz(*)   
      character*120       filename

      integer fjcol, ljcol, ixstrt, ixstop, ipnt, jpnt, ix, i 
      integer jsup, jcol
  
      open(unit = 10, file=filename(1:8))
      fjcol = xsuper(1)
      do jsup = 1, nsuper
         ljcol = xsuper(jsup+1) - 1
         ixstrt = xlnz(fjcol)
         jpnt   = xlindx(jsup)
         do jcol = fjcol, ljcol
            ixstop    = xlnz(jcol+1) - 1
            ipnt      = jpnt
            do  ix = ixstrt, ixstop
                 i      = lindx(ipnt)
                 write(10,222) i,jcol,lnz(ix)
                 ipnt   = ipnt + 1
            enddo
            ixstrt = ixstop + 1
            jpnt   = jpnt + 1
         enddo
         fjcol = ljcol + 1
      enddo
 222  format('L(',I6,',',I6,')=',1pe22.15,';')
      write(10,*) 'd = diag(L);'
      write(10,*) 'n = size(L,1);'
      write(10,*) 'for j = 1:n'
      write(10,*) '   L(j,j) = 1.0;'
      write(10,*) 'end;'
      write(10,*) 'A = L*diag(d)*L'';'
      write(10,*) 'S = zeros(n);'
      write(10,*) 'S(n,n) = 1/d(n);'
      write(10,*) 'for j = n-1:-1:1'
      write(10,*) 'inz = j + find(L(j+1:n,j)~=0);'
      write(10,*) 'S(inz,j) = -S(inz,inz)*L(inz,j);'
      write(10,*) 'S(j,inz) =  S(inz,j)'';'
      write(10,*) 'S(j,j) = 1/d(j) - S(j,inz)*L(inz,j);'
      write(10,*) 'end;'
      write(10,*) 'di0 = diag(inv(A));'
      write(10,*) 'di = diag(S);'
      write(10,*) 'err = norm(di0-di)'
      close(10)
      end subroutine
