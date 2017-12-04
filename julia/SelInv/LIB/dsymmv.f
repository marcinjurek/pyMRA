      SUBROUTINE DSYMMV( N, A, pointers, indices, X, Y )
c
c     Purpose:
c
c     This subroutine performs sparse symmetric matrix-vector
c     multiplication
c
c                 y <---- A*x
c
c     Arguements:
c
c     N          INTEGER (INPUT).
c                Dimension of the matrix
c
c     A          DOUBLE PRECISION array of size NNZ (INPUT).
c                Nonzero entries of a sparse matrix A.
c
c     POINTERS   INTEGER array of size N+1. (INPUT)
c                Column pointers of A.  POINTERS(j) gives the 
c                index of the starting point of column j in A.
c                NNZ = COLPTR(j+1)-1
c
c     INDICES    INTEGER array of size NNZ. (INPUT)
c                The row indices of each nonzero entry in A.
c
c     X          INTEGER array of size N. (INPUT)
c                The input vector.
c
c     Y          INTEGER array of size N. (OUTPUT)
c                The output vector.
*
*     .. Scalar Arguments .. 
      INTEGER            N
*     ..
*     .. Array Arguments .. 
      INTEGER            pointers( * ), indices( * )       
      DOUBLE PRECISION   A( * ) 
      DOUBLE PRECISION   X( * ), Y( * ) 
*
*==============================================================
*
*  does a sparse matrix * vector multiplication to compute 
*
*                   y = A * x  
*
*  where A is stored in the sparse column format
*
*==============================================================
*
*     .. Parameter .. 
*
*     .. Local Scalars .. 
      INTEGER            I, J, K, K1, K2
*
*     Initialization 
*
      DO 10 J = 1, N
	 Y( J ) = 0.0d0
  10  CONTINUE
*
*     Compute y = A*x
*
      DO 100 J = 1,N
         K1 = pointers( J )
         K2 = pointers( J+1 ) -1
         DO 110 K = K1, K2
            I = indices( K )
*
*           Compute y(i) = y(i) + a(i,j) * x(j)
*
            Y( I ) = Y( I ) + A( K )*X( J )
*    
	    if (i .ne. j) then
               y(j) = y(j) + a(k)*x(i)		
            endif 
*
 110     CONTINUE
 100  CONTINUE
*
*     End of DSYMMV
*
      END
