/**********************************************************************
***********************************************************************

   Version:        0.1
   Last modified:  October 27, 2009
   Authors:        
     Chao Yang,
     Computer Research Division, Lawrence Berkeley National Lab
     Lin Lin,
     Program in Applied and Computational Mathematics, Princeton
     University

***********************************************************************
***********************************************************************

     lap2d5pt.c is the driver routine that computes the selected
     inversion of a 2D Laplacian operator 
          -\partial^2/\partial x^2-\partial ^2/\partial y^2
     with Dirichlet boundary condition, discretized on a 2D square domain 
     [0, 1] * [0, 1] with grid size nx * ny.

     
     Usage: lap2d5pt.x -nx=<nx> -ny=<ny> -order=<-1|2|3> -chkerr=<0|1> -printa=<0|1>\n");
     -nx        :   The number of grid points in x direction. (default: nx=10)
     -ny        :   The number of grid points in y direction. (default: ny=10)
     -order     :   Reordering strategy:
		    order = -1 (default) : Multiple Minimum Degree Reordering.
		    If METIS is supported, then the following choices
		    are also available:
		    order = 2 : Node Nested Dissection
		    order = 3 : Edge Nested Dissection

     -chkerr    :   Whether full inversion is to be computed for
		    efficiency comparison.  
		    chkerr = 0 (default). Only perform selected inversion.
		    chkerr = 1. Both selected inversion and full
    		      inversion are computed.

     -printa    :   Output the sparse matrix A. 
                    printa = 0 (default). Do not output A.
		    printa = 1 . Output A.
***********************************************************************
**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include <unistd.h>
#include <string.h>

#define mesh(i,j) mesh[nx*((j)-1)+(i)-1]
#define nzvals(i) nzvals[(i)-1]
#define rowind(i) rowind[(i)-1]
#define colptr(i) colptr[(i)-1]
#define xsol(i)   xsol[(i)-1]
#define rhs(i)    rhs[(i)-1]
#define diag(i)   diag[(i)-1]
#define diag2(i)  diag2[(i)-1]

extern int ldlt_preprocess__(int *, int *, int *, int *, int *, int *, int *);
extern int ldlt_fact__(int *, int *, int *, double *);
extern int ldlt_solve__(int *, double *, double *);
extern int ldlt_free__(int *);
extern int ldlt_selinv__(int *, double *, int *);

#ifdef TIMING
extern double getime(void);
#endif

int main(int argc, char ** argv)
{
   int i, j, nnodes,nedges, ia, nnz, ibeg, iend ;
   int nx= -1, ny = -1, count, node;
   int *mesh, *perm;
   int *rowind, *colptr;
   double *nzvals, *rhs, *xsol, *diag, *diag2;
   int token, order=-1;
   int Lnnz;
   double t0,t1, errmax; 
   long myclktck;
   double dval, errabs;
   int chkerr = 0, printa = 0, dumpL=0;
   int ierr = 0;
   double hx, hy;  /* mesh size */




   ia = 1;
   while (ia < argc) {
      if ( !strncmp(argv[ia],"-nx",3) ) {
 	 nx = atoi(&argv[ia][4]);
      }
      else if ( !strncmp(argv[ia],"-ny",3) ) {
    	 ny = atoi(&argv[ia][4]);
      }
      else if ( !strncmp(argv[ia],"-order",6) ) {
    	 order = atoi(&argv[ia][7]);
      }
      else if ( !strncmp(argv[ia],"-chkerr",7) ) {
    	 chkerr = atoi(&argv[ia][8]);
         if (chkerr != 0) chkerr = 1;
      }
      else if ( !strncmp(argv[ia],"-printa",7) ) {
    	 printa = atoi(&argv[ia][8]);
         if (printa != 0) printa = 1;
      }
      else if ( !strncmp(argv[ia],"-dumpL",6) ) {
    	 dumpL = atoi(&argv[ia][7]);
         if (dumpL != 0) dumpL = 1;
      }
      else {
 	 fprintf(stderr, "invalide argument!\n");
	 fprintf(stderr, "Usage: lap2d5pt.x -nx=<nx> -ny=<ny> -order=<-1|2|3> -chkerr=<0|1> -printa=<0|1>\n");
         return 1;
      }
      ia++;
   }

   if (nx == -1 && ny > 0) nx = ny;
   if (ny == -1 && nx > 0) ny = nx;

   if (nx == -1) nx = 10;
   if (ny == -1) ny = 10;

   hx = 1.0 / (nx+1);
   hy = 1.0 / (ny+1);

   nnodes = nx*ny;
   mesh = (int*)malloc(nnodes*sizeof(int));

   for (i = 0; i<nnodes; i++) mesh[i]=i+1;
   
   /* first pass to count the number of edges */
   /* Dirichlet BC */ 
   nedges = 0; 
   for (j = 1; j <= ny; j++) {
     for (i = 1; i <= nx; i++) {
       if (j < ny) nedges++;
       if (j > 1)  nedges++;
       if (i < nx) nedges++;
       if (i > 1)  nedges++;   
     }
   }

   /* print the matrix dimension and number of nonzeros */
   nnz = nedges/2 + nnodes;
   printf("%d  %d\n", nnodes, nnz);

   colptr = (int*)malloc((nnodes+1)*sizeof(int));
   rowind = (int*)malloc(nnz*sizeof(int));
   nzvals = (double*)malloc(nnz*sizeof(double));
   colptr(1) = 1;
   count = 0;
   node  = 0;

   dval = 2.0/(hx*hx) + 2.0/(hy*hy);
   printf(" Dirichlet boundary condition\n");

   for (j = 1; j <= ny; j++) {
      for (i = 1; i <= nx; i++) {
	 /* diagonal */
         if (printa)
            printf("%d %d  %8.2e\n", mesh(i,j), mesh(i,j), dval); 

         rowind[count] = mesh(i,j);
         nzvals[count] = dval;
         count++;

         /* lower */
         if (i < nx) {
            if (printa) 
	       printf("%d %d -1.0\n", mesh(i+1,j), mesh(i,j));

            rowind[count] = mesh(i+1,j);
            nzvals[count] = -1.0/(hx*hx);
            count++;
         }


         /* right */
         if (j < ny) {
            if (printa) 
	       printf("%d %d -1.0\n", mesh(i,j+1), mesh(i,j)); 

            rowind[count] = mesh(i,j+1);
            nzvals[count] = -1.0/(hy*hy);
            count++;
         }        

         node++;
         colptr(node+1) = count+1; 
      } 
   }
   if (count != nnz) {
       printf(" count = %d, nnz = %d\n", count, nnz);  
       return 1;
   }

   token = 0;
  

   ldlt_preprocess__(&token, &nnodes, colptr, rowind, &Lnnz, &order, perm);   
   ldlt_fact__(&token, colptr, rowind, nzvals);

   rhs = (double*)calloc(nnodes,sizeof(double));
   xsol = (double*)calloc(nnodes,sizeof(double));
   diag = (double*)malloc(nnodes*sizeof(double));

   if (chkerr) {
#ifdef TIMING
      t0 = getime();
#endif
/*
      for (i=1; i<=nnodes; i++) {
         rhs(i) = 1.0;
         ldlt_solve__(&token, xsol, rhs);
         diag(i) = xsol(i);
         rhs(i) = 0.0;
      } 
*/
      ldlt_diaginv__(&token, diag); 
      /*for (i=1; i<=nnodes;i++) printf("diag = %15.7e\n", diag(i));*/
#ifdef TIMING
      t1 = getime();
      printf(" TIME FOR FULL INVERSION = %11.3e\n", t1-t0);
#endif
   }

   /* selected inversion */
   diag2 = (double*)calloc(nnodes,sizeof(double));
   ldlt_selinv__(&token, diag2, &dumpL);
     /* for (i=1; i<=nnodes;i++) printf("diag2 = %15.7e\n", diag2(i)); */

   if (chkerr) {
      errmax = 0.0;
      for (i=1; i<=nnodes; i++) {
          errabs = fabs(diag(i)-diag2(i));
          /* printf("err(%d) = %11.3e\n", i, errabs); */
          if ( errabs > errmax )
              errmax = fabs(diag(i)-diag2(i)); 
      }
      printf(" ERROR OF THE DIAGONAL ELEMENTS BETWEEN FULL INVERSION AND SELECTED INVERSION = %11.3e\n", errmax);
   }

   ldlt_free__(&token); 

   if (order == 0) free(perm);
   free(mesh);
   free(colptr); 
   free(rowind); 
   free(nzvals); 
   free(rhs); 
   free(xsol); 
   free(diag); 
   free(diag2); 
}
