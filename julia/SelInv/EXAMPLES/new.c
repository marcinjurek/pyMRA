#include <stdlib.h>

void alloc( int n, int * vec)
{ 
  vec=(int *)malloc(n*sizeof(int));
  int i;
  for (i=0;i<n;i++) 
      vec[i]=1;
}
