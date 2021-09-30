#include <stdio.h> 

int main()
{
  const int N=100000;
  const int dim=100;
  double a [dim][dim];
  double b [dim][dim];
  
  for (int n=0; n<N; n++)
  {
    for  (int i=0; i<dim; i++)
    {
      for (int j=0; j<dim; j++)
      {
    	b[i][j] = a[i][j];
      }
    }
  }
  printf("a[0] = %lf ",a[0][0]);
  return 0;
}