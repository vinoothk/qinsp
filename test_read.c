#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
double frobeniusNorm(const gsl_matrix *matrix, int size1, int size2)
{
    double result = 0.0;
    for(int i = 0; i < size1; ++i)
    {
        for(int j = 0; j < size2; ++j)
        {
            double value = (gsl_matrix_get(matrix,i,j));
            // printf("fffff norm %f",gsl_matrix_get(matrix,i,j));
            result += value * value;
            // return;
        }
    }
    // gsl_matrix_fprintf(stdout,matrix, "%g");
   //  printf("fffff norm %f",result);
    return sqrt(result);
}

int main()
{
//    char buffer[4096*] ;
   char *record;
   double output;
   double dnorm[500];
// *line;





     char * line = NULL;
    size_t len = 0;
    ssize_t read;
   int i=0,j=0;
//    int mat[100][100];
    gsl_matrix *a = gsl_matrix_alloc (500, 250);
    gsl_vector *v = gsl_vector_alloc(250);
   FILE *fstream = fopen("A.csv","r");
   if(fstream == NULL)
   {
      printf("\n file opening failed ");
      return -1 ;
   }
   while(((read = getline(&line, &len, fstream)) != -1))
    //    (line=fgets(buffer,sizeof(buffer),fstream))!=NULL)
   {
    // 
    j=0;
    
     record = strtok(line,",");
    //  printf("line : %s",line) ;
     
     while(record != NULL)
     {
   //  printf("record : %f",atof(record)) ;
    //here you can put the record into the array as per your requirement.
    //  a[i][j++] = atof(record) ;
    gsl_matrix_set (a, i, j,atof(record));
    //  printf("%d \t  %d\n",i,j);
      ++j;
     record = strtok(NULL,",");
       
     }
     ++i ;
   }
  
   output = frobeniusNorm(a,500,250);
   printf("frobenius norm %f",output);
   for (i=0; i<500;i++)
   {
      gsl_matrix_get_row(v,a,i);
      dnorm[i] = pow(gsl_blas_dnrm2(v), 2 );

      printf("%f",dnorm [i]);

   }
   // _gslblas_dnrm2(v);

   return 0 ;
 }