#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include "vose_algo.c"
#include <time.h>
#include<gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>


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
void sample_C(const gsl_matrix *matrix, int m, int n,int r,int c, double row_norms[], double LS_prob_rows[], double LS_prob_columns[500][250],double A_Frobenius)
{
  int  *rows = malloc(sizeof(int)*r);
  int a[1];int rng_choose;
  int *columns_tmp =  malloc(sizeof(int)*1);
  int columns[c];
  clock_t start, end;
  double rt_sampling_C;
  // const gsl_rng_type * T;
  // gsl_rng * rg;
  gsl_rng * rg = gsl_rng_alloc (gsl_rng_taus);
  for(int i = 0; i < m; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            //Testing
            double value = (gsl_matrix_get(matrix,i,j));
            // printf("%f \n",gsl_matrix_get(matrix,i,j));
            // printf("Vinooth");
            // printf(" %f\n ", LS_prob_columns[i][j]);

            // result += value * value;
            // return;
        }
    }  
    // for(int i=0;i<m;i++)
    // {
    // Testing
    //   // LS_prob_rows[i] = row_norms[i] / pow(A_Frobenius, 2);
    //   printf("%f\n",LS_prob_rows[i]);
      
    // }
    start = clock();
    rows = vose(LS_prob_rows,m,r);

    for(int i=0;i<c;i++)
    {
     gsl_ran_choose (rg, a, 1, rows,m, sizeof (int));
    //  printf("%d\n",a[0]);
    //  for(j=0;j<n;j++)
    //  {
    //    LS_prob_columns_tmp[j] = 
    //  }
     columns_tmp = vose(LS_prob_columns[a[0]],n,1);
     columns[i] = columns_tmp[0];
    //  printf("%d\n",columns[i]);
    }
    end = clock();
    rt_sampling_C = ((double) (end - start)) / CLOCKS_PER_SEC;
  // for(int j=0;j<c;j++)
  //    {
  //      printf("columns %d \t s%f\n",columns[j],LS_prob_columns[a[0]][j]);
  //    }
    
  printf("time elapsed %f",rt_sampling_C );
}
int main()
{
//    char buffer[4096*] ;
   int m = 500, n =250;
   int r = 200,c = 200;
   char *record;
   double A_Frobenius;
   double row_norms[m];
// *line;
   double LS_prob_rows[m];
   double LS_prob_columns[m][n];
    char * line = NULL;
    size_t len = 0;
    size_t read;
    int i=0,j=0;
//    int mat[100][100];
    gsl_matrix *a = gsl_matrix_alloc (m, n);
    gsl_vector *v = gsl_vector_alloc(n);
   FILE *fstream = fopen("A.csv","r");
   srand(time(NULL));
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
  
   A_Frobenius = frobeniusNorm(a,m,n);
  //  printf("frobenius norm %f",A_Frobenius);
   for (i=0; i<m;i++)
   {
      gsl_matrix_get_row(v,a,i);
      row_norms[i] = pow(gsl_blas_dnrm2(v), 2 );

      // printf("%f\n",row_no rms[i]);

   }

  //normalized length-square row probability distribution

    for(i=0;i<m;i++)
    {

      LS_prob_rows[i] = row_norms[i] / pow(A_Frobenius, 2);
      // printf("%f\n",LS_prob_rows[i]);
      
    }

    for(i=0;i<m;i++)
    {
      for (j=0; j<n; j++)
      {

        LS_prob_columns[i][j] = (pow(gsl_matrix_get(a,i,j),2) / row_norms[i] );
        // printf("%f \n",LS_prob_columns[i][j]);
      }
      // printf("\n");
    }
    sample_C(a,m,n,r,c,row_norms,LS_prob_rows,LS_prob_columns,A_Frobenius);
   return 0 ;
 }
