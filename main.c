#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include "vose_gsl.c"
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
            // printf("al %f\n",value);
            result += value * value;
            // return;
        }
    }
    // gsl_matrix_fprintf(stdout,matrix, "%g");
    printf("fffff norm %f",result);
    return sqrt(result);
}
void sample_C(const gsl_matrix *A, int m, int n,int r,int c, double row_norms[], gsl_vector *LS_prob_rows, const gsl_matrix *LS_prob_columns,double A_Frobenius)
{
    int  *rows = malloc(sizeof(int)*r);
    // double *R_row = malloc(sizeof(double)*n);
    int a[1];
    int *columns_tmp =  malloc(sizeof(int)*1);
    int columns[c];
    clock_t start, end;
    double rt_sampling_C;
    double tmp;
    gsl_vector *v = gsl_vector_alloc(n);
    gsl_vector *R_row = gsl_vector_alloc(n);
    gsl_matrix *R_C = gsl_matrix_alloc (r, c);
    gsl_vector *vtmp = gsl_vector_alloc(c);
    gsl_vector *column_norms = gsl_vector_alloc(c);
    
    double R_row_norm;
    gsl_matrix *LS_prob_columns_R = gsl_matrix_alloc (r, n);
    // const gsl_rng_type * T;
    // gsl_rng * rg;
    gsl_rng * rg = gsl_rng_alloc (gsl_rng_taus);  
    
    start = clock();
    rows = vose(LS_prob_rows,m,r);
    // for(int s=0; s < r;s++)
    // {
    //   printf(" %d\n",rows[s]);
    // }
    // exit(0);
    
    // end = clock();
    // rt_sampling_C = ((double) (end - start)) / CLOCKS_PER_SEC;
    for(int i=0;i<c;i++)
    {
     gsl_ran_choose (rg, a, 1, rows,r, sizeof (int));
     gsl_matrix_get_row(v,LS_prob_columns,a[0]);
     columns_tmp = vose(v,n,1);
    //  printf("%d aaaaaaa-\t -- \n",a[0]);
        // printf("%f\n",LS_prob_columns[a[0]][i]); 
     columns[i] = columns_tmp[0];
    //  printf("%d\n",columns[i]);
    }
    // exit(0);
    end = clock();
    rt_sampling_C = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("time elapsed %f",rt_sampling_C );
   
    for(int s=0; s < r;s++)
    {
      
      tmp = A_Frobenius/ sqrt(r) * sqrt(row_norms[rows[s]]);
      gsl_matrix_get_row(R_row,A,rows[s]);
      gsl_vector_scale(R_row,tmp);
      
      
      R_row_norm = pow(gsl_blas_dnrm2(R_row),2);
      printf("%f\n",R_row_norm);
      gsl_vector_mul(R_row,R_row);
      gsl_vector_scale(R_row,1/R_row_norm);
      gsl_matrix_set_row(LS_prob_columns_R,s,R_row);
    }
    for(int s=0;s<r;s++)
    {
      for(int t=0;t<c;t++)
      {
        gsl_matrix_set(R_C,s,t,gsl_matrix_get(A,rows[s],columns[t]));
        
      }
      tmp = A_Frobenius/ sqrt(r) * sqrt(row_norms[rows[s]]);
      // R_C[s,:] = R_C[s,:] * s
      gsl_matrix_get_row(vtmp,R_C,s);

      gsl_vector_scale(vtmp,tmp);
      gsl_matrix_set_row(R_C,s,vtmp);
      
    }
    // for(int s=0;s<r;s++)
    // {
    //   for(int t=0;t<c;t++)
    //   {
    //     printf("%f\n",gsl_matrix_get(R_C,s,t));
    //    }
    //   }


    gsl_vector_set_zero(column_norms);
    for(int t=0;t<c;t++)
    {
      for(int s=0;s<r;s++)
      {
        // column_norms[t] += np.abs(R_C[s, t])**2
        pow(gsl_matrix_get(R_C,t,s),2);
        gsl_vector_set(column_norms,t,gsl_vector_get(column_norms,t) +pow(gsl_matrix_get(R_C,t,s),2));                                                            
      }

    }
    //   for(int s=0; s < c;s++)
    // {
    //   printf(" column_norms%f\n",gsl_vector_get(column_norms,s));
    // }
}
int main()
{
//    char buffer[4096*] ;
   int m = 472, n =472;//610, 9724
   int r = 400,c = 400;
   char *record;
   double A_Frobenius;
   double row_norms[m];
// *line;
   double LS_prob_row_tmp;
   double LS_prob_column_tmp;
  //  double LS_prob_columns[m][n];
    // double **LS_prob_columns =  (int **)malloc(m * sizeof(double *));
    // for (int i=0; i<m; i++)
        //  LS_prob_columns[i] = (int *)malloc(n * sizeof(double));
    char * line = NULL;
    size_t len = 0;
    size_t read;
    int i=0,j=0;
//    int mat[100][100];
    gsl_matrix *a = gsl_matrix_alloc (m, n);
    gsl_matrix *LS_prob_columns = gsl_matrix_alloc (m, n);
    gsl_vector *v = gsl_vector_alloc(n);
    gsl_vector *LS_prob_rows = gsl_vector_alloc(m);
   FILE *fstream = fopen("Data/Snp.csv","r");
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
    // printf("record : %f",atof(record)) ;
    //here you can put the record into the array as per your requirement.
    //  a[i][j++] = atof(record) ;
    gsl_matrix_set (a, i, j,atof(record));
    // printf("%f \t  \n",record);
      ++j;
     record = strtok(NULL,",");
       
     }
     ++i ;
   }
  
   A_Frobenius = frobeniusNorm(a,m,n);
   printf("frobenius norm %f",A_Frobenius);
   for (i=0; i<m;i++)
   {
      gsl_matrix_get_row(v,a,i);
      row_norms[i] = pow(gsl_blas_dnrm2(v), 2 );

      // printf("%f\n",row_norms[i]);

   }
  // exit(0);
  //normalized length-square row probability distribution

    for(i=0;i<m;i++)
    {

      LS_prob_row_tmp = row_norms[i] / pow(A_Frobenius, 2);
      gsl_vector_set(LS_prob_rows,i,LS_prob_row_tmp);
      // printf("%f\n",gsl_vector_get(LS_prob_rows,i));
      
    }
    // exit(0);
    for(int i=0;i<m;i++)
    {
      for (int j=0; j<n; j++)
      {
        // printf("%f\n",gsl_matrix_get(a,i,j));
        // LS_prob_columns[i][j] 
        LS_prob_column_tmp = (pow(gsl_matrix_get(a,i,j),2) / row_norms[i] );
        gsl_matrix_set(LS_prob_columns,i,j,LS_prob_column_tmp);
        // gsl_matrix_get()
        // printf("%f\n",gsl_matrix_get(a,i,j));
        // printf("LS_prob_columns%f \n",gsl_matrix_get(LS_prob_columns,i,j));
      }
      // printf("\n");
    }
    sample_C(a,m,n,r,c,row_norms,LS_prob_rows,LS_prob_columns,A_Frobenius);
    gsl_matrix_free(a);
   return 0 ;
 }
