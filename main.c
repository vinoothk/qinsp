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
#include <gsl/gsl_linalg.h>


double frobeniusNorm(const gsl_matrix *matrix, int size1, int size2)
{
    double result = 0.0;
    for(int i = 0; i < size1; ++i)
    {
        for(int j = 0; j < size2; ++j)
        {
            double value = (gsl_matrix_get(matrix,i,j));
            result += value * value;
        }
    }
    printf("fffff norm %f",result);
    return sqrt(result);
}

// void sample_me_lsyst(A, b, m, n, samples, rank, r, w, rows, sigma, row_norms, LS_prob_rows, LS_prob_columns, A_Frobenius)
// {


// }                                            
void linear_eqs(const gsl_matrix *A, const gsl_vector *b, int m, int n, int r, int c, int rank, int Nsamples, int NcompX, double row_norms[], gsl_vector *LS_prob_rows, const gsl_matrix *LS_prob_columns,double A_Frobenius)
{
  // m_rows, n_cols = np.shape(A)
    int  *rows = malloc(sizeof(int)*r);
    gsl_vector *sigma = gsl_vector_alloc(c); 
    gsl_matrix *vh = gsl_matrix_alloc (c, c);
    gsl_matrix *w = gsl_matrix_alloc(r,c);
    gsl_vector *u_approx = gsl_vector_alloc(m);
    gsl_vector *v_approx = gsl_vector_alloc(n);
    double *lambdas = (int *)calloc(rank, sizeof(double));
    clock_t start,end;
    double rt_ls_prob;
    
    # 1- Generating LS probability distributions to sample from matrix A
    
    start = clock();
    LS = ls_probs(m, n, A);
    end = clock();
    

    rt_ls_prob = ((double) (end - start)) / CLOCKS_PER_SEC;

    # 2- Building matrix C by sampling "r" rows and "c" columns from matrix A and computing SVD of matrix C
    
    sample_C(A, m, n, r, c, row_norms, LS_prob_rows, LS_prob_columns, A_Frobenius,w,vh,sigma);
    
    ul_approx = gsl_matrix_alloc (m, rank);
    vl_approx = gsl_matrix_alloc (n, rank);
    gsl_matrix_Set_zero(ul_approx);
    gsl_matrix_set_zero(vl_approx);
    // for l in range(rank):
    //     ul_approx[:, l], vl_approx[:, l] = 
    for (int l=0;l<rank;l++)
    {
      uvl_vector(l, A, r, w, rows, sigma, row_norms, A_Frobenius);
    }
    start = clock();
    sample_me_lsyst(A, b, m, n, Nsamples, rank, r, w, rows,sigma,row_norms,LS_prob_rows,LS_prob_columns,A_Frobenius,lambdas);
    end = clock();
    rt_sampling_me = ((double) (end - start)) / CLOCKS_PER_SEC;

    w_vector = np.zeros(r)
    for l in range(rank):
        w_vector[:] += (lambdas[l] / sigma[l] ** 3) * w[:, l]

    w_norm = la.norm(w_vector)

    gsl_vector *w_vector = gsl_vector_alloc(r);
    gsl_vector *vtmp2 = gsl_vector_alloc();
    for( int l=0; l<rank; l++)
    {
      pow(lambdas[l] / sigma[l], 3) * gsl_matrix_get_column(w,)
      w_vector
    }


}
void sample_me_lsyst(const gsl_matrix *A, const gsl_vector* b, int m, int n, int samples, int rank, int r, gsl_vector *w, int *rows, gsl_vector *sigma, double row_norms[], gsl_vector *LS_prob_rows, const gsl_matrix *LS_prob_columns, double A_Frobenius, double *lambdas)
{
  int reps = 10;
  int 
  if(reps % 2 == 0)
  {
    m = reps / 2;
  }
  else
  {
    m = reps / 2;
  }
  gsl_matrix *matrix_elements = gsl_matrix_alloc(reps,rank);

    for(int i=0; i<reps; i++)
    {
      for(int l=0; l<rank; l++)
      {
        gsl_vector *X = gsl_vector_alloc(samples);
        
        gsl_vector_set_zero(X);

        for(int k=0; k<samples; 
        {
          int *sample_i;
          gsl_vector *vtmp = gsl_vector_alloc(n);

            sample_i = vose(LS_prob_rows,m,1);
            
            gsl_matrix_get_row(vtmp,LS_prob_columns,*sample_i);
            
            sample_j = vose(LS_prob_columns,n,1);

            # calculates v_j
            
            int v_j = 0; 
            
            for (int s=0; s<r; s++)
            {
              v_j+= gsl_matrix_get(rows[s],sample_j) * gsl_matrix_get(w,s,l) / (sqrt(row_norms[rows[s]]);

            }

            v_j = v_j * A_Frobenius / (sqrt(r) * gsl_vector_get(sigma,l));
            
            X[k] = ((A_Frobenius ** 2 * b[sample_i]) / (A[sample_i, sample_j])) * v_j
            
            double tmp = pow(A_Frobenius,2) * gsl_vector_get(b,sample_i) / gsl_matrix_get(sample_i,sample_j) * v_j;
            
            gsl_vector_set(X,k,tmp);
        }

        gsl_matrix_set(matrix_elements,i,l,(gsl_vector_sum(X)/samples));

      }
    }

    # take median of all repeated estimates    
    double tmp1;
    gsl_vector *vtmp = gsl_vector_alloc(reps);
    for (int l=0; l<rank; l++)
    {
      gsl_stats_median(lambdas,1,rank); 
      gsl_matrix_get_column(vtmp,matrix_elements,l);
      gsl_sort_vector(vtmp);
      if(m==0)
        lambdas[l] = gsl_vector_get(vtmp,reps/2);
      else 
        lambdas[l] = (gsl_vector_get(vtmp,reps+1/2) + gsl_vector_get(vtmp,reps-1/2))/2;
    }

}




void uvl_vector(int m, int n, int l, const gsl_matrix* A, int r, gsl_vector *w, int *rows, gsl_vector *sigma, double row_norms[], A_Frobenius, gsl_vector *u_approx, gsl_vector *v_approx)
{
    
    # building approximated v^l vector
    gsl_vector *tmp_vector = gsl_vector_alloc(n);
    double tmp_val;
    factor = A_Frobenius / ( np.sqrt(r) * gsl_vector_get(sigma,l) )
    for s in range(r):
        v_approx[:] += ( A[rows[s], :] / np.sqrt(row_norms[rows[s]]) ) * w[s, l]

    for (int s=0; s<r; s++)
    {

      gsl_matrix_get_row(tmp_vector,A,rows[s]);
      tmp_val = np.sqrt(row_norms[rows[s]]) ) * gsl_matrix_get(w,s,l);
      gsl_vector_scale(tmp_vector,tmp_val);
      gsl_vector_add(v_approx,tmp_vector);

    }
    gsl_vector_scale(v_approx,factor);
    
    dgemm('N', 'N', 1, A, v_approx, 0.0, u_approx);

    gsl_vector_scale(u_approx, 1/gsl_vector_get(sigma,l));


}


void sample_C(const gsl_matrix *A, int m, int n,int r,int c, double row_norms[], gsl_vector *LS_prob_rows, const gsl_matrix *LS_prob_columns,double A_Frobenius,gsl_matrix *C,gsl_matrix *vh,gsl_vector* sigma)
{
    // int  *rows = malloc(sizeof(int)*r); 
    int a[1];
    int *columns_tmp =  malloc(sizeof(int)*1);
    int columns[c];
    clock_t start, end;
    double rt_sampling_C,rt_building_C,rt_svd_C;
    double tmp;
    gsl_vector *v = gsl_vector_alloc(n);
    gsl_vector *R_row = gsl_vector_alloc(n);
    gsl_matrix *R_C = gsl_matrix_alloc (r, c);
    // gsl_matrix *C = gsl_matrix_alloc(r,c);
    gsl_vector *vtmp = gsl_vector_alloc(c);
    gsl_vector *column_norms = gsl_vector_alloc(r);
    gsl_vector *vtmp2 = gsl_vector_alloc(r);
    double R_row_norm;
    gsl_matrix *LS_prob_columns_R = gsl_matrix_alloc (r, n);
    // gsl_matrix *Vh = gsl_matrix_alloc (c, c);
    gsl_rng * rg = gsl_rng_alloc (gsl_rng_taus); 
    // gsl_vector *S = gsl_vector_alloc(c); 
    // gsl_vector *work = gsl_vector_alloc(c);
    start = clock();
    rows = vose(LS_prob_rows,m,r);
    
    for(int i=0;i<c;i++)
    {
    
     gsl_ran_choose (rg, a, 1, rows,r, sizeof (int));
    
     gsl_matrix_get_row(v,LS_prob_columns,a[0]);
    
     columns_tmp = vose(v,n,1);
    
     columns[i] = columns_tmp[0];
    }

    end = clock();
    
    rt_sampling_C = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    printf("time elapsed %f",rt_sampling_C );

    start = clock();
    
    for(int s=0; s < r;s++)
    {
      
        tmp = A_Frobenius/ sqrt(r) * sqrt(row_norms[rows[s]]);
        
        gsl_matrix_get_row(R_row,A,rows[s]);
        
        gsl_vector_scale(R_row,tmp);
      
        
        R_row_norm = pow(gsl_blas_dnrm2(R_row),2);
        // printf("%f\n",R_row_norm);
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
        
        gsl_matrix_get_row(vtmp,R_C,s);

        gsl_vector_scale(vtmp,tmp);
        
        gsl_matrix_set_row(R_C,s,vtmp);
      
    }
    gsl_vector_set_zero(column_norms);

    for(int t=0;t<c;t++)
    {
      for(int s=0;s<r;s++)
      {
        pow(gsl_matrix_get(R_C,t,s),2);
        
        gsl_vector_set(column_norms, t, gsl_vector_get(column_norms,t) + pow(gsl_matrix_get(R_C,t,s),2));                                                            
      }
    }

    for(int t=0; t<c; t++)
    {
      gsl_matrix_get_col(vtmp2,C,t);
      
      gsl_matrix_get_col(vtmp,R_C,t);
      
      tmp = (A_Frobenius )/ sqrt(gsl_vector_get(column_norms,t)) / sqrt(c);
      
      gsl_vector_scale(vtmp,tmp);
      
      gsl_matrix_set_col(C,t,vtmp);
    }
    end = clock();
    
    rt_building_C = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    printf("rt_building_C %f",rt_building_C );
  
    start = clock();
    
    gsl_linalg_SV_decomp(C,vh,sigma,work);
    
    end = clock();
    
    rt_svd_C = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    printf("rt_svd_C %f",rt_svd_C );
    // for(int i =0;i<r;i++)
    // {
    //   for(int j=0;j<c;j++)
    //   {
    //     printf("%f\t,",gsl_matrix_get(Vh,i,j));
    //   }
    //   printf("\n");
    // }
    // for(int j=0;j<c;j++)
    //   {
    //     printf("%f\n,",gsl_vector_get(S,j));
    //   }
    // sample_me_lsyst(A, b, m, n, samples, rank, r, w, rows, sigma, row_norms, LS_prob_rows, LS_prob_columns, A_Frobenius);
}
int main()
{
   int m = 500, n =250;//610, 9724
   int r = 400,c = 100;
   int Nsamples =50,NcompX = 50;
   int rank = 3;
   char *record;

   char *recordb;
   double A_Frobenius;
   double row_norms[m];
   double LS_prob_row_tmp;
   double LS_prob_column_tmp;
    char * line = NULL;
    char *lineb = NULL;
    size_t len = 0;
    size_t lenb = 0;
    size_t read;
    size_t readb;
    int i=0,j=0;
    gsl_matrix *a = gsl_matrix_alloc (m, n);
    gsl_matrix *LS_prob_columns = gsl_matrix_alloc (m, n);
    gsl_vector *v = gsl_vector_alloc(n);
    gsl_vector *v = gsl_vector_alloc(b);
    gsl_vector *LS_prob_rows = gsl_vector_alloc(m);
   FILE *fstream = fopen("Data/A.csv","r");
   FILE *fstreamb = fopen("Data/B.csv","r");
   srand(time(NULL));
   if(fstream == NULL)
   {
      printf("\n file-1 opening failed  ");
      return -1 ;
   }
   if(fstreamb == NULL)
   {
      printf("\n file-2 opening failed  ");
      return -1 ;
   }
   while(((read = getline(&line, &len, fstream)) != -1))
   {
    // 
    j=0;
     record = strtok(line,",");
         
     while(record != NULL)
     {
    gsl_matrix_set (a, i, j,atof(record));
      ++j;
     record = strtok(NULL,",");
       
     }
     ++i ;
   }
   int ib = 0;

    while(((readb = getline(&lineb, &lenb, fstreamb)) != -1))
   {
    gsl_vector_set(b,atof(lineb));
   }


   A_Frobenius = frobeniusNorm(a,m,n);
   printf("frobenius norm %f",A_Frobenius);
   for (i=0; i<m;i++)
   {
      gsl_matrix_get_row(v,a,i);
      row_norms[i] = pow(gsl_blas_dnrm2(v), 2 );
   }
    for(i=0;i<m;i++)
    {

      LS_prob_row_tmp = row_norms[i] / pow(A_Frobenius, 2);
      gsl_vector_set(LS_prob_rows,i,LS_prob_row_tmp);
    }
    for(int i=0;i<m;i++)
    {
      for (int j=0; j<n; j++)
      {
        LS_prob_column_tmp = (pow(gsl_matrix_get(a,i,j),2) / row_norms[i] );
        gsl_matrix_set(LS_prob_columns,i,j,LS_prob_column_tmp);
      }
    }
    // sample_C(a,m,n,r,c,row_norms,LS_prob_rows,LS_prob_columns,A_Frobenius);
    linear_eqs(A, b, r, c, rank, Nsamples, NcompX,row_norms, LS_prob_rows, LS_prob_columns, A_Frobenius);
    gsl_matrix_free(a);
   return 0 ;
 }
