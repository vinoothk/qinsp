#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include "vose_gsl.c"
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_vector.h>

double frobeniusNorm(const gsl_matrix *matrix, int size1, int size2)
{ 
  double result = 0.0;
  for (int i = 0; i < size1; ++i)
  {
    for (int j = 0; j < size2; ++j)
    {
      double value = (gsl_matrix_get(matrix, i, j));
      result += value * value;
    }
  }
  printf("fffff norm %f", result);
  return sqrt(result);
}

// void sample_me_lsyst(A, b, m, n, samples, rank, r, w, rows, sigma, row_norms, LS_prob_rows, LS_prob_columns, A_Frobenius)
// {

// }


double approx_solution(const gsl_matrix *A, int rank, int r, gsl_matrix *w, int* rows,gsl_vector *sigma,double row_norms[], double A_Frobenius,double* lambdas, int comp)
{
  
    double approx_value = 0;
    
    for (int l=0; l< rank; l++)
    {
       double v_comp=0;

       for (int s=0; s<r; s++)
       {
          v_comp +=  gsl_matrix_get(A,rows[s],comp) * gsl_matrix_get(w,s,l) / sqrt(row_norms[rows[s]]);

       }
       v_comp = v_comp * A_Frobenius / sqrt(r) * gsl_vector_get(sigma,l);

       approx_value += v_comp * lambdas[l] * pow(gsl_vector_get(sigma,l), 2);

    }
    // printf("Approximate solution %f\n", approx_value);
    return approx_value;

}

int* sample_from_x(const gsl_matrix *A, int r, int c, int n, int *rows, double row_norms[], gsl_matrix *LS_prob_columns_R, double A_Frobenius, gsl_vector *w_vector, double w_norm)
{
    int keep_going = 1;
    int out_j = 0;
    int counter = 0;
    int *out_j_conter = malloc(sizeof(int) * 2);
    int i_sample[1];
    gsl_rng *rg = gsl_rng_alloc(gsl_rng_taus);
    // gsl_rng *rg1= gsl_rng_alloc(gsl_rng_taus);
    const gsl_rng_type * T;
    gsl_rng *rg1;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    rg1 = gsl_rng_alloc (T);

    // gsl_rng * rg1;
    int b[r];
    int *j_sample;
    unsigned int coin;
    for (int i=0; i<r; i++)
    {
      b[i] = i;
    }
    gsl_vector *vtmp = gsl_vector_alloc(n);
    double R_j_norm,prob;
    double *Rw_dot = (double*)malloc(1*sizeof(double));
    // printf("sample from X inside");
    while (keep_going)
    {

        counter += 1;
        // # sample row index uniformly at random
        // i_sample = np.random.choice(r)
        gsl_ran_choose(rg,i_sample,1,b,r,sizeof(int));
        //(rg, a, 1, rows, r, sizeof(int));

        // # sample column index from length-square distribution of corresponding row
        gsl_matrix_get_row(vtmp, LS_prob_columns_R , i_sample[0]);
        j_sample = vose(vtmp,n,1);
        printf("isample %d\n",i_sample[0] );
        // # column j_sample of matrix R

        gsl_vector *R_j = gsl_vector_alloc(r);
        gsl_vector_set_zero(R_j);

        // # compute entries of R_j
        // R_j = (A_Frobenius/sqrt(r)) * R_j
        // gsl_vector_scale(R_j,(A_Frobenius/sqrt(r)));
        
        for (int s=0; s<r; s++)
        {
          gsl_vector_set(R_j,s,gsl_matrix_get(A,rows[s],j_sample[0])/sqrt(row_norms[rows[s]]));
        }
        gsl_vector_scale(R_j,A_Frobenius/sqrt(r));
        // for(int i=0;i<r;i++)
        // {
        //   printf("R_j%f\n", gsl_vector_get(R_j,i));
        // }
        break;
        // # norm of column vector R_j 
        R_j_norm = gsl_blas_dnrm2(R_j);

        // # inner product of R_j and w
        // Rw_dot = np.dot(R_j, w_vector)
        gsl_blas_ddot(R_j,w_vector,Rw_dot);
        // # probability to select j_sample as output
        printf("\n\n\n rw_dot,wnorm, rjnorm------ %f \t %f \t %f\n",*Rw_dot,w_norm,R_j_norm);
        prob = pow(*Rw_dot /( w_norm * R_j_norm),2);
        printf("prob \t %f\n",prob);

        // # determine if we output j_sample given above probability

        coin = gsl_ran_binomial(rg1,prob,1);

        if (coin == 1)
        {
              
              out_j = *j_sample;
            
            // # if we get heads from coin, then stop while loop
              keep_going = 0;
        }
        printf("out_j counter\t %d %d \t%d ", out_j,counter,coin);
    }
    gsl_rng_free (rg1);
    out_j_conter[0] = (out_j);
    out_j_conter[1] = counter;
    // printf("out_j counter\t %d  \t%d ", out_j_conter[0],out_j_conter[1]);
    return out_j_conter;

}
                  

void sample_me_lsyst(const gsl_matrix *A, const gsl_vector *b, int m, int n, int samples, int rank, int r, gsl_matrix *w, int *rows, gsl_vector *sigma, double row_norms[], gsl_vector *LS_prob_rows, const gsl_matrix *LS_prob_columns, double A_Frobenius, double *lambdas)
{
  int reps = 10;
  int *sample_j;
  int med;
  if (reps % 2 == 0)
  {
    med = 0;
  }
  else
  {
    med = 1;
  }
  gsl_matrix *matrix_elements = gsl_matrix_alloc(reps, rank);

  double vector_sum=0;

  for (int i = 0; i < reps; i++)
  {
    for (int l = 0; l < rank; l++)
    {
      gsl_vector *X = gsl_vector_alloc(samples);

      gsl_vector_set_zero(X);

        for(int k=0; k<samples;k++) 
        {
        int *sample_i;
        gsl_vector *vtmp = gsl_vector_alloc(n);

        sample_i = vose(LS_prob_rows, m, 1);

        gsl_matrix_get_row(vtmp, LS_prob_columns, sample_i[0]);

        sample_j = vose(vtmp, n, 1);

// #calculates v_j

        double v_j = 0;

        for (int s = 0; s < r; s++)
        {
              v_j+= gsl_matrix_get(A, rows[s],sample_j[0]) * gsl_matrix_get(w,s,l) / (sqrt(row_norms[rows[s]]));
              // printf("tmp  %f\n" ,v_j);
        }

        v_j = v_j * A_Frobenius / (sqrt(r) * gsl_vector_get(sigma, l));

        // X[k] = ((A_Frobenius * *2 * b[sample_i]) / (A[sample_i, sample_j])) *v_j

        double tmp = pow(A_Frobenius, 2) * gsl_vector_get(b, *sample_i) / gsl_matrix_get(A,*sample_i, *sample_j) * v_j;
        
        gsl_vector_set(X, k, tmp);
        }
        for (int i=0; i<samples; i++)
        {
          
         vector_sum += gsl_vector_get(X,i);
        //  printf("vector sum %f \n", vector_sum);
        
        }

        gsl_matrix_set(matrix_elements,i,l,(vector_sum) / samples);
    }
  }

// #take median of all repeated estimates
    // double tmp1;
    gsl_vector *vtmp1 = gsl_vector_alloc(reps);

    for (int l = 0; l < rank; l++)
    {
      // gsl_stats_median(lambdas, 1, rank);
      gsl_matrix_get_col(vtmp1, matrix_elements, l);
     
      gsl_sort_vector(vtmp1);
      
      if (med == 0)
        lambdas[l] = gsl_vector_get(vtmp1, reps / 2);
      else
        lambdas[l] = (gsl_vector_get(vtmp1, (reps + 1 )/ 2) + gsl_vector_get(vtmp1, (reps - 1) / 2  )) / 2;
      // printf("lambads %f \n",gsl_vector_get(vtmp1,l));
    }
    // printf("sample_me_lsyst");
    
}

void uvl_vector(int m, int n, int l, const gsl_matrix *A, int r, gsl_matrix *w, int *rows, gsl_vector *sigma, double row_norms[],double A_Frobenius, gsl_matrix *ul_approx, gsl_matrix *vl_approx)
{

// #building approximated v ^ l vector
  gsl_vector *tmp_vector = gsl_vector_alloc(n);
  double tmp_val;
  double factor;

  factor = A_Frobenius / ( sqrt(r) * gsl_vector_get(sigma,l) );
  // printf("%f",factor);
    // for s in range(r):
    //     v_approx[:] += ( A[rows[s], :] / np.sqrt(row_norms[rows[s]]) ) * w[s, l]
    gsl_vector *u_approx = gsl_vector_alloc(m);
    gsl_vector *v_approx = gsl_vector_alloc(n);

    for (int s=0; s<r; s++)
    {

      gsl_matrix_get_row(tmp_vector, A, rows[s]);
      tmp_val = sqrt(row_norms[rows[s]]) * gsl_matrix_get(w,s,l);
      gsl_vector_scale(tmp_vector, tmp_val);
      gsl_vector_add(v_approx, tmp_vector);
    }
    gsl_vector_scale(v_approx, factor);
    gsl_matrix_set_col(vl_approx,l,v_approx);
    gsl_blas_dgemv(CblasNoTrans, 1, A, v_approx, 0.0, u_approx);

    gsl_vector_scale(u_approx, 1 / gsl_vector_get(sigma, l));
    gsl_matrix_set_col(ul_approx,l,u_approx);
    // printf("Uvl function");
}

void sample_C(const gsl_matrix *A, int m, int n, int r, int c, int *rows, double row_norms[], gsl_vector *LS_prob_rows,const gsl_matrix *LS_prob_columns, double A_Frobenius, gsl_matrix *C, gsl_matrix *vh, gsl_vector *sigma,gsl_matrix *LS_prob_columns_R)
{
  // int  *rows = malloc(sizeof(int)*r);
  int a[1];
  int *columns_tmp = malloc(sizeof(int) * 1);
  int columns[c];
  clock_t start, end;
  double rt_sampling_C, rt_building_C, rt_svd_C;
  double tmp;
  gsl_vector *v = gsl_vector_alloc(n);
  gsl_vector *R_row = gsl_vector_alloc(n);
  gsl_matrix *R_C = gsl_matrix_alloc(r, c);
  // gsl_matrix *C = gsl_matrix_alloc(r,c);
  gsl_vector *vtmp = gsl_vector_alloc(c);
  gsl_vector *column_norms = gsl_vector_alloc(c);
  gsl_vector *vtmp2 = gsl_vector_alloc(r);
  double R_row_norm;
  // gsl_matrix *LS_prob_columns_R = gsl_matrix_alloc(r, n);
  // gsl_matrix *Vh = gsl_matrix_alloc (c, c);
  gsl_rng *rg = gsl_rng_alloc(gsl_rng_taus);
  // gsl_vector *S = gsl_vector_alloc(c);
  gsl_vector *work = gsl_vector_alloc(c);
  start = clock();
  rows = vose(LS_prob_rows, m, r);
  
  for (int i = 0; i < c; i++)
  {

    gsl_ran_sample(rg, a, 1, rows, r, sizeof(int));

    gsl_matrix_get_row(v, LS_prob_columns, a[0]);

    columns_tmp = vose(v, n, 1);

    columns[i] = columns_tmp[0];
    // printf("%d\n",columns[i]);
  }

  end = clock();

  rt_sampling_C = ((double)(end - start)) / CLOCKS_PER_SEC;

  printf("time elapsed rt_sampling_C %f", rt_sampling_C);

  start = clock();

  for (int s = 0; s < r; s++)
  {

    tmp = A_Frobenius / (sqrt(r) * sqrt(row_norms[rows[s]]));
    gsl_matrix_get_row(R_row, A, rows[s]);
    
    

    gsl_vector_scale(R_row, tmp);
    
    R_row_norm = fabs(pow(gsl_blas_dnrm2(R_row), 2));

    // printf("R_row_norm \t %f\n",R_row_norm);    
    gsl_vector_mul(R_row, R_row);

    gsl_vector_scale(R_row, 1 / R_row_norm);

    gsl_matrix_set_row(LS_prob_columns_R, s, R_row);
  }
  

  for (int s = 0; s < r; s++)
  {
    for (int t = 0; t < c; t++)
    {
      gsl_matrix_set(R_C, s, t, gsl_matrix_get(A, rows[s], columns[t]));
    }
    tmp = A_Frobenius / (sqrt(r) * sqrt(row_norms[rows[s]]));

    gsl_matrix_get_row(vtmp, R_C, s);

    gsl_vector_scale(vtmp, tmp);

    gsl_matrix_set_row(R_C, s, vtmp);

  }
  
  gsl_vector_set_zero(column_norms);
  
  for (int t = 0; t < c; t++)
  {
    for (int s = 0; s < r; s++)
    {

      gsl_vector_set(column_norms, t, gsl_vector_get(column_norms, t) + pow(gsl_matrix_get(R_C, s, t), 2));
    }
  }
   
 
  for (int t = 0; t < c; t++)
  {
    // gsl_matrix_get_col(vtmp2, C, t);

    gsl_matrix_get_col(vtmp2, R_C, t);
    
  
    tmp = (A_Frobenius) / sqrt(gsl_vector_get(column_norms, t)) / sqrt(c);

    gsl_vector_scale(vtmp2, tmp);

    gsl_matrix_set_col(C, t, vtmp2);
  }

//  for (int i=0; i<r; i++)
//   {
//     for (int j=0; j<c; j++)
//     {
//       printf("%f\n",gsl_matrix_get(C,i,j));
//     }
//   }
  end = clock();

  rt_building_C = ((double)(end - start)) / CLOCKS_PER_SEC;

  printf("rt_building_C %f", rt_building_C);

  start = clock();

  gsl_linalg_SV_decomp(C, vh, sigma, work);

  end = clock();

  rt_svd_C = ((double)(end - start)) / CLOCKS_PER_SEC;

  printf("rt_svd_C %f", rt_svd_C);
  // for(int i =0;i<r;i++)
  // {
  //   for(int j=0;j<c;j++)
  //   {
  //     printf("%f\t,",gsl_matrix_get(C,i,j));
  //   }
  //   printf("\n");
  // }
  // for(int j=0;j<c;j++)
  //   {
  //     printf("%f\n,",gsl_vector_get(S,j));
  //   }
  // sample_me_lsyst(A, b, m, n, samples, rank, r, w, rows, sigma, row_norms, LS_prob_rows, LS_prob_columns, A_Frobenius);
}
void linear_eqs(const gsl_matrix *A, const gsl_vector *b, int m, int n, int r, int c, int rank, int Nsamples, int NcompX, double row_norms[], gsl_vector *LS_prob_rows, const gsl_matrix *LS_prob_columns, double A_Frobenius)
{
  // m_rows, n_cols = np.shape(A)
  // int t; 
  
  
  int *rows = malloc(sizeof(int) * r);
  gsl_vector *sigma = gsl_vector_alloc(c);
  gsl_matrix *vh = gsl_matrix_alloc(c, c);
  gsl_matrix *w = gsl_matrix_alloc(r, c);
  // gsl_vector *u_approx = gsl_vector_alloc(m);
  // gsl_vector *v_approx = gsl_vector_alloc(n);
  gsl_matrix *LS_prob_columns_R = gsl_matrix_alloc(r, n);

  double *lambdas = (double *)calloc(rank, sizeof(double));
  clock_t start, end;
  double  rt_sample_me, rt_sampling_sol;
  int *sample_return_from_x = (int*)malloc(sizeof(int) * 2) ;
  gsl_vector *w_vector = gsl_vector_alloc(r);
  gsl_vector_set_zero(w_vector);
  
  
  gsl_vector *vtmp2 = gsl_vector_alloc(r);
   

  // # 1 - Generating LS probability distributions to sample from matrix A

  // start = clock();
  // LS = ls_probs(m, n, A);
  // end = clock(); 

  // rt_ls_prob = ((double)(end - start)) / CLOCKS_PER_SEC;

  // # 2 - Building matrix C by sampling "r" rows and "c" columns from matrix A and computing SVD of matrix C
  printf("Linear EQ");
  sample_C(A, m, n, r, c, rows,row_norms, LS_prob_rows, LS_prob_columns, A_Frobenius, w, vh, sigma,LS_prob_columns_R);
  // exit(0);
  
  gsl_matrix *ul_approx = gsl_matrix_alloc(m, rank);
  gsl_matrix *vl_approx = gsl_matrix_alloc(n, rank);
  gsl_matrix_set_zero(ul_approx);
  gsl_matrix_set_zero(vl_approx);
  // for l in range(rank):
  //     ul_approx[:, l], vl_approx[:, l] =
  for (int l = 0; l < rank; l++)
  {
    uvl_vector(m,n,l, A, r, w, rows, sigma, row_norms, A_Frobenius,ul_approx,vl_approx);
  }
  // for (int i=0; i<n; i++)
  // {
  //   for (int j=0; j<rank; j++)
  //   {
  //     printf("%f ,",gsl_matrix_get(vl_approx,i,j));
  //   }
  //   printf("\n");
  // }

  // printf("\n AFTER UVL ---------=-=-===-=-=== %f\n",gsl_matrix_get(LS_prob_columns_R,10,10));
  // start = clock();
  sample_me_lsyst(A, b, m, n, Nsamples, rank, r, w, rows, sigma, row_norms, LS_prob_rows, LS_prob_columns, A_Frobenius, lambdas);
  // exit(0);
  // end = clock();
  // rt_sample_me = ((double)(end - start)) / CLOCKS_PER_SEC;
  // printf("\n\nrt_sample_me vinooth-------------\n\n %f",rt_sample_me);

    // start = clock();
    // for(t=0;t<NcompX; t++)
    // {
    //   printf("samplex %d",NcompX);
    // }

    
    
    for (int k = 0; k < rank; k++)
    {
      gsl_matrix_get_col(vtmp2, w, k);
      // printf("lambdas %f\n",lambdas[k]);
      gsl_vector_scale(vtmp2, pow( lambdas[k] / gsl_vector_get(sigma,k),3));

      gsl_vector_add(w_vector, vtmp2);
      
    }
    
    // for(int i=0;i<r;i++ )
    // {
    //   printf(" w_vector %f\n",gsl_vector_get(w_vector,i));
    // }
    double w_norm = gsl_blas_dnrm2(w_vector);

    // sampled_comp = np.zeros(NcompX, dtype = np.uint
    printf("wnorm \n\n\n\n\t%f" ,w_norm);
    gsl_vector *sampled_comp = gsl_vector_alloc(NcompX);
    gsl_vector_set_zero(sampled_comp);

    

    gsl_vector *no_of_rejected_samples = gsl_vector_alloc(NcompX);
    gsl_vector_set_zero(no_of_rejected_samples);
  //  printf("\n after saample beffore X \n");
    for (int t=0; t<NcompX; t++)
    {
      // printf("vinooth before samplex -- %d",NcompX);
      // break;
      // sample_return_from_x =
      sample_return_from_x = sample_from_x(A, r, c, n, rows, row_norms, LS_prob_columns_R, A_Frobenius, w_vector, w_norm);
      gsl_vector_set(sampled_comp,t,sample_return_from_x[0]);
      gsl_vector_set(no_of_rejected_samples,t,sample_return_from_x[1]); 
    }
    // printf("\n after saample from X \n");
    // end = clock();
   
    // for(int i=0;i<NcompX;i++ )
    // {
    //   // printf(" sampled_comp %f\n",gsl_vector_get(sampled_comp,i));
    //   printf(" no_of_rej %f\n",gsl_vector_get(no_of_rejected_samples,i));
    // }
    // rt_sampling_sol = ((double)(end - start)) / CLOCKS_PER_SEC;
    // printf("rt_svd_C %f", rt_sampling_sol);
    // for t in range(NcompX):
    //     x_tilde[t] = approx_solution(A, rank, r, w, svd_C[1], svd_C[2],
    //                                  LS[0], LS[3], lambdas, sampled_comp[t])
    

    gsl_vector *x_tilde = gsl_vector_alloc(NcompX);

    
    for (int t=0; t<NcompX; t++)
    {

      gsl_vector_set(x_tilde,t, approx_solution(A, rank, r, w, rows, sigma, row_norms, A_Frobenius, lambdas, gsl_vector_get(sampled_comp,t)));

    }


}

int main()
{
  int m = 500, n = 250; //610, 9724
  int r = 200, c = 200;
  int Nsamples = 50, NcompX = 50;
  int rank = 3;
  char *record;

  // char *recordb;
  double A_Frobenius;
  double row_norms[m];
  double LS_prob_row_tmp;
  double LS_prob_column_tmp;
  char *line = NULL;
  char *lineb = NULL;
  size_t len = 0;
  size_t lenb = 0;
  size_t read;
  size_t readb;
  int i = 0, j = 0;
  gsl_matrix *a = gsl_matrix_alloc(m, n);
  gsl_matrix *LS_prob_columns = gsl_matrix_alloc(m, n);
  gsl_vector *v = gsl_vector_alloc(n);
  // gsl_vector *v = gsl_vector_alloc(b);
  gsl_vector *LS_prob_rows = gsl_vector_alloc(m);
  FILE *fstream = fopen("Data/A.csv", "r");
  gsl_vector *b = gsl_vector_alloc(m);
  FILE *fstreamb = fopen("Data/B.csv", "r");
  srand(time(NULL));
  if (fstream == NULL)
  {
    printf("\n file-1 opening failed  ");
    return -1;
  }
  if (fstreamb == NULL)
  {
    printf("\n file-2 opening failed  ");
    return -1;
  }
  while (((read = getline(&line, &len, fstream)) != -1))
  {
    //
    j = 0;
    record = strtok(line, ",");

    while (record != NULL)
    {
      gsl_matrix_set(a, i, j, atof(record));
      ++j;
      record = strtok(NULL, ",");
    }
    ++i;
  }
  int ib = 0;

  while (((readb = getline(&lineb, &lenb, fstreamb)) != -1))
  {
    gsl_vector_set(b,ib, atof(lineb));
    ib++;
  }

  A_Frobenius = frobeniusNorm(a, m, n);
  printf("frobenius norm %f", A_Frobenius);
  for (i = 0; i < m; i++)
  {
    gsl_matrix_get_row(v, a, i);
    row_norms[i] = pow(gsl_blas_dnrm2(v), 2);
  }
  for (i = 0; i < m; i++)
  {

    LS_prob_row_tmp = row_norms[i] / pow(A_Frobenius, 2);
    gsl_vector_set(LS_prob_rows, i, LS_prob_row_tmp);
  }
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      LS_prob_column_tmp = (pow(gsl_matrix_get(a, i, j), 2) / row_norms[i]);
      gsl_matrix_set(LS_prob_columns, i, j, LS_prob_column_tmp);
    }
  }
  // sample_C(a,m,n,r,c,row_norms,LS_prob_rows,LS_prob_columns,A_Frobenius);
  linear_eqs(a, b,m,n, r, c, rank, Nsamples, NcompX, row_norms, LS_prob_rows, LS_prob_columns, A_Frobenius);
  gsl_matrix_free(a);
  return 0;
}
