#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <math.h>

double frobeniusNorm(const gsl_matrix *matrix, int size1, int size2)
{
    double result = 0.0;
    for(int i = 0; i < size1; ++i)
    {
        for(int j = 0; j < size2; ++j)
        {
            double value = (gsl_matrix_get(matrix,i,j));
            printf("fffff norm %f",gsl_matrix_get(matrix,i,j));
            result += value * value;
            return;
        }
    }
    // gsl_matrix_fprintf(stdout,matrix, "%g");
    printf("fffff norm %f",result);
    return sqrt(result);
}
int
main (void)
{
  double a_data[] = { 0.18, 0.60, 0.57, 0.96,
                      0.41, 0.24, 0.99, 0.58,
                      0.14, 0.30, 0.97, 0.66,
                      0.51, 0.13, 0.19, 0.85 };

  gsl_matrix *a = gsl_matrix_alloc (500, 250);                   

//   const double b_data[] = { 0.18, 0.60, 0.57, 0.96 };
  double x;
  // int i;
  double y;
  int i,j;
//   gsl_matrix_view m
//     = gsl_matrix_view_array (a_data, 4, 4);

//   gsl_vector b =  { 1.0, 2.0, 3.0, 4.0 };
  // gsl_vector * v = gsl_vector_alloc (4);
 
  {
     FILE * f = fopen ("B.dat", "rb");
     gsl_matrix_fread (f, a);
     fclose (f);
  }

  for (i = 0; i < 500; i++)
    for (j = 0; j < 250; j++)
      {
        // double mij = gsl_matrix_get (m, i, j);
        double aij = gsl_matrix_get (a, i, j);
        // if (mij != aij) k++;
      }
  // gsl_vector_fscanf(stdin,v);
       
      //  for (i = 0; i < 4000000; i++)
      //    {
      //      gsl_vector_set (v, i, i+1);
      //    }
       
    //    for (i = 0; i < 4000000; i++)
    //      {
    //        printf ("v_%d = %g\n", i, gsl_vector_get (v, i));
    //      }
     
    
    // = gsl_vector_view_array (b_data, 4);

//   gsl_vector *x = gsl_vector_alloc (4);

//   int s;

//   gsl_permutation * p = gsl_permutation_alloc (4);

//   gsl_linalg_LU_decomp (&m.matrix, p, &s);

//   gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

//   printf ("x = \n");
//   gsl_vector_fprintf (stdout, x, "%g");

//   gsl_permutation_free (p);
//   gsl_vector_free (x);

    // x = gsl_blas_dnrm2(v);
    // printf("euclidian norm %f",pow(x,2));
    // gsl_matrix_fprintf(stdout,a , "%g");
    y = frobeniusNorm(a,500,250);
    printf("frobenius norm %f",y);

  return 0;
}
