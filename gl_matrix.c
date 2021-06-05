#include <stdio.h>
#include <gsl/gsl_matrix.h>

int
main (void)
{
  int i, j, k = 0;
  // FILE *f;
//   gsl_matrix * m = gsl_matrix_alloc (100, 100);
  gsl_matrix *a = gsl_matrix_alloc (500, 250);
  // gsl_matrix *a = gsl_matrix_alloc (500, 250);

//   for (i = 0; i < 100; i++)
//     for (j = 0; j < 100; j++)
//       gsl_matrix_set (m, i, j, 0.23 + i + j);

//   {
//      FILE * f = fopen ("test.dat", "wb");
//      gsl_matrix_fwrite (f, m);
//      fclose (f);
//   }

  // {
  //    FILE * f = fopen ("A.csv", "rb");
  //    gsl_matrix_fread (f, a);
  //    fclose (f);
  // }

  // {
  //    FILE * f = fopen ("B.dat", "wb");
  //    gsl_matrix_fwrite (f, a);
  //    fclose (f);
  // }
  // f=fopen("A.dat","rb");
  // gsl_matrix_fscanf(f,a);
  // fclose(f);

  // /* print matrix the easy way */
  // printf("Matrix mm\n");
  // gsl_matrix_fprintf(stdout,a,"%f");

  for (i = 0; i < 500; i++)
    for (j = 0; j < 250; j++)
      {
        // double mij = gsl_matrix_get (m, i, j);
        double aij = gsl_matrix_get (a, i, j);
        printf("%g \n",aij);
        // return;s
        // if (mij != aij) k++;
      }
  
 
//   gsl_matrix_fprintf(stdout,a , "%g");
// //   gsl_matrix_free (m);
  // gsl_matrix_free (a);

//   printf ("differences = %d (should be zero)\n", k);
//   return (k > 0);
} 