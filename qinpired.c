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
  char *record;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  int i=0,j=0;
  int x,y;
  // // int i;
double output;
  
  gsl_matrix *a = gsl_matrix_alloc (500, 250);
//    int mat[100][100];
   FILE *fstream = fopen("A.csv","r");
   if(fstream == NULL)
   {
      printf("\n file opening failed ");
      return -1 ;
   }
   while(((read = getline(&line, &len, fstream)) != -1))
    // (line=fgets(buffer,sizeof(buffer),fstream))!=NULL)
   {
    // 
    j=0;
    
     record = strtok(line,",");
    //  printf("line : %s",line) ;
     
     while(record != NULL)
     {
    printf("record : %f",atof(record)) ;
    //here you can put the record into the array as per your requirement.
    //  a[i][j++] = atof(record) ;
    gsl_matrix_set (a, i, j,atof(record));
    //  printf("%d \t  %d\n",i,j);
      ++j;
     record = strtok(NULL,",");
       
     }
     ++i ;   
   }          
 
  // {
  //    FILE * f = fopen ("B.dat", "rb");
  //    gsl_matrix_fread (f, a);
  //    fclose (f);
  // }

  for (x = 0; x < 500; x++)
    for (y = 0; y < 250; y++)
      {
        // double mij = gsl_matrix_get (m, i, j);
        double aij = gsl_matrix_get (a, x, y);
        printf("aij = %f\n",aij);
        // if (mij != aij) k++;
      }
  
      output = frobeniusNorm(a,500,250);
    printf("frobenius norm %f",output);

  return 0;
}
