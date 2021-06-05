#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
int main()
{
//    char buffer[4096*] ;
   char *record;
// *line;
     char * line = NULL;
    size_t len = 0;
    ssize_t read;
   int i=0,j=0;
//    int mat[100][100];
    gsl_matrix *a = gsl_matrix_alloc (500, 250);
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
//    for (i = 0; i < 500; i++)
//     for (j = 0; j < 250; j++)
//       {
//         // double mij = gsl_matrix_get (m, i, j);
//         double aij = gsl_matrix_get (a, i, j);
//         printf("%f \n",aij);
//         // return;s
//         // if (mij != aij) k++;
//       }
   return 0 ;
 }