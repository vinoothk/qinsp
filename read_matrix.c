#include <stdio.h>          // Needed for printf() and feof()
#include <math.h>           // Needed for pow().   
#include <stdlib.h>         // Needed for exit() and atof()
#include <string.h>         // Needed for strcmp()
#include <gsl/gsl_matrix.h> // Needed for matrix manipulations

/*------------ Defines -----------------------------------------------------------------*/
#define MAX_SIZE  10000000   // Maximum size of the time series array
#define NUM_LAG   1000       // Number of lags to calculate for


/*
 *----------------------------------------------------------------------------------------
 * Main program
 *----------------------------------------------------------------------------------------
 */
int main(int argc, char* argv[]) {

//--------Initialization--------------------------------------------------------------
    double ac_value;      // computed autocorrelation value
    int i,j;              // Loop counter
    long int N;
    double mean, variance;
    gsl_matrix * chains;

    char filename[100];
    FILE* in_file;      // input file
    FILE* out_file;     // output file

    int no_params;      // number of parameters to calculate autocorrelation for
    int first_column;   // Which column first corresponds to a chain
    int ch;             // to determine number of samples in file

    printf("-------------------------------------------------------------------------\n");
    //--------Check that there are the correct number of arguments passed-----------------
    if(argc != 4) { 
        printf("usage: ./auto_corr chainfile no_params first_column \n");
        exit(1); // 0 means success typically, non-zero indicates an error
}

    //--------Extract arguments-----------------------------------------------------------
    sprintf(filename,"%s",argv[1]); // convert input file to string
    in_file = fopen(filename,"rb");  // open input file for reading

    no_params = atoi(argv[2]);      
    first_column = atoi(argv[3]);

    //--------What is the number of samples in chain file?--------------------------------
    N = 0; // Initialize count
    while(!feof(in_file)) {
        ch = fgetc(in_file);
        if(ch == '\n'){
            N++;
        }
    }
    printf("Number of samples: %li\n", N); // print number of samples
    if (N > MAX_SIZE) { // throw error if there are too many samples
        printf("ERROR - Too many samples! MAX_SIZE = %i", MAX_SIZE);
        exit(2);
    }

    //--------Generate a gsl matrix from the chains---------------------------------------
    printf("%i\n", no_params);
    chains = gsl_matrix_alloc(N, no_params); // allocate memory for gsl_matrix(rows, cols)
    // print the matrix (for testing)
    printf("Chain matrix \n");
    for (int m=0;m<N;m++) { //rows
        for (int n=0; n<no_params;n++) { // columns
            printf("%f  ",gsl_matrix_get(chains,m,n));
        }
        printf("\n");
    } 
    // gsl_matrix_fprintf(stdout,chains,"%f"); // easy way to print, no formatting though
    gsl_matrix_fscanf(in_file, chains);        // read in chains to the gsl_matrix
    fclose(in_file);