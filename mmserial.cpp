#include <iostream>
#include <fstream>
#include <cstdlib>
#include <mpi.h>

using namespace std;

int main(){
	MPI_Init(NULL, NULL);
	double start, end;
  float A[128][128];
  float B[128][128];
  float C[128][128];
  srand(-1337); 
  //Popuate the matrixes with random numbers between 0 and 1.
  for(int i = 0; i < 128; i++){
    for(int j = 0; j <128; j++){
      A[i][j] = (float) (rand() / ((float) RAND_MAX / 2)) - 1;
      B[i][j] = (float) (rand() / ((float) RAND_MAX / 2)) - 1;
      C[i][j] = 0;
    }
  }
  start = MPI_Wtime();
  //Matrix multiply
  for(int i = 0; i < 128; i++){
    for(int j = 0; j < 128; j++){
      for(int k = 0; k < 128; k++){
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
	end = MPI_Wtime();
	cout << "Elapsed time: " << end - start << endl;
  ofstream output;

  output.open("serial.output");
  // Print the output to a file for testing
  for(int i = 0; i < 128; i++) {
  	  for(int j = 0; j < 128; j++) {
  	  	   output << C[i][j] << " ";
  	  }
  	  output << endl;
  }
  MPI_Finalize();
  return 0;
}
