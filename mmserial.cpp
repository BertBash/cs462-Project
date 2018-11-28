#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(){
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
  
  //Matrix multiply
  for(int i = 0; i < 128; i++){
    for(int j = 0; j < 128; j++){
      for(int k = 0; k < 128; k++){
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }

  ofstream output;

  output.open("serial.output");
  // Print the output to a file for testing
  for(int i = 0; i < 128; i++) {
  	  for(int j = 0; j < 128; j++) {
  	  	   output << C[i][j];
  	  }
  	  output << endl;
  }
  
  return 0;
}
