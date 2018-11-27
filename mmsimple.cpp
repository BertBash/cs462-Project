#include <iostream>
#include <cmath>
#include <mpi.h>
#include <ctime>
#include <cstdllib>

using namespace std;

int main(){
  int size, rank;
  double time;
  
  MPI_Init(NULL,NULL);
  
  time = MPI_Wtime();
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int row_color = rank / sqrt(size);
  int col_color = rank % sqrt(size);
  
  int num_row, num_col, i, j;
  num_col = num_row = (128 / sqrt(size));
  srand(time(NULL)+rank);
  
  MPI_Comm row_comm;
  MPI_Comm col_comm;
  
  float **A; // Entire row beginning to end
  float **B; // Entire column beginning to end
  float **C; // Subsection result
  
  MPI_Comm_split(MPI_COMM_WORLD, row_comm, rank, &row_comm);
  MPI_Comm_split(MPI_COMM_WORLD, col_comm, rank, &col_comm);
  
  
  //Generate Matricies
  A = (float**) malloc(sizeof(float*) * num_row);
  for(i = 0; i < num_row; i++) {
    A[i] = (float*) malloc(sizeof(float) * 128);
    for(j = 0; j < 128; j++)
      A[i][j] = (float) (rand/ (float) RAND_MAX);
  }
  B = (float**) malloc(sizeof(float*) * 128);
  for(i = 0; i < num_row; i++) {
    B[i] = (float*) malloc(sizeof(float) * num_col);
    for(j = 0; j < num_col; j++)
      B[i][j] = (float) (rand/ (float) RAND_MAX);
  }
  C = (float**) malloc(sizeof(float*) * num_row);
  for(i = 0; i < num_row; i++) {
    C[i] = (float*) malloc(sizeof(float) * num_col);
  }

  if(rank == 0)
    cout << "Matricies Set" << endl;
  
  //Communicate
  
  //Crunch Numbers  
  
  
  if(rank == 0)
    cout << "Calculated" << endl;

  //Check
    if(rank == 0)
      cout << "Checked" << endl;
    
  //Finish
  if(rank == 0)
    cout << "Exiting" << endl;
    
  MPI_Finalize();
  
  return 0;
}
