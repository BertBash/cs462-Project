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
  
  float *A_start, *B_start; // Initialize data
  float *A_chunk; // Entire row chunk beginning to end
  float *B_chunk; // Entire column chunk beginning to end
  float *C_chunk; // Subsection result
  
  MPI_Comm_split(MPI_COMM_WORLD, row_comm, rank, &row_comm);
  MPI_Comm_split(MPI_COMM_WORLD, col_comm, rank, &col_comm);
  
  
  //Generate Matricies
  A_chunk = (float*) malloc(sizeof(float) * num_row * 128);
  B_chunk = (float*) malloc(sizeof(float) * num_col * 128);
  C_chunk = (float*) malloc(sizeof(float) * num_row * num_col);

  A_start = (float*) malloc(sizeof(float) * num_row * num_col);
  B_start = (float*) malloc(sizeof(float) * num_row * num_col);
 
  for(i = 0; i < num_row * num_col; i++) {
    A_start[i] = (float) (rand/ (float) RAND_MAX);
    B_start[i] = (float) (rand/ (float) RAND_MAX);
  }
  


  if(rank == 0)
    cout << "Matricies Set" << endl;
  
  //Communicate
 
  MPI_Alltoall(A_start, num_row * num_col, MPI_FLOAT, A_chunk, num_row * 128, MPI_FLOAT, row_comm);
  MPI_Alltoall(B_start, num_row * num_col, MPI_FLOAT, B_chunk, num_col * 128, MPI_FLOAT, col_comm);

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
