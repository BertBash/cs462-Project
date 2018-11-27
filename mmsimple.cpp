#include <iostream>
#include <cmath>
#include <mpi.h>
#include <ctime>
#include <cstdlib>
#include <cstdio>

using namespace std;

int main(){
  int size, rank;
  int row_size, row_rank;
  int col_size, col_rank;
  double start_time;
  
  MPI_Init(NULL,NULL);
  
  start_time = MPI_Wtime();
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int row_color = rank / (int)sqrt((double)size);
  int col_color = rank % (int)sqrt((double)size);
  
  int num_row, num_col, i, j;
  num_col = num_row = (128 / sqrt(size));
  srand(time(NULL)+rank);
  
  MPI_Comm row_comm;
  MPI_Comm col_comm;
  
  float *A_start, *B_start; // Initialize data
  float *A_chunk; // Entire row chunk beginning to end
  float *B_chunk; // Entire column chunk beginning to end
  float *C_chunk; // Subsection result
  
  MPI_Comm_split(MPI_COMM_WORLD, row_color, rank, &row_comm);
  MPI_Comm_split(MPI_COMM_WORLD, col_color, rank, &col_comm);
  
  MPI_Comm_rank(row_comm, &row_rank);
  MPI_Comm_rank(col_comm, &col_rank);
  MPI_Comm_size(row_comm, &row_size);
  MPI_Comm_size(col_comm, &col_size);


  
  //Generate Matricies
  A_chunk = (float*) malloc(sizeof(float) * num_row * 128);
  B_chunk = (float*) malloc(sizeof(float) * num_col * 128);
  C_chunk = (float*) malloc(sizeof(float) * num_row * num_col);

  A_start = (float*) malloc(sizeof(float) * num_row * num_col);
  B_start = (float*) malloc(sizeof(float) * num_row * num_col);
 
  for(i = 0; i < num_row * num_col; i++) {
    A_start[i] = (float) (rand()/ (float) RAND_MAX);
    B_start[i] = (float) (rand()/ (float) RAND_MAX);
  }
  


  if(rank == 0)
    cout << "Matricies Set" << endl;
  
  //Communicate
 
  MPI_Alltoall(A_start, num_row * num_col, MPI_FLOAT, A_chunk, num_row * 128, MPI_FLOAT, row_comm);
  MPI_Alltoall(B_start, num_row * num_col, MPI_FLOAT, B_chunk, num_col * 128, MPI_FLOAT, col_comm);

  //Crunch Numbers  
	
	printf("(%d/%d)Processor[%d][%d] is on range %d x %d\n", rank, size, row_color, col_color, row_color * num_row, col_color * num_col);

  
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
