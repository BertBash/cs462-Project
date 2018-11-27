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
  double start_time, end_time;
  
  MPI_Init(NULL,NULL);
  
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int row_color = rank / (int)sqrt((double)size);
  int col_color = rank % (int)sqrt((double)size);
  
  int num_row, num_col, i, j, k, n, m, num_elements;
  num_col = num_row = (128 / sqrt(size));
  srand(time(NULL)+rank);
  
  MPI_Comm row_comm;
  MPI_Comm col_comm;
  
  float *A_start, *B_start; // Initialize data
  float *A_chunk; // Entire row chunk beginning to end
  float *B_chunk; // Entire column chunk beginning to end
  float *C_chunk; // Subsection result
  float *C; // Final matrix
  
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
  if(rank==0) C = (float*) malloc(sizeof(float) * 128 * 128);
 
  for(i = 0; i < num_row * num_col; i++) {
    A_start[i] = (float) (rand()/ (float) RAND_MAX);
    B_start[i] = (float) (rand()/ (float) RAND_MAX);
  }
  
  // Get all processors sync'd
	MPI_Barrier(MPI_COMM_WORLD);

	// Grab timestamp and continue with communication and number crunching.
	start_time = MPI_Wtime();
//  if(rank == 0) cout << "Start time: " << start_time << endl;

//  if(rank == 0)
//    cout << "Matricies Set" << endl;
  
  //Communicate
 
  // Each process is sending their A and B (horizontal for A, vertical for B
  MPI_Allgather(A_start, num_row * num_col, MPI_FLOAT, A_chunk, num_row * num_col, MPI_FLOAT, row_comm);
  MPI_Allgather(B_start, num_row * num_col, MPI_FLOAT, B_chunk, num_col * num_row, MPI_FLOAT, col_comm);

  //Crunch Numbers  
  // size * (k / length) + row * length + k % length
  num_elements = num_row * num_col;
	for(i = num_row * row_rank, n = 0; i < num_row + num_row * row_rank; i++, n++) {
		for(j = num_col * col_rank, m = 0; j < num_col + num_col * col_rank; j++, m++) {
			C_chunk[m+num_col*n] = 0;
			for(k = 0; k < 128; k++) {
				C_chunk[m+num_col*n] += A_chunk[num_elements * (k / num_row) + i * num_row + k % num_row] 
											 * B_chunk[num_elements * (k / num_col) + j * num_col + k % num_col];
			}
		}
	}
  
//  if(rank == 0)
 //   cout << "Calculated" << endl;
//	if(rank == 0) cout << "End time: " << (end_time = MPI_Wtime()) << endl;
	MPI_Gather(C_chunk, num_elements, MPI_FLOAT, C, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);
 	end_time = MPI_Wtime();
	if(rank == 0) cout << "Elapsed time: " << end_time - start_time << endl;
  //Check
//    if(rank == 0)
//      cout << "Checked" << endl;
    
  //Finish
//  if(rank == 0)
//    cout << "Exiting" << endl;
    
  MPI_Finalize();
  
  return 0;
}
