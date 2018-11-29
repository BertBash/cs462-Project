#include <iostream>
#include <cmath>
#include <mpi.h>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <fstream>

using namespace std;

int main(){
  
  MPI_Init(NULL,NULL);

  int size, rank;              // Number of processors, Rank of current processor
  int row_size, row_rank;      // Number of processors in row comm, rank in that comm
  int col_size, col_rank;  	 // Same as above for columns
  float *A, *B, *C; 				 // Matrix A, B, C initialized on Processor 0
  float *A_scatter, *B_scatter, *C_gather;
  double start_time, end_time; // Book keeping time stamps
  float *A_start, *B_start;    // Initialize data
  float *A_chunk;              // Entire row chunk beginning to end
  float *B_chunk;              // Entire column chunk beginning to end
  float *C_chunk;              // Subsection result
  int num_row, num_col;			 // Number of elements on a row and column of a processor's section
  int i, j, k, n, m;				 // Loop vars
  int sqrt_p;									// Number of processors in one dimension
  int num_elements;            // num_row * num_col
  int row_color, col_color;	 // Row and column id
  MPI_Comm row_comm;				 // Communicator for its row
  MPI_Comm col_comm;				 // Communicator for its col
  
  
  // Initialize the matrices
  if(rank == 0){
    A = (float *) malloc(sizeof(float)*128*128);
    B = (float *) malloc(sizeof(float)*128*128);
    C = (float *) malloc(sizeof(float)*128*128);
		C_gather = (float*) malloc(sizeof(float)*128*128);
    srand(-1337); 
    //Popuate the matrixes with random numbers between 0 and 1.
    for(int i = 0; i < 128; i++){
      for(int j = 0; j <128; j++){
        A[i*128+j] = (float) (rand() / ((float) RAND_MAX / 2)) - 1;
        B[i*128+j] = (float) (rand() / ((float) RAND_MAX / 2)) - 1;
        C[i*128+j] = 0;
      }
    }
  }

  // Grab timestamp.
  if(rank == 0) start_time = MPI_Wtime();
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Get the communicator ID that the processor will be grouped by with comms
  row_color = rank / (int)sqrt((double)size);
  col_color = rank % (int)sqrt((double)size);
 
  num_col = num_row = (128 / sqrt(size));
  num_elements = num_row * num_col;
  
  sqrt_p = sqrt((double)size);

  MPI_Comm_split(MPI_COMM_WORLD, row_color, rank, &row_comm);
  MPI_Comm_split(MPI_COMM_WORLD, col_color, rank, &col_comm);
  
  MPI_Comm_rank(row_comm, &row_rank);
  MPI_Comm_rank(col_comm, &col_rank);
  MPI_Comm_size(row_comm, &row_size);
  MPI_Comm_size(col_comm, &col_size);

  // Allocate Matrix sections
  A_chunk = (float*) malloc(sizeof(float) * num_row * 128);
  B_chunk = (float*) malloc(sizeof(float) * num_col * 128);
  C_chunk = (float*) malloc(sizeof(float) * num_elements);

  A_start = (float*) malloc(sizeof(float) * num_row * num_col);
  B_start = (float*) malloc(sizeof(float) * num_row * num_col);

  
  // Fill in A_start, B_start (which is the subsection of the matrix before gathers)
//cerr << "sqrt_p: " << sqrt_p << endl;
  // Fill in A_scatter and B_scatter so that the entire matrix will be sent for each chunk when using MPI_Scatter
  if(rank == 0) {
  	cerr << rank << " is here" << endl;
  	A_scatter = (float*)malloc(sizeof(float) * 128 * 128);
  	B_scatter = (float*)malloc(sizeof(float) * 128 * 128);
		for(i = 0; i < size; i++) {
			for(j = 0; j < num_elements; j++) {
				
//				cerr << "Size: " << size << ", j: " << j << ", ";
//				cerr << "A_scatter[" << j + i * num_elements << "] = A[" << j % num_row + (j / num_row) * 128 + (i % sqrt_p) * num_row + (i / sqrt_p) * num_col * 128 << "]" << endl;
//				cerr << "j: " << j << endl;
//				cerr << "num_row: " << num_row << endl;
//				cerr << "i: " << i << endl;
//				cerr << "sqrt_p: " << sqrt_p << endl;
//				cerr.flush();
				A_scatter[j + i * num_elements] = A[j % num_row + (j / num_row) * 128 + (i % sqrt_p) * num_row + (i / sqrt_p) * num_col * 128];
				B_scatter[j + i * num_elements] = B[j % num_row + (j / num_row) * 128 + (i % sqrt_p) * num_row + (i / sqrt_p) * num_col * 128];
			}
		}
  	cerr << rank << " is here after" << endl;
		/*
    0for(i = 0; i < sqrt_p; i++) {
    	for(j = 0; j < sqrt_p; j++) {
    		if(i == j && i == 0) continue;
				// Using A_start and B_start as tmp holders for what's being sent
				for(n = 0; n < num_row; n++) {
					for(m = 0; m < num_col; m++) {
						A_start[m + n * num_col] = A[128 * (n + i * num_row) + (m + j * num_col)];
						B_start[m + n * num_col] = B[128 * (n + i * num_row) + (m + j * num_col)];
					}
				}
				MPI_Send();
				MPI_Send();
 	   }
 	  }*/
 	  // Fill A_start and B_start for proc 0
 	  //
 	  //
 	  //
  }
  MPI_Barrier(MPI_COMM_WORLD);
	// Scatter out the "Adjusted" matrices
	MPI_Scatter(A_scatter, num_elements, MPI_FLOAT, A_start, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Scatter(B_scatter, num_elements, MPI_FLOAT, B_start, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);

  
  // Communicate
 
  // Each process is sending their A and B (horizontal for A, vertical for B
  MPI_Allgather(A_start, num_elements, MPI_FLOAT, A_chunk, num_elements, MPI_FLOAT, row_comm);
  MPI_Allgather(B_start, num_elements, MPI_FLOAT, B_chunk, num_elements, MPI_FLOAT, col_comm);


  //Crunch Numbers  
  for(i = 0; i < num_row; i++) {
    for(j = 0; j < num_col; j++) {
      C_chunk[j+num_col*i] = 0;
      for(k = 0; k < 128; k++) {
        C_chunk[j+num_col*i] += A_chunk[num_elements * (k / num_row) + i * num_row + k % num_row] 
                              * B_chunk[k * num_col + j];
			}
	 	}
  }
  
  MPI_Gather(C_chunk, num_elements, MPI_FLOAT, C_gather, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);

	// Shift C back to standard matrix format

	if(rank == 0) {
		for(i = 0; i < size; i++) {
			for(j = 0; j < num_elements; j++) {
				C[j % num_row + 128 * (j / num_col + (i / sqrt_p) * num_row)] = C_gather[j + i * num_elements];
			}
		}
	}

  end_time = MPI_Wtime();
  if(rank == 0) cout << "Elapsed time: " << end_time - start_time << endl;
  if(rank == 0) {
		ofstream output;
		output.open("simple.output");
		for(int i = 0; i < 128; i++) {
  	  for(int j = 0; j < 128; j++) {
				output << C[j + i * 128];
			}
			output << endl;
  	}
  }
    
  //Finish
  MPI_Finalize();
  
  return 0;
}
