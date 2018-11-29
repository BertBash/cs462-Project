#include <iostream>
#include <cmath>
#include <mpi.h>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>

using namespace std;

int main(){

//  int dim[1], period[1], reorder;
//  int up, down, right, left;
//  MPI_Comm row_cart_comm, col_cart_comm;
  
  MPI_Init(NULL,NULL);

  int print = 1; // Set to 1 to print an output file matrix


  float *A_tmp, *B_tmp;    // Initialize data
  float *A_chunk;              // Entire row chunk beginning to end
  float *B_chunk;              // Entire column chunk beginning to end
  float *C_chunk;              // Subsection result
  float *C_gather; // Final matrix
  float *A, *B, *C;
  float *A_scatter, *B_scatter;
  int size, rank;
  int row_size, row_rank;
  int col_size, col_rank;
  double start_time, end_time;
  int sqrt_p = sqrt(size);
  int row_color = rank / sqrt_p;
  int col_color = rank % sqrt_p;
  int num_row, num_col;
  int i, j, k, n, m;
  int num_elements;
  MPI_Comm row_comm;
  MPI_Comm col_comm;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(rank == 0){
    A = (float*) malloc(sizeof(float)*128*128);
    B = (float*) malloc(sizeof(float)*128*128);
    C = (float*) malloc(sizeof(float)*128*128);
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
  sqrt_p = sqrt(size);


  row_color = rank / sqrt_p;
  col_color = rank % sqrt_p;
  
  num_col = num_row = (128 / sqrt_p);
  num_elements = num_row * num_col;
  
  MPI_Comm_split(MPI_COMM_WORLD, row_color, rank, &row_comm);
  MPI_Comm_split(MPI_COMM_WORLD, col_color, rank, &col_comm);
  

  MPI_Comm_rank(row_comm, &row_rank);
  MPI_Comm_rank(col_comm, &col_rank);
  MPI_Comm_size(row_comm, &row_size);
  MPI_Comm_size(col_comm, &col_size);

  // Allocate Matricies
  A_chunk = (float*) malloc(sizeof(float) * num_row * num_col);
  B_chunk = (float*) malloc(sizeof(float) * num_row * num_col);
  A_tmp = (float*) malloc(sizeof(float) * num_row * num_col);
  B_tmp = (float*) malloc(sizeof(float) * num_row * num_col);
  C_chunk = (float*) malloc(sizeof(float) * num_row * num_col);

 
  for(i = 0; i < num_row * num_col; i++) {
    C_chunk[i] = 0;
  }

	if(rank == 0) {
		A_scatter = (float*)malloc(sizeof(float) * 128 * 128);
    B_scatter = (float*)malloc(sizeof(float) * 128 * 128);
    for(i = 0; i < size; i++) {
      for(j = 0; j < num_elements; j++) {
        A_scatter[j + i * num_elements] =
          A[j % num_row + (j / num_row) * 128 + (i % sqrt_p) * num_row + (i / sqrt_p) * num_col * 128];
        B_scatter[j + i * num_elements] =
          B[j % num_row + (j / num_row) * 128 + (i % sqrt_p) * num_row + (i / sqrt_p) * num_col * 128];
      }
    }
  }
  // Scatter out the "Adjusted" matrices
  MPI_Scatter(A_scatter, num_elements, MPI_FLOAT, A_chunk, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatter(B_scatter, num_elements, MPI_FLOAT, B_chunk, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);

  
  // Get all processors sync'd
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Status junk;
  
  int sending_to = ((col_color < row_color) ? rank - row_color + sqrt_p : rank - row_color);

  // Initial Shift
  if(row_color == 0); // Do nothing
  /*else if(row_color % 2 == 1) { // Odd rows easily able to shift
    // Even rank sends first
    if(col_color % 2 == 0) {
      MPI_Send(A_chunk, num_elements, MPI_FLOAT, sending_to, 0, MPI_COMM_WORLD);
      MPI_Recv(A_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
    }
    else {
      MPI_Recv(A_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
    	MPI_Send(A_chunk, num_elements, MPI_FLOAT, sending_to, 0, MPI_COMM_WORLD);
    }
		memcpy(A_chunk, A_tmp, sizeof(float) * num_elements);
  }*/
  else { // Even rows shift incrementally by 1 (Could be sped up)
    for(i = 0; i < row_color; i++) {
      if(col_color % 2 == 0) {
				MPI_Send(A_chunk, num_elements, MPI_FLOAT, ((col_color == 0) ? rank - 1 + sqrt_p : rank - 1), 0, MPI_COMM_WORLD);
				MPI_Recv(A_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
      }
      else {
        MPI_Recv(A_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
				MPI_Send(A_chunk, num_elements, MPI_FLOAT, ((col_color == 0) ? rank - 1 + sqrt_p : rank - 1), 0, MPI_COMM_WORLD);
      }
      memcpy(A_chunk, A_tmp, sizeof(float) * num_elements);
    }
  }
  sending_to = rank - col_color * sqrt_p;
  if(sending_to < 0) sending_to += size;
  if(col_color == 0); // Do nothing
  /*else if(col_color % 2 == 1) { // Odd rows easily able to shift
    // Even rank sends first
    if(row_color % 2 == 0) {
      MPI_Send(B_chunk, num_elements, MPI_FLOAT, sending_to, 0, MPI_COMM_WORLD);
      MPI_Recv(B_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
    }
    else {
      MPI_Recv(B_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
    	MPI_Send(B_chunk, num_elements, MPI_FLOAT, sending_to, 0, MPI_COMM_WORLD);
    }
		memcpy(B_chunk, B_tmp, sizeof(float) * num_elements);
  }*/
  else { // Even rows shift incrementally by 1 (Could be sped up)
    for(i = 0; i < col_color; i++) {
			if(row_color % 2 == 0) {
        MPI_Send(B_chunk, num_elements, MPI_FLOAT,((row_color == 0) ? rank - sqrt_p + size : rank - sqrt_p), 0, MPI_COMM_WORLD);
        MPI_Recv(B_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
			}
			else {
        MPI_Recv(A_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
    	  MPI_Send(A_chunk, num_elements, MPI_FLOAT, ((row_color == 0) ? rank - sqrt_p + size : rank - sqrt_p), 0, MPI_COMM_WORLD);
			}
			memcpy(B_chunk, B_tmp, sizeof(float) * num_elements);
    }
  }
  // Calculations
  for(i = 0; i < num_row; i++) {
    for(j = 0; j < num_col; j++) {
      for(k = 0; k < num_row; k++) {
          C_chunk[j + num_col * i] += A_chunk[k + num_col * i] * B_chunk[j + num_col * k];
      }
    }
  }
  // Remaining shifts
  for(n = 1; n < sqrt(size); n++) {
    // Shift A
    sending_to = ((col_color < row_color) ? rank - row_color + sqrt_p : rank - row_color);
		if(col_color % 2 == 0) {
        MPI_Send(A_chunk, num_elements, MPI_FLOAT, ((col_color == 0) ? rank - 1 + sqrt_p : rank - 1), 0, MPI_COMM_WORLD);
        MPI_Recv(A_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
		}
		else {
        MPI_Recv(A_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
    	  MPI_Send(A_chunk, num_elements, MPI_FLOAT, ((col_color == 0) ? rank - 1 + sqrt_p : rank - 1), 0, MPI_COMM_WORLD);
		}
    memcpy(A_chunk, A_tmp, sizeof(float) * num_elements);
    // Shift B
    sending_to = rank - col_color * sqrt_p;
		if(row_color % 2 == 0) {

        MPI_Send(B_chunk, num_elements, MPI_FLOAT,((row_color == 0) ? rank - sqrt_p + size : rank - sqrt_p), 0, MPI_COMM_WORLD);
        MPI_Recv(B_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
		}
		else {
        MPI_Recv(A_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
    	  MPI_Send(A_chunk, num_elements, MPI_FLOAT, ((row_color == 0) ? rank - sqrt_p + size : rank - sqrt_p), 0, MPI_COMM_WORLD);
		}
		memcpy(B_chunk, B_tmp, sizeof(float) * num_elements);
		// More calculations
    for(i = 0; i < num_row; i++) {
      for(j = 0; j < num_col; j++) {
        for(k = 0; k < num_row; k++) {
          C_chunk[j + num_col * i] += A_chunk[k + num_col * i] * B_chunk[j + num_col * k];
        }
      }
    }
  }

  MPI_Gather(C_chunk, num_elements, MPI_FLOAT, C_gather, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);
 	// Shift C back to standard matrix format (This doesn't work, we aren't getting the right numbers at all)
  if(rank == 0) {
    for(i = 0; i < size; i++) {
      for(j = 0; j < num_elements; j++) {
        C[j % num_row + 128 * (j / num_col + (i / sqrt_p) * num_row) + (i % sqrt_p) * num_row] = C_gather[j + i * num_elements];
      }
    }
  }
  if(rank == 0) end_time = MPI_Wtime();
  if(rank == 0) cout << "Elapsed time: " << end_time - start_time << endl;
  // Print output file when desired
  if(print) {
    if(rank == 0) {
      ofstream output;
      output.open("cannon.output");
      for(int i = 0; i < 128; i++) {
        for(int j = 0; j < 128; j++) {
          output << C_gather[j + i * 128] << " ";
        }
        output << endl;
      }
    }
  } 
    
  //Finish
  MPI_Finalize();
  
  return 0;
}
