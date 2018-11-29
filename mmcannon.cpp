#include <iostream>
#include <cmath>
#include <mpi.h>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cstring>

using namespace std;

int main(){
  float *a, *b, *c;
  int size, rank;
  int row_size, row_rank;
  int col_size, col_rank;
  double start_time, end_time;

//  int dim[1], period[1], reorder;
//  int up, down, right, left;
//  MPI_Comm row_cart_comm, col_cart_comm;
  
  MPI_Init(NULL,NULL);

  float *A_tmp, *B_tmp;    // Initialize data
  float *A_chunk;              // Entire row chunk beginning to end
  float *B_chunk;              // Entire column chunk beginning to end
  float *C_chunk;              // Subsection result
  float *C; // Final matrix
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(rank == 0){
    a = (float *) malloc(sizeof(float)*128*128);
    b = (float *) malloc(sizeof(float)*128*128);
    c = (float *) malloc(sizeof(float)*128*128);
      
    srand(-1337); 
    //Popuate the matrixes with random numbers between 0 and 1.
    for(int i = 0; i < 128; i++){
      for(int j = 0; j <128; j++){
        a[i*128+j] = (float) (rand() / ((float) RAND_MAX / 2)) - 1;
        b[i*128+j] = (float) (rand() / ((float) RAND_MAX / 2)) - 1;
        c[i*128+j] = 0;
      }
    }
  }
  
  int sqrt_p = sqrt(size);
  int row_color = rank / sqrt_p;
  int col_color = rank % sqrt_p;
  int num_row, num_col, i, j, k, n, m, num_elements;
  num_col = num_row = (128 / sqrt_p);
  srand(time(NULL)+rank);
  
  MPI_Comm row_comm;
  MPI_Comm col_comm;
  
  MPI_Comm_split(MPI_COMM_WORLD, row_color, rank, &row_comm);
  MPI_Comm_split(MPI_COMM_WORLD, col_color, rank, &col_comm);
  

  MPI_Comm_rank(row_comm, &row_rank);
  MPI_Comm_rank(col_comm, &col_rank);
  MPI_Comm_size(row_comm, &row_size);
  MPI_Comm_size(col_comm, &col_size);

//  dim[0] = sqrt(size);
//  period[0]=1;
//  reorder=0;

//  MPI_Cart_create(row_comm, 1, dim, period, reorder, &row_cart_comm);
//  MPI_Cart_create(col_comm, 1, dim, period, reorder, &col_cart_comm);

  //Generate Matricies
  A_chunk = (float*) malloc(sizeof(float) * num_row * num_col);
  B_chunk = (float*) malloc(sizeof(float) * num_row * num_col);
  A_tmp = (float*) malloc(sizeof(float) * num_row * num_col);
  B_tmp = (float*) malloc(sizeof(float) * num_row * num_col);
  C_chunk = (float*) malloc(sizeof(float) * num_row * num_col);

//  A_start = (float*) malloc(sizeof(float) * num_row * num_col);
//  B_start = (float*) malloc(sizeof(float) * num_row * num_col);
  if(rank==0) C = (float*) malloc(sizeof(float) * 128 * 128);
 
  for(i = 0; i < num_row * num_col; i++) {
    A_chunk[i] = (float) (rand()/ ((float) RAND_MAX / 2)) - 1;
    B_chunk[i] = (float) (rand()/ ((float) RAND_MAX / 2)) - 1;
    C_chunk[i] = 0;
  }
  num_elements = num_row * num_col;
  
  // Get all processors sync'd
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Status junk;
  
  // Grab timestamp and continue with communication and number crunching.
  start_time = MPI_Wtime();
  int sending_to = ((col_color < row_color) ? rank - row_color + sqrt_p : rank - row_color);

  // Initial Shift
  if(row_color == 0); // Do nothing
  if(row_color % 2 == 1) { // Odd rows easily able to shift
    // Even rank sends first
    if(col_color % 2 == 0) {
      MPI_Send(A_chunk, num_elements, MPI_FLOAT, sending_to, 0, MPI_COMM_WORLD);
      MPI_Recv(A_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
    }
    else {
      MPI_Recv(A_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
    	MPI_Send(A_chunk, num_elements, MPI_FLOAT, sending_to, 0, MPI_COMM_WORLD);
    }
  }
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
  if(col_color % 2 == 1) { // Odd rows easily able to shift
    // Even rank sends first
    if(row_color % 2 == 0) {
      MPI_Send(B_chunk, num_elements, MPI_FLOAT, sending_to, 0, MPI_COMM_WORLD);
      MPI_Recv(B_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
    }
    else {
      MPI_Recv(B_tmp, num_elements, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &junk);
    	MPI_Send(B_chunk, num_elements, MPI_FLOAT, sending_to, 0, MPI_COMM_WORLD);
    }
  }
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
    for(i = 0; i < num_row; i++) {
      for(j = 0; j < num_col; j++) {
        for(k = 0; k < num_row; k++) {
          C_chunk[j + num_col * i] += A_chunk[k + num_col * i] * B_chunk[j + num_col * k];
        }
      }
    }
  }

//  MPI_Cart_shift(row_cart_comm, 1, -1 * row_color, 
  
  //Communicate
 
  // Each process is sending their A and B (horizontal for A, vertical for B
/*  MPI_Allgather(A_start, num_row * num_col, MPI_FLOAT, A_chunk, num_row * num_col, MPI_FLOAT, row_comm);
  MPI_Allgather(B_start, num_row * num_col, MPI_FLOAT, B_chunk, num_col * num_row, MPI_FLOAT, col_comm);

  //Crunch Numbers  
  num_elements = num_row * num_col;
  for(i = num_row * row_rank, n = 0; i < num_row + num_row * row_rank; i++, n++) {
    for(j = num_col * col_rank, m = 0; j < num_col + num_col * col_rank; j++, m++) {
      C_chunk[m+num_col*n] = 0;
      for(k = 0; k < 128; k++) {
        C_chunk[m+num_col*n] += A_chunk[num_elements * (k / num_row) + i *	num_row + k	% num_row] 
                              * B_chunk[num_elements * (k / num_col) + j * num_col + k % num_col];
		}
	 }
  }
*/
  MPI_Gather(C_chunk, num_elements, MPI_FLOAT, C, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);
  end_time = MPI_Wtime();
  if(rank == 0) cout << "Elapsed time: " << end_time - start_time << endl;
  
    
  //Finish
  MPI_Finalize();
  
  return 0;
}
