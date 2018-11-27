#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;

int main(){
//  float m1[128*128];
//  float m2[128*128];
//  float a[128*128];
  int size, rank;
  double time;
  
  MPI_Init(NULL,NULL);
  
  time = MPI_Wtime();
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

	int row_color = rank / sqrt(size);
	int col_color = rank % sqrt(size);

	int num_row, num_col, i;
	num_col = num_row = (128 / sqrt(size));

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
  	}
  	B = (float**) malloc(sizeof(float*) * 128);
  	for(i = 0; i < num_row; i++) {
  		B[i] = (float*) malloc(sizeof(float) * num_col);
  	}
  	C = (float**) malloc(sizeof(float*) * num_row);
  	for(i = 0; i < num_row; i++) {
		C[i] = (float*) malloc(sizeof(float) * num_col);
	}
/*  for(int i = 0; i < 128; i++){
    m1[i*128+i] = 1;
    m2[i*128+i] = 1;
  }
*/  
  if(rank == 0)
      cout << "Matricies Set" << endl;  
  //Algorithm

  
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
