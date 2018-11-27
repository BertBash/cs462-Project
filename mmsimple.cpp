#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;

int main(){
  float m1[128*128];
  float m2[128*128];
  float a[128*128];
  int size, rank;
  double time;
  
  MPI_Init(NULL,NULL);
  
  time = MPI_Wtime();
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

	int row_comm = rank / sqrt(size);

  //Generate Matricies
  for(int i = 0; i < 128; i++){
    m1[i*128+i] = 1;
    m2[i*128+i] = 1;
  }
  
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
