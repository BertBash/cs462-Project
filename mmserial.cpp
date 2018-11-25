#include <iostream>

using namespace std;

int main(){
  float m1[128][128];
  float m2[128][128];
  float a[128][128];
  
  //Popuate the matrixes with random numbers between 0 and 1.
  for(int i = 0; i < 128; i++){
    for(int j = 0; j <128; j++){
      m1[i][j] = 0;
      m2[i][j] = 0;
      a[i][j] = 0;
    }
    m1[i][i] = 1;
    m2[i][i] = 1;
  }
  
  //Matrix multiply
  for(int i = 0; i < 128; i++){
      for(int j = 0; j < 128; j++){
	for(int k = 0; k < 128; k++){
	  a[i][j] += m1[i][k] * m2[k][j];
	}
      }
  }
  
  //Test
  for(int i = 0; i < 128; i++)
    for(int j = 0; j < 128; j++)
      if( i == j && a[i][j] != 1)
	cout << "Failure." << a[i][j] << endl;
      else if(i != j && a[i][j] != 0)
	cout << "No way." << a[i][j] << endl;

  
  return 0;
}
