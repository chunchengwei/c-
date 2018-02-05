//******************************************************************************
// File Name: main.cc
// Author: Chuncheng Wei
// Mail: weicc1989@gmail.com
// Created Time : Sun 04 Feb 2018 10:54:20 PM CST
// Last Modified: Mon 05 Feb 2018 11:09:37 PM CST
//******************************************************************************

#include "paralle_hdf5.h"
#include <iostream>


using namespace std;

void sum(double *Para, int n_par, double *Result, void *context)
{
  // mpi set
  int mpi_size, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  double sum = 0;
  for (int i=0; i<n_par; i++)
    sum += Para[i];

  Result[0] = sum;
  Result[1] = mpi_size;
  Result[2] = mpi_rank;
}



int main (int argc, char * argv[])
{
  vector<double> Ekin = {1,3,5};

  MPI_Init(&argc, &argv);

  h5_loop_run("dset.h5", "dset", "result.h5", "sum", sum, Ekin, 5, NULL);

  MPI_Finalize();

  return 0;
}

