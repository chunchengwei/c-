//******************************************************************************
// File Name: test.cc
// Author: Chuncheng Wei
// Mail: weicc1989@gmail.com
// Created Time : Sun 04 Feb 2018 10:54:20 PM CST
// Last Modified: Mon 05 Feb 2018 06:26:39 PM CST
//******************************************************************************

#include <iostream>
#include "paralle_hdf5.h"


using namespace std;

int main (int argc, char * argv[])
{
  /*
   * Initialize MPI
   */
  int mpi_size, mpi_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  double * rdata;
  hsize_t rcount[2];
  hsize_t dims  [2];

  h5_load_data("dset.h5", "dset", rdata, rcount, dims);

  cout << "(" << mpi_rank << "/" << mpi_size << ") "
       << rcount[0] << "," << rcount[1] << " | "
       << dims[0]   << "," << dims[1] << endl;

  int sum = 0;
  for (int i=0; i<3; i++) {
    sum = 0;
    for (int j=0; j<rcount[1]; j++)
      sum += rdata[i * rcount[1] + j];
    cout << "(" << mpi_rank << "/" << mpi_size << ") " << i << ": sum = " << sum << endl;
  }

  double ene[3] = {1,3,5};
  h5_creat_file("result.h5", "sum", ene, 3);

	double wdata[5*3];
  hsize_t wcount[2] = {3, 3};
	sum = 0;
	  for (int i=0; i<5; i++) {
	    sum = 0;
	    for (int j=0; j<rcount[1]; j++)
	      sum += rdata[i * rcount[1] + j];
	    wdata[i*3] = sum;
	    wdata[i*3+1] = i;
	    wdata[i*3+2] = mpi_rank;
	  }

  h5_append_data("result.h5", "sum", wdata, wcount, mpi_size * 5);

  return 0;
}

