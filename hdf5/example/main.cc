//******************************************************************************
// File Name: main.cc
// Author: Chuncheng Wei
// Mail: weicc1989@gmail.com
// Created Time : Sun 04 Feb 2018 10:54:20 PM CST
// Last Modified: Tue 27 Mar 2018 11:05:42 AM CST
//******************************************************************************

#include <iostream>
#include <vector>
#include <string>

#include "paralle_hdf5.h"

#define RANK 2

using namespace std;

void sum(double *Para, int n_par, double *Result[], void *context)
{
  // mpi set
  int mpi_size, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  double sum = 0;
  for (int i=0; i<n_par; i++)
    sum += Para[i];

  Result[0][0] = sum;
  Result[0][1] = mpi_size;
  Result[0][2] = mpi_rank;

  double times = 1;
  for (int i=0; i<n_par; i++)
    times *= Para[i];

  Result[1][0] = times;
  Result[1][1] = mpi_size;
}



int main (int argc, char * argv[])
{
  // data set
  vector<string> dname = {"Ekin", "Ekin2"};
  double Ekin[] = {1, 3, 5};
  double Ekin2[] = {1, 3, 5};
  double * data[] = {Ekin, Ekin2};
  int n_col[] = {3, 3};

  // ext data_set
  vector<string> ext_dname = {"Sum", "Times"};
  // double extd1[] = {1,2,3,4,5,6};
  // double extd2[] = {7,8,9,0};
  // double * ext_data[] = {extd1, extd2};
  int ext_n_col[] = {3, 2};


  // MPI_Init(&argc, &argv);
  // int mpi_size, mpi_rank;
  // MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  // MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  /******************************************************************
   * parameter loading ...
   ******************************************************************/
  // double * para;
  // int pdims[RANK];  // dimensions in current process
  // int pdimsf[RANK]; // dimensions in file

  // h5_load_param("dset.h5", "dset", para, pdims, pdimsf);

  /******************************************************************
   * creat output file, and save Ekin
   ******************************************************************/
  // h5_creat_file("result.h5", dname, data, n_col, ext_dname, ext_n_col);


  // h5_append_ext_data("result.h5", ext_dname, ext_data, ext_n_col, 2, 2*mpi_size);



  // MPI_Finalize();


  // MPI & run h5_loop_run()
  MPI_Init(&argc, &argv);

  h5_loop_run("dset.h5", "dset", "result.h5", dname, data, n_col, ext_dname, ext_n_col, sum,
      2, NULL);

  MPI_Finalize();

  return 0;
}

