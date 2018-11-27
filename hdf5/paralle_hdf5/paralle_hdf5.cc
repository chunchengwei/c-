//******************************************************************************
// File Name: paralle_hdf5.cc
// Author: Chuncheng Wei
// Mail: weicc1989@gmail.com
// Created Time : Sun 04 Feb 2018 09:47:16 PM CST
// Last Modified: Wed Nov 28 01:17:50 2018
//******************************************************************************

#include <iostream>
#include "paralle_hdf5.h"

#define RANK 2

using std::string;
using std::vector;
using std::cout;
using std::endl;


void h5_load_param(string fname, string dname,
    double*& data_CP, int *dims_CP, int *dims_TP) {

  /*****************
   * MPI variables *
   *****************/
  int mpi_size, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


  /******************************************************************
   * h5 file reading ...
   ******************************************************************/

  // Set up file access property list with parallel I/O access.
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  // open file
  hid_t file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);

  // open dataset
  hid_t dset_id = H5Dopen(file_id, dname.c_str(), H5P_DEFAULT);

  // get dataspace
  hid_t filespace = H5Dget_space(dset_id);

  // get dataset dimensions
  int ndim = H5Sget_simple_extent_ndims(filespace);
  if (ndim != RANK) {
    cout << "ERROR: only deal RANK=2 dataset" << endl;
    return;
  }
  hsize_t	dims[RANK];
  H5Sget_simple_extent_dims(filespace, dims, NULL);

  /*
   * Each process defines dataset in memory and read it from the hyperslab
   * in the file.
   */
  hsize_t	offset[RANK] = {mpi_rank, 0};
  hsize_t	stride[RANK] = {mpi_size, 1};
  hsize_t	rcount[RANK] = {dims[0] / mpi_size, dims[1]};

  if (mpi_rank < dims[0] % mpi_size) rcount[0] += 1;

  // creat memspace
  hid_t memspace = H5Screate_simple(RANK, rcount, NULL);

  // select_hyperslab
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, rcount, NULL);

  // read data
  data_CP = new double[rcount[0] * rcount[1]];
  H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, data_CP);

  // set dims_XX
  dims_CP[0] = rcount[0];
  dims_CP[1] = rcount[1];
  dims_TP[0] = dims[0];
  dims_TP[1] = dims[1];

  // close
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset_id);
  H5Fclose(file_id);
}

void h5_creat_file(string fname,
    const vector<std::string> & dname, double * data[], const int * n_col,
    const vector<std::string> & ext_dname, const int * ext_n_col) {

  /*****************
   * MPI variables *
   *****************/
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


  /******************************************************************
   * h5 file creating ...
   ******************************************************************/

  // Set up file access property list with parallel I/O access.
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  // creat file
  hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  /***********************
   * dataset writing ...
   ***********************/
  int nd = dname.size();
  hsize_t dims[RANK] = {1};

  for (int i = 0; i < nd; i++) {

    // creat dataspace
    dims[1] = n_col[i];
    hid_t filespace = H5Screate_simple (RANK, dims, NULL);

    // creat dataset
    hid_t dset_id = H5Dcreate(file_id, dname[i].c_str(), H5T_NATIVE_DOUBLE, filespace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // write dataset
    if (mpi_rank == 0)
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
          H5P_DEFAULT, data[i]);

    // close
    H5Sclose(filespace);
    H5Dclose(dset_id);
  }

  /**************************
   * ext_dataset writing ...
   **************************/
  nd = ext_dname.size();
  dims[0] = 0;
  hsize_t maxdims[RANK] = {H5S_UNLIMITED};
  hsize_t chunk_dims[RANK] = {1};

  for (int i = 0; i < nd; i++) {

    // creat dataspace with unlimited dimensions.
    dims[1] = ext_n_col[i];
    maxdims[1] = ext_n_col[i];
    chunk_dims[1] = ext_n_col[i];
    hid_t filespace = H5Screate_simple (RANK, dims, maxdims);

    // set chunk to dataset creat property list.
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, RANK, chunk_dims);

    // creat dataset
    hid_t dset_id = H5Dcreate(file_id, ext_dname[i].c_str(), H5T_NATIVE_DOUBLE, filespace,
        H5P_DEFAULT, plist_id, H5P_DEFAULT);

    // close
    H5Pclose(plist_id);
    H5Sclose(filespace);
    H5Dclose(dset_id);
  }

  // close
  H5Fclose(file_id);
}

void h5_append_ext_data(string fname,
    const vector<string> & dname, double * data[], const int * n_col,
    int n_row_CP, int n_row_TP) {


  /*****************
   * MPI variables *
   *****************/
  int mpi_size, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


  /******************************************************************
   * h5 file reading ...
   ******************************************************************/

  // Set up file access property list with parallel I/O access.
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  // open file
  hid_t file_id = H5Fopen(fname.c_str(), H5F_ACC_RDWR, plist_id);
  H5Pclose(plist_id);

  /***********************
   * dataset writing ...
   ***********************/
  int nd = dname.size();

  for (int i = 0; i < nd; i++) {

    // open dataset
    hid_t dset_id = H5Dopen(file_id, dname[i].c_str(), H5P_DEFAULT);

    // get dataspace
    hid_t filespace = H5Dget_space(dset_id);

    // get dataset dimensions
    hsize_t old_dimsf[RANK];
	  H5Sget_simple_extent_dims(filespace, old_dimsf, NULL);
	  H5Sclose(filespace);

    // extend dataset
	  hsize_t new_dimsf[RANK] = {old_dimsf[0] + n_row_TP, old_dimsf[1]};
	  H5Dset_extent(dset_id, new_dimsf);

    /*
	   * Select a hyperslab in extended portion of dataset.
     */
    hsize_t	offset[RANK] = {mpi_rank + old_dimsf[0], 0};
    hsize_t	stride[RANK] = {mpi_size, 1};
    hsize_t	count[RANK]  = {n_row_CP, n_col[i]};

    // get dataspace
    filespace = H5Dget_space(dset_id);

	  // Create property list for collective dataset write.
	  plist_id = H5Pcreate(H5P_DATASET_XFER);
	  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    if (n_row_CP > 0) {

      // creat memspace
      hid_t memspace = H5Screate_simple(RANK, count, NULL);

      // select_hyperslab
	    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, NULL);

      // write data[i]
	    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data[i]);

      // close
	    H5Sclose(memspace);

    } else {

      // select_hyperslab
      H5Sselect_none(filespace);

      // write data[i]
	    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, filespace, plist_id, data[i]);
    }

    // close
	  H5Pclose(plist_id);
	  H5Sclose(filespace);
	  H5Dclose(dset_id);
  }

  // close
	H5Fclose(file_id);
}

void h5_loop_run(
    std::string pfname, std::string pdname, std::string fname,
    const std::vector<std::string> & dname, double * data[], const int * n_col,
    const std::vector<std::string> & ext_dname, const int * ext_n_col,
    void (*f)(double *Para, int n_par, double * Result[], void *),
    int minibatch, void *context) {

  /*****************
   * MPI variables *
   *****************/
  int mpi_size, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


  /******************************************************************
   * parameter loading ...
   ******************************************************************/
  double * para;
  int pdims[RANK];  // dimensions in current process
  int pdimsf[RANK]; // dimensions in file

  h5_load_param(pfname, pdname, para, pdims, pdimsf);

  /******************************************************************
   * creat output file, and save Ekin
   ******************************************************************/
  h5_creat_file(fname, dname, data, n_col, ext_dname, ext_n_col);

  /******************************************************************
   * looping ...
   * append ext_data to ext_dataset in output file, after dealling
   * minibatch rows of parameter.
   ******************************************************************/
  int payload_remain_CP = pdims[0];   // for current process
  int payload_remain_TP = pdimsf[0];  // for total process
  int steps = pdimsf[0] / (mpi_size * minibatch);

  while (payload_remain_TP) {

    // set payload_CP, payload_TP, which represent payloads in this epochs.
    int payload_CP;
    int payload_TP;

    if (steps > 0) {
      payload_CP = minibatch;
      payload_TP = mpi_size * payload_CP;
    }
    else if (steps == 0) {
      payload_CP = payload_remain_TP / mpi_size;
      payload_TP = mpi_size * payload_CP;
    }
    else {
      payload_CP = (mpi_rank < payload_remain_TP ? 1 : 0);
      payload_TP = payload_remain_TP;
    }
    steps--;

    /*********************
     * calculate ext_data
     *********************/
    int nd = ext_dname.size();
    double * ext_data[nd];
    for (int i = 0; i < nd; i++)
      ext_data[i] = new double[payload_CP * ext_n_col[i]];

    for (int j = 0; j < payload_CP; j++) {

      // pacr: pointer to current para row
      int  row_in_p = j + pdims[0] - payload_remain_CP;
      double * pacr = para + row_in_p * pdims[1];

      // edcr[i] pointer to current ext_data[i] row
      double * edcr[nd];
      for (int i = 0; i < nd; i++)
        edcr[i] = ext_data[i] + j * ext_n_col[i];

      // run f()
      f(pacr, pdims[1], edcr, context);
    }

    /****************************************************************
     * append ext_data to file
     ****************************************************************/
    h5_append_ext_data(fname, ext_dname, ext_data, ext_n_col, payload_CP, payload_TP);

    // free ext_data
    for (int i = 0; i < nd; i++)
      delete[] ext_data[i];

    // update payload_remain_XX
    payload_remain_CP -= payload_CP;
    payload_remain_TP -= payload_TP;
  }

  // free para
  delete[] para;
}
