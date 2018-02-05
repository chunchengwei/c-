//******************************************************************************
// File Name: paralle_hdf5.cc
// Author: Chuncheng Wei
// Mail: weicc1989@gmail.com
// Created Time : Sun 04 Feb 2018 09:47:16 PM CST
// Last Modified: Mon 05 Feb 2018 05:30:30 PM CST
//******************************************************************************

#include <iostream>
#include "paralle_hdf5.h"

#define RANK 2

using std::string;
using std::vector;
using std::cout;
using std::endl;


void h5_load_data(string fname, string dname,
    double*& rdata, hsize_t *rcount, hsize_t *dims) {

  /******************************************************************
   * load data from h5 file.
   *
   * current process:
   *    rdata []:     save data for each process.
   *    rcount[RANK]: dimensions of rdata.
   *
   * all process:
   *    dims  [RANK]: dimensions of data in file.
   *
   * MPI: mpi_size, mpi_rank.
   ******************************************************************/



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
  H5Sget_simple_extent_dims(filespace, dims, NULL);

  /*
   * Each process defines dataset in memory and read it from the hyperslab
   * in the file.
   */
  hsize_t	offset[RANK] = {mpi_rank, 0};
  hsize_t	stride[RANK] = {mpi_size, 1};
  rcount[0] = dims[0] / mpi_size;
  rcount[1] = dims[1];

  if (mpi_rank < dims[0] % mpi_size) rcount[0] += 1;

  // creat memspace
  hid_t memspace = H5Screate_simple(RANK, rcount, NULL);

  // select_hyperslab
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, rcount, NULL);

  // read data
  rdata = new double[rcount[0] * rcount[1]];

  H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, rdata);

  // close
  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset_id);
  H5Fclose(file_id);
}

void h5_creat_file(string fname, string dname,
    double *ene, int n_col) {

  /******************************************************************
   * creat h5 file with unlimit dataset
   * file name:     fname
   *
   * dataset name:  dname
   * dataset shape: [H5S_UNLIMITED, n_col]
   *
   * dataset name:  "Ekin"
   * dataset shape: [1,             n_col]
   ******************************************************************/



  /*****************
   * MPI variables *
   *****************/
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


  // Set up file access property list with parallel I/O access.
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  // creat file
  hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);


  /*
   * Create the dataspace with unlimited dimensions.
   */
  hsize_t dims[RANK]       = {0,             n_col};
  hsize_t maxdims[RANK]    = {H5S_UNLIMITED, n_col};
  hsize_t chunk_dims[RANK] = {1,             n_col};

  hid_t filespace = H5Screate_simple (RANK, dims, maxdims);

  // set chunk to dataset creat property list.
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(plist_id, RANK, chunk_dims);

  // creat dataset
  hid_t dset_id = H5Dcreate(file_id, dname.c_str(), H5T_NATIVE_DOUBLE, filespace,
      H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  // close
  H5Sclose(filespace);
  H5Dclose(dset_id);


  /*
   * Create the "Ekin" dataspace.
   */
  dims[0] = 1;
  filespace = H5Screate_simple (RANK, dims, NULL);

  // creat dataset
  dset_id = H5Dcreate(file_id, "Ekin", H5T_NATIVE_DOUBLE, filespace,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(filespace);

  // write dataset
  if (mpi_rank == 0)
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, ene);

  // close
  H5Dclose(dset_id);
  H5Fclose(file_id);
}

void h5_append_data(string fname, string dname,
    double *data, hsize_t *count, int ntot_row) {

  /******************************************************************
   * append data to h5 file.
   *
   * current process:
   *    data[]:      raw data (length = count[0] * count[1]).
   *    count[RANK]: elements to be write in each dimension.
   *
   * all process:
   *    ntot_row:    total rows to be append by all process
   *
   * MPI: mpi_size, mpi_rank.
   ******************************************************************/



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

  // open dataset
  hid_t dset_id = H5Dopen(file_id, dname.c_str(), H5P_DEFAULT);

  // get dataspace
  hid_t filespace = H5Dget_space(dset_id);

  // get dataset dimensions
  hsize_t old_dimsf[RANK];
	H5Sget_simple_extent_dims(filespace, old_dimsf, NULL);
	H5Sclose(filespace);

  // extend dataset
	hsize_t new_dimsf[RANK] = {old_dimsf[0] + ntot_row, old_dimsf[1]};
	H5Dset_extent(dset_id, new_dimsf);

  /*
	 * Select a hyperslab in extended portion of dataset.
   */
  hsize_t	offset[RANK] = {mpi_rank + old_dimsf[0], 0};
  hsize_t	stride[RANK] = {mpi_size, 1};

  // creat memspace
  hid_t memspace = H5Screate_simple(RANK, count, NULL);

  // get dataspace
  filespace = H5Dget_space(dset_id);

  // select_hyperslab
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, NULL);

	// Create property list for collective dataset write.
	plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  // write data
	H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
      plist_id, data);
	H5Pclose(plist_id);

  // close
	H5Sclose(memspace);
	H5Sclose(filespace);
	H5Dclose(dset_id);
	H5Fclose(file_id);
}

void h5_loop_run(string pfname, string pdname, string ofname, string odname,
    void (*f)(double *Para, int n_par, double *Result, void *),
    vector<double> &Ekin, int minibatch, void *context) {

  /******************************************************************
   * h5 file contain parameter matrix
   *    pfname: h5 file name
   *    pdname: dataset name
   *
   * h5 file save output
   *    ofname: h5 file name
   *    odname: dataset name
   *
   * void f(Para, n_par, Result, void*)
   *    Para:    parameter array pointer
   *    n_par:   parameter size
   *    Result:  save output
   *    context: any additional information user wants to pass
   *
   * Ekin:
   *
   * minibatch: each process save output to h5 file, after
   *            dealling minibatch rows of parameter.
   *
   * MPI: mpi_size, mpi_rank.
   ******************************************************************/



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
  hsize_t pdims[RANK];  // dimensions in current process
  hsize_t pdimsf[RANK]; // dimensions in file

  h5_load_data(pfname, pdname, para, pdims, pdimsf);

  /******************************************************************
   * creat output file, and save Ekin
   ******************************************************************/
  h5_creat_file(ofname, odname, Ekin.data(), Ekin.size());

  /******************************************************************
   * looping ...
   * append output to dataset in output file, after dealling
   * minibatch rows of parameter.
   ******************************************************************/
  int payload_remain = pdims[0];  // for current process

  while (payload_remain) {

    // set n_payload, ntot_payload
    int n_payload;                // current process
    int ntot_payload;             // total process
    if (payload_remain < minibatch) {
      n_payload       = payload_remain;
      ntot_payload    = pdimsf[0] % (mpi_size * minibatch);
      payload_remain  = 0;
    } else {
      n_payload       = minibatch;
      ntot_payload    = mpi_size * minibatch;
      payload_remain -= minibatch;
    }

    // calculate output
    hsize_t odims[RANK] = {n_payload, Ekin.size()};
    double output[odims[0] * odims[1]];
    for (int i=0; i<n_payload; i++) {
      int  row_in_p = i + pdims[0] - payload_remain - n_payload;
      double * pacr = para   + row_in_p * pdims[1];   // pointer to current para row
      double * opcr = output +        i * odims[1];   // pointer to current output row
      f(pacr, pdims[1], opcr, context);
    }

    /****************************************************************
     * append output to file
     ****************************************************************/
    h5_append_data(ofname, odname, output, odims, ntot_payload);
  }
}
