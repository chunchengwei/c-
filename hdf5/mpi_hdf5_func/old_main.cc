//******************************************************************************
// File Name: main.cc
// Author: Chuncheng Wei
// Mail: weicc1989@gmail.com
// Created Time : Fri 02 Feb 2018 10:18:39 PM CST
// Last Modified: Sun 04 Feb 2018 05:36:15 PM CST
//******************************************************************************

#include "hdf5.h"
#include <iostream>

using namespace std;

int main (int argc, char * argv[])
{
  /*
   * HDF5 APIs definitions
   */
  hid_t   file_id, dset_id; /* file and dataset identifiers */
  hid_t   filespace;        /* file and memory filespace identifiers */
  hid_t	  plist_id;         /* property list identifier */
  herr_t	status;

  /*
   * MPI variables
   */
  int mpi_size, mpi_rank;
  MPI_Comm comm  = MPI_COMM_WORLD;
  MPI_Info info  = MPI_INFO_NULL;

  /*
   * Initialize MPI
   */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  /*
   * Set up file access property list with parallel I/O access
   */
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);
  file_id = H5Fopen("dset.h5", H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);

  /*
   * Create property list for collective dataset read.
   */
  dset_id = H5Dopen(file_id, "dset", H5P_DEFAULT);
  filespace = H5Dget_space(dset_id);

  int ndim = H5Sget_simple_extent_ndims(filespace);
  hsize_t dims[ndim];
  ndim = H5Sget_simple_extent_dims(filespace, dims, NULL);

  /*
   * Each process defines dataset in memory and read it from the hyperslab
   * in the file.
   */
  int n_payload = dims[0] / mpi_size;

  hsize_t	offset[2] = {mpi_rank, 0};
  hsize_t	stride[2] = {mpi_size, 1};
  hsize_t	count [2] = {n_payload, dims[1]};

  if (mpi_rank < dims[0] % mpi_size)
    count[0] = n_payload + 1;

  hid_t memspace = H5Screate_simple(2, count, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, NULL);

  int data[count[0] * count[1]];

  status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace,
           H5P_DEFAULT, data);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset_id);
  H5Fclose(file_id);

  int sum = 0;
  for (int i=0; i<3; i++) {
    sum = 0;
    for (int j=0; j<count[1]; j++)
      sum += data[i * count[1] + j];
    cout << "(" << mpi_rank << "/" << mpi_size << ") " << i << ": sum = " << sum << endl;
  }


  /****************************************
   * creat unlimit dataset with new file. *
   ****************************************/
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);
  file_id = H5Fcreate("result.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  /* Create the data space with unlimited dimensions. */
  hsize_t dimsf[2]       = {0, 3};
  hsize_t maxdimsf[2]    = {H5S_UNLIMITED, 3};
  hsize_t chunk_dimsf[2] = {1, 3};
  filespace = H5Screate_simple (2, dimsf, maxdimsf);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(plist_id, 2, chunk_dimsf);
  dset_id = H5Dcreate(file_id, "sum", H5T_NATIVE_INT, filespace,
                      H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Dclose(dset_id);
  H5Fclose(file_id);






  /*
   * raw data: data
   * dimansion: count[2]
   * loop the with minibatch
   */

  int minibatch = 5;
  int payload_remain = count[0];

  while (payload_remain) {

    // set n_row to write.
    int n_row, n_row_all_process;
    if (payload_remain < minibatch) {
      n_row = payload_remain;
      n_row_all_process = dims[0] % (mpi_size * minibatch);
    } else {
      n_row = minibatch;
      n_row_all_process = mpi_size * n_row;
    }

	  hsize_t old_dimsf[2];
	  hsize_t new_dimsf[2];
	  /*************************************
	   * reopen file and write data to it. *
	   *************************************/
	  plist_id = H5Pcreate(H5P_FILE_ACCESS);
	  H5Pset_fapl_mpio(plist_id, comm, info);
	  file_id = H5Fopen("result.h5", H5F_ACC_RDWR, plist_id);
	  H5Pclose(plist_id);

	  dset_id = H5Dopen(file_id, "sum", H5P_DEFAULT);
	  filespace = H5Dget_space(dset_id);
	  H5Sget_simple_extent_dims(filespace, old_dimsf, NULL);
	  H5Sclose(filespace);

	  new_dimsf[0] = old_dimsf[0] + n_row_all_process;
	  new_dimsf[1] = old_dimsf[1];
	  status = H5Dset_extent(dset_id, new_dimsf);

	  /* Select a hyperslab in extended portion of dataset  */
	  // hsize_t	stride[2] = {mpi_size, 1};
	  hsize_t	wcount [2] = {n_row, new_dimsf[1]};
	  offset[0] = mpi_rank + old_dimsf[0];

	  memspace = H5Screate_simple(2, wcount, NULL);
	  filespace = H5Dget_space(dset_id);
	  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, wcount, NULL);

	  int wdata[wcount[0] * wcount[1]];
	  sum = 0;
	  for (int i=0; i<wcount[0]; i++) {
	    sum = 0;
	    for (int j=0; j<count[1]; j++)
	      sum += data[(i+count[0]-payload_remain) * count[1] + j];
	    wdata[i*3] = sum;
	    wdata[i*3+1] = i;
	    wdata[i*3+2] = mpi_rank;
	  }

	  /* Create property list for collective dataset write. */
	  plist_id = H5Pcreate(H5P_DATASET_XFER);
	  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	  status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
	           plist_id, wdata);
	  H5Pclose(plist_id);

	  H5Sclose(memspace);
	  H5Sclose(filespace);
	  H5Dclose(dset_id);
	  H5Fclose(file_id);

    payload_remain = payload_remain < minibatch ? 0 : payload_remain - minibatch;
  }


  MPI_Finalize();

  return 0;

}
