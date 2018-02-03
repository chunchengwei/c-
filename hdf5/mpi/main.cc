//******************************************************************************
// File Name: main.cc
// Author: Chuncheng Wei
// Mail: weicc1989@gmail.com
// Created Time : Fri 02 Feb 2018 10:18:39 PM CST
// Last Modified: Sat 03 Feb 2018 11:40:05 PM CST
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
  hsize_t *dims;  /* dataset dimensions */
  int ndim;

  hid_t	plist_id;                 /* property list identifier */
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

  ndim = H5Sget_simple_extent_ndims(filespace);
  dims = new hsize_t[ndim];
  ndim = H5Sget_simple_extent_dims(filespace, dims, NULL);

  cout << "[" << mpi_rank << "/" << mpi_size << "]: ";
  for (int i=0; i<ndim; i++)
    cout << dims[i] << " ";
  cout << endl;

  /*
   * Each process defines dataset in memory and read it from the hyperslab
   * in the file.
   */
  int n_payload = dims[0] / mpi_size;
  int r_bk = dims[0] % mpi_size;

  hsize_t	count[2];	          /* hyperslab selection parameters */
  count[1] = dims[1];

  hsize_t	offset[2];
  offset[1] = 0;

  if (mpi_rank < r_bk) {
    count[0] = n_payload + 1;
    offset[0] = mpi_rank * count[0];
  } else {
    count[0] = n_payload;
    offset[0] = r_bk * (n_payload + 1) + (mpi_rank - r_bk) * count[0];
  }

  hid_t memspace = H5Screate_simple(2, count, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

  int data[count[0] * count[1]];
  /*
   * Create property list for collective dataset write.
   */
  // plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  // plist_id = H5Pcreate(H5P_DATASET_XFER);
  // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  // status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace,
					 // plist_id, data);
  // H5Pclose(plist_id);

  status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace,
           H5P_DEFAULT, data);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset_id);
  H5Fclose(file_id);


  int sum = 0;
  for (int i=0; i<count[0]; i++) {
    sum = 0;
    for (int j=0; j<count[1]; j++)
      sum += data[i * count[1] + j];
    cout << "(" << mpi_rank << "/" << mpi_size << ") " << i << ": sum = " << sum << endl;
  }

  /*
   * open a new h5 file and write sum to it.
   */
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);

  file_id = H5Fcreate("result.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  /*
   * Create the dataspace for the dataset.
   */
  hsize_t dimsf[2] = {dims[0], 3};                 /* dataset dimensions */
  filespace = H5Screate_simple(2, dimsf, NULL);

  /*
   * Create the dataset with default properties and close filespace.
   */
  dset_id = H5Dcreate(file_id, "sum", H5T_NATIVE_INT, filespace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(filespace);

  hsize_t	wcount[2] = {count[0], 3};	          /* hyperslab selection parameters */
  hsize_t	woffset[2] = {offset[0], 0};
  memspace = H5Screate_simple(2, wcount, NULL);

  /*
   * Select hyperslab in the file.
   */
  filespace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, woffset, NULL, wcount, NULL);

  int wdata[wcount[0] * wcount[1]];
  sum = 0;
  for (int i=0; i<count[0]; i++) {
    sum = 0;
    for (int j=0; j<count[1]; j++)
      sum += data[i * count[1] + j];
    wdata[i*3] = sum;
    wdata[i*3+1] = i;
    wdata[i*3+2] = mpi_rank;
  }

  /*
   * Create property list for collective dataset write.
   */
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
		       plist_id, wdata);
  H5Pclose(plist_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Dclose(dset_id);
  H5Fclose(file_id);

  MPI_Finalize();

  return 0;

}
