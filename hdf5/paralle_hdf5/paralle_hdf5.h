//******************************************************************************
// File Name: paralle_hdf5.h
// Author: Chuncheng Wei
// Mail: weicc1989@gmail.com
// Created Time : Sun 04 Feb 2018 09:47:22 PM CST
// Last Modified: Mon 05 Feb 2018 05:13:08 PM CST
//******************************************************************************

#ifndef _PARALLE_HDF5_H
#define _PARALLE_HDF5_H

#include "hdf5.h"
#include <vector>


/******************************************************************
 * load data from h5 file.
 *
 * current process:
 *    rdata []:     save data for each process.
 *    rcount[RANK]: dimensions of rdata.
 *
 * all process:
 *    dims  [RANK]: dimensions of data in file.
 ******************************************************************/
void h5_load_data(std::string fname, std::string dname,
    double*& rdata, hsize_t *rcount, hsize_t *dims);

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
void h5_creat_file(std::string fname, std::string dname,
    double *ene, int n_col);

/******************************************************************
 * append data to h5 file.
 *
 * current process:
 *    data[]:      raw data (length = count[0] * count[1]).
 *    count[RANK]: elements to be write in each dimension.
 *
 * all process:
 *    ntot_row:    total rows to be append by all process
 ******************************************************************/
void h5_append_data(std::string fname, std::string dname,
    double *data, hsize_t *count, int ntot_row);

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
 ******************************************************************/
void h5_loop_run(
    std::string pfname, std::string pdname,
    std::string ofname, std::string odname,
    void (*f)(double *Para, int n_par, double *Result, void *),
    std::vector<double> &Ekin, int minibatch, void *context);

#endif
