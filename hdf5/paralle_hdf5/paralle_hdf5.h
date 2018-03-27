//******************************************************************************
// File Name: paralle_hdf5.h
// Author: Chuncheng Wei
// Mail: weicc1989@gmail.com
// Created Time : Sun 04 Feb 2018 09:47:22 PM CST
// Last Modified: Sat 24 Mar 2018 04:41:59 PM CST
//******************************************************************************

#ifndef _PARALLE_HDF5_H
#define _PARALLE_HDF5_H

#include "hdf5.h"
#include <vector>
#include <string>


/******************************************************************
 * global naming rule:
 *
 *    xx_CP: for current process
 *    xx_TP: for total process
 *
 ******************************************************************/


/******************************************************************
 * load data from h5 file.
 *
 *    data_CP[]:     save param, read from h5 file
 *    dims_CP[RANK]: dimensions of data_CP.
 *    dims_TP[RANK]: dimensions of data in file.
 *
 * ps. when not use data_CP, run "delete[] data_CP;"
 *
 ******************************************************************/
void h5_load_param(std::string fname, std::string dname,
    double*& data_CP, int *dims_CP, int *dims_TP);

/******************************************************************
 * creat h5 file with unlimit dataset
 * file name:     fname
 *
 * dataset name:  dname
 * dataset shape: [1,             n_col]
 *
 * dataset name:  ext_dname
 * dataset shape: [H5S_UNLIMITED, ext_n_col]
 ******************************************************************/
void h5_creat_file(std::string fname,
    const std::vector<std::string> & dname, double * data[], const int * n_col,
    const std::vector<std::string> & ext_dname, const int * ext_n_col);

/******************************************************************
 * append data to h5 file.
 *
 *    data[i][]:  data write to h5 file. (length = n_row_CP * n_col[i]).
 *
 *    n_row_CP, n_col[i]: dimension of data[i][].
 *    n_row_TP:   total rows to be append by all process
 *
 ******************************************************************/
void h5_append_ext_data(std::string fname,
    const std::vector<std::string> & dname, double * data[], const int * n_col,
    int n_row_CP, int n_row_TP);

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
 * minibatch: each process save output to h5 file, after
 *            dealling minibatch rows of parameter.
 ******************************************************************/
void h5_loop_run(
    std::string pfname, std::string pdname, std::string fname,
    const std::vector<std::string> & dname, double * data[], const int * n_col,
    const std::vector<std::string> & ext_dname, const int * ext_n_col,
    void (*f)(double *Para, int n_par, double * Result[], void *),
    int minibatch, void *context);

#endif
