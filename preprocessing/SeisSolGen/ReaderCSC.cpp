/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
*                                                                             *
* This file is part of the SeisSolGen project. For conditions of              *
* distribution and use, please see the copyright notice                       *
* in its root directory.                                                      *
******************************************************************************/
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <immintrin.h>

#include "ReaderCSC.hpp"

namespace seissolgen {

  void ReaderCSC::parse_file(std::string tFilename, int*& ptr_rowidx, int*& ptr_colidx, double*& ptr_values, int& numRows, int& numCols, int& numElems) {
    std::ifstream file;
    char a[256];
    int i = 0;

    file.open(tFilename.c_str());

    if (!file.is_open()) {
      throw new std::runtime_error("GeneratorCSC::pharse_file cannot open file");
      return;
    }

    // Determine how many lines of comments
    // are in the header of the file
    file.getline(a, 256);

    while (a[0] == '%') {
      file.getline(a, 256);
      i++;
    }

    file.close();

    // Now skip these lines
    // and read the next line
    // This one contains
    // rows columns No.Elems
    file.open(tFilename.c_str());

    for (int j = 0; j < i; j++)
      file.getline(a, 256);

    file >> numRows >> numCols >> numElems;
    file.close();

#ifdef DEBUG
    std::cout << numRows << " " << numCols << " " << numElems << std::endl;
#endif

    // allocate memory for CSR structure
    ptr_rowidx = (int*)_mm_malloc(sizeof(int) * numElems, 64);
    ptr_values = (double*)_mm_malloc(sizeof(double) * numElems, 64);
    ptr_colidx = (int*)_mm_malloc(sizeof(int) * (numCols + 1), 64);
    int* colidx_id = (int*)_mm_malloc(sizeof(int) * (numCols), 64);

    double* values = ptr_values;
    int* rowidx = ptr_rowidx;
    int* colidx = ptr_colidx;

    memset(values, 0, sizeof(double)*numElems);
    memset(rowidx, 0, sizeof(int)*numElems);
    memset(colidx_id, 0, sizeof(int)*numCols);

    for (int j = 0; j < (numCols + 1); j++)
      colidx[j] = numElems;

#ifdef DEBUG
    std::cout << "Alloc done!" << std::endl;
#endif

    // read element data into data structure
    file.open(tFilename.c_str());

    // skip comment + frist line
    for (int j = 0; j < i + 1; j++)
      file.getline(a, 256);

    i = 0;
    ptr_colidx[0] = 0;
    int last_col = 0;

    while (i < (numElems)) {
      int r, c;
      double val;
      file >> r >> c >> val;
      r--;
      c--;
      rowidx[i] = r;
      values[i] = val;
      i++;
      colidx_id[c] = 1;
      colidx[c + 1] = i;
#ifdef DEBUG
      std::cout << r << " " << c << " " << i << " " << colidx[c + 1] << std::endl;
#endif
    }

    file.close();

    for (int j = 0; j < numCols; j++) {
      if (colidx_id[j] == 0) {
        colidx[j + 1] = colidx[j];
      }
    }

    _mm_free(colidx_id);
#ifdef DEBUG
    std::cout << "Parsing done!" << std::endl;
    std::cout << "Printing CSC Matrix we just read:" << std::endl;

    double* tmp = new double[numRows * numCols];

    for (int l = 0; l < (numRows * numCols); l++)
      tmp[l] = 0.0;

    for (int t = 0; t < numCols; t++) {
      int lcl_colElems = colidx[t + 1] - colidx[t];

      for (int z = 0; z < lcl_colElems; z++) {
        tmp[(t * numRows) + rowidx[colidx[t] + z]] = values[colidx[t] + z];
      }
    }

    for (int l = 0; l < numRows; l++) {
      for (int k = 0; k < numCols; k++) {
        std::cout << tmp[(k * numRows) + l] << " ";
      }

      std::cout << std::endl;
    }

    delete[] tmp;
#endif
  }

}

