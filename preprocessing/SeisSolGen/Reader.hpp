/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
*                                                                             *
* This file is part of the SeisSolGen project. For conditions of              *
* distribution and use, please see the copyright notice                       *
* in its root directory.                                                      *
******************************************************************************/
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#ifndef READER_HPP
#define READER_HPP

#include <string>

namespace seissolgen {

  /**
   * Abstract definition of the a file Reader for Matrix Market format for CSC and CSR.
   */
  class Reader {
    public:
      /**
       * reads sparse matrix definition from file into arbitray sparse representations
       * based on row and column indecies.
       *
       * @param tFilename file to open
       * @param ptr_rowidx pointer to row indecies which is set by this routine
       * @param ptr_colidx pointer to col indecies which is set by this routine
       * @param ptr_values pointer to values which is set by this routine
       * @param numRows number of rows which is set by this routine
       * @param numCols number of cols which is set by this routine
       * @param numElems number of elements which is set by this routine
       */
      virtual void parse_file(std::string tFilename, int*& ptr_rowidx, int*& ptr_colidx, double*& ptr_values, int& numRows, int& numCols, int& numElems) = 0;
  };

}

#endif /* READER_HPP */

