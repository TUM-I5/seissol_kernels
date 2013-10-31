/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
*                                                                             *
* This file is part of the SeisSolGen project. For conditions of              *
* distribution and use, please see the copyright notice                       *
* in its root directory.                                                      *
******************************************************************************/
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <string>

namespace seissolgen {

  /**
   * Abstract definition of the a code generator for sparse-dense and dense-sparse
   * matrix multiplications.
   */
  class Generator {
    public:

      /**
       * Generates dense x sparse multiplication code
       *
       * @param tFilename file that contains sparse matrix definition
       * @param nRows number of rows of the left dense matrix A
       * @param nCols number of cols of the left dense matrix A
       * @param bIsRowMajor specifies whethter A and C are row-major or column-major
       * @param lda lda of matrix A
       * @param ldc ldc of matrix C
       */
      virtual std::string generate_code_right(std::string tFilename, int nRows, int nCols, bool bIsRowMajor, int lda, int ldc) = 0;

      /**
       * Generates sparse x dense multiplication code
       *
       * @param tFilename file that contains sparse matrix definition
       * @param nRows number of rows of the right dense matrix B
       * @param nCols number of cols of the right dense matrix B
       * @param bIsRowMajor specifies whethter B and C are row-major or column-major
       * @param ldb ldb of matrix B
       * @param ldc ldc of matrix C
       */
      virtual std::string generate_code_left(std::string tFilename, int nRows, int nCols, bool bIsRowMajor, int ldb, int ldc) = 0;
  };

}

#endif /* GENERATOR_HPP */

