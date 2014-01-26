/******************************************************************************
* Copyright (C) 2012-2014 Technische Universitaet Muenchen                    *
*                                                                             *
* This file is part of the SeisSolGen project. For conditions of              *
* distribution and use, please see the copyright notice                       *
* in its root directory.                                                      *
******************************************************************************/
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#ifndef GENERATORCSC_HPP
#define GENERATORCSC_HPP

#include <vector>

#include "Generator.hpp"

namespace seissolgen {

  /**
   * FIXME add description
   */
  class GeneratorCSC : public Generator {
    private:
      bool bGenerateExitForCK_;
      bool bAdd_;
      std::vector<int> BasisfunctionsCounter_;

      void generate_code_left_innerloop_scalar(std::stringstream& codestream, int ldc, int l, int z, int* rowidx, int* colidx);

      void generate_code_left_innerloop_2vector(std::stringstream& codestream, int ldc, int l, int z, int* rowidx, int* colidx);

      void generate_code_left_innerloop_3vector(std::stringstream& codestream, int ldc, int l, int z, int* rowidx, int* colidx);

      void generate_code_left_innerloop_4vector(std::stringstream& codestream, int ldc, int l, int z, int* rowidx, int* colidx);

    public:
      GeneratorCSC();

      GeneratorCSC(bool bGenerateExitForCK, int nMaxOrder, bool bAdd);

      virtual std::string generate_code_right(std::string tFilename, int nRows, int nCols, bool bIsRowMajor, int lda, int ldc);

      virtual std::string generate_code_left(std::string tFilename, int nRows, int nCols, bool bIsRowMajor, int ldb, int ldc);
  };

}

#endif /* GENERATORCSC_HPP */

