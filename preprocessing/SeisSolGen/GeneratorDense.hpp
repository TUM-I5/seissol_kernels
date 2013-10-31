/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
*                                                                             *
* This file is part of the SeisSolGen project. For conditions of              *
* distribution and use, please see the copyright notice                       *
* in its root directory.                                                      *
******************************************************************************/
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#ifndef GENERATORDENSE_HPP
#define GENERATORDENSE_HPP

namespace seissolgen {

  /**
   * This class generates a dense matrix multiply kernel
   */
  class GeneratorDense {
    private:
      bool bGenerateExitForCK_;
      bool bAdd_;

    public:
      GeneratorDense();

      GeneratorDense(bool bGenerateExitForCK, bool bAdd);

      std::string generate_dense(bool bIsColMajor, int M, int N, int K, int lda, int ldb, int ldc);
  };

}

#endif /* GENERATORDENSE_HPP */

