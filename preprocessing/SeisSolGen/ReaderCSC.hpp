/******************************************************************************
* Copyright (C) 2012-2014 Technische Universitaet Muenchen                    *
*                                                                             *
* This file is part of the SeisSolGen project. For conditions of              *
* distribution and use, please see the copyright notice                       *
* in its root directory.                                                      *
******************************************************************************/
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#ifndef READERCSC_HPP
#define READERCSC_HPP

#include "Reader.hpp"

namespace seissolgen {

  /**
   * FIXME add description
   */
  class ReaderCSC : public Reader {
    public:
      virtual void parse_file(std::string tFilename, int*& ptr_rowidx, int*& ptr_colidx, double*& ptr_values, int& numRows, int& numCols, int& numElems);
  };

}

#endif /* READERCSC_HPP */

