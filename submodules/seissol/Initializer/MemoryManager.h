/** @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
 *
 * According to good scientific practice, publications on results achieved in whole or in part due to this software should cite at least one paper or referring to an URL presenting this software.
 *
 * The owner wishes to make the software available to all users to use, reproduce, modify, distribute and redistribute also for commercial purposes under the following conditions of the original BSD license. Linking this software module statically or dynamically with other modules is making a combined work based on this software. Thus, the terms and conditions of this license cover the whole combination. As a special exception, the copyright holders of this software give you permission to link it with independent modules or to instantiate templates and macros from this software's source files to produce an executable, regardless of the license terms of these independent modules, and to copy and distribute the resulting executable under terms of your choice, provided that you also meet, for each linked independent module, the terms and conditions of this license of that module.
 *
 * Copyright (c) 2013
 * Technische Universitaet Muenchen
 * Department of Informatics
 * Chair of Scientific Computing
 * http://www5.in.tum.de/
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * All advertising materials mentioning features or use of this software must display the following acknowledgement: This product includes software developed by the Technische Universitaet Muenchen (TUM), Germany, and its contributors.
 * Neither the name of the Technische Universitaet Muenchen, Munich, Germany nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Memory management of SeisSol.
 **/

#ifndef MEMORYMANAGER_H_
#define MEMORYMANAGER_H_

#include <utils/logger.h>

#include "XmlParser.hpp"
#include "MemoryAllocator.h"

namespace seissol {
  namespace initializers {
    class MemoryManager;
  }
}

/**
 * Memory manager of SeisSol.
 **/
class seissol::initializers::MemoryManager {
  private: // explicit private for unit tests
    //! memory allocator
    seissol::MemoryAllocator m_memoryAllocator;

    /**
     * Addresses of the global flux matrices (multiplied by the inverse diagonal mass matrix):
     *
     *    0:  \f$ M^{-1} F^{-, 1} \f$
     *    1 : \f$ M^{-1} F^{-, 2} \f$
     *    2:  \f$ M^{-1} F^{-, 3} \f$
     *    3 : \f$ M^{-1} F^{-, 4} \f$
     *    4:  \f$ M^{-1} F^+{+, 1, 1, 1} \f$
     *    5:  \f$ M^{-1} F^+{+, 1, 1, 2} \f$
     *    6:  \f$ M^{-1} F^+{+, 1, 1, 3} \f$
     *    7:  \f$ M^{-1} F^+{+, 1, 2, 1} \f$
     *    8:  \f$ M^{-1} F^+{+, 1, 1, 2} \f$
     *    9:  \f$ M^{-1} F^+{+, 1, 1, 3} \f$
     *    [..]
     *    51: \f$ M^{-1} F^+{+, 4, 4, 3} \f$
     *    52: \f$ N_{k,i} A_k^+ N_{k,i}^{-1}\f$ or \f$ N_{k,i} A_{k(i)}^- N_{k,i}^{-1} \f$
     *    53: \f$ M^{-1} K^\xi \f$ (for prefetches)
     *
     *   Remark: The ordering of the pointers is given as above, however the chunks in memory are allowed to have a different orderning as given in the XML-file.
     **/ 
    double** m_fluxMatrixPointers;

    /**
     * Addresses of the global stiffness matrices (multiplied by the inverse diagonal mass matrix):
     *
     *    0:  \f$ M^{-1} K^\xi \f$
     *    1:  \f$ M^{-1} K^\eta \f$
     *    2:  \f$ M^{-1} K^\zeta f$
     *    3:  \f$ M^{-1} ( K^\xi )^T \f$
     *    4:  \f$ M^{-1} ( K^\eta )^T \f$
     *    5:  \f$ M^{-1} ( K^\zeta )^T \f$
     *    6:  \f$ M^{-1} F^{-, 1} \f$ (for prefetches)
     *
     *   Remark: The ordering of the pointers is identical to the ordering of the memory chunks (except for the additional flux matrix).
     **/ 
    double** m_stiffnessMatrixPointers;

    /**
     * Initializes a global matrix with the given values.
     *
     * @param i_sparse true if sparse, false if dense.
     * @param i_numberOfRows number of rows.
     * @param i_numberOfColumns number of columns.
     * @param i_rows rows in sparse coordinate format.
     * @param i_columns columns in sparse coordinate format.
     * @param i_values values in sparse coordinate format.
     * @param io_matrix global matrix, which values are set.
     **/
    void initializeGlobalMatrix(       bool                       i_sparse,
                                       unsigned int               i_numberOfRows,
                                       unsigned int               i_numberOfColumns,
                                 const std::vector<unsigned int> &i_rows,
                                 const std::vector<unsigned int> &i_columns,
                                 const std::vector<double>       &i_values,
                                       double*                    o_matrix );

    /**
     * Allocates memory for the global matrices and initializes it.
     *
     * @param i_matrixReader XML matrix reader.
     **/
    void initializeGlobalMatrices( const seissol::XmlParser &i_matrixReader );

  public:
    /**
     * Constructor, which allocates memory for the global matrices and initializes them.
     *
     * @param i_matrixReader XML matrix reader.
     **/
    MemoryManager( const seissol::XmlParser &i_matrixReader );

    /**
     * Destructor, which frees all allocated memory.
     **/
    ~MemoryManager();

    /**
     * Gets the pointers to the 52 memory chunks of the flux matrices.
     * Additionally the 53rd pointer gives the address of the stiffness matrix \f$ M^{-1} K^\xi \f$. 
     *
     * @return 52 pointers to the memory chunks of the flux matrices (and as 53rd pointer \f$ M^{-1} K^\xi \f$).
     **/
    double** getFluxMatrixPointers() const;

    /**
     * Get the pointers to the 2*3 (non-transposed and transposed) memory chunks of the stiffness matrices.
     * Additionally the seventh pointer gives the address of the flux matrix \f$ M^{-1} F^{-, 1} \f$.
     *
     * @return 6 pointers to the memory chunks of the flux matrices (and as 7th pointer \f$ M^{-1} F^{-, 1} \f$).
     **/    
    double** getStiffnessMatrixPointers() const;
};

#endif
