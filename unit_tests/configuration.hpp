/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Configuration for the kernel generation unit tests.
 **/

#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include <cstdio>
#include <Initializer/typedefs.hpp>
#include <cxxtest/TestSuite.h>
#include <cxxtest/GlobalFixture.h>

//! dummy variable for the FLOP counter
static unsigned long long libxsmm_num_total_flops = 0;

namespace unit_test {
  class Configuration;
}

/**
 * Holds the kernel generation unit tests.
 */
class unit_test::Configuration: public CxxTest::TestSuite {
  private:
    //! how often should a test with random numebers be repeated?
    int m_numberOfRepeats;

    //! path to the source directory of the unit tests
    std::string m_unitTestsSrcDirectory;

    //! path to the matrices directory
    std::string m_matricesDirectory;

    //! minimum order to test
    unsigned int m_minimumOrder;
    unsigned int m_maximumOrder;

    //! zero tolerance
    real m_zeroTolerance;

    //! sparse dense switch
    int m_sparseSwitch[60];

  public:
    /**
     * Constructor, which sets up the configuration.
     */
    Configuration(){
      // setup the configuration
      setup();
    }

    /**
     * Set up the configuration.
     **/
    void setup(){
      m_numberOfRepeats = 500;
    
      // get current directory
      m_unitTestsSrcDirectory = __FILE__;
      m_unitTestsSrcDirectory.erase(  m_unitTestsSrcDirectory.end() - 18,  m_unitTestsSrcDirectory.end() );

      // derive matrices directory
      m_matricesDirectory = m_unitTestsSrcDirectory  + "/../preprocessing/matrices";

      m_minimumOrder = 2;
      m_maximumOrder = 6;

      m_zeroTolerance = 10e-13;

      // set up sparse switch
#define SPARSE_SWITCH
#include <initialization/bind.h>
#undef SPARSE_SWITCH
    }

    /**
     * Get the number of repeats.
     * @return number of repeats.
     **/
    int getNumberOfRepeats() const {
      return m_numberOfRepeats;
    }

   /**
    * Get the path of the unit test source directory.
    * @return path to the unit test src directory.
    **/
    std::string getUnitTestsSrcDirectory() const {
      return m_unitTestsSrcDirectory;
    }

   /**
    * Get the path of the matrices directory.
    * @return path to the matrices directory.
    **/
    std::string getMatricesDirectory() const {
      return m_matricesDirectory;
    }

   /**
    * Get the path of the matrices files.
    * @return path to the matrices file.
    **/
    std::string getMatricesFile() const {
      std::stringstream l_matricesFile;
      l_matricesFile << m_matricesDirectory << "/matrices_" << NUMBER_OF_BASIS_FUNCTIONS << ".xml";
      return l_matricesFile.str();
    }

    /**
     * Get the minimum order to test against.
     * @return minimum order.
     **/
    unsigned int getMinimumOrder() const {
      return m_minimumOrder;
    }

    /**
     * Get the maximum order to test against.
     * @return maximum order.
     **/
    unsigned int getMaximumOrder() const {
      return m_maximumOrder;
    }

    /**
     * Get the numerical zero toleratnce to test again.
     * @return numerical zero tolerance.
     **/
    real getNumericalZeroTolerance() const {
      return m_zeroTolerance;
    }

    /**
     * Get the sparse switch of a specific matrix.
     * @param i_matrixId matrix id.
     * @return true if sparse, false if dense.
     **/
    bool getSparseSwitch( unsigned int i_matrixId ) {
      if( m_sparseSwitch[i_matrixId] != -1 ) return true;
      else return false;
    }
};

#endif
