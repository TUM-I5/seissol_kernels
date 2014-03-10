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
 * Configuration for the kernel generation unit tests.
 **/

#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include <cstdio>
#include <cxxtest/TestSuite.h>
#include <cxxtest/GlobalFixture.h>

//! dummy variable for the FLOP counter
static unsigned long long num_flops = 0;

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
    double m_zeroTolerance;

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
    double getNumericalZeroTolerance() const {
      return m_zeroTolerance;
    }
};

#endif
