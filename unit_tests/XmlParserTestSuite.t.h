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
 * Test suite, which tests the XML parser.
 **/

#include "configuration.hpp"
#include "Initializer/XmlParser.hpp"

namespace unit_test{
  class XmlParserTestSuite;
}

class unit_test::XmlParserTestSuite: public CxxTest::TestSuite {
  //private
  public:
    //! Configuration of the unit tests
    unit_test::Configuration m_configuration;

    /**
     * Tests the xml parser against a hardcoded XML-file.
     **/
    void testXmlParser() {
      // set up test values
      const unsigned int l_testIds[]             = { 18, 2, 4, 6, 8, 13, 14, 5, 10, 16,
                                                     19, 11, 12, 9, 3, 15, 7, 1, 17, 20    };

      const std::string l_testNames[]            = { "223", "112", "121", "123", "132",
                                                     "211", "212", "122", "141", "221",
                                                     "232", "142", "143", "133", "113",
                                                     "213", "131", "111", "222", "233"     };

      const unsigned int l_testNumberOfRows[]    = { 20, 4, 10, 84, 56, 84, 56, 20, 20, 10,
                                                     84, 20, 20, 10, 4, 10, 20, 84, 56, 20 };

      const unsigned int l_testNumberOfColumns[] = { 20, 4, 10, 84, 56, 84, 56, 20, 20, 10,
                                                     84, 20, 20, 10, 4, 10, 20, 5, 2, 1    };

      const bool l_testSparsities[]              = { true, true,  false, true,  false,
                                                     true, true,  true,  true,  true,
                                                     true, false, true,  true,  true,
                                                     true, true,  true,  false, true       };
      
      std::vector< std::vector<unsigned int> > l_testRows;
      std::vector< std::vector<unsigned int> > l_testColumns;
      std::vector< std::vector<real>       > l_testValues;

      // add first speficied matrix
      std::vector<unsigned int> l_rows;
      std::vector<unsigned int> l_columns;
      std::vector<real        > l_values;
      l_rows.push_back(     1            );
      l_columns.push_back( 10            ); 
      l_values.push_back( 421.2215141241 );
      
      l_rows.push_back(     3             );
      l_columns.push_back( 10             );
      l_values.push_back(  42.29146891241 );

      l_rows.push_back(     6              );
      l_columns.push_back( 11              );
      l_values.push_back(   5.221445141241 );

      l_rows.push_back(     5              );
      l_columns.push_back( 12              );
      l_values.push_back(   3.295146891241 );

      l_testRows.push_back(    l_rows    );
      l_testColumns.push_back( l_columns );
      l_testValues.push_back(  l_values  );
      l_rows.clear();
      l_columns.clear();
      l_values.clear();

      // add first 11 matrices without entries
      for( int l_i = 0; l_i < 11; l_i++) {
        l_testRows.push_back(    std::vector<unsigned int>() );
        l_testColumns.push_back( std::vector<unsigned int>() );
        l_testValues.push_back(  std::vector<real>()         );
      }

      // add second specified matrix
      l_rows.push_back(     5             );
      l_columns.push_back(  3             ); 
      l_values.push_back(  91.42215141241 );
      
      l_rows.push_back(     5             );
      l_columns.push_back(  4             );
      l_values.push_back(  21.29146891241 );

      l_rows.push_back(    13            );
      l_columns.push_back( 12            );
      l_values.push_back(  86.21445141241);

      l_testRows.push_back(    l_rows    );
      l_testColumns.push_back( l_columns );
      l_testValues.push_back(  l_values  );

      // add remaining 7 matrices without entries
      for( int l_i = 0; l_i < 7; l_i++) {
        l_testRows.push_back(    std::vector<unsigned int>() );
        l_testColumns.push_back( std::vector<unsigned int>() );
        l_testValues.push_back(  std::vector<real>()         );
      }

      // add test xml file
      std::string l_matricesPath = m_configuration.getUnitTestsSrcDirectory() + "/globalTestMatrices.xml";

      // setup the xml-parser
      seissol::XmlParser l_matrixReader( l_matricesPath );

      //! vectors, which hold information about our matrices
      std::vector< unsigned int > l_matrixIds;
      std::vector< std::string  > l_matrixNames;
      std::vector< bool         > l_matrixSparsities;
      std::vector< unsigned int > l_matrixNumberOfRows;
      std::vector< unsigned int > l_matrixNumberOfColumns;

      std::vector< std::vector<unsigned int> > l_matrixRows;
      std::vector< std::vector<unsigned int> > l_matrixColumns;
      std::vector< std::vector<real>       >   l_matrixValues;
      
      // read the flux matrices
      l_matrixReader.readGlobalMatrices( "flux",
                                         l_matrixIds, l_matrixNames,
                                         l_matrixNumberOfRows, l_matrixNumberOfColumns, l_matrixSparsities,
                                         l_matrixRows, l_matrixColumns, l_matrixValues );

      // assert we have same dimensions
      TS_ASSERT_EQUALS( l_matrixIds.size(), l_matrixNames.size()      );
      TS_ASSERT_EQUALS( l_matrixIds.size(), l_matrixSparsities.size() );
      TS_ASSERT_EQUALS( l_matrixIds.size(), l_matrixRows.size()       );
      TS_ASSERT_EQUALS( l_matrixIds.size(), l_matrixColumns.size()    );
      TS_ASSERT_EQUALS( l_matrixIds.size(), l_matrixValues.size()     );
      TS_ASSERT_EQUALS( l_matrixIds.size(), 20                        );

      // compare values of the parser against test values
      for( int l_matrix = 0; l_matrix < l_matrixIds.size(); l_matrix++) {
        // assert meta information matches
        TS_ASSERT_EQUALS( l_testIds            [l_matrix], l_matrixIds            [l_matrix] );
        TS_ASSERT_EQUALS( l_testNames          [l_matrix], l_matrixNames          [l_matrix] );
        TS_ASSERT_EQUALS( l_testNumberOfRows   [l_matrix], l_matrixNumberOfRows   [l_matrix] );
        TS_ASSERT_EQUALS( l_testNumberOfColumns[l_matrix], l_matrixNumberOfColumns[l_matrix] );
        TS_ASSERT_EQUALS( l_testSparsities     [l_matrix], l_matrixSparsities     [l_matrix] );

        // assert number of entries matches
        TS_ASSERT_EQUALS( l_matrixRows[l_matrix].size(), l_matrixColumns[l_matrix].size() );
        TS_ASSERT_EQUALS( l_matrixRows[l_matrix].size(), l_matrixValues[l_matrix].size()  );
        TS_ASSERT_EQUALS( l_matrixRows[l_matrix].size(), l_testRows[l_matrix].size()      );

        for( int l_entry = 0; l_entry < l_matrixRows[l_matrix].size(); l_entry++) {
           TS_ASSERT_DELTA( l_matrixValues[l_matrix][l_entry],
                            l_testValues[l_matrix][l_entry],
                            m_configuration.getNumericalZeroTolerance() );
        }
      }
    }
};
