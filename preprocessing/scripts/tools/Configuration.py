#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2015, SeisSol Group
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
#
# Configuration
#
import logging
import os.path

class Configuration():
  # ids of the matrices
  m_matrixIds = { "kXiDivM":    53,
                  "kEtaDivM":   54,
                  "kZetaDivM":  55,

                  "kXiDivMT":   56,
                  "kEtaDivMT":  57,
                  "kZetaDivMT": 58,

                  "fM1DivM":     0,
                  "fM2DivM":     1,
                  "fM3DivM":     2,
                  "fM4DivM":     3,

                  "fP111DivM":   4,
                  "fP112DivM":   5,
                  "fP113DivM":   6,
                  "fP121DivM":   7,
                  "fP122DivM":   8,
                  "fP123DivM":   9,
                  "fP131DivM":  10,
                  "fP132DivM":  11,
                  "fP133DivM":  12,
                  "fP141DivM":  13,
                  "fP142DivM":  14,
                  "fP143DivM":  15,

                  "fP211DivM":  16,
                  "fP212DivM":  17,
                  "fP213DivM":  18,
                  "fP221DivM":  19,
                  "fP222DivM":  20,
                  "fP223DivM":  21,
                  "fP231DivM":  22,
                  "fP232DivM":  23,
                  "fP233DivM":  24,
                  "fP241DivM":  25,
                  "fP242DivM":  26,
                  "fP243DivM":  27,

                  "fP311DivM":  28,
                  "fP312DivM":  29,
                  "fP313DivM":  30,
                  "fP321DivM":  31,
                  "fP322DivM":  32,
                  "fP323DivM":  33,
                  "fP331DivM":  34,
                  "fP332DivM":  35,
                  "fP333DivM":  36,
                  "fP341DivM":  37,
                  "fP342DivM":  38,
                  "fP343DivM":  39,

                  "fP411DivM":  40,
                  "fP412DivM":  41,
                  "fP413DivM":  42,
                  "fP421DivM":  43,
                  "fP422DivM":  44,
                  "fP423DivM":  45,
                  "fP431DivM":  46,
                  "fP432DivM":  47,
                  "fP433DivM":  48,
                  "fP441DivM":  49,
                  "fP442DivM":  50,
                  "fP443DivM":  51,

                  "fluxSolver": 52,
                  "starMatrix": 59  }

  # names of matrix ids
  m_matrixNames = {}

  # matrix bindings in the corresponding integration kernels
  m_matrixBinds = {
                    'time': {            "kXiDivMT":    0,
                                         "kEtaDivMT":   1,
                                         "kZetaDivMT":  2,
                                         "starMatrix":  3  },
                    'volume': {
                                         "kXiDivM":     0,
                                         "kEtaDivM":    1,
                                         "kZetaDivM":   2,
                                         "starMatrix":  3  },

                    'boundary': {
                                         "fM1DivM":     0,
                                         "fM2DivM":     1,
                                         "fM3DivM":     2,
                                         "fM4DivM":     3,

                                         "fP111DivM":   4,
                                         "fP112DivM":   5,
                                         "fP113DivM":   6,
                                         "fP121DivM":   7,
                                         "fP122DivM":   8,
                                         "fP123DivM":   9,
                                         "fP131DivM":  10,
                                         "fP132DivM":  11,
                                         "fP133DivM":  12,
                                         "fP141DivM":  13,
                                         "fP142DivM":  14,
                                         "fP143DivM":  15,

                                         "fP211DivM":  16,
                                         "fP212DivM":  17,
                                         "fP213DivM":  18,
                                         "fP221DivM":  19,
                                         "fP222DivM":  20,
                                         "fP223DivM":  21,
                                         "fP231DivM":  22,
                                         "fP232DivM":  23,
                                         "fP233DivM":  24,
                                         "fP241DivM":  25,
                                         "fP242DivM":  26,
                                         "fP243DivM":  27,

                                         "fP311DivM":  28,
                                         "fP312DivM":  29,
                                         "fP313DivM":  30,
                                         "fP321DivM":  31,
                                         "fP322DivM":  32,
                                         "fP323DivM":  33,
                                         "fP331DivM":  34,
                                         "fP332DivM":  35,
                                         "fP333DivM":  36,
                                         "fP341DivM":  37,
                                         "fP342DivM":  38,
                                         "fP343DivM":  39,

                                         "fP411DivM":  40,
                                         "fP412DivM":  41,
                                         "fP413DivM":  42,
                                         "fP421DivM":  43,
                                         "fP422DivM":  44,
                                         "fP423DivM":  45,
                                         "fP431DivM":  46,
                                         "fP432DivM":  47,
                                         "fP433DivM":  48,
                                         "fP441DivM":  49,
                                         "fP442DivM":  50,
                                         "fP443DivM":  51,

                                         "fluxSolver": 52 }
                  }

  m_globalMatrices = { 'time': 3, 'volume': 3, 'boundary': 52 }

  m_matricesDir = ""
  m_matrixMarketFiles = {}

  m_architectures = ['wsm', 'snb', 'hsw', 'knc', 'noarch']


  m_bytesPerReal = { 's': 4,
                     'd': 8 }

  m_alignments = {  'wsm': 16,
                    'snb': 32,
                    'hsw': 32,
                    'knc': 64,
                    'noarch': 16 }

  m_maximumDegree = 0

  m_pathToSparseDenseConfigs = '-1'
  m_pathToGemmCodeGenerator = '-1'
  m_pathToGeneratedCodeDir = '-1'

  # Constuctor
  def __init__( self,
                i_matricesDir              = "matrices",
                i_maximumOrder             = 8,
                i_pathToSparseDenseConfigs = 'sparse_dense',
                i_pathToGemmCodeGenerator  = 'SeisSolGen/generator.exe',
                i_pathToGeneratedCodeDir   = 'generated_code'):
    logging.debug( "Constructed a Configuration()" )

    # set up directories
    self.m_pathToSparseDenseConfigs = i_pathToSparseDenseConfigs
    self.m_pathToGemmCodeGenerator  = i_pathToGemmCodeGenerator
    self.m_pathToGeneratedCodeDir   = i_pathToGeneratedCodeDir

    # set maximum degree
    self.m_maximumDegree = i_maximumOrder - 1

    # inverse name -> key dict
    self.m_matrixNames = {v: k for k, v in self.m_matrixIds.items()}

    # check that the matrices dir exists
    if( not os.path.exists( i_matricesDir ) ):
      logging.error( "matrices directory does not exist: " + i_matricesDir )
      exit(1)

    # setup paths to global matrix market files
    for l_order in range(2, i_maximumOrder+1):
      # create a new dict for this order
      self.m_matrixMarketFiles[l_order] = {}

      for l_matrix in self.m_matrixIds.keys():
        if( l_matrix is "fluxSolver" or l_matrix is "starMatrix" ):
          continue

        self.m_matrixMarketFiles[l_order][l_matrix] = i_matricesDir + "/" + l_matrix + "_3D_" + str(l_order-1) + "_maple.mtx"

        if( not os.path.exists( self.m_matrixMarketFiles[l_order][l_matrix] ) ):
          logging.error( "matrix market files does not exist: " + self.m_matrixMarketFiles[l_order][l_matrix] )
          exit(1)

    # setup path to local matrix market files
    self.m_matrixMarketFiles["starMatrix"] = i_matricesDir + "/starMatrix_3D_maple.mtx"
    if( not os.path.exists( self.m_matrixMarketFiles["starMatrix"] ) ):
      logging.error( "matrix market files does not exist: " + self.m_matrixMarketFiles["starMatrix"] )
      exit(1)
