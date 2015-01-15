#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2014, SeisSol Group
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
# Sets up the sparse-dense switches.
#
import argparse
import logging
import pprint
import numpy
import scipy.io.mmio
import scipy.sparse
import operator
import xml.etree.ElementTree as etree
import sys
import os

# matrices
l_matrices = {
  'time':                  [ "kXiDivMT",
                             "kEtaDivMT",
                             "kZetaDivMT"
                           ],

  'volume':                [ "kXiDivM",
                             "kEtaDivM",
                             "kZetaDivM"
                           ],

  'boundary_local':        [ "fM1DivM",
                             "fM2DivM",
                             "fM3DivM",
                             "fM4DivM"
                           ],

  'boundary_neighboring': [ "fP111DivM",
                            "fP112DivM",
                            "fP113DivM",
                            "fP121DivM",
                            "fP122DivM",
                            "fP123DivM",
                            "fP131DivM",
                            "fP132DivM",
                            "fP133DivM",
                            "fP141DivM",
                            "fP142DivM",
                            "fP143DivM",

                            "fP211DivM",
                            "fP212DivM",
                            "fP213DivM",
                            "fP221DivM",
                            "fP222DivM",
                            "fP223DivM",
                            "fP231DivM",
                            "fP232DivM",
                            "fP233DivM",
                            "fP241DivM",
                            "fP242DivM",
                            "fP243DivM",

                            "fP311DivM",
                            "fP312DivM",
                            "fP313DivM",
                            "fP321DivM",
                            "fP322DivM",
                            "fP323DivM",
                            "fP331DivM",
                            "fP332DivM",
                            "fP333DivM",
                            "fP341DivM",
                            "fP342DivM",
                            "fP343DivM",

                            "fP411DivM",
                            "fP412DivM",
                            "fP413DivM",
                            "fP421DivM",
                            "fP422DivM",
                            "fP423DivM",
                            "fP431DivM",
                            "fP432DivM",
                            "fP433DivM",
                            "fP441DivM",
                            "fP442DivM",
                            "fP443DivM" ]
}

l_genKernels = [ 'swsm','dwsm', 'ssnb', 'dsnb', 'sknc', 'dknc', 'shsw', 'dhsw', 'snoarch', 'dnoarch' ]

# add a command line parser
l_parser    = argparse.ArgumentParser( description='Sets up the sparse-dense switches.' )
l_parser.add_argument( '--matrices_dir',
                       dest     = "matrices_dir",
                       required = True,
                       help     = "Path to matrices directory.",
                       metavar  = "MATRICES_DIR" )

l_parser.add_argument( '--output_dir',
                       dest     = "output_dir",
                       required = True,
                       help     = "Path to output directory.",
                       metavar  = "OUTPUT_DIR" )

l_arguments = vars(l_parser.parse_args())


# create a logger
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
l_logger = logging.getLogger('simple_example')
l_logger.setLevel(logging.INFO)

l_logger.info('**************************************************')
l_logger.info('| Welcome to the sparse-dense generation script! |')
l_logger.info('**************************************************')

# split flux matrices into blocks for nnz
l_nnzFlux = {}
for l_fluxMatrix in l_matrices['boundary_neighboring']:
  l_pathToMatrix = l_arguments['matrices_dir'] + '/' + l_fluxMatrix + '_3D_7_maple.mtx'

  # read the matrix
  l_denseMatrix = scipy.io.mmio.mmread( l_pathToMatrix )

  # convert to a list of lists
  l_denseMatrix = scipy.sparse.lil_matrix(l_denseMatrix)

  l_nnz = len(numpy.nonzero(l_denseMatrix)[0])
  l_nnzFlux[l_fluxMatrix] = l_nnz

# sort the dictionary by non-zeros
l_nnzFlux = sorted(l_nnzFlux.items(), key=operator.itemgetter(1))

# only work on names
l_nnzFlux = [e[0] for e in l_nnzFlux]

# create the local auto-tuning runs
for l_starSparse in [ [], ['starMatrix'] ]:
  for l_boundaryNeighboring in range( -1, len(l_nnzFlux) ):
    l_boundaryNeighboringSparse = l_nnzFlux[:l_boundaryNeighboring+1]

    # setup xml
    l_root = etree.Element("sparse_matrices")

    for l_order in range(2, 9):
      l_order = etree.SubElement( l_root, "O"+str(l_order) )

      # write dummies for time and volume
      etree.SubElement( l_order, 'time' )
      etree.SubElement( l_order, 'volume' )

      # write star configuration
      l_star = etree.SubElement( l_order, 'local' )
      for l_matrix in l_starSparse:
        etree.SubElement( l_star, l_matrix )

      # write boundary configuration
      l_boundary = etree.SubElement( l_order, 'boundary' )
      for l_matrix in l_boundaryNeighboringSparse:
        etree.SubElement( l_boundary, l_matrix )

      # warp into a tree and write
      l_tree = etree.ElementTree( l_root )
      for l_kernelConfig in l_genKernels:
        l_outDir = l_arguments['output_dir'] + '/neighboring' + '_' + str(len(l_starSparse)) \
                                                              + '_' + str(len(l_boundaryNeighboringSparse))
        if not os.path.exists(l_outDir):
          os.makedirs(l_outDir)
        l_tree.write( l_outDir + '/' + l_kernelConfig + '.xml' )

# create the local auto-tuning runs
for l_starSparse in [ [], ['starMatrix'] ]:
  for l_stiffTime in range( -1, len(l_matrices['time']) ):
    l_stiffTimeSparse = l_matrices['time'][:l_stiffTime+1]

    for l_stiffVolume in range( -1, len(l_matrices['volume']) ):
      l_stiffVolumeSparse = l_matrices['volume'][:l_stiffVolume+1]

      for l_boundaryLocal in range( -1, len(l_matrices['boundary_local']) ):
        l_boundaryLocalSparse = l_matrices['boundary_local'][:l_boundaryLocal+1]

        # setup xml
        l_root = etree.Element("sparse_matrices")

        for l_order in range(2, 9):
          l_order = etree.SubElement( l_root, "O"+str(l_order) )

          # write star configuration
          l_star = etree.SubElement( l_order, 'local' )
          for l_matrix in l_starSparse:
            etree.SubElement( l_star, l_matrix )

          # write time configuration
          l_time = etree.SubElement( l_order, 'time' )
          for l_matrix in l_stiffTimeSparse:
            etree.SubElement( l_time, l_matrix )

          # write volume configuration
          l_volume = etree.SubElement( l_order, 'volume' )
          for l_matrix in l_stiffVolumeSparse:
            etree.SubElement( l_volume, l_matrix )

          # write boundary configuration
          l_boundary = etree.SubElement( l_order, 'boundary' )
          for l_matrix in l_boundaryLocalSparse:
            etree.SubElement( l_boundary, l_matrix )


        # warp into a tree and write
        l_tree = etree.ElementTree( l_root )
        for l_kernelConfig in l_genKernels:
          l_outDir = l_arguments['output_dir'] + '/local'                               +\
                                                  '_'  + str(len(l_starSparse))         +\
                                                  '_'  + str(len(l_stiffTimeSparse))    +\
                                                  '_'  + str(len(l_stiffVolumeSparse))  +\
                                                  '_'  + str(len(l_boundaryLocalSparse))
          if not os.path.exists(l_outDir):
            os.makedirs(l_outDir)                    
          l_tree.write( l_outDir + '/' + l_kernelConfig + '.xml' )

l_logger.info('********************************************')
l_logger.info('| sparse-dense generation script finished! |')
l_logger.info('********************************************')
