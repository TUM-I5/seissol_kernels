#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
# Analyzes the results of the tuning runs.

import logging
import argparse
import glob
import os
import re
import csv
import statistics
import dicttoxml
from xml.dom.minidom import parseString

# seissol_proxy regular expressions
l_proxyExp = { 'float'    : '\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?',
               'time'      : 'time for seissol proxy  :',
}

# command line interface
logging.info( "parsing command line arguments" )

l_parser    = argparse.ArgumentParser( description='Performs a convergence analysis.' )
l_parser.add_argument( '--log_dirs',
                       dest     = "log_dirs",
                       required = True,
                       nargs    = "+",
                       help     = "Directories where the log files are located.",
                       metavar  = "LOG_DIR" )

l_parser.add_argument( '--output_dir',
                       dest     = "output_dir",
                       required = True,
                       help     = "Path to the output directory where the PDFs go.",
                       metavar  = "OUTPUT_DIR" )

l_parser.add_argument( '--raw',
                       action  = "store_true",
                       help    = "Generate raw csv-files: One line per run." )

l_parser.add_argument( '--type',
                        choices=['all', 'local', 'neighboring'],
                        default='all' )

l_arguments = vars(l_parser.parse_args())

if( not l_arguments['raw'] and len( l_arguments['log_dirs'] ) ):
  print 'non-raw output not supported for multiple directories'
  exit()

# get all files matching the regexp
l_files = []
for l_dir in l_arguments['log_dirs']:
  l_files = l_files + glob.glob( l_dir+'/*tune*.log' )

l_runtime = {}

for l_file in l_files:
  # cut off directories
  l_base = os.path.basename(l_file)

  # local integration
  if 'local' in l_base and l_arguments['type'] in ['all', 'local']:
    l_match = re.search( "tune_local_([0-9])+_([0-9])+_([0-9])+_([0-9])+_([0-9])+.log" , l_base )

    # create a new dictionary for this base if required
    if( not l_base in l_runtime ):
      l_runtime[l_base] = {}

      # save the information of the file name
      l_runtime[l_base]['star_matrix']                 = int(l_match.group(1))
      l_runtime[l_base]['time_kernel']                 = int(l_match.group(2))
      l_runtime[l_base]['volume_kernel']               = int(l_match.group(3))
      l_runtime[l_base]['boundary_kernel_local']       = int(l_match.group(4))
      l_runtime[l_base]['boundary_kernel_neighboring'] = 0
      l_runtime[l_base]['order']                       = int(l_match.group(5))
      l_runtime[l_base]['time']                        = ''
  elif 'neighboring' in l_base and l_arguments['type'] in ['all', 'neighboring']:
    l_match = re.search( "tune_neighboring_0_([0-9])+_([0-9])+.log" , l_base )

    # create a new dictionary for this base if required
    if( not l_base in l_runtime ):
      l_runtime[l_base] = {}

      # save the information of the file name
      l_runtime[l_base]['star_matrix']                 = 0
      l_runtime[l_base]['time_kernel']                 = 0
      l_runtime[l_base]['volume_kernel']               = 0
      l_runtime[l_base]['boundary_kernel_local']       = 0
      l_runtime[l_base]['boundary_kernel_neighboring'] = int(l_match.group(1))
      l_runtime[l_base]['order']                       = int(l_match.group(2))
      l_runtime[l_base]['time']                        = ''
  else:
    continue

  # create an empty list for the runtime
  l_measurements = []

  # get the runtime
  l_fileContents = open( l_file, "r" )

  # iterate over log file
  for l_line in l_fileContents:
    l_match = re.search('(?<='+l_proxyExp['time']+')'+l_proxyExp['float'], l_line )

    if l_match:
      l_measurements = l_measurements + [ float(l_match.group(0)) ]

  # compute statistic values
  l_runtime[l_base]['repeats'] = len( l_measurements )
  if len( l_measurements ) > 0:
    if l_arguments['raw']:
      l_runtime[l_base]['time']    = l_runtime[l_base]['time']  + str(l_measurements)
      l_runtime[l_base]['time']    = l_runtime[l_base]['time'].replace("][",",")
    else:
      l_runtime[l_base]['min']     = min( l_measurements )
      l_runtime[l_base]['mean']    = statistics.mean( l_measurements )
      l_runtime[l_base]['max']     = max( l_measurements )
      l_runtime[l_base]['stddev']  = statistics.stdev( l_measurements )
      l_runtime[l_base]['var']     = statistics.variance( l_measurements )

if not l_arguments['raw']:
  # save to csv file
  with open( l_arguments['output_dir']+'/timings.csv', 'w' ) as l_csvFile:
    l_fieldNames = ['order', 'star_matrix', 'time_kernel', 'volume_kernel', 'boundary_kernel_local', 'boundary_kernel_neighboring']
    l_fieldNames = l_fieldNames + ['repeats', 'min', 'mean', 'max', 'stddev', 'var']
    l_writer = csv.DictWriter(l_csvFile, fieldnames=l_fieldNames)
    l_writer.writeheader()

    for l_result in l_runtime.keys():
      l_writer.writerow( l_runtime[l_result] )

else:
  # save raw results to xml
  with open( l_arguments['output_dir']+'/timings_'+l_arguments['type']+'.xml', 'w' ) as l_csvFile:
    l_csvFile.write( parseString( dicttoxml.dicttoxml(l_runtime) ).toprettyxml() )
