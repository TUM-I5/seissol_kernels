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

# seissol_proxy regular expressions
l_proxyExp = { 'float'    : '\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?',
               'time'      : 'time for seissol proxy  :',
}

# command line interface
logging.info( "parsing command line arguments" )

l_parser    = argparse.ArgumentParser( description='Performs a convergence analysis.' )
l_parser.add_argument( '--log_dir',
                       dest     = "log_dir",
                       required = True,
                       help     = "Directory where the log files are located.",
                       metavar  = "LOG_DIR" )

l_parser.add_argument( '--output_dir',
                       dest     = "output_dir",
                       required = True,
                       help     = "Path to the output directory where the PDFs go.",
                       metavar  = "OUTPUT_DIR" )

l_arguments = vars(l_parser.parse_args())

# get all files matching the regexp
l_files = glob.glob( l_arguments['log_dir']+'/*tune*.log' )

# cut off directories
for l_file in range(len(l_files)):
  l_files[l_file] = os.path.basename(l_files[l_file])

l_runtime = {}
for l_file in l_files: 
  l_match = re.search( "tune_local_([0-9])+_([0-9])+_([0-9])+_([0-9])+_([0-9])+.log" , l_file )
  # create a new dictionary for this file
  l_runtime[l_file] = {}

  # save the information of the file name
  l_runtime[l_file]['star_matrix']           = l_match.group(1)
  l_runtime[l_file]['time_kernel']           = l_match.group(2)
  l_runtime[l_file]['volume_kernel']         = l_match.group(3)
  l_runtime[l_file]['boundary_kernel_local'] = l_match.group(4)
  l_runtime[l_file]['order']                 = l_match.group(5)

  # create an empty list for the runtime
  l_measurements = []

  # get the runtime
  l_fileContents = open( l_arguments['log_dir']+'/'+l_file, "r" )

  # iterate over log file
  for l_line in l_fileContents:
    l_match = re.search('(?<='+l_proxyExp['time']+')'+l_proxyExp['float'], l_line )

    if l_match:
      l_measurements = l_measurements + [ float(l_match.group(0)) ]

  # compute statistic values
  l_runtime[l_file]['repeats'] = len( l_measurements )
  if len( l_measurements ) > 0:
    l_runtime[l_file]['min']     = min( l_measurements )
    l_runtime[l_file]['mean']    = statistics.mean( l_measurements )
    l_runtime[l_file]['max']     = max( l_measurements )
    l_runtime[l_file]['stddev']  = statistics.stdev( l_measurements )
    l_runtime[l_file]['var']     = statistics.variance( l_measurements )

# save to csv file
with open( l_arguments['output_dir']+'timings.csv', 'w' ) as l_csvFile:
  l_fieldNames = ['order', 'star_matrix', 'time_kernel', 'volume_kernel', 'boundary_kernel_local']
  l_fieldNames = l_fieldNames + ['repeats', 'min', 'mean', 'max', 'stddev', 'var']
  l_writer = csv.DictWriter(l_csvFile, fieldnames=l_fieldNames)
  l_writer.writeheader()

  for l_result in l_runtime.keys():
    l_writer.writerow( l_runtime[l_result] )
