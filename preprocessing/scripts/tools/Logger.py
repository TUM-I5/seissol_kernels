#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
# Copyright (c) 2013, SeisSol Group
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
# Logs information of the scripts.
#
from __future__ import print_function

from datetime import datetime


identation = '  ';

# Log a given string.
#
# \param i_string string, which will be logged.
def log( i_string,
         i_level=1,
         i_returnString = False ):
  # resulting string
  l_completeString = ''

  # add time stamp
  l_identation = str(datetime.now());

  # indent to the given level
  for l_level in range(1, i_level):
    l_identation = l_identation + identation;

  l_lines = i_string.splitlines()

  # add the given string line by line
  for l_line in i_string.splitlines():
    l_completeString = l_completeString + l_identation  + '  ' + l_line + '\n'

  if( i_returnString == True ):
    return l_completeString
  else:
     print(l_completeString, end='')

# Prints the welome message.
#
def printWelcomeMessage():
  # define tliccense
  l_license = getTumBsdLicense()

  print('')
  log('*******************************************************')
  log('** Welcome to the offline assembly script of SeisSol **')
  log('*******************************************************')
  log('')
  log(getTumBsdLicense())
  log('*******************************************************')
  print('')

# Prints the finish message.
#
def printFinishMessage():
  print('')
  log('******************************************************************')
  log('** The offline assembly script of SeisSol finished successfully **')
  log('******************************************************************')
  print('')

# Get TUMs BSD license
def getTumBsdLicense():
  return '''Copyright (c) 2012-2014, SeisSol Group
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
'''

def getFileHeader():
  l_string = '''@file
This file is part of SeisSol.

@author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
@author Alexander Heinecke (alexander.heinecke AT mytum.de, http://www5.in.tum.de/wiki/index.php/Alexander_Heinecke,_M.Sc.,_M.Sc._with_honors)

'''
  l_string = l_string + '@date ' + str(datetime.now()) + '\n'

  l_string = l_string + '\n@section LICENSE\n' + getTumBsdLicense()

  l_string = l_string + '''\n@section DESCRIPTION
Remark: This file was generated.
''' 

  return l_string

# Writes an file header (@file, @author @license) to an file.
#
# @param i_pathToOutputFile path to the file, where the license will appended to.
# @param i_identation ident each line by this string.
def writeFileHeader( i_pathToOutputFile,
                     i_identation = '' ):

  with open( i_pathToOutputFile, 'a') as l_file:
    for l_line in getFileHeader().splitlines():
      l_file.write( i_identation  + l_line + '\n')
