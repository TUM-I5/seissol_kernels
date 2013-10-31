#!/bin/env python
# @file
# This file is part of SeisSol.
# 
# @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
# 
# 
# @section LICENSE
# This software was developed at Technische Universitaet Muenchen, who is the owner of the software.
# 
# According to good scientific practice, publications on results achieved in whole or in part due to this software should cite at least one paper or referring to an URL presenting this software.
# 
# The owner wishes to make the software available to all users to use, reproduce, modify, distribute and redistribute also for commercial purposes under the following conditions of the original BSD license. Linking this software module statically or dynamically with other modules is making a combined work based on this software. Thus, the terms and conditions of this license cover the whole combination. As a special exception, the copyright holders of this software give you permission to link it with independent modules or to instantiate templates and macros from this software's source files to produce an executable, regardless of the license terms of these independent modules, and to copy and distribute the resulting executable under terms of your choice, provided that you also meet, for each linked independent module, the terms and conditions of this license of that module.
# 
# Copyright (c) 2013
# Technische Universitaet Muenchen
# Department of Informatics
# Chair of Scientific Computing
# http://www5.in.tum.de/
# 
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# All advertising materials mentioning features or use of this software must display the following acknowledgement: This product includes software developed by the Technische Universitaet Muenchen (TUM), Germany, and its contributors.
# Neither the name of the Technische Universitaet Muenchen, Munich, Germany nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
#
# @section DESCRIPTION
#
# Generates batch scripts.
#

import Logger as l_logger
import numpy
import re

class BatchGenerator():
  #! setup
  m_setup = {}

  #! holds a batch script
  m_batchScript = {}
  
  # Constructor
  #
  def __init__( self ):
    self.m_setup['name']                        = 'performance_generated_kernels'
    self.m_setup['tmp_dir_path']                = ''
#    self.m_setup['build_environment']           = 'slurm'
    self.m_setup['#ranks']                      = 16
    self.m_setup['#cores per rank']             = 1
#    self.m_setup['cluster']                     = 'mac_cluster'
    self.m_setup['build_environment']           = 'loadleveler'
    self.m_setup['cluster']                     = 'supermuc'
    self.m_setup['seissol_dir']                 = '${HOME}/workspace/seissol_3d/'
    self.m_setup['benchmark_dir']               = '${HOME}/workspace/seissol_benchmarks/performance/Tetra_Elastic'
    self.m_setup['mail_address']                = 'breuera@in.tum.de'
    self.m_setup['preparatory_command']         = ''
    
    self.m_setup['run_setup']                   = [ { 'binary':     '${HOME}/workspace/seissol_3d/build/SeisSol_release_generatedKernels_mpi_9_56',
                                                      'config':     'PARAMETERS.par',
                                                      'output_dir': 'generated',
                                                    },
                                                    { 'binary':     '${HOME}/workspace/seissol_3d/build/SeisSol_release_classic_mpi_9_56',
                                                      'config':     'PARAMETERS.par',
                                                      'output_dir': 'classic'          
                                                    }
                                                  ]

    self.setup()

  def setup( self ):
    l_logger.log('setting up')

    # set up 
    self.m_setup['tmp_dir_name']                = self.m_setup['tmp_dir_path'] + '/' + self.m_setup['name'] + '_tmp'
    self.m_setup['file_name']                   = self.m_setup['name'] + '.' + self.m_setup['build_environment']

    # print the setup
    for l_key in self.m_setup:
      l_logger.log( str(l_key)+': '+str(self.m_setup[l_key]), 2 )
  
  # Generates the batch scripts
  #
  # @return file name of the written script
  def generateBatchScript( self ):
    l_logger.log('generating batch script ' + self.m_setup['file_name'], 1)

    # generate slurm header
    if( self.m_setup['build_environment'] == 'slurm' and self.m_setup['cluster'] == 'mac_cluster' ):
      self.m_batchScript['header'] =                                    \
        '#!/bin/bash \n'                                               +\
        '#SBATCH -o ' + self.m_setup['name'] + '.%j.%N.out \n'         +\
        '#SBATCH -J ' + self.m_setup['name'] +'\n'                     +\
        '#SBATCH --get-user-env \n'

      if( self.m_setup['cluster'] == 'mac_cluster' ):
        # use sandy bridge partition on mac cluster and disable TURBO
        self.m_batchScript['header'] = self.m_batchScript['header'] + '#SBATCH --partition=snb \n'        +\
                                                                      '#SBATCH --constraint=turbo_off \n'
      else:
        assert(false)
      self.m_batchScript['header'] = self.m_batchScript['header']                +\
        '#SBATCH --ntasks=' + str(self.m_setup['#ranks'])                 + '\n' +\
        '#SBATCH --cpus-per-task=' + str(self.m_setup['#cores per rank']*2) + '\n' +\
        '#SBATCH --mail-type=end                                             \n' +\
        '#SBATCH --mail-user=' + self.m_setup['mail_address']             + '\n' +\
        '#SBATCH --export=NONE                                               \n' +\
        '#SBATCH --time=24:00:00                                             \n' +\
        '#SBATCH --nice=500'
    elif(  self.m_setup['build_environment'] == 'loadleveler' and self.m_setup['cluster'] == 'supermuc' ):
      self.m_batchScript['header'] = \
'''#!/bin/bash
## we are optimizing not saving power
#@ energy_policy_tag = NONE
#@ job_type = parallel
#@ class = test
#@ island_count = 1,1
#@ node = 1
#@ tasks_per_node=16
#@ wall_clock_limit = 00:10:00
##                    0 h 10 min 00 secs
'''+\
      '#@ job_name = ' + self.m_setup['name'] +             '\n' +\
      '#@ network.MPI = sn_all,not_shared,us                 \n' +\
      '#@ initialdir = $(home)                               \n' +\
      '#@ output = ' + self.m_setup['name'] + '_$(jobid).out \n' +\
      '#@ error = ' + self.m_setup['name'] + '_$(jobid).err  \n' +\
      '#@ notification=always                                \n' +\
      '#@ notify_user=breuera@in.tum.de                      \n' +\
      '#@ queue                                              \n'
    else:
      assert( False )

    # load environment
    if( self.m_setup['cluster'] == 'mac_cluster' or  self.m_setup['cluster'] == 'supermuc' ):
      # lrz module system
      self.m_batchScript['environment'] =   \
        '\n'                               +\
        '. /etc/profile \n'                +\
        '. /etc/profile.d/modules.sh \n'

      # tools and seissol 
      self.m_batchScript['environment'] = self.m_batchScript['environment']       +\
        'NUMBER_OF_MPI_RANKS=' + str(self.m_setup['#ranks'])             + '\n'   +\
        'SEISSOL_DIR='   + self.m_setup['seissol_dir']                   + '\n'   +\
        'BENCHMARK_DIR=' + self.m_setup['benchmark_dir']                 + '\n'   +\
        'TIME=$(date +"%y_%m_%d-%H_%M_%S")'                              + '\n'   +\
        'echo \'start time: \'$(date)'                                   + '\n'   +\
        'echo \'loading intel compiler suite 13.1\''                     + '\n'   +\
        'module unload ccomp/intel'                                      + '\n'   +\
        'module unload fortran/intel'                                    + '\n'   +\
        'source /lrz/sys/intel/ifort_131_163/bin/compilervars.sh intel64'+ '\n'   +\
        'source /lrz/sys/intel/mpi_41_0_030/bin64/mpivars.sh'            + '\n'   +\
        'source /lrz/sys/intel/icc_131_163/bin/compilervars.sh intel64'  + '\n'
        #'module unload ccomp/intel'                                      + '\n'   +\
        #'module load ccomp/intel/13.1'                                   + '\n'   +\
        #'module unload fortran/intel'                                    + '\n'   +\
        #'module load fortran/intel/13.1'                                 + '\n'   +\
        #'export OMP_NUM_THREADS=' + str(self.m_setup['#cores per rank']) + '\n \n'
    else:
      assert( False )

    # executable part of the script
    self.m_batchScript['executable'] =\
      'echo \'Measuring kernel performance, Benchmark: ' + self.m_setup['benchmark_dir'] + '\' \n\n' +\
      'echo \'Removing temporary folder: ' + self.m_setup['tmp_dir_name'] + '\'                \n'  +\
      'rm -R ' + self.m_setup['tmp_dir_name'] +                                               '\n'  +\
      'echo \'copying benchmark folder to: ' + self.m_setup['tmp_dir_name'] + '\'              \n'  +\
      'cp -R $BENCHMARK_DIR ' + self.m_setup['tmp_dir_name'] +                                '\n'  +\
      'echo \' moving to temporary benchmark directory\'                                       \n'  +\
      'cd ' + self.m_setup['tmp_dir_name'] +                                                  '\n'  +\
      'pwd                                                                                     \n'  +\
       self.m_setup['preparatory_command'] +                                                  '\n'

    self.m_batchScript['executable'] = self.m_batchScript['executable'] + \
'''
echo 'creating DGPATH-file'
echo ${SEISSOL_DIR}/Maple/ > DGPATH

echo 'creating output directory'
mkdir output

echo 'executing code'
'''
    # itearate over setups
    for l_runSetup in self.m_setup['run_setup']:
      # get the name of the executable by removing the path
      l_executableName = l_runSetup['binary'].split( '/' )[-1]

      # add machine specific run command
      if not 'launcher' in l_runSetup:
        if( self.m_setup['cluster'] == 'mac_cluster' ):
          l_runCommand = 'mpiexec.hydra ./' + l_executableName + ' ' + l_runSetup['config']
        elif( self.m_setup['cluster'] == 'supermuc' ):
          l_runCommand ='poe ./' + l_executableName + ' ' + l_runSetup['config'] + ' -procs ${NUMBER_OF_MPI_RANKS}'
      else:
        l_runCommand = l_runSetup['launcher'] + ' ./' + l_executableName + ' ' + l_runSetup['config']

      # create a symbolic link and run the code
      self.m_batchScript['executable'] = self.m_batchScript['executable'] +\
        'ln -s ' + l_runSetup['binary'] + ' ' + l_executableName   + '\n' +\
        'echo executing ' + l_runCommand                           + '\n' +\
        'date                                                         \n' +\
        l_runCommand                                               + '\n'
     
      # save snapshots locally
      self.m_batchScript['executable'] = self.m_batchScript['executable']             +\
        'date                                                                   \n'   +\
        'echo saving snapshots                                                  \n'   +\
        'mkdir output/' + l_runSetup['output_dir']                           + '\n'   +\
        'mv output/data*receiver* output/' +  l_runSetup['output_dir']       + '\n'   +\
        'mv epik* output/' +  l_runSetup['output_dir']       + '\n'                   +\
        'echo saving config                                                     \n'   +\
        'cp ' + l_runSetup['config'] + ' output/' + l_runSetup['output_dir'] + '\n' +\
         'mv ' + self.m_setup['tmp_dir_name'] + '/output/' + l_runSetup['output_dir'] + ' ' + self.m_setup['benchmark_dir'] + '/../' + l_runSetup['output_dir'] + '_${TIME}' + '\n\n'

    # save the output and remove tmpdir 
    self.m_batchScript['executable'] = self.m_batchScript['executable']    +\
    'rm -R ' + self.m_setup['tmp_dir_name'] + '\n' 

    # open file
    l_file = open( self.m_setup['file_name'] ,'w')
    
    # write the script

    l_file.write( str(self.m_batchScript['header'])      )
    l_file.write( str(self.m_batchScript['environment']) )
    l_file.write( str(self.m_batchScript['executable'])  )

    return self.m_setup['file_name'] 


l_batchGenerator = BatchGenerator()

# script to submit all batch jobs
l_submitScript = open( 'submit_all.sh', 'w' )
l_submitScript.write( '#!/bin/bash\n' )


# generate scripts to test the sparse/dense-switch
l_sparseDenseScripts = True

if l_sparseDenseScripts:
  l_sparseDenseRatios = []

  # define different basis functions
  l_basisFunctions = [4, 10, 20, 35, 56]

  l_batchGenerator.m_setup['#ranks']                      = 16
  l_batchGenerator.m_setup['#cores per rank']             = 1

  # create sparse dense search space
  for l_sparseDenseRatio in numpy.arange(0, 1.01, 0.01).tolist():
    l_sparseDenseRatios = l_sparseDenseRatios + ["%0.2f" % l_sparseDenseRatio]

  # generate a script for each sparse dense ratio
  for l_sparseDenseRatio in l_sparseDenseRatios:
    # reset setup
    l_batchGenerator.m_setup['run_setup'] = []
    l_batchGenerator.m_setup['preparatory_command'] = 'export PATH=/home/hpc/pr63so/di56dok/software/scalasca/scalasca-1.4.3/bin:$PATH \n'

    for l_basisFunction in l_basisFunctions:
      l_batchGenerator.m_setup['name'] = 'performance_generated_kernels_' + l_sparseDenseRatio 
      l_batchGenerator.m_setup['tmp_dir_path'] = '$TMPDIR'
      l_batchGenerator.m_setup['preparatory_command'] = l_batchGenerator.m_setup['preparatory_command'] + 'ln -s matrices/matrices_' + l_sparseDenseRatio + '_' + str(l_basisFunction) + '.xml ' + 'matrices_' + str(l_basisFunction) + '.xml\n'
      l_batchGenerator.m_setup['benchmark_dir'] = '${HOME}/workspace/seissol_benchmarks/performance/LOH1_success'
      l_batchGenerator.m_setup['build_environment'] = 'slurm'
      l_batchGenerator.m_setup['cluster'] = 'mac_cluster'

      l_binary = '${HOME}/workspace/seissol_3d/build/' +\
                   'SeisSol_release_generatedKernels_mpi_scalasca_9_' +  str(l_basisFunction)

      l_outputDir = 'sparse_dense_' + str(l_basisFunction) + '_' + l_sparseDenseRatio

      l_launcher = 'scalasca -analyze mpiexec.hydra -n 16'

      l_parameters = 'parameters_' + str(l_basisFunction) + '.par'

      l_batchGenerator.m_setup['run_setup']  = l_batchGenerator.m_setup['run_setup'] +\
        [{ 'binary': l_binary,
           'config': l_parameters,
           'output_dir': l_outputDir,
           'launcher': l_launcher
          } ]
    l_batchGenerator.setup()
    l_submitScript.write( 'sbatch ' + l_batchGenerator.generateBatchScript()+'\n' )

# generate scripts to test the kernel performance
l_kernelPerformance = False

if l_kernelPerformance:

  # define different run types
  l_runs = ['cycles', 'flops_cycles', 'flops', 'loop', 'overall']

  # define different basis functions
  l_basisFunctions = [4, 10, 20, 35, 56]

  for l_run in l_runs:
    l_binaryDir = l_batchGenerator.m_setup['seissol_dir'] + 'build_' + l_run

    l_batchGenerator.m_setup['tmp_dir_path'] = '$TMPDIR'
    l_batchGenerator.m_setup['benchmark_dir'] = '${HOME}/workspace/seissol_benchmarks/performance/13_07_23_parco'
    l_batchGenerator.m_setup['build_environment'] = 'slurm'
    l_batchGenerator.m_setup['cluster'] = 'mac_cluster'

    for l_basisFunction in l_basisFunctions:
      # iterate over different code versions
      if l_run in ['cycles', 'flops_cycles' 'flops']:
        l_types = ['generatedKernels']
      else:
        l_types = ['generatedKernels', 'classic']

      for l_type in l_types:
        # reset run setup
        l_batchGenerator.m_setup['run_setup'] = []

        l_binary = '${HOME}/workspace/seissol_3d/build_' + l_run + '/' +\
                   'SeisSol_release_' + l_type + '_mpi_9_' +  str(l_basisFunction)
        if( l_run in 'cycles', 'flops_cycles' ):
          l_config = 'PARAMETERS_cycles.par'
        else:
          l_config = 'PARAMETERS_' + str(l_basisFunction) + '.par'

        l_outputDir = l_run + '_'  + l_type + '_' + str(l_basisFunction)

        l_batchGenerator.m_setup['run_setup']  = l_batchGenerator.m_setup['run_setup'] +\
          [{ 'binary': l_binary,
             'config': l_config,
             'output_dir': l_outputDir          
            } ]

        l_batchGenerator.m_setup['name'] = l_run + '_' + l_type + '_' + str(l_basisFunction)

        l_batchGenerator.setup()
        l_submitScript.write( 'sbatch ' + l_batchGenerator.generateBatchScript()+'\nsleep 15\n' )

# generate verification runs
l_verification = False

if l_verification:

  # define different benchmarks
  l_benchmarks = ['Tetra_Elastic', 'LOH1_success', 'Tetra_LOH4_generated_kernels']

  # define different basis functions
  l_basisFunctions = [4, 10, 20, 35, 56]

  for l_benchmark in l_benchmarks:
    l_batchGenerator.m_setup['benchmark_dir'] = '${HOME}/workspace/seissol_benchmarks/verification/'+l_benchmark
    l_batchGenerator.m_setup['build_environment'] = 'slurm'
    l_batchGenerator.m_setup['cluster'] = 'mac_cluster'

    for l_basisFunction in l_basisFunctions:
      for l_type in ['classic', 'generatedKernels']:
        # reset run setup
        l_batchGenerator.m_setup['run_setup'] = []

        if( l_benchmark == 'Tetra_Elastic' ):
          l_batchGenerator.m_setup['tmp_dir_path'] = '$TMPDIR'
          if( l_type == 'generatedKernels' ):
            l_batchGenerator.m_setup['#ranks']                      = 2
            l_batchGenerator.m_setup['#cores per rank']             = 8
          else:
            l_batchGenerator.m_setup['#ranks']                      = 16
            l_batchGenerator.m_setup['#cores per rank']             = 1
        else:
          l_batchGenerator.m_setup['tmp_dir_path'] = '$HOME'
          if( l_type == 'generatedKernels' ):
            l_batchGenerator.m_setup['#ranks']                      = 16
            l_batchGenerator.m_setup['#cores per rank']             = 8
          else:
            l_batchGenerator.m_setup['#ranks']                      = 128
            l_batchGenerator.m_setup['#cores per rank']             = 1

        if( l_type == 'generatedKernels' ):
          l_binary = '${HOME}/workspace/seissol_3d/build/'+\
                     'SeisSol_release_' + l_type + '_hybrid_9_' +  str(l_basisFunction)
        else:
          l_binary = '${HOME}/workspace/seissol_3d/build/'+\
                     'SeisSol_release_' + l_type + '_mpi_9_' +  str(l_basisFunction)

        l_batchGenerator.m_setup['run_setup']  = l_batchGenerator.m_setup['run_setup'] +\
          [{ 'binary': l_binary,
             'config': 'parameters_' + str(l_basisFunction) + '.par',
             'output_dir': l_benchmark + '_'  + l_type + '_' + str(l_basisFunction)       
            } ]

        l_batchGenerator.m_setup['name'] = l_benchmark  + '_' + l_type + '_' + str(l_basisFunction)

        l_batchGenerator.setup()
        l_submitScript.write( 'sbatch ' + l_batchGenerator.generateBatchScript()+'\nsleep 15\n' )
