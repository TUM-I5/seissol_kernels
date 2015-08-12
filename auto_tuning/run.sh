#!/bin/bash
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
# Entry point for the auto-tuning runs.

# Global config
TUNE_CELLS[2]=50000
TUNE_CELLS[3]=50000
TUNE_CELLS[4]=50000
TUNE_CELLS[5]=50000
TUNE_CELLS[6]=50000
TUNE_CELLS[7]=50000

TUNE_TS[2]=25000
TUNE_TS[3]=15000
TUNE_TS[4]=8500
TUNE_TS[5]=4000
TUNE_TS[6]=1750
TUNE_TS[7]=750

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-m MODE -c CONFIGURATIONS -a GENERATED_KERNELS -k SEISSOL_KERNELS_DIR -g GEMM_CODE_GEN_DIR -o OUTPUT_DIR -t WORK_DIR -r REPEATS -d]
Runs the auto-tuning of the kernels.
     -h                     display this help and exit.
     -m MODE                mode to run:
                             1: generate code
                             2: build code
                             3: run unit test
                             4: auto-tune
                             5: 2 + 3
                             6: 2 + 4
                             7: 1 + 2 + 3
                             8: 1 + 2 + 3 + 4
     -c CONFIGURATIONS      path to the directory containing the sparse-dense configurations.
     -a GENERATED_KERNELS   Kernels to tune for:  swsm, dwsm, ssnb, dsnb, sknc, dknc, shsw, dhsw
     -k SEISSOL_KERNELS_DIR path to the seissol_kernels module.
     -g GEMM_CODE_GEN_EXE   path to the GemmCodeGenerator executable.
     -o OUTPUT_DIR          path to the dictory where the output goes.
     -w WORK_DIR            path to the working directory (optional).
     -r REPEATS             number of repeats for the tuning runs (default: 1)
     -d                     delete working directories after usage: generated code (1, 7, 8), builds (2, 5, 6, 7, 8).
EOF
exit 1
}

# Print with data
echo_date() {
    echo -e `date +%y/%m/%d\ %H:%M:%S` $*
}

# Generates the code for a single configuration
#
# @param $1 ${SEISSOL_KERNELS_DIR}
# @param $2 {MATRICES_DIR}
# @param $3 ${CONFIG}
# @param $4 ${GEMM_CODE_DIR}
# @param $5 ${CODE_DIR}
#

generate_code() {
  echo_date "generating code for $3"
  if ! python $1/preprocessing/scripts/offlineAssembly.py --generateMatrixKernels $2 $3 $4 $5
  then
    echo_date "offline assembly failed, exiting: $1 $2 $3 $4 $5" >&2
    exit 1
  fi
}

# Builds the code.
#
# @param $1 ${GENERATED_KERNELS}
# @param $2 ${ORDER}
# @param $3 ${SEISSOL_KERNELS_DIR}
# @param $4 ${CODE_DIR}
# @param $5 ${BUILD_DIR}
# @param $6 ${SEISSOL_PROXY} (optional)
build() {
  # store call dir
  CALL_DIR=$(pwd)

  if [ -z $6 ]; then
    PROXY=OFF
  else
    PROXY=$6
  fi

  echo_date "building code for $4"

  # jump to build dir
  cd $5

  # call cmake
  if ! cmake $3 -DGENERATED_KERNELS=$1 -DCONVERGENCE_ORDER=$2 -DGENERATED_CODE=$4 -DSEISSOL_PROXY=$PROXY
  then
    echo_date "cmake failed, exiting: $1 $2 $3 $4 $5 $6" >&2
    exit 1
  fi

  if ! make -j
  then
    echo_date "make failed, exiting: $1 $2 $3 $4 $5 $6" >&2
    exit 1
  fi

  cd ${CALL_DIR}
}

# Runs the unit tests.
#
# @param $1 ${BUILD_DIR}
unit_tests() {
  if ! $1/time_unit_tests; then
    echo_date "time unit tests failed, exiting: $1" >&2
    exit 1
  fi

  if ! $1/volume_unit_tests; then
    echo_date "volume unit tests failed, exiting: $1" >&2
    exit 1
  fi

  if ! $1/flux_unit_tests; then
    echo_date "flux unit tests failed, exiting: $1" >&2
    exit 1
  fi
}

# Runs the tuning
#
# @param $1 ${BUILD_DIR}
# @param $2 #cells
# @param $3 #timesteps
# @param $4 kernel (all, local, neigh, ader, vol, bndlocal)
# @param $5 number of tuning repeats
tune() {
  CURRENT_REPEAT=0
  while [ $CURRENT_REPEAT -lt $5 ]; do
    if ! $1/seissol_proxy $2 $3 $4; then
      echo_date "seissol proxy failed, exiting: $1 $2 $3 $4" >&2
      exit
    fi
    let CURRENT_REPEAT=CURRENT_REPEAT+1 
  done
}

#
# Parse command line arguments.
#

OPTIND=1
while getopts "hm:c:a:k:g:o:w:r:d" opt; do
    case "$opt" in
        m) MODE=$OPTARG
            # set up the switches
            if [ $MODE == 1 ] || [ $MODE == 7 ] || [ $MODE == 8 ];
            then
              GENERATE_SWITCH=true
            fi

            if [ $MODE == 2 ] || [ $MODE == 5 ] || [ $MODE == 6 ] || [ $MODE == 7 ] || [ $MODE == 8 ];
            then
              BUILD_SWITCH=true
            fi

            if [ $MODE == 3 ] || [ $MODE == 5 ] || [ $MODE == 7 ] || [ $MODE == 8 ];
            then
              UT_SWITCH=true
            fi

            if [ $MODE == 4 ] || [ $MODE == 6 ] || [ $MODE == 8 ];
            then
              TUNE_SWITCH=true
            fi
            ;;
        c) CONFIGURATIONS=$OPTARG
            ;;
        a) GENERATED_KERNELS=$OPTARG
            ;;
        k) SEISSOL_KERNELS_DIR=$OPTARG
            ;;
        g) GEMM_CODE_GEN_EXE=$OPTARG
            ;;
        o) OUTPUT_DIR=$OPTARG
            ;;
        w) WORK_DIR=$OPTARG
            ;;
        r) TUNE_REPEATS=$OPTARG
            ;;
        d) DELETE=true
            ;;
        *)
            show_help >&2
            ;;
    esac
done

shift $((OPTIND-1))

if [ -z "${TUNE_REPEATS}" ];
then
  TUNE_REPEATS=1
fi

if [ -z "${MODE}"                ] ||
   [ -z "${CONFIGURATIONS}"      ] ||
   [ -z "${GENERATED_KERNELS}"   ] ||
   [ -z "${SEISSOL_KERNELS_DIR}" ] ||
   [ -z "${GEMM_CODE_GEN_EXE}"   ] ||
   [ -z "${OUTPUT_DIR}"          ];
then
    show_help >&2
else
  echo_date '---------------------------------------------------------'
  echo_date '| Welcome to the auto-tuning script of seissol_kernels! |'
  echo_date '---------------------------------------------------------'
  echo
  echo_date 'We are using the following configuration:'
  echo_date "\t Code generation enabled?  ${GENERATE_SWITCH}"
  echo_date "\t Build enabled?            ${BUILD_SWITCH}"
  echo_date "\t Unit tests enabled?       ${UT_SWITCH}"
  echo_date "\t Auto tuning enabled?      ${TUNE_SWITCH}"
  if [ $TUNE_SWITCH ]; then
  echo_date "\t Number of tuning repeats: ${TUNE_REPEATS}"
  fi
  echo_date "\t CONFIGURATIONS:           ${CONFIGURATIONS}"
  echo_date "\t GENERATED_KERNELS:        ${GENERATED_KERNELS}"
  echo_date "\t SEISSOL_KERNELS_DIR:      ${SEISSOL_KERNELS_DIR}"
  echo_date "\t GEMM_CODE_GEN_EXE:        ${GEMM_CODE_GEN_EXE}"
  echo_date "\t OUTPUT_DIR:               ${OUTPUT_DIR}"
  if [ "$WORK_DIR" != "" ]; then
    echo_date "\t WORK_DIR:               ${WORK_DIR}\n"
  fi
  echo_date "\t Deleting enabled?         ${DELETE}"
fi

#
# Set up the infrastructure for running the code.
#
if [[ -z "$WORK_DIR" ]]; then
  WORK_DIR=$(mktemp -d)
fi
echo_date "using working directory: ${WORK_DIR}"

# use Intel compilers
export CC=icc
export CXX=icpc
export CPP=icpc

# get full path of everything
CONFIGURATIONS=$(readlink -f $CONFIGURATIONS)
SEISSOL_KERNELS_DIR=$(readlink -f $SEISSOL_KERNELS_DIR)
GEMM_CODE_GEN_EXE=$(readlink -f $GEMM_CODE_GEN_EXE)
OUTPUT_DIR=$(readlink -f $OUTPUT_DIR)
WORK_DIR=$(readlink -f $WORK_DIR)

LOG_DIR=${OUTPUT_DIR}/logs
mkdir $LOG_DIR
echo_date "using log directory: ${LOG_DIR}"

GENERATED_DIR=${WORK_DIR}/generated_code
if [ $GENERATE_SWITCH ]; then
  mkdir ${GENERATED_DIR}
fi
BUILDS_DIR=${WORK_DIR}/builds
if [ $BUILD_SWITCH ]; then
  mkdir ${BUILDS_DIR}
fi

#
# Iterate over the configurations
#
echo

# Repeat the tuning runs if requested
CURRENT_REPEAT=0
while [ $CURRENT_REPEAT -lt $TUNE_REPEATS ]; do
  echo_date "iteration: ${CURRENT_REPEAT}"

  echo_date "iterating over configurations.."
  for CONFIG in ${CONFIGURATIONS}/*/; do
    # get base name
    CONFIG_BASE=$(basename ${CONFIG})

    echo_date "config: ${CONFIG_BASE}"

    CODE_DIR=${GENERATED_DIR}/gen_$CONFIG_BASE

    if [ $GENERATE_SWITCH ] && [ $CURRENT_REPEAT -lt 1 ]; then
      # generate code for this configuration
      mkdir -p ${CODE_DIR}
      echo_date "generating code"
      generate_code ${SEISSOL_KERNELS_DIR} ${SEISSOL_KERNELS_DIR}/preprocessing/matrices ${CONFIG} ${GEMM_CODE_GEN_EXE} ${CODE_DIR} >> ${LOG_DIR}/gen_${CONFIG_BASE}.log  2>> ${LOG_DIR}/gen_${CONFIG_BASE}.err
    fi

    if [ $BUILD_SWITCH ] || [ $UT_SWITCH ] || [ $TUNE_SWITCH ]; then
      echo -n -e `date +%y/%m/%d\ %H:%M:%S` $* 'iterating over orders: '
      for ORDER in {2..7}; do
        echo -n "${ORDER}.."

        BUILD_DIR=${BUILDS_DIR}/${CONFIG_BASE}/${ORDER}
        if [ $BUILD_SWITCH ] && [ $CURRENT_REPEAT -lt 1 ]; then
          # build this setting
          mkdir -p $BUILD_DIR

          build ${GENERATED_KERNELS} ${ORDER} ${SEISSOL_KERNELS_DIR} ${CODE_DIR} ${BUILD_DIR} ${SEISSOL_KERNELS_DIR}/proxy >> ${LOG_DIR}/build_${CONFIG_BASE}_${ORDER}.log 2>> ${LOG_DIR}/build_${CONFIG_BASE}_${ORDER}.err
        fi

        if [ $UT_SWITCH ]  && [ $CURRENT_REPEAT -lt 1 ]; then
          # run the unit tests
          unit_tests ${BUILD_DIR}  >> ${LOG_DIR}/ut_${CONFIG_BASE}_${ORDER}.log 2>> ${LOG_DIR}/ut_${CONFIG_BASE}_${ORDER}.err
        fi

        if [ $TUNE_SWITCH ]; then
          echo "#cells ${TUNE_CELLS[${ORDER}]}" >> ${LOG_DIR}/tune_${CONFIG_BASE}_${ORDER}.log
          echo "#ts ${TUNE_TS[${ORDER}]}"       >> ${LOG_DIR}/tune_${CONFIG_BASE}_${ORDER}.log
          # run the tuning runs
          tune ${BUILD_DIR} ${TUNE_CELLS[${ORDER}]} ${TUNE_TS[${ORDER}]} all 1  >> ${LOG_DIR}/tune_${CONFIG_BASE}_${ORDER}.log 2>> ${LOG_DIR}/tune_${CONFIG_BASE}_${ORDER}.err
        fi

        # clean up build
        if [ $DELETE ] && [ $BUILD_SWITCH ]  && [ $CURRENT_REPEAT -lt 1 ]; then
          rm -r $BUILD_DIR
        fi
      done
      echo

      # clean up generated code
      if [ $DELETE ] && [ $GENERATE_SWITCH ]  && [ $CURRENT_REPEAT -lt 1 ]; then
        rm -r ${CODE_DIR}
      fi
    fi

  done

let CURRENT_REPEAT=CURRENT_REPEAT+1 
done

echo_date '-------------------------'
echo_date '| Auto-tuning finished. |'
echo_date '-------------------------'
