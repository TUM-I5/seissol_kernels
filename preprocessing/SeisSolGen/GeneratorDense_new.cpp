/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
*                                                                             *
* This file is part of the SeisSolGen project. For conditions of              *
* distribution and use, please see the copyright notice                       *
* in its root directory.                                                      *
******************************************************************************/
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <assert.h>

#include "GeneratorDense.hpp"

namespace seissolgen {

  GeneratorDense::GeneratorDense() : bGenerateExitForCK_(false), bAdd_(true) {
  }

  GeneratorDense::GeneratorDense(bool bGenerateExitForCK, bool bAdd) : bGenerateExitForCK_(bGenerateExitForCK), bAdd_(bAdd) {
  }

  void sse_inner_blocked_kernel_6(std::stringstream& codestream, int lda, bool alignA) {
    codestream << "    b_0 = _mm_loaddup_pd(b0);" << std::endl;
    codestream << "    b_1 = _mm_loaddup_pd(b1);" << std::endl;
    codestream << "    b_2 = _mm_loaddup_pd(b2);" << std::endl << std::endl;
    codestream << "    b0++; b1++; b2++;" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_0 = _mm_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_0 = _mm_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << 2 << ";" << std::endl;
    codestream << "    c_0_0 = _mm_add_pd(c_0_0, _mm_mul_pd(a_0, b_0));" << std::endl;
    codestream << "    c_0_1 = _mm_add_pd(c_0_1, _mm_mul_pd(a_0, b_1));" << std::endl;
    codestream << "    c_0_2 = _mm_add_pd(c_0_2, _mm_mul_pd(a_0, b_2));" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_1 = _mm_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_1 = _mm_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << 2 << ";" << std::endl;
    codestream << "    c_1_0 = _mm_add_pd(c_1_0, _mm_mul_pd(a_1, b_0));" << std::endl;
    codestream << "    c_1_1 = _mm_add_pd(c_1_1, _mm_mul_pd(a_1, b_1));" << std::endl;
    codestream << "    c_1_2 = _mm_add_pd(c_1_2, _mm_mul_pd(a_1, b_2));" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_2 = _mm_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_2 = _mm_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << lda - 4 << ";" << std::endl;
    codestream << "    c_2_0 = _mm_add_pd(c_2_0, _mm_mul_pd(a_2, b_0));" << std::endl;
    codestream << "    c_2_1 = _mm_add_pd(c_2_1, _mm_mul_pd(a_2, b_1));" << std::endl;
    codestream << "    c_2_2 = _mm_add_pd(c_2_2, _mm_mul_pd(a_2, b_2));" << std::endl << std::endl;
  }

  void sse_inner_blocked_kernel_4(std::stringstream& codestream, int lda, bool alignA) {
    codestream << "    b_0 = _mm_loaddup_pd(b0);" << std::endl;
    codestream << "    b_1 = _mm_loaddup_pd(b1);" << std::endl;
    codestream << "    b_2 = _mm_loaddup_pd(b2);" << std::endl << std::endl;
    codestream << "    b0++; b1++; b2++;" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_0 = _mm_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_0 = _mm_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << 2 << ";" << std::endl;
    codestream << "    c_0_0 = _mm_add_pd(c_0_0, _mm_mul_pd(a_0, b_0));" << std::endl;
    codestream << "    c_0_1 = _mm_add_pd(c_0_1, _mm_mul_pd(a_0, b_1));" << std::endl;
    codestream << "    c_0_2 = _mm_add_pd(c_0_2, _mm_mul_pd(a_0, b_2));" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_1 = _mm_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_1 = _mm_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << lda - 2 << ";" << std::endl;
    codestream << "    c_1_0 = _mm_add_pd(c_1_0, _mm_mul_pd(a_1, b_0));" << std::endl;
    codestream << "    c_1_1 = _mm_add_pd(c_1_1, _mm_mul_pd(a_1, b_1));" << std::endl;
    codestream << "    c_1_2 = _mm_add_pd(c_1_2, _mm_mul_pd(a_1, b_2));" << std::endl << std::endl;
  }

  void sse_inner_blocked_kernel_2(std::stringstream& codestream, int lda, bool alignA) {
    codestream << "    b_0 = _mm_loaddup_pd(b0);" << std::endl;
    codestream << "    b_1 = _mm_loaddup_pd(b1);" << std::endl;
    codestream << "    b_2 = _mm_loaddup_pd(b2);" << std::endl << std::endl;
    codestream << "    b0++; b1++; b2++;" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_0 = _mm_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_0 = _mm_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << lda << ";" << std::endl;
    codestream << "    c_0_0 = _mm_add_pd(c_0_0, _mm_mul_pd(a_0, b_0));" << std::endl;
    codestream << "    c_0_1 = _mm_add_pd(c_0_1, _mm_mul_pd(a_0, b_1));" << std::endl;
    codestream << "    c_0_2 = _mm_add_pd(c_0_2, _mm_mul_pd(a_0, b_2));" << std::endl << std::endl;
  }

  void sse_inner_blocked_kernel_1(std::stringstream& codestream, int lda, bool alignA) {
    codestream << "    b_0 = _mm_load_sd(b0);" << std::endl;
    codestream << "    b_1 = _mm_load_sd(b1);" << std::endl;
    codestream << "    b_2 = _mm_load_sd(b2);" << std::endl << std::endl;
    codestream << "    b0++; b1++; b2++;" << std::endl << std::endl;

    codestream << "    a_0 = _mm_load_sd(a0);" << std::endl;

    codestream << "    a0+=" << lda << ";" << std::endl;
    codestream << "    c_0_0 = _mm_add_sd(c_0_0, _mm_mul_sd(a_0, b_0));" << std::endl;
    codestream << "    c_0_1 = _mm_add_sd(c_0_1, _mm_mul_sd(a_0, b_1));" << std::endl;
    codestream << "    c_0_2 = _mm_add_sd(c_0_2, _mm_mul_sd(a_0, b_2));" << std::endl << std::endl;
  }

  void avx_inner_blocked_kernel_12(std::stringstream& codestream, int lda, bool alignA) {
    codestream << "    b_0 = _mm256_broadcast_sd(b0);" << std::endl;
    codestream << "    b_1 = _mm256_broadcast_sd(b1);" << std::endl;
    codestream << "    b_2 = _mm256_broadcast_sd(b2);" << std::endl << std::endl;
    codestream << "    b0++; b1++; b2++;" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_0 = _mm256_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_0 = _mm256_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << 4 << ";" << std::endl;
    codestream << "#ifdef __AVX2__" << std::endl;
    codestream << "    c_0_0 = _mm256_fmadd_pd(a_0, b_0, c_0_0);" << std::endl;
    codestream << "    c_0_1 = _mm256_fmadd_pd(a_0, b_1, c_0_1);" << std::endl;
    codestream << "    c_0_2 = _mm256_fmadd_pd(a_0, b_2, c_0_2);" << std::endl;
    codestream << "#else" << std::endl;
    codestream << "    c_0_0 = _mm256_add_pd(c_0_0, _mm256_mul_pd(a_0, b_0));" << std::endl;
    codestream << "    c_0_1 = _mm256_add_pd(c_0_1, _mm256_mul_pd(a_0, b_1));" << std::endl;
    codestream << "    c_0_2 = _mm256_add_pd(c_0_2, _mm256_mul_pd(a_0, b_2));" << std::endl;
    codestream << "#endif" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_1 = _mm256_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_1 = _mm256_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << 4 << ";" << std::endl;
    codestream << "#ifdef __AVX2__" << std::endl;
    codestream << "    c_0_0 = _mm256_fmadd_pd(a_1, b_0, c_1_0);" << std::endl;
    codestream << "    c_0_1 = _mm256_fmadd_pd(a_1, b_1, c_1_1);" << std::endl;
    codestream << "    c_0_2 = _mm256_fmadd_pd(a_1, b_2, c_1_2);" << std::endl;
    codestream << "#else" << std::endl;
    codestream << "    c_1_0 = _mm256_add_pd(c_1_0, _mm256_mul_pd(a_1, b_0));" << std::endl;
    codestream << "    c_1_1 = _mm256_add_pd(c_1_1, _mm256_mul_pd(a_1, b_1));" << std::endl;
    codestream << "    c_1_2 = _mm256_add_pd(c_1_2, _mm256_mul_pd(a_1, b_2));" << std::endl << std::endl;
    codestream << "#endif" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_2 = _mm256_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_2 = _mm256_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << lda - 8 << ";" << std::endl;
    codestream << "#ifdef __AVX2__" << std::endl;
    codestream << "    c_0_0 = _mm256_fmadd_pd(a_2, b_0, c_2_0);" << std::endl;
    codestream << "    c_0_1 = _mm256_fmadd_pd(a_2, b_1, c_2_1);" << std::endl;
    codestream << "    c_0_2 = _mm256_fmadd_pd(a_2, b_2, c_2_2);" << std::endl;
    codestream << "#else" << std::endl;
    codestream << "    c_2_0 = _mm256_add_pd(c_2_0, _mm256_mul_pd(a_2, b_0));" << std::endl;
    codestream << "    c_2_1 = _mm256_add_pd(c_2_1, _mm256_mul_pd(a_2, b_1));" << std::endl;
    codestream << "    c_2_2 = _mm256_add_pd(c_2_2, _mm256_mul_pd(a_2, b_2));" << std::endl << std::endl;
    codestream << "#endif" << std::endl << std::endl;
  }

  void avx_inner_blocked_kernel_8(std::stringstream& codestream, int lda, bool alignA) {
    codestream << "    b_0 = _mm256_broadcast_sd(b0);" << std::endl;
    codestream << "    b_1 = _mm256_broadcast_sd(b1);" << std::endl;
    codestream << "    b_2 = _mm256_broadcast_sd(b2);" << std::endl << std::endl;
    codestream << "    b0++; b1++; b2++;" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_0 = _mm256_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_0 = _mm256_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << 4 << ";" << std::endl;
    codestream << "#ifdef __AVX2__" << std::endl;
    codestream << "    c_0_0 = _mm256_fmadd_pd(a_0, b_0, c_0_0);" << std::endl;
    codestream << "    c_0_1 = _mm256_fmadd_pd(a_0, b_1, c_0_1);" << std::endl;
    codestream << "    c_0_2 = _mm256_fmadd_pd(a_0, b_2, c_0_2);" << std::endl;
    codestream << "#else" << std::endl;
    codestream << "    c_0_0 = _mm256_add_pd(c_0_0, _mm256_mul_pd(a_0, b_0));" << std::endl;
    codestream << "    c_0_1 = _mm256_add_pd(c_0_1, _mm256_mul_pd(a_0, b_1));" << std::endl;
    codestream << "    c_0_2 = _mm256_add_pd(c_0_2, _mm256_mul_pd(a_0, b_2));" << std::endl;
    codestream << "#endif" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_1 = _mm256_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_1 = _mm256_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << lda - 4  << ";" << std::endl;
    codestream << "#ifdef __AVX2__" << std::endl;
    codestream << "    c_0_0 = _mm256_fmadd_pd(a_1, b_0, c_1_0);" << std::endl;
    codestream << "    c_0_1 = _mm256_fmadd_pd(a_1, b_1, c_1_1);" << std::endl;
    codestream << "    c_0_2 = _mm256_fmadd_pd(a_1, b_2, c_1_2);" << std::endl;
    codestream << "#else" << std::endl;
    codestream << "    c_1_0 = _mm256_add_pd(c_1_0, _mm256_mul_pd(a_1, b_0));" << std::endl;
    codestream << "    c_1_1 = _mm256_add_pd(c_1_1, _mm256_mul_pd(a_1, b_1));" << std::endl;
    codestream << "    c_1_2 = _mm256_add_pd(c_1_2, _mm256_mul_pd(a_1, b_2));" << std::endl << std::endl;
    codestream << "#endif" << std::endl << std::endl;
  }

  void avx_inner_blocked_kernel_4(std::stringstream& codestream, int lda, bool alignA) {
    codestream << "    b_0 = _mm256_broadcast_sd(b0);" << std::endl;
    codestream << "    b_1 = _mm256_broadcast_sd(b1);" << std::endl;
    codestream << "    b_2 = _mm256_broadcast_sd(b2);" << std::endl << std::endl;
    codestream << "    b0++; b1++; b2++;" << std::endl << std::endl;

    if (alignA == true) {
      codestream << "    a_0 = _mm256_load_pd(a0);" << std::endl;
    } else {
      codestream << "    a_0 = _mm256_loadu_pd(a0);" << std::endl;
    }

    codestream << "    a0+=" << lda << ";" << std::endl;
    codestream << "#ifdef __AVX2__" << std::endl;
    codestream << "    c_0_0 = _mm256_fmadd_pd(a_0, b_0, c_0_0);" << std::endl;
    codestream << "    c_0_1 = _mm256_fmadd_pd(a_0, b_1, c_0_1);" << std::endl;
    codestream << "    c_0_2 = _mm256_fmadd_pd(a_0, b_2, c_0_2);" << std::endl;
    codestream << "#else" << std::endl;
    codestream << "    c_0_0 = _mm256_add_pd(c_0_0, _mm256_mul_pd(a_0, b_0));" << std::endl;
    codestream << "    c_0_1 = _mm256_add_pd(c_0_1, _mm256_mul_pd(a_0, b_1));" << std::endl;
    codestream << "    c_0_2 = _mm256_add_pd(c_0_2, _mm256_mul_pd(a_0, b_2));" << std::endl;
    codestream << "#endif" << std::endl << std::endl;
  }

  void avx_inner_blocked_kernel_1(std::stringstream& codestream, int lda, bool alignA) {
    codestream << "    b_0 = _mm256_load_sd(b0);" << std::endl;
    codestream << "    b_1 = _mm256_load_sd(b1);" << std::endl;
    codestream << "    b_2 = _mm256_load_sd(b2);" << std::endl << std::endl;
    codestream << "    b0++; b1++; b2++;" << std::endl << std::endl;

    codestream << "    a_0 = _mm256_load_sd(a0);" << std::endl;

    codestream << "    a0+=" << lda << ";" << std::endl;
    codestream << "#ifdef __AVX2__" << std::endl;
    codestream << "    c_0_0 = _mm256_fmadd_sd(a_0, b_0, c_0_0);" << std::endl;
    codestream << "    c_0_1 = _mm256_fmadd_sd(a_0, b_1, c_0_1);" << std::endl;
    codestream << "    c_0_2 = _mm256_fmadd_sd(a_0, b_2, c_0_2);" << std::endl;
    codestream << "#else" << std::endl;
    codestream << "    c_0_0 = _mm256_add_sd(c_0_0, _mm256_mul_sd(a_0, b_0));" << std::endl;
    codestream << "    c_0_1 = _mm256_add_sd(c_0_1, _mm256_mul_sd(a_0, b_1));" << std::endl;
    codestream << "    c_0_2 = _mm256_add_sd(c_0_2, _mm256_mul_sd(a_0, b_2));" << std::endl;
    codestream << "#endif" << std::endl << std::endl;
  }

  void mic_inner_blocked_kernel_for_56(std::stringstream& codestream, int lda) {
    codestream << "    b_0 = _mm512_extload_pd(b0, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);" << std::endl;
    codestream << "    b_1 = _mm512_extload_pd(b1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);" << std::endl;
    codestream << "    b_2 = _mm512_extload_pd(b2, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);" << std::endl << std::endl;

    //codestream << "    _mm_prefetch((const char*)a0+16, _MM_HINT_T0);" << std::endl;
    //codestream << "    _mm_prefetch((const char*)a0+" << 2*lda << ", _MM_HINT_T1);" << std::endl;
    codestream << "    a_0 = _mm512_load_pd(a0);" << std::endl;
    codestream << "    c_0_0 = _mm512_fmadd_pd(a_0, b_0, c_0_0);" << std::endl;
    codestream << "    c_0_1 = _mm512_fmadd_pd(a_0, b_1, c_0_1);" << std::endl;
    codestream << "    c_0_2 = _mm512_fmadd_pd(a_0, b_2, c_0_2);" << std::endl;
    codestream << "    a0 += 8;" << std::endl;

    //codestream << "    _mm_prefetch((const char*)a0+16, _MM_HINT_T0);" << std::endl;
    //codestream << "    _mm_prefetch((const char*)a0+" << 2*lda << ", _MM_HINT_T1);" << std::endl;
    codestream << "    a_1 = _mm512_load_pd(a0);" << std::endl;
    codestream << "    c_1_0 = _mm512_fmadd_pd(a_1, b_0, c_1_0);" << std::endl;
    codestream << "    c_1_1 = _mm512_fmadd_pd(a_1, b_1, c_1_1);" << std::endl;
    codestream << "    c_1_2 = _mm512_fmadd_pd(a_1, b_2, c_1_2);" << std::endl;
    codestream << "    a0 += 8;" << std::endl;

    //codestream << "    _mm_prefetch((const char*)a0+16, _MM_HINT_T0);" << std::endl;
    //codestream << "    _mm_prefetch((const char*)a0+" << 2*lda << ", _MM_HINT_T1);" << std::endl;
    codestream << "    a_2 = _mm512_load_pd(a0);" << std::endl;
    codestream << "    c_2_0 = _mm512_fmadd_pd(a_2, b_0, c_2_0);" << std::endl;
    codestream << "    c_2_1 = _mm512_fmadd_pd(a_2, b_1, c_2_1);" << std::endl;
    codestream << "    c_2_2 = _mm512_fmadd_pd(a_2, b_2, c_2_2);" << std::endl;
    codestream << "    a0 += 8;" << std::endl;

    //codestream << "    _mm_prefetch((const char*)a0+16, _MM_HINT_T0);" << std::endl;
    //codestream << "    _mm_prefetch((const char*)a0+" << 2*lda << ", _MM_HINT_T1);" << std::endl;
    codestream << "    a_3 = _mm512_load_pd(a0);" << std::endl;
    codestream << "    c_3_0 = _mm512_fmadd_pd(a_3, b_0, c_3_0);" << std::endl;
    codestream << "    c_3_1 = _mm512_fmadd_pd(a_3, b_1, c_3_1);" << std::endl;
    codestream << "    c_3_2 = _mm512_fmadd_pd(a_3, b_2, c_3_2);" << std::endl;
    codestream << "    a0 += 8;" << std::endl;

    //codestream << "    _mm_prefetch((const char*)a0+16, _MM_HINT_T0);" << std::endl;
    //codestream << "    _mm_prefetch((const char*)a0+" << 2*lda << ", _MM_HINT_T1);" << std::endl;
    codestream << "    a_4 = _mm512_load_pd(a0);" << std::endl;
    codestream << "    c_4_0 = _mm512_fmadd_pd(a_4, b_0, c_4_0);" << std::endl;
    codestream << "    c_4_1 = _mm512_fmadd_pd(a_4, b_1, c_4_1);" << std::endl;
    codestream << "    c_4_2 = _mm512_fmadd_pd(a_4, b_2, c_4_2);" << std::endl;
    codestream << "    a0 += 8;" << std::endl;

    //codestream << "    _mm_prefetch((const char*)a0+16, _MM_HINT_T0);" << std::endl;
    //codestream << "    _mm_prefetch((const char*)a0+" << 2*lda << ", _MM_HINT_T1);" << std::endl;
    codestream << "    a_5 = _mm512_load_pd(a0);" << std::endl;
    codestream << "    c_5_0 = _mm512_fmadd_pd(a_5, b_0, c_5_0);" << std::endl;
    codestream << "    c_5_1 = _mm512_fmadd_pd(a_5, b_1, c_5_1);" << std::endl;
    codestream << "    c_5_2 = _mm512_fmadd_pd(a_5, b_2, c_5_2);" << std::endl;
    codestream << "    a0 += 8;" << std::endl;

    //codestream << "    _mm_prefetch((const char*)a0+16, _MM_HINT_T0);" << std::endl;
    //codestream << "    _mm_prefetch((const char*)a0+" << 2*lda << ", _MM_HINT_T1);" << std::endl;
    codestream << "    a_6 = _mm512_load_pd(a0);" << std::endl;
    codestream << "    c_6_0 = _mm512_fmadd_pd(a_6, b_0, c_6_0);" << std::endl;
    codestream << "    c_6_1 = _mm512_fmadd_pd(a_6, b_1, c_6_1);" << std::endl;
    codestream << "    c_6_2 = _mm512_fmadd_pd(a_6, b_2, c_6_2);" << std::endl;
    codestream << "    a0 += " << lda - 48 << ";" << std::endl << std::endl;

    codestream << "    b0++;" << std::endl;
    codestream << "    b1++;" << std::endl;
    codestream << "    b2++;" << std::endl;
  }

  std::string GeneratorDense::generate_dense(bool bIsColMajor, int M, int N, int K, int lda, int ldb, int ldc) {
    std::stringstream codestream;
    int maxM = 0;
    int remainder = 0;
    int block_end = 0;
    int kb = 0;
    bool alignA = false;
    bool alignB = false;
    bool alignC = false;

    // @TODO add unrolling factors depending on K
    int unroll_factor = 8;

#ifdef DEBUG
    std::cout << "Generating dense matrix multiplication" << std::endl;
    std::cout << "M=" << M << " N=" << N << " K=" << K << std::endl;
    std::cout << "lda=" << lda << " ldb=" << ldb << " ldc=" << ldc << std::endl;
#endif

    assert(N % 3 == 0);
    assert(bIsColMajor == true);

    // generating SSE3 code
    codestream << "#if defined(__SSE3__) && !defined(__AVX__)" << std::endl;
    codestream << "__m128d c_0_0;" << std::endl;
    codestream << "__m128d c_0_1;" << std::endl;
    codestream << "__m128d c_0_2;" << std::endl;
    codestream << "__m128d c_1_0;" << std::endl;
    codestream << "__m128d c_1_1;" << std::endl;
    codestream << "__m128d c_1_2;" << std::endl;
    codestream << "__m128d c_2_0;" << std::endl;
    codestream << "__m128d c_2_1;" << std::endl;
    codestream << "__m128d c_2_2;" << std::endl;
    codestream << "__m128d b_0;" << std::endl;
    codestream << "__m128d b_1;" << std::endl;
    codestream << "__m128d b_2;" << std::endl;
    codestream << "__m128d a_0;" << std::endl;
    codestream << "__m128d a_1;" << std::endl;
    codestream << "__m128d a_2;" << std::endl;
    codestream << "#endif" << std::endl << std::endl;
    codestream << "#if defined(__SSE3__) && defined(__AVX__)" << std::endl;
    codestream << "__m256d c_0_0;" << std::endl;
    codestream << "__m256d c_0_1;" << std::endl;
    codestream << "__m256d c_0_2;" << std::endl;
    codestream << "__m256d c_1_0;" << std::endl;
    codestream << "__m256d c_1_1;" << std::endl;
    codestream << "__m256d c_1_2;" << std::endl;
    codestream << "__m256d c_2_0;" << std::endl;
    codestream << "__m256d c_2_1;" << std::endl;
    codestream << "__m256d c_2_2;" << std::endl;
    codestream << "__m256d b_0;" << std::endl;
    codestream << "__m256d b_1;" << std::endl;
    codestream << "__m256d b_2;" << std::endl;
    codestream << "__m256d a_0;" << std::endl;
    codestream << "__m256d a_1;" << std::endl;
    codestream << "__m256d a_2;" << std::endl;
    codestream << "#endif" << std::endl << std::endl;

    /////////////////////////
    /////////////////////////
    // generating SSE code //
    /////////////////////////
    /////////////////////////

    // calculate the maximum number of row
    // we can process with max. blocking
    int mSix = (M / 6) * 6;

    if (lda % 2 == 0)
      alignA = true;

    if (ldb % 2 == 0)
      alignB = true;

    if (ldc % 2 == 0)
      alignC = true;

    codestream << "#if defined(__SSE3__) && !defined(__AVX__)" << std::endl;
    codestream << "double* c0 = C;" << std::endl;
    codestream << "double* c1 = C+" << ldc << ";" << std::endl;
    codestream << "double* c2 = C+" << 2 * ldc << ";" << std::endl;
    codestream << "for(int n = 0; n < " << N << "; n+=3)" << std::endl;
    codestream << "{" << std::endl;
    codestream << "  for(int m = 0; m < " << mSix << "; m+=6)" << std::endl;
    codestream << "  {" << std::endl;
    codestream << "    double* b0 = B+(n*" << ldb << ");" << std::endl;
    codestream << "    double* b1 = B+((n+1)*" << ldb << ");" << std::endl;
    codestream << "    double* b2 = B+((n+2)*" << ldb << ");" << std::endl;
    codestream << "    double* a0 = A+m;" << std::endl;

    if (bAdd_) {
      if (alignC == true) {
        codestream << "    c_0_0 = _mm_load_pd(c0);" << std::endl;
        codestream << "    c_0_1 = _mm_load_pd(c1);" << std::endl;
        codestream << "    c_0_2 = _mm_load_pd(c2);" << std::endl;
        codestream << "    c_1_0 = _mm_load_pd(c0+2);" << std::endl;
        codestream << "    c_1_1 = _mm_load_pd(c1+2);" << std::endl;
        codestream << "    c_1_2 = _mm_load_pd(c2+2);" << std::endl;
        codestream << "    c_2_0 = _mm_load_pd(c0+4);" << std::endl;
        codestream << "    c_2_1 = _mm_load_pd(c1+4);" << std::endl;
        codestream << "    c_2_2 = _mm_load_pd(c2+4);" << std::endl << std::endl;
      } else {
        codestream << "    c_0_0 = _mm_loadu_pd(c0);" << std::endl;
        codestream << "    c_0_1 = _mm_loadu_pd(c1);" << std::endl;
        codestream << "    c_0_2 = _mm_loadu_pd(c2);" << std::endl;
        codestream << "    c_1_0 = _mm_loadu_pd(c0+2);" << std::endl;
        codestream << "    c_1_1 = _mm_loadu_pd(c1+2);" << std::endl;
        codestream << "    c_1_2 = _mm_loadu_pd(c2+2);" << std::endl;
        codestream << "    c_2_0 = _mm_loadu_pd(c0+4);" << std::endl;
        codestream << "    c_2_1 = _mm_loadu_pd(c1+4);" << std::endl;
        codestream << "    c_2_2 = _mm_loadu_pd(c2+4);" << std::endl << std::endl;
      }
    } else {
      codestream << "    c_0_0 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_0_1 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_0_2 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_1_0 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_1_1 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_1_2 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_2_0 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_2_1 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_2_2 = _mm_setzero_pd();" << std::endl << std::endl;
    }

#ifndef FULL_UNROLL
    if (this->bGenerateExitForCK_ == true) {
      codestream << "    for (int k = 0; k < exit_col; k++)" << std::endl;
    } else {
      codestream << "    for (int k = 0; k < " << K << "; k++)" << std::endl;
    }
    codestream << "    {" << std::endl;
#else
    for(int k = 0; k < K; k++) {
      if (this->bGenerateExitForCK_ == true) {
        for (int o = 0; o < this->BasisfunctionsCounter_.size(); o++) {
          if (l == this->BasisfunctionsCounter_[o]) {
            codestream << "if ( __builtin_expect(exit_col == " << k << ", false) ) { goto sse_six_end; }" << std::endl;
          }
        }
      }
#endif
    sse_inner_blocked_kernel_6(codestream, lda, alignA);
#ifndef FULL_UNROLL
    codestream << "    }" << std::endl;
#else
    }
    if (this->bGenerateExitForCK_ == true) {
      codestream << "sse_six_end:" << std::endl;
    }
#endif

    if (alignC == true) {
      codestream << "    _mm_store_pd(c0, c_0_0);" << std::endl;
      codestream << "    _mm_store_pd(c1, c_0_1);" << std::endl;
      codestream << "    _mm_store_pd(c2, c_0_2);" << std::endl;
      codestream << "    _mm_store_pd(c0+2, c_1_0);" << std::endl;
      codestream << "    _mm_store_pd(c1+2, c_1_1);" << std::endl;
      codestream << "    _mm_store_pd(c2+2, c_1_2);" << std::endl;
      codestream << "    _mm_store_pd(c0+4, c_2_0);" << std::endl;
      codestream << "    _mm_store_pd(c1+4, c_2_1);" << std::endl;
      codestream << "    _mm_store_pd(c2+4, c_2_2);" << std::endl;
    } else {
      codestream << "    _mm_storeu_pd(c0, c_0_0);" << std::endl;
      codestream << "    _mm_storeu_pd(c1, c_0_1);" << std::endl;
      codestream << "    _mm_storeu_pd(c2, c_0_2);" << std::endl;
      codestream << "    _mm_storeu_pd(c0+2, c_1_0);" << std::endl;
      codestream << "    _mm_storeu_pd(c1+2, c_1_1);" << std::endl;
      codestream << "    _mm_storeu_pd(c2+2, c_1_2);" << std::endl;
      codestream << "    _mm_storeu_pd(c0+4, c_2_0);" << std::endl;
      codestream << "    _mm_storeu_pd(c1+4, c_2_1);" << std::endl;
      codestream << "    _mm_storeu_pd(c2+4, c_2_2);" << std::endl;
    }

    codestream << "    c0+=6;" << std::endl;
    codestream << "    c1+=6;" << std::endl;
    codestream << "    c2+=6;" << std::endl;
    codestream << "  }" << std::endl;

    int mFour = (M/4)*4;
    codestream << "  for(int m = " << mSix << "; m < " << mFour << "; m+=4)" << std::endl;
    codestream << "  {" << std::endl;
    codestream << "    double* b0 = B+(n*" << ldb << ");" << std::endl;
    codestream << "    double* b1 = B+((n+1)*" << ldb << ");" << std::endl;
    codestream << "    double* b2 = B+((n+2)*" << ldb << ");" << std::endl;
    codestream << "    double* a0 = A+m;" << std::endl;

    if (bAdd_) {
      if (alignC == true) {
        codestream << "    c_0_0 = _mm_load_pd(c0);" << std::endl;
        codestream << "    c_0_1 = _mm_load_pd(c1);" << std::endl;
        codestream << "    c_0_2 = _mm_load_pd(c2);" << std::endl;
        codestream << "    c_1_0 = _mm_load_pd(c0+2);" << std::endl;
        codestream << "    c_1_1 = _mm_load_pd(c1+2);" << std::endl;
        codestream << "    c_1_2 = _mm_load_pd(c2+2);" << std::endl;
      } else {
        codestream << "    c_0_0 = _mm_loadu_pd(c0);" << std::endl;
        codestream << "    c_0_1 = _mm_loadu_pd(c1);" << std::endl;
        codestream << "    c_0_2 = _mm_loadu_pd(c2);" << std::endl;
        codestream << "    c_1_0 = _mm_loadu_pd(c0+2);" << std::endl;
        codestream << "    c_1_1 = _mm_loadu_pd(c1+2);" << std::endl;
        codestream << "    c_1_2 = _mm_loadu_pd(c2+2);" << std::endl;
      }
    } else {
      codestream << "    c_0_0 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_0_1 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_0_2 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_1_0 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_1_1 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_1_2 = _mm_setzero_pd();" << std::endl;
    }

#ifndef FULL_UNROLL
    if (this->bGenerateExitForCK_ == true) {
      codestream << "    for (int k = 0; k < exit_col; k++)" << std::endl;
    } else {
      codestream << "    for (int k = 0; k < " << K << "; k++)" << std::endl;
    }
    codestream << "    {" << std::endl;
#else
    for(int k = 0; k < K; k++) {
      if (this->bGenerateExitForCK_ == true) {
        for (int o = 0; o < this->BasisfunctionsCounter_.size(); o++) {
          if (l == this->BasisfunctionsCounter_[o]) {
            codestream << "if ( __builtin_expect(exit_col == " << k << ", false) ) { goto sse_four_end; }" << std::endl;
          }
        }
      }
#endif
    sse_inner_blocked_kernel_4(codestream, lda, alignA);
#ifndef FULL_UNROLL
    codestream << "    }" << std::endl;
#else
    }
    if (this->bGenerateExitForCK_ == true) {
      codestream << "sse_four_end:" << std::endl;
    }
#endif

    if (alignC == true) {
      codestream << "    _mm_store_pd(c0, c_0_0);" << std::endl;
      codestream << "    _mm_store_pd(c1, c_0_1);" << std::endl;
      codestream << "    _mm_store_pd(c2, c_0_2);" << std::endl;
      codestream << "    _mm_store_pd(c0+2, c_1_0);" << std::endl;
      codestream << "    _mm_store_pd(c1+2, c_1_1);" << std::endl;
      codestream << "    _mm_store_pd(c2+2, c_1_2);" << std::endl;
    } else {
      codestream << "    _mm_storeu_pd(c0, c_0_0);" << std::endl;
      codestream << "    _mm_storeu_pd(c1, c_0_1);" << std::endl;
      codestream << "    _mm_storeu_pd(c2, c_0_2);" << std::endl;
      codestream << "    _mm_storeu_pd(c0+2, c_1_0);" << std::endl;
      codestream << "    _mm_storeu_pd(c1+2, c_1_1);" << std::endl;
      codestream << "    _mm_storeu_pd(c2+2, c_1_2);" << std::endl;
    }

    codestream << "    c0+=4;" << std::endl;
    codestream << "    c1+=4;" << std::endl;
    codestream << "    c2+=4;" << std::endl;
    codestream << "  }" << std::endl;

    int mTwo = (M/2)*2;
    codestream << "  for(int m = " << mFour << "; m < " << mTwo << "; m+=2)" << std::endl;
    codestream << "  {" << std::endl;
    codestream << "    double* b0 = B+(n*" << ldb << ");" << std::endl;
    codestream << "    double* b1 = B+((n+1)*" << ldb << ");" << std::endl;
    codestream << "    double* b2 = B+((n+2)*" << ldb << ");" << std::endl;
    codestream << "    double* a0 = A+m;" << std::endl;

    if (bAdd_) {
      if (alignC == true) {
        codestream << "    c_0_0 = _mm_load_pd(c0);" << std::endl;
        codestream << "    c_0_1 = _mm_load_pd(c1);" << std::endl;
        codestream << "    c_0_2 = _mm_load_pd(c2);" << std::endl;
      } else {
        codestream << "    c_0_0 = _mm_loadu_pd(c0);" << std::endl;
        codestream << "    c_0_1 = _mm_loadu_pd(c1);" << std::endl;
        codestream << "    c_0_2 = _mm_loadu_pd(c2);" << std::endl;
      }
    } else {
      codestream << "    c_0_0 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_0_1 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_0_2 = _mm_setzero_pd();" << std::endl;
    }

#ifndef FULL_UNROLL
    if (this->bGenerateExitForCK_ == true) {
      codestream << "    for (int k = 0; k < exit_col; k++)" << std::endl;
    } else {
      codestream << "    for (int k = 0; k < " << K << "; k++)" << std::endl;
    }
    codestream << "    {" << std::endl;
#else
    for(int k = 0; k < K; k++) {
      if (this->bGenerateExitForCK_ == true) {
        for (int o = 0; o < this->BasisfunctionsCounter_.size(); o++) {
          if (l == this->BasisfunctionsCounter_[o]) {
            codestream << "if ( __builtin_expect(exit_col == " << k << ", false) ) { goto sse_two_end; }" << std::endl;
          }
        }
      }
#endif
    sse_inner_blocked_kernel_2(codestream, lda, alignA);
#ifndef FULL_UNROLL
    codestream << "    }" << std::endl;
#else
    }
    if (this->bGenerateExitForCK_ == true) {
      codestream << "sse_two_end:" << std::endl;
    }
#endif

    if (alignC == true) {
      codestream << "    _mm_store_pd(c0, c_0_0);" << std::endl;
      codestream << "    _mm_store_pd(c1, c_0_1);" << std::endl;
      codestream << "    _mm_store_pd(c2, c_0_2);" << std::endl;
    } else {
      codestream << "    _mm_storeu_pd(c0, c_0_0);" << std::endl;
      codestream << "    _mm_storeu_pd(c1, c_0_1);" << std::endl;
      codestream << "    _mm_storeu_pd(c2, c_0_2);" << std::endl;
    }

    codestream << "    c0+=2;" << std::endl;
    codestream << "    c1+=2;" << std::endl;
    codestream << "    c2+=2;" << std::endl;
    codestream << "  }" << std::endl;

    codestream << "  for(int m = " << mTwo << "; m < " << M << "; m+=2)" << std::endl;
    codestream << "  {" << std::endl;
    codestream << "    double* b0 = B+(n*" << ldb << ");" << std::endl;
    codestream << "    double* b1 = B+((n+1)*" << ldb << ");" << std::endl;
    codestream << "    double* b2 = B+((n+2)*" << ldb << ");" << std::endl;
    codestream << "    double* a0 = A+m;" << std::endl;

    if (bAdd_) {
      codestream << "    c_0_0 = _mm_load_sd(c0);" << std::endl;
      codestream << "    c_0_1 = _mm_load_sd(c1);" << std::endl;
      codestream << "    c_0_2 = _mm_load_sd(c2);" << std::endl;
    } else {
      codestream << "    c_0_0 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_0_1 = _mm_setzero_pd();" << std::endl;
      codestream << "    c_0_2 = _mm_setzero_pd();" << std::endl;
    }

#ifndef FULL_UNROLL
    if (this->bGenerateExitForCK_ == true) {
      codestream << "    for (int k = 0; k < exit_col; k++)" << std::endl;
    } else {
      codestream << "    for (int k = 0; k < " << K << "; k++)" << std::endl;
    }
    codestream << "    {" << std::endl;
#else
    for(int k = 0; k < K; k++) {
      if (this->bGenerateExitForCK_ == true) {
        for (int o = 0; o < this->BasisfunctionsCounter_.size(); o++) {
          if (l == this->BasisfunctionsCounter_[o]) {
            codestream << "if ( __builtin_expect(exit_col == " << k << ", false) ) { goto sse_one_end; }" << std::endl;
          }
        }
      }
#endif
    sse_inner_blocked_kernel_2(codestream, lda, alignA);
#ifndef FULL_UNROLL
    codestream << "    }" << std::endl;
#else
    }
    if (this->bGenerateExitForCK_ == true) {
      codestream << "sse_one_end:" << std::endl;
    }
#endif

    codestream << "    _mm_store_sd(c0, c_0_0);" << std::endl;
    codestream << "    _mm_store_sd(c1, c_0_1);" << std::endl;
    codestream << "    _mm_store_sd(c2, c_0_2);" << std::endl;

    codestream << "    c0+=1;" << std::endl;
    codestream << "    c1+=1;" << std::endl;
    codestream << "    c2+=1;" << std::endl;
    codestream << "  }" << std::endl;

    codestream << "  c0+=" << (2 * ldc) << ";" << std::endl;
    codestream << "  c1+=" << (2 * ldc) << ";" << std::endl;
    codestream << "  c2+=" << (2 * ldc) << ";" << std::endl;
    codestream << "}" << std::endl << std::endl;

    codestream << "#endif" << std::endl << std::endl;

    /////////////////////////
    // generating AVX code //
    /////////////////////////

    // calculate the maximum number of row
    // we can process with max. blocking
    // SSE case
    maxM = (M / 12) * 12;
    remainder = M - maxM;

    if (lda % 4 == 0)
      alignA = true;

    if (ldb % 4 == 0)
      alignB = true;

    if (ldc % 4 == 0)
      alignC = true;

    codestream << "#if defined(__SSE3__) && defined(__AVX__)" << std::endl;
    codestream << "double* c0 = C;" << std::endl;
    codestream << "double* c1 = C+" << ldc << ";" << std::endl;
    codestream << "double* c2 = C+" << 2 * ldc << ";" << std::endl;
    codestream << "for(int n = 0; n < " << N << "; n+=3)" << std::endl;
    codestream << "{" << std::endl;
    codestream << "  for(int m = 0; m < " << maxM << "; m+=12)" << std::endl;
    codestream << "  {" << std::endl;
    codestream << "    double* b0 = B+(n*" << ldb << ");" << std::endl;
    codestream << "    double* b1 = B+((n+1)*" << ldb << ");" << std::endl;
    codestream << "    double* b2 = B+((n+2)*" << ldb << ");" << std::endl;
    codestream << "    double* a0 = A+m;" << std::endl;

    if (bAdd_) {
      if (alignC == true) {
        codestream << "    c_0_0 = _mm256_load_pd(c0);" << std::endl;
        codestream << "    c_0_1 = _mm256_load_pd(c1);" << std::endl;
        codestream << "    c_0_2 = _mm256_load_pd(c2);" << std::endl;
        codestream << "    c_1_0 = _mm256_load_pd(c0+4);" << std::endl;
        codestream << "    c_1_1 = _mm256_load_pd(c1+4);" << std::endl;
        codestream << "    c_1_2 = _mm256_load_pd(c2+4);" << std::endl;
        codestream << "    c_2_0 = _mm256_load_pd(c0+8);" << std::endl;
        codestream << "    c_2_1 = _mm256_load_pd(c1+8);" << std::endl;
        codestream << "    c_2_2 = _mm256_load_pd(c2+8);" << std::endl << std::endl;
      } else {
        codestream << "    c_0_0 = _mm256_loadu_pd(c0);" << std::endl;
        codestream << "    c_0_1 = _mm256_loadu_pd(c1);" << std::endl;
        codestream << "    c_0_2 = _mm256_loadu_pd(c2);" << std::endl;
        codestream << "    c_1_0 = _mm256_loadu_pd(c0+4);" << std::endl;
        codestream << "    c_1_1 = _mm256_loadu_pd(c1+4);" << std::endl;
        codestream << "    c_1_2 = _mm256_loadu_pd(c2+4);" << std::endl;
        codestream << "    c_2_0 = _mm256_loadu_pd(c0+8);" << std::endl;
        codestream << "    c_2_1 = _mm256_loadu_pd(c1+8);" << std::endl;
        codestream << "    c_2_2 = _mm256_loadu_pd(c2+8);" << std::endl << std::endl;
      }
    } else {
      codestream << "    c_0_0 = _mm256_setzero_pd();" << std::endl;
      codestream << "    c_0_1 = _mm256_setzero_pd();" << std::endl;
      codestream << "    c_0_2 = _mm256_setzero_pd();" << std::endl;
      codestream << "    c_1_0 = _mm256_setzero_pd();" << std::endl;
      codestream << "    c_1_1 = _mm256_setzero_pd();" << std::endl;
      codestream << "    c_1_2 = _mm256_setzero_pd();" << std::endl;
      codestream << "    c_2_0 = _mm256_setzero_pd();" << std::endl;
      codestream << "    c_2_1 = _mm256_setzero_pd();" << std::endl;
      codestream << "    c_2_2 = _mm256_setzero_pd();" << std::endl << std::endl;
    }

    //codestream << "    #pragma unroll(" << unroll_factor << ")" << std::endl;
    codestream << "    for (int k = 0; k < " << K << "; k++)" << std::endl;
    codestream << "    {" << std::endl;

    if (this->bGenerateExitForCK_ == true) {
      codestream << "    if ( __builtin_expect(exit_col == k, false) ) { break; }" << std::endl;
    }

    avx_inner_blocked_kernel(codestream, lda, alignA);

    codestream << "    }" << std::endl;

    if (alignC == true) {
      codestream << "    _mm256_store_pd(c0, c_0_0);" << std::endl;
      codestream << "    _mm256_store_pd(c1, c_0_1);" << std::endl;
      codestream << "    _mm256_store_pd(c2, c_0_2);" << std::endl;
      codestream << "    _mm256_store_pd(c0+4, c_1_0);" << std::endl;
      codestream << "    _mm256_store_pd(c1+4, c_1_1);" << std::endl;
      codestream << "    _mm256_store_pd(c2+4, c_1_2);" << std::endl;
      codestream << "    _mm256_store_pd(c0+8, c_2_0);" << std::endl;
      codestream << "    _mm256_store_pd(c1+8, c_2_1);" << std::endl;
      codestream << "    _mm256_store_pd(c2+8, c_2_2);" << std::endl;
    } else {
      codestream << "    _mm256_storeu_pd(c0, c_0_0);" << std::endl;
      codestream << "    _mm256_storeu_pd(c1, c_0_1);" << std::endl;
      codestream << "    _mm256_storeu_pd(c2, c_0_2);" << std::endl;
      codestream << "    _mm256_storeu_pd(c0+4, c_1_0);" << std::endl;
      codestream << "    _mm256_storeu_pd(c1+4, c_1_1);" << std::endl;
      codestream << "    _mm256_storeu_pd(c2+4, c_1_2);" << std::endl;
      codestream << "    _mm256_storeu_pd(c0+8, c_2_0);" << std::endl;
      codestream << "    _mm256_storeu_pd(c1+8, c_2_1);" << std::endl;
      codestream << "    _mm256_storeu_pd(c2+8, c_2_2);" << std::endl;
    }

    codestream << "    c0+=12;" << std::endl;
    codestream << "    c1+=12;" << std::endl;
    codestream << "    c2+=12;" << std::endl;
    codestream << "  }" << std::endl;
    codestream << "  c0+=" << (2 * ldc) + (ldc - maxM) << ";" << std::endl;
    codestream << "  c1+=" << (2 * ldc) + (ldc - maxM) << ";" << std::endl;
    codestream << "  c2+=" << (2 * ldc) + (ldc - maxM) << ";" << std::endl;
    codestream << "}" << std::endl << std::endl;

    // do the remainder
    int A_M = maxM;
    //codestream << "int m = " << maxM << ";" << std::endl;
    codestream << "double* c0_base = C+" << maxM << ";" << std::endl;
    codestream << "double* c1_base = C+" << maxM + ldc << ";" << std::endl;
    codestream << "double* c2_base = C+" << maxM + (2 * ldc) << ";" << std::endl;
    codestream << "for(int n = 0; n < " << N << "; n+=3)" << std::endl;
    codestream << "{" << std::endl;

    codestream << "  c0 = c0_base;" << std::endl;
    codestream << "  c1 = c1_base;" << std::endl;
    codestream << "  c2 = c2_base;" << std::endl << std::endl;

    codestream << "  double* b0;" << std::endl;
    codestream << "  double* b1;" << std::endl;
    codestream << "  double* b2;" << std::endl;
    codestream << "  double* a0;" << std::endl;

    if (remainder >= 8) {
      remainder -= 8;
      codestream << "  b0 = B+(n*" << ldb << ");" << std::endl;
      codestream << "  b1 = B+((n+1)*" << ldb << ");" << std::endl;
      codestream << "  b2 = B+((n+2)*" << ldb << ");" << std::endl;
      codestream << "  a0 = A+" << A_M << ";" << std::endl;
      A_M += 8;

      if (bAdd_) {
        if (alignC == true) {
          codestream << "  c_0_0 = _mm256_load_pd(c0);" << std::endl;
          codestream << "  c_0_1 = _mm256_load_pd(c1);" << std::endl;
          codestream << "  c_0_2 = _mm256_load_pd(c2);" << std::endl;
          codestream << "  c_1_0 = _mm256_load_pd(c0+4);" << std::endl;
          codestream << "  c_1_1 = _mm256_load_pd(c1+4);" << std::endl;
          codestream << "  c_1_2 = _mm256_load_pd(c2+4);" << std::endl << std::endl;
        } else {
          codestream << "  c_0_0 = _mm256_loadu_pd(c0);" << std::endl;
          codestream << "  c_0_1 = _mm256_loadu_pd(c1);" << std::endl;
          codestream << "  c_0_2 = _mm256_loadu_pd(c2);" << std::endl;
          codestream << "  c_1_0 = _mm256_loadu_pd(c0+4);" << std::endl;
          codestream << "  c_1_1 = _mm256_loadu_pd(c1+4);" << std::endl;
          codestream << "  c_1_2 = _mm256_loadu_pd(c2+4);" << std::endl << std::endl;
        }
      } else {
        codestream << "  c_0_0 = _mm256_setzero_pd();" << std::endl;
        codestream << "  c_0_1 = _mm256_setzero_pd();" << std::endl;
        codestream << "  c_0_2 = _mm256_setzero_pd();" << std::endl;
        codestream << "  c_1_0 = _mm256_setzero_pd();" << std::endl;
        codestream << "  c_1_1 = _mm256_setzero_pd();" << std::endl;
        codestream << "  c_1_2 = _mm256_setzero_pd();" << std::endl << std::endl;
      }

      //            codestream << "  #pragma unroll(" << unroll_factor << ")" << std::endl;
      codestream << "  for (int k = 0; k < " << K << "; k++)" << std::endl;
      codestream << "  {" << std::endl;

      if (this->bGenerateExitForCK_ == true) {
        codestream << "  if ( __builtin_expect(exit_col == k, false) ) { break; }" << std::endl;
      }

      codestream << "  b_0 = _mm256_broadcast_sd(b0);" << std::endl;
      codestream << "  b_1 = _mm256_broadcast_sd(b1);" << std::endl;
      codestream << "  b_2 = _mm256_broadcast_sd(b2);" << std::endl << std::endl;
      codestream << "  b0++; b1++; b2++;" << std::endl << std::endl;

      if (alignA == true) {
        codestream << "  a_0 = _mm256_load_pd(a0);" << std::endl;
      } else {
        codestream << "  a_0 = _mm256_loadu_pd(a0);" << std::endl;
      }

      codestream << "  a0+=" << 4 << ";" << std::endl;
      codestream << "  c_0_0 = _mm256_add_pd(c_0_0, _mm256_mul_pd(a_0, b_0));" << std::endl;
      codestream << "  c_0_1 = _mm256_add_pd(c_0_1, _mm256_mul_pd(a_0, b_1));" << std::endl;
      codestream << "  c_0_2 = _mm256_add_pd(c_0_2, _mm256_mul_pd(a_0, b_2));" << std::endl << std::endl;

      if (alignA == true) {
        codestream << "  a_1 = _mm256_load_pd(a0);" << std::endl;
      } else {
        codestream << "  a_1 = _mm256_loadu_pd(a0);" << std::endl;
      }

      codestream << "  a0+=" << lda - 4 << ";" << std::endl;
      codestream << "  c_1_0 = _mm256_add_pd(c_1_0, _mm256_mul_pd(a_1, b_0));" << std::endl;
      codestream << "  c_1_1 = _mm256_add_pd(c_1_1, _mm256_mul_pd(a_1, b_1));" << std::endl;
      codestream << "  c_1_2 = _mm256_add_pd(c_1_2, _mm256_mul_pd(a_1, b_2));" << std::endl << std::endl;

      codestream << "  }" << std::endl;

      if (alignC == true) {
        codestream << "  _mm256_store_pd(c0, c_0_0);" << std::endl;
        codestream << "  _mm256_store_pd(c1, c_0_1);" << std::endl;
        codestream << "  _mm256_store_pd(c2, c_0_2);" << std::endl;
        codestream << "  _mm256_store_pd(c0+4, c_1_0);" << std::endl;
        codestream << "  _mm256_store_pd(c1+4, c_1_1);" << std::endl;
        codestream << "  _mm256_store_pd(c2+4, c_1_2);" << std::endl << std::endl;
      } else {
        codestream << "  _mm256_storeu_pd(c0, c_0_0);" << std::endl;
        codestream << "  _mm256_storeu_pd(c1, c_0_1);" << std::endl;
        codestream << "  _mm256_storeu_pd(c2, c_0_2);" << std::endl;
        codestream << "  _mm256_storeu_pd(c0+4, c_1_0);" << std::endl;
        codestream << "  _mm256_storeu_pd(c1+4, c_1_1);" << std::endl;
        codestream << "  _mm256_storeu_pd(c2+4, c_1_2);" << std::endl << std::endl;
      }

      codestream << "  c0+=8;" << std::endl;
      codestream << "  c1+=8;" << std::endl;
      codestream << "  c2+=8;" << std::endl;
    }

    if (remainder >= 4) {
      remainder -= 4;
      codestream << "  b0 = B+(n*" << ldb << ");" << std::endl;
      codestream << "  b1 = B+((n+1)*" << ldb << ");" << std::endl;
      codestream << "  b2 = B+((n+2)*" << ldb << ");" << std::endl;
      codestream << "  a0 = A+" << A_M << ";" << std::endl;
      A_M += 4;

      if (bAdd_) {
        if (alignC == true) {
          codestream << "  c_0_0 = _mm256_load_pd(c0);" << std::endl;
          codestream << "  c_0_1 = _mm256_load_pd(c1);" << std::endl;
          codestream << "  c_0_2 = _mm256_load_pd(c2);" << std::endl;
        } else {
          codestream << "  c_0_0 = _mm256_loadu_pd(c0);" << std::endl;
          codestream << "  c_0_1 = _mm256_loadu_pd(c1);" << std::endl;
          codestream << "  c_0_2 = _mm256_loadu_pd(c2);" << std::endl;
        }
      } else {
        codestream << "  c_0_0 = _mm256_setzero_pd();" << std::endl;
        codestream << "  c_0_1 = _mm256_setzero_pd();" << std::endl;
        codestream << "  c_0_2 = _mm256_setzero_pd();" << std::endl;

      }

      //            codestream << "  #pragma unroll(" << unroll_factor << ")" << std::endl;
      codestream << "  for (int k = 0; k < " << K << "; k++)" << std::endl;
      codestream << "  {" << std::endl;

      if (this->bGenerateExitForCK_ == true) {
        codestream << "  if ( __builtin_expect(exit_col == k, false) ) { break; }" << std::endl;
      }

      codestream << "  b_0 = _mm256_broadcast_sd(b0);" << std::endl;
      codestream << "  b_1 = _mm256_broadcast_sd(b1);" << std::endl;
      codestream << "  b_2 = _mm256_broadcast_sd(b2);" << std::endl << std::endl;
      codestream << "  b0++; b1++; b2++;" << std::endl << std::endl;

      if (alignA == true) {
        codestream << "  a_0 = _mm256_load_pd(a0);" << std::endl;
      } else {
        codestream << "  a_0 = _mm256_loadu_pd(a0);" << std::endl;
      }

      codestream << "  a0+=" << lda << ";" << std::endl;
      codestream << "  c_0_0 = _mm256_add_pd(c_0_0, _mm256_mul_pd(a_0, b_0));" << std::endl;
      codestream << "  c_0_1 = _mm256_add_pd(c_0_1, _mm256_mul_pd(a_0, b_1));" << std::endl;
      codestream << "  c_0_2 = _mm256_add_pd(c_0_2, _mm256_mul_pd(a_0, b_2));" << std::endl << std::endl;

      codestream << "  }" << std::endl;

      if (alignC == true) {
        codestream << "  _mm256_store_pd(c0, c_0_0);" << std::endl;
        codestream << "  _mm256_store_pd(c1, c_0_1);" << std::endl;
        codestream << "  _mm256_store_pd(c2, c_0_2);" << std::endl;
      } else {
        codestream << "  _mm256_storeu_pd(c0, c_0_0);" << std::endl;
        codestream << "  _mm256_storeu_pd(c1, c_0_1);" << std::endl;
        codestream << "  _mm256_storeu_pd(c2, c_0_2);" << std::endl;
      }

      codestream << "  c0+=4;" << std::endl;
      codestream << "  c1+=4;" << std::endl;
      codestream << "  c2+=4;" << std::endl;
    }

    if (remainder >= 2) {
      remainder -= 2;
      codestream << "  b0 = B+(n*" << ldb << ");" << std::endl;
      codestream << "  b1 = B+((n+1)*" << ldb << ");" << std::endl;
      codestream << "  b2 = B+((n+2)*" << ldb << ");" << std::endl;
      codestream << "  a0 = A+" << A_M << ";" << std::endl;
      A_M += 2;

      if (bAdd_) {
        if (alignC == true) {
          codestream << "  c_0_0_128 = _mm_load_pd(c0);" << std::endl;
          codestream << "  c_0_1_128 = _mm_load_pd(c1);" << std::endl;
          codestream << "  c_0_2_128 = _mm_load_pd(c2);" << std::endl << std::endl;
        } else {
          codestream << "  c_0_0_128 = _mm_loadu_pd(c0);" << std::endl;
          codestream << "  c_0_1_128 = _mm_loadu_pd(c1);" << std::endl;
          codestream << "  c_0_2_128 = _mm_loadu_pd(c2);" << std::endl << std::endl;
        }
      } else {
        codestream << "  c_0_0_128 = _mm_setzero_pd();" << std::endl;
        codestream << "  c_0_1_128 = _mm_setzero_pd();" << std::endl;
        codestream << "  c_0_2_128 = _mm_setzero_pd();" << std::endl << std::endl;
      }

      //            codestream << "  #pragma unroll(" << unroll_factor << ")" << std::endl;
      codestream << "  for (int k = 0; k < " << K << "; k++)" << std::endl;
      codestream << "  {" << std::endl;

      if (this->bGenerateExitForCK_ == true) {
        codestream << "  if ( __builtin_expect(exit_col == k, false) ) { break; }" << std::endl;
      }

      codestream << "  b_0_128 = _mm_loaddup_pd(b0);" << std::endl;
      codestream << "  b_1_128 = _mm_loaddup_pd(b1);" << std::endl;
      codestream << "  b_2_128 = _mm_loaddup_pd(b2);" << std::endl << std::endl;
      codestream << "  b0++; b1++; b2++;" << std::endl << std::endl;

      if (alignA == true) {
        codestream << "  a_0_128 = _mm_load_pd(a0);" << std::endl;
      } else {
        codestream << "  a_0_128 = _mm_loadu_pd(a0);" << std::endl;
      }

      codestream << "  a0+=" << lda << ";" << std::endl;
      codestream << "  c_0_0_128 = _mm_add_pd(c_0_0_128, _mm_mul_pd(a_0_128, b_0_128));" << std::endl;
      codestream << "  c_0_1_128 = _mm_add_pd(c_0_1_128, _mm_mul_pd(a_0_128, b_1_128));" << std::endl;
      codestream << "  c_0_2_128 = _mm_add_pd(c_0_2_128, _mm_mul_pd(a_0_128, b_2_128));" << std::endl << std::endl;

      codestream << "  }" << std::endl;

      if (alignC == true) {
        codestream << "  _mm_store_pd(c0, c_0_0_128);" << std::endl;
        codestream << "  _mm_store_pd(c1, c_0_1_128);" << std::endl;
        codestream << "  _mm_store_pd(c2, c_0_2_128);" << std::endl << std::endl;
      } else {
        codestream << "  _mm_storeu_pd(c0, c_0_0_128);" << std::endl;
        codestream << "  _mm_storeu_pd(c1, c_0_1_128);" << std::endl;
        codestream << "  _mm_storeu_pd(c2, c_0_2_128);" << std::endl << std::endl;
      }

      codestream << "  c0+=2;" << std::endl;
      codestream << "  c1+=2;" << std::endl;
      codestream << "  c2+=2;" << std::endl;
    }

    if (remainder == 1) {
      codestream << "  b0 = B+(n*" << ldb << ");" << std::endl;
      codestream << "  b1 = B+((n+1)*" << ldb << ");" << std::endl;
      codestream << "  b2 = B+((n+2)*" << ldb << ");" << std::endl;
      codestream << "  a0 = A+" << A_M << ";" << std::endl;
      A_M++;

      if (bAdd_) {
        codestream << "  c_0_0_128 = _mm_load_sd(c0);" << std::endl;
        codestream << "  c_0_1_128 = _mm_load_sd(c1);" << std::endl;
        codestream << "  c_0_2_128 = _mm_load_sd(c2);" << std::endl << std::endl;
      } else {
        codestream << "  c_0_0_128 = _mm_setzero_pd();" << std::endl;
        codestream << "  c_0_1_128 = _mm_setzero_pd();" << std::endl;
        codestream << "  c_0_2_128 = _mm_setzero_pd();" << std::endl << std::endl;
      }

      //            codestream << "  #pragma unroll(" << unroll_factor << ")" << std::endl;
      codestream << "  for (int k = 0; k < " << K << "; k++)" << std::endl;
      codestream << "  {" << std::endl;

      if (this->bGenerateExitForCK_ == true) {
        codestream << "  if ( __builtin_expect(exit_col == k, false) ) { break; }" << std::endl;
      }

      codestream << "  b_0_128 = _mm_load_sd(b0);" << std::endl;
      codestream << "  b_1_128 = _mm_load_sd(b1);" << std::endl;
      codestream << "  b_2_128 = _mm_load_sd(b2);" << std::endl << std::endl;
      codestream << "  b0++; b1++; b2++;" << std::endl << std::endl;
      codestream << "  a_0_128 = _mm_load_sd(a0);" << std::endl;
      codestream << "  a0+=" << lda << ";" << std::endl;
      codestream << "  c_0_0_128 = _mm_add_sd(c_0_0_128, _mm_mul_sd(a_0_128, b_0_128));" << std::endl;
      codestream << "  c_0_1_128 = _mm_add_sd(c_0_1_128, _mm_mul_sd(a_0_128, b_1_128));" << std::endl;
      codestream << "  c_0_2_128 = _mm_add_sd(c_0_2_128, _mm_mul_sd(a_0_128, b_2_128));" << std::endl << std::endl;

      codestream << "  }" << std::endl;
      codestream << "  _mm_store_sd(c0, c_0_0_128);" << std::endl;
      codestream << "  _mm_store_sd(c1, c_0_1_128);" << std::endl;
      codestream << "  _mm_store_sd(c2, c_0_2_128);" << std::endl << std::endl;

      codestream << "  c0+=1;" << std::endl;
      codestream << "  c1+=1;" << std::endl;
      codestream << "  c2+=1;" << std::endl;
    }

    codestream << "  c0_base+=" << (3 * ldc) << ";" << std::endl;
    codestream << "  c1_base+=" << (3 * ldc) << ";" << std::endl;
    codestream << "  c2_base+=" << (3 * ldc) << ";" << std::endl;

    codestream << "}" << std::endl;
    codestream << "#endif" << std::endl << std::endl;

    /////////////////////////
    // generating MIC code //
    /////////////////////////

    if ( (M == 56) ) {
      //if ( (M == 56) && (K == 56) ) {
      codestream << "#if defined(__MIC__)" << std::endl;
      codestream << "__m512d c_0_0;" << std::endl;
      codestream << "__m512d c_1_0;" << std::endl;
      codestream << "__m512d c_2_0;" << std::endl;
      codestream << "__m512d c_3_0;" << std::endl;
      codestream << "__m512d c_4_0;" << std::endl;
      codestream << "__m512d c_5_0;" << std::endl;
      codestream << "__m512d c_6_0;" << std::endl << std::endl;

      codestream << "__m512d c_0_1;" << std::endl;
      codestream << "__m512d c_1_1;" << std::endl;
      codestream << "__m512d c_2_1;" << std::endl;
      codestream << "__m512d c_3_1;" << std::endl;
      codestream << "__m512d c_4_1;" << std::endl;
      codestream << "__m512d c_5_1;" << std::endl;
      codestream << "__m512d c_6_1;" << std::endl << std::endl;

      codestream << "__m512d c_0_2;" << std::endl;
      codestream << "__m512d c_1_2;" << std::endl;
      codestream << "__m512d c_2_2;" << std::endl;
      codestream << "__m512d c_3_2;" << std::endl;
      codestream << "__m512d c_4_2;" << std::endl;
      codestream << "__m512d c_5_2;" << std::endl;
      codestream << "__m512d c_6_2;" << std::endl << std::endl;

      codestream << "__m512d b_0;" << std::endl;
      codestream << "__m512d b_1;" << std::endl;
      codestream << "__m512d b_2;" << std::endl << std::endl;

      codestream << "__m512d a_0;" << std::endl;
      codestream << "__m512d a_1;" << std::endl;
      codestream << "__m512d a_2;" << std::endl;
      codestream << "__m512d a_3;" << std::endl;
      codestream << "__m512d a_4;" << std::endl;
      codestream << "__m512d a_5;" << std::endl;
      codestream << "__m512d a_6;" << std::endl << std::endl;

      codestream << "double* c0 = C;" << std::endl;
      codestream << "double* c1 = C + " << ldc << ";" << std::endl;
      codestream << "double* c2 = C + " << 2 * ldc << ";" << std::endl;
      codestream << "#pragma prefetch c0,c1,c2" << std::endl;
      codestream << "for(int n = 0; n < " << N << "; n+=3)" << std::endl;
      codestream << "{" << std::endl;
      codestream << "  double* b0 = B+(n*" << ldb << ");" << std::endl;
      codestream << "  double* b1 = B+((n+1)*" << ldb << ");" << std::endl;
      codestream << "  double* b2 = B+((n+2)*" << ldb << ");" << std::endl;
      codestream << "  double* a0 = A;" << std::endl << std::endl;

      if (bAdd_) {
        codestream << "  c_0_0 = _mm512_load_pd(c0);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_1_0 = _mm512_load_pd(c0+8);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+8+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_2_0 = _mm512_load_pd(c0+16);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+16+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_3_0 = _mm512_load_pd(c0+24);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+24+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_4_0 = _mm512_load_pd(c0+32);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+32+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_5_0 = _mm512_load_pd(c0+40);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+40+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_6_0 = _mm512_load_pd(c0+48);" << std::endl << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+48+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;

        codestream << "  c_0_1 = _mm512_load_pd(c1);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_1_1 = _mm512_load_pd(c1+8);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+8+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_2_1 = _mm512_load_pd(c1+16);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+16+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_3_1 = _mm512_load_pd(c1+24);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+24+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_4_1 = _mm512_load_pd(c1+32);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+32+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_5_1 = _mm512_load_pd(c1+40);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+40+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_6_1 = _mm512_load_pd(c1+48);" << std::endl << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+48+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;

        codestream << "  c_0_2 = _mm512_load_pd(c2);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_1_2 = _mm512_load_pd(c2+8);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+8+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_2_2 = _mm512_load_pd(c2+16);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+16+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_3_2 = _mm512_load_pd(c2+24);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+24+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_4_2 = _mm512_load_pd(c2+32);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+32+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_5_2 = _mm512_load_pd(c2+40);" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+40+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_6_2 = _mm512_load_pd(c2+48);" << std::endl << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+48+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
      } else {
        codestream << "  c_0_0 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_1_0 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+8+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_2_0 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+16+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_3_0 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+24+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_4_0 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+32+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_5_0 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+40+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_6_0 = _mm512_setzero_pd();" << std::endl << std::endl;
        //codestream << "  _mm_prefetch((const char*)c0+48+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;

        codestream << "  c_0_1 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_1_1 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+8+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_2_1 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+16+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_3_1 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+24+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_4_1 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+32+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_5_1 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+40+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_6_1 = _mm512_setzero_pd();" << std::endl << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+48+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;

        codestream << "  c_0_2 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_1_2 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+8+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_2_2 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+16+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_3_2 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+24+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_4_2 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+32+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_5_2 = _mm512_setzero_pd();" << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+40+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
        codestream << "  c_6_2 = _mm512_setzero_pd();" << std::endl << std::endl;
        //codestream << "  _mm_prefetch((const char*)c1+48+" << 3*ldc << ", _MM_HINT_T1);" << std::endl;
      }

      /// loop
      codestream << "  #pragma prefetch b0,b1,b2,a0" << std::endl;
      codestream << "  for(int k = 0; k < " << K << "; k++)" << std::endl;
      codestream << "  {" << std::endl;

      if (this->bGenerateExitForCK_ == true) {
        codestream << "    if ( __builtin_expect(exit_col == k, false) ) { break; }" << std::endl;
      }

      //codestream << "    _mm_prefetch((const char*)b0+" << 3*ldb << ", _MM_HINT_T1);" << std::endl;
      //codestream << "    _mm_prefetch((const char*)b1+" << 3*ldb << ", _MM_HINT_T1);" << std::endl;
      //codestream << "    _mm_prefetch((const char*)b2+" << 3*ldb << ", _MM_HINT_T1);" << std::endl;
      mic_inner_blocked_kernel_for_56(codestream, lda);

      codestream << "  }" << std::endl << std::endl;

      codestream << "  _mm512_store_pd(c0, c_0_0);" << std::endl;
      codestream << "  _mm512_store_pd(c0+8, c_1_0);" << std::endl;
      codestream << "  _mm512_store_pd(c0+16, c_2_0);" << std::endl;
      codestream << "  _mm512_store_pd(c0+24, c_3_0);" << std::endl;
      codestream << "  _mm512_store_pd(c0+32, c_4_0);" << std::endl;
      codestream << "  _mm512_store_pd(c0+40, c_5_0);" << std::endl;
      codestream << "  _mm512_store_pd(c0+48, c_6_0);" << std::endl << std::endl;

      codestream << "  _mm512_store_pd(c1, c_0_1);" << std::endl;
      codestream << "  _mm512_store_pd(c1+8, c_1_1);" << std::endl;
      codestream << "  _mm512_store_pd(c1+16, c_2_1);" << std::endl;
      codestream << "  _mm512_store_pd(c1+24, c_3_1);" << std::endl;
      codestream << "  _mm512_store_pd(c1+32, c_4_1);" << std::endl;
      codestream << "  _mm512_store_pd(c1+40, c_5_1);" << std::endl;
      codestream << "  _mm512_store_pd(c1+48, c_6_1);" << std::endl << std::endl;

      codestream << "  _mm512_store_pd(c2, c_0_2);" << std::endl;
      codestream << "  _mm512_store_pd(c2+8, c_1_2);" << std::endl;
      codestream << "  _mm512_store_pd(c2+16, c_2_2);" << std::endl;
      codestream << "  _mm512_store_pd(c2+24, c_3_2);" << std::endl;
      codestream << "  _mm512_store_pd(c2+32, c_4_2);" << std::endl;
      codestream << "  _mm512_store_pd(c2+40, c_5_2);" << std::endl;
      codestream << "  _mm512_store_pd(c2+48, c_6_2);" << std::endl << std::endl;

      codestream << "  c0 += " << 3 * ldc << ";" << std::endl;
      codestream << "  c1 += " << 3 * ldc << ";" << std::endl;
      codestream << "  c2 += " << 3 * ldc << ";" << std::endl;

      codestream << "}" << std::endl;
      codestream << "#endif" << std::endl << std::endl;
    }

    // generate fallback c code
    if ( (M == 56) ) {
      //if ( (M == 56) && (K == 56) ) {
      codestream << "#if !defined(__SSE3__) && !defined(__AVX__) && !defined(__MIC__)" << std::endl;
    } else {
      codestream << "#if !defined(__SSE3__) && !defined(__AVX__)" << std::endl;
    }

    codestream << "for (int n = 0; n < " << N << "; n++)" << std::endl;
    codestream << "{" << std::endl;
    codestream << "  for (int k = 0; k < " << K << "; k++)" << std::endl;
    codestream << "  {" << std::endl;

    if (this->bGenerateExitForCK_ == true) {
      codestream << "    if ( __builtin_expect(exit_col == k, false) ) { break; }" << std::endl;
    }

    codestream << "    for(int m = 0; m < " << M << "; m++)" << std::endl;
    codestream << "    {" << std::endl;
    codestream << "      C[(n*" << ldc << ")+m] += A[(k*" << lda << ")+m] * B[(n*" << ldb << ")+k];" << std::endl;
    codestream << "    }" << std::endl;
    codestream << "  }" << std::endl;
    codestream << "}" << std::endl;
    codestream << "#endif" << std::endl << std::endl;

    codestream << "#ifndef NDEBUG" << std::endl;

    if (this->bGenerateExitForCK_ == true) {
      codestream << "num_flops += " << 2 * N* M << "*exit_col;" << std::endl;
    } else {
      codestream << "num_flops += " << 2 * N* M* K << ";" << std::endl;
    }

    codestream << "#endif" << std::endl << std::endl;

    return codestream.str();
  }

}

