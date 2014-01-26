/******************************************************************************
* Copyright (C) 2012-2014 Technische Universitaet Muenchen                    *
*                                                                             *
* This file is part of the SeisSolGen project. For conditions of              *
* distribution and use, please see the copyright notice                       *
* in its root directory.                                                      *
******************************************************************************/
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <immintrin.h>
#include <sys/time.h>

#include "gen_matmul_dense.hpp"

#ifndef ORDER_NUMBER
#define ORDER_NUMBER 20
#endif

#define REPS 100000

inline double sec(struct timeval start, struct timeval end) {
  return ((double)(((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)))) / 1.0e6;
}

int main(int argc, char* argv[]) {
  // allocate
  double* a = (double*)_mm_malloc(ORDER_NUMBER * ORDER_NUMBER * sizeof(double), 64);
  double* b = (double*)_mm_malloc(ORDER_NUMBER * 9 * sizeof(double), 64);
  double* c = (double*)_mm_malloc(ORDER_NUMBER * 9 * sizeof(double), 64);
  double* c_gold = (double*)_mm_malloc(ORDER_NUMBER * 9 * sizeof(double), 64);

  // touch
  for (int i = 0; i < ORDER_NUMBER; i++) {
    for (int j = 0; j < ORDER_NUMBER; j++) {
      a[(j * ORDER_NUMBER) + i] = (double)(i + (j * ORDER_NUMBER));

      if (j < 9) {
        b[(j * 9) + i] = (double)(i + (j * 9));
        c[(j * 9) + i] = 0.0;
        c_gold[(j * 9) + i] = 0.0;
      }
    }
  }

  // C routine
  struct timeval start, end;
  gettimeofday(&start, NULL);

  for (int t = 0; t < REPS; t++) {
    for (int n = 0; n < 9; n++) {
      for (int k = 0; k < ORDER_NUMBER; k++) {
#pragma simd

        for (int m = 0; m < ORDER_NUMBER; m++) {
          c_gold[(n * ORDER_NUMBER) + m] += a[(k * ORDER_NUMBER) + m] * b[(n * ORDER_NUMBER) + k];
        }
      }
    }
  }

  gettimeofday(&end, NULL);
  double total = sec(start, end);

  std::cout << total << "s for C" << std::endl;
  std::cout << ((double)(REPS * ORDER_NUMBER * ORDER_NUMBER) * 9.0 * 2.0) / (total * 1.0e9) << " GFLOPS for C" << std::endl;

  gettimeofday(&start, NULL);

  for (int t = 0; t < REPS; t++) {
    dense_test(a, b, c, 84);
#if 0
    __m128d a_0, a_1, a_2, a_3;
    __m128d b_0, b_1, b_2;
    __m128d c_0_0, c_0_1, c_0_2;
    __m128d c_1_0, c_1_1, c_1_2;
    __m128d c_2_0, c_2_1, c_2_2;
    __m128d c_3_0, c_3_1, c_3_2;

    int saveend;
    int newrange;
    int end = (ORDER_NUMBER / 6) * 6;
    int start = 0;

    //std::cout << start << " " << end << std::endl;

    for (int n = 0; n < 9; n += 3) {
      for (int m = start; m < end; m += 6) {
        c_0_0 = _mm_load_pd(&c[((n + 0) * ORDER_NUMBER) + m]);
        c_0_1 = _mm_load_pd(&c[((n + 1) * ORDER_NUMBER) + m]);
        c_0_2 = _mm_load_pd(&c[((n + 2) * ORDER_NUMBER) + m]);
        c_1_0 = _mm_load_pd(&c[((n + 0) * ORDER_NUMBER) + m + 2]);
        c_1_1 = _mm_load_pd(&c[((n + 1) * ORDER_NUMBER) + m + 2]);
        c_1_2 = _mm_load_pd(&c[((n + 2) * ORDER_NUMBER) + m + 2]);
        c_2_0 = _mm_load_pd(&c[((n + 0) * ORDER_NUMBER) + m + 4]);
        c_2_1 = _mm_load_pd(&c[((n + 1) * ORDER_NUMBER) + m + 4]);
        c_2_2 = _mm_load_pd(&c[((n + 2) * ORDER_NUMBER) + m + 4]);

#pragma unroll (ORDER_NUMBER)

        for (int k = 0; k < ORDER_NUMBER; k++) {
          b_0 = _mm_loaddup_pd(&b[((n + 0) * ORDER_NUMBER) + k]);
          b_1 = _mm_loaddup_pd(&b[((n + 1) * ORDER_NUMBER) + k]);
          b_2 = _mm_loaddup_pd(&b[((n + 2) * ORDER_NUMBER) + k]);

          a_0 = _mm_load_pd(&a[(k * ORDER_NUMBER) + m]);
          c_0_0 = _mm_add_pd(c_0_0, _mm_mul_pd(a_0, b_0));
          c_0_1 = _mm_add_pd(c_0_1, _mm_mul_pd(a_0, b_1));
          c_0_2 = _mm_add_pd(c_0_2, _mm_mul_pd(a_0, b_2));

          a_1 = _mm_load_pd(&a[(k * ORDER_NUMBER) + m + 2]);
          c_1_0 = _mm_add_pd(c_1_0, _mm_mul_pd(a_1, b_0));
          c_1_1 = _mm_add_pd(c_1_1, _mm_mul_pd(a_1, b_1));
          c_1_2 = _mm_add_pd(c_1_2, _mm_mul_pd(a_1, b_2));

          a_2 = _mm_load_pd(&a[(k * ORDER_NUMBER) + m + 4]);
          c_2_0 = _mm_add_pd(c_2_0, _mm_mul_pd(a_2, b_0));
          c_2_1 = _mm_add_pd(c_2_1, _mm_mul_pd(a_2, b_1));
          c_2_2 = _mm_add_pd(c_2_2, _mm_mul_pd(a_2, b_2));
        }

        _mm_store_pd(&c[((n + 0)*ORDER_NUMBER) + m], c_0_0);
        _mm_store_pd(&c[((n + 1)*ORDER_NUMBER) + m], c_0_1);
        _mm_store_pd(&c[((n + 2)*ORDER_NUMBER) + m], c_0_2);
        _mm_store_pd(&c[((n + 0)*ORDER_NUMBER) + m + 2], c_1_0);
        _mm_store_pd(&c[((n + 1)*ORDER_NUMBER) + m + 2], c_1_1);
        _mm_store_pd(&c[((n + 2)*ORDER_NUMBER) + m + 2], c_1_2);
        _mm_store_pd(&c[((n + 0)*ORDER_NUMBER) + m + 4], c_2_0);
        _mm_store_pd(&c[((n + 1)*ORDER_NUMBER) + m + 4], c_2_1);
        _mm_store_pd(&c[((n + 2)*ORDER_NUMBER) + m + 4], c_2_2);
      }
    }

    int m = end;

    for (int n = 0; n < 9; n += 3) {
      c_0_0 = _mm_load_pd(&c[((n + 0) * ORDER_NUMBER) + m]);
      c_0_1 = _mm_load_pd(&c[((n + 1) * ORDER_NUMBER) + m]);
      c_0_2 = _mm_load_pd(&c[((n + 2) * ORDER_NUMBER) + m]);
      c_1_0 = _mm_load_pd(&c[((n + 0) * ORDER_NUMBER) + m + 2]);
      c_1_1 = _mm_load_pd(&c[((n + 1) * ORDER_NUMBER) + m + 2]);
      c_1_2 = _mm_load_pd(&c[((n + 2) * ORDER_NUMBER) + m + 2]);
      c_2_0 = _mm_load_sd(&c[((n + 0) * ORDER_NUMBER) + m + 4]);
      c_2_1 = _mm_load_sd(&c[((n + 1) * ORDER_NUMBER) + m + 4]);
      c_2_2 = _mm_load_sd(&c[((n + 2) * ORDER_NUMBER) + m + 4]);

#pragma unroll (ORDER_NUMBER)

      for (int k = 0; k < ORDER_NUMBER; k++) {
        b_0 = _mm_loaddup_pd(&b[((n + 0) * ORDER_NUMBER) + k]);
        b_1 = _mm_loaddup_pd(&b[((n + 1) * ORDER_NUMBER) + k]);
        b_2 = _mm_loaddup_pd(&b[((n + 2) * ORDER_NUMBER) + k]);

        a_0 = _mm_load_pd(&a[(k * ORDER_NUMBER) + m]);
        c_0_0 = _mm_add_pd(c_0_0, _mm_mul_pd(a_0, b_0));
        c_0_1 = _mm_add_pd(c_0_1, _mm_mul_pd(a_0, b_1));
        c_0_2 = _mm_add_pd(c_0_2, _mm_mul_pd(a_0, b_2));

        a_1 = _mm_load_pd(&a[(k * ORDER_NUMBER) + m + 2]);
        c_1_0 = _mm_add_pd(c_1_0, _mm_mul_pd(a_1, b_0));
        c_1_1 = _mm_add_pd(c_1_1, _mm_mul_pd(a_1, b_1));
        c_1_2 = _mm_add_pd(c_1_2, _mm_mul_pd(a_1, b_2));

        a_2 = _mm_load_sd(&a[(k * ORDER_NUMBER) + m + 4]);
        c_2_0 = _mm_add_sd(c_2_0, _mm_mul_sd(a_2, b_0));
        c_2_1 = _mm_add_sd(c_2_1, _mm_mul_sd(a_2, b_1));
        c_2_2 = _mm_add_sd(c_2_2, _mm_mul_sd(a_2, b_2));
      }

      _mm_store_pd(&c[((n + 0)*ORDER_NUMBER) + m], c_0_0);
      _mm_store_pd(&c[((n + 1)*ORDER_NUMBER) + m], c_0_1);
      _mm_store_pd(&c[((n + 2)*ORDER_NUMBER) + m], c_0_2);
      _mm_store_pd(&c[((n + 0)*ORDER_NUMBER) + m + 2], c_1_0);
      _mm_store_pd(&c[((n + 1)*ORDER_NUMBER) + m + 2], c_1_1);
      _mm_store_pd(&c[((n + 2)*ORDER_NUMBER) + m + 2], c_1_2);
      _mm_store_sd(&c[((n + 0)*ORDER_NUMBER) + m + 4], c_2_0);
      _mm_store_sd(&c[((n + 1)*ORDER_NUMBER) + m + 4], c_2_1);
      _mm_store_sd(&c[((n + 2)*ORDER_NUMBER) + m + 4], c_2_2);
    }

#endif

#if 0
    saveend = end;
    newrange = ORDER_NUMBER - end;
    start = end;
    end = saveend + (newrange / 4) * 4;

    //std::cout << start << " " << end << std::endl;

    for (int n = 0; n < 9; n += 3) {
      for (int m = start; m < end; m += 4) {
        c_0_0 = _mm_load_pd(&c[((n + 0) * ORDER_NUMBER) + m]);
        c_0_1 = _mm_load_pd(&c[((n + 1) * ORDER_NUMBER) + m]);
        c_0_2 = _mm_load_pd(&c[((n + 2) * ORDER_NUMBER) + m]);
        c_1_0 = _mm_load_pd(&c[((n + 0) * ORDER_NUMBER) + m + 2]);
        c_1_1 = _mm_load_pd(&c[((n + 1) * ORDER_NUMBER) + m + 2]);
        c_1_2 = _mm_load_pd(&c[((n + 2) * ORDER_NUMBER) + m + 2]);

#pragma unroll (ORDER_NUMBER)

        for (int k = 0; k < ORDER_NUMBER; k++) {
          b_0 = _mm_loaddup_pd(&b[((n + 0) * ORDER_NUMBER) + k]);
          b_1 = _mm_loaddup_pd(&b[((n + 1) * ORDER_NUMBER) + k]);
          b_2 = _mm_loaddup_pd(&b[((n + 2) * ORDER_NUMBER) + k]);

          a_0 = _mm_load_pd(&a[(k * ORDER_NUMBER) + m]);
          c_0_0 = _mm_add_pd(c_0_0, _mm_mul_pd(a_0, b_0));
          c_0_1 = _mm_add_pd(c_0_1, _mm_mul_pd(a_0, b_1));
          c_0_2 = _mm_add_pd(c_0_2, _mm_mul_pd(a_0, b_2));

          a_1 = _mm_load_pd(&a[(k * ORDER_NUMBER) + m + 2]);
          c_1_0 = _mm_add_pd(c_1_0, _mm_mul_pd(a_1, b_0));
          c_1_1 = _mm_add_pd(c_1_1, _mm_mul_pd(a_1, b_1));
          c_1_2 = _mm_add_pd(c_1_2, _mm_mul_pd(a_1, b_2));

        }

        _mm_store_pd(&c[((n + 0)*ORDER_NUMBER) + m], c_0_0);
        _mm_store_pd(&c[((n + 1)*ORDER_NUMBER) + m], c_0_1);
        _mm_store_pd(&c[((n + 2)*ORDER_NUMBER) + m], c_0_2);
        _mm_store_pd(&c[((n + 0)*ORDER_NUMBER) + m + 2], c_1_0);
        _mm_store_pd(&c[((n + 1)*ORDER_NUMBER) + m + 2], c_1_1);
        _mm_store_pd(&c[((n + 2)*ORDER_NUMBER) + m + 2], c_1_2);
      }
    }


    saveend = end;
    newrange = ORDER_NUMBER - end;
    start = end;
    end = saveend + (newrange / 2) * 2;
    int rem = (ORDER_NUMBER - end);

    //std::cout << start << " " << end << std::endl;
    for (int n = 0; n < 9; n += 3) {
      for (int m = start; m < end; m += 2) {
        c_0_0 = _mm_load_pd(&c[((n + 0) * ORDER_NUMBER) + m]);
        c_0_1 = _mm_load_pd(&c[((n + 1) * ORDER_NUMBER) + m]);
        c_0_2 = _mm_load_pd(&c[((n + 2) * ORDER_NUMBER) + m]);

#pragma unroll (ORDER_NUMBER)

        for (int k = 0; k < ORDER_NUMBER; k++) {
          b_0 = _mm_loaddup_pd(&b[((n + 0) * ORDER_NUMBER) + k]);
          b_1 = _mm_loaddup_pd(&b[((n + 1) * ORDER_NUMBER) + k]);
          b_2 = _mm_loaddup_pd(&b[((n + 2) * ORDER_NUMBER) + k]);

          a_0 = _mm_load_pd(&a[(k * ORDER_NUMBER) + m]);
          c_0_0 = _mm_add_pd(c_0_0, _mm_mul_pd(a_0, b_0));
          c_0_1 = _mm_add_pd(c_0_1, _mm_mul_pd(a_0, b_1));
          c_0_2 = _mm_add_pd(c_0_2, _mm_mul_pd(a_0, b_2));

        }

        _mm_store_pd(&c[((n + 0)*ORDER_NUMBER) + m], c_0_0);
        _mm_store_pd(&c[((n + 1)*ORDER_NUMBER) + m], c_0_1);
        _mm_store_pd(&c[((n + 2)*ORDER_NUMBER) + m], c_0_2);
      }
    }


    if (rem > 0) {
      //std::cout << "rem 1 " << std::endl;
      int m = ORDER_NUMBER - 1;

      for (int n = 0; n < 9; n += 3) {
        c_0_0 = _mm_load_sd(&c[((n + 0) * ORDER_NUMBER) + m]);
        c_0_1 = _mm_load_sd(&c[((n + 1) * ORDER_NUMBER) + m]);
        c_0_2 = _mm_load_sd(&c[((n + 2) * ORDER_NUMBER) + m]);

#pragma unroll (ORDER_NUMBER)

        for (int k = 0; k < ORDER_NUMBER; k++) {
          b_0 = _mm_load_sd(&b[((n + 0) * ORDER_NUMBER) + k]);
          b_1 = _mm_load_sd(&b[((n + 1) * ORDER_NUMBER) + k]);
          b_2 = _mm_load_sd(&b[((n + 2) * ORDER_NUMBER) + k]);

          a_0 = _mm_load_sd(&a[(k * ORDER_NUMBER) + m]);
          c_0_0 = _mm_add_sd(c_0_0, _mm_mul_sd(a_0, b_0));
          c_0_1 = _mm_add_sd(c_0_1, _mm_mul_sd(a_0, b_1));
          c_0_2 = _mm_add_sd(c_0_2, _mm_mul_sd(a_0, b_2));

        }

        _mm_store_sd(&c[((n + 0)*ORDER_NUMBER) + m], c_0_0);
        _mm_store_sd(&c[((n + 1)*ORDER_NUMBER) + m], c_0_1);
        _mm_store_sd(&c[((n + 2)*ORDER_NUMBER) + m], c_0_2);

      }
    }

#endif

  }

  gettimeofday(&end, NULL);
  total = sec(start, end);

  std::cout << total << "s for initrinsics" << std::endl;
  std::cout << ((double)(REPS * ORDER_NUMBER * ORDER_NUMBER) * 9.0 * 2.0) / (total * 1.0e9) << " GFLOPS for initrinsics" << std::endl;

  // check result
  double max_error = 0.0;

  for (int i = 0; i < ORDER_NUMBER; i++) {
    for (int j = 0; j < 9; j++) {
      //std::cout << c_gold[(j*ORDER_NUMBER)+i] << " " << c[(j*ORDER_NUMBER)+i] << std::endl;
      if (max_error < fabs( c_gold[(j * ORDER_NUMBER) + i] - c[(j * ORDER_NUMBER) + i]))
        max_error = fabs( c_gold[(j * ORDER_NUMBER) + i] - c[(j * ORDER_NUMBER) + i]);
    }
  }

  std::cout << "max. error: " << max_error << std::endl;

  // free
  _mm_free(a);
  _mm_free(b);
  _mm_free(c);
  _mm_free(c_gold);

  return 0;
}
