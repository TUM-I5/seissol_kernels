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
#define ORDER_NUMBER 56
#endif

#ifndef GAUSS_NUMBER
#define GAUSS_NUMBER 52
#endif

#define N_ONE 36

#define REPS 100000

inline double sec(struct timeval start, struct timeval end) {
  return ((double)(((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)))) / 1.0e6;
}


void run_A_square() {
  // allocate
  double* a = (double*)_mm_malloc(GAUSS_NUMBER * ORDER_NUMBER * sizeof(double), 64);
  double* b = (double*)_mm_malloc(ORDER_NUMBER * N_ONE * sizeof(double), 64);
  double* c = (double*)_mm_malloc(GAUSS_NUMBER * N_ONE * sizeof(double), 64);
  double* c_gold = (double*)_mm_malloc(GAUSS_NUMBER * N_ONE * sizeof(double), 64);

  // touch
  for (int i = 0; i < GAUSS_NUMBER; i++) {
    for (int j = 0; j < ORDER_NUMBER; j++) {
      a[(j * GAUSS_NUMBER) + i] = (double)(i + (j * GAUSS_NUMBER));
    }
  }

  for (int i = 0; i < ORDER_NUMBER; i++) {
    for (int j = 0; j < N_ONE; j++) {
      b[(j * ORDER_NUMBER) + i] = (double)(i + (j * 9));
    }
  }

  for (int i = 0; i < GAUSS_NUMBER; i++) {
    for (int j = 0; j < N_ONE; j++) {
      c[(j * GAUSS_NUMBER) + i] = 0.0;
      c_gold[(j * GAUSS_NUMBER) + i] = 0.0;
    }
  }

  // C routine
  struct timeval start, end;
  gettimeofday(&start, NULL);

  for (int t = 0; t < REPS; t++) {
    for (int n = 0; n < N_ONE; n++) {
      for (int k = 0; k < ORDER_NUMBER; k++) {
#pragma simd

        for (int m = 0; m < GAUSS_NUMBER; m++) {
          c_gold[(n * GAUSS_NUMBER) + m] += a[(k * GAUSS_NUMBER) + m] * b[(n * ORDER_NUMBER) + k];
        }
      }
    }
  }

  gettimeofday(&end, NULL);
  double total = sec(start, end);

  std::cout << total << "s for C" << std::endl;
  std::cout << ((double)(REPS * GAUSS_NUMBER * ORDER_NUMBER) * N_ONE * 2.0) / (total * 1.0e9) << " GFLOPS for C" << std::endl;

  gettimeofday(&start, NULL);

  for (int t = 0; t < REPS; t++) {
    dense_test_square(a, b, c);
  }

  gettimeofday(&end, NULL);
  total = sec(start, end);

  std::cout << total << "s for initrinsics" << std::endl;
  std::cout << ((double)(REPS * GAUSS_NUMBER * ORDER_NUMBER) * N_ONE * 2.0) / (total * 1.0e9) << " GFLOPS for initrinsics" << std::endl;

  // check result
  double max_error = 0.0;

  for (int i = 0; i < GAUSS_NUMBER; i++) {
    for (int j = 0; j < N_ONE; j++) {
      //std::cout << c_gold[(j*GAUSS_NUMBER)+i] << " " << c[(j*GAUSS_NUMBER)+i] << std::endl;
      if (max_error < fabs( c_gold[(j * GAUSS_NUMBER) + i] - c[(j * GAUSS_NUMBER) + i]))
        max_error = fabs( c_gold[(j * GAUSS_NUMBER) + i] - c[(j * GAUSS_NUMBER) + i]);
    }
  }

  std::cout << "max. error: " << max_error << std::endl;

  // free
  _mm_free(a);
  _mm_free(b);
  _mm_free(c);
  _mm_free(c_gold);
}

void run_A_rect() {
  // allocate
  double* a = (double*)_mm_malloc(ORDER_NUMBER * 9 * sizeof(double), 64);
  double* b = (double*)_mm_malloc(9 * 6 * sizeof(double), 64);
  double* c = (double*)_mm_malloc(ORDER_NUMBER * 6 * sizeof(double), 64);
  double* c_gold = (double*)_mm_malloc(ORDER_NUMBER * 6 * sizeof(double), 64);

  // touch
  for (int i = 0; i < ORDER_NUMBER; i++) {
    for (int j = 0; j < 9; j++) {
      a[(j * ORDER_NUMBER) + i] = (double)(i + (j * ORDER_NUMBER));
      if (j < 6) c[(j * ORDER_NUMBER) + i] = 0.0;
      if (j < 6) c_gold[(j * ORDER_NUMBER) + i] = 0.0;

      if (i < 9 && j < 6) b[(j * 9) + i] = (double)(i + (j * 9));
    }
  }

  // C routine
  struct timeval start, end;
  gettimeofday(&start, NULL);

  for (int t = 0; t < REPS; t++) {
    for (int n = 0; n < 6; n++) {
      for (int k = 0; k < 9; k++) {
#pragma simd

        for (int m = 0; m < ORDER_NUMBER; m++) {
          c_gold[(n * ORDER_NUMBER) + m] += a[(k * ORDER_NUMBER) + m] * b[(n * 9) + k];
        }
      }
    }
  }

  gettimeofday(&end, NULL);
  double total = sec(start, end);

  std::cout << total << "s for C" << std::endl;
  std::cout << ((double)(REPS * ORDER_NUMBER * 9.0) * 6.0 * 2.0) / (total * 1.0e9) << " GFLOPS for C" << std::endl;

  gettimeofday(&start, NULL);

  for (int t = 0; t < REPS; t++) {
    dense_test_rect(a, b, c);
  }

  gettimeofday(&end, NULL);
  total = sec(start, end);

  std::cout << total << "s for initrinsics" << std::endl;
  std::cout << ((double)(REPS * ORDER_NUMBER * 9.0) * 6.0 * 2.0) / (total * 1.0e9) << " GFLOPS for initrinsics" << std::endl;

  // check result
  double max_error = 0.0;

  for (int i = 0; i < ORDER_NUMBER; i++) {
    for (int j = 0; j < 6; j++) {
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
}


int main(int argc, char* argv[]) {
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "RUNNING (" << GAUSS_NUMBER << "x" << ORDER_NUMBER << ") X (" << ORDER_NUMBER << "x" << N_ONE << ") = (" << GAUSS_NUMBER << "x" << N_ONE << ")"  << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
  run_A_square();
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "RUNNING (" << ORDER_NUMBER << "x" << 9 << ") X (" << 9 << "x" << 6 << ") = (" << ORDER_NUMBER << "x" << 6 << ")"  << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
  run_A_rect();
  std::cout << "------------------------------------------------" << std::endl;
  return 0;
}
