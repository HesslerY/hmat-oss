/*
  HMat-OSS (HMatrix library, open source software)

  Copyright (C) 2014-2015 Airbus Group SAS

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

  http://github.com/jeromerobert/hmat-oss
*/

/*! \file
  \ingroup HMatrix
  \brief Dense Matrix implementation.
*/
#include "config.h"

#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <cmath>
#endif

#include "scalar_array.hpp"

#include "data_types.hpp"
#include "lapack_overloads.hpp"
#include "blas_overloads.hpp"
#include "lapack_exception.hpp"
#include "common/memory_instrumentation.hpp"
#include "system_types.h"
#include "common/my_assert.h"
#include "common/context.hpp"

#include <cstring> // memset
#include <algorithm> // swap
#include <iostream>
#include <fstream>
#include <cmath>
#include <fcntl.h>

#ifndef _WIN32
#include <sys/mman.h> // mmap
#endif

#include <sys/stat.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <stdlib.h>

#ifdef HAVE_JEMALLOC
#define JEMALLOC_NO_DEMANGLE
#include <jemalloc/jemalloc.h>
#endif

#ifdef _MSC_VER
// Intel compiler defines isnan in global namespace
// MSVC defines _isnan
# ifndef __INTEL_COMPILER
#  define isnan _isnan
# endif
#elif __GLIBC__ == 2 && __GLIBC_MINOR__ < 23
// https://sourceware.org/bugzilla/show_bug.cgi?id=19439
#elif __cplusplus >= 201103L || !defined(__GLIBC__)
using std::isnan;
#endif

namespace hmat {

/** ScalarArray */
template<typename T>
ScalarArray<T>::ScalarArray(T* _m, int _rows, int _cols, int _lda)
  : ownsMemory(false), m(_m), rows(_rows), cols(_cols), lda(_lda) {
  if (lda == -1) {
    lda = rows;
  }
  assert(lda >= rows);
}

// #define POISON_ALLOCATION
#ifdef POISON_ALLOCATION
/*! \brief Fill an array with NaNs.

  The purpose of this function is to help spotting initialized memory
  earlier by making sure that any code using uninitialized memory
  encounters NaNs.
 */
template<typename T> void poisonArray(T* array, size_t n);

template<> static void poisonArray(S_t* array, size_t n) {
  const float nanFloat = nanf("");
  for (size_t i = 0; i < n; i++) {
    array[i] = nanFloat;
  }
}
template<> static void poisonArray(D_t* array, size_t n) {
  const double nanDouble = nan("");
  for (size_t i = 0; i < n; i++) {
    array[i] = nanDouble;
  }
}
template<> static void poisonArray(C_t* array, size_t n) {
  poisonArray<S_t>((S_t*) array, 2 * n);
}
template<> static void poisonArray(Z_t* array, size_t n) {
  poisonArray<D_t>((D_t*) array, 2 * n);
}
#endif

template<typename T>
ScalarArray<T>::ScalarArray(int _rows, int _cols)
  : ownsMemory(true), rows(_rows), cols(_cols), lda(_rows) {
  size_t size = ((size_t) rows) * cols * sizeof(T);
#ifdef HAVE_JEMALLOC
  m = (T*) je_calloc(size, 1);
#else
  m = (T*) calloc(size, 1);
#endif
  HMAT_ASSERT_MSG(m, "Trying to allocate %ldb of memory failed (rows=%d cols=%d sizeof(T)=%d)", size, rows, cols, sizeof(T));
  MemoryInstrumenter::instance().alloc(size, MemoryInstrumenter::FULL_MATRIX);
#ifdef POISON_ALLOCATION
  // This memory is not initialized, fill it with NaNs to force a
  // crash when using it.
  poisonArray<T>(m, ((size_t) rows) * cols);
#endif
}

template<typename T>
ScalarArray<T>* ScalarArray<T>::Zero(int rows, int cols) {
  ScalarArray<T>* result = new ScalarArray<T>(rows, cols);
#ifdef POISON_ALLOCATION
  // The memory was poisoned in ScalarArray<T>::ScalarArray(), set it back to 0;
  result->clear();
#endif
  return result;
}

template<typename T> ScalarArray<T>::~ScalarArray() {
  if (ownsMemory) {
    size_t size = ((size_t) rows) * cols * sizeof(T);
    MemoryInstrumenter::instance().free(size, MemoryInstrumenter::FULL_MATRIX);
#ifdef HAVE_JEMALLOC
    je_free(m);
#else
    free(m);
#endif
    m = NULL;
  }
}

template<typename T> void ScalarArray<T>::clear() {
  assert(lda == rows);
  size_t size = ((size_t) rows) * cols * sizeof(T);
  memset(m, 0, size);
}

template<typename T> size_t ScalarArray<T>::storedZeros() {
  size_t result = 0;
  for (int col = 0; col < cols; col++) {
    for (int row = 0; row < rows; row++) {
      if (std::abs(get(row, col)) < 1e-16) {
        result++;
      }
    }
  }
  return result;
}

template<typename T> void ScalarArray<T>::scale(T alpha) {
  increment_flops(Multipliers<T>::mul * ((size_t) rows) * cols);
  if (lda == rows) {
    if (alpha == Constants<T>::zero) {
      this->clear();
    } else {
      // Warning: check for overflow
      size_t nm = ((size_t) rows) * cols;
      const size_t block_size_blas = 1 << 30;
      while (nm > block_size_blas) {
        proxy_cblas::scal(block_size_blas, alpha, m + nm - block_size_blas, 1);
        nm -= block_size_blas;
      }
      proxy_cblas::scal(nm, alpha, m, 1);
    }
  } else {
    T* x = m;
    if (alpha == Constants<T>::zero) {
      for (int col = 0; col < cols; col++) {
        memset(x, 0, sizeof(T) * rows);
        x += lda;
      }
    } else {
      for (int col = 0; col < cols; col++) {
        proxy_cblas::scal(rows, alpha, x, 1);
        x += lda;
      }
    }
  }
}

template<typename T> void ScalarArray<T>::transpose() {
  assert(lda == rows);
  assert(m);
#ifdef HAVE_MKL_IMATCOPY
  proxy_mkl::imatcopy(rows, cols, m);
  std::swap(rows, cols);
  lda = rows;
#else
  if (rows == cols) {
    // "Fast" path
    for (int col = 0; col < cols; col++) {
      for (int row = 0; row < col; row++) {
        T tmp = get(row, col);
        get(row, col) = get(col, row);
        get(col, row) = tmp;
      }
    }
  } else {
    ScalarArray<T> tmp(rows, cols);
    tmp.copyMatrixAtOffset(this, 0, 0);
    std::swap(rows, cols);
    lda = rows;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        get(i, j) = tmp.get(j, i);
      }
    }
  }
#endif
}

template<typename T> ScalarArray<T>* ScalarArray<T>::copy(ScalarArray<T>* result) const {
  if(result == NULL)
    result = new ScalarArray<T>(rows, cols);

  if (lda == rows) {
    size_t size = ((size_t) rows) * cols * sizeof(T);
    memcpy(result->m, m, size);
  } else {
    for (int col = 0; col < cols; col++) {
      size_t resultOffset = ((size_t) result->rows) * col;
      size_t offset = ((size_t) lda) * col;
      memcpy(result->m + resultOffset, m + offset, rows * sizeof(T));
    }
  }

  return result;
}

template<typename T> ScalarArray<T>* ScalarArray<T>::copyAndTranspose() const {
  ScalarArray<T>* result = new ScalarArray<T>(cols, rows);
  //result->clear();
#ifdef HAVE_MKL_IMATCOPY
  if (lda == rows) {
    proxy_mkl::omatcopy(rows, cols, m, result->m);
  } else {
#endif
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      result->get(j, i) = get(i, j);
    }
  }
#ifdef HAVE_MKL_IMATCOPY
  }
#endif
  return result;
}


template<typename T>
void ScalarArray<T>::gemm(char transA, char transB, T alpha,
                         const ScalarArray<T>* a, const ScalarArray<T>* b,
                         T beta) {
  int m = (transA == 'N' ? a->rows : a->cols);
  int n = (transB == 'N' ? b->cols : b->rows);
  int k = (transA == 'N' ? a->cols : a->rows);
  assert(a->lda >= (transA == 'N' ? m : k));
  assert(b->lda >= (transB == 'N' ? k : n));
  assert(rows == m);
  assert(cols == n);
  {
    const size_t _m = m, _n = n, _k = k;
    const size_t adds = _m * _n * _k;
    const size_t muls = _m * _n * _k;
    increment_flops(Multipliers<T>::add * adds + Multipliers<T>::mul * muls);
  }
  proxy_cblas::gemm(transA, transB, m, n, k, alpha, a->m, a->lda, b->m, b->lda,
                    beta, this->m, this->lda);
}

template<typename T>
void ScalarArray<T>::copyMatrixAtOffset(const ScalarArray<T>* a,
                                       int rowOffset, int colOffset) {
  assert(rowOffset + a->rows <= rows);
  assert(colOffset + a->cols <= cols);


  // Use memcpy when copying the whole matrix. This avoids BLAS calls.
  if ((rowOffset == 0) && (colOffset == 0)
      && (a->rows == rows) && (a->cols == cols)
      && (a->lda == a->rows) && (lda == rows)) {
    size_t size = ((size_t) rows) * cols;
    memcpy(m, a->m, size * sizeof(T));
    return;
  }

  for (int col = 0; col < a->cols; col++) {
    proxy_cblas::copy(a->rows, a->m + col * a->lda, 1,
                m + rowOffset + ((colOffset + col) * lda), 1);
  }
}

template<typename T>
void ScalarArray<T>::copyMatrixAtOffset(const ScalarArray<T>* a,
                                       int rowOffset, int colOffset,
                                       int rowsToCopy, int colsToCopy) {
  assert(rowOffset + rowsToCopy <= rows);
  assert(colOffset + colsToCopy <= cols);
  for (int col = 0; col < colsToCopy; col++) {
    proxy_cblas::copy(rowsToCopy, a->m + col * a->lda, 1,
                (m + rowOffset + ((colOffset + col) * lda)), 1);
  }
}

template<typename T>
void ScalarArray<T>::axpy(T alpha, const ScalarArray<T>* a) {
  assert(rows == a->rows);
  assert(cols == a->cols);
  size_t size = ((size_t) rows) * cols;

  increment_flops(Multipliers<T>::add * size
		  + (alpha == Constants<T>::pone ? 0 : Multipliers<T>::mul * size));
  // Fast path
  if ((lda == rows) && (a->lda == a->rows) && (size < 1000000000)) {
    proxy_cblas::axpy(size, alpha, a->m, 1, m, 1);
    return;
  }

  for (int col = 0; col < cols; col++) {
    proxy_cblas::axpy(rows, alpha, a->m + ((size_t) col) * a->lda, 1, m + ((size_t) col) * lda, 1);
  }
}

template<typename T>
double ScalarArray<T>::normSqr() const {
  size_t size = ((size_t) rows) * cols;
  T result = Constants<T>::zero;

  // Fast path
  if ((size < 1000000000) && (lda == rows)) {
    result += proxy_cblas_convenience::dot_c(size, m, 1, m, 1);
    return hmat::real(result);
  }
  for (int col = 0; col < cols; col++) {
    result += proxy_cblas_convenience::dot_c(rows, m + col * lda, 1, m + col * lda, 1);
  }
  return hmat::real(result);
}

template<typename T> double ScalarArray<T>::norm() const {
  return sqrt(normSqr());
}

template<typename T> void ScalarArray<T>::fromFile(const char * filename) {
  FILE * f = fopen(filename, "rb");
  /* Read the header before data : [stype, rows, cols, sieof(T), 0] */
  int code;
  int r = fread(&code, sizeof(int), 1, f);
  HMAT_ASSERT(r == 1);
  HMAT_ASSERT(code == Constants<T>::code);
  r = fread(&rows, sizeof(int), 1, f);
  lda = rows;
  HMAT_ASSERT(r == 1);
  r = fread(&cols, sizeof(int), 1, f);
  HMAT_ASSERT(r == 1);
  r = fseek(f, 2 * sizeof(int), SEEK_CUR);
  HMAT_ASSERT(r == 0);
  if(m)
      free(m);
  size_t size = ((size_t) rows) * cols * sizeof(T);
  m = (T*) calloc(size, 1);
  r = fread(m, size, 1, f);
  fclose(f);
  HMAT_ASSERT(r == 1);
}

template<typename T> void ScalarArray<T>::toFile(const char *filename) const {
  int ierr;
  int fd;
  size_t size = ((size_t) rows) * cols * sizeof(T) + 5 * sizeof(int);

  HMAT_ASSERT(lda == rows);

  fd = open(filename, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
  HMAT_ASSERT(fd != -1);
  ierr = lseek(fd, size - 1, SEEK_SET);
  HMAT_ASSERT(ierr != -1);
  ierr = write(fd, "", 1);
  HMAT_ASSERT(ierr == 1);
#ifndef _WIN32
  void* mmapedFile = mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  ierr = (mmapedFile == MAP_FAILED) ? 1 : 0;
  HMAT_ASSERT(!ierr);
  /* Write the header before data : [stype, rows, cols, sieof(T), 0] */
  int *asIntArray = (int*) mmapedFile;
  asIntArray[0] = Constants<T>::code;
  asIntArray[1] = rows;
  asIntArray[2] = cols;
  asIntArray[3] = sizeof(T);
  asIntArray[4] = 0;
  asIntArray += 5;
  T* mat = (T*) asIntArray;
  memcpy(mat, m, size - 5 * sizeof(int));
  close(fd);
  munmap(mmapedFile, size);
#else
  HMAT_ASSERT_MSG(false, "mmap not available on this platform");
#endif

}

template<typename T> size_t ScalarArray<T>::memorySize() const {
   return ((size_t) rows) * cols * sizeof(T);
}

template<typename T> void checkNanReal(const ScalarArray<T>* m) {
  for (int col = 0; col < m->cols; col++) {
    for (int row = 0; row < m->rows; row++) {
      HMAT_ASSERT(!isnan(m->get(row, col)));
    }
  }
}

template<typename T> void checkNanComplex(const ScalarArray<T>* m) {
  for (int col = 0; col < m->cols; col++) {
    for (int row = 0; row < m->rows; row++) {
      HMAT_ASSERT(!isnan(m->get(row, col).real()));
      HMAT_ASSERT(!isnan(m->get(row, col).imag()));
    }
  }
}

template<> void ScalarArray<S_t>::checkNan() const {
  checkNanReal(this);
}
template<> void ScalarArray<D_t>::checkNan() const {
  checkNanReal(this);
}
template<> void ScalarArray<C_t>::checkNan() const {
  checkNanComplex(this);
}
template<> void ScalarArray<Z_t>::checkNan() const {
  checkNanComplex(this);
}

// the classes declaration
template class ScalarArray<S_t>;
template class ScalarArray<D_t>;
template class ScalarArray<C_t>;
template class ScalarArray<Z_t>;

}  // end namespace hmat
