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

#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <cmath>
#endif

#include "compression.hpp"

#include <vector>
#include <cfloat>
#include <cstring>

#include "cluster_tree.hpp"
#include "assembly.hpp"
#include "rk_matrix.hpp"
#include "lapack_operations.hpp"
#include "lapack_overloads.hpp"
#include "blas_overloads.hpp"
#include "full_matrix.hpp"
#include "common/context.hpp"
#include "common/my_assert.h"
#include "cluster_assembly_function.hpp"
#include "random_pivot_manager.hpp"

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

using std::vector;
using std::min;
using std::max;

namespace hmat {

/** \brief Updates a row to reflect its current value in the matrix.

    \param rowVec the row to update
    \param row the index of this row
    \param rows the rows that have already been computed
    \param cols the cols that have already been computed
    \params k the number of rows that have already been computed
 */
template<typename T>
static void updateRow(Vector<T>& rowVec, int row, const vector<Vector<T>*>& rows,
                      const vector<Vector<T>*>& cols, int k) {
  for (int l = 0; l < k; l++) {
    rowVec.axpy(Constants<T>::mone * (*cols[l])[row], rows[l]);
  }
}

template<typename T>
static void updateCol(Vector<T>& colVec, int col, const vector<Vector<T>*>& cols,
                      const vector<Vector<T>*>& rows, int k) {
  updateRow<T>(colVec, col, cols, rows, k);
}


template<typename T> static void findMax(const FullMatrix<T>* m, int& i, int& j) {
  if (m->data.lda == m->data.rows) {
    // quick path
    int k = proxy_cblas::i_amax(m->data.rows*m->data.cols, m->data.const_ptr(), 1);
    i = k % m->data.rows ;
    j = (k - i) / m->data.rows ;
  } else {
    i = 0;
    j = 0;
    double maxNorm = 0.;
    for (int col = 0; col < m->cols(); col++) {
      int row = proxy_cblas::i_amax(m->data.rows, m->data.const_ptr(0,col), 1);
      const double norm = squaredNorm<T>(m->get(row, col));
      if (norm > maxNorm) {
        i = row;
        j = col;
        maxNorm = norm;
      }
    }
  }
}


/*! \brief Find a column that is free and not null, or return an error.

  \param block The block assembly function
  \param colFree an array of the columns that are free to choose
                 from. This array is updated by this function.
  \param col the column to be returned. It doesn't need to be zeroed beforehand.'
  \return the index of the chosen column, or -1 if no column can be found.
 */
template<typename T>
static int findCol(const ClusterAssemblyFunction<T>& block, vector<bool>& colFree,
                   Vector<typename Types<T>::dp>& col) {
  int colCount = colFree.size();
  bool found = false;
  int i;
  for (i = 0; i < colCount; i++) {
    if (colFree[i]) {
      col.clear();
      block.getCol(i, col);
      colFree[i] = false;
      if (!col.isZero()) {
        found = true;
        break;
      }
    }
  }
  return (found ? i : -1);
}


template<typename T>
static int findMinRow(const ClusterAssemblyFunction<T>& block,
                      vector<bool>& rowFree,
                      const vector<Vector<typename Types<T>::dp>*>& aCols,
                      const vector<Vector<typename Types<T>::dp>*>& bCols,
                      const Vector<typename Types<T>::dp>& aRef,
                      Vector<typename Types<T>::dp>& row) {

  int rowCount = aRef.rows;
  double minNorm2;
  int i_ref;
  bool found = false;

  while (!found) {
    i_ref = -1;
    minNorm2 = DBL_MAX;
    for (int i = 0; i < rowCount; i++) {
      if (rowFree[i]) {
        double norm2 = squaredNorm<typename Types<T>::dp>(aRef[i]);
        if (norm2 < minNorm2) {
          i_ref = i;
          minNorm2 = norm2;
        }
      }
    }
    if (i_ref == -1) {
      return i_ref;
    }
    row.clear();
    block.getRow(i_ref, row);
    updateRow<typename Types<T>::dp>(row, i_ref, bCols, aCols, aCols.size());
    found = !row.isZero();
    rowFree[i_ref] = false;
  }
  return i_ref;
}

template<typename T>
static int findMinCol(const ClusterAssemblyFunction<T>& block,
                      vector<bool>& colFree,
                      const vector<Vector<typename Types<T>::dp>*>& aCols,
                      const vector<Vector<typename Types<T>::dp>*>& bCols,
                      const Vector<typename Types<T>::dp>& bRef,
                      Vector<typename Types<T>::dp>& col) {
  int colCount = bRef.rows;
  double minNorm2;
  int j_ref;
  bool found = false;

  while (!found) {
    j_ref = -1;
    minNorm2 = DBL_MAX;
    for (int j = 0; j < colCount; j++) {
      if (colFree[j]) {
        double norm2 = squaredNorm<typename Types<T>::dp>(bRef[j]);
        if (norm2 < minNorm2) {
          j_ref = j;
          minNorm2 = norm2;
        }
      }
    }
    if (j_ref == -1) {
      return j_ref;
    }
    col.clear();
    block.getCol(j_ref, col);
    updateCol<typename Types<T>::dp>(col, j_ref, aCols, bCols, bCols.size());
    found = !col.isZero();
    colFree[j_ref] = false;
  }
  return j_ref;
}


template<typename T>
RkMatrix<T>* truncatedSvd(FullMatrix<T>* m, double epsilon) {
  DECLARE_CONTEXT;

  if (m->isZero()) {
    return new RkMatrix<T>(NULL, m->rows_, NULL, m->cols_);
  }
  // In the case of non-square matrix, we don't calculate singular vectors
  // bigger than the minimum dimension of the matrix. However this is not
  // necessary here, since k < min (n, p) for M matrix (nxp).
  ScalarArray<T> *u = NULL, *v = NULL;

  // TODO compress with something else than SVD
  m->data.truncatedSvdDecomposition(&u, &v, epsilon);

  return new RkMatrix<T>(u, m->rows_, v, m->cols_);
}


template<typename T>
RkMatrix<typename Types<T>::dp>*
doCompressionSVD(const ClusterAssemblyFunction<T>& block, double compressionEpsilon) {
  DECLARE_CONTEXT;
  typedef typename Types<T>::dp dp_t;
  FullMatrix<dp_t>* m = block.assemble();
  RkMatrix<dp_t>* result = truncatedSvd(m, compressionEpsilon);
  delete m;
  return result;
}

RkMatrix<Types<S_t>::dp>*
CompressionSVD::compress(const ClusterAssemblyFunction<S_t>& block) const {
    return doCompressionSVD<S_t>(block, epsilon_);
}
RkMatrix<Types<D_t>::dp>*
CompressionSVD::compress(const ClusterAssemblyFunction<D_t>& block) const {
    return doCompressionSVD<D_t>(block, epsilon_);
}
RkMatrix<Types<C_t>::dp>*
CompressionSVD::compress(const ClusterAssemblyFunction<C_t>& block) const {
    return doCompressionSVD<C_t>(block, epsilon_);
}
RkMatrix<Types<Z_t>::dp>*
CompressionSVD::compress(const ClusterAssemblyFunction<Z_t>& block) const {
    return doCompressionSVD<Z_t>(block, epsilon_);
}


template<typename T>
RkMatrix<typename Types<T>::dp>*
doCompressionAcaFull(const ClusterAssemblyFunction<T>& block, double compressionEpsilon) {
  DECLARE_CONTEXT;
  typedef typename Types<T>::dp dp_t;
  // TODO: use ClusterAssemblyFunction to optimize with blockinfo_t
  FullMatrix<dp_t>* m = block.assemble();

  double estimateSquaredNorm = 0;
  int maxK = min(m->rows(), m->cols());

  ScalarArray<dp_t> tmpA(m->rows(), maxK);
  ScalarArray<dp_t> tmpB(m->cols(), maxK);
  int nu;

  for (nu = 0; nu < maxK; nu++) {
    int i_nu, j_nu;
    findMax(m, i_nu, j_nu);
    const dp_t delta = m->get(i_nu, j_nu);
    if (squaredNorm(delta) == 0.) {
      break;
    }

    // Creation of the vectors A_i_nu and B_j_nu
    {
      Vector<dp_t> va_nu(tmpA, nu);
      Vector<dp_t> vb_nu(tmpB, nu);

      for (int i = 0; i < m->rows(); i++)
        va_nu[i] = m->get(i, j_nu);
      for (int j = 0; j < m->cols(); j++)
        vb_nu[j] = m->get(i_nu, j) / delta;

      // performs the rank 1 operation m := m - va_nu*vb_nu^T
      // in order to nullify m->get(i_nu, j_nu) (the previous maximum value)
      m->data.rankOneUpdate(Constants<dp_t>::mone, va_nu, vb_nu);

      // Update the estimate norm
      // Let S_{k-1} be the previous estimate. We have (for the Frobenius norm):
      //  ||S_k||^2 = ||S_{k-1}||^2 + \sum_{l = 0}^{nu-1} (<a_k, a_l> <b_k, b_l> + <a_l, a_k> <b_l, b_k>))
      //              + ||a_k||^2 ||b_k||^2
      // The sum
      double newEstimate = 0.0;
      for (int l = 0; l < nu - 1; l++) {
        Vector<dp_t> a_l(tmpA, l);
        Vector<dp_t> b_l(tmpB, l);
        newEstimate += hmat::real(Vector<dp_t>::dot(&va_nu, &a_l) * Vector<dp_t>::dot(&vb_nu, &b_l));
      }
      estimateSquaredNorm += 2.0 * newEstimate;
      const double a_nu_norm_2 = va_nu.normSqr();
      const double b_nu_norm_2 = vb_nu.normSqr();
      const double ab_norm_2 = a_nu_norm_2 * b_nu_norm_2;
      estimateSquaredNorm += ab_norm_2;

      // Evaluate the stopping criterion
      // ||a_nu|| ||b_nu|| < compressionEpsilon * ||S_nu||
      // <=> ||a_nu||^2 ||b_nu||^2 < compressionEpsilon^2 ||S_nu||^2
      if (ab_norm_2 < compressionEpsilon * compressionEpsilon * estimateSquaredNorm) {
        break;
      }
    }
  }
  delete m;

  if (nu == 0) {
    return new RkMatrix<dp_t>(NULL, block.rows, NULL, block.cols);
  }

  tmpA.cols=nu; // resize tmpA and tmpB for the copy
  tmpB.cols=nu;

  return new RkMatrix<dp_t>(tmpA.copy(), block.rows, tmpB.copy(), block.cols);
}

RkMatrix<Types<S_t>::dp>*
CompressionAcaFull::compress(const ClusterAssemblyFunction<S_t>& block) const {
    return doCompressionAcaFull<S_t>(block, epsilon_);
}
RkMatrix<Types<D_t>::dp>*
CompressionAcaFull::compress(const ClusterAssemblyFunction<D_t>& block) const {
    return doCompressionAcaFull<D_t>(block, epsilon_);
}
RkMatrix<Types<C_t>::dp>*
CompressionAcaFull::compress(const ClusterAssemblyFunction<C_t>& block) const {
    return doCompressionAcaFull<C_t>(block, epsilon_);
}
RkMatrix<Types<Z_t>::dp>*
CompressionAcaFull::compress(const ClusterAssemblyFunction<Z_t>& block) const {
    return doCompressionAcaFull<Z_t>(block, epsilon_);
}

template<typename T>
RkMatrix<typename Types<T>::dp>*
doCompressionAcaPartial(const ClusterAssemblyFunction<T>& block, double compressionEpsilon, bool useRandomPivots) {
  typedef typename Types<T>::dp dp_t;

  double estimateSquaredNorm = 0;

  const int rowCount = block.rows->size();
  const int colCount = block.cols->size();
  int maxK = min(rowCount, colCount);
  // Contains false for the rows that were already used as pivot
  vector<bool> rowFree(rowCount, true);
  int rowPivotCount = 0;
  // idem for columns
  vector<bool> colFree(colCount, true);
  vector<Vector<dp_t>*> aCols;
  vector<Vector<dp_t>*> bCols;

  int I = 0;
  int J = 0;
  int k = 0;

  RandomPivotManager<T> randomPivotManager(block, useRandomPivots ? max(rowCount, colCount) : 0);

  do {
    Vector<dp_t>* bCol = new Vector<dp_t>(block.cols->size());
    // Calculation of row I and its residue
    block.getRow(I, *bCol);
    updateRow(*bCol, I, bCols, aCols, k);
    rowFree[I] = false;

    // Find max and argmax of the residue
    double maxNorm2 = 0.;
    for (int j = 0; j < colCount; j++) {
      const double norm2 = squaredNorm<dp_t>((*bCol)[j]);
      if (colFree[j] && norm2 > maxNorm2) {
        maxNorm2 = norm2;
        J = j;
      }
    }

    Pivot<dp_t > randomOrDefaultPivot = randomPivotManager.GetPivot();
    if(I!=randomOrDefaultPivot.row && squaredNorm(randomOrDefaultPivot.value) > maxNorm2){
      I = randomOrDefaultPivot.row;
      continue;
    }

    if ((*bCol)[J] == Constants<dp_t>::zero) {
      delete bCol;
      // We look for another row which has not already been used.
      I = 0;
      while (!rowFree[I]) {
        I++;
      }
    } else {
      // Find pivot and scale column B
      dp_t pivot = Constants<dp_t>::pone / (*bCol)[J];
      bCol->scale(pivot);
      bCols.push_back(bCol);

      // Compute column J and residue
      Vector<dp_t>* aCol = new Vector<dp_t>(block.rows->size());
      block.getCol(J, *aCol);
      updateCol(*aCol, J, aCols, bCols, k);
      randomPivotManager.AddUsedPivot(bCol, aCol, I, J);
      colFree[J] = false;
      aCols.push_back(aCol);

      // Find max and argmax of the residue
      maxNorm2 = 0.;
      for (int i = 0; i < rowCount; i++) {
        const double norm2 = squaredNorm<dp_t>((*aCol)[i]);
        if (rowFree[i] && norm2 > maxNorm2) {
          maxNorm2 = norm2;
          I = i;
        }
      }

      // Update the estimated norm
      // Let S_{k-1} be the previous estimate. We have (for the Frobenius norm):
      //  ||S_k||^2 = ||S_{k-1}||^2 + \sum_{l = 0}^{nu-1} (<a_k, a_l> <b_k, b_l> + <a_l, a_k> <b_l, b_k>))
      //              + ||a_k||^2 ||b_k||^2
      double newEstimate = 0.0;
      for (int l = 0; l < k; l++) {
        newEstimate += hmat::real(Vector<dp_t>::dot(aCol, aCols[l]) * Vector<dp_t>::dot(bCol, bCols[l]));
      }
      estimateSquaredNorm += 2.0 * newEstimate;
      const double aColNorm_2 = aCol->normSqr();
      const double bColNorm_2 = bCol->normSqr();
      const double ab_norm_2 = aColNorm_2 * bColNorm_2;
      estimateSquaredNorm += ab_norm_2;
      k++;

      // Evaluate the stopping criterion
      // ||a_nu|| ||b_nu|| < compressionEpsilon * ||S_nu||
      // <=> ||a_nu||^2 ||b_nu||^2 < compressionEpsilon^2 ||S_nu||^2
      if (ab_norm_2 < compressionEpsilon * compressionEpsilon * estimateSquaredNorm) {
        break;
      }
    }
    rowPivotCount++;
  } while (rowPivotCount < maxK);

  ScalarArray<dp_t> *newA, *newB;
  if (k != 0) {
    newA = new ScalarArray<dp_t>(block.rows->size(), k);
    for (int i = 0; i < k; i++) {
      newA->copyMatrixAtOffset(aCols[i], 0, i);
      delete aCols[i];
      aCols[i] = NULL;
    }
    newB = new ScalarArray<dp_t>(block.cols->size(), k);
    for (int i = 0; i < k; i++) {
      newB->copyMatrixAtOffset(bCols[i], 0, i);
      delete bCols[i];
      bCols[i] = NULL;
    }
  } else {
    // If k == 0, block is only made of zeros.
    return new RkMatrix<dp_t>(NULL, block.rows, NULL, block.cols);
  }

  return new RkMatrix<dp_t>(newA, block.rows, newB, block.cols);
}

RkMatrix<Types<S_t>::dp>*
CompressionAcaPartial::compress(const ClusterAssemblyFunction<S_t>& block) const {
    return doCompressionAcaPartial<S_t>(block, epsilon_, useRandomPivots_);
}
RkMatrix<Types<D_t>::dp>*
CompressionAcaPartial::compress(const ClusterAssemblyFunction<D_t>& block) const {
    return doCompressionAcaPartial<D_t>(block, epsilon_, useRandomPivots_);
}
RkMatrix<Types<C_t>::dp>*
CompressionAcaPartial::compress(const ClusterAssemblyFunction<C_t>& block) const {
    return doCompressionAcaPartial<C_t>(block, epsilon_, useRandomPivots_);
}
RkMatrix<Types<Z_t>::dp>*
CompressionAcaPartial::compress(const ClusterAssemblyFunction<Z_t>& block) const {
    return doCompressionAcaPartial<Z_t>(block, epsilon_, useRandomPivots_);
}


template<typename T>
RkMatrix<typename Types<T>::dp>*
doCompressionAcaPlus(const ClusterAssemblyFunction<T>& block, double compressionEpsilon) {

  typedef typename Types<T>::dp dp_t;
  double estimateSquaredNorm = 0;
  int i_ref, j_ref;
  int rowCount = block.rows->size(), colCount = block.cols->size();
  int maxK = min(rowCount, colCount);
  Vector<dp_t> bRef(colCount), aRef(rowCount);
  vector<bool> rowFree(rowCount, true), colFree(colCount, true);
  vector<Vector<dp_t>*> aCols, bCols;

  j_ref = findCol(block, colFree, aRef);
  if (j_ref == -1) {
	// The block is completely zero.
    return new RkMatrix<dp_t>(NULL, block.rows, NULL, block.cols);
  }

  // The reference row is chosen such that it intersects the reference
  // column at its argmin index.
  i_ref = findMinRow(block, rowFree, aCols, bCols, aRef, bRef);

  int k = 0;
  do {
    Vector<dp_t>* bVec = new Vector<dp_t>(colCount);
    Vector<dp_t>* aVec = new Vector<dp_t>(rowCount);
    int i_star, j_star;
    dp_t i_star_value, j_star_value;

    i_star = aRef.absoluteMaxIndex();
    i_star_value = aRef[i_star];

    j_star = bRef.absoluteMaxIndex();
    j_star_value = bRef[j_star];

    if (squaredNorm<dp_t>(i_star_value) > squaredNorm<dp_t>(j_star_value)) {
      // i_star is fixed, we look for j_star
      block.getRow(i_star, *bVec);
      // Calculate the residue
      updateRow<dp_t>(*bVec, i_star, bCols, aCols, k);
      j_star = bVec->absoluteMaxIndex();
      dp_t pivot = (*bVec)[j_star];
      // Calculate a
      block.getCol(j_star, *aVec);
      updateCol<dp_t>(*aVec, j_star, aCols, bCols, k);
      if(pivot != Constants<dp_t>::zero) aVec->scale(Constants<dp_t>::pone / pivot);
    } else {
      // j_star is fixed, we look for i_star
      block.getCol(j_star, *aVec);
      updateCol<dp_t>(*aVec, j_star, aCols, bCols, k);
      i_star = aVec->absoluteMaxIndex();
      dp_t pivot = (*aVec)[i_star];
      // Calculate b
      block.getRow(i_star, *bVec);
      updateRow<dp_t>(*bVec, i_star, bCols, aCols, k);
      if(pivot != Constants<dp_t>::zero) bVec->scale(Constants<dp_t>::pone / pivot);
    }

    rowFree[i_star] = false;
    colFree[j_star] = false;

    aCols.push_back(aVec);
    bCols.push_back(bVec);
    // Update the estimate norm
    // Let S_{k-1} be the previous estimate. We have (for the Frobenius norm):
    //  ||S_k||^2 = ||S_{k-1}||^2 + \sum_{l = 0}^{nu-1} (<a_k, a_l> <b_k, b_l> + <u_l, u_k> <b_l, b_k>))
    //              + ||a_k||^2 ||b_k||^2
    double newEstimate = 0.0;
    for (int l = 0; l < k; l++) {
      newEstimate += hmat::real(Vector<dp_t>::dot(aVec, aCols[l]) * Vector<dp_t>::dot(bVec, bCols[l]));
    }
    estimateSquaredNorm += 2.0 * newEstimate;
    const double aVecNorm_2 = aVec->normSqr();
    const double bVecNorm_2 = bVec->normSqr();
    const double ab_norm_2 = aVecNorm_2 * bVecNorm_2;
    estimateSquaredNorm += ab_norm_2;
    k++;

    // Evaluate the stopping criterion
    // ||a_nu|| ||b_nu|| < compressionEpsilon * ||S_nu||
    // <=> ||a_nu||^2 ||b_nu||^2 < compressionEpsilon^2 ||S_nu||^2
    if (ab_norm_2 < compressionEpsilon * compressionEpsilon * estimateSquaredNorm) {
      break;
    }

    // Update of a_ref and b_ref
    aRef.axpy(Constants<dp_t>::mone * (*bCols[k - 1])[j_ref], aCols[k - 1]);
    bRef.axpy(Constants<dp_t>::mone * (*aCols[k - 1])[i_ref], bCols[k - 1]);
    const bool needNewA = aRef.isZero() || (j_star == j_ref);
    const bool needNewB = bRef.isZero() || (i_star == i_ref);

    // If the row or the column of reference have been already chosen as pivot,
    // we can not keep it and then we take one or two others.
    if (needNewA && needNewB) {
      bool found = false;
      while (!found) {
        aRef.clear();
        j_ref = findCol(block, colFree, aRef);
         // We can not find non-zero column anymore, done!
        if (j_ref == -1) {
          break;
        }
        updateCol<dp_t>(aRef, j_ref, aCols, bCols, k);
        found = !aRef.isZero();
      }
      if (!found) {
        break;
      }
      bRef.clear();
      i_ref = findMinRow(block, rowFree, aCols, bCols, aRef, bRef);
      // We can not find non-zero line anymore, done!
      if (i_ref == -1) {
        break;
      }
    } else if (needNewB) {
      bRef.clear();
      i_ref = findMinRow(block, rowFree, aCols, bCols, aRef, bRef);
      // We can not find non-zero line anymore, done!
      if (i_ref == -1) {
        break;
      }
    } else if (needNewA) {
      aRef.clear();
      j_ref = findMinCol(block, colFree, aCols, bCols, bRef, aRef);
      // We can not find non-zero column anymore, done!
      if (j_ref == -1) {
        break;
      }
    }
  } while (k < maxK);

  assert(k > 0);
  ScalarArray<dp_t>* newA = new ScalarArray<dp_t>(block.rows->size(), k);
  for (int i = 0; i < k; i++) {
    newA->copyMatrixAtOffset(aCols[i], 0, i);
    delete aCols[i];
    aCols[i] = NULL;
  }
  ScalarArray<dp_t>* newB = new ScalarArray<dp_t>(block.cols->size(), k);
  for (int i = 0; i < k; i++) {
    newB->copyMatrixAtOffset(bCols[i], 0, i);
    delete bCols[i];
    bCols[i] = NULL;
  }
  return new RkMatrix<dp_t>(newA, block.rows, newB, block.cols);
}

RkMatrix<Types<S_t>::dp>*
CompressionAcaPlus::compress(const ClusterAssemblyFunction<S_t>& block) const {
    return doCompressionAcaPlus<S_t>(block, epsilon_);
}
RkMatrix<Types<D_t>::dp>*
CompressionAcaPlus::compress(const ClusterAssemblyFunction<D_t>& block) const {
    return doCompressionAcaPlus<D_t>(block, epsilon_);
}
RkMatrix<Types<C_t>::dp>*
CompressionAcaPlus::compress(const ClusterAssemblyFunction<C_t>& block) const {
    return doCompressionAcaPlus<C_t>(block, epsilon_);
}
RkMatrix<Types<Z_t>::dp>*
CompressionAcaPlus::compress(const ClusterAssemblyFunction<Z_t>& block) const {
    return doCompressionAcaPlus<Z_t>(block, epsilon_);
}

#include <iostream>

template<typename T>
RkMatrix<typename Types<T>::dp>* compress(
    const CompressionAlgorithm* method, const Function<T>& f,
    const ClusterData* rows, const ClusterData* cols,
    double epsilon, const AllocationObserver & ao) {
    typedef typename Types<T>::dp dp_t;
    ClusterAssemblyFunction<T> block(f, rows, cols, ao);
    int nloop=-1; // so we assemble only one strata
    if(block.info.number_of_strata > 1 && method->isIncremental()) {
        block.stratum = 0;
        // enable strata assembling for AcaPartial & AcaPlus only
        nloop = block.info.number_of_strata;
    }
    RkMatrix<dp_t>* rk = compressOneStratum(method, block);
    rk->truncate(epsilon);
    for(block.stratum = 1; block.stratum < nloop; block.stratum++) {
        assert(method->isIncremental());
        RkMatrix<dp_t>* stratumRk = compressOneStratum(method, block);
        if(stratumRk->rank() > 0) {
            // Pass a negative value to tell formattedAddParts to not call truncate()
            // FIXME: investigate why calling truncate from formattedAddParts or here
            //        gives different results
            rk->formattedAddParts(-epsilon, &Constants<dp_t>::pone, &stratumRk, 1);
            rk->truncate(epsilon);
        }
        delete stratumRk;
    }
    return rk;
}

template<typename T> RkMatrix<typename Types<T>::dp>* compressOneStratum(
    const CompressionAlgorithm* method, const ClusterAssemblyFunction<T> & block) {

  typedef typename Types<T>::dp dp_t;
  RkMatrix<dp_t>* rk = method->compress(block);

  if (HMatrix<T>::validateCompression) {
    FullMatrix<dp_t>* full = block.assemble();
    rk->checkNan();
    FullMatrix<dp_t>* rkFull = rk->eval();
    const double approxNorm = rkFull->norm();
    const double fullNorm = full->norm();

    // If I meet a NaN, I save & leave
    // TODO : improve this behaviour
    if (isnan(approxNorm)) {
      rkFull->toFile("Rk");
      full->toFile("Full");
      HMAT_ASSERT(false);
    }

    rkFull->axpy(Constants<T>::mone, full);
    double diffNorm = rkFull->norm();
    if (diffNorm > HMatrix<T>::validationErrorThreshold * fullNorm ) {
      std::cout << block.rows->description() << "x" << block.cols->description() << std::endl
           << std::scientific
           << "|M|  = " << fullNorm << std::endl
           << "|Rk| = " << approxNorm << std::endl
           << "|M - Rk| / |M| = " << diffNorm / fullNorm << std::endl
           << "Rank = " << rk->rank() << " / " << min(full->rows(), full->cols()) << std::endl << std::endl;

      if (HMatrix<T>::validationReRun) {
        // Call compression a 2nd time, for debugging with gdb the work of the compression algorithm...
        RkMatrix<dp_t>* rk_bis = NULL;

        rk_bis = method->compress(block);
        delete rk_bis ;
      }

      if (HMatrix<T>::validationDump) {
        std::string filename;
        std::ostringstream convert;   // stream used for the conversion
        convert << block.stratum << "_"<< block.rows->description() << "x" << block.cols->description();

        filename = "Rk_";
        filename += convert.str(); // set 'Result' to the contents of the stream
        delete rkFull;
        rkFull = rk->eval();
        rkFull->toFile(filename.c_str());
        filename = "Full_"+convert.str(); // set 'Result' to the contents of the stream
        full->toFile(filename.c_str());
      }
    }

    delete rkFull;
    delete full;
  }
  return rk;
}

// Declaration of the used templates
template RkMatrix<S_t>* truncatedSvd(FullMatrix<S_t>* m, double eps);
template RkMatrix<D_t>* truncatedSvd(FullMatrix<D_t>* m, double eps);
template RkMatrix<C_t>* truncatedSvd(FullMatrix<C_t>* m, double eps);
template RkMatrix<Z_t>* truncatedSvd(FullMatrix<Z_t>* m, double eps);

template RkMatrix<Types<S_t>::dp>* compress<S_t>(const CompressionAlgorithm* method, const Function<S_t>& f, const ClusterData* rows, const ClusterData* cols, double epsilon, const AllocationObserver &);
template RkMatrix<Types<D_t>::dp>* compress<D_t>(const CompressionAlgorithm* method, const Function<D_t>& f, const ClusterData* rows, const ClusterData* cols, double epsilon, const AllocationObserver &);
template RkMatrix<Types<C_t>::dp>* compress<C_t>(const CompressionAlgorithm* method, const Function<C_t>& f, const ClusterData* rows, const ClusterData* cols, double epsilon, const AllocationObserver &);
template RkMatrix<Types<Z_t>::dp>* compress<Z_t>(const CompressionAlgorithm* method, const Function<Z_t>& f, const ClusterData* rows, const ClusterData* cols, double epsilon, const AllocationObserver &);

}  // end namespace hmat

