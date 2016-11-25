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

#include "hmat_cpp_interface.hpp"
#include "h_matrix.hpp"
#include "rk_matrix.hpp"
#include "cluster_tree.hpp"
#include "common/context.hpp"
#include "disable_threading.hpp"
#include "json.hpp"
#include "iengine.hpp"
#include "serialization.hpp"
#include <cstring>
#include <fstream>

namespace hmat {

// HMatInterface
template<typename T>
HMatInterface<T>::HMatInterface(IEngine<T>* engine, const ClusterTree* _rows, const ClusterTree* _cols,
                                SymmetryFlag sym, AdmissibilityCondition * admissibilityCondition) :
  engine_(engine),factorizationType(hmat_factorization_none)
{
  DECLARE_CONTEXT;
  engine_->hmat = new HMatrix<T>(_rows, _cols, &HMatSettings::getInstance(), 0, sym, admissibilityCondition);
}

template<typename T>
HMatInterface<T>::~HMatInterface() {
  engine_->destroy();
  delete engine_->hmat;
  delete engine_;
}

template<typename T>
HMatInterface<T>::HMatInterface(IEngine<T>* engine, HMatrix<T>* h, hmat_factorization_t factorization):
  engine_(engine)
{
  engine_->setHMatrix(h);
      factorizationType = factorization;
}

static void io_write(void * buffer, size_t size, void * user_data) {
    size_t r = fwrite(buffer, size, 1, (FILE*)user_data);
    HMAT_ASSERT(r == 1);
}

template<typename T>
void HMatInterface<T>::assemble(Assembly<T>& f, SymmetryFlag sym, bool,
                                   hmat_progress_t * progress, bool ownAssembly) {
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  engine_->progress(progress);
  engine_->assembly(f, sym, ownAssembly);
  FILE * file = fopen("assembled.hmat", "wb");
  HMAT_ASSERT(file != NULL);
  MatrixStructMarshaller<T>(io_write, file).write(engine_.hmat, factorizationType);
  fclose(file);
  exit(0);
}

template<typename T>
void HMatInterface<T>::factorize(hmat_factorization_t t, hmat_progress_t * progress) {
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  engine_->progress(progress);
  if(progress != NULL)
    progress->max = engine_->hmat->rows()->size();
  engine_->factorization(t);
  factorizationType = t;
  engine_->hmat->checkStructure();
}

template<typename T>
void HMatInterface<T>::inverse(hmat_progress_t * progress) {
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  engine_->progress(progress);
  engine_->inverse();
}

template<typename T>
void HMatInterface<T>::gemv(char trans, T alpha, ScalarArray<T>& x, T beta,
                            ScalarArray<T>& y) const {
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  reorderVector(&x, trans == 'N' ? engine_->hmat->cols()->indices() : engine_->hmat->rows()->indices());
  reorderVector(&y, trans == 'N' ? engine_->hmat->rows()->indices() : engine_->hmat->cols()->indices());
  engine_->gemv(trans, alpha, x, beta, y);
  restoreVectorOrder(&x, trans == 'N' ? engine_->hmat->cols()->indices() : engine_->hmat->rows()->indices());
  restoreVectorOrder(&y, trans == 'N' ? engine_->hmat->rows()->indices() : engine_->hmat->cols()->indices());
}

template<typename T>
void HMatInterface<T>::gemm(char transA, char transB, T alpha,
                            const HMatInterface<T>* a,
                            const HMatInterface<T>* b, T beta) {
    DISABLE_THREADING_IN_BLOCK;
    DECLARE_CONTEXT;
    engine_->gemm(transA, transB, alpha, *a->engine_, *b->engine_, beta);
    engine_->hmat->checkStructure();
}

template<typename T>
void HMatInterface<T>::gemm(ScalarArray<T>& c, char transA, char transB, T alpha,
                            ScalarArray<T>& a, const HMatInterface<T>& b,
                            T beta) {
  DECLARE_CONTEXT;
  // C <- AB + C  <=>  C^t <- B^t A^t + C^t
  // On fait les operations dans ce sens pour etre dans le bon sens
  // pour la memoire, et pour reordonner correctement les "vecteurs" A
  // et C.
  if (transA == 'N') {
    a.transpose();
  }
  c.transpose();
  b.gemv(transB == 'N' ? 'T' : 'N', alpha, a, beta, c);
  c.transpose();
  if (transA == 'N') {
    a.transpose();
  }
}


template<typename T>
void HMatInterface<T>::solve(ScalarArray<T>& b) const {
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  reorderVector<T>(&b, engine_->hmat->cols()->indices());
  engine_->solve(b, factorizationType);
  restoreVectorOrder<T>(&b, engine_->hmat->cols()->indices());
}

template<typename T>
void HMatInterface<T>::solve(HMatInterface<T>& b) const {
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  engine_->solve(*b.engine_, factorizationType);
}

template<typename T>
void HMatInterface<T>::solveLower(ScalarArray<T>& b, bool transpose) const {
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  if (transpose)
    reorderVector<T>(&b, engine_->hmat->rows()->indices());
  else
    reorderVector<T>(&b, engine_->hmat->cols()->indices());
  engine_->solveLower(b, factorizationType, transpose);
  if (transpose)
    restoreVectorOrder<T>(&b, engine_->hmat->rows()->indices());
  else
    restoreVectorOrder<T>(&b, engine_->hmat->cols()->indices());
}

template<typename T>
HMatInterface<T>* HMatInterface<T>::copy(bool structOnly) const {
  DECLARE_CONTEXT;
  HMatInterface<T>* result = new HMatInterface<T>(engine_->clone(), NULL);
  engine_->copy(*(result->engine_), structOnly);
  assert(result->engine_->hmat);
  result->engine_->hmat->checkStructure();
  return result;
}

template<typename T>
void HMatInterface<T>::transpose() {
  DECLARE_CONTEXT;
  engine_->transpose();
  engine_->hmat->checkStructure();
}


template<typename T>
double HMatInterface<T>::norm() const {
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  return engine_->hmat->norm();
}

template<typename T>
void HMatInterface<T>::scale(T alpha) {
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  engine_->scale(alpha);
}

template<typename T>
void HMatInterface<T>::addIdentity(T alpha) {
  DECLARE_CONTEXT;
  engine_->hmat->addIdentity(alpha);
}

template<typename T>
void HMatInterface<T>::addRand(double epsilon) {
  DECLARE_CONTEXT;
  engine_->addRand(epsilon);
}

template<typename T>
void HMatInterface<T>::info(hmat_info_t & result) const {
  DECLARE_CONTEXT;
    memset(&result, 0, sizeof(hmat_info_t));
    engine_->info(result);
}

template<typename T>
void HMatInterface<T>::dumpTreeToFile(const std::string& filename) const {
  DECLARE_CONTEXT;
  std::ofstream out(filename.c_str());
  HMatrixJSONDumper<T>(engine_->hmat, out).dump();
}

template<typename T>
int HMatInterface<T>::nodesCount() const {
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  return engine_->hmat->nodesCount();
}

template<typename T>
void HMatInterface<T>::walk(TreeProcedure<HMatrix<T> > *proc){
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  return engine_->hmat->walk(proc);
}

template<typename T>
void HMatInterface<T>::apply_on_leaf(const LeafProcedure<HMatrix<T> >& proc){
  DISABLE_THREADING_IN_BLOCK;
  DECLARE_CONTEXT;
  engine_->applyOnLeaf(proc);
}

// Explicit template instantiation
template class HMatInterface<S_t>;
template class HMatInterface<D_t>;
template class HMatInterface<C_t>;
template class HMatInterface<Z_t>;


} // end namespace hmat


