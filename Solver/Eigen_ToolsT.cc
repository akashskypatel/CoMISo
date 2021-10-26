/*===========================================================================*\
 *                                                                           *
 *                               CoMISo                                      *
 *      Copyright (C) 2008-2009 by Computer Graphics Group, RWTH Aachen      *
 *                           www.rwth-graphics.de                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of CoMISo.                                             *
 *                                                                           *
 *  CoMISo is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  CoMISo is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with CoMISo.  If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                           *
\*===========================================================================*/



//=============================================================================
//
//  CLASS Eigen_Tools - IMPLEMENTATION
//
//=============================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_EIGEN3_AVAILABLE


#define COMISO_Eigen_TOOLS_C

//== INCLUDES =================================================================

#include "Eigen_Tools.hh"
#include <CoMISo/Utils/gmm.hh>
#include <CoMISo/Utils/CoMISoError.hh>
#include <unsupported/Eigen/SparseExtra>
#include <queue>

//== NAMESPACES ===============================================================

namespace COMISO_EIGEN
{

//== IMPLEMENTATION ==========================================================

//-----------------------------------------------------------------------------

template <class MatrixT, class REALT, class INTT>
void get_ccs_symmetric_data(const MatrixT& _mat, const char _uplo,
    std::vector<REALT>& _values, std::vector<INTT>& _rowind,
    std::vector<INTT>& _colptr)
{
  // Assumes col major

  int m = _mat.innerSize();
  int n = _mat.outerSize();

  _values.resize(0);
  _rowind.resize(0);
  _colptr.resize(0);

  INTT iv(0);

  typedef typename MatrixT::InnerIterator It;

  switch (_uplo)
  {
  case 'l':
  case 'L':
    // for all columns
    for (int i = 0; i < n; ++i)
    {
      _colptr.push_back(iv);

      // row it
      It it(_mat, i);

      for (; it; ++it)
      {
        if (it.index() >= i)
        {
          _values.push_back(it.value());
          _rowind.push_back(it.index());
          ++iv;
        }
      }
    }
    _colptr.push_back(iv);
    break;

  case 'u':
  case 'U':
    // for all columns
    for (int i = 0; i < n; ++i)
    {
      _colptr.push_back(iv);

      // row it
      It it(_mat, i);

      for (; it; ++it)
      {
        if (it.index() <= i)
        {
          _values.push_back(it.value());
          _rowind.push_back(it.index());
          ++iv;
        }
      }
    }
    _colptr.push_back(iv);
    break;

  case 'c':
  case 'C':
    // for all columns
    for (int i = 0; i < n; ++i)
    {
      _colptr.push_back(iv);

      // row it
      It it(_mat, i);

      for (; it; ++it)
      {
        _values.push_back(it.value());
        _rowind.push_back(it.index());
        ++iv;
      }
    }
    _colptr.push_back(iv);
    break;

  default:
    std::cerr << "ERROR: parameter uplo must bei either 'U' or 'L' or 'C'!!!\n";
    break;
  }
}

//-----------------------------------------------------------------------------

// inspect the matrix: dimension, symmetry, zero_rows, zero_cols, nnz, max, min,
// max_abs, min_abs, NAN, INF
template <class MatrixT> void inspect_matrix(const MatrixT& _A)
{

  std::cerr << "################### INSPECT MATRIX ##################\n";
  std::cerr << "#outer size  : " << _A.outerSize() << std::endl;
  std::cerr << "#inner size  : " << _A.innerSize() << std::endl;
  std::cerr << "#rows        : " << _A.rows() << std::endl;
  std::cerr << "#cols        : " << _A.cols() << std::endl;
  std::cerr << "#nonzeros    : " << _A.nonZeros() << std::endl;
  std::cerr << "#nonzeros/row: " << (double(_A.nonZeros()) / double(_A.rows()))
            << std::endl;
  std::cerr << "symmetric    : " << is_symmetric(_A) << std::endl;

  MatrixT trans(_A.transpose());

  int zero_rows = 0;
  int zero_cols = 0;

  for (int i = 0; i < _A.outerSize(); ++i)
  {
    typename MatrixT::InnerIterator it(_A, i);
    if (!it)
      ++zero_rows;
  }

  for (int i = 0; i < trans.outerSize(); ++i)
  {
    typename MatrixT::InnerIterator it(trans, i);
    if (!it)
      ++zero_cols;
  }

  std::cerr << "zero rows    : " << zero_rows << std::endl;
  std::cerr << "zero cols    : " << zero_cols << std::endl;

  typedef typename MatrixT::Scalar Scalar;
  Scalar vmin = std::numeric_limits<Scalar>::max();
  Scalar vmax = std::numeric_limits<Scalar>::min();
  Scalar vmin_abs = std::numeric_limits<Scalar>::max();
  Scalar vmax_abs = 0;

  int n_nan = 0;
  int n_inf = 0;

  // inspect elements
  for (int i = 0; i < _A.outerSize(); ++i)
  {
    typename MatrixT::InnerIterator it(_A, i);

    for (; it; ++it)
    {
      if (it.value() < vmin)
        vmin = it.value();
      if (it.value() > vmax)
        vmax = it.value();

      if (fabs(it.value()) < vmin_abs)
        vmin_abs = fabs(it.value());
      if (fabs(it.value()) > vmax_abs)
        vmax_abs = fabs(it.value());

      if (std::isnan(it.value()))
        ++n_nan;
      if (std::isinf(it.value()))
        ++n_inf;
    }
  }

  std::cerr << "min  val     : " << vmin << std::endl;
  std::cerr << "max  val     : " << vmax << std::endl;
  std::cerr << "min |val|    : " << vmin_abs << std::endl;
  std::cerr << "max |val|    : " << vmax_abs << std::endl;
  std::cerr << "#nan         : " << n_nan << std::endl;
  std::cerr << "#inf         : " << n_inf << std::endl;

  std::cerr << "min eval     : "
            << "..." << std::endl;
  std::cerr << "max eval     : "
            << "..." << std::endl;
  std::cerr << "min|eval|    : "
            << "..." << std::endl;
  std::cerr << "max|eval|    : "
            << "..." << std::endl;
}

//-----------------------------------------------------------------------------

// symmetric ?
template <class MatrixT> bool is_symmetric(const MatrixT& _A)
{
  typedef typename MatrixT::InnerIterator It;
  typedef typename MatrixT::Scalar Scalar;

  int nouter(_A.outerSize());
  int ninner(_A.innerSize());

  if (nouter != ninner)
    return false;

  bool symmetric(true);

  for (int c = 0; c < nouter; ++c)
  {
    for (It it(_A, c); it; ++it)
    {
      int r(it.index());

      Scalar val(it.value());

      // find diagonal partner element
      bool found(false);
      for (It dit(_A, r); dit; ++dit)
      {
        if (dit.index() < c)
        {
        }
        else if (dit.index() == c)
        {
          if (dit.value() == val)
            found = true;
          break;
        }
        else
        {
          break;
        }
      }
      if (!found)
      {
        symmetric = false;
        break;
      }
    }
  }
  return symmetric;
}

//-----------------------------------------------------------------------------

template <class Eigen_MatrixT, class IntT>
void permute(
    const Eigen_MatrixT& _QR, const std::vector<IntT>& _Pvec, Eigen_MatrixT& _A)
{
  typedef typename Eigen_MatrixT::Scalar Scalar;

  int m = _QR.innerSize();
  int n = _QR.outerSize();

  if (_Pvec.size() == 0)
  {
    _A = _QR;
    return;
  }

  if (_Pvec.size() != (size_t)_QR.cols() && _Pvec.size() != 0)
  {
    std::cerr << __FUNCTION__
              << " wrong size of permutation vector, should have #cols length "
                 "(or zero)"
              << std::endl;
  }

  // build sparse permutation matrix
  typedef Eigen::Triplet<Scalar> Triplet;
  std::vector<Triplet> triplets;
  triplets.reserve(_QR.nonZeros());
  _A = Eigen_MatrixT(m, n);

  typedef typename Eigen_MatrixT::InnerIterator It;

  for (int c = 0; c < n; ++c) // cols
  {
    for (It it(_QR, c); it; ++it) // rows
    {
      int r(it.index());

      Scalar val(it.value());

      int newcol(_Pvec[c]);

      triplets.push_back(Triplet(r, newcol, val));
    }
  }
  _A.setFromTriplets(triplets.begin(), triplets.end());
}

//-----------------------------------------------------------------------------

#if COMISO_SUITESPARSE_AVAILABLE

/// Eigen to Cholmod_sparse interface
template <class MatrixT>
void cholmod_to_eigen(const cholmod_sparse& _AC, MatrixT& _A)
{
  // initialize dimensions
  typedef typename MatrixT::Scalar Scalar;
  typedef Eigen::Triplet<Scalar> Triplet;
  size_t nzmax(_AC.nzmax);
  std::cerr << __FUNCTION__ << " row " << _AC.nrow << " col " << _AC.ncol
            << " stype " << _AC.stype << std::endl;
  _A = MatrixT((int)_AC.nrow, (int)_AC.ncol);
  std::vector<Triplet> triplets;
  triplets.reserve(nzmax);

  if (!_AC.packed)
  {
    std::cerr << "Warning: " << __FUNCTION__
              << " does not support unpacked matrices yet!!!" << std::endl;
    return;
  }

  // Pointer to data
  double* X((double*)_AC.x);

  // complete matrix stored
  if (_AC.stype == 0)
  {
    // which integer type?
    if (_AC.itype == CHOLMOD_LONG)
    {
      SuiteSparse_long* P((SuiteSparse_long*)_AC.p);
      SuiteSparse_long* I((SuiteSparse_long*)_AC.i);

      for (SuiteSparse_long i = 0; i < (SuiteSparse_long)_AC.ncol; ++i)
        for (SuiteSparse_long j = P[i]; j < P[i + 1]; ++j)
          //_A( I[j], i) += X[j]; // += really needed?
          triplets.push_back(Triplet((int)I[j], (int)i, X[j]));
    }
    else
    {
      int* P((int*)_AC.p);
      int* I((int*)_AC.i);

      for (int i = 0; i < (int)_AC.ncol; ++i)
        for (int j = P[i]; j < P[i + 1]; ++j)
          triplets.push_back(Triplet((int)I[j], (int)i, X[j]));
      //_A( I[j], i) += X[j];
    }
  }
  else // only upper or lower diagonal stored
  {
    // which integer type?
    if (_AC.itype == CHOLMOD_LONG)
    {
      SuiteSparse_long* P((SuiteSparse_long*)_AC.p);
      SuiteSparse_long* I((SuiteSparse_long*)_AC.i);

      for (SuiteSparse_long i = 0; i < (SuiteSparse_long)_AC.ncol; ++i)
        for (SuiteSparse_long j = P[i]; j < P[i + 1]; ++j)
        {
          //_A(I[j], i) += X[j];
          triplets.push_back(Triplet((int)I[j], (int)i, X[j]));

          // add up symmetric part
          if (I[j] != i)
            triplets.push_back(Triplet((int)i, (int)I[j], X[j]));
          //_A(i,I[j]) += X[j];
        }
    }
    else
    {
      int* P((int*)_AC.p);
      int* I((int*)_AC.i);

      for (int i = 0; i < (int)_AC.ncol; ++i)
        for (int j = P[i]; j < P[i + 1]; ++j)
        {
          //_A(I[j], i) += X[j];
          triplets.push_back(Triplet(I[j], i, X[j]));

          // add up symmetric part
          if (I[j] != i)
            //  _A(i,I[j]) += X[j];
            triplets.push_back(Triplet(i, I[j], X[j]));
        }
    }
  }
  _A.setFromTriplets(triplets.begin(), triplets.end());
}

/// GMM to Cholmod_sparse interface
template <class MatrixT>
void eigen_to_cholmod(const MatrixT& _A, cholmod_sparse*& _AC,
    cholmod_common* _common, int _sparsity_type, bool _long_int)
{
  /* _sparsity_type
   * 0:  matrix is "unsymmetric": use both upper and lower triangular parts
   *     (the matrix may actually be symmetric in pattern and value, but
   *     both parts are explicitly stored and used).  May be square or
   *     rectangular.
   * >0: matrix is square and symmetric, use upper triangular part.
   *     Entries in the lower triangular part are ignored.
   * <0: matrix is square and symmetric, use lower triangular part.
   *     Entries in the upper triangular part are ignored. */

  int m = _A.innerSize();
  int n = _A.outerSize();

  // get upper or lower
  char uplo = 'c';
  if (_sparsity_type < 0)
    uplo = 'l';
  if (_sparsity_type > 0)
    uplo = 'u';

  if (_long_int) // long int version
  {
    std::vector<double> values;
    std::vector<SuiteSparse_long> rowind;
    std::vector<SuiteSparse_long> colptr;

    // get data of gmm matrix
    COMISO_EIGEN::get_ccs_symmetric_data(_A, uplo, values, rowind, colptr);

    // allocate cholmod matrix
    _AC = cholmod_l_allocate_sparse(
        m, n, values.size(), true, true, _sparsity_type, CHOLMOD_REAL, _common);

    // copy data to cholmod matrix
    for (SuiteSparse_long i = 0; i < (SuiteSparse_long)values.size(); ++i)
    {
      ((double*)(_AC->x))[i] = values[i];
      ((SuiteSparse_long*)(_AC->i))[i] = rowind[i];
    }

    for (SuiteSparse_long i = 0; i < (SuiteSparse_long)colptr.size(); ++i)
      ((SuiteSparse_long*)(_AC->p))[i] = colptr[i];
  }
  else // int version
  {
    std::vector<double> values;
    std::vector<int> rowind;
    std::vector<int> colptr;

    // get data of gmm matrix
    COMISO_EIGEN::get_ccs_symmetric_data(_A, uplo, values, rowind, colptr);

    // allocate cholmod matrix
    _AC = cholmod_allocate_sparse(
        m, n, values.size(), true, true, _sparsity_type, CHOLMOD_REAL, _common);

    // copy data to cholmod matrix
    for (unsigned int i = 0; i < values.size(); ++i)
    {
      ((double*)(_AC->x))[i] = values[i];
      ((int*)(_AC->i))[i] = rowind[i];
    }
    for (unsigned int i = 0; i < colptr.size(); ++i)
      ((int*)(_AC->p))[i] = colptr[i];
  }
}
#endif

/*
/// Eigen to Cholmod_dense interface
template<class MatrixT>
void cholmod_to_eigen_dense( const cholmod_dense& _AC, MatrixT& _A)
{
  // initialize dimensions
  typedef typename MatrixT::Scalar Scalar;
  typedef Eigen::Triplet< Scalar > Triplet;
  size_t nzmax( _AC.nzmax);
  _A = MatrixT(_AC.nrow, _AC.ncol);
  std::vector< Triplet > triplets;
  triplets.reserve(nzmax);

  if(!_AC.packed)
  {
    std::cerr << "Warning: " << __FUNCTION__ << " does not support unpacked
    matrices yet!!!" << std::endl; return;
  }

  // Pointer to data
  double* X((double*)_AC.x);

  // complete matrix stored
  if(_AC.stype == 0)
  {
    // which integer type?
    if(_AC.itype == CHOLMOD_LONG)
    {
      SuiteSparse_long* P((SuiteSparse_long*)_AC.p);
      SuiteSparse_long* I((SuiteSparse_long*)_AC.i);

      for(SuiteSparse_long i=0; i<(SuiteSparse_long)_AC.ncol; ++i)
        for(SuiteSparse_long j= P[i]; j< P[i+1]; ++j)
          //_A( I[j], i) += X[j]; // += really needed?
          triplets.push_back( Triplet( I[j], i, X[j]));
    }
    else
    {
      int* P((int*)_AC.p);
      int* I((int*)_AC.i);

      for(int i=0; i<(int)_AC.ncol; ++i)
        for(int j= P[i]; j< P[i+1]; ++j)
          triplets.push_back( Triplet( I[j], i, X[j]));
      //_A( I[j], i) += X[j];
    }

  }
  else // only upper or lower diagonal stored
  {
    // which integer type?
    if(_AC.itype == CHOLMOD_LONG)
    {
      SuiteSparse_long* P((SuiteSparse_long*)_AC.p);
      SuiteSparse_long* I((SuiteSparse_long*)_AC.i);

      for(SuiteSparse_long i=0; i<(SuiteSparse_long)_AC.ncol; ++i)
        for(SuiteSparse_long j=P[i]; j<P[i+1]; ++j)
        {
          //_A(I[j], i) += X[j];
          triplets.push_back( Triplet( I[j], i, X[j]));

          // add up symmetric part
          if( I[j] != i)
            triplets.push_back( Triplet( i, I[j], X[j]));
          //_A(i,I[j]) += X[j];
        }
    }
    else
    {
      int* P((int*)_AC.p);
      int* I((int*)_AC.i);

      for(int i=0; i<(int)_AC.ncol; ++i)
        for(int j=P[i]; j<P[i+1]; ++j)
        {
          //_A(I[j], i) += X[j];
          triplets.push_back( Triplet( I[j], i, X[j]));

          // add up symmetric part
          if( I[j] != i)
            //  _A(i,I[j]) += X[j];
            triplets.push_back( Triplet( i, I[j], X[j]));
        }
    }
  }
  _A.setFromTriplets( triplets.begin(), triplets.end());
}


/// GMM to Cholmod_sparse interface
template<class MatrixT>
void eigen_to_cholmod_dense( const MatrixT& _A, cholmod_dense* &_AC,
cholmod_common* _common, bool _long_int)
{
  int m = _A.innerSize();
  int n = _A.outerSize();

  // allocate cholmod matrix
  _AC = cholmod_l_allocate_sparse(m,n,values.size(),true,true,_sparsity_type,
      CHOLMOD_REAL, _common);
  _AC = cholmod_l_allocate_dense (m,n,n, xtype, cc)

  if( _long_int) // long int version
  {
    std::vector<double> values;
    std::vector<SuiteSparse_long> rowind;
    std::vector<SuiteSparse_long> colptr;

    // get data of gmm matrix
    COMISO_EIGEN::get_ccs_symmetric_data( _A, uplo, values, rowind, colptr);

    // allocate cholmod matrix
    _AC = cholmod_l_allocate_sparse(m,n,values.size(),true,true,_sparsity_type,
        CHOLMOD_REAL,_common);

    // copy data to cholmod matrix
    for(SuiteSparse_long i=0; i<(SuiteSparse_long)values.size(); ++i)
    {
      ((double*) (_AC->x))[i] = values[i];
      ((SuiteSparse_long*)(_AC->i))[i] = rowind[i];
    }

    for(SuiteSparse_long i=0; i<(SuiteSparse_long)colptr.size(); ++i)
      ((SuiteSparse_long*)(_AC->p))[i] = colptr[i];
  }
  else // int version
  {
     std::vector<double> values;
     std::vector<int> rowind;
     std::vector<int> colptr;

     // get data of gmm matrix
     COMISO_EIGEN::get_ccs_symmetric_data( _A, uplo, values, rowind, colptr);

     // allocate cholmod matrix
     _AC = cholmod_allocate_sparse(m,n,values.size(),true,true,_sparsity_type,
         CHOLMOD_REAL, _common);

     // copy data to cholmod matrix
     for(unsigned int i=0; i<values.size(); ++i)
     {
       ((double*)(_AC->x))[i] = values[i];
       ((int*)   (_AC->i))[i] = rowind[i];
     }
     for(unsigned int i=0; i<colptr.size(); ++i)
       ((int*)(_AC->p))[i] = colptr[i];
  }

}*/

/*
// convert a gmm col-sparse matrix into an eigen sparse matrix
template<class GMM_MatrixT, class EIGEN_MatrixT>
void gmm_to_eigen( const GMM_MatrixT& _G, EIGEN_MatrixT& _E)
{
  typedef typename EIGEN_MatrixT::Scalar Scalar;

  typedef typename gmm::linalg_traits<GMM_MatrixT>::const_sub_col_type ColT;
  typedef typename gmm::linalg_traits<ColT>::const_iterator CIter;

  // build matrix triplets
  typedef Eigen::Triplet< Scalar > Triplet;
  std::vector< Triplet > triplets;
  triplets.reserve(gmm::nnz(_G));

  for(unsigned int i=0; i<gmm::mat_ncols(_G); ++i)
  {
     ColT col = mat_const_col( _G, i );

     CIter it  = gmm::vect_const_begin( col );
     CIter ite = gmm::vect_const_end( col );
     for ( ; it!=ite; ++it )
       triplets.push_back( Triplet( it.index(), i, *it));

  }

  // generate eigen matrix
  _E = EIGEN_MatrixT( gmm::mat_nrows(_G), gmm::mat_ncols(_G));
  _E.setFromTriplets( triplets.begin(), triplets.end());
}
*/

/*!
This unusual approach to implement the different partial specializations for the
template gmm_to_eigen() function was forced by the MS VC2015 migration: The VC
linker is unable to match the explicitly instantiated template partial
specializations. In previous versions (VS2012, VS2013) we were able to work
around this by partially specializing the function for the template arguments
we needed. This no longer works.

Using a separate function name (gmm_to_eigen_impl::f()) for the partial
specializations we work around the MSVC linker limitations. Now the interface
function gmm_to_eigen no longer has any partial specializations, instead the
partial specializations are confined to gmm_to_eigen_impl::f().
*/
namespace gmm_to_eigen_impl
{

// convert a gmm col-sparse matrix into an eigen sparse matrix
template <class GMM_VectorT, class EIGEN_MatrixT>
void f(const gmm::col_matrix<GMM_VectorT>& _G, EIGEN_MatrixT& _E)
{
  typedef typename EIGEN_MatrixT::Scalar Scalar;

  typedef typename gmm::col_matrix<GMM_VectorT> GMM_MatrixT;
  typedef typename gmm::linalg_traits<GMM_MatrixT>::const_sub_col_type ColT;
  typedef typename gmm::linalg_traits<ColT>::const_iterator CIter;

  // build matrix triplets
  typedef Eigen::Triplet<Scalar> Triplet;
  std::vector<Triplet> triplets;
  triplets.reserve(gmm::nnz(_G));

  for (unsigned int i = 0; i < gmm::mat_ncols(_G); ++i)
  {
    ColT col = mat_const_col(_G, i);

    auto it = gmm::vect_const_begin(col);
    auto ite = gmm::vect_const_end(col);
    for (; it != ite; ++it)
    {
      triplets.push_back(
          Triplet(static_cast<int>(it.index()), static_cast<int>(i), *it));
    }
  }

  // generate eigen matrix
  _E = EIGEN_MatrixT(gmm::mat_nrows(_G), gmm::mat_ncols(_G));
  _E.setFromTriplets(triplets.begin(), triplets.end());
}

// convert a gmm row-sparse matrix into an eigen sparse matrix
template <class GMM_VectorT, class EIGEN_MatrixT>
void f(const gmm::row_matrix<GMM_VectorT>& _G, EIGEN_MatrixT& _E)
{
  typedef typename EIGEN_MatrixT::Scalar Scalar;

  typedef typename gmm::row_matrix<GMM_VectorT> GMM_MatrixT;
  typedef typename gmm::linalg_traits<GMM_MatrixT>::const_sub_row_type RowT;
  typedef typename gmm::linalg_traits<RowT>::const_iterator CIter;

  // build matrix triplets
  typedef Eigen::Triplet<Scalar> Triplet;
  std::vector<Triplet> triplets;
  triplets.reserve(gmm::nnz(_G));

  for (gmm::size_type i = 0; i < gmm::mat_nrows(_G); ++i)
  {
    RowT row = mat_const_row(_G, i);

    auto it = gmm::vect_const_begin(row);
    auto ite = gmm::vect_const_end(row);
    for (; it != ite; ++it)
      triplets.push_back(
          Triplet(static_cast<int>(i), static_cast<int>(it.index()), *it));
  }

  // generate eigen matrix
  _E = EIGEN_MatrixT(static_cast<int>(gmm::mat_nrows(_G)),
      static_cast<int>(gmm::mat_ncols(_G)));
  _E.setFromTriplets(triplets.begin(), triplets.end());
}

// convert a gmm col-sparse matrix into an eigen sparse matrix
template <class GMM_RealT, class EIGEN_MatrixT>
void f(const gmm::csc_matrix<GMM_RealT>& _G, EIGEN_MatrixT& _E)
{
  typedef typename EIGEN_MatrixT::Scalar Scalar;

  typedef typename gmm::csc_matrix<GMM_RealT> GMM_MatrixT;
  typedef typename gmm::linalg_traits<GMM_MatrixT>::const_sub_col_type ColT;
  typedef typename gmm::linalg_traits<ColT>::const_iterator CIter;

  // build matrix triplets
  typedef Eigen::Triplet<Scalar> Triplet;
  std::vector<Triplet> triplets;
  triplets.reserve(gmm::nnz(_G));

  for (unsigned int i = 0; i < gmm::mat_ncols(_G); ++i)
  {
    ColT col = mat_const_col(_G, i);

    CIter it = gmm::vect_const_begin(col);
    CIter ite = gmm::vect_const_end(col);
    for (; it != ite; ++it)
      triplets.push_back(Triplet(static_cast<int>(it.index()), i, *it));
  }

  // generate eigen matrix
  _E = EIGEN_MatrixT(static_cast<int>(gmm::mat_nrows(_G)),
      static_cast<int>(gmm::mat_ncols(_G)));
  _E.setFromTriplets(triplets.begin(), triplets.end());
}

} // namespace gmm_to_eigen_impl

template <class GMM_MatrixT, class EIGEN_MatrixT>
void gmm_to_eigen(const GMM_MatrixT& _G, EIGEN_MatrixT& _E)
{
  gmm_to_eigen_impl::f(_G, _E);
}

template <typename MatrixT>
void write_eigen_matrix_ascii(const std::string& _filename, const MatrixT& _m)
{
  Eigen::saveMarket(_m, _filename);
}

template <typename MatrixT>
void read_eigen_matrix_ascii(const std::string& _filename, MatrixT& _m)
{
  Eigen::SparseMatrix<double> m;
  Eigen::loadMarket(m, _filename);
  _m = m;
}

template <typename ScalarT, int OPTIONS, typename StorageIndexT>
void read_eigen_matrix_ascii(const std::string& _filename,
    Eigen::SparseMatrix<ScalarT, OPTIONS, StorageIndexT>& _m)
{
  Eigen::loadMarket(_m, _filename);
}

template <typename VectorT>
void write_eigen_vector_ascii(const std::string& _filename, const VectorT& _v)
{
  Eigen::saveMarketVector(_v, _filename);
}

template <typename VectorT>
void read_eigen_vector_ascii(const std::string& _filename, VectorT& _v)
{
  Eigen::loadMarketVector(_v, _filename);
}

namespace detail
{

using TypeIndex = char;
template <typename T> struct ScalarTypeIndexT;
template <> struct ScalarTypeIndexT<double> { enum : TypeIndex { VALUE = 1 }; };
template <> struct ScalarTypeIndexT<float>  { enum : TypeIndex { VALUE = 2 }; };
template <> struct ScalarTypeIndexT<int>    { enum : TypeIndex { VALUE = 3 }; };

// Write a simple type to a binary stream
template <typename T> void write_simple(std::ostream& _os, const T& _t)
{
  _os.write((char*)&_t, sizeof(_t));
}

// File format version, in case we ever want to update the format
using FileFormatType = int;
static constexpr FileFormatType SPARSE_MATRIX_FILE_FORMAT_VERSION = 1;
static constexpr FileFormatType DENSE_MATRIX_FILE_FORMAT_VERSION = 1;

using DensityId = char;
static constexpr DensityId SPARSE_ID = 'S'; // used to indicate a sparse matrix is written
static constexpr DensityId DENSE_ID = 'D'; // used to indicate a dense matrix is written

} // namespace

template <typename ScalarT, int OPTIONS, typename StorageIndexT>
void write_eigen_matrix(const std::string& _filename,
    const Eigen::SparseMatrix<ScalarT, OPTIONS, StorageIndexT>& _m)
{
  using Matrix = Eigen::SparseMatrix<ScalarT, OPTIONS, StorageIndexT>;

  Matrix m_compressed;
  if (!_m.isCompressed())
  {
    m_compressed = _m;
    m_compressed.makeCompressed();
  }

  const Matrix& m = _m.isCompressed() ? _m : m_compressed;

  // open file
  std::ofstream out_file(_filename, std::ios::binary);
  COMISO_THROW_TODO_if(
      !out_file.is_open(), "Could not open file: " << _filename);

  // write that we are writing a sparse matrix
  detail::write_simple(out_file, detail::SPARSE_ID);

  // write file format version
  detail::write_simple(out_file, detail::SPARSE_MATRIX_FILE_FORMAT_VERSION);

  // write template info
  detail::write_simple(out_file, detail::ScalarTypeIndexT<ScalarT>::VALUE);
  detail::write_simple(out_file, OPTIONS);
  detail::write_simple(out_file, detail::ScalarTypeIndexT<StorageIndexT>::VALUE);

  // write matrix dimensions
  detail::write_simple(out_file, m.outerSize());
  detail::write_simple(out_file, m.innerSize());
  detail::write_simple(out_file, m.nonZeros());

  // write matrix data
  auto value_ptr = m.valuePtr();
  auto inner_ptr = m.innerIndexPtr();
  auto outer_ptr = m.outerIndexPtr();
  out_file.write((char*)value_ptr, m.nonZeros() * sizeof(ScalarT));
  out_file.write((char*)inner_ptr, m.nonZeros() * sizeof(StorageIndexT));
  out_file.write((char*)outer_ptr, (m.outerSize()+1) * sizeof(StorageIndexT));
}

namespace detail
{

// read a simple type from binary stream
template <typename T> void read_simple(std::istream& _is, T& _t)
{
  _is.read((char*)&_t, sizeof(T));
}

// Actual method that loads a matrix from a file
template <typename CallScalarT, int CALL_OPTIONS, typename CallStorageIndexT,
          typename FileScalarT, int FILE_OPTIONS, typename FileStorageIndexT>
void read_eigen_matrix_storage(std::istream& _is,
    Eigen::SparseMatrix<CallScalarT, CALL_OPTIONS, CallStorageIndexT>& _m)
{
  using Matrix = Eigen::SparseMatrix<FileScalarT, FILE_OPTIONS, FileStorageIndexT>;
  // read matrix dimensions
  typename Matrix::Index outerSize;
  typename Matrix::Index innerSize;
  typename Matrix::Index nonZeros;
  read_simple(_is, outerSize);
  read_simple(_is, innerSize);
  read_simple(_is, nonZeros);

  // read matrix data
  std::vector<FileScalarT> values(nonZeros);
  std::vector<FileStorageIndexT> inner_indices(nonZeros);
  std::vector<FileStorageIndexT> outer_indices(outerSize+1);
  _is.read((char*)values.data(), values.size() * sizeof(FileScalarT));
  _is.read((char*)inner_indices.data(),
      inner_indices.size() * sizeof(FileStorageIndexT));
  _is.read((char*)outer_indices.data(),
      outer_indices.size() * sizeof(FileStorageIndexT));

  // create triplets from matrix data
  std::vector<Eigen::Triplet<CallScalarT, CallStorageIndexT>> trips;
  size_t outer_id = 0;
  for (size_t i = 0; i < values.size(); ++i)
  {
    while (i >= outer_indices[outer_id + 1])
      ++outer_id;
    auto val = static_cast<CallScalarT>(values[i]);
    auto col = static_cast<CallStorageIndexT>(inner_indices[i]);
    auto row = static_cast<CallStorageIndexT>(outer_id);
    auto COL_MAJOR = Eigen::ColMajor;
    if (FILE_OPTIONS == COL_MAJOR) // TODO: if constexpr
      std::swap(col, row);
    trips.emplace_back(row, col, val);
  }

  // initialize _m from triplets
  auto rows = static_cast<CallStorageIndexT>(outerSize);
  auto cols = static_cast<CallStorageIndexT>(innerSize);
  auto COL_MAJOR = Eigen::ColMajor;
  if (FILE_OPTIONS == COL_MAJOR) // TODO: if constexpr
    std::swap(cols, rows);
  _m.resize(rows, cols);
  _m.setFromTriplets(trips.begin(), trips.end());
}

// Helper method to get correct storage index type from EigenTypeID argument
template <typename CallScalarT, int CALL_OPTIONS, typename CallStorageIndexT,
          typename FileScalarT, int FILE_OPTIONS>
void read_eigen_matrix_options(std::istream& _is,
    Eigen::SparseMatrix<CallScalarT, CALL_OPTIONS, CallStorageIndexT>& _m,
    TypeIndex _storage_index_id)
{
  // Dispatch to load method with correct storage index type
  switch (_storage_index_id)
  {
  case ScalarTypeIndexT<int>::VALUE:
    read_eigen_matrix_storage<CallScalarT, CALL_OPTIONS, CallStorageIndexT,
        FileScalarT, FILE_OPTIONS, int>(_is, _m);
    break;
  default:
    COMISO_THROW_TODO("Undexpected storage id: " << _storage_index_id);
  }
}


// Helper method to get correct ordering type from int argument
template <typename CallScalarT, int CALL_OPTIONS, typename CallStorageIndexT,
          typename FileScalarT>
void read_eigen_matrix_scalar(std::istream& _is,
    Eigen::SparseMatrix<CallScalarT, CALL_OPTIONS, CallStorageIndexT>& _m,
    int _options, TypeIndex _storage_index_is)
{
  // dispatch to load method with correct storage ordering
  if (_options == Eigen::RowMajor)
  {
    read_eigen_matrix_options<CallScalarT, CALL_OPTIONS, CallStorageIndexT,
        FileScalarT, Eigen::RowMajor>(_is, _m, _storage_index_is);
  }
  else if (_options == Eigen::ColMajor)
  {
    read_eigen_matrix_options<CallScalarT, CALL_OPTIONS, CallStorageIndexT,
        FileScalarT, Eigen::ColMajor>(_is, _m, _storage_index_is);
  }
  else
    COMISO_THROW_TODO("Read unexpected options: " << _options);
}

}

template <typename ScalarT, int OPTIONS, typename StorageIndexT>
void read_eigen_matrix(const std::string& _filename,
    Eigen::SparseMatrix<ScalarT, OPTIONS, StorageIndexT>& _m)
{
  // open file
  std::ifstream in_file(_filename, std::ios::binary);
  COMISO_THROW_TODO_if(
      !in_file.is_open(), "Could not open file: " << _filename);

  // Read whether the file actually contains a sparse matrix
  detail::DensityId density_id;
  detail::read_simple(in_file, density_id);
  COMISO_THROW_TODO_if(density_id != detail::SPARSE_ID,
      "Trying to read a non-sparse matrix into a sparse matrix");

  // Read file format version (ignored so far)
  detail::FileFormatType ff_version;
  detail::read_simple(in_file, ff_version);

  // read template type info
  detail::TypeIndex sclr_id;
  int options;
  detail::TypeIndex strg_index_id;
  detail::read_simple(in_file, sclr_id);
  detail::read_simple(in_file, options);
  detail::read_simple(in_file, strg_index_id);

  // dispatch to load method with correct scalar type
  switch (sclr_id)
  {
  case detail::ScalarTypeIndexT<double>::VALUE:
    detail::read_eigen_matrix_scalar<ScalarT, OPTIONS, StorageIndexT, double>(
        in_file, _m, options, strg_index_id);
    break;
  case detail::ScalarTypeIndexT<float>::VALUE:
    detail::read_eigen_matrix_scalar<ScalarT, OPTIONS, StorageIndexT, float>(
        in_file, _m, options, strg_index_id);
    break;
  case detail::ScalarTypeIndexT<int>::VALUE:
    detail::read_eigen_matrix_scalar<ScalarT, OPTIONS, StorageIndexT, int>(
        in_file, _m, options, strg_index_id);
    break;
  default:
    COMISO_THROW_TODO("Unexpected scalar id: " << sclr_id);
  }
}



// Write a dense matrix in our custom file format
template <typename ScalarT, int ROWS, int COLS>
void write_eigen_matrix(
    const std::string& _filename, const Eigen::Matrix<ScalarT, ROWS, COLS>& _m)
{
  using Matrix = Eigen::Matrix<ScalarT, ROWS, COLS>;

  // open file
  std::ofstream out_file(_filename, std::ios::binary);
  COMISO_THROW_TODO_if(
      !out_file.is_open(), "Could not open file: " << _filename);

  // write that we are writing a sparse matrix
  detail::write_simple(out_file, detail::DENSE_ID);

  // write file format version
  detail::write_simple(out_file, detail::DENSE_MATRIX_FILE_FORMAT_VERSION);

  // write template info
  detail::write_simple(out_file, detail::ScalarTypeIndexT<ScalarT>::VALUE);

  // write matrix dimensions
  detail::write_simple(out_file, _m.rows());
  detail::write_simple(out_file, _m.cols());

  // write matrix data (always row major)
  for (typename Matrix::Index row = 0; row < _m.rows(); ++row)
  {
    for (typename Matrix::Index col = 0; col < _m.cols(); ++col)
    {
      ScalarT val = _m(row, col);
      detail::write_simple(out_file, val);
    }
  }
}


namespace detail
{

template <typename CallScalarT, int ROWS, int COLS,
    typename FileScalarT>
void read_eigen_matrix_scalar(
    std::istream& _is, Eigen::Matrix<CallScalarT, ROWS, COLS>& _m)
{
  using Matrix = Eigen::Matrix<FileScalarT, Eigen::Dynamic, Eigen::Dynamic>;
  // read matrix dimensions
  typename Matrix::Index rows;
  typename Matrix::Index cols;
  detail::read_simple(_is, rows);
  detail::read_simple(_is, cols);

  int D = Eigen::Dynamic; // Needed to get rid of "Warning: if expr is const"

  COMISO_THROW_TODO_if(ROWS != D && ROWS != rows,
      "Dimension missmatch. Call asks for "
          << ROWS << " rows, but file stores a matrix with "
          << (int)rows << " rows.");

  COMISO_THROW_TODO_if(COLS != D && COLS != cols,
      "Dimension missmatch. Call asks for "
          << COLS << " cols, but file stores a matrix with "
          << (int)cols << " cols.");

  _m.resize(rows, cols);

  // read matrix data (always row major)
  for (typename Matrix::Index row = 0; row < _m.rows(); ++row)
  {
    for (typename Matrix::Index col = 0; col < _m.cols(); ++col)
    {
      FileScalarT val;
      detail::read_simple(_is, val);
      _m(row, col) = static_cast<CallScalarT>(val);
    }
  }

}

} // namespace detail

// Load a dense matrix in our custom file format
template <typename ScalarT, int ROWS, int COLS>
void read_eigen_matrix(
    const std::string& _filename, Eigen::Matrix<ScalarT, ROWS, COLS>& _m)
{
  // open file
  std::ifstream in_file(_filename, std::ios::binary);
  COMISO_THROW_TODO_if(
      !in_file.is_open(), "Could not open file: " << _filename);

  // Read whether the file actually contains a sparse matrix
  detail::DensityId density_id;
  detail::read_simple(in_file, density_id);
  COMISO_THROW_TODO_if(density_id != detail::DENSE_ID,
      "Trying to read a non-dense matrix into a dense matrix");

  // Read file format version (ignored so far)
  detail::FileFormatType ff_version;
  detail::read_simple(in_file, ff_version);

  // read template type info
  detail::TypeIndex sclr_id;
  detail::read_simple(in_file, sclr_id);

  // dispatch to load method with correct scalar type
  switch (sclr_id)
  {
  case detail::ScalarTypeIndexT<double>::VALUE:
    detail::read_eigen_matrix_scalar<ScalarT, ROWS, COLS, double>(in_file, _m);
    break;
  case detail::ScalarTypeIndexT<float>::VALUE:
    detail::read_eigen_matrix_scalar<ScalarT, ROWS, COLS, float>(in_file, _m);
    break;
  case detail::ScalarTypeIndexT<int>::VALUE:
    detail::read_eigen_matrix_scalar<ScalarT, ROWS, COLS, int>(in_file, _m);
    break;
  default:
    COMISO_THROW_TODO("Unexpected scalar id: " << sclr_id);
  }
}
//=============================================================================
} // namespace COMISO_EIGEN
//=============================================================================

//=============================================================================
#endif // COMISO_EIGEN3_AVAILABLE
//=============================================================================
