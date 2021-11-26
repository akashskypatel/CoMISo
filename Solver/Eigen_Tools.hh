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


#ifndef COMISO_Eigen_TOOLS_HH
#define COMISO_Eigen_TOOLS_HH


//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_EIGEN3_AVAILABLE

//== INCLUDES =================================================================

#include <Base/Code/Quality.hh>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

LOW_CODE_QUALITY_SECTION_BEGIN
#include <Eigen/Dense>
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>
LOW_CODE_QUALITY_SECTION_END

#if COMISO_SUITESPARSE_AVAILABLE
#include <cholmod.h>
#endif


//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO_EIGEN
{

/** \class EigenTools Eigen_Tools.hh

    A collection of helper functions for manipulating (Eigen) matrices.
*/


//== FUNCTION DEFINITION ======================================================

/// Get matrix data (CSC matrix format) from matrix
/** Used by Cholmod wrapper
 *  @param _mat matrix
 *  @param _c uplo parameter (l, L, u, U, c, C)
 *  @param _values values vector
 *  @param _rowind row indices
 *  @param _colptr column pointer  */
template<class MatrixT, class REALT, class INTT>
void get_ccs_symmetric_data( const MatrixT&      _mat,
                             const char          _c,
                             std::vector<REALT>& _values,
                             std::vector<INTT>&  _rowind,
                             std::vector<INTT>&  _colptr );

/// Inspect the matrix (print)
/** Prints useful matrix informations such as, dimension, symmetry, zero_rows, zero_cols, nnz, max, min, max_abs, min_abs, NAN, INF
  * @param _A matrix */
template<class MatrixT>
void inspect_matrix( const MatrixT& _A);

/** checks for symmetry
  * @param _A matrix
  * @return symmetric? (bool)*/
template<class MatrixT>
bool is_symmetric( const MatrixT& _A);

template< class Eigen_MatrixT, class IntT >
void permute( const Eigen_MatrixT& _QR, const std::vector< IntT>& _Pvec, Eigen_MatrixT& _A);


/// Residuum norm of linear system
/** residuum = Ax-b
 * @param _A Matrix
 * @param _x Variables
 * @param _rhs right hand side
 * @return norm Ax-rhs */
template <class MatrixT, class VectorT>
double residuum_norm(const MatrixT& _A, const VectorT& _x, const VectorT& _rhs);

/// Convert factored LSE to quadratic representation
/** Conversion is done by computing _F^t _F where the last column is the _rhs
 * @param _F Factored Matrix (input)
 * @param _Q Quadratic Matrix (output)
 * @param _rhs right hand side (output) */
template <class ScalarT, int OPTIONS1, class Storage1T, int OPTIONS2,
    class Storage2T>
void factored_to_quadratic(
    const Eigen::SparseMatrix<ScalarT, OPTIONS1, Storage1T>& _F,
          Eigen::SparseMatrix<ScalarT, OPTIONS2, Storage2T>& _Q,
          Eigen::Matrix<ScalarT, Eigen::Dynamic, 1>& _rhs);

/// Eliminate multiple variables from a CSC matrix.
/**
 *  @param _elmn_vars indices of variables to be eliminated
 *  @param _elmn_vals values c_i of x_i to be eliminated, x_i = c_i
 *  @param _A CSC Matrix of the equation system
 *  @param _x variable vector of equation system
 *  @param _rhs right-hand side vector of equation system  */
template <class ScalarT, class IntegerT>
void eliminate_csc_vars(const std::vector<IntegerT>& _elmn_vars,
    const std::vector<ScalarT>& _elmn_vals,
    Eigen::SparseMatrix<ScalarT, Eigen::ColMajor>& _A,
    Eigen::Matrix<ScalarT, Eigen::Dynamic, 1>& _x,
    Eigen::Matrix<ScalarT, Eigen::Dynamic, 1>& _rhs);

/// Same as above but operate on data buffers
template <class ScalarT, class Integer1T, class Integer2T>
void eliminate_csc_vars(const std::vector<Integer1T>& _elmn_vars,
    const int _rows, const ScalarT* const _val_src,
    const Integer2T* const _rows_src, const Integer2T* const _cols_src,
    ScalarT* const _val_dst, Integer2T* const _rows_dst,
    Integer2T* const _cols_dst);

/// Same as above but input and output buffers are the same, i.e. work in-place.
template <class ScalarT, class Integer1T, class Integer2T>
void eliminate_csc_vars(const std::vector<Integer1T>& _elmn_vars,
    const int _n_rows, ScalarT* const _val, Integer2T* const _rows,
    Integer2T* const _cols);

/// Create a map for _n_vars variables to their new position if variables in
/// _elmn_vars are eliminated
template <class IntegerT>
std::vector<int> make_new_index_map(
    const std::vector<IntegerT>& _elmn_vars, int _n_vars);


/// do in-place elimination in CSC format by setting row and column to zero and
/// diagonal entry to one
/**
 *  @param _i index of variable to be eliminated
 *  @param _xi value the eliminated variable to set to
 *  @param _A CSC Matrix of the equation system
 *  @param _x variable vector of equation system
 *  @param _rhs right-hand side vector of equation system */
template <class ScalarT, class RealT>
void fix_var_csc_symmetric(const unsigned int _i, const ScalarT _xi,
    Eigen::SparseMatrix<RealT, Eigen::ColMajor>& _A,
    Eigen::Matrix<ScalarT, Eigen::Dynamic, 1>& _x,
    Eigen::Matrix<ScalarT, Eigen::Dynamic, 1>& _rhs);

// same as above but operate directly on csc storage buffers
// See https://en.wikipedia.org/wiki/Sparse_matrix
// for description of csc format.
template <class ScalarT, class IntegerT, class RealT>
void fix_var_csc_symmetric(const unsigned int _i, const ScalarT _xi,
    RealT* const _val, IntegerT* const _rows, IntegerT* const _cols,
    ScalarT* const _x, ScalarT* const _rhs);


#if COMISO_SUITESPARSE_AVAILABLE

/// Eigen to Cholmod_sparse interface
template<class MatrixT>
void cholmod_to_eigen( const cholmod_sparse& _AC, MatrixT& _A);

template<class MatrixT>
void eigen_to_cholmod( const MatrixT&  _A,
                     cholmod_sparse* &_AC,
                     cholmod_common* _common,
                     int             _sparsity_type = 0,
                     bool            _long_int      = false);
#endif

// convert a gmm column-sparse matrix into an Eigen sparse matrix
template <class GMM_MatrixT, class EIGEN_MatrixT>
void gmm_to_eigen(const GMM_MatrixT& _G, EIGEN_MatrixT& _E);

template <class GMM_VectorT, class EIGEN_VectorT>
void to_eigen_vec(const GMM_VectorT& _G, EIGEN_VectorT& _E);

template <class ScalarT, class EIGEN_VectorT>
void to_eigen_vec(const std::vector<ScalarT>& _G, EIGEN_VectorT& _E);

// convert an Eigen sparse matrix into a gmm sparse matrix
template <class EIGEN_MatrixT, class GMM_MatrixT>
void eigen_to_gmm(const EIGEN_MatrixT& _E, GMM_MatrixT& _G);

// convert an Eigen sparse matrix into a gmm csc matrix
template <class EIGEN_MatrixT, class GMM_CSC_MatrixT>
void eigen_to_gmm_csc(const EIGEN_MatrixT& _E, GMM_CSC_MatrixT& _G);

// convert an Eigen csc matrix into a gmm csc matrix
// Expects that the buffers of the gmm csc matrix are already big enough
template <class ScalarT, class GMM_CSC_MatrixT>
void eigen_to_gmm_csc(const Eigen::SparseMatrix<ScalarT, Eigen::ColMajor>& _E, GMM_CSC_MatrixT& _G);

template <class EIGEN_VectorT, class GMM_VectorT>
void from_eigen_vec(const EIGEN_VectorT& _E, const GMM_VectorT& _G);

template <class EIGEN_VectorT, class ScalarT>
void from_eigen_vec(const EIGEN_VectorT& _E, std::vector<ScalarT>& _v);


// Write a matrix in MatrixMarket format
template <typename MatrixT>
void write_matrix_ascii(const std::string& _filename, const MatrixT& _m);

// Load a matrix from MatrixMarket format
template <typename MatrixT>
void read_matrix_ascii(const std::string& _filename, MatrixT& _m);

// Load a sparse matrix from MatrixMarket format
template <typename ScalarT, int OPTIONS, typename StorageIndexT>
void read_matrix_ascii(const std::string& _filename,
    Eigen::SparseMatrix<ScalarT, OPTIONS, StorageIndexT>& _m);

// Write a vector in MatrixMarket format
template <typename VectorT>
void write_vector_ascii(const std::string& _filename, const VectorT& _v);

// Load a vector from MatrixMarket format
template <typename VectorT>
void read_vector_ascii(const std::string& _filename, VectorT& _v);

// Write a sparse matrix in our custom compact storage file format
template <typename ScalarT, int OPTIONS, typename StorageIndexT>
void write_matrix(const std::string& _filename,
    const Eigen::SparseMatrix<ScalarT, OPTIONS, StorageIndexT>& _m);

// Load a sparse matrix from our custom compact storage file format
template <typename ScalarT, int OPTIONS, typename StorageIndexT>
void read_matrix(const std::string& _filename,
    Eigen::SparseMatrix<ScalarT, OPTIONS, StorageIndexT>& _m);

// Write a dense matrix in our custom file format
template <typename ScalarT, int ROWS, int COLS>
void write_matrix(const std::string& _filename,
    const Eigen::Matrix<ScalarT, ROWS, COLS>& _m);

// Load a dense matrix in our custom file format
template <typename ScalarT, int ROWS, int COLS>
void read_matrix(const std::string& _filename,
    Eigen::Matrix<ScalarT, ROWS, COLS>& _m);


// Write a matrix in MatrixMarket format
template <typename MatrixT>
void write_matrix_ascii(const std::string& _filename, const MatrixT& _m);

// Load a matrix from MatrixMarket format
template <typename MatrixT>
void read_matrix_ascii(const std::string& _filename, MatrixT& _m);

// Load a sparse matrix from MatrixMarket format
template <typename ScalarT, int OPTIONS, typename StorageIndexT>
void read_matrix_ascii(const std::string& _filename,
    Eigen::SparseMatrix<ScalarT, OPTIONS, StorageIndexT>& _m);

// Write a vector in MatrixMarket format
template <typename VectorT>
void write_vector_ascii(const std::string& _filename, const VectorT& _v);

// Load a vector from MatrixMarket format
template <typename VectorT>
void read_vector_ascii(const std::string& _filename, VectorT& _v);

// Write a sparse matrix in our custom compact storage file format
template <typename ScalarT, int OPTIONS, typename StorageIndexT>
void write_matrix(const std::string& _filename,
    const Eigen::SparseMatrix<ScalarT, OPTIONS, StorageIndexT>& _m);

// Load a sparse matrix from our custom compact storage file format
template <typename ScalarT, int OPTIONS, typename StorageIndexT>
void read_matrix(const std::string& _filename,
    Eigen::SparseMatrix<ScalarT, OPTIONS, StorageIndexT>& _m);

// Write a dense matrix in our custom file format
template <typename ScalarT, int ROWS, int COLS>
void write_matrix(const std::string& _filename,
    const Eigen::Matrix<ScalarT, ROWS, COLS>& _m);

// Load a dense matrix in our custom file format
template <typename ScalarT, int ROWS, int COLS>
void read_matrix(const std::string& _filename,
    Eigen::Matrix<ScalarT, ROWS, COLS>& _m);


//=============================================================================
} // namespace COMISO_Eigen
//=============================================================================
#define COMISO_Eigen_TOOLS_TEMPLATES
#include "Eigen_ToolsT.cc"

//=============================================================================
#endif // COMISO_EIGEN3_AVAILABLE
//=============================================================================//=============================================================================
#endif // Eigen_TOOLS_HH defined
//=============================================================================

