/*===========================================================================*\
 *                                                                           *
 *                        ExactConstraintProjection                          *
 *      Copyright (C) 2025 by Computer Graphics Group, University of Bern    *
 *                           http://cgg.unibe.ch                             *
 *                                                                           *
 *      Author: David Bommes                                                 *
 *                                                                           *
\*===========================================================================*/

#pragma once

#include <CoMISo/Config/config.hh>
#include <CoMISo/Config/CoMISoDefines.hh>
#include <CoMISo/Solver/Eigen_Tools.hh>
#include <Eigen/Sparse>

#include <numeric>

namespace COMISO {

class COMISODLLEXPORT ExactConstraintProjection
{
public:

  // Sparse vector/matrix types with integer scalars
  using DVectorInt  = Eigen::VectorXi;
  using SVectorInt  = Eigen::SparseVector<int>;

  using SMatrixIntC = COMISO_EIGEN::HalfSparseColMatrix<int>;
  using SMatrixIntR = COMISO_EIGEN::HalfSparseRowMatrix<int>;

  using PairII = std::pair<int,int>;
  using PairDD = std::pair<double,double>;

  // default constructor
  ExactConstraintProjection() {}

  // transform linear system _A*x=_b into integer-reduced row echelon form (IRREF)
  // Return value: true upon success, otherwise false
  // Note I : _A and _b are expected to have integer coefficients only
  template<class SMatrixEigen,class DVectorEigen>
  bool initialize(const SMatrixEigen& _A, const DVectorEigen& _b)
  {
    COMISO::StopWatch sw;
    sw.start();

    assert(_A.rows() == _b.size());

    // init empty row matrix and col matrix
    A_IRREF_R_ = SMatrixIntR(_A.rows(),_A.cols());
    A_IRREF_C_ = SMatrixIntC(_A.rows(),_A.cols());

    for (int k=0; k<_A.outerSize(); ++k)
      for (typename SMatrixEigen::InnerIterator it(_A,k); it; ++it)
      {
        // verify the requirement of integer coefficients
        const int val_int = static_cast<int>(it.value());
        assert(it.value() == val_int);

        // TODO: can be done faster? (probably needs assumption on row/col major input)
        A_IRREF_R_.coeffRef(it.row(),it.col()) = val_int;
        A_IRREF_C_.coeffRef(it.row(),it.col()) = val_int;
      }

    // init rhs
    b_IRREF_.resize(_b.size());
    for(int i=0; i<_b.size(); ++i)
    {
      // check requirement of integer coefficients
      assert(_b[i] == static_cast<int>(_b[i]));
      b_IRREF_[i] = _b[i];
    }

    // perform transformation
    transform_to_IRREF();

    std::cerr << "IRREF transformation took " << sw.stop()/1000.0 << " seconds" << std::endl;
  }

  // modify _x such that _A*x=b is exactly satisfied without any numerical error
  // Return value: true upon success, otherwise false
  // Note: As described in the addendum of "Exact Constraint Satisfaction for Truly Seamless Parametrization", the projection might be impossible for inhomogenous systems with _b not equal to 0
  template<class DVectorEigen>
  bool project(DVectorEigen& _x)
  {
    // 1. determine K=max_i ceil(log2(|x_i|) + 1 and delta = 2^K
    double max_abs = 0.0;
    for(int i=0; i<_x.size(); ++i)
      max_abs = std::max(max_abs, std::abs(_x[i]));

    K_ = std::ceil(std::log2(max_abs)) + static_cast<double>(K_margin_);
    delta_ = std::pow(2,K_);

    // 2. truncate free variables (collect divisors)
    double max_diff_free_variables = 0.0;
    int    max_lcm = 1;
    std::vector<size_t> dependent_variables;
    for(size_t i=0; i<is_free_variable_.size(); ++i)
      if(is_free_variable_[i])
      {
        int lcm_i = 1;
        // collect pivots
        for(SVectorInt::InnerIterator it_col(A_IRREF_C_.col(i)); it_col; ++it_col)
        {
          const int row_idx = it_col.index();
          const int col_idx = pivot_[it_col.index()];

          lcm_i = std::lcm(lcm_i,A_IRREF_C_.coeffRef(row_idx,col_idx));
        }

        max_lcm = std::max(max_lcm, lcm_i);

        double x_old = _x[i];
        _x[i] = truncate_to_F_delta(_x[i]/lcm_i)*lcm_i;
        max_diff_free_variables = std::max(max_diff_free_variables, std::abs(x_old-_x[i]));
      }
      else
        dependent_variables.push_back(i);

    // 4. compute dependent variables (use safe_dot)
    double max_diff_dependent_variables = 0.0;
    for(const auto pivot_i : dependent_variables)
    {
      assert(A_IRREF_C_.col(pivot_i).nonZeros() == 1); // a dependent variable has a single 1 in its column
      auto col_it = SVectorInt::InnerIterator(A_IRREF_C_.col(pivot_i));
      int row_idx = col_it.index();
      int C_pivot = col_it.value();

      // setup coefficients for safe_dot(a,b)
      std::vector<PairDD> dp;
      for(SVectorInt ::InnerIterator row_it(A_IRREF_R_.row(row_idx)); row_it; ++row_it)
        if(row_it.index() != pivot_i)
        {
          assert((_x[row_it.index()]/C_pivot)*C_pivot == _x[row_it.index()]); // assume divisibility
          const double a = static_cast<double>(row_it.value());
          const double b = _x[row_it.index()]/C_pivot;
          dp.emplace_back( PairDD(a,b) );
        }

      double b_div = b_IRREF_[row_idx]/C_pivot;

      // check divisibility of rhs
      if(static_cast<int>(b_div * C_pivot) != b_IRREF_[row_idx])
        std::cerr << "Warning: rhs value " << b_IRREF_[row_idx] << " is not exactly divisible by row pivot value " << C_pivot << std::endl;

      double x_old = _x[pivot_i];
      _x[pivot_i]  = b_div - safe_dot(dp);
      max_diff_dependent_variables = std::max(max_diff_dependent_variables, std::abs(x_old-_x[pivot_i]));
    }

    // output statistics on max/avg change
    std::cerr << "max_diff_free_variables      = " << max_diff_free_variables      << std::endl;
    std::cerr << "max_diff_dependent_variables = " << max_diff_dependent_variables << std::endl;
    std::cerr << "max_lcm                      = " << max_lcm << std::endl;
  }

  private:

  // transform linear system A_IRREF_R_*x=b_IRREF_ into integer-reduced row echelon form (IRREF)
  // simultaneously transfrom A_IRREF_C_*x=b_IRREF_
  bool transform_to_IRREF()
  {
    int n = A_IRREF_R_.rows();
    int m = A_IRREF_R_.cols();

    // results
    int n_pivots = 0;
    int pivot_max_abs_val = 0;

    // initialize pivots to uninitialized, i.e. -1
    pivot_.clear();
    pivot_.resize(n,-1);
    // initialize set of free variables
    is_free_variable_.clear();
    is_free_variable_.resize(m,true);


    // current number of nonzeros in row
    std::vector<int> row_nnz(n);
    // position of minimal non-zero abs(value) in row, -1 if there is none
    std::vector<int> row_min_pivot_col_idx(n);
    std::vector<int> row_min_pivot_col_val(n);

    // create two priority queues
    // one for rows with unit pivot, one for rows without unit pivot
    // priority is [number of nonzeros in row, row_idx]
    // update are lazy, i.e. modified elements are not removed but filtered on-the-fly
    std::set<PairII> to_process_with_unit_pivot;
    std::set<PairII> to_process;

    // lambda to add/update rows in the queues
    auto enqueue_row = [&](const int _i)
    {
      // reset default values
      row_nnz[_i] = 0;
      row_min_pivot_col_idx[_i] = -1;
      row_min_pivot_col_val[_i] = std::numeric_limits<int>::max();

      for(SMatrixIntR::SparseVector::InnerIterator it(A_IRREF_R_.row(_i)); it; ++it)
      {
        if(it.value() != 0)
        {
          row_nnz[_i] += 1;

          // minimal abs element in row?
          if(std::abs(it.value()) < std::abs(row_min_pivot_col_val[_i]) )
          {
            row_min_pivot_col_idx[_i] = it.index();
            row_min_pivot_col_val[_i] = it.value();
          }
        }
      }

      // only add to queue if nonzeros exist
      if(row_nnz[_i] > 0)
      {
        // has unit pivot?
        if( std::abs(row_min_pivot_col_val[_i]) == 1)
          to_process_with_unit_pivot.insert(PairII(row_nnz[_i], _i));
        else
          to_process.insert(PairII(row_nnz[_i], _i));
      }
      else
      {
        if(b_IRREF_[_i] != 0)
          std::cerr << "Warning: infeasible linear condition with zero coefficients but non-zero rhs = " << b_IRREF_[_i] << " detected during elimination" << std::endl;
      }
    };

    // enqueue all initial rows
    for(int i=0; i<n; ++i)
      enqueue_row(i);

    while(!to_process_with_unit_pivot.empty() || !to_process.empty())
    {
      // get next row
      // prefer those with unit pivot if available
      PairII cur;
      int row_cur   = -1;
      if(!to_process_with_unit_pivot.empty())
      {
        cur = *(to_process_with_unit_pivot.begin());
        to_process_with_unit_pivot.erase(to_process_with_unit_pivot.begin());
        row_cur = cur.second;

        // no unit pivot anymore (element in queue can be outdated)
        if(std::abs(row_min_pivot_col_val[row_cur]) != 1)
          continue;
      }
      else
      {
        cur = *(to_process.begin());
        to_process.erase(to_process.begin());
        row_cur = cur.second;
      }

      // outdated, or already processed, or zero row?
      if(cur.first != row_nnz[row_cur] || pivot_[row_cur] != -1 || row_nnz[row_cur] == 0)
        continue;

      // choose  pivot element to be of minimal magnitude
      const int pivot_cur = row_min_pivot_col_idx[row_cur];
      pivot_[row_cur] = pivot_cur;
      int pivot_val = row_min_pivot_col_val[row_cur];
      pivot_max_abs_val = std::max(pivot_max_abs_val, pivot_val);
      ++n_pivots;
      // mark pivot as dependent variable
      is_free_variable_[pivot_cur] = false;

      // verify data consistency
      assert(A_IRREF_R_.coeff(row_cur,pivot_cur) == pivot_val);

      // copy pivot column since it will be modified
      SVectorInt pivot_col = A_IRREF_C_.col(pivot_cur);
      for(SVectorInt::InnerIterator it_col(pivot_col); it_col; ++it_col)
        if(it_col.index() != row_cur && it_col.value() != 0) // skip current row and zero coefficients
        {
          // index of row where pivot colum will be zeroed
          const int row_elim = it_col.index();
          const int val_elim = it_col.value();

          // scale current row if pivot != 1
          if(pivot_val != 1)
          {
            for (SMatrixIntR::SparseVector::InnerIterator it_row( A_IRREF_R_.row(row_elim)); it_row; ++it_row)
            {
              // scale current row in row matrix
              it_row.valueRef() *= pivot_val;
              // also update in column matrix
              A_IRREF_C_.coeffRef(row_elim,it_row.index()) *= pivot_val;
            }
            // update rhs
            b_IRREF_[row_elim] *= pivot_val;
          }

          // subtract scaled row_cur (with pivot) from row_elim (in row and col matrix)
          for (SVectorInt::InnerIterator it_row( A_IRREF_R_.row(row_cur)); it_row; ++it_row)
          {
            int delta = val_elim*it_row.value();
            A_IRREF_R_.coeffRef(row_elim,it_row.index()) -= delta;
            A_IRREF_C_.coeffRef(row_elim,it_row.index()) -= delta;
          }
          // update rhs
          b_IRREF_[row_elim] -= val_elim*b_IRREF_[row_cur];

          A_IRREF_R_.row(row_elim).prune(0);

          // check consistency, i.e. zeroing of pivot column
          assert(A_IRREF_R_.coeff(row_elim,pivot_cur) == 0);
          assert(A_IRREF_C_.coeff(row_elim,pivot_cur) == 0);

          // determine gcd
          int gcd_row = b_IRREF_[row_elim];
          for (SVectorInt::InnerIterator it_row( A_IRREF_R_.row(row_elim)); it_row; ++it_row)
          {
            gcd_row = std::gcd(gcd_row,it_row.value());
            if(gcd_row == 1) // early termination
              break;
          }
          // divide row if gcd larger than 1
          if(gcd_row > 1)
          {
            b_IRREF_[row_elim] /= gcd_row;

            for (SVectorInt::InnerIterator it_row( A_IRREF_R_.row(row_elim)); it_row; ++it_row)
            {
              it_row.valueRef() /= gcd_row; // in row matrix
              A_IRREF_C_.coeffRef(row_elim,it_row.index()) /= gcd_row; // in column matrix
            }
          }

          // update priority queue for row_elim
          enqueue_row(row_elim);
        }
      // prune pivot column
      A_IRREF_C_.col(pivot_cur).prune(0);
    }

    std::cerr << "#independent      conditions in IRREF = " << n_pivots << std::endl;
    std::cerr << "#linear dependent conditions in IRREF = " << n-n_pivots << std::endl;
    std::cerr << "value |pivot_max| = " << pivot_max_abs_val << std::endl;
    if(!row_nnz.empty())
      std::cerr << "max #nnz in row = " << *std::max_element(row_nnz.begin(),row_nnz.end()) << std::endl;


    // verify consistency
    if(1)
      check_consistency();
  }

  void check_consistency()
  {
    A_IRREF_R_.prune(0);
    A_IRREF_C_.prune(0);

    // first verify consistency of col and row matrices
    for (int k=0; k<A_IRREF_R_.outerSize(); ++k)
      for (SVectorInt::InnerIterator it(A_IRREF_R_.row(k)); it; ++it)
      {
        int val_r = it.value();
        int val_c = A_IRREF_C_.coeff(k,it.index());
        if( val_r != val_c)
          std::cerr << "ERROR: inconsistent row and col matrix at (i,j)=(" << it.row() << "," << it.col() << ") and values "
                    << val_r << " vs. " << val_c << " detected in row matrix" << std::endl;
      }

    for (int k=0; k<A_IRREF_C_.outerSize(); ++k)
      for (SVectorInt::InnerIterator it(A_IRREF_C_.col(k)); it; ++it)
      {
        int val_r = it.value();
        int val_c = A_IRREF_R_.coeff(it.index(),k);
        if( val_r != val_c)
          std::cerr << "ERROR: inconsistent row and col matrix at (i,j)=(" << it.row() << "," << it.col() << ") and values "
                    << val_r << " vs. " << val_c << " detected in col matrix" << std::endl;
      }

    // verify that all rows without pivot element are zero
    // and pivot columns have only one non-zero
    for(int i=0; i< static_cast<int>(pivot_.size()); ++i)
      if(pivot_[i] != -1)
      {
        int nnz_c = A_IRREF_C_.col(pivot_[i]).nonZeros();
        if(nnz_c != 1)
          std::cerr << "ERROR: pivot column has #nonzeros = " << nnz_c << " but should have only one" << std::endl;
      }
      else
      {
        int nnz_r = A_IRREF_R_.row(i).nonZeros();
        if(nnz_r != 0)
          std::cerr << "ERROR: non-pivot row has #nonzeros = " << nnz_r << " but should have zero" << std::endl;
        if(b_IRREF_[i] != 0)
          std::cerr << "ERROR: zero row with non-zero rhs = " << b_IRREF_[i] << std::endl;
      }

    int nnz_max_in_row = 0;
    for (int k=0; k<A_IRREF_R_.outerSize(); ++k)
      nnz_max_in_row = std::max(nnz_max_in_row, static_cast<int>(A_IRREF_R_.row(k).nonZeros()));
    std::cerr << "max #nnz in row checked = " << nnz_max_in_row << std::endl;
  }

  // zero all mantissa bits in conflict with F_delta
  double truncate_to_F_delta(const double _d) const
  {
    if(_d >= 0.0)
      return (_d+delta_)-delta_;
    else
      return (_d-delta_)+delta_;
  }

  // assume pairs where p.first is an integer coefficient and p.second is a double coefficient
  double safe_dot(const std::vector<PairDD>& _dp) const
  {
    std::queue<PairDD> pos;
    std::queue<PairDD> neg;

    // construct queues of positive and negative terms
    // assure that p.first is always positive, i.e. for negative terms p.second must be negative
    for(const auto p : _dp)
    {
      if(p.first*p.second >= 0.0)
        pos.push(PairDD(std::abs(p.first), std::abs(p.second)));
      else
        neg.push(PairDD(std::abs(p.first), -std::abs(p.second)));
    }

    double r = 0.0;

    while(!pos.empty() || !neg.empty())
    {
      if(!pos.empty() && (r <= 0.0 || neg.empty()))
      {
        // get next pair
        const PairDD p = pos.front();
        pos.pop();

        double k = std::min(p.first, std::floor((delta_ - r) / p.second));

        // catch and handle infeasible case
        if (k == 0.0)
        {
          std::cerr << "ERROR: safe dot ended up in infeasible case ---> numerical precision loss cannot be guaranteed" << std::endl;
          std::cerr << "r=" << r << ", delta=" << delta_ << std::endl;
          // perform full update and ignore lost precision
          k = p.first;
        }

        // update r
        r += k * p.second;
        // re-add remainder
        if (k < p.first)
          pos.push(PairDD(p.first - k, p.second));
      }
      else
      {
        // get next pair
        const PairDD p = neg.front();
        neg.pop();

        double k = std::min(p.first, std::floor((-delta_ - r) / p.second));

        // catch and handle infeasible case
        if (k == 0.0)
        {
          std::cerr << "ERROR: safe dot ended up in infeasible case ---> numerical precision loss cannot be guaranteed" << std::endl;
          std::cerr << "r=" << r << ", delta=" << delta_ << std::endl;
          // perform full update and ignore lost precision
          k = p.first;
        }

        // update r
        r += k * p.second;
        // re-add remainder
        if (k < p.first)
          neg.push(PairDD(p.first - k, p.second));
     }
    }
    return r;
  }

private:
  // constraint matrix in IRREF both in row and column storage
  SMatrixIntR A_IRREF_R_;
  SMatrixIntC A_IRREF_C_;
  // rhs of constraint system in IRREF
  DVectorInt b_IRREF_;

  // colum of pivot of row i is stored in pivot_[i]
  // -1 if row has no pivot ---> zero row
  std::vector<int> pivot_;

  // distinguish between free and dependent variables. dependent are those chosen as pivots.
  std::vector<bool> is_free_variable_;

  // largest required exponent
  double K_ = 0;
  int    K_margin_ = 1;
  // largest number 2^K_
  double delta_ = 0.0;
};

}