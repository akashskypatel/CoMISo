/*===========================================================================*\
 *                                                                           *
 *                        ExactConstraintProjection                          *
 *      Copyright (C) 2025 by Computer Graphics Group, University of Bern    *
 *                           http://cgg.unibe.ch                             *
 *                                                                           *
 *      Author: David Bommes                                                 *
 *                                                                           *
\*===========================================================================*/

#include "ExactConstraintProjection.hh"
#include <Base/Debug/DebOut.hh>

#include <numeric>
#include <set>
#include <cassert>

namespace COMISO {

bool
ExactConstraintProjection::
initialize(const SMatrixIntR& _A, const DVectorInt& _b)
{
  assert(_A.rows() == _b.size());
  // copy matrix and rhs
  A_IRREF_R_ = _A;
  b_IRREF_ = _b;

  // initialize column version
  A_IRREF_C_ = SMatrixIntC(_A.rows(), _A.cols());
  for(int i=0; i<A_IRREF_R_.rows(); ++i)
    for(SVectorInt::InnerIterator row_it(A_IRREF_R_.row(i)); row_it; ++row_it)
      A_IRREF_C_.coeffRef(i,row_it.index()) = row_it.value();

  // perform transformation
  bool valid = transform_to_IRREF();
  return valid;
}


//-----------------------------------------------------------------------------


bool
ExactConstraintProjection::
transform_to_IRREF()
{
  DEB_enter_func
  COMISO::StopWatch sw; sw.start();

  int n = A_IRREF_R_.rows();
  int m = A_IRREF_R_.cols();

  // results
  int  n_pivots = 0;
  int  pivot_max_abs_val = 0;
  bool valid = true;

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
      if(std::abs(it.value()) != 0)
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
      {
        DEB_line(2, "Warning: infeasible linear condition with zero coefficients but non-zero rhs = "
                << b_IRREF_[_i] << " detected during elimination");
        valid = false;
      }
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

    // choose  pivot element to be of minimal magnitude and with minimal number of nonzeros in column
    int pivot_cur = row_min_pivot_col_idx[row_cur];
    int pivot_val = row_min_pivot_col_val[row_cur];

    // determine best pivot element ---> minimal number of nonzeros in column
    int nnz_pivot_col = A_IRREF_C_.col(pivot_cur).nonZeros();
    for (SMatrixIntR::SparseVector::InnerIterator it_row( A_IRREF_R_.row(row_cur)); it_row; ++it_row)
      if(std::abs(it_row.value()) != 0)
      {
        int nnz = A_IRREF_C_.col(it_row.index()).nonZeros();
        if (nnz < nnz_pivot_col &&
            (!prioritize_unit_pivots_ || std::abs(it_row.value()) <= std::abs(pivot_val)))
        {
          nnz_pivot_col = nnz;
          pivot_cur = it_row.index();
          pivot_val = it_row.value();

        }
      }

    // set pivot
    pivot_[row_cur] = pivot_cur;
    pivot_max_abs_val = std::max(pivot_max_abs_val, pivot_val);
    ++n_pivots;
    // mark pivot as dependent variable
    is_free_variable_[pivot_cur] = false;

    if(enable_detailed_logging_)
    {
      int nnz_min   = INT_MAX;
      int nnz_min_1 = INT_MAX;
      for (SMatrixIntR::SparseVector::InnerIterator it_row( A_IRREF_R_.row(row_cur)); it_row; ++it_row)
      {
        int nnz = A_IRREF_C_.col(it_row.index()).nonZeros();
        nnz_min = std::min(nnz,nnz_min);
        if(std::abs(it_row.value()) == 1)
          nnz_min_1 = std::min(nnz, nnz_min_1);
      }

      DEB_line(2, "*** process row " << row_cur << ", remaining "
                << to_process_with_unit_pivot.size() + to_process.size()
                << ", pivot_idx = " << pivot_cur
                << ", pivot_val = " << pivot_val
                << ", #nnz in pivot col = " << A_IRREF_C_.col(pivot_cur).nonZeros()
                << ", min #nnz in col = " << nnz_min
                << ", min #nnz in col with |1| coeff = " << nnz_min_1);
    }

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

        A_IRREF_R_.prune_row(row_elim, 0);

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
    A_IRREF_C_.prune_col(pivot_cur, 0);
  }

  DEB_line(2, "#independent      conditions in IRREF = " << n_pivots);
  DEB_line(2, "#linear dependent conditions in IRREF = " << n-n_pivots);
  DEB_line(2, "value |pivot_max| = " << pivot_max_abs_val);
  if(!row_nnz.empty())
  {
    DEB_line(2, "max #nnz in row = " << *std::max_element(row_nnz.begin(), row_nnz.end()));
    DEB_line(2, "#nnz in IRREF   = " << std::accumulate(row_nnz.begin(), row_nnz.end(), 0));
  }

  DEB_line(2, "transform_to_IRREF took " << sw.stop()/1000.0 << "seconds");

  // verify consistency
  if(1)
    check_consistency();

  return valid;
}

//-----------------------------------------------------------------------------


double
ExactConstraintProjection::
round_to_F_delta(const double _d) const
{
  // the floating point operation automatically rounds to the closest number (it's not only a truncation!!!)
  if(_d >= 0.0)
    return (_d+delta_)-delta_;
  else
    return (_d-delta_)+delta_;
}


//-----------------------------------------------------------------------------


bool
ExactConstraintProjection::
is_in_F_delta(const double _d) const
{
  if(_d >= 0.0)
    return (((_d+delta_)-delta_) == _d);
  else
    return (((_d-delta_)+delta_) == _d);
}



//-----------------------------------------------------------------------------


double
ExactConstraintProjection::
safe_dot(const std::vector<PairDD>& _dp) const
{
  std::queue<PairDD> pos;
  std::queue<PairDD> neg;

  // construct queues of positive and negative terms
  // assure that p.first is always positive, i.e. for negative terms p.second must be negative
  for(const auto& p : _dp)
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
        DEB_error("safe_dot: ended up in infeasible case for pos ---> numerical precision loss might occur: "
                << "r=" << r << ", delta=" << delta_ << ", pos.size()= " << pos.size() << ", neg.size()= " << neg.size());
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
        DEB_error("safe_dot: ended up in infeasible case for neg ---> numerical precision loss might occur: "
            << "r=" << r << ", delta=" << delta_ << ", pos.size()= " << pos.size() << ", neg.size()= " << neg.size());
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


//-----------------------------------------------------------------------------


void
ExactConstraintProjection::
check_consistency()
{
  DEB_enter_func
  A_IRREF_R_.prune(0.);
  A_IRREF_C_.prune(0.);

  // first verify consistency of col and row matrices
  for (int k=0; k<A_IRREF_R_.outerSize(); ++k)
    for (SVectorInt::InnerIterator it(A_IRREF_R_.row(k)); it; ++it)
    {
      int val_r = it.value();
      int val_c = A_IRREF_C_.coeff(k,it.index());
      if( val_r != val_c)
        DEB_error("ExactConstraintProjection: inconsistent row and col matrix at (i,j)=(" << it.row() << "," << it.col()
                << ") and values " << val_r << " vs. " << val_c << " detected in row matrix");
    }

  for (int k=0; k<A_IRREF_C_.outerSize(); ++k)
    for (SVectorInt::InnerIterator it(A_IRREF_C_.col(k)); it; ++it)
    {
      int val_r = it.value();
      int val_c = A_IRREF_R_.coeff(it.index(),k);
      if( val_r != val_c)
        DEB_error("ExactConstraintProjection: inconsistent row and col matrix at (i,j)=(" << it.row() << "," << it.col()
                << ") and values " << val_r << " vs. " << val_c << " detected in col matrix");
    }

  // verify that all rows without pivot element are zero
  // and pivot columns have only one non-zero
  for(int i=0; i< static_cast<int>(pivot_.size()); ++i)
    if(pivot_[i] != -1)
    {
      int nnz_c = A_IRREF_C_.col(pivot_[i]).nonZeros();
      if(nnz_c != 1)
        DEB_error("ExactConstraintProjection: pivot column has #nonzeros = " << nnz_c << " but should have only one");
    }
    else
    {
      int nnz_r = A_IRREF_R_.row(i).nonZeros();
      if(nnz_r != 0)
        DEB_error("ExactConstraintProjection: non-pivot row has #nonzeros = " << nnz_r << " but should have zero");
      if(b_IRREF_[i] != 0)
        DEB_error("ExactConstraintProjection: zero row with non-zero rhs = " << b_IRREF_[i]);
    }

  int nnz_max_in_row = 0;
  for (int k=0; k<A_IRREF_R_.outerSize(); ++k)
    nnz_max_in_row = std::max(nnz_max_in_row, static_cast<int>(A_IRREF_R_.row(k).nonZeros()));
  DEB_line(2, "max #nnz in row checked = " << nnz_max_in_row);
}

} // NAMESPACE COMISO
