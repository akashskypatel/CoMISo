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


#include "ConstrainedSolverT.cc"
#include <Base/Debug/DebTime.hh>

namespace COMISO {

//-----------------------------------------------------------------------------

int
ConstrainedSolver::
find_gcd(std::vector<int>& _v_gcd, int& _n_ints)
{
  bool found_gcd = false;
  bool done      = false;
  bool all_same  = true;
  int i_gcd = -1;
  int prev_val   = -1;
  // check integer coefficient pairwise
  while( !done)
  {
    // assume gcd of all pairs is the same
    all_same = true;
    for( int k=0; k<_n_ints-1 && !done; ++k)
    {
      // use abs(.) to get same sign needed for all_same
      _v_gcd[k] = abs(gcd(_v_gcd[k],_v_gcd[k+1]));

      if( k>0 && prev_val != _v_gcd[k])
        all_same = false;

      prev_val = _v_gcd[k];

      // if a 2 was found, all other entries have to be divisible by 2
      if(_v_gcd[k] == 2)
      {
        bool all_ok=true;
        for( int l=0; l<_n_ints; ++l)
          if( abs(_v_gcd[l]%2) != 0)
          {
            all_ok = false;
            break;
          }
        done      = true;
        if( all_ok )
        {
          found_gcd = true;
          i_gcd     = 2;
        }
        else
        {
          found_gcd = false;
        }
      }
    }
    // already done (by successfull "2"-test)?
    if(!done)
    {
      // all gcds the same?
      // we just need to check one final gcd between first 2 elements
      if( all_same && _n_ints >1)
      {
        _v_gcd[0]  = abs(gcd(_v_gcd[0],_v_gcd[1]));
        // we are done
        _n_ints = 1;
      }

      // only one value left, either +-1 or gcd
      if( _n_ints == 1)
      {
        done = true;
        if( (_v_gcd[0])*(_v_gcd[0]) == 1)
        {
          found_gcd = false;
          //std::cerr << __FUNCTION__ << " Info: No gcd found!" << std::endl;
        }
        else
        {
          i_gcd = _v_gcd[0];
          found_gcd = true;
          //std::cerr << __FUNCTION__ << " Info: Found gcd = " << i_gcd << std::endl;
        }
      }
    }
    // we now have n_ints-1 gcds to check next iteration
    --_n_ints;
  }
  return i_gcd;
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::solve(
  const RowMatrix&   _constraints,
  const ColMatrix&   _A,
  Vector&            _x,
  Vector&            _rhs,
  std::vector<int>&  _idx_to_round,
  double             _reg_factor,
  bool               _show_miso_settings)
{
  HalfSparseRowMatrix constraints(_constraints);
  HalfSparseColMatrix A(_A);
  solve(constraints, A, _x, _rhs, _idx_to_round, _reg_factor,
      _show_miso_settings);
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::solve(
  HalfSparseRowMatrix& _constraints,
  HalfSparseColMatrix& _A,
  Vector&              _x,
  Vector&              _rhs,
  std::vector<int>&    _idx_to_round,
  double               _reg_factor,
  bool                 _show_miso_settings)
{
  DEB_time_func_def;

#ifdef COMISO_CONSTRAINEDSOLVER_DUMP_SYSTEMS
  {
    RowMatrix constraints = _constraints;
    ColMatrix A = _A;
    write_matrix("ConstrainedSolver_constraints.mtx", constraints);
    write_matrix("ConstrainedSolver_system.mtx", A);
    write_matrix("ConstrainedSolver_rhs.vec", _rhs);
    Eigen::VectorXi idx_to_round;
    COMISO_EIGEN::to_eigen_vec(_idx_to_round, idx_to_round);
    write_matrix("ConstrainedSolver_idx_to_round.vec", idx_to_round);
  }
#endif

  // show options dialog
  if (_show_miso_settings)
    miso_.show_options_dialog();

  const size_t nrows = _A.rows();
  const size_t ncols = _A.cols();
  const size_t ncons = _constraints.rows();

  DEB_line(2, "Initital dimension: "
                  << nrows << " x " << ncols
                  << ", number of constraints: " << ncons
                  << ", number of integer variables: " << _idx_to_round.size()
                  << ", use reordering: "
                  << (use_constraint_reordering_ ? "yes" : "no"));

  // c_elim[i] = index of variable which is eliminated in condition i
  // or -1 if condition is invalid
  std::vector<int> c_elim(ncons);
  make_constraints_independent(_constraints, _idx_to_round, c_elim);

  // re-indexing vector
  std::vector<int> new_idx;

  ColMatrix Acsc;
  eliminate_constraints(
        _constraints, _A, _x, _rhs, _idx_to_round, c_elim, new_idx, Acsc);

  DEB_line(2,
      "Eliminated dimension: " << (int)Acsc.rows() << " x " << (int)Acsc.cols()
                               << "\n#nonzeros: " << (int)Acsc.nonZeros());

  miso_.solve(Acsc, _x, _rhs, _idx_to_round);

  rhs_update_table_.store(_constraints, c_elim, new_idx);
  // restore eliminated vars to fulfill the given conditions
  restore_eliminated_vars(_constraints, _x, c_elim, new_idx);

#ifdef COMISO_CONSTRAINEDSOLVER_DUMP_SYSTEMS
  write_vector("ConstrainedSolver_solution.vec", _x);
#endif
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::solve_const(
  const RowMatrix&        _constraints,
  const ColMatrix&        _A,
        Vector&           _x,
  const Vector&           _rhs,
  const std::vector<int>& _idx_to_round,
  const double            _reg_factor,
  const bool              _show_miso_settings)
{
  // copy matrices
  HalfSparseRowMatrix constraints(_constraints);
  HalfSparseColMatrix A(_A);

  // ... and vectors
  auto rhs = _rhs;
  auto idx_to_round = _idx_to_round;

  // call non-const function
  solve(constraints,
    A,
    _x,
    rhs,
    idx_to_round,
    _reg_factor,
    _show_miso_settings);
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::resolve(
        Vector& _x,
  const Vector* _constraint_rhs,
  const Vector* _rhs)
{
  DEB_time_func_def;
  // apply stored updates and eliminations to exchanged rhs
  if (_constraint_rhs)
  {
    // apply linear transformation of Gaussian elimination
    rhs_update_table_.cur_constraint_rhs_.resize(
        rhs_update_table_.D_.rows());
    rhs_update_table_.cur_constraint_rhs_ =
        rhs_update_table_.D_ * (*_constraint_rhs);

    // update rhs of stored constraints
    const size_t nc = rhs_update_table_.constraints_.cols();
    const size_t nr = rhs_update_table_.constraints_.rows();
    for (size_t i = 0; i < nr; ++i)
    {
      rhs_update_table_.constraints_.coeffRef(i, nc - 1) =
          -rhs_update_table_.cur_constraint_rhs_[i];
    }
  }
  if (_rhs)
    rhs_update_table_.cur_rhs_ = *_rhs;

  auto rhs_red = rhs_update_table_.cur_rhs_;

  rhs_update_table_.apply(rhs_update_table_.cur_constraint_rhs_, rhs_red);
  rhs_update_table_.eliminate(rhs_red);

  miso_.resolve(_x, rhs_red);

  // restore eliminated vars to fulfill the given conditions
  restore_eliminated_vars(rhs_update_table_.constraints_, _x, rhs_update_table_.c_elim_, rhs_update_table_.new_idx_);

#ifdef COMISO_CONSTRAINEDSOLVER_DUMP_SYSTEMS
  write_vector("ConstrainedSolver_solution_resolve_" +
                       std::to_string(n_resolves_++) + ".vec", _x);
#endif
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::solve(
  const RowMatrix&  _constraints,
  const RowMatrix&  _B,
  Vector&           _x,
  std::vector<int>& _idx_to_round,
  double            _reg_factor,
  bool              _show_miso_settings)
{
  DEB_enter_func;

  Vector rhs;
  ColMatrix A;
  COMISO_EIGEN::factored_to_quadratic(_B, A, rhs);
  HalfSparseRowMatrix constraints(_constraints);
  HalfSparseColMatrix A_halfsparse(A);
  // solve
  solve(constraints, A_halfsparse, _x, rhs, _idx_to_round,
      _reg_factor, _show_miso_settings);
}



//-----------------------------------------------------------------------------


void ConstrainedSolver::solve(
  HalfSparseRowMatrix& _constraints,
  HalfSparseColMatrix& _B,
  Vector&              _x,
  std::vector<int>&    _idx_to_round,
  double               _reg_factor,
  bool                 _show_miso_settings)
{
  DEB_time_func_def;

  // show options dialog
  if (_show_miso_settings)
    miso_.show_options_dialog();

  const size_t nrows = _B.rows();
  const size_t ncols = _B.cols();
  const size_t ncons = _constraints.rows();

  DEB_line(2, "Initital dimension: "
                  << nrows << " x " << ncols
                  << ", number of constraints: " << ncons
                  << ", number of integer variables: " << _idx_to_round.size()
                  << ", use reordering: "
                  << (use_constraint_reordering_ ? "yes" : "no"));

  // c_elim[i] = index of variable which is eliminated in condition i
  // or -1 if condition is invalid
  std::vector<int> c_elim(ncons);
  make_constraints_independent(_constraints, _idx_to_round, c_elim);

  // re-indexing vector
  std::vector<int> new_idx;

  eliminate_constraints(_constraints, _B, _idx_to_round, c_elim, new_idx);

  Vector rhs;
  ColMatrix Bcol;
  _B.prune(0.0);
  Bcol = _B;
  ColMatrix Acsc;
  COMISO_EIGEN::factored_to_quadratic(Bcol, Acsc, rhs);
  miso_.solve(Acsc, _x, rhs, _idx_to_round);

  rhs_update_table_.store(_constraints, c_elim, new_idx);

  // restore eliminated vars to fulfill the given conditions
  restore_eliminated_vars(_constraints, _x, c_elim, new_idx);

#ifdef COMISO_CONSTRAINEDSOLVER_DUMP_SYSTEMS
  write_vector("ConstrainedSolver_solution.vec", _x);
#endif
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::solve_const(
  const RowMatrix&        _constraints,
  const RowMatrix&        _B,
        Vector&           _x,
  const std::vector<int>& _idx_to_round,
  const double            _reg_factor,
  const bool              _show_miso_settings)
{
  // copy and vectors
  auto idx_to_round = _idx_to_round;

  // call non-const function
  solve(_constraints,
    _B,
    _x,
    idx_to_round,
    _reg_factor,
    _show_miso_settings);
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::resolve(
  const RowMatrix& _B,
        Vector&    _x,
  const Vector*    _constraint_rhs)
{
  // rhs of quadratic system is the negative of the last column of BtB,
  // or, since BtB is symmetric, the last row of BtB, without the last element.
  // Thus, rhs is equal to the last row of Bt multiplied by B.
  // The last row of Bt is the last col of B.
  // Let r_i be the ith row of B and c_i the ith element of the last colum of B.
  // rhs is sum c_i * r_i, for i = 0..n_rows
  const size_t nc = _B.cols();
  const size_t nr = _B.rows();
  Vector rhs(nc);
  rhs.setZero();
  for (size_t i = 0; i < nr; ++i)
    rhs -= _B.coeff(i, nc - 1) * _B.row(i);
  rhs.conservativeResize(nc - 1);

  // solve
  resolve(_x, _constraint_rhs, &rhs);
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::make_constraints_independent(
        HalfSparseRowMatrix& _constraints,
  const std::vector<int>&    _idx_to_round,
        std::vector<int>&    _c_elim)
{
  DEB_time_func_def;

  // setup linear transformation for rhs, start with identity
  const auto nr = _constraints.rows();
  rhs_update_table_.D_.resize(nr, nr);
  rhs_update_table_.D_.setIdentity();
  rhs_update_table_.update_D_ = rhs_update_table_.D_;

  _c_elim.clear();
  _c_elim.resize(nr, -1);

  if (use_constraint_reordering_)
  {
    make_constraints_independent_reordering(
        _constraints, _idx_to_round, _c_elim);
  }
  else
  {
    // number of variables
    const size_t n_vars = _constraints.cols();

    // build round map
    std::vector<bool> roundmap(n_vars, false);
    for (unsigned int i = 0; i < _idx_to_round.size(); ++i)
      roundmap[_idx_to_round[i]] = true;

    // copy constraints into column matrix (for faster update via iterators)
    HalfSparseColMatrix constraints_c(_constraints);

    std::vector<bool> visited(nr, false);
    // for all constraints
    for (int i = 0; i < nr; ++i)
    {
      visited[i] = true;
      make_constraint_independent(
          _constraints, constraints_c, i, _c_elim[i], roundmap, visited);
    }
    _constraints.prune(0.0);
  }

  rhs_update_table_.D_ = rhs_update_table_.update_D_;
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::make_constraints_independent_reordering(
        HalfSparseRowMatrix& _constraints,
  const std::vector<int>&    _idx_to_round,
        std::vector<int>&    _c_elim)
{
  DEB_enter_func;

  const auto n_vars = _constraints.cols();
  const auto n_rows = _constraints.rows();

  // build round map
  std::vector<bool> roundmap(n_vars, false);
  for (size_t i = 0; i < _idx_to_round.size(); ++i)
    roundmap[_idx_to_round[i]] = true;


  HalfSparseColMatrix _constraints_c(_constraints);

  // init priority queue
  MutablePriorityQueueT<int, int> queue;
  queue.clear(n_rows);
  for (int i = 0; i < n_rows; ++i)
    queue.update(i, COMISO_EIGEN::count_non_zeros(_constraints.row(i), true));

  std::vector<bool> row_visited(n_rows, false);
  std::vector<gmm::size_type> row_ordering;
  row_ordering.reserve(n_rows);

  std::vector<int> changed_rows;

  while (!queue.empty())
  {
    // get next row
    int row = queue.get_next();
    row_ordering.push_back(row);
    row_visited[row] = true;

    make_constraint_independent(_constraints, _constraints_c, row, _c_elim[row],
        roundmap, row_visited, changed_rows);

    for (int row : changed_rows)
    {
      auto cur_nnz = COMISO_EIGEN::count_non_zeros(_constraints.row(row), true);
      queue.update(row, cur_nnz);
    }
  }

  _constraints.prune(0.0);

  // correct ordering
  auto c_tmp = std::move(_constraints);
  _constraints = RowMatrix(n_rows, n_vars);
  auto d_tmp = rhs_update_table_.D_;

  std::vector<int> elim_temp(_c_elim);
  _c_elim.resize(0); _c_elim.resize(elim_temp.size(), -1);

  auto& update_D = rhs_update_table_.update_D_;
  for (int i = 0; i < n_rows; ++i)
  {
    _constraints.row(i) = std::move(c_tmp.row(row_ordering[i]));
    update_D.row(i) = d_tmp.row(row_ordering[i]);

    _c_elim[i] = elim_temp[row_ordering[i]];
  }

  rhs_update_table_.D_ = update_D;
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::make_constraint_independent(
        HalfSparseRowMatrix& _constraints,
        HalfSparseColMatrix& _constraints_c,
        int                  _row,
        int&                 _col,
  const std::vector<bool>&   _integer,
  const std::vector<bool>&   _ignore,
        std::vector<int>&    _changed_rows)
{
  DEB_enter_func;

  _changed_rows.clear();

  auto& row_id = _row;
  int n_vars = (int)_constraints.cols();

  // get elimination variable
  int elim_j = -1;
  int elim_int_j = -1;

  // iterate over current row, until variable found
  // first search for real valued variable
  // if not found for integers with value +-1
  // and finally take the smallest integer variable

  double elim_val = FLT_MAX;
  double max_elim_val = -FLT_MAX;

  // new: gcd
  std::vector<int> v_gcd;
  //TODO: why do i need .get_matrix? should be implicitly convertible.
  v_gcd.resize(COMISO_EIGEN::count_non_zeros(
                   _constraints.row(row_id), true),-1);
  int n_ints(0);
  bool gcd_update_valid(true);

  const SparseVector& row = _constraints.row(row_id);
  for (SparseVector::InnerIterator row_it(row); row_it; ++row_it)
  {
    int cur_j = static_cast<int>(row_it.index());
    if (cur_j == (int)n_vars - 1 || row_it.value() == 0)
      continue; // do not use the constant part and ignore zero values
    // found real valued var? -> finished (UPDATE: no not any more, find biggest
    // real value to avoid x/1e-13)
    if (!_integer[cur_j])
    {
      if (std::abs(row_it.value()) > max_elim_val)
      {
        elim_j = (int)cur_j;
        max_elim_val = std::abs(row_it.value());
      }
      // break;
    }
    else
    {
      double cur_row_val(std::abs(row_it.value()));
      // gcd
      // If the coefficient of an integer variable is not an integer, then
      // the variable most probably will not be. This is expected if all
      // coeffs are the same, e.g. 0.5).
      // This happens quite often in some ReForm test cases, so downgrading
      // the warning below to DEB_line at high verbosity.
      if (double(int(cur_row_val)) != cur_row_val)
      {
        DEB_line(11,
            "coefficient of integer variable is NOT integer : " << cur_row_val);
        gcd_update_valid = false;
      }

      v_gcd[n_ints] = static_cast<int>(cur_row_val);
      ++n_ints;

      // store integer closest to 1, must be greater than epsilon_
      if (std::abs(cur_row_val - 1.0) < elim_val && cur_row_val > epsilon_)
      {
        elim_int_j = (int)cur_j;
        elim_val = std::abs(cur_row_val - 1.0);
      }
    }
  }

  // first try to eliminate a valid (>epsilon_) real valued variable (safer)
  if (max_elim_val <= epsilon_)
    elim_j = elim_int_j; // use the best found integer

  _col = elim_j;

  // if no integer or real valued variable greater than epsilon_ existed, then
  // elim_j is now -1 and this row is not considered as a valid constraint

  // error check result
  if (elim_j == -1)
  {
    DEB_warning_if( // redundant or incompatible?
        std::abs(_constraints.coeff(row_id, n_vars - 1)) > epsilon_, 1,
        "incompatible condition: " << std::abs(
            _constraints.coeff(row_id, n_vars - 1)))
  }
  else if (_integer[elim_j] && elim_val > 1e-6)
  {
    if (do_gcd_ && gcd_update_valid)
    {
      // perform gcd update
      DEB_only(bool gcd_ok =) update_constraint_gcd(
                   _constraints.row(row_id), elim_j, v_gcd, n_ints);
      DEB_warning_if(!gcd_ok, 1,
          " GCD update failed! " << DEB_os_str(_constraints.row(row_id)));
    }
    else
    {
      DEB_warning_if(!do_gcd_, 1,
          "NO +-1 coefficient found, integer rounding cannot be guaranteed. "
          "Try using the GCD option! "
              << DEB_os_str(_constraints.row(row_id)));
      DEB_warning_if(do_gcd_, 1,
          "GCD of non-integer cannot be computed! "
              << DEB_os_str(_constraints.row(row_id)))
    }
  }

  // is this condition dependent?
  if (elim_j != -1)
  {
    // get elim variable value
    double elim_val_cur = _constraints.coeff(row_id, elim_j);

    // iterate over column
    const SparseVector& col = _constraints_c.col(elim_j);
    for (SparseVector::InnerIterator c_it(col); c_it; ++c_it)
    {
      if (c_it.value() == 0.0)
        continue;
      //        if( c_it.index() > i)
      if (!_ignore[c_it.index()])
      {
        //          sw.start();
        double val = -c_it.value() / elim_val_cur;
        add_row_simultaneously((int)c_it.index(), val,
            _constraints, row_id, _constraints,
            _constraints_c, elim_j);

        _changed_rows.push_back(c_it.index());

        // update linear transition of rhs
        rhs_update_table_.update_D_.row(c_it.index()) +=
            val * rhs_update_table_.update_D_.row(row_id);
      }
    }
  }

}


//-----------------------------------------------------------------------------


void ConstrainedSolver::make_constraint_independent(
        HalfSparseRowMatrix& _constraints,
        HalfSparseColMatrix& _constraints_c,
        int                  _row,
        int&                 _col,
  const std::vector<bool>&   _integer,
  const std::vector<bool>&   _ignore)
{
  std::vector<int> changed_rows;
  make_constraint_independent(_constraints, _constraints_c, _row, _col,
      _integer, _ignore, changed_rows);
}

//-----------------------------------------------------------------------------


void ConstrainedSolver::eliminate_constraints(
  const HalfSparseRowMatrix& _constraints,
        HalfSparseColMatrix& _Bcol,
        std::vector<int>&    _idx_to_round,
  const std::vector<int>&    _c_elim,
        std::vector<int>&    _new_idx)
{
  DEB_time_func_def;

  // store columns which should be eliminated
  std::vector<int> elim_cols;
  elim_cols.reserve(_c_elim.size());

  for (unsigned int i = 0; i < _c_elim.size(); ++i)
  {
    int cur_j = _c_elim[i];

    if (cur_j != -1)
    {
      double cur_val = _constraints.coeff(i, cur_j);

      // store index
      elim_cols.push_back(_c_elim[i]);

      SparseVector col = _Bcol.col(cur_j);
      for (SparseVector::InnerIterator it(col); it; ++it)
      {
        auto coeff = it.value();
        for (SparseVector::InnerIterator con_it(_constraints.row(i)); con_it; ++con_it)
          _Bcol.coeffRef(it.index(), con_it.index()) -= con_it.value() * coeff / cur_val;
      }
    }
  }

  // eliminate columns
  eliminate_columns(_Bcol, elim_cols);

  _new_idx = COMISO_EIGEN::make_new_index_map(
      elim_cols, static_cast<int>(_constraints.cols()));

  // update _idx_to_round (in place)
  std::vector<int> round_old(_idx_to_round);
  unsigned int wi = 0;
  for (unsigned int i = 0; i < _idx_to_round.size(); ++i)
  {
    if (_new_idx[_idx_to_round[i]] != -1)
    {
      _idx_to_round[wi] = _new_idx[_idx_to_round[i]];
      ++wi;
    }
  }

  // resize, sort and make unique
  _idx_to_round.resize(wi);
  sort_unique(_idx_to_round);

  DEB_line(4, "remaining         variables: " << (int)_Bcol.cols());
  DEB_line(4, "remaining integer variables: " << _idx_to_round.size());
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::eliminate_constraints(
    const HalfSparseRowMatrix& _constraints,
          HalfSparseColMatrix& _A,
          Vector&              _x,
          Vector&              _rhs,
          std::vector<int>&    _idx_to_round,
    const std::vector<int>&    _v_elim,
          std::vector<int>&    _new_idx,
          ColMatrix&           _Acsc)
{
  DEB_time_func_def;

  // store variable indices to be eliminated
  std::vector<int> elim_varids;
  elim_varids.reserve(_v_elim.size());

  rhs_update_table_.clear();

  {
    DEB_time_session_def("Constraint integration");
    SparseVector constraint_k;
    for (unsigned int i = 0; i < _v_elim.size(); ++i)
    {
      int cur_j = _v_elim[i];

      if (cur_j != -1)
      {
        const double cur_val_inv = 1.0 / _constraints.coeff(i, cur_j);

        // store index
        elim_varids.push_back(cur_j);
        rhs_update_table_.add_elim_id(cur_j);

        // copy col
        const SparseVector col = _A.col(cur_j);

        // add cur_j-th row multiplied with constraint[k] to each row k

        // loop over all constraint entries and over all column entries.
        // Ignore last element of constraint as it corresponds to rhs.
        constraint_k = _constraints.row(i);
        constraint_k.conservativeResize(_A.cols());

        for (SparseVector::InnerIterator con_it(constraint_k); con_it; ++con_it)
          rhs_update_table_.append(static_cast<int>(con_it.index()), -1.0 * (con_it.value() * cur_val_inv), cur_j, false);

        for (SparseVector::InnerIterator col_it(col); col_it; ++col_it)
        {
          for (SparseVector::InnerIterator con_it(constraint_k); con_it; ++con_it)
            _A.coeffRef(con_it.index(), col_it.index()) -=  col_it.value() * cur_val_inv * con_it.value();
        }

        // add cur_j-th col multiplied with condition[k] to each col k

        const SparseVector col_updated = _A.col(cur_j);

        for (SparseVector::InnerIterator col_it(col_updated); col_it; ++col_it)
          rhs_update_table_.append(static_cast<int>(col_it.index()), -col_it.value() * cur_val_inv, i, true);

        for (SparseVector::InnerIterator con_it(constraint_k); con_it; ++con_it)
        {
          for (SparseVector::InnerIterator col_it(col_updated); col_it; ++col_it)
            _A.coeffRef(col_it.index(), con_it.index()) -= col_it.value() * (con_it.value() * cur_val_inv);
        }
      }
    }

    // cache current rhs's
    Vector constraint_rhs_vec = -_constraints.col(_constraints.cols() - 1);
    rhs_update_table_.cur_constraint_rhs_ = constraint_rhs_vec;
    rhs_update_table_.cur_rhs_ = _rhs;
    // apply transformation due to elimination
    rhs_update_table_.apply(constraint_rhs_vec, _rhs);
  }

  {
    DEB_time_session_def("Constraint elimination");
    // eliminate vars
    _Acsc = _A;
    std::vector<double> elim_varvals(elim_varids.size(), 0);
    COMISO_EIGEN::eliminate_csc_vars(
        elim_varids, elim_varvals, _Acsc, _x, _rhs);
  }

  {
    DEB_time_session_def("Reindexing");
    _new_idx = COMISO_EIGEN::make_new_index_map(
        elim_varids, static_cast<int>(_A.cols() + 1));

    // update _idx_to_round (in place)
    unsigned int wi = 0;
    for (unsigned int i = 0; i < _idx_to_round.size(); ++i)
    {
      if (_new_idx[_idx_to_round[i]] != -1)
      {
        _idx_to_round[wi] = _new_idx[_idx_to_round[i]];
        ++wi;
      }
    }

    // resize, sort and make unique
    _idx_to_round.resize(wi);
    sort_unique(_idx_to_round);
  }
}


//-----------------------------------------------------------------------------

void ConstrainedSolver::restore_eliminated_vars(
  const HalfSparseRowMatrix& _constraints,
        Vector&               _x,
  const std::vector<int>&    _c_elim,
  const std::vector<int>&    _new_idx)
{
  DEB_time_func_def;
  // restore original ordering of _x
  Vector zeros(_new_idx.size());
  zeros.setZero(_new_idx.size());
  _x.conservativeResizeLike(zeros);
  // last variable is the constant term 1.0 for the multiplication with a
  // constraint row which contains -rhs as last column.
  _x[_x.size() - 1] = 1.0;

  // reverse iterate from prelast element
  for (int i = static_cast<int>(_new_idx.size()) - 2; i >= 0; --i)
  {
    if (_new_idx[i] != -1)
    {
      // error handling
      DEB_warning_if((i < _new_idx[i]), 1, "UNSAFE Ordering!!!")
      _x[i] = _x[_new_idx[i]];
    }
  }

  // reverse iterate
  for (int i = static_cast<int>(_c_elim.size()) - 1; i >= 0; --i)
  {
    const int cur_var = _c_elim[i];

    if (cur_var != -1)
    {
      // get variable value and set to zero
      double cur_val = _constraints.coeff(i, cur_var);

      _x[cur_var] -= _constraints.row(i).dot(_x) / cur_val;
    }
  }

  // resize
  _x.conservativeResize(_x.size() - 1);
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::add_row(
          Eigen::Index _target_row,
          double       _coeff,
    const RowMatrix&   _source_mat,
          Eigen::Index _source_row,
          ColMatrix&   _target_mat )
{
  // TODO: this might be slow
  for (RowMatrix::InnerIterator it(_source_mat, _source_row); it; ++it)
    _target_mat.coeffRef(_target_row, it.col()) += _coeff * it.value();
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::add_row_simultaneously(
    const Eigen::Index         _target_row,
    const double               _coeff,
    const HalfSparseRowMatrix& _source_mat,
    const Eigen::Index         _source_row,
          HalfSparseRowMatrix& _target_rmat,
          HalfSparseColMatrix& _target_cmat,
    const Eigen::Index         _zero_col )
{
  const SparseVector& row = _source_mat.row(_source_row);
  for (SparseVector::InnerIterator it(row); it; ++it)
  {
    if (it.value() == 0.0)
      continue;
    if (it.index() == _zero_col)
    {
      _target_rmat.coeffRef(_target_row, it.index()) = 0.0;
      _target_cmat.coeffRef(_target_row, it.index()) = 0.0;
    }
    else
    {
      _target_rmat.coeffRef(_target_row, it.index()) += _coeff * it.value();
      _target_cmat.coeffRef(_target_row, it.index()) += _coeff * it.value();
      //    if( _rmat(_row_i, r_it.index())*_rmat(_row_i, r_it.index()) <
      //    epsilon_squared_ )
      if (std::abs(_target_rmat.coeff(_target_row, it.index())) < epsilon_)
      {
        _target_rmat.coeffRef(_target_row, it.index()) = 0.0;
        _target_cmat.coeffRef(_target_row, it.index()) = 0.0;
      }
    }
  }
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::eliminate_columns(
  HalfSparseColMatrix& _M,
  const std::vector<int>& _columns)
{
  // nothing to do?
  if (_columns.size() == 0)
    return;

  const auto columns = make_sorted_unique(_columns);
  const auto n_cols_after = _M.cols() - columns.size();
  int offset = 0;
  for (int i = 0; i < n_cols_after; ++i)
  {
    while (offset < columns.size() && i + offset == columns[offset])
      ++offset;
    if (offset > 0) // prevent self assignment which sets col to zero
      _M.col(i) = _M.col(i + offset);
  }
  _M.outerResize(n_cols_after);
}


//-----------------------------------------------------------------------------


bool ConstrainedSolver::update_constraint_gcd(
          RowMatrix&        _constraints,
    const Eigen::Index      _row_i,
    const int               _elim_j,
          std::vector<int>& _v_gcd,
          int&              _n_ints)
{
  DEB_enter_func;
  // find gcd
  double i_gcd = find_gcd(_v_gcd, _n_ints);

  if (std::abs(i_gcd) == 1.0)
    return false;

  _constraints.row(_row_i) *= 1.0 / i_gcd;

  // TODO: really size_t? used to be gmm::size_type, but does that make sense?
  DEB_only(auto elim_coeff = static_cast<size_t>(std::abs(_constraints.coeff(_row_i, _elim_j)));)
  DEB_error_if(elim_coeff != 1, "elimination coefficient " << elim_coeff
    << " will (most probably) NOT lead to an integer solution!");
  return true;
}


//-----------------------------------------------------------------------------


bool ConstrainedSolver::update_constraint_gcd(
          SparseVector&     _row,
    const int               _elim_j,
          std::vector<int>& _v_gcd,
          int&              _n_ints)
{
  DEB_enter_func;
  // find gcd
  double i_gcd = find_gcd(_v_gcd, _n_ints);

  if (std::abs(i_gcd) == 1.0)
    return false;

  _row *= 1.0 / i_gcd;

  LOW_CODE_QUALITY_VARIABLE_ALLOW(_elim_j);
  // TODO: really size_t? used to be gmm::size_type, but does that make sense?
  DEB_only(auto elim_coeff = static_cast<size_t>(std::abs(_row.coeff(_elim_j)));)
  DEB_error_if(elim_coeff != 1, "elimination coefficient " << elim_coeff
    << " will (most probably) NOT lead to an integer solution!");
  return true;
}


//-----------------------------------------------------------------------------


void ConstrainedSolver::RHSUpdateTable::apply(
    Vector& _constraint_rhs, Vector& _rhs)
{
  std::vector<RHSUpdateTableEntry>::const_iterator t_it, t_end;
  t_end = table_.end();
  int cur_j = -1;
  double cur_rhs = 0.0;
  for (t_it = table_.begin(); t_it != t_end; ++t_it)
  {
    if (t_it->rhs_flag)
      _rhs[t_it->i] += t_it->f * _constraint_rhs[t_it->j];
    else
    {
      if (t_it->j != cur_j)
      {
        cur_j = t_it->j;
        cur_rhs = _rhs[cur_j];
      }
      _rhs[t_it->i] += t_it->f * cur_rhs;
    }
  }
}

void ConstrainedSolver::RHSUpdateTable::eliminate(Vector& _rhs)
{
  std::vector<int> evar(elim_var_ids_);
  std::sort(evar.begin(), evar.end());
  evar.push_back(std::numeric_limits<int>::max());

  int cur_evar_idx = 0;
  size_t nc = _rhs.size();
  for (size_t i = 0; i < nc; ++i)
  {
    size_t next_i = evar[cur_evar_idx];

    if (i != next_i)
      _rhs[i - cur_evar_idx] = _rhs[i];
    else
      ++cur_evar_idx;
  }
  _rhs.conservativeResize(nc - cur_evar_idx);
}

} // namespace COMISO

 // explicit instantiation

#include <CoMISo/Config/GmmTypes.hh>
#include <CoMISo/Config/StdTypes.hh>

namespace COMISO
{
using namespace COMISO_GMM;
using namespace COMISO_STD;

template void ConstrainedSolver::solve(const WSRowMatrix&, const WSColMatrix&,
    DoubleVector&, const DoubleVector&, IntVector&, const double, const bool);

template void ConstrainedSolver::solve(const WSRowMatrix&, const WSRowMatrix&,
    DoubleVector&, IntVector&, const double, const bool);

template void ConstrainedSolver::solve_const(const WSRowMatrix&,
    const WSColMatrix&, DoubleVector&, const DoubleVector&, const IntVector&,
    const double, const bool);

template void ConstrainedSolver::resolve(
    const WSRowMatrix&, DoubleVector&, const DoubleVector*);

template void ConstrainedSolver::solve_const(const WSRowMatrix&,
    const WSRowMatrix&, DoubleVector&, const IntVector&, const double,
    const bool);

template void ConstrainedSolver::solve_const(const RSRowMatrix&,
    const RSRowMatrix&, DoubleVector&, const IntVector&, const double,
    const bool);

}//namespace COMISO
