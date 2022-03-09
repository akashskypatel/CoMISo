//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_EIGEN3_AVAILABLE

//== INCLUDES =================================================================
#include "ConstraintTools.hh"

#include <CoMISo/Utils/MutablePriorityQueueT.hh>
#include <CoMISo/Solver/Eigen_Tools.hh>

#include <Base/Debug/DebOut.hh>

#include <limits>
#include <numeric>

namespace COMISO 
{
namespace ConstraintTools
{

//-----------------------------------------------------------------------------

void
remove_dependent_linear_constraints(
    std::vector<NConstraintInterface*>& _constraints, const double _eps)
{
  // split into linear and nonlinear
  std::vector<NConstraintInterface*> lin_const, nonlin_const;

  for (unsigned int i = 0; i < _constraints.size(); ++i)
  {
    if (_constraints[i]->is_linear() &&
        _constraints[i]->constraint_type() == NConstraintInterface::NC_EQUAL)
    {
      lin_const.push_back(_constraints[i]);
    }
    else
      nonlin_const.push_back(_constraints[i]);
  }

  remove_dependent_linear_constraints_only_linear_equality(lin_const);

  for (unsigned int i = 0; i < lin_const.size(); ++i)
    nonlin_const.push_back(lin_const[i]);

  // return filtered constraints
  _constraints.swap(nonlin_const);
}


//-----------------------------------------------------------------------------

void 
remove_dependent_linear_constraints_only_linear_equality(
    std::vector<NConstraintInterface*>& _constraints, const double _eps)
{
  DEB_enter_func;
  // make sure that constraints are available
  if (_constraints.empty())
    return;

  // 1. copy (normalized) data into gmm dynamic sparse matrix
  size_t n(_constraints[0]->n_unknowns());
  size_t m(_constraints.size());
  std::vector<double> x(n, 0.0);
  NConstraintInterface::SVectorNC g;
  HalfSparseRowMatrix A(m, n + 1);
  for (unsigned int i = 0; i < _constraints.size(); ++i)
  {
    // store rhs in last column
    A.coeffRef(i, n) = _constraints[i]->eval_constraint(x.data());
    // get and store coefficients
    _constraints[i]->eval_gradient(x.data(), g);
    double v_max(0.0);
    for (NConstraintInterface::SVectorNC::InnerIterator it(g); it; ++it)
    {
      A.coeffRef(i, it.index()) = it.value();
      v_max = std::max(v_max, std::abs(it.value()));
    }
    // normalize row
    if (v_max != 0.0)
      A.row(i) *= 1.0 / v_max;
  }

  std::vector<int> _elmn_clmn(A.rows(), -1);
  gauss_elimination(A, _elmn_clmn, IntVector(), nullptr, _eps, FL_DO_GCD);

  std::vector<size_t> keep;
  for (size_t i = 0; i < _elmn_clmn.size(); ++i)
  {
    if (_elmn_clmn[i] >= 0)
    {
      keep.push_back(i); // this rows was used to eliminate a variable, so it is
                         // is independent from the others
    }
  }

  DEB_line(2, "removed " << _constraints.size() - keep.size()
                         << " dependent linear constraints out of "
                         << _constraints.size());

  // 4. store result
  std::vector<NConstraintInterface*> new_constraints;
  for (unsigned int i = 0; i < keep.size(); ++i)
    new_constraints.push_back(_constraints[keep[i]]);

  // return linearly independent ones
  _constraints.swap(new_constraints);
}


//-----------------------------------------------------------------------------

// Dependent linear constraint removal implementation
class GaussElimination
{
public:
  GaussElimination(HalfSparseRowMatrix& _constraints, IntVector& _c_elim,
      const IntVector& _indcs_to_round, HalfSparseRowMatrix* _update_D,
      const double _eps, const uint _flags)
      : constraints_(_constraints), c_elim_(_c_elim),
        indcs_to_round_(_indcs_to_round), update_D_(_update_D), epsilon_(_eps),
        flags_(_flags)
  {}

  void run()
  {
    if (update_D_ != nullptr)
    {// setup linear transformation for rhs, start with identity
      const auto row_nmbr = constraints_.rows();
      update_D_->innerResize(row_nmbr);
      update_D_->outerResize(row_nmbr);
      for (int i = 0; i < row_nmbr; ++i)
        update_D_->coeffRef(i, i) = 1.0;
    }

    if (reorder())
      make_constraints_independent_reordering();
    else
      make_constraints_independent_no_reordering();
  }

private:
  HalfSparseRowMatrix& constraints_;
  IntVector& c_elim_;
  const IntVector& indcs_to_round_;
  double epsilon_;
  HalfSparseRowMatrix* update_D_;
  const uint flags_;

private:

  bool do_gcd() const { return (flags_ & FL_DO_GCD) == FL_DO_GCD; }
  bool reorder() const { return (flags_ & FL_REORDER) == FL_REORDER; }

  // adjust _constraints such that col _col is zero except in row _row and
  // and those rows that should be ignored.
  // col will be chosen automatically.
  // _changed_rows returns the rows that were changed.
  void make_constraint_independent(
          HalfSparseRowMatrix& _constraints,
          HalfSparseColMatrix& _constraints_c,
          int                  _row,
          int&                 _col,
    const std::vector<bool>&   _integer,
    const std::vector<bool>&   _ignore,
          std::vector<int>&    _changed_rows);

  // same as above but without returning changed rows
  void make_constraint_independent(
            HalfSparseRowMatrix& _constraints,
            HalfSparseColMatrix& _constraints_c,
      const int                  _row,
            int&                 _col,
      const std::vector<bool>&   _integer,
      const std::vector<bool>&   _ignore);

  // add _coeff * (_source_row of _source_mat)  to _target_row of _target_rmat
  // and target_cmat. set element in _zero_col to 0 if it exists
  void add_row_simultaneously(
    const Eigen::Index         _target_row,
    const double               _coeff,
    const HalfSparseRowMatrix& _source_mat,
    const Eigen::Index         _source_row,
          HalfSparseRowMatrix& _target_rmat,
          HalfSparseColMatrix& _target_cmat,
    const Eigen::Index         _zero_col = -1);

  // TODO if no gcd correction was possible, at least use a variable divisible
  // by 2 as new elim_j (to avoid in-exactness e.g. 1/3)
  bool update_constraint_gcd(
          SparseVector&     _row,
    const int               _elim_j,
          std::vector<int>& _v_gcd,
          int&              _n_ints);

  void make_constraints_independent_reordering();

  void make_constraints_independent_no_reordering();

  static int find_gcd(IntVector& _v_gcd, int& _n_ints);
};



//-----------------------------------------------------------------------------


void GaussElimination::make_constraint_independent(
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

  double elim_val = std::numeric_limits<double>::max();
  double max_elim_val = -std::numeric_limits<double>::max();

  // new: gcd
  std::vector<int> v_gcd;
  v_gcd.resize(
      COMISO_EIGEN::count_non_zeros(_constraints.row(row_id), true), -1);
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
  else if (_integer[elim_j] && elim_val > 1e-6)// TODO: why not use epsion_?
  {
    if (do_gcd() && gcd_update_valid)
    {
      // perform gcd update
      DEB_only(bool gcd_ok =) update_constraint_gcd(
          _constraints.row(row_id), elim_j, v_gcd, n_ints);
      DEB_warning_if(!gcd_ok, 1,
          " GCD update failed! " << DEB_os_str(_constraints.row(row_id)));
    }
    else
    {
      DEB_warning_if(do_gcd(), 1,
          "NO +-1 coefficient found, integer rounding cannot be guaranteed. "
          "Try using the GCD option! "
              << DEB_os_str(_constraints.row(row_id)));
      DEB_warning_if(do_gcd(), 1,
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
        add_row_simultaneously((int)c_it.index(), val, _constraints, row_id,
            _constraints, _constraints_c, elim_j);

        _changed_rows.push_back(c_it.index());

        // update linear transition of rhs
        if (update_D_ != nullptr)
          update_D_->row(c_it.index()) += val * update_D_->row(row_id);
      }
    }
  }
}

void GaussElimination::make_constraint_independent(
        HalfSparseRowMatrix& _constraints,
        HalfSparseColMatrix& _constraints_c,
        int                  _row,
        int&                 _col,
  const std::vector<bool>&   _integer,
  const std::vector<bool>&   _ignore)
{
  IntVector changed_rows;
  make_constraint_independent(_constraints, _constraints_c, _row, _col,
      _integer, _ignore, changed_rows);
}


void GaussElimination::make_constraints_independent_reordering()
{
  DEB_enter_func;

  const auto n_vars = constraints_.cols();
  const auto n_rows = constraints_.rows();

  // build round map
  std::vector<bool> roundmap(n_vars, false);
  for (size_t i = 0; i < indcs_to_round_.size(); ++i)
    roundmap[indcs_to_round_[i]] = true;

  HalfSparseColMatrix _constraints_c(constraints_);

  // init priority queue
  MutablePriorityQueueT<int, int> queue;
  queue.clear(n_rows);
  for (int i = 0; i < n_rows; ++i)
    queue.update(i, COMISO_EIGEN::count_non_zeros(constraints_.row(i), true));

  std::vector<bool> row_visited(n_rows, false);
  std::vector<size_t> row_ordering;
  row_ordering.reserve(n_rows);

  std::vector<int> changed_rows;

  while (!queue.empty())
  {
    // get next row
    int row = queue.get_next();
    row_ordering.push_back(row);
    row_visited[row] = true;

    make_constraint_independent(constraints_, _constraints_c, row, c_elim_[row],
        roundmap, row_visited, changed_rows);

    for (int row : changed_rows)
    {
      auto cur_nnz = COMISO_EIGEN::count_non_zeros(constraints_.row(row), true);
      queue.update(row, cur_nnz);
    }
  }

  constraints_.prune(0.0);

  // correct ordering
  auto c_tmp = std::move(constraints_);
  constraints_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(n_rows, n_vars);
  HalfSparseRowMatrix d_tmp;
  if (update_D_ != nullptr)
    d_tmp = *update_D_;

  std::vector<int> elim_temp(c_elim_);
  c_elim_.resize(0);
  c_elim_.resize(elim_temp.size(), -1);

  for (int i = 0; i < n_rows; ++i)
  {
    constraints_.row(i) = std::move(c_tmp.row(row_ordering[i]));
    if (update_D_ != nullptr)
      update_D_->row(i) = d_tmp.row(row_ordering[i]);

    c_elim_[i] = elim_temp[row_ordering[i]];
  }
}

void GaussElimination::make_constraints_independent_no_reordering()
{
  const size_t row_nmbr = constraints_.rows();
  const size_t var_nmbr = constraints_.cols();

  // build round map
  std::vector<bool> roundmap(var_nmbr, false);
  for (unsigned int i = 0; i < indcs_to_round_.size(); ++i)
    roundmap[indcs_to_round_[i]] = true;

  // copy constraints into column matrix (for faster update via iterators)
  HalfSparseColMatrix constraints_c(constraints_);

  std::vector<bool> visited(row_nmbr, false);
  // for all constraints
  for (int i = 0; i < row_nmbr; ++i)
  {
    visited[i] = true;
    make_constraint_independent(constraints_, constraints_c, i, c_elim_[i],
        roundmap, visited);
  }
  constraints_.prune(0.0);
}

void GaussElimination::add_row_simultaneously(
    const Eigen::Index         _target_row,
    const double               _coeff,
    const HalfSparseRowMatrix& _source_mat,
    const Eigen::Index         _source_row,
          HalfSparseRowMatrix& _target_rmat,
          HalfSparseColMatrix& _target_cmat,
    const Eigen::Index         _zero_col)
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

bool GaussElimination::update_constraint_gcd(
  SparseVector& _row,
  const int _elim_j,
  std::vector<int>& _v_gcd,
  int& _n_ints)
{
  DEB_enter_func;
  // find gcd
  double i_gcd = find_gcd(_v_gcd, _n_ints);

  if (std::abs(i_gcd) == 1.0)
    return false;

  _row *= 1.0 / i_gcd;

  LOW_CODE_QUALITY_VARIABLE_ALLOW(_elim_j);
  // TODO: really size_t? used to be gmm::size_type, but does that make sense?
  DEB_only(
      auto elim_coeff = static_cast<size_t>(std::abs(_row.coeff(_elim_j))));
  DEB_error_if(elim_coeff != 1,
      "elimination coefficient "
          << elim_coeff
          << " will (most probably) NOT lead to an integer solution!");
  return true;
}

int GaussElimination::find_gcd(std::vector<int>& _v_gcd, int& _n_ints)
{
  bool found_gcd = false;
  bool done = false;
  bool all_same = true;
  int i_gcd = -1;
  int prev_val = -1;
  // check integer coefficient pairwise
  while (!done)
  {
    // assume gcd of all pairs is the same
    all_same = true;
    for (int k = 0; k < _n_ints - 1 && !done; ++k)
    {
      // use abs(.) to get same sign needed for all_same
      _v_gcd[k] = std::abs(std::gcd(_v_gcd[k], _v_gcd[k + 1]));

      if (k > 0 && prev_val != _v_gcd[k])
        all_same = false;

      prev_val = _v_gcd[k];

      // if a 2 was found, all other entries have to be divisible by 2
      if (_v_gcd[k] == 2)
      {
        bool all_ok = true;
        for (int l = 0; l < _n_ints; ++l)
          if (abs(_v_gcd[l] % 2) != 0)
          {
            all_ok = false;
            break;
          }
        done = true;
        if (all_ok)
        {
          found_gcd = true;
          i_gcd = 2;
        }
        else
        {
          found_gcd = false;
        }
      }
    }
    // already done (by successful "2"-test)?
    if (!done)
    {
      // all gcds the same?
      // we just need to check one final gcd between first 2 elements
      if (all_same && _n_ints > 1)
      {
        _v_gcd[0] = std::abs(std::gcd(_v_gcd[0], _v_gcd[1]));
        // we are done
        _n_ints = 1;
      }

      // only one value left, either +-1 or gcd
      if (_n_ints == 1)
      {
        done = true;
        if ((_v_gcd[0]) * (_v_gcd[0]) == 1)
        {
          found_gcd = false;
          // std::cerr << __FUNCTION__ << " Info: No gcd found!" << std::endl;
        }
        else
        {
          i_gcd = _v_gcd[0];
          found_gcd = true;
          // std::cerr << __FUNCTION__ << " Info: Found gcd = " << i_gcd <<
          // std::endl;
        }
      }
    }
    // we now have n_ints-1 gcds to check next iteration
    --_n_ints;
  }
  return i_gcd;
}

void gauss_elimination(HalfSparseRowMatrix& _constraints, IntVector& _c_elim,
    const IntVector& _indcs_to_round, HalfSparseRowMatrix* _update_D,
    const double _eps, const uint _flags)
{
  GaussElimination(
      _constraints, _c_elim, _indcs_to_round, _update_D, _eps, _flags).run();
}

} // namespace ConstraintTools

//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_EIGEN3_AVAILABLE
//=============================================================================

