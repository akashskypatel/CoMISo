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

using BoolVector = std::vector<bool>;

//-----------------------------------------------------------------------------

ConstraintRemovalResult remove_dependent_linear_constraints(
        ConstraintVector& _constraints, const double _eps, const EliminationMethod _elim_method)
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

  ConstraintRemovalResult crr = remove_dependent_linear_constraints_only_linear_equality(lin_const, _eps, _elim_method);

  for (unsigned int i = 0; i < lin_const.size(); ++i)
    nonlin_const.push_back(lin_const[i]);

  // return filtered constraints
  _constraints.swap(nonlin_const);

  return crr;
}

//-----------------------------------------------------------------------------

ConstraintRemovalResult remove_dependent_linear_constraints_only_linear_equality(
        ConstraintVector& _constraints, const double _eps, const EliminationMethod _elim_method)
{
  DEB_enter_func;

  switch(_elim_method)
  {
    case ELIMINATION_EIGEN:
#if COMISO_EIGEN3_AVAILABLE
      return remove_dependent_linear_constraints_only_linear_equality_eigen(_constraints, _eps);
#endif
    case ELIMINATION_GMM:
#if COMISO_GMM_AVAILABLE
      return remove_dependent_linear_constraints_only_linear_equality_gmm(_constraints, _eps);
#endif
    default:
    DEB_warning( 1, "remove_dependent_linear_constraints_only_linear_equality was called with invalid EliminationMethod, or neither Eigen3 or GMM are available");
      return ConstraintRemovalResult();
  }
}

//-----------------------------------------------------------------------------


// TODO: replace with std::gcd when C++17 is available
int gcd(int _a, int _b)
{
  while (_b != 0)
  {
    int t(_b);
    _b = _a % _b;
    _a = t;
  }
  return _a;
}

//-----------------------------------------------------------------------------


ConstraintRemovalResult remove_dependent_linear_constraints_only_linear_equality_eigen(
    ConstraintVector& _constraints, const double _eps)
{
  DEB_enter_func;

  ConstraintRemovalResult result;

  // make sure that constraints are available
  if (_constraints.empty())
    return result;

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

  IntVector elmn_clmn_indcs;
  gauss_elimination(A, elmn_clmn_indcs, IntVector(), nullptr, _eps, FL_DO_GCD);

  std::vector<size_t> keep;
  for (size_t i = 0; i < elmn_clmn_indcs.size(); ++i)
  {
    if (elmn_clmn_indcs[i] >= 0)
    {
      keep.push_back(i); // this rows was used to eliminate a variable, so it is
                         // is independent from the others
    }
  }

  DEB_line(2, "removed " << _constraints.size() - keep.size()
                         << " dependent linear constraints out of "
                         << _constraints.size());

  //  number of removed constraints
  result.n_constraints_eliminated = _constraints.size()-keep.size();

  // only update _constraints if at least one constraint has been removed, otherwise preserve order and leave _constraints untouched!!!
  if(result.n_constraints_eliminated > 0)
  {
    DEB_line(2, "removed " << result.n_constraints_eliminated <<
                           " dependent linear constraints out of " << _constraints.size());
    // 4. store result
    std::vector<NConstraintInterface *> new_constraints;
    for (unsigned int i = 0; i < keep.size(); ++i)
      new_constraints.push_back(_constraints[keep[i]]);

    // return linearly independent ones
    _constraints.swap(new_constraints);
  }

  return result;
}


//-----------------------------------------------------------------------------


class GaussElimination
{
public:
  GaussElimination(HalfSparseRowMatrix& _constraints,
      IntVector& _elmn_clmn_indcs, const IntVector& _indcs_to_round,
      HalfSparseRowMatrix* _update_D, const double _eps, const uint _flags)
      : constraints_(_constraints), constraints_clmn_(constraints_),
        elmn_clmn_indcs_(_elmn_clmn_indcs), indcs_to_round_(_indcs_to_round),
        round_map_(constraints_.cols(), false), epsilon_(_eps),
        update_D_(_update_D), flags_(_flags),
        visited_(constraints_.rows(), false),
        n_rows_linearly_dependent_(0),
        n_rows_contradicting_(0)

  {
    for (const auto indx_to_round : indcs_to_round_)
      round_map_[indx_to_round] = true; // build round map
  }

  GaussEliminationResult run()
  {
    const auto row_nmbr = constraints_.rows();
    elmn_clmn_indcs_.clear();
    elmn_clmn_indcs_.resize(row_nmbr, -1);

    // reset
    n_rows_linearly_dependent_=0;
    n_rows_contradicting_=0;

    if (update_D_ != nullptr)
    {// setup linear transformation for rhs, start with identity
      update_D_->innerResize(row_nmbr);
      update_D_->outerResize(row_nmbr);
      for (int i = 0; i < row_nmbr; ++i)
        update_D_->coeffRef(i, i) = 1.0;
    }

    if (reorder())
      make_independent_reordering();
    else
      make_independent_no_reordering();

    GaussEliminationResult result;
    result.n_rows_linearly_dependent = n_rows_linearly_dependent_;
    result.n_rows_contradicting = n_rows_contradicting_;

    return result;
  }

private:
  HalfSparseRowMatrix& constraints_;
  // constraints copy into column matrix (for faster update via iterators)
  HalfSparseColMatrix constraints_clmn_;
  IntVector& elmn_clmn_indcs_;
  const IntVector& indcs_to_round_;
  BoolVector round_map_;
  double epsilon_;
  HalfSparseRowMatrix* update_D_;
  const uint flags_;
  BoolVector visited_;
  IntVector chng_row_indcs_; // storage for the changed rows in make_independent
  int       n_rows_linearly_dependent_;
  int       n_rows_contradicting_;

private:

  bool do_gcd() const { return (flags_ & FL_DO_GCD) == FL_DO_GCD; }
  bool reorder() const { return (flags_ & FL_REORDER) == FL_REORDER; }

  // Chooses a non-zero column in the specified row of the constraints matrix.
  // The choice is stored elmn_clmn_indcs_[_row_indx]. By adding multiples of
  // the selected row, the constraint matrix is adjusted so that the chosen
  // column is zero in all rows except the input row and those that have already
  // been visited. chng_row_indcs_ contains the rows that have been changed.
  void make_independent(const int _row_indx);

  // Add _coeff * the _source_row of the constraints to _target_row of the
  // constraints. Set the element in _zero_col to 0 if it is non-zero. The
  // operation is performed simultaneously to both copies of the constraints
  // matrix (row and column ordered).
  void add_row_simultaneously(const int _target_row, const double _coeff,
      const int _source_row, const int _zero_col = -1);

  void make_independent_reordering();
  void make_independent_no_reordering();

  // TODO if no gcd correction was possible, at least use a variable divisible
  // by 2 as new elim_j (to avoid in-exactness e.g. 1/3)
  static bool update_constraint_gcd(
      SparseVector& _row, const int _elim_j, IntVector& _v_gcd, int& _n_ints);

  static int find_gcd(IntVector& _v_gcd, int& _n_ints);
};


//-----------------------------------------------------------------------------


void GaussElimination::make_independent(int _row_indx)
{
  DEB_enter_func;

  visited_[_row_indx] = true;
  chng_row_indcs_.clear();

  const int n_vars = (int)constraints_.cols();

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
  IntVector v_gcd;
  v_gcd.resize(
      COMISO_EIGEN::count_non_zeros(constraints_.row(_row_indx), true), -1);
  int n_ints(0);
  bool gcd_update_valid(true);

  const SparseVector& row = constraints_.row(_row_indx);
  for (SparseVector::InnerIterator row_it(row); row_it; ++row_it)
  {
    int cur_j = static_cast<int>(row_it.index());
    if (cur_j == (int)n_vars - 1 || row_it.value() == 0)
      continue; // do not use the constant part and ignore zero values
    // found real valued var? -> finished (UPDATE: no not any more, find biggest
    // real value to avoid x/1e-13)
    if (!round_map_[cur_j])
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

  elmn_clmn_indcs_[_row_indx] = elim_j;

  // if no integer or real valued variable greater than epsilon_ existed, then
  // elim_j is now -1 and this row is not considered as a valid constraint

  // error check result
  if (elim_j == -1)
  {
    DEB_warning_if( // redundant or incompatible?
        std::abs(constraints_.coeff(_row_indx, n_vars - 1)) > epsilon_, 1,
        "incompatible condition: " << std::abs(
            constraints_.coeff(_row_indx, n_vars - 1)))
  }
  else if (round_map_[elim_j] && elim_val > 1e-6)// TODO: why not use epsilon_?
  {
    if (do_gcd() && gcd_update_valid)
    {
      // perform gcd update
      DEB_only(bool gcd_ok =) update_constraint_gcd(
          constraints_.row(_row_indx), elim_j, v_gcd, n_ints);
      DEB_warning_if(!gcd_ok, 1,
          " GCD update failed! " << DEB_os_str(constraints_.row(_row_indx)));
    }
    else
    {
      DEB_warning_if(do_gcd(), 1,
          "NO +-1 coefficient found, integer rounding cannot be guaranteed. "
          "Try using the GCD option! "
              << DEB_os_str(constraints_.row(_row_indx)));
      DEB_warning_if(do_gcd(), 1,
          "GCD of non-integer cannot be computed! "
              << DEB_os_str(constraints_.row(_row_indx)))
    }
  }

  if (elim_j == -1) // is this condition dependent?
  {
    // record number of linearly dependent rows
    ++n_rows_linearly_dependent_;
    // record number of contradicting rows
    if(std::abs(constraints_.coeff(_row_indx, n_vars - 1)) > epsilon_)
      ++n_rows_contradicting_;

    return;
  }

  // get elim variable value
  double elim_val_cur = constraints_.coeff(_row_indx, elim_j);

  // iterate over column
  const SparseVector& col = constraints_clmn_.col(elim_j);
  for (SparseVector::InnerIterator c_it(col); c_it; ++c_it)
  {
    if (c_it.value() == 0.0)
      continue;
    //        if( c_it.index() > i)
    if (!visited_[c_it.index()])
    {
      double val = -c_it.value() / elim_val_cur;
      add_row_simultaneously((int)c_it.index(), val, _row_indx, elim_j);

      chng_row_indcs_.push_back(c_it.index());

      // update linear transition of rhs
      if (update_D_ != nullptr)
        update_D_->row(c_it.index()) += val * update_D_->row(_row_indx);
    }
  }
}


void GaussElimination::make_independent_reordering()
{
  DEB_enter_func;

  const auto n_vars = constraints_.cols();
  const auto n_rows = constraints_.rows();

  // init priority queue
  MutablePriorityQueueT<int, int> queue;
  queue.clear(n_rows);

  const auto queue_update_row = [this, &queue](const int _i)
  {
    queue.update(_i, COMISO_EIGEN::count_non_zeros(constraints_.row(_i), true));
  };

  for (int i = 0; i < n_rows; ++i)
    queue_update_row(i);

  IntVector row_ordering;
  row_ordering.reserve(n_rows);

  while (!queue.empty())
  {
    row_ordering.push_back(queue.get_next());
    make_independent(row_ordering.back());

    for (int i : chng_row_indcs_)
      queue_update_row(i);
  }

  constraints_.prune(0.0);

  // correct ordering
  auto c_tmp = std::move(constraints_);
  constraints_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(n_rows, n_vars);
  HalfSparseRowMatrix d_tmp;
  if (update_D_ != nullptr)
    d_tmp = *update_D_;

  IntVector elim_temp(elmn_clmn_indcs_);
  elmn_clmn_indcs_.resize(0);
  elmn_clmn_indcs_.resize(elim_temp.size(), -1);

  for (int i = 0; i < n_rows; ++i)
  {
    constraints_.row(i) = std::move(c_tmp.row(row_ordering[i]));
    if (update_D_ != nullptr)
      update_D_->row(i) = d_tmp.row(row_ordering[i]);

    elmn_clmn_indcs_[i] = elim_temp[row_ordering[i]];
  }
}

void GaussElimination::make_independent_no_reordering()
{
  // for all constraints
  for (int i = 0, n = constraints_.rows(); i < n; ++i)
    make_independent(i);

  constraints_.prune(0.0);
}

void GaussElimination::add_row_simultaneously(const int _target_row,
    const double _coeff, const int _source_row,
    const int _zero_col)
{
  const SparseVector& row = constraints_.row(_source_row);
  for (SparseVector::InnerIterator it(row); it; ++it)
  {
    if (it.value() == 0.0)
      continue;
    if (it.index() == _zero_col)
    {
      constraints_.coeffRef(_target_row, it.index()) = 0.0;
      constraints_clmn_.coeffRef(_target_row, it.index()) = 0.0;
    }
    else
    {
      constraints_.coeffRef(_target_row, it.index()) += _coeff * it.value();
      constraints_clmn_.coeffRef(_target_row, it.index()) += _coeff * it.value();
      //    if( _rmat(_row_i, r_it.index())*_rmat(_row_i, r_it.index()) <
      //    epsilon_squared_ )
      if (std::abs(constraints_.coeff(_target_row, it.index())) < epsilon_)
      {
        constraints_.coeffRef(_target_row, it.index()) = 0.0;
        constraints_clmn_.coeffRef(_target_row, it.index()) = 0.0;
      }
    }
  }
}

bool GaussElimination::update_constraint_gcd(
    SparseVector& _row, const int _elim_j, IntVector& _v_gcd, int& _n_ints)
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

int GaussElimination::find_gcd(IntVector& _v_gcd, int& _n_ints)
{
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
      _v_gcd[k] = std::abs(gcd(_v_gcd[k], _v_gcd[k + 1]));

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
          i_gcd = 2;
      }
    }
    // already done (by successful "2"-test)?
    if (!done)
    {
      // all gcds the same?
      // we just need to check one final gcd between first 2 elements
      if (all_same && _n_ints > 1)
      {
        _v_gcd[0] = std::abs(gcd(_v_gcd[0], _v_gcd[1]));
        // we are done
        _n_ints = 1;
      }

      // only one value left, either +-1 or gcd
      if (_n_ints == 1)
      {
        done = true;
        if ((_v_gcd[0]) * (_v_gcd[0]) != 1)
          i_gcd = _v_gcd[0];
      }
    }
    // we now have n_ints-1 gcds to check next iteration
    --_n_ints;
  }
  return i_gcd;
}

GaussEliminationResult gauss_elimination(HalfSparseRowMatrix& _constraints,
    IntVector& _elmn_clmn_indcs, const IntVector& _indcs_to_round,
    HalfSparseRowMatrix* _update_D, const double _eps, const uint _flags)
{
  return GaussElimination(
      _constraints, _elmn_clmn_indcs, _indcs_to_round, _update_D, _eps, _flags)
      .run();
}


//-----------------------------------------------------------------------------

#if COMISO_GMM_AVAILABLE

//-----------------------------------------------------------------------------

static gmm::size_type
find_max_abs_coeff(SVectorGMM& _v)
{
  size_t n = _v.size();
  gmm::size_type imax(0);
  double       vmax(-1.0);

  gmm::linalg_traits<SVectorGMM>::const_iterator c_it   = gmm::vect_const_begin(_v);
  gmm::linalg_traits<SVectorGMM>::const_iterator c_end  = gmm::vect_const_end(_v);

  for(; c_it != c_end; ++c_it)
    if(c_it.index() != n-1)
      if(std::abs(*c_it) > vmax)
      {
        imax = c_it.index();
        vmax = std::abs(*c_it);
      }

  return imax;
}


static void
add_row_simultaneously( gmm::size_type _row_i,
                        double      _coeff,
                        SVectorGMM& _row,
                        RMatrixGMM& _rmat,
                        CMatrixGMM& _cmat,
                        const double _eps )
{
  typedef gmm::linalg_traits<SVectorGMM>::const_iterator RIter;
  RIter r_it  = gmm::vect_const_begin(_row);
  RIter r_end = gmm::vect_const_end(_row);

  for(; r_it!=r_end; ++r_it)
  {
    _rmat(_row_i, r_it.index()) += _coeff*(*r_it);
    _cmat(_row_i, r_it.index()) += _coeff*(*r_it);
    if( std::abs(_rmat(_row_i, r_it.index())) < _eps )
    {
      _rmat(_row_i, r_it.index()) = 0.0;
      _cmat(_row_i, r_it.index()) = 0.0;
    }
  }
}

ConstraintTools::ConstraintRemovalResult
remove_dependent_linear_constraints_only_linear_equality_gmm( std::vector<NConstraintInterface*>& _constraints, const double _eps)
{
  DEB_enter_func;

  ConstraintRemovalResult result;

  // make sure that constraints are available
  if(_constraints.empty()) return result;

  // 1. copy (normalized) data into gmm dynamic sparse matrix
  size_t n(_constraints[0]->n_unknowns());
  size_t m(_constraints.size());
  std::vector<double> x(n, 0.0);
  NConstraintInterface::SVectorNC g;
  RMatrixGMM A;
  gmm::resize(A,m, n+1);
  for(unsigned int i=0; i<_constraints.size(); ++i)
  {
    // store rhs in last column
    A(i,n) = _constraints[i]->eval_constraint(x.data());
    // get and store coefficients
    _constraints[i]->eval_gradient(x.data(), g);
    double v_max(0.0);
    for (NConstraintInterface::SVectorNC::InnerIterator it(g); it; ++it)
    {
      A(i,it.index()) = it.value();
      v_max = std::max(v_max, std::abs(it.value()));
    }
    // normalize row
    if(v_max != 0.0)
      gmm::scale(A.row(i), 1.0/v_max);
  }

  // 2. get additionally column matrix to exploit column iterators
  CMatrixGMM Ac;
  gmm::resize(Ac, gmm::mat_nrows(A), gmm::mat_ncols(A));
  gmm::copy(A, Ac);

  // 3. initialize priorityqueue for sorting
  // init priority queue
  MutablePriorityQueueT<gmm::size_type, gmm::size_type> queue;
  queue.clear(m);
  for (gmm::size_type i = 0; i<m; ++i)
  {
    gmm::size_type cur_nnz = gmm::nnz( gmm::mat_row(A,i));
    if (A(i,n) != 0.0)
      --cur_nnz;

    queue.update(i, cur_nnz);
  }

  // track row status -1=undecided, 0=remove, 1=keep
  std::vector<int> row_status(m, -1);
  std::vector<gmm::size_type> keep;
//  std::vector<int> remove;

//  std::vector<double> val_keep, val_remove;

  // for all conditions
  while(!queue.empty())
  {
    // get next row
    gmm::size_type i = queue.get_next();
    gmm::size_type j = find_max_abs_coeff(A.row(i));
    double aij = A(i,j);
    if(std::abs(aij) <= _eps)
    {
//      std::cerr << "drop " << aij << "in row " << i << "and column " << j << std::endl;
      // constraint is linearly dependent
      row_status[i] = 0;
      if(std::abs(A(i,n)) > _eps) {
        ++result.n_infeasible_detected;
        DEB_warning( 1, "Warning: found dependent constraint with nonzero rhs " << A(i,n));
      }

//      val_remove.push_back(std::abs(aij));
    }
    else
    {
//      std::cerr << "keep " << aij << "in row " << i << "and column " << j << std::endl;
//      val_keep.push_back(std::abs(aij));

      // constraint is linearly independent
      row_status[i] = 1;
      keep.push_back(i);

      // update undecided constraints
      // copy col
      SVectorGMM col = Ac.col(j);

      // copy row
      SVectorGMM row = A.row(i);

      // iterate over column
      gmm::linalg_traits<SVectorGMM>::const_iterator c_it   = gmm::vect_const_begin(col);
      gmm::linalg_traits<SVectorGMM>::const_iterator c_end  = gmm::vect_const_end(col);

      for(; c_it != c_end; ++c_it)
        if( row_status[c_it.index()] == -1) // only process unvisited rows
        {
          // row idx
          gmm::size_type k = c_it.index();

          double s = -(*c_it)/aij;

//          add_row_simultaneously( k, s, row, A, Ac, _eps); // DO NOT TRUNCATE DURING COMPUTATION TO PRESERVE DOUBLE PRECISION UNTIL DECISION!!!
          add_row_simultaneously( k, s, row, A, Ac, 0.0);
          // make sure the eliminated entry is numerically 0 on all other rows
          A( k, j) = 0;
          Ac(k, j) = 0;

          gmm::size_type cur_nnz = gmm::nnz( gmm::mat_row(A,k));
          if( A(k,n) != 0.0)
            --cur_nnz;

          queue.update(k, cur_nnz);
        }
    }
  }

  //  number of removed constraints
  result.n_constraints_eliminated = _constraints.size()-keep.size();

  // only update _constraints if at least one constraint has been removed, otherwise preserve order and leave _constraints untouched!!!
  if(result.n_constraints_eliminated > 0)
  {
    DEB_line(2, "removed " << result.n_constraints_eliminated <<
                           " dependent linear constraints out of " << _constraints.size());

    // 4. store updated constraints
    std::vector<NConstraintInterface *> new_constraints;
    for (unsigned int i = 0; i < keep.size(); ++i)
      new_constraints.push_back(_constraints[keep[i]]);

    // return linearly independent ones
    _constraints.swap(new_constraints);
  }

  return result;

  // debug output
//  std::sort(val_keep.begin(),val_keep.end());
//  std::sort(val_remove.begin(),val_remove.end());
//  std::cerr << "KEEP VALUES" << std::endl;
//  for(auto v : val_keep)
//    std::cerr << std::scientific << v << std::endl;
//  std::cerr << "REMOVE VALUES" << std::endl;
//  for(auto v : val_remove)
//    std::cerr << std::scientific << v << std::endl;
}



#endif // COMISO_GMM_AVAILABLE

} // namespace ConstraintTools

//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_EIGEN3_AVAILABLE
//=============================================================================

