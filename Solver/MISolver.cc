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

#include <CoMISo/Config/config.hh>
#include "MISolver.hh"
#include "CoMISo/Utils/CoMISoError.hh"

#if (COMISO_QT_AVAILABLE)
#include <CoMISo/QtWidgets/MISolverDialogUI.hh>
#endif

#if COMISO_CPLEX_AVAILABLE
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
#endif

#if COMISO_GUROBI_AVAILABLE
#include <gurobi_c++.h>
#endif

//#define COMISO_MISOLVER_PERFORMANCE_TEST
#ifdef COMISO_MISOLVER_PERFORMANCE_TEST
#include "SparseQRSolver.hh"
#include "UMFPACKSolver.hh"
#include "EigenLDLTSolver.hh"
#endif

#if COMISO_SUITESPARSE_AVAILABLE
#include "CholmodSolver.hh"
#elif COMISO_EIGEN3_AVAILABLE
#include "EigenLDLTSolver.hh"
#else
#error "MISolver requires Suitesparse or Eigen3 support"
#endif

#include "IterativeSolverT.hh"

#include <CoMISo/Utils/gmm.hh>
#include <CoMISo/Utils/Tools.hh>

#include <Base/Debug/DebTime.hh>
#include <Base/Utils/StopWatch.hh>

#include <queue>
#include <float.h>

namespace COMISO
{
namespace
{
typedef unsigned int uint;

// gmm Column and ColumnIterator of CSC matrix
typedef gmm::linalg_traits<MISolver::CSCMatrix>::const_sub_col_type Col;
typedef gmm::linalg_traits<Col>::const_iterator ColIter;

using ToRoundSet = std::set<std::pair<double, uint>>;
using ToRoundSetIter = ToRoundSet::iterator;

// extend ToRoundSetIter with a flag to tell if it has been set or not
class ToRoundSetIterExt
{ 
public:
  void set(ToRoundSetIter _iter)
  {
    iter_ = _iter;
    null_ = false;
  }

  void clear() { null_ = true; }

  bool is_null() const { return null_; }

  const ToRoundSetIter& get() const
  {
    DEB_error_if(is_null(), "accessing null to-round iterator");
    return iter_;
  }

  ToRoundSetIter::pointer operator->() const { return get().operator->(); }
  ToRoundSetIter::reference operator*() const { return get().operator*(); }
  operator ToRoundSetIter() { return get(); }

private:
  ToRoundSetIter iter_;
  bool null_ = true;
};

} // namespace

// base class selected based on the available packages
class MISolver::DirectSolver : public
#if COMISO_SUITESPARSE_AVAILABLE
                               CholmodSolver
#elif COMISO_EIGEN3_AVAILABLE
                               EigenLDLTSolver
#else
#error "MISolver requires Suitesparse or Eigen3 support"
#endif
{
};

class MISolver::IterativeSolver : public IterativeSolverT<double>
{
};

// Constructor
MISolver::MISolver()
{// default parameters
  rounding_type_ = RoundingType::DEFAULT;

  initial_full_solution_ = true;
  iter_full_solution_ = true;
  final_full_solution_ = true;

  max_local_iters_ = 100000;
  max_local_error_ = 1e-3;
  max_cg_iters_ = 50;
  max_cg_error_ = 1e-3;

  multiple_rounding_threshold_ = 0.5;
  max_time_ = 60;

  direct_solver_ = new DirectSolver;
  iter_solver_ = new IterativeSolver;
}

MISolver::~MISolver()
{
  delete direct_solver_;
  delete iter_solver_;
}
//-----------------------------------------------------------------------------

void MISolver::solve(CSCMatrix& _A, Vecd& _x, Vecd& _rhs, const Veci& _to_round)
{
  DEB_enter_func;
  DEB_line(2, "integer variables #: " << _to_round.size());
  DEB_line(2, "continuous variables #: " << _x.size() - _to_round.size());

  // nothing to solve?
  if (gmm::mat_ncols(_A) == 0 || gmm::mat_nrows(_A) == 0)
    return;

  if (_to_round.empty())
    return solve_no_rounding(_A, _x, _rhs);

  switch (rounding_type_)
  {
  case RoundingType::NONE:
    return solve_no_rounding(_A, _x, _rhs);
  case RoundingType::DIRECT: 
    return solve_direct_rounding(_A, _x, _rhs, _to_round);
  case RoundingType::MULTIPLE: 
    return solve_multiple_rounding(_A, _x, _rhs, _to_round);
  case RoundingType::GUROBI:
    return solve_gurobi(_A, _x, _rhs, _to_round);
  case RoundingType::CPLEX:
    return solve_cplex(_A, _x, _rhs, _to_round);
  };
}

//-----------------------------------------------------------------------------

// inline function...
void MISolver::solve_cplex(
    CSCMatrix& _A, Vecd& _x, Vecd& _rhs, const Veci& _to_round)
{
  DEB_enter_func;
  DEB_out(2, "gurobi_max_time_: " << max_time_ << "\n");

#if COMISO_CPLEX_AVAILABLE

  // get round-indices in set
  std::set<int> to_round;
  for (unsigned int i = 0; i < _to_round.size(); ++i)
    to_round.insert(_to_round[i]);

  try
  {

    IloEnv env_;
    IloModel model(env_);

    // set time limite
    //    model.getEnv().set(GRB_DoubleParam_TimeLimit, gurobi_max_time_);

    unsigned int n = _rhs.size();

    // 1. allocate variables
    std::vector<IloNumVar> vars;
    for (unsigned int i = 0; i < n; ++i)
      if (to_round.count(i))
        vars.push_back(IloNumVar(env_, -IloIntMax, IloIntMax, IloNumVar::Int));
      else
        vars.push_back(
            IloNumVar(env_, -IloInfinity, IloInfinity, IloNumVar::Float));

    // 2. setup_energy

    // build objective function from linear system E = x^tAx - 2x^t*rhs
    IloExpr objective(env_);

    for (unsigned int i = 0; i < _A.nc; ++i)
      for (unsigned int j = _A.jc[i]; j < _A.jc[i + 1]; ++j)
      {
        objective += _A.pr[j] * vars[_A.ir[j]] * vars[i];
      }
    for (unsigned int i = 0; i < n; ++i)
      objective -= 2 * _rhs[i] * vars[i];

    // ToDo: objective correction!!!

    //    _A.jc[c+1]
    //    _A.pr[write]
    //    _A.ir[write]
    //    _A.nc
    //    _A.nr

    // minimize
    model.add(IloMinimize(env_, objective));

    // 4. solve
    IloCplex cplex(model);
    cplex.setParam(IloCplex::TiLim, max_time_);

#ifdef 0
    // set parameters comparable to CoMISo
    {
      cplex.setParam(IloCplex::MIPSearch  , 1);  // Traditional Branch-and-Cut
      cplex.setParam(IloCplex::NodeSel    , 0);  // Depth-First
      cplex.setParam(IloCplex::VarSel     , -1); // closest to integer
      cplex.setParam(IloCplex::MIPEmphasis, 1);  // concentrate on feasibility
    }
#endif

    cplex.solve();

    // 5. store result
    _x.resize(n);
    for (unsigned int i = 0; i < n; ++i)
      _x[i] = cplex.getValue(vars[i]);

    DEB_out(2, "CPLEX objective: " << cplex.getObjValue() << "\n");
  }
  catch (IloException& e)
  {
    PROGRESS_RESUME_ABORT; // resume a processed abort request
    DEB_warning(2, "CPLEX Concert exception caught: " << e.getMessage())
  }
  catch (...)
  {
    PROGRESS_RESUME_ABORT; // resume a processed abort request
    DEB_warning(1, "CPLEX Unknown exception caught")
  }

#else
  DEB_warning(1, "CPLEX solver is not available, please install it")
#endif
}

//-----------------------------------------------------------------------------

void MISolver::solve_no_rounding(CSCMatrix& _A, Vecd& _x, Vecd& _rhs)
{
  COMISO_THROW_if(
      !direct_solver_->calc_system_gmm(_A), UNSPECIFIED_EIGEN_FAILURE);
  COMISO_THROW_if(!direct_solver_->solve(_x, _rhs), UNSPECIFIED_EIGEN_FAILURE);
}

//-----------------------------------------------------------------------------

void MISolver::resolve(Vecd& _x, Vecd& _rhs)
{
  direct_solver_->solve(_x, _rhs);
}

//-----------------------------------------------------------------------------

void MISolver::solve_direct_rounding(
    CSCMatrix& _A, Vecd& _x, Vecd& _rhs, const Veci& _to_round)
{
  DEB_enter_func;
  Veci to_round(_to_round);
  // copy to round vector and make it unique
  std::sort(to_round.begin(), to_round.end());
  Veci::iterator last_unique;
  last_unique = std::unique(to_round.begin(), to_round.end());
  size_t r = (size_t)(last_unique - to_round.begin());
  to_round.resize(r);

  // initalize old indices
  Veci old_idx(_rhs.size());
  for (unsigned int i = 0; i < old_idx.size(); ++i)
    old_idx[i] = i;
  direct_solver_->calc_system_gmm(_A);
  direct_solver_->solve(_x, _rhs);

#ifdef COMISO_MISOLVER_PERFORMANCE_TEST
  // check solver performance (only for testing!!!)
  {
    Base::StopWatch sw;

    // performance comparison code
#if (COMISO_SUITESPARSE_SPQR_AVAILABLE)
    {
      sw.start();
      COMISO::SparseQRSolver spqr;
      spqr.calc_system_gmm(_A);
      std::cerr << "SparseQR factor took: " << sw.stop() / 1000.0 << "s\n";
      Vecd x2(_x);
      sw.start();
      spqr.solve(x2, _rhs);
      std::cerr << "SparseQR solve took: " << sw.stop() / 1000.0 << "s\n";
      Vecd res(_x);
      gmm::add(_x, gmm::scaled(x2, -1.0), res);
      std::cerr << "DIFFERENCE IN RESULT: " << gmm::vect_norm2(res)
                << std::endl;
    }
#endif

    // performance comparison code
#if (COMISO_SUITESPARSE_AVAILABLE)
    {
      sw.start();
      COMISO::UMFPACKSolver umf;
      umf.calc_system_gmm(_A);
      std::cerr << "UMFPack factor took: " << sw.stop() / 1000.0 << "s\n";
      Vecd x3(_x);
      sw.start();
      umf.solve(x3, _rhs);
      std::cerr << "UMFPack solve took: " << sw.stop() / 1000.0 << "s\n";
      Vecd res2(_x);
      gmm::add(_x, gmm::scaled(x3, -1.0), res2);
      std::cerr << "UMFPACK DIFFERENCE IN RESULT: " << gmm::vect_norm2(res2)
                << std::endl;
    }

    // performance comparison code
    {
      sw.start();
      COMISO::CholmodSolver chol;
      chol.calc_system_gmm(_A);
      std::cerr << "Choldmod factor took: " << sw.stop() / 1000.0 << "s\n";
      Vecd x4(_x);
      sw.start();
      chol.solve(x4, _rhs);
      std::cerr << "Choldmod solve took: " << sw.stop() / 1000.0 << "s\n";
      Vecd res(_x);
      gmm::add(_x, gmm::scaled(x4, -1.0), res);
      std::cerr << "DIFFERENCE IN RESULT: " << gmm::vect_norm2(res)
                << std::endl;
    }
#endif

#if (COMISO_EIGEN3_AVAILABLE)
    // performance comparison code
    {
      sw.start();
      COMISO::EigenLDLTSolver ldlt;
      ldlt.calc_system_gmm(_A);
      DEB_out(2, "Eigen LDLT factor took: " << sw.stop() / 1000.0 << "s\n")
          Vecd x5(_x);
      sw.start();
      ldlt.solve(x5, _rhs);
      DEB_out(2, "Eigen LDLT solve took: " << sw.stop() / 1000.0 << "s\n")
          Vecd res(_x);
      gmm::add(_x, gmm::scaled(x5, -1.0), res);
      DEB_warning(2, "DIFFERENCE IN RESULT: " << gmm::vect_norm2(res))
    }
#endif
  }
#endif

  // round and eliminate variables
  Vecui elim_i;
  Vecd elim_v;
  for (unsigned int i = 0; i < to_round.size(); ++i)
  {
    _x[to_round[i]] = double_round(_x[to_round[i]]);
    elim_i.push_back(to_round[i]);
    elim_v.push_back(_x[to_round[i]]);
    // update old idx
    old_idx[to_round[i]] = -1;
  }

  Veci::iterator new_end = std::remove(old_idx.begin(), old_idx.end(), -1);
  old_idx.resize(new_end - old_idx.begin());
  // eliminate vars from linear system
  Vecd xr(_x);
  COMISO_GMM::eliminate_csc_vars2(elim_i, elim_v, _A, xr, _rhs);

  // std::cerr << "size A: " << gmm::mat_nrows(_A) << " " << gmm::mat_ncols(_A)
  // 	    << std::endl;
  // std::cerr << "size x  : " << xr.size() << std::endl;
  // std::cerr << "size rhs: " << _rhs.size() << std::endl;

  // final full solution
  if (gmm::mat_ncols(_A) > 0)
  {
    //    direct_solver_->update_system_gmm(_A);
    direct_solver_->calc_system_gmm(_A);
    direct_solver_->solve(xr, _rhs);
  }

  // store solution values to result vector
  for (unsigned int i = 0; i < old_idx.size(); ++i)
    _x[old_idx[i]] = xr[i];
}

//-----------------------------------------------------------------------------

bool MISolver::update_solution_is_local(
    const CSCMatrix& _A, Vecd& _x, const Vecd& _rhs, const Vecui& _neigh_i)
{
  DEB_enter_func;
  bool converged = false; // set to not converged

  if (max_local_iters_ > 0) // compute new solution
  {
    DEB_out(11, "use local iteration ");
    converged = iter_solver_->gauss_seidel_local(
        _A, _x, _rhs, _neigh_i, max_local_iters_, max_local_error_);

    ++n_local_;
  }

  if (converged)
  {
    DEB_line(11, "Local iteration converged");
    return true;
  }

  bool only_local_updates = true; // set to false if we do any global updates

  DEB_line(6, "Local iterations failed to converge, needs global update");
  if (max_cg_iters_ > 0) 
  {// conjugate gradient iterations
    DEB_out(6, ", cg ");

    int max_cg_iters = max_cg_iters_;
    double tolerance = max_cg_error_;
    converged =
        iter_solver_->conjugate_gradient(_A, _x, _rhs, max_cg_iters, tolerance);

    DEB_out(6, "( converged " << converged << " "
                              << " iters " << max_cg_iters << " "
                              << " res_norm " << tolerance << "\n");
    only_local_updates = false;
    ++n_cg_;
  }

  if (!converged && iter_full_solution_ && gmm::mat_ncols(_A) > 0)
  {
    DEB_out(6, ", full ");
    if (cholmod_step_done_)
      direct_solver_->update_system_gmm(_A);
    else
    {
      direct_solver_->calc_system_gmm(_A);
      cholmod_step_done_ = true;
    }
    // const_cast<> to work around broken const correctness in Eigen
    direct_solver_->solve(_x, const_cast<Vecd&>(_rhs)); 

    only_local_updates = false;
    ++n_full_;
  }

  DEB_line(6, "");
  return only_local_updates;
}

//-----------------------------------------------------------------------------

void MISolver::solve_multiple_rounding(
    CSCMatrix& _A, Vecd& _x, Vecd& _rhs, const Veci& _to_round)
{
  DEB_time_func_def;
  // some statistics
  n_local_ = 0;
  n_cg_ = 0;
  n_full_ = 0;

  // reset cholmod step flag
  cholmod_step_done_ = false;

  if (initial_full_solution_)
  {
    DEB_line(2, "initial full solution");
    // TODO: we can throw more specific outcomes in the body of the
    // functions below
    COMISO_THROW_if(
        !direct_solver_->calc_system_gmm(_A), UNSPECIFIED_EIGEN_FAILURE);
    COMISO_THROW_if(!direct_solver_->solve(_x, _rhs), UNSPECIFIED_EIGEN_FAILURE);

    cholmod_step_done_ = true;

    ++n_full_;
  }

  DEB_warning_if(max_local_iters_ == 0, 1,
      "No local iterations available, this method will not perform well");

  // use a doubly-indexed data structure to be able to access a "to_round"
  // variable by both its index and round residue
  ToRoundSet to_round; // variables yet to round sorted by round residue 
  // store iterators to the to_round variable in the residue sorted container
  std::vector<ToRoundSetIterExt> to_round_indx(_x.size());

  const auto add_to_round_index = 
      [&_x, &to_round, &to_round_indx](const uint _tr_indx)
  {
    const auto x_rr = round_residue(_x[_tr_indx]);
    to_round_indx[_tr_indx].set(to_round.emplace(x_rr, _tr_indx).first);
  };

  for (const auto tr_indx : _to_round)
  {// initialize the doubly-indexed data
    if (to_round_indx[tr_indx].is_null()) // not added already?
      add_to_round_index(tr_indx);
  }

  Vecui neigh_i; // neighbors for local optimization
  while (!to_round.empty()) // loop until solution computed
  {
    DEB_line(11, "Integer DOF's left: " << to_round.size());
    DEB_line(11, "residuum_norm: " << COMISO_GMM::residuum_norm(_A, _x, _rhs));

    neigh_i.clear(); // clear neigh for local update

    double rnd_err_sum = 0;
    for (auto tr_it = to_round.begin(); // start at the front
         tr_it != to_round.end() && // if there are still variables to round ..
         (rnd_err_sum == 0 || // ... at least one variable should be rounded
             rnd_err_sum + tr_it->first <= multiple_rounding_threshold_);
         tr_it = to_round.erase(tr_it))
    {
      rnd_err_sum += tr_it->first;
      const auto tr_indx = tr_it->second;
      const auto rnd_x = double_round(_x[tr_indx]); // store rounded value
      to_round_indx[tr_indx].clear(); // clear pointer

      // compute neighbors
      const Col col = gmm::mat_const_col(_A, tr_indx);
      for (auto it = gmm::vect_const_begin(col), ite = gmm::vect_const_end(col);
           it != ite; ++it)
      {
        if (it.index() != tr_indx)
          neigh_i.push_back(static_cast<int>(it.index()));
      }
      // eliminate x_i from _A and _rhs, and set _x[tr_indx] = rnd_x
      COMISO_GMM::fix_var_csc_symmetric(tr_indx, rnd_x, _A, _x, _rhs);
    }

    // 3-stage solution update w.r.t. roundings: local GS / CG / SparseCholesky
    if (update_solution_is_local(_A, _x, _rhs, neigh_i)) // only local updates?
    {// re-sort only the updated variables 
      for (const auto updt_indx : iter_solver_->updated_variable_indices())
      { // refresh the rounding residues for the updated variables only
        if (to_round_indx[updt_indx].is_null())
          continue; // _x[updt_indx] does not need to be rounded
        to_round.erase(to_round_indx[updt_indx]); // erase from the set
        add_to_round_index(updt_indx);
      }
      continue;
    }

    // global iteration(s) performed, so update all
    ToRoundSet to_round_temp; // temporary storage for the to_round entries
    std::swap(to_round_temp, to_round); // swap() is fast, and clears to_round
    for (const auto& tr : to_round_temp) // re-add all variables still to round
      add_to_round_index(tr.second);
  }

  if (final_full_solution_ && gmm::mat_ncols(_A) > 0) // final full solution?
  {
    DEB_line(3, "final full solution");
    if (cholmod_step_done_)
      direct_solver_->update_system_gmm(_A);
    else
      direct_solver_->calc_system_gmm(_A);

    direct_solver_->solve(_x, _rhs);
    ++n_full_;
  }

  // output statistics
  DEB_line(2, " *** Statistics of MiSo Solver ***");
  DEB_line(2, "Number of LOCAL iterations  = " << n_local_);
  DEB_line(2, "Number of CG    iterations  = " << n_cg_);
  DEB_line(2, "Number of FULL  iterations  = " << n_full_);
  DEB_line(2, "Number of ROUNDING          = " << _to_round.size());
}

//-----------------------------------------------------------------------------

void MISolver::solve_gurobi(
    CSCMatrix& _A, Vecd& _x, Vecd& _rhs, const Veci& _to_round)
{
  DEB_enter_func;
#if COMISO_GUROBI_AVAILABLE

  // get round-indices in set
  std::set<int> to_round;
  for (unsigned int i = 0; i < _to_round.size(); ++i)
    to_round.insert(_to_round[i]);

  try
  {
    GRBEnv env = GRBEnv();

    GRBModel model = GRBModel(env);

    // set time limite
    model.getEnv().set(GRB_DoubleParam_TimeLimit, max_time_);

    unsigned int n = _rhs.size();

    // 1. allocate variables
    std::vector<GRBVar> vars;
    for (unsigned int i = 0; i < n; ++i)
      if (to_round.count(i))
        vars.push_back(
            model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER));
      else
        vars.push_back(
            model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS));

    // Integrate new variables
    model.update();

    // 2. setup_energy

    // build objective function from linear system E = x^tAx - 2x^t*rhs
    GRBQuadExpr objective;

    for (unsigned int i = 0; i < _A.nc; ++i)
      for (unsigned int j = _A.jc[i]; j < _A.jc[i + 1]; ++j)
      {
        objective += _A.pr[j] * vars[_A.ir[j]] * vars[i];
      }
    for (unsigned int i = 0; i < n; ++i)
      objective -= 2 * _rhs[i] * vars[i];

    //    _A.jc[c+1]
    //    _A.pr[write]
    //    _A.ir[write]
    //    _A.nc
    //    _A.nr

    // minimize
    model.set(GRB_IntAttr_ModelSense, 1);
    model.setObjective(objective);

    // 4. solve
    model.optimize();

    // 5. store result
    _x.resize(n);
    for (unsigned int i = 0; i < n; ++i)
      _x[i] = vars[i].get(GRB_DoubleAttr_X);

    DEB_out(
        2, "GUROBI objective: " << model.get(GRB_DoubleAttr_ObjVal) << "\n");
  }
  catch (GRBException& e)
  {
    PROGRESS_RESUME_ABORT; // resume a processed abort request
    DEB_warning(2,
        "Error code = " << e.getErrorCode() << "[" << e.getMessage() << "]\n");
  }
  catch (...)
  {
    PROGRESS_RESUME_ABORT; // resume a processed abort request
    DEB_warning(1, "Exception during optimization");
  }

#else
  DEB_warning(1, "GUROBI solver is not available, please install it");
#endif
}

//----------------------------------------------------------------------------

void MISolver::show_options_dialog()
{
  DEB_enter_func;
#if (COMISO_QT_AVAILABLE)
  MISolverDialog* pd = new MISolverDialog(*this);
  pd->exec();
#else
  DEB_warning(1, "Qt not available to show solver dialog!!!");
#endif
}

// end namespace COMISO
} // namespace COMISO
