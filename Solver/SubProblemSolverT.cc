// Copyright 2022 Autodesk, Inc. All rights reserved.

#define COMISO_SUBSYSTEMSOLVERT_C

#include <CoMISo/Solver/SubProblemSolverT.hh>
#include <CoMISo/Utils/ProblemSubsetMapT.hh>


namespace COMISO
{

template <int DIM>
SubProblemSolverT<DIM>::SubProblemSolverT(
    size_t _max_var_num, const ValueVector& _fixed_values)
    : sbst_map_(_max_var_num), solver_()
{
  reset(_fixed_values);
}


template <int DIM>
void
SubProblemSolverT<DIM>::add_equation(LinearEquation _eq)
{
  sbst_map_.map(_eq);
  solver_.add_equation(_eq);
}

template <int DIM>
void
SubProblemSolverT<DIM>::add_constraint(LinearEquation _eq)
{
  sbst_map_.map(_eq);
  solver_.add_constraint(_eq);
}

template <int DIM>
void SubProblemSolverT<DIM>::add_constraints(std::vector<LinearEquation> _eqs)
{
  for (auto& eq : _eqs)
    add_constraint(std::move(eq));
}

template <int DIM>
void
SubProblemSolverT<DIM>::reset(ValueVector _fixed_values)
{
  sbst_map_.reset_fixed_values(std::move(_fixed_values));
}

template <int DIM>
void
SubProblemSolverT<DIM>::set_integers(IndexVector _int_var_indcs)
{
  sbst_map_.map(_int_var_indcs);
  solver_.set_integers(std::move(_int_var_indcs));
}

template <int DIM>
void
SubProblemSolverT<DIM>::solve(Result& _result)
{
  // Compute dense result
  typename MultiDimConstrainedSolverT<DIM>::PointVector res;
  solver_.solve(res);

  // Transform result into value vector with indices mapped back to original
  // problem
  const auto size = res.size();
  _result.clear();
  _result.reserve(size);
  for (size_t i = 0; i < size; ++i)
    _result.emplace_back(sbst_map_.mapped_back(i), res[i]);

  // reset map
  sbst_map_.reset();
}

template <int DIM>
const typename SubProblemSolverT<DIM>::ValueVector&
SubProblemSolverT<DIM>::fixed_values() const
{
  return sbst_map_.fixed_values();
}

}//namespace COMISO


