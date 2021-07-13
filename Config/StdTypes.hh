// (C) Copyright 2021 by Autodesk, Inc.

#ifndef STDTYPES_HH_INCLUDED
#define STDTYPES_HH_INCLUDED

#include <array>
#include <vector>

namespace COMISO_STD
{  

// Vector Types
typedef std::vector<double> DoubleVector;
typedef std::vector<int> IntVector;
typedef std::vector<unsigned int> UIntVector;

template <size_t DIM> struct SolverBaseT
{
  using Point = std::array<double, DIM>; // Object to minimize

  using PointVector = std::vector<Point>;

  struct Value // It means that X_name_ = val_.
               // Used to get the results and to set the fixed variables
  {
    size_t var_name; // Variable name.
    Point point;     // Value of the variable
    bool operator<(const Value& _vl) const { return var_name < _vl.var_name; }
    bool operator==(const Value& _vl) const { return var_name == _vl.var_name; }
  }; // struct Value

  struct LinearTerm // It means that coeff_ * X_name_
  {
    size_t var_name; // Variable name.
    double coeff;
    bool operator<(const LinearTerm& _lt) const
    {
      return var_name < _lt.var_name;
    }
  }; // struct LinearTerm

  using LinearTermVector = std::vector<LinearTerm>;

  struct LinearEquation // a linear equation in the form Sum(c_i * x_i) =
                        // constant_term
  {
    LinearTermVector linear_terms;
    // Construction must create an equation 'nothing' = 0, so const_term is
    // initialized with zeros.
    Point const_term{};
  }; // struct LinearEquation

}; // struct Types


}//namespace COMISO_STD

#endif//STDTYPES_HH_INCLUDED
