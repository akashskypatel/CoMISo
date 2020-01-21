/*===========================================================================*\
 *                                                                           *
 *                               CoMISo                                      *
 *      Copyright (C) 2008-2019 by Computer Graphics Group, RWTH Aachen      *
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
#include <iostream>

#include <CoMISo/NSolver/NewtonSolver.hh>
#include <CoMISo/NSolver/NProblemInterface.hh>
#include <vector>

#include <CoMISo/Utils/ExactConstraintSatisfaction.hh>
//------------------------------------------------------------------------------------------------------

class SmallNProblem : public COMISO::NProblemInterface
{
public:
  // specify a function which has several local minima
  // f(x,y) = x^4 + y^4

  // number of unknown variables, here x and y = 2
  virtual int    n_unknowns   (                                )
  {
    return 2;
  }

  // initial value where the optimization should start from
  virtual void   initial_x    (       double* _x               )
  {
    _x[0] = 4.0;
    _x[1] = 2.0;
  }

  // function evaluation at location _x
  virtual double eval_f       ( const double* _x               )
  {
    return std::pow(_x[0], 4) + std::pow(_x[1], 4);
  }

  // gradient evaluation at location _x
  virtual void   eval_gradient( const double* _x, double*    _g)
  {
    _g[0] =  4.0*std::pow(_x[0], 3);
    _g[1] =  4.0*std::pow(_x[1], 3);
   }

  // hessian matrix evaluation at location _x
  virtual void   eval_hessian ( const double* _x, SMatrixNP& _H)
  {
    _H.resize(n_unknowns(), n_unknowns());
    _H.setZero();

    _H.coeffRef(0,0) =  12.0*std::pow(_x[0], 2);
    _H.coeffRef(1,0) =  0.0;
    _H.coeffRef(0,1) =  0.0;
    _H.coeffRef(1,1) =  12.0*std::pow(_x[1], 2);
  }

  // print result
  virtual void   store_result ( const double* _x               )
  {
    solution.resize(n_unknowns());
    for (int i = 0; i < n_unknowns(); ++i)
      solution[i] = _x[i];
  }

  // advanced properties
  virtual bool   constant_hessian() const { return false; }

  std::vector<double> solution;
};

// Example main
int main(void)
{
  std::cout << "---------- 1) Get an instance of a problem..." << std::endl;
  SmallNProblem problem;

  std::cout << "---------- 2) Set up constraints..." << std::endl;

/*
  int n_constraints = 2; // there will be one constraints
  Eigen::VectorXd b;
  Eigen::SparseMatrix<double> A(n_constraints, problem.n_unknowns());
  //A.resize(n_constraints, problem.n_unknowns());
  //A.setZero();
  b.resize(n_constraints);
  b.setZero();

  // first constraint: first variable equals three times second
  //different number of constraints :
  if(n_constraints == 1){
    A.coeffRef(0,0) =  2;
    A.coeffRef(0,1) = -6;
    b.coeffRef(0)   =  0;
  }else if(n_constraints == 3){
    A.coeffRef(0,0) =  2;
    A.coeffRef(0,1) = -6;
    A.coeffRef(1,0) =  5;
    A.coeffRef(1,1) =  0;
    A.coeffRef(2,0) =  0;
    A.coeffRef(2,1) = -3;
    b.coeffRef(0)   =  0;
    b.coeffRef(1)   =  15;
    b.coeffRef(2)   =  -3;
  }else if(n_constraints == 2){
    //A.coeffRef(0,0) =  2;
    //A.coeffRef(0,1) =  -1;
    b.coeffRef(0)   =  0;
    A.coeffRef(1,0) =  -2;
    //A.coeffRef(1,1) =  8;
    b.coeffRef(1)   =  21;
  }else if(n_constraints == 4){
    A.coeffRef(0,0) =  2;
    A.coeffRef(0,1) = -6;
    b.coeffRef(0)   =  0;
    A.coeffRef(1,0) =  2;
    A.coeffRef(1,1) = -6;
    b.coeffRef(1)   =  0;
    A.coeffRef(2,0) =  2;
    A.coeffRef(2,1) = -6;
    b.coeffRef(2)   =  0;
    A.coeffRef(3,0) =  2;
    A.coeffRef(3,1) = -6;
    b.coeffRef(3)   =  0;
  }
*/


  int n_constraints = 20;
  int rows = n_constraints;
  int cols = 100;
  double error = 0.00000000001;

  std::vector<Eigen::Triplet<int>> triplets;
  Eigen::SparseMatrix<int> A(rows, cols);
  Eigen::VectorXd x(cols);
  Eigen::VectorXi b(rows);

  //set triplets

  triplets.clear();
  for(int i = 0; i < rows; i++){
    if(i == 0){
      triplets.push_back({i,  0        , 1});
      triplets.push_back({i, (cols - 2), 1});
      triplets.push_back({i, (cols - 1), -1});
    }else{
      triplets.push_back({i,  0     ,  1});
      triplets.push_back({i,  i     ,  1});
      triplets.push_back({i, (i + 1), -1});
    }
  }
  A.setFromTriplets(triplets.begin(), triplets.end());

  //set x_vec 1,2,3,4...

  for(int i = 0; i < cols; i++){
    x.coeffRef(i) = i+1 + error;
  }

  //set b_vec = 0;

  b.setZero(rows);

/*

  std::cout << A << std::endl << b << std::endl;

  std::cout << "---------- 3) Solve with Newton Solver..." << std::endl;
  COMISO::NewtonSolver nsolver;
  nsolver.set_verbosity(15);
  nsolver.solve(&problem, A, b);

  std::cout << "---------- 4) Print solution..." << std::endl;
  std::cout << std::setprecision(100);
  for (unsigned int i = 0; i < problem.n_unknowns(); ++i)
    std::cout << "x[" << i << "] = " << problem.solution[i] << std::endl;


  std::cout << "---------- 5) Check constraint violation..." << std::endl;
  Eigen::VectorXd x;
  x.resize(problem.n_unknowns());
  for (unsigned int i = 0; i < problem.n_unknowns(); ++i)
    x.coeffRef(i) = problem.solution[i];

  std::cout << "Constraint violation: " << (A.cast<double>() *x - b.cast<double>()).squaredNorm() << std::endl;
*/

  std::cout << "Constraint violation: " << (A.cast<double>() *x - b.cast<double>()).squaredNorm() << std::endl;

  std::cout << "---------- 6) Try to exactly fulfill constraints..." << std::endl;

  ExactConstraintSatisfaction satisfy;
  satisfy.printMatrix(A);
  satisfy.evaluation(A, b, x);


  std::cout << "values of vector x : " << x << std::endl;

  std::cout << "---------- 7) Check constraint violation again..." << std::endl;

  std::cout << "Constraint violation: " << (A.cast<double>() *x - b.cast<double>()).squaredNorm() << std::endl;


  return 0;
}
