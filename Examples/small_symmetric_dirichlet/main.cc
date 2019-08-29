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
#include <CoMISo/Utils/StopWatch.hh>
#include <vector>
#include <CoMISo/NSolver/NewtonSolver.hh>
#include <CoMISo/NSolver/SymmetricDirichletProblem.hh>

//------------------------------------------------------------------------------------------------------

// Example main
int main(void)
{
  std::cout << "---------- 1) Get an instance of a SymmetricDirichletProblem..." << std::endl;

  // then create finite element problem and add sets
  unsigned int n_vertices = 4;
  COMISO::SymmetricDirichletProblem sd_problem(n_vertices);
  COMISO::SymmetricDirichletProblem::IndexVector indices(0,1,2);
  COMISO::SymmetricDirichletProblem::ReferencePositionVector2D positions;
  positions << 0, 0,
               1, 0,
               0, 1;
  sd_problem.add_triangle(indices, positions);
  COMISO::SymmetricDirichletProblem::IndexVector indices2(3,2,1);
  sd_problem.add_triangle(indices2, positions); // same reference positions can be used because we want both triangles to be isosceles

  std::vector<double> initial_solution{0.0,0.0,2,0.0,0.0,2.0,3.0,4.0};
  sd_problem.x() = initial_solution;

  auto f = sd_problem.eval_f(initial_solution.data());
  std::vector<double> grad(8);
  sd_problem.eval_gradient(initial_solution.data(), grad.data());

  std::cout << "energy: " << f << std::endl;
  std::cout << "grad: " << grad[0] << ", " << grad[1] << ", " << grad[2] << ", " << grad[3] << ", " << grad[4] << ", " << grad[5] << ", " << grad[6] << ", " << grad[7] << std::endl;

  std::cout << "---------- 2) Set up constraints..." << std::endl;
  // fix first vertex to origin to fix translational degree of freedom
  sd_problem.add_fix_point_constraint(0, 0.0, 0.0);
  // fix v coordinate of second vertex to 0 to fix rotational degree of freedom
  sd_problem.add_fix_coordinate_constraint(1, 1, 0.0);

  std::cout << "---------- 3) Solve with Newton Solver..." << std::endl;
  COMISO::SymmetricDirichletProblem::VectorD b;
  COMISO::SymmetricDirichletProblem::SMatrixD A;
  sd_problem.get_constraints(A, b);
  COMISO::NewtonSolver nsolver(1e-6,1e-9,20000);
  nsolver.solve(&sd_problem, A, b);

  // print result
  for(unsigned int i=0; i<sd_problem.x().size(); ++i)
    std::cerr << "x[" << i << "] = " << sd_problem.x()[i] << std::endl;

  return 0;
}

