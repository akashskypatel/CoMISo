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
#include <CoMISo/Utils/StopWatch.hh>
#include <vector>
#include <CoMISo/NSolver/FiniteElementProblem.hh>
#include <CoMISo/NSolver/NPDerivativeChecker.hh>
#include <CoMISo/NSolver/IPOPTSolver.hh>
//#include <CoMISo/NSolver/AcceleratedQuadraticProxy.hh>


// create log barrier term for triangle orientation
class TriangleOrientationLogBarrierElement
{
public:

  // define dimensions
  const static int NV = 6;
  const static int NC = 1;

  typedef Eigen::Matrix<size_t,NV,1> VecI;
  typedef Eigen::Matrix<double,NV,1> VecV;
  typedef Eigen::Matrix<double,NC,1> VecC;
  typedef Eigen::Triplet<double> Triplet;

  inline double eval_f       (const VecV& _x, const VecC& _c) const
  {
    return -_c[0]*std::log((-_x[0]+_x[2])*(-_x[1]+_x[5])-(-_x[0]+_x[4])*(-_x[1]+_x[3]));
  }

  inline void   eval_gradient(const VecV& _x, const VecC& _c, VecV& _g) const
  {
    // get c*determinant
    double d = _c[0]/((-_x[0]+_x[2])*(-_x[1]+_x[5])-(-_x[0]+_x[4])*(-_x[1]+_x[3]));

    _g[0] = -d*(-_x[5]+_x[3]);
    _g[1] = -d*(-_x[2]+_x[4]);
    _g[2] = -d*(-_x[1]+_x[5]);
    _g[3] = -d*( _x[0]-_x[4]);
    _g[4] = -d*( _x[1]-_x[3]);
    _g[5] = -d*(-_x[0]+_x[2]);
  }

  inline void   eval_hessian (const VecV& _x, const VecC& _c, std::vector<Triplet>& _triplets) const
  {
    double d = 1.0/((-_x[0]+_x[2])*(-_x[1]+_x[5])-(-_x[0]+_x[4])*(-_x[1]+_x[3]));

    const VecV g;
    g[0] = d*(-_x[5]+_x[3]);
    g[1] = d*(-_x[2]+_x[4]);
    g[2] = d*(-_x[1]+_x[5]);
    g[3] = d*( _x[0]-_x[4]);
    g[4] = d*( _x[1]-_x[3]);
    g[5] = d*(-_x[0]+_x[2]);

    _triplets.clear();
    for(unsigned int i=0; i<NV; ++i)
      for(unsigned int j=0; j<NV; ++j)
      {
        double v = _c[0]*g[i]*g[j];
        _triplets.push_back(Triplet(i,j,v));
        _triplets.push_back(Triplet(j,i,v));
      }

    d *= _c[0];
    _triplets.push_back(Triplet(0,5, d));
    _triplets.push_back(Triplet(0,3,-d));
    _triplets.push_back(Triplet(5,0, d));
    _triplets.push_back(Triplet(3,0,-d));

    _triplets.push_back(Triplet(1,2, d));
    _triplets.push_back(Triplet(1,4,-d));
    _triplets.push_back(Triplet(2,1, d));
    _triplets.push_back(Triplet(4,1,-d));

    _triplets.push_back(Triplet(2,1, d));
    _triplets.push_back(Triplet(2,5,-d));
    _triplets.push_back(Triplet(1,2, d));
    _triplets.push_back(Triplet(5,2,-d));

    _triplets.push_back(Triplet(3,4, d));
    _triplets.push_back(Triplet(3,0,-d));
    _triplets.push_back(Triplet(4,3, d));
    _triplets.push_back(Triplet(0,3,-d));

    _triplets.push_back(Triplet(4,3, d));
    _triplets.push_back(Triplet(4,1,-d));
    _triplets.push_back(Triplet(3,4, d));
    _triplets.push_back(Triplet(1,4,-d));

    _triplets.push_back(Triplet(5,0, d));
    _triplets.push_back(Triplet(5,2,-d));
    _triplets.push_back(Triplet(0,5, d));
    _triplets.push_back(Triplet(2,5,-d));
  }
};

// create a simple finite element (x-c)^2
class SimpleElement2
{
public:

  // define dimensions
  const static int NV = 1;
  const static int NC = 1;

  typedef Eigen::Matrix<size_t,NV,1> VecI;
  typedef Eigen::Matrix<double,NV,1> VecV;
  typedef Eigen::Matrix<double,NC,1> VecC;
  typedef Eigen::Triplet<double> Triplet;

  inline double eval_f(const VecV& _x, const VecC& _c) const
  {
    return std::pow(_x[0]-_c[0],2);
  }

  inline void   eval_gradient(const VecV& _x, const VecC& _c, VecV& _g) const
  {
    _g[0] = 2.0*(_x[0]-_c[0]);
  }

  inline void   eval_hessian (const VecV& _x, const VecC& _c, std::vector<Triplet>& _triplets) const
  {
    _triplets.clear();
    _triplets.push_back(Triplet(0,0,2));
  }
};


//------------------------------------------------------------------------------------------------------

// Example main
int main(void)
{
  std::cout << "---------- 1) Get an instance of a FiniteElementProblem..." << std::endl;

  // first create sets of different finite elements
  COMISO::FiniteElementSet<TriangleOrientationLogBarrierElement> fe_set;

  TriangleOrientationLogBarrierElement::VecI tidx;
  tidx << 0,1,2,3,4,5;
  fe_set.instances().add_element(tidx, SimpleElement::VecC(1.0));

  // second set for boundary conditions
  COMISO::FiniteElementSet<SimpleElement2> fe_set2;
  fe_set2.instances().add_element(SimpleElement2::VecI(0), SimpleElement2VecV(0));
  fe_set2.instances().add_element(SimpleElement2::VecI(1), SimpleElement2VecV(0));
  fe_set2.instances().add_element(SimpleElement2::VecI(2), SimpleElement2VecV(2));
  fe_set2.instances().add_element(SimpleElement2::VecI(3), SimpleElement2VecV(0));
  fe_set2.instances().add_element(SimpleElement2::VecI(4), SimpleElement2VecV(1));
  fe_set2.instances().add_element(SimpleElement2::VecI(5), SimpleElement2VecV(1));

  // then create finite element problem and add sets
  COMISO::FiniteElementProblem fe_problem(3);
  fe_problem.add_set(&fe_set );
  fe_problem.add_set(&fe_set2);

  std::cout << "---------- 2) (optional for debugging) Check derivatives of problem..." << std::endl;
  COMISO::NPDerivativeChecker npd;
  npd.check_all(&fe_problem);

  // check if IPOPT solver available in current configuration
  #if( COMISO_IPOPT_AVAILABLE)
    std::cout << "---------- 3) Get IPOPT solver... " << std::endl;
    COMISO::IPOPTSolver ipsol;

    std::cout << "---------- 4) Solve..." << std::endl;
    // there are no constraints -> provide an empty vector
    std::vector<COMISO::NConstraintInterface*> constraints;
    ipsol.solve(&fe_problem, constraints);
  #endif

  // print result
  for(unsigned int i=0; i<fe_problem.x().size(); ++i)
    std::cerr << "x[" << i << "] = " << fe_problem.x()[i] << std::endl;

  return 0;
}

