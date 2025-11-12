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
#include <iostream>
#include <CoMISo/NSolver/FiniteElementTinyAD.hh>

#if COMISO_TINYAD_AVAILABLE

const Eigen::Index N = Eigen::Dynamic;
using TADG = TinyAD::Scalar<N,double,false>;
using TADH = TinyAD::Scalar<N,double,true>;

using VecVd = Eigen::Matrix<double,N,1>;


//------------------------------------------------------------------------------------------------------

template <class ScalarT>
ScalarT eval_f(const Eigen::Matrix<ScalarT,N,1>& _x)
{
  ScalarT t0 = 0.0;
  ScalarT t1 = 0.0;
  t0 += t1;

  for(Eigen::Index i=0; i<_x.size(); i+=2)
    t0 += (i+1)*_x[i];

  for(Eigen::Index i=1; i<_x.size(); i+=2)
    t1 += 2.0*_x[i];

  ScalarT f = pow(t0,4)+pow(t1,4);

  return f;
}


//------------------------------------------------------------------------------------------------------

// Example main
int main(void)
{
  std::cerr << "----------- Test TinyAD Dynamic ----------- " << std::endl;
  VecVd x;
  x.resize(3);
  x << 1, -2, 3;

  double fx = eval_f(x);

  std::cerr << "f(x) = " << fx << std::endl;

  auto x_ad_g = TADG::make_active(x);
  auto f_g = eval_f(x_ad_g);

  std::cerr << "grad f(x) = " << f_g.grad.transpose() << std::endl;

  auto x_ad_H = TADH::make_active(x);
  auto f_H = eval_f(x_ad_H);

  std::cerr << "Hess f(x) = " << std::endl << f_H.Hess << std::endl << std::endl;

  std::cerr << "----------- Test Runtime of TinyAD Dynamic----------- " << std::endl;
  COMISO::StopWatch sw;
  sw.start();
  Eigen::Matrix<double,N,N> H(3,3);
  H.setZero();


  int n_iters = 1000000;
  for(int i=1; i<n_iters+1; ++i)
  {
    VecVd xs = 1.0/double(i)*x;

    auto x_ad_H2 = TADH::make_active(xs);
    auto f_H2 = eval_f(x_ad_H2);

    H += f_H2.Hess;
  }

  auto t = sw.stop();

  std::cerr << "#terms  = " << n_iters << std::endl;
  std::cerr << "runtime = " << t/1000.0 << "s" << std::endl << std::endl;

  std::cerr << "accumulated Hess f(x) = " << std::endl
            << H << std::endl;

  return 0;
}

#else

int main(void)
{
  std::cerr << "TinyAD not available..." << std::endl;
}

#endif
