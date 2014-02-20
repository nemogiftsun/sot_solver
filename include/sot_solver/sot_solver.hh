/*
 * Copyright 2014, Nirmal Giftsun, LAAS-CNRS
 *
 * This file is part of sot_solver.
 * sot_solver is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * sot_solver is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.  You should
 * have received a copy of the GNU Lesser General Public License along
 * with sot_solver.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef DG_SOT_SOLVER_HH
#define DG_SOT_SOLVER_HH


#include <dynamic-graph/command-setter.h>
#include <dynamic-graph/command-getter.h>
#include <dynamic-graph/factory.h>
#include <dynamic-graph/all-commands.h>
#include <dynamic-graph/entity.h>
#include <dynamic-graph/pool.h>
#include <dynamic-graph/signal-ptr.h>
#include <dynamic-graph/signal-time-dependent.h>
#include <dynamic-graph/linear-algebra.h>

#include <sot/core/exception-dynamic.hh>
#include <sot/core/matrix-homogeneous.hh>
#include <sot/core/flags.hh>

#include <Eigen/Householder>
#include <Eigen/QR>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;
using namespace dynamicgraph;
using namespace dynamicgraph::sot;


namespace dynamicgraph {
  namespace sot {
			namespace sotsolver{

class SotSolver : public Entity{
public: /* --- CONSTRUCTOR ---- */
      SotSolver(const std::string& inName);
      virtual ~SotSolver(void);
	  static const std::string CLASS_NAME;
	  virtual const std::string& getClassName( void ) const { return CLASS_NAME;}
      /// Header documentation of the python class
      virtual std::string getDocString () const {return "sot-collision\n";}
    /* --- SIGNALS --- */
    /* --- COMMANDS --- */
	  void testsvd(void);
	  MatrixXd testhouseholder(MatrixXd& A);
	  MatrixXd givensmatrix(MatrixXd& P,MatrixXd& m);
    void testgivensmatrix(void);
	  void solve(void);
    void testsolver(void);
    void createRandomPIMatrixOfRank(int desired_rank, int rows, int cols, MatrixXd& m);

	private: /* --- INTERNAL COMPUTATIONS --- */
	  std::vector< Eigen::MatrixXd > Ctasks;
	  Eigen::VectorXd dummy;

	};


#endif



























    }
  }
}
