#include "sot_solver/sot_solver.hh"

#include <dynamic-graph/all-commands.h>

#include <Eigen/Dense>
#include <Eigen/Jacobi>
#include <Eigen/QR>


using namespace std;
using namespace dynamicgraph::sot;
using namespace dynamicgraph::sot::sotsolver;
using namespace dynamicgraph;
using namespace Eigen;
using namespace timer;



DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN(SotSolver,"SotSolver");

SotSolver::SotSolver(const std::string& inName):Entity(inName)
{

			using namespace ::dynamicgraph::command;
      std::string docstring;
      docstring =     "\n"
    "    tests\n"
    "\n";

addCommand(std::string("testsolver"),
		     makeCommandVoid0(*this,&SotSolver::testsolver,docstring));
addCommand(std::string("testsomething"),
		     makeCommandVoid0(*this,&SotSolver::testsomething,docstring));

}
SotSolver::~SotSolver()
{
   
}

/* create partial isometry matrix*/
void SotSolver::createRandomPIMatrixOfRank(int desired_rank, int rows, int cols, MatrixXd& m)
{
  if(desired_rank == 0)
  {
    m.setZero(rows,cols);
    return;
  }
  if(desired_rank == 1)
  {
    // here we normalize the vectors to get a partial isometry
    m = VectorXd::Random(rows).normalized() * VectorXd::Random(cols).normalized().transpose();
    return;
  }

  MatrixXd a = MatrixXd::Random(rows,rows);
  MatrixXd d = MatrixXd::Identity(rows,cols);
  MatrixXd  b = MatrixXd::Random(cols,cols);

  // set the diagonal such that only desired_rank non-zero entries reamain
  int diag_size = (std::min)(d.rows(),d.cols());
  if(diag_size != desired_rank)
    d.diagonal().segment(desired_rank, diag_size-desired_rank) = VectorXd::Zero(diag_size-desired_rank);

  HouseholderQR<MatrixXd> qra(a);
  HouseholderQR<MatrixXd> qrb(b);
  m = qra.householderQ() * d * qrb.householderQ();
}


void SotSolver::decomposeLQ(MatrixXd &A,MatrixXd &L,MatrixXd &Q,MatrixXd &P,const string& type)
{

     if(type == "full")
     {
         // full pivoting
				 FullPivHouseholderQR<MatrixXd> qr2(A.rows(), A.cols());
				 qr2.compute(A.transpose());
				 L = qr2.matrixQR().triangularView<Upper>().transpose();
				 P = qr2.colsPermutation();
         Q = qr2.matrixQ().transpose();    
     }
     else
     {
                
         // Column pivoting Householder QR decomposition
				 ColPivHouseholderQR<MatrixXd> qr1(A.transpose().rows(), A.transpose().cols());
				 qr1.compute(A.transpose());
				 L = qr1.matrixQR().triangularView<Upper>().transpose();
				 P = qr1.colsPermutation();
				 Q = qr1.matrixQ();

     }        

}

void SotSolver::pinv(MatrixXd &L)
{
     JacobiSVD<MatrixXd> svd(L, ComputeThinU | ComputeThinV);
     MatrixXd D;D.setIdentity(svd.singularValues().rows(),svd.singularValues().rows());
     for(int k = 0; k < svd.singularValues().rows(); k++)
     {
      D.col(k)(k) = svd.singularValues().col(0)(k);
 
     }
          std::cout << "singular" << D<<endl;
     L= (svd.matrixV()*(D.inverse())*svd.matrixU().transpose());

}

void SotSolver::testsolver(void)
{
   // timer
   Timer timer;
 
   MatrixXd A;
   MatrixXd x = VectorXd::Random(50);
   A = MatrixXd::Random(6, 27) * MatrixXd::Random(27,50);
   MatrixXd b;
   b = A * x;
   std::cout << "Here is the right hand side b:\n" << b << endl;
   std::cout << "The b  svd is:\n";
   timer.start();
   std::cout << A*A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b) << endl;  
   timer.stop();
   std::cout <<"svd time is "<<timer.getElapsedTimeInMicroSec()<< endl;
   timer.start();
   solve(A,b);
   timer.stop();
   std::cout <<"the solver side b is "<<A*b<< endl;
   std::cout <<"cod solver time is "<<timer.getElapsedTimeInMicroSec()<< endl;

}
void SotSolver::testsomething(void)
{
		MatrixXd m = MatrixXd::Random(4,4);
   
}
/* solve lower triangular matrix by forward substitution */
void SotSolver::solveLowerTriangular(MatrixXd& A, VectorXd& b)
{
    const int nrows = A.rows();
    for (int i=0; i<nrows-1; i++)
    {
      b[i] = b[i]/A(i,i);
      b.tail(nrows-i-1) -= b[i]* A.col(i).tail(nrows-i-1);
    }
    b[nrows-1] /= A(nrows-1,nrows-1);
}

/* triangularize lower triangular matrix using givens transformation */
MatrixXd SotSolver::makeGivensTransformation(MatrixXd& m, const int rank)
{

  MatrixXd givens_transformation, intermediate_transformation;
  givens_transformation.setIdentity(m.rows(),m.rows());
  intermediate_transformation.setIdentity(m.rows(),m.rows());

  // compute orthogonal givens transformations     
   for(int i = rank ; i < m.rows(); i++)
   {
       for(int j = rank-1 ; j >= 0; j--)
      {           
          JacobiRotation<float> G;
          Vector2f v(m.col(j)(j),m.col(j)(i));
          G.makeGivens(v.x(), v.y());
          intermediate_transformation.setIdentity(m.rows(),m.rows());
          intermediate_transformation(j,j) = G.c() ;
          intermediate_transformation(i,i) = G.c() ;
          intermediate_transformation(j,i) = G.s();
          intermediate_transformation(i,j) = -G.s();    
          m.applyOnTheLeft(j, i, G.adjoint());
          givens_transformation = intermediate_transformation.transpose()*givens_transformation;
      }    
   }

  return givens_transformation;
  
}

void SotSolver::solve(MatrixXd &A, MatrixXd &b)
{

// do complete orthogonal decomposition 
  MatrixXd L,Q,P;
  // LQ decomposition
  // select column pivoting
  const string type = "column" ;
  decomposeLQ(A,L,Q,P,type);
  // check rank defeciency
	int rank = 0;
	const int rowsize = L.rows();
	const int colsize = L.cols();
	for(int k = 0; k < rowsize; k++)
	{
		 if(abs(L.col(k)(k)) > 0.00001)
		 (rank = rank+1);
		 else
		 break;
      
	}
  // triangularize L-Lower triangular matrix to match up with the rank of the system
  MatrixXd GT; 
  // rank not deficient
  if (rank == rowsize)
  {
     GT.setIdentity(rowsize,rowsize);
  }
  else //rank deficient
  {
     // find Givens transformation to triangularize L
     GT = makeGivensTransformation(L,rank);
  }
	MatrixXd q =  Q.block(0,0,colsize,rank);
	MatrixXd l = L.block(0,0,rank,rank);
	MatrixXd g = GT.block(0,0,rank,rowsize);

  VectorXd b_sol = g*P.transpose()*b;
  solveLowerTriangular(l, b_sol);
	b = q*b_sol;
  // solve L = Gt*Pt*b

  // 


}
