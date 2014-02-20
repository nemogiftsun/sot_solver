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
addCommand(std::string("testgivensmatrix"),
		     makeCommandVoid0(*this,&SotSolver::testgivensmatrix,docstring));
addCommand(std::string("solve"),
		     makeCommandVoid0(*this,&SotSolver::solve,docstring));

}
SotSolver::~SotSolver()
{
   
}

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


MatrixXd SotSolver::testhouseholder(MatrixXd &A)
{
     MatrixXd R;

     std::cout << "------------------------------------------------\n";
     std::cout << "------Normal one--------------------------------\n";
     /*
     // normal one 
     HouseholderQR<MatrixXd> qr(A.rows(), A.cols());
     qr.compute(A);
     std::cout << "The matrix q is:\n";
     MatrixXd hq = qr.householderQ();
     R = qr.matrixQR().triangularView<Upper>();
     std::cout << (hq) << endl;
     std::cout << "The matrix r is:\n";
     std::cout << R << endl;
     std::cout << "The matrix product is:\n";
     std::cout << (qr.householderQ()* R)<< endl;*/
     std::cout << "-----------------------------------------------------\n";
     std::cout << "------Column pivoting--------------------------------\n";
     // column pivoting
     ColPivHouseholderQR<MatrixXd> qr1(A.transpose().rows(), A.transpose().cols());
     qr1.compute(A.transpose());
     std::cout << "The matrix q is:\n";
     MatrixXd hq = qr1.matrixQ();
     std::cout << (hq) << endl;
     std::cout << "The matrix qr is:\n";
     R = qr1.matrixQR().triangularView<Upper>();
     std::cout << (R) << endl;
     std::cout << "The matrix p is:\n";
     MatrixXd P = qr1.colsPermutation();
     std::cout << (P) << endl;
     //std::cout << "The matrix product is:\n";
     //std::cout << (qr1.matrixQ()* R*P.inverse())<< endl;
     std::cout << "------------------------------------------------\n";
     std::cout << "------Full pivoting--------------------------------\n";
     // full pivoting
     /*
     FullPivHouseholderQR<MatrixXd> qr2(A.rows(), A.cols());
     qr2.compute(A.transpose());
     std::cout << "The matrix q is:\n";
     std::cout << (qr2.matrixQ()) << endl;
     std::cout << "The matrix qr is:\n";
     R = qr1.matrixQR().triangularView<Upper>();
     std::cout << (R) << endl;
     std::cout << "------------------------------------------------\n";
     std::cout << "The matrix p is:\n";
      P = qr1.colsPermutation();
     std::cout << (P) << endl;
     std::cout << "The matrix A is:\n";
     std::cout << A << endl;
     */
     
     MatrixXd L = R.transpose();
     MatrixXd Q = qr1.matrixQ();
     MatrixXd Qt = qr1.matrixQ().transpose();
     std::cout << "The matrix L is:\n";
     std::cout << L<< endl;
     std::cout << "The matrix Qt is:\n";
  int rank = 0;
  const int rowsize = L.rows();
  const int colsize = L.cols();
  for(int k = 0; k < colsize; k++)
  {
     if(abs(L.col(k)(k)) > 0.00001)
     (rank = rank+1);
     else
     break;
  }
     std::cout << Qt<< endl;
     std::cout << "The matrix pinv is:\n";
     std::cout << P.inverse()<< endl;
     std::cout << "The matrix A is:\n";
     std::cout << A << endl;
     std::cout << "The matrix product is:\n";
     std::cout << (P*L*Qt)<< endl;
     MatrixXd GT = givensmatrix(P,L);
     std::cout << "The matrix product is:\n";
     std::cout << (P*GT*L*Qt)<< endl;
     //GT = (GT);
     
     
     MatrixXd inverse;
     MatrixXd q =  Qt.block(0,0,rank,colsize);
     MatrixXd l = L.block(0,0,rank,rank);
     MatrixXd g = GT.block(0,0,rowsize,rank);
     std::cout << "The matrix prod is is:\n"<<endl;
     std::cout << g*l*q<<endl;
     inverse = q.transpose()*l.inverse()*g.transpose();
     std::cout << "orthonormality 1:\n";
     std::cout << g*g.transpose()<< endl;
     std::cout << "orthonormality 2:\n";
     std::cout << q*q.transpose()<< endl;
     //std::cout << (Q*Qt*(P*GT).transpose())<< endl;
     return inverse;
 

}

void SotSolver::testsolver(void)
{

   MatrixXd A = MatrixXd::Random(6, 8);
     //Index rows = internal::random<Index>(20,200), cols = internal::random<int>(20,200);
     //Index rank = internal::random<Index>(1, (std::min)(rows, cols)-1);
     //createRandomPIMatrixOfRank(3,6,8,A);

   A = MatrixXd::Random(6, 4) * MatrixXd::Random(4,8);
   //MatrixXf A = MatrixXf::Random(3, 2);
   VectorXd b = VectorXd::Random(6);
   //std::cout << "The least-squares solution using svd is:\n";
   //std::cout << A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b) << endl;  
   
   std::cout << A*(testhouseholder(A)*b);
   std::cout << "Here is the right hand side b:\n" << b << endl;


}
void SotSolver::testgivensmatrix(void)
{
MatrixXd m = MatrixXd::Random(4,4);
std::cout << "Here is the matrix m:\n" << m << endl;
JacobiRotation<double> J;
//cout << "Here is the m elements" << endl << m(0,1) << endl;
Vector2f v(m.col(0)(2),m.col(0)(3));
J.makeGivens(v.x(), v.y());
//J.makeJacobi(m.col(0)(2),m.col(0)(3));
//J.makeJacobi(m,2,3);
cout << "Here is the J.m_c" << endl << J.c() << endl;
cout << "Here is the J.m_s" << endl << J.s() << endl;
m.applyOnTheLeft(2, 3, J.adjoint());
//m.applyOnTheRight(0, 1, J);
cout << "Here is the m" << endl << m << endl;

   

}
MatrixXd SotSolver::givensmatrix(MatrixXd& P,MatrixXd& m)
{
  int rowsize = m.rows();
  int colsize = m.cols();
  MatrixXd L = m;
  
  MatrixXd temp_transformation;
  
  temp_transformation.setIdentity(rowsize,rowsize);
  int rank = 0;
  for(int k = 0; k < colsize; k++)
  {
     if(abs(m.col(k)(k)) > 0.00001)
     (rank = rank+1);
     else
     break;
  }
  //cout << "Here is the L matrix" << endl << L << endl;
  //std::cout << "rank of the matrix:\n" << rank << endl;

  for(int j = rank-1 ; j >= 0; j--)
  {
     for(int i = rowsize -1; i >= rank; i--)
     {
        MatrixXd givens_transformation = P;
        JacobiRotation<float> G;
        Vector2f v(m.col(j)(j),m.col(j)(i));
        G.makeGivens(v.x(), v.y());
        temp_transformation.setIdentity(rowsize,rowsize);
        temp_transformation(j,j) = G.c() ;
        temp_transformation(i,i) = G.c() ;
        temp_transformation(j,i) = G.s();
        temp_transformation(i,j) = -G.s();    

        m.applyOnTheLeft(j, i, G.adjoint());
        //cout << "Here is the givenstransformedmatrix :" << endl << m << endl;
        //givens_transformation = temp_transformation.transpose()*givens_transformation;
        givens_transformation = givens_transformation*temp_transformation.transpose();
      
     }
     
  }
  m = givens_transformation * L;
  cout << "Here is the GT * L" << endl << m << endl;
  return givens_transformation.transpose();
  
}

void SotSolver::solve(void)
{


}
