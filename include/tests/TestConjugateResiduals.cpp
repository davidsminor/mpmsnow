#include "tests/TestConjugateResiduals.h"

#include "MpmSim/ConjugateResiduals.h"

#include <iostream>

#include <Eigen/Eigenvalues>

using namespace MpmSim;
using namespace Eigen;

// class for applying a garden variety matrix to a vector:
class DenseMatrix : public ProceduralMatrix
{
public:
	DenseMatrix( const Eigen::MatrixXf& mat ) : m_mat( mat )
	{
	}

	virtual void multVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const
	{
		result = m_mat * x;
	}

private:

	const Eigen::MatrixXf& m_mat;

};

namespace MpmSimTest
{

void testConjugateResiduals()
{
	const int matrixSize = 12;

	// create random symmetric indefinite matrix:
	MatrixXf m = Eigen::MatrixXf::Random(matrixSize,matrixSize);
	MatrixXf mt = m.transpose();
	m += mt;
	
	// solver should converge in "matrixSize" steps, but it doesn't for numerical reasons
	// (I guess??), so we run it for twice that number of steps:
	ConjugateResiduals solver( 2 * matrixSize, 0 );
	VectorXf v = Eigen::VectorXf::Random(matrixSize);
	VectorXf result(matrixSize);
	
	// test on symmetric indefinite matrix:
	solver( DenseMatrix( m ), v, result );
	
	assert( ( m * result - v ).norm() < 1.e-4 );

}

}
