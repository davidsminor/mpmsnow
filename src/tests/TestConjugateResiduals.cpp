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
	DenseMatrix( const Eigen::MatrixXf& mat ) : m_mat( mat ), m_inv( m_mat.inverse() )
	{
	}

	virtual void multVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const
	{
		result = m_mat * x;
	}
	virtual void multInverseVector( const Eigen::VectorXf& x, Eigen::VectorXf& result ) const
	{
		result = m_inv * x;
	}
	virtual void subspaceProject( Eigen::VectorXf& x ) const
	{
	}

private:

	const Eigen::MatrixXf& m_mat;
	Eigen::MatrixXf m_inv;

};

namespace MpmSimTest
{

void testSolve()
{
	std::cerr << "testSolve()" << std::endl;
	const int matrixSize = 6;

	// create random symmetric indefinite matrix:
	MatrixXf A = Eigen::MatrixXf::Random(matrixSize,matrixSize);
	MatrixXf At = A.transpose();
	A += At;
	
	// solver should converge in "matrixSize" steps. "true" is for turning on logging:
	ConjugateResiduals solver( matrixSize, 0, 0, true );
	VectorXf v = Eigen::VectorXf::Random(matrixSize);
	VectorXf result(matrixSize);
	result.setZero();
	
	// test on symmetric indefinite matrix:
	solver( DenseMatrix( A ), v, result );
	std::cerr << A * result << std::endl << std::endl;
	std::cerr << v << std::endl << std::endl;
	// did that work? (Should the accuracy really suck as much as this? 1.e-3? Really?)
	std::cerr << ( A * result - v ).norm() << std::endl;
	assert( ( A * result - v ).norm() < 1.e-3 );

	// test residuals always decrease in norm:
	for( int i=1; i < matrixSize; ++i )
	{
		assert( solver.residuals[i].squaredNorm() < solver.residuals[i-1].squaredNorm() );
	}
	
	// test A orthogonality of residuals:
	for( int i=0; i < matrixSize; ++i )
	{
		for( int j=0; j < matrixSize; ++j )
		{
			VectorXf test( matrixSize );
			float f = fabs( solver.residuals[ i ].dot( A * solver.residuals[ j ] ) );
			if( i != j )
			{
				assert( f < 1.e-4 );
			}
		}
	}
	
	// test A^2 orthogonality of search directions:
	for( int i=0; i < matrixSize; ++i )
	{
		for( int j=0; j < matrixSize; ++j )
		{
			VectorXf test( matrixSize );
			float f = fabs( ( A * solver.searchDirections[ i ] ).dot( A * solver.searchDirections[ j ] ) );
			if( i != j )
			{
				assert( f < 1.e-4 );
			}
		}
	}

}

void testPreconditioner()
{
	std::cerr << "testPreconditioner()" << std::endl;
	const bool iKnowWhatImDoing( false );
	assert( iKnowWhatImDoing );
}


void testConjugateResiduals()
{
	std::cerr << "testConjugateResiduals()" << std::endl;
	testSolve();
	//testPreconditioner();
}

}
