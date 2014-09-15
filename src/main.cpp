#ifdef WIN32
#include <windows.h>
#endif

#include <algorithm>
#include <vector>
#include <iostream>

#include <GL/gl.h>
#include <GL/glut.h>

#include "MpmSim/Sim.h"
#include "MpmSim/SquareMagnitudeTermination.h"
#include "MpmSim/CollisionPlane.h"
#include "MpmSim/GravityField.h"
#include "MpmSim/CubicBsplineShapeFunction.h"
#include "MpmSim/SnowConstitutiveModel.h"
#include "MpmSim/ConjugateResiduals.h"

using namespace Eigen;
using namespace MpmSim;

////////////////////////////////////////////////////////////////////////////////
// constants
const unsigned int window_width = 512;
const unsigned int window_height = 512;

float g_x(0);
float g_y(0);
float g_r(0.5);

// GL functionality
bool initGL(int *argc, char** argv);

// rendering callbacks
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);

int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;

const float g_timeStep = 0.0025f;
const float g_gridSize = 0.025f;

std::auto_ptr< Sim > g_sim;

CubicBsplineShapeFunction g_shapeFunction;
SnowConstitutiveModel g_snowModel(
		1.4e5f, // young's modulus
		0.2f, // poisson ratio
		10, // hardening
		2.5e-2f, // compressive strength
		7.5e-3f	// tensile strength
);

// collision objects:
Sim::CollisionObjectSet g_collisionObjects;
Sim::ForceFieldSet g_fields;

int g_time(0);

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv)
{
	// collisions and fields:
	Eigen::Vector4f plane( 0, 1, 0, 0.5f );
	Eigen::Vector4f plane2( 0, -1, 0, -1.5f );
	Eigen::Vector3f gravity( 0, -9.8, 0 );
	g_collisionObjects.objects.push_back( new CollisionPlane( plane, 0.7f, true ) );
	//g_collisionObjects.objects.push_back( new CollisionPlane( plane2, 0.7f, true ) );
	//((CollisionPlane*)g_collisionObjects.objects.back())->setV( Eigen::Vector3f(0,-1,0) );
	g_fields.fields.push_back( new GravityField( gravity ) );
	
	// initial configuration:
	float particleSpacing = 0.5f * g_gridSize;
	float initialDensity = 400;

	std::vector<Eigen::Vector3f> x;
	std::vector<float> m;
	for( int i=0; i < 40; ++i )
	{
		for( int j=0; j < 40; ++j )
		{

			float xr = particleSpacing * ( i - 20 + 0.5f );
			float yr = particleSpacing * ( j - 20 + 0.5f );
			Vector3f pos( xr, yr, 0 );
			if( pos.norm() < 0.2f )
			{
				x.push_back( pos );
				m.push_back( initialDensity * particleSpacing * particleSpacing * particleSpacing );
			}
		}
	}
	g_sim.reset(
		new MpmSim::Sim(
			x, m, g_gridSize, g_shapeFunction, g_snowModel, g_collisionObjects, g_fields, 2
		)
	);

	initGL( &argc, argv );
	
	// start rendering mainloop
	glutMainLoop();

}

////////////////////////////////////////////////////////////////////////////////
//! Initialize GL
////////////////////////////////////////////////////////////////////////////////
bool initGL(int *argc, char **argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(window_width, window_height);
	glutCreateWindow("Thnowing on our roatht!");
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);

	// default initialization
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glEnable( GL_CULL_FACE );
	glEnable( GL_DEPTH_TEST );

	// viewport
	glViewport(0, 0, window_width, window_height);

	// projection
	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();

	return true;
}

/*
class GridDrawWotsit : public LinearSolver::Debug
{
public:
	GridDrawWotsit( Grid& g, ParticleData& d, float timeStep, const std::vector< CollisionObject* >& collisionObjects ) :
		m_g( g ),
		m_d( d ),
		m_timeStep( timeStep ),
		m_collisionObjects( collisionObjects )
	{
	}
	
	Grid& m_g;
	ParticleData& m_d;
	float m_timeStep;
	const std::vector< CollisionObject* >& m_collisionObjects;

	virtual void operator()( Eigen::VectorXf& velocities )
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glDisable( GL_DEPTH_TEST );
		
		float maxMass(0);
		for( int i=0; i < m_g.masses().size(); ++i )
		{
			if( m_g.masses()[i] > maxMass ) maxMass = m_g.masses()[i];
		}

		glBegin( GL_QUADS );
		for( int i=0; i < m_g.dimensions()[0]; ++i )
		{
			for( int j=0; j < m_g.dimensions()[1]; ++j )
			{
				int index = m_g.coordsToIndex( Vector3i( i, j, m_g.dimensions()[2]/2 ) );
				float mass = m_g.masses()[ index ] / maxMass;
				glColor3f( mass, mass, mass);
				glVertex3f( m_g.minCoord()[0] + (i-0.5) * m_g.gridH(), m_g.minCoord()[1] + (j-0.5) * m_g.gridH(), 0 );
				glVertex3f( m_g.minCoord()[0] + (i+0.5) * m_g.gridH(), m_g.minCoord()[1] + (j-0.5) * m_g.gridH(), 0 );
				glVertex3f( m_g.minCoord()[0] + (i+0.5) * m_g.gridH(), m_g.minCoord()[1] + (j+0.5) * m_g.gridH(), 0 );
				glVertex3f( m_g.minCoord()[0] + (i-0.5) * m_g.gridH(), m_g.minCoord()[1] + (j+0.5) * m_g.gridH(), 0 );
			}
		}
		glEnd();

		for( size_t p = 0; p < m_d.particleX.size(); ++p )
		{
			float r = 2 * pow( m_d.particleVolumes[p], 1.0f/3 ) / ( 4 * 3.1415926 / 3 );
			Eigen::Vector3f x = m_d.particleF[p] * Eigen::Vector3f(1,0,0);
			Eigen::Vector3f y = m_d.particleF[p] * Eigen::Vector3f(0,1,0);
			Eigen::Vector3f z = m_d.particleF[p] * Eigen::Vector3f(0,0,1);
			
			glBegin( GL_QUADS );
			glColor3f( 1,0,1 );
			glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( x[0] + y[0] ), m_d.particleX[p][1] + 0.5f * r * ( x[1] + y[1] ), 0 );
			glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( -x[0] + y[0] ), m_d.particleX[p][1] + 0.5f * r * ( -x[1] + y[1] ), 0 );
			glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( -x[0] - y[0] ), m_d.particleX[p][1] + 0.5f * r * ( -x[1] - y[1] ), 0 );
			glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( x[0] - y[0] ), m_d.particleX[p][1] + 0.5f * r * ( x[1] - y[1] ), 0 );
			glEnd();
			
			glBegin( GL_LINE_LOOP );
			glColor3f( 1,1,1 );
			glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( x[0] + y[0] ), m_d.particleX[p][1] + 0.5f * r * ( x[1] + y[1] ), 0 );
			glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( -x[0] + y[0] ), m_d.particleX[p][1] + 0.5f * r * ( -x[1] + y[1] ), 0 );
			glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( -x[0] - y[0] ), m_d.particleX[p][1] + 0.5f * r * ( -x[1] - y[1] ), 0 );
			glVertex3f( m_d.particleX[p][0] + 0.5f * r * ( x[0] - y[0] ), m_d.particleX[p][1] + 0.5f * r * ( x[1] - y[1] ), 0 );
			glEnd();

		}

		glBegin( GL_LINES );
		glColor3f( 1, 0, 0 );
		for( int i=0; i < m_g.dimensions()[0]; ++i )
		{
			for( int j=0; j < m_g.dimensions()[1]; ++j )
			{
				int index = m_g.coordsToIndex( Vector3i( i, j, m_g.dimensions()[2]/2 ) );
				Vector3f v = velocities.segment<3>( 3*index ) / sqrt( m_g.masses()[ index ] );
				glVertex3f( m_g.minCoord()[0] + i * m_g.gridH(), m_g.minCoord()[1] + j * m_g.gridH(), 0 );
				glVertex3f( m_g.minCoord()[0] + i * m_g.gridH() + v[0]*m_timeStep * 5, m_g.minCoord()[1] + j * m_g.gridH() + v[1]*m_timeStep * 5, 0 );
			}
		}
		glEnd();
		
		glEnable( GL_LIGHTING );
		for( size_t i=0; i < m_collisionObjects.size(); ++i )
		{
			m_collisionObjects[i]->draw();
		}
		glDisable( GL_LIGHTING );
	
		glutSwapBuffers();
	}

};
*/

////////////////////////////////////////////////////////////////////////////////
//! Display callback
////////////////////////////////////////////////////////////////////////////////
void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glOrtho( g_x - g_r, g_x + g_r, g_y - g_r, g_y + g_r, -1, 1);
	
	//static int count=0;
	//((CollisionPlane*)g_collisionObjects.objects.back())->setCoeffs( Eigen::Vector4f(0,-1,0,-1.5f + g_timeStep * count++  ) );
	
	MpmSim::SquareMagnitudeTermination t( 60, 0.01f );
	g_sim->advance( g_timeStep, t );
	glDisable( GL_DEPTH_TEST );

	const std::vector<Eigen::Vector3f>& particleX = g_sim->particleVariable<MpmSim::VectorData>( "p" )->m_data;
	const std::vector<Eigen::Matrix3f>& F = g_sim->particleVariable<MpmSim::MatrixData>( "F" )->m_data;
	const std::vector<float>& volume = g_sim->particleVariable<MpmSim::ScalarData>( "volume" )->m_data;
	
	for( size_t p = 0; p < particleX.size(); ++p )
	{
		float r = 2 * pow( volume[p], 1.0f/3 ) / ( 4 * 3.1415926 / 3 );
		Eigen::Vector3f x = F[p] * Eigen::Vector3f(1,0,0);
		Eigen::Vector3f y = F[p] * Eigen::Vector3f(0,1,0);
		Eigen::Vector3f z = F[p] * Eigen::Vector3f(0,0,1);
		float thingy = 50 * ( F[p].determinant() - 1 );
		
		glBegin( GL_QUADS );
		glColor3f( 1 - thingy, 0, 1 + thingy );
		glVertex3f( particleX[p][0] + 0.5f * r * ( x[0] + y[0] ), particleX[p][1] + 0.5f * r * ( x[1] + y[1] ), 0 );
		glVertex3f( particleX[p][0] + 0.5f * r * ( -x[0] + y[0] ), particleX[p][1] + 0.5f * r * ( -x[1] + y[1] ), 0 );
		glVertex3f( particleX[p][0] + 0.5f * r * ( -x[0] - y[0] ), particleX[p][1] + 0.5f * r * ( -x[1] - y[1] ), 0 );
		glVertex3f( particleX[p][0] + 0.5f * r * ( x[0] - y[0] ), particleX[p][1] + 0.5f * r * ( x[1] - y[1] ), 0 );
		glEnd();
		
		glBegin( GL_LINE_LOOP );
		glColor3f( 1,1,1 );
		glVertex3f( particleX[p][0] + 0.5f * r * ( x[0] + y[0] ), particleX[p][1] + 0.5f * r * ( x[1] + y[1] ), 0 );
		glVertex3f( particleX[p][0] + 0.5f * r * ( -x[0] + y[0] ), particleX[p][1] + 0.5f * r * ( -x[1] + y[1] ), 0 );
		glVertex3f( particleX[p][0] + 0.5f * r * ( -x[0] - y[0] ), particleX[p][1] + 0.5f * r * ( -x[1] - y[1] ), 0 );
		glVertex3f( particleX[p][0] + 0.5f * r * ( x[0] - y[0] ), particleX[p][1] + 0.5f * r * ( x[1] - y[1] ), 0 );
		glEnd();

	}
	
	glutSwapBuffers();

	glutPostRedisplay();
}


////////////////////////////////////////////////////////////////////////////////
//! Keyboard events handler
////////////////////////////////////////////////////////////////////////////////
void keyboard(unsigned char key, int /*x*/, int /*y*/)
{
	std::cerr << key << std::endl;
    switch(key)
	{
    case(27) :
        exit(0);
        break;
	case 'p':
		break;
    }
}

////////////////////////////////////////////////////////////////////////////////
//! Mouse event handlers
////////////////////////////////////////////////////////////////////////////////
void mouse(int button, int state, int x, int y)
{
    if (state == GLUT_DOWN)
	{
        mouse_buttons |= 1<<button;
    }
	else if (state == GLUT_UP)
	{
        mouse_buttons &= ~( 1<<button );
    }

    mouse_old_x = x;
    mouse_old_y = y;
    glutPostRedisplay();
}

void motion(int x, int y)
{
	float dx, dy;
	dx = x - mouse_old_x;
	dy = y - mouse_old_y;

	if (mouse_buttons & 1)
	{
		g_x -= 0.001 * dx * g_r;
		g_y += 0.001 * dy * g_r;
	}
	else if (mouse_buttons & 4)
	{
		g_r *= exp( 0.01 * dy );
	}

	mouse_old_x = x;
	mouse_old_y = y;
}
