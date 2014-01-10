#ifdef WIN32
#include <windows.h>
#endif

#include <algorithm>
#include <vector>
#include <iostream>

#include <GL/gl.h>
#include <GL/glut.h>

#include "Grid.h"

using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// constants
const unsigned int window_width = 512;
const unsigned int window_height = 512;

float g_theta(0);
float g_phi(0);
float g_r(5);

// GL functionality
bool initGL(int *argc, char** argv);

// rendering callbacks
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);

int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;


ParticleData g_particles;

int g_time(0);







////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	
	// test code for force differential:
	Vector3f w(1,1,1);
	w.normalize();
	Matrix3f Rorig = AngleAxisf( 0.25*M_PI, w ).toRotationMatrix();
	Matrix3f Sorig = Matrix3f::Zero();
	Sorig(0,0) = 1;
	Sorig(1,1) = 2;
	Sorig(2,2) = 3;

	Matrix3f F = Matrix3f::Random();//Rorig * Sorig;
	Matrix3f R;
	Matrix3f S;
	Affine3f trans( F );	
	trans.computeRotationScaling( &R, &S );
	Matrix3f Finv;
	Matrix3f FinvTrans;
	float J;
	bool invertible;
	
	F.computeInverseAndDetWithCheck( Finv, J, invertible );
	FinvTrans = Finv.transpose();
	
	Matrix3f dEdF = 2 * MU * ( F - R ) + LAMBDA * ( J - 1 ) * J * FinvTrans;
	
	Matrix3f dF = Matrix3f::Random() * 0.0001;
	Matrix3f F_= F + dF;
	Matrix3f R_;
	Matrix3f S_;
	Affine3f trans2( F_ );	
	trans2.computeRotationScaling( &R_, &S_ );
	Matrix3f Finv_;
	Matrix3f FinvTrans_;
	float J_;
	
	F_.computeInverseAndDetWithCheck( Finv_, J_, invertible );
	FinvTrans_ = Finv_.transpose();
	
	Matrix3f dEdF_ = 2 * MU * ( F_ - R_ ) + LAMBDA * ( J_ - 1 ) * J_ * FinvTrans_;
	
	
	double dJ = J * Grid::matrixDoubleDot( FinvTrans, dF );
	
	
	Matrix3f M = R.transpose() * dF - dF.transpose() * R;
	
	Matrix3f G;
	G(0,0) = S(0,0) + S(1,1);
	G(1,1) = S(0,0) + S(2,2);
	G(2,2) = S(1,1) + S(2,2);
	
	G(0,1) = G(1,0) = S(1,2);
	G(0,2) = G(2,0) = -S(0,2);
	G(1,2) = G(2,1) = S(0,1);
	
	Vector3f m( M(0,1), M(0,2), M(1,2) );
	
	w = G.inverse() * m;
	
	Matrix3f RtdR;
	RtdR(0,0) = RtdR(1,1) = RtdR(2,2) = 0;
	
	RtdR(0,1) = w[0];
	RtdR(1,0) = -w[0];
	
	RtdR(0,2) = w[1];
	RtdR(2,0) = -w[1];
	
	RtdR(1,2) = w[2];
	RtdR(2,1) = -w[2];
	
	
	Matrix3f dR = R * RtdR;
	
	Matrix3f dFinvTrans = - FinvTrans * dF.transpose() * FinvTrans;
	
	Matrix3f Ap = 2 * MU * ( dF - dR ) + LAMBDA * ( dJ * J * FinvTrans + ( J - 1 ) * ( dJ * FinvTrans + J * dFinvTrans ) );
	
	std::cerr << dEdF_ - dEdF << std::endl << std::endl;
	std::cerr << Ap << std::endl << std::endl;
	
	
	
	
	
	// initial configuration:
	Vector3f rotVector( 1, 9, 3 );
	rotVector.normalize();
	float particleSpacing( 0.1 );
	float particleVolume = particleSpacing * particleSpacing * particleSpacing;

	for( int i=3; i <= 7; ++i )
	{
		for( int j=3; j <= 7; ++j )
		{
			for( int k=3; k <= 7; ++k )
			{
				Vector3f pos( 0.1 * float( i ) - 0.5, 0.1 * float( j ) - 0.5, 0.1 * float( k ) - 0.5 );
				
				g_particles.particleX.push_back( pos );
				g_particles.particleV.push_back( rotVector.cross( pos ) + 2 * rotVector.dot( pos ) * rotVector );
				g_particles.particleM.push_back( INITIALDENSITY * particleVolume );

				g_particles.particleF.push_back( Matrix3f::Identity() );
				g_particles.particleR.push_back( Matrix3f::Identity() );
				g_particles.particleS.push_back( Matrix3f::Identity() );
				g_particles.particleFinvTrans.push_back( Matrix3f::Identity() );
				g_particles.particleJ.push_back( 1.0f );

			}
		}
	}

	srand( 10 );
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
    glDisable(GL_DEPTH_TEST);
	glEnable( GL_CULL_FACE );

    // viewport
    glViewport(0, 0, window_width, window_height);

    // projection
    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();
	gluPerspective(	30,
 					double( window_width ) / window_height,
 					0.01,
 					20);
	
    return true;
}



////////////////////////////////////////////////////////////////////////////////
//! Display callback
////////////////////////////////////////////////////////////////////////////////
void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(	g_r * cos( g_theta ) * sin( g_phi ),
				g_r * sin( g_theta ),
 				g_r * cos( g_theta ) * cos( g_phi ),
 				0.0, 0.0, 0.0,
 				0.0, 1.0, 0.0
	);

	// display!
	glPointSize(2);
	glBegin( GL_POINTS );
	glColor3f( 1,1,1 );
	for( size_t i = 0; i < g_particles.particleX.size(); ++i )
	{
		glVertex3f( g_particles.particleX[i][0], g_particles.particleX[i][1], g_particles.particleX[i][2] );
	}

	glEnd();
	

	// instantiate grid and rasterize!
	Grid g( g_particles );
	g.draw();
	if( g_particles.particleVolumes.empty() )
	{
		g.computeDensities( g_particles );
		g_particles.particleVolumes.resize( g_particles.particleX.size() );
		for( size_t i = 0; i < g_particles.particleDensities.size(); ++i )
		{
			g_particles.particleVolumes[i] = g_particles.particleM[i] / g_particles.particleDensities[i];
		}
	}
	
	if( g_time == 6 )
	{
		g.testForces( g_particles );
		exit(0);
	}
	
	
	// update grid velocities using internal stresses...
	g.updateGridVelocities( g_particles );
	
	// transfer the grid velocities back onto the particles:
	g.updateParticleVelocities( g_particles );

	// update particle deformation gradients:
	g.updateDeformationGradients( g_particles );

	// update positions...
	for( size_t i = 0; i < g_particles.particleX.size(); ++i )
	{
		g_particles.particleX[i] += g_particles.particleV[i] * TIME_STEP;
	}
	++g_time;
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
		g_theta += dy * 0.2 * M_PI / 180.0;
		g_phi -= dx * 0.2 * M_PI / 180.0;
	}
	else if (mouse_buttons & 4)
	{
		g_r += dy * 0.01;
	}

	mouse_old_x = x;
	mouse_old_y = y;
}
