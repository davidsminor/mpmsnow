#ifdef WIN32
#include <windows.h>
#endif

#include <algorithm>
#include <vector>
#include <iostream>

#include <GL/gl.h>
#include <GL/glut.h>

#include "Grid.h"

#include "CollisionPlane.h"

using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// constants
const unsigned int window_width = 512;
const unsigned int window_height = 512;

float g_theta(0.1);
float g_phi(0.1);
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

// collision objects:
// ground plane at y = -1.5
CollisionPlane g_groundPlane( Vector4f( 0, 1, 0, 0 ) );
std::vector< CollisionObject* > g_collisionObjects;

int g_time(0);







////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv)
{
	
	// initial configuration:
	Vector3f rotVector( 1, 9, 3 );
	rotVector.normalize();
	float particleSpacing( 0.1 );
	float particleVolume = particleSpacing * particleSpacing * particleSpacing;

	srand( 10 );
	for( int i=3; i <= 7; ++i )
	{
		for( int j=3; j <= 7; ++j )
		{
			for( int k=3; k <= 7; ++k )
			{
				for( int n=0; n < 3; ++n )
				{
					float xr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
					float yr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
					float zr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
					Vector3f pos( particleSpacing * float( i + xr - 0.5f ) - 0.5, 1.0f + particleSpacing * float( j + yr - 0.5f ) - 0.5, particleSpacing * float( k + zr - 0.5f ) - 0.5 );
					
					g_particles.particleX.push_back( pos );
					g_particles.particleV.push_back( 0.2 * rotVector.cross( pos ) );
					g_particles.particleM.push_back( INITIALDENSITY * particleVolume );

					g_particles.particleF.push_back( Matrix3f::Identity() );
					g_particles.particleR.push_back( Matrix3f::Identity() );
					g_particles.particleS.push_back( Matrix3f::Identity() );
					g_particles.particleFinvTrans.push_back( Matrix3f::Identity() );
					g_particles.particleJ.push_back( 1.0f );
				}
			}
		}
	}
	
	// set up collision objects:
	g_collisionObjects.push_back( &g_groundPlane );

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
	
	// draw collision objects:
	for( size_t i=0; i < g_collisionObjects.size(); ++i )
	{
		g_collisionObjects[i]->draw();
	}

	// instantiate grid and rasterize!
	Grid g( g_particles );
	//g.draw();
	if( g_particles.particleVolumes.empty() )
	{
		g.computeDensities( g_particles );
		g_particles.particleVolumes.resize( g_particles.particleX.size() );
		for( size_t i = 0; i < g_particles.particleDensities.size(); ++i )
		{
			g_particles.particleVolumes[i] = g_particles.particleM[i] / g_particles.particleDensities[i];
		}
	}
	
	/*
	if( g_time == 6 )
	{
		g.testForceDifferentials( g_particles );
		exit(0);
	}
	*/
	
	
	// update grid velocities using internal stresses...
	g.updateGridVelocities( g_particles, g_collisionObjects );
	
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
