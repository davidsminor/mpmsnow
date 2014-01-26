#ifdef WIN32
#include <windows.h>
#endif

#ifdef HAVE_CORTEX
#include "IECore/ImagePrimitive.h"
#include "IECore/Writer.h"
#include <boost/format.hpp>
#endif



#include <algorithm>
#include <vector>
#include <iostream>

#include <GL/gl.h>
#include <GL/glut.h>

#include "Grid.h"

#include "CollisionPlane.h"
#include "ConjugateGradients.h"
#include "ConjugateResiduals.h"

using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// constants
const unsigned int window_width = 512;
const unsigned int window_height = 512;

float g_theta(0.2);
float g_phi(-0.2);
float g_r(6);

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
GLuint g_matrixTexture;

//#define DISPLAYMATRIX 1

// collision objects:
// ground plane at y = 0
CollisionPlane g_groundPlane( Vector4f( 0, 1, 0, 0 ) );
std::vector< CollisionObject* > g_collisionObjects;

#define LAT_DIVISIONS 4
#define LONG_DIVISIONS 6
std::vector< Vector3f > g_spherePoints;

int g_time(0);







////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv)
{
	
	// initial configuration:
	Vector3f rotVector( 1, 9, 3 );
	rotVector.normalize();
	float particleSpacing( 0.05 );
	float particleVolume = particleSpacing * particleSpacing * particleSpacing;

	srand( 10 );
	const int particlesPerCell = 3;
	for( int i=-3; i <= 3; ++i )
	{
		for( int j=-3; j <= 3; ++j )
		{
			for( int k=-3; k <= 3; ++k )
			{
				for( int n=0; n < particlesPerCell; ++n )
				{
					float xr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
					float yr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
					float zr = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
					Vector3f pos( particleSpacing * float( i + xr ), 0.5f + particleSpacing * float( j + yr ), particleSpacing * float( k + zr ) );
					
					g_particles.particleX.push_back( pos );
					g_particles.particleV.push_back( 0.6 * rotVector.cross( pos ) );
					g_particles.particleM.push_back( INITIALDENSITY * particleVolume / particlesPerCell );

					g_particles.particleF.push_back( Matrix3f::Identity() );
					g_particles.particleFplastic.push_back( Matrix3f::Identity() );
					g_particles.particleR.push_back( Matrix3f::Identity() );
					g_particles.particleS.push_back( Matrix3f::Identity() );
					g_particles.particleFinvTrans.push_back( Matrix3f::Identity() );
					g_particles.particleJ.push_back( 1.0f );
					
					g_particles.particleMu.push_back( MU );
					g_particles.particleLambda.push_back( LAMBDA );
				}
			}
		}
	}
	
	// set up collision objects:
	g_collisionObjects.push_back( &g_groundPlane );
	
	// make a little sphere to draw particles with:
	for( int i=0; i<LONG_DIVISIONS; ++i )
	{
		float latAngle = 3.1415926f / LAT_DIVISIONS;
		float longAngle0 = i * 2 * 3.1415926f / LONG_DIVISIONS;
		float longAngle1 = ( i + 1 ) * 2 * 3.1415926f / LONG_DIVISIONS;
		g_spherePoints.push_back( Vector3f( sin( latAngle ) * cos( longAngle0 ), cos( latAngle ), sin( latAngle ) * sin( longAngle0 ) ) );
		g_spherePoints.push_back( Vector3f( 0,1,0 ) );
		g_spherePoints.push_back( Vector3f( sin( latAngle ) * cos( longAngle1 ), cos( latAngle ), sin( latAngle ) * sin( longAngle1 ) ) );
		
		
		latAngle = 3.1415926f - latAngle;
		g_spherePoints.push_back( Vector3f( 0,-1,0 ) );
		g_spherePoints.push_back( Vector3f( sin( latAngle ) * cos( longAngle0 ), cos( latAngle ), sin( latAngle ) * sin( longAngle0 ) ) );
		g_spherePoints.push_back( Vector3f( sin( latAngle ) * cos( longAngle1 ), cos( latAngle ), sin( latAngle ) * sin( longAngle1 ) ) );
		
		for( int j=1; j < LAT_DIVISIONS-1; ++j )
		{
			float latAngle0 = j * 3.1415926f / LAT_DIVISIONS;
			float latAngle1 = ( j + 1 ) * 3.1415926f / LAT_DIVISIONS;

			g_spherePoints.push_back( Vector3f( sin( latAngle0 ) * cos( longAngle0 ), cos( latAngle0 ), sin( latAngle0 ) * sin( longAngle0 ) ) );
			g_spherePoints.push_back( Vector3f( sin( latAngle0 ) * cos( longAngle1 ), cos( latAngle0 ), sin( latAngle0 ) * sin( longAngle1 ) ) );
			g_spherePoints.push_back( Vector3f( sin( latAngle1 ) * cos( longAngle1 ), cos( latAngle1 ), sin( latAngle1 ) * sin( longAngle1 ) ) );
		
			g_spherePoints.push_back( Vector3f( sin( latAngle1 ) * cos( longAngle1 ), cos( latAngle1 ), sin( latAngle1 ) * sin( longAngle1 ) ) );
			g_spherePoints.push_back( Vector3f( sin( latAngle1 ) * cos( longAngle0 ), cos( latAngle1 ), sin( latAngle1 ) * sin( longAngle0 ) ) );
			g_spherePoints.push_back( Vector3f( sin( latAngle0 ) * cos( longAngle0 ), cos( latAngle0 ), sin( latAngle0 ) * sin( longAngle0 ) ) );
		}
	}
	
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

#ifdef DISPLAYMATRIX
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	if( g_matrixTexture == 0 )
	{
		// instantiate grid and rasterize!
		Grid g( g_particles );
		if( g_particles.particleVolumes.empty() )
		{
			g.computeDensities( g_particles );
			g_particles.particleVolumes.resize( g_particles.particleX.size() );
			for( size_t i = 0; i < g_particles.particleDensities.size(); ++i )
			{
				g_particles.particleVolumes[i] = g_particles.particleM[i] / g_particles.particleDensities[i];
			}
		}
		g_matrixTexture = g.matrixTexture( g_particles, g_collisionObjects );
	}
	
	glEnable( GL_TEXTURE_2D );
	glBindTexture( GL_TEXTURE_2D, g_matrixTexture );
	
	glBegin( GL_QUADS );
	glTexCoord2f(1,1);
	glVertex2f(1,-1);
	glTexCoord2f(1,0);
	glVertex2f(1,1);
	glTexCoord2f(0,0);
	glVertex2f(-1,1);
	glTexCoord2f(0,1);
	glVertex2f(-1,-1);
	glEnd();
	
	glDisable( GL_TEXTURE_2D );

    glutSwapBuffers();
    glutPostRedisplay();
	return;

#endif

	gluLookAt(	g_r * cos( g_theta ) * sin( g_phi ),
				g_r * sin( g_theta ),
 				g_r * cos( g_theta ) * cos( g_phi ),
 				0.0, 0.0, 0.0,
 				0.0, 1.0, 0.0
	);

	// display!
	/*
	glPointSize(2);
	glBegin( GL_POINTS );
	glColor3f( 1,1,1 );
	for( size_t i = 0; i < g_particles.particleX.size(); ++i )
	{
		glVertex3f( g_particles.particleX[i][0], g_particles.particleX[i][1], g_particles.particleX[i][2] );
	}

	glEnd();
	*/
	

	// lights:
	
	GLfloat light_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
	
	glEnable( GL_LIGHT0 );
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	
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
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnable( GL_LIGHTING );
	
    glVertexPointer(3, GL_FLOAT, 0, (float*)&g_spherePoints[0]);
    glNormalPointer(GL_FLOAT, 0, (float*)&g_spherePoints[0]);
	
	for( size_t p = 0; p < g_particles.particleX.size(); ++p )
	{
		float r = 2 * pow( g_particles.particleVolumes[p], 1.0f/3 ) / ( 4 * 3.1415926 / 3 );
		glPushMatrix();
		glTranslatef( g_particles.particleX[p][0],g_particles.particleX[p][1],g_particles.particleX[p][2] );
		glScalef( r, r, r );
		glDrawArrays(GL_TRIANGLES, 0, (int)g_spherePoints.size());
		glPopMatrix();
	}
	
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisable( GL_LIGHTING );
	
	// update grid velocities using internal stresses...
	g.updateGridVelocities( g_particles, g_collisionObjects, ConjugateGradients( 30, 1.e-10 ) );
	
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
	
#ifdef HAVE_CORTEX
	Imath::Box2i screenWindow( Imath::V2i( 0, 0 ), Imath::V2i( window_width-1, window_height-1 ) );
	IECore::ImagePrimitivePtr image = IECore::ImagePrimitive::createGreyscale( 0.5f, screenWindow, screenWindow );
	glReadPixels(1, 1, window_width, window_height, GL_RED, GL_FLOAT, &image->variableData< IECore::FloatVectorData >( "Y" )->writable()[0]);
	IECore::Writer::create( image, boost::str( boost::format( "./frame.%04d.tif" ) % g_time ) )->write();
#endif
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
