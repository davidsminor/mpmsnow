#ifdef WIN32
#include <windows.h>
#endif

#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>

#include <GL/gl.h>
#include <GL/glut.h>

////////////////////////////////////////////////////////////////////////////////
// constants
const unsigned int window_width = 512;
const unsigned int window_height = 512;

int g_frame(0);
float g_vscale(1.0f);

float g_theta(0);
float g_phi(0);
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

int g_gridSize[3];
float g_gridH;
float g_gridMin[3];
std::vector<float> g_masses;
std::vector<float> g_explicitVelocities;
std::vector< std::vector<float> > g_iterates;

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
	//if( argc == 2 )
	{

		//std::ifstream f( argv[1], std::ofstream::binary );
		std::ifstream f( "debug.dat", std::ofstream::binary );
		
		// read grid size:
		f.read( (char*)&g_gridH, sizeof(float) );

		// read grid origin
		f.read( (char*)g_gridMin, 3 * sizeof(float) );
		
		// read grid dimensions:
		f.read( (char*)g_gridSize, 3 * sizeof(int) );

		// read masses:
		g_masses.resize( g_gridSize[0] * g_gridSize[1] * g_gridSize[2] );
		f.read( (char*)&g_masses.front(), g_masses.size() * sizeof(float) );
		
		// read explicit velocities:
		g_explicitVelocities.resize( 3 * g_gridSize[0] * g_gridSize[1] * g_gridSize[2] );
		f.read( (char*)&g_explicitVelocities.front(), g_explicitVelocities.size() * sizeof(float) );
		
		while( f )
		{
			g_iterates.resize( g_iterates.size() + 1 );
			g_iterates.back().resize( 3 * g_gridSize[0] * g_gridSize[1] * g_gridSize[2] );
			f.read( (char*)&g_iterates.back().front(), 3 * g_masses.size() * sizeof(float) );
		}
		std::cerr << "grid size: " << g_gridSize[0] << " "  << g_gridSize[1] << " "  << g_gridSize[2] << std::endl;
		std::cerr << g_iterates.size() << " iterates" << std::endl;
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
	glutCreateWindow("iterate viewer");
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


////////////////////////////////////////////////////////////////////////////////
//! Display callback
////////////////////////////////////////////////////////////////////////////////
void display()
{
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective( 45,
 		1.0f,
 		0.02,
 		100);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(
		g_r * cos( g_theta ) * cos( g_phi ),
 		g_r * sin( g_theta ),
 		g_r * cos( g_theta ) * sin( g_phi ),
 		0,0,0,
 		0,1,0);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDisable( GL_DEPTH_TEST );
	
	// draw axes:
	glBegin( GL_LINES );

	glColor3f( 1,0,0 );
	glVertex3f( 0,0,0 );
	glVertex3f( 1,0,0 );
	glColor3f( 0,1,0 );
	glVertex3f( 0,0,0 );
	glVertex3f( 0,1,0 );
	glColor3f( 0,0,1 );
	glVertex3f( 0,0,0 );
	glVertex3f( 0,0,1 );

	glEnd();
	
	int maxDim( g_gridSize[0] );
	if( g_gridSize[1] > maxDim ) maxDim = g_gridSize[1];
	if( g_gridSize[2] > maxDim ) maxDim = g_gridSize[2];
	
	// draw masses:
	glEnable( GL_BLEND );
	glBlendFunc( GL_ONE, GL_ONE );
	glDisable(GL_DEPTH_TEST);
	glDepthMask( GL_FALSE );
	
	glPointSize(4);
	glBegin( GL_POINTS );
	float x(g_gridMin[0]);
	for( int i=0; i < g_gridSize[0]; ++i, x += g_gridH )
	{
		float y(g_gridMin[1]);
		for( int j=0; j < g_gridSize[1]; ++j, y += g_gridH )
		{
			float z(g_gridMin[2]);
			for( int k=0; k < g_gridSize[2]; ++k, z += g_gridH )
			{
				int offset = i + g_gridSize[1] * ( j + g_gridSize[2] * k );
				glColor3f( g_masses[offset] == 0 ? 0.1f : 0.0f,0,g_masses[offset] );
				glVertex3f( x,y,z );

			}
		}

	}

	glEnd();

	glDisable( GL_BLEND );
	glEnable(GL_DEPTH_TEST);
	glDepthMask( GL_TRUE );


	// draw explicit velocities:
	glBegin( GL_LINES );
	
	glColor3f(0,1,0);

	x = (g_gridMin[0]);
	for( int i=0; i < g_gridSize[0]; ++i, x += g_gridH )
	{
		float y(g_gridMin[1]);
		for( int j=0; j < g_gridSize[1]; ++j, y += g_gridH )
		{
			float z(g_gridMin[2]);
			for( int k=0; k < g_gridSize[2]; ++k, z += g_gridH )
			{
				int offset = 3 * (i + g_gridSize[1] * ( j + g_gridSize[2] * k ));
				glVertex3f( x,y,z );
				glVertex3f(
					x + g_vscale * g_explicitVelocities[offset],
					y + g_vscale * g_explicitVelocities[offset+1],
					z + g_vscale * g_explicitVelocities[offset+2]
				);

			}
		}

	}

	// draw velocities:
	glBegin( GL_LINES );
	
	glColor3f(1,1,1);

	x = (g_gridMin[0]);
	for( int i=0; i < g_gridSize[0]; ++i, x += g_gridH )
	{
		float y(g_gridMin[1]);
		for( int j=0; j < g_gridSize[1]; ++j, y += g_gridH )
		{
			float z(g_gridMin[2]);
			for( int k=0; k < g_gridSize[2]; ++k, z += g_gridH )
			{
				int offset = 3 * (i + g_gridSize[1] * ( j + g_gridSize[2] * k ));
				glVertex3f( x,y,z );
				glVertex3f(
					x + g_vscale * g_iterates[g_frame][offset],
					y + g_vscale * g_iterates[g_frame][offset+1],
					z + g_vscale * g_iterates[g_frame][offset+2]
				);

			}
		}

	}

	glEnd();

	glutSwapBuffers();
}


////////////////////////////////////////////////////////////////////////////////
//! Keyboard events handler
////////////////////////////////////////////////////////////////////////////////
void keyboard(unsigned char key, int /*x*/, int /*y*/)
{
    switch(key)
	{
    case(27) :
        exit(0);
        break;
	case 's':
		++g_frame;
		if( g_frame == g_iterates.size() )
		{
			g_frame = g_iterates.size() - 1;
		}
		glutPostRedisplay();
		break;
	case 'a':
		--g_frame;
		if( g_frame < 0 )
		{
			g_frame = 0;
		}
		glutPostRedisplay();
		break;
	case 'q':
		g_vscale *= 0.9;
		glutPostRedisplay();
		break;
	case 'w':
		g_vscale /= 0.9;
		glutPostRedisplay();
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
}

void motion(int x, int y)
{
	float dx, dy;
	dx = x - mouse_old_x;
	dy = y - mouse_old_y;

	if (mouse_buttons & 1)
	{
		g_phi += 0.01 * dx;
		g_theta += 0.01 * dy;
	}
	else if (mouse_buttons & 4)
	{
		g_r *= exp( 0.01 * dy );
	}

	mouse_old_x = x;
	mouse_old_y = y;
    glutPostRedisplay();
}
