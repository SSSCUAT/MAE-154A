
  

#include <iostream>
using std::cerr; 
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <cstdlib> // for exit function
 
//#include <stdlib.h>

//#include <stdio.h>       /* Std. I/O routines */

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
 
#include "ac3d.h"

#include "imageloader.h"

//#include <gl/glut.h>     /* Graphics library  */
#include <time.h>        /* Timing routines   */
#include <math.h>        /* Math library      */ 

#include "flight.h"      /* Misc. definitions for flight routines                          */
#include "graphics.h"    /* Color & drawing definitions                                    */
#include "vehicle.h"     /* Vehicle Parameters
/* declare functions used by the simulation */
void initial_states(void);
void setup_graphics(void);
clock_t gettime(void);
void display_and_dynamics(void);
void reshape(int w, int h);
void mousebutton(int button, int state, int x, int y);
void keyboard(unsigned char key, int x, int y);
void arrowkeys(int key, int x, int y);
void mousepos(int x, int y);
void draw_view(void);
void read_state(void);
double limit(double val, double low, double hi);
void perspective_projection(void);

// declare functions used in graphics.c
float getgnd(int i, int j);
void setup_terrain(void);
void getcolr(int i, int j);
void draw_rect(int i, int j);
void draw_runway(void);
void draw_tower(void);
void draw_waypoint(void);
//void draw_plane(void);
void draw_new_plane(void);
void draw_own_vehicle(void);
void draw_ground(void);
void draw_airport(void);
void draw_ground_plane(void);
void draw_cockpit(void);
void draw_text_readout(void);
void print_string1(char *s, void *font);
void print_string2(char *s, void *font);
void print_string3(char *s, void *font);
void print_string4(char *s, void *font);
void draw_wings(int right_or_left);
void draw_canards(int right_or_left);
//void draw_fuselage(void);
void draw_propeller(void);
void draw_elevator(int right_or_left);

double get_time_step(void);
//clock_t gettime(void);
double difftime(clock_t nowtime, clock_t lasttime);

ifstream indata; // indata is like cin
double airspeed_read;
double alpha_read;
double beta_read;
double TARGET_N;
double TARGET_E;
double TARGET_D;
double cg_loc_x;
double cg_loc_z;

//added for ac3d rendering
//-------------------------

int display_list;
 

void init_gfx(void)
{
    ac_prepare_render();
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
	
//    set_projection(WINDOW_WIDTH,WINDOW_HEIGHT);
}
//----------------------

//double limit(double val, double low, double hi);
GLuint loadMipmappedTexture(Image *image) {
	GLuint textureId_dam;
	glGenTextures(1, &textureId_dam);
	glBindTexture(GL_TEXTURE_2D, textureId_dam);
	gluBuild2DMipmaps(GL_TEXTURE_2D,
					  GL_RGB,
					  image->width, image->height,
					  GL_RGB,
					  GL_UNSIGNED_BYTE,
					  image->pixels);
	return textureId_dam;
}


GLuint _textureId_dam; //The id of the texture
GLuint _textureId_dam2; //The id of the texture
GLuint _textureId_dam3; //The id of the texture
GLuint _textureId_dam4; //The id of the texture
GLuint _textureId_dam5; //The id of the texture




/*---------------------*/
/*----- MAIN LOOP -----*/
/*---------------------*/

int main(int argc, char** argv) {
	
	ACObject *ob;
    static char *acFileName="pioneer_body_tex.ac";


	
	/* set up graphics window */
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(200, 200);
	glutInitWindowPosition(10, 10);
	glutCreateWindow("Flight");
	
	
	/* initialize states and controls of the aircraft    */
	initial_states();

    glutFullScreen();

	setup_graphics();
	
	/* save the current time so we can compute the time step later */
	lasttime = gettime();
	
	/* define functions to be called during the sim loop                            */
	glutDisplayFunc(display_and_dynamics);   /* main simulation loop                */
	glutReshapeFunc(reshape);                /* called if the window is resized     */
	glutMouseFunc(mousebutton);              /* called if a mouse button is pressed */
	glutKeyboardFunc(keyboard);              /* called if a key is pressed          */
	glutSpecialFunc(arrowkeys);              /* called if an arrow key is pressed   */
	glutPassiveMotionFunc(mousepos);         /* called if the mouse moves           */
	
	//added for ac3d rendering
	//----------------------------------------
	//init_gfx();
	
	if ( argc > 1 ) acFileName = argv[argc-1];
	ob = ac_load_ac3d(acFileName);
	if (ob != NULL)
	{
		display_list = ac_display_list_render_object(ob);
		glutMainLoop();
	}
	else
	{
		printf("exiting.\n");
	}
    return 0;
	//-------------------------------------------

// commented out for ac3d rendering
//----------------------------------
//	glutMainLoop();                          /* indicates the loop can begin        */

//	return 0;
//-----------------------------------
	
}

/*-----------------------------*/
/* Initialize the state vector */
/*-----------------------------*/

void initial_states(void) {
	
	/* defined in <flight.h> */
	
	
	zoom = 1.0;                   /* variables for chase-plane view */
	pan = 0.0;
	sway = 0.0;
	view = OUT_THE_WINDOW;
	cg_loc_x = 13.0;
	cg_loc_z = 1.0;
}

/*-----------------------------------------*/
/* Function Calls for Display and Dynamics */
/*-----------------------------------------*/


void display_and_dynamics(void) {
	//ifstream indata; // indata is like cin
	double num; // variable for input value
	double num1; // variable for input value
	double num2; // variable for input value
	
	// render the graphics 
	draw_view();

    //going every 4th line to speed it up
	read_state();
    read_state();
    read_state();
    read_state();
	
		glutSwapBuffers();
	
	// refresh the display and continue the sim loop
	glutPostRedisplay();
}


/*-----------------------------*/
/* Reshape the graphics window */
/*-----------------------------*/

void reshape(int w, int h) {
	
	XMAXSCREEN = (double)(w);
	YMAXSCREEN = (double)(h);
	XCTR = XMAXSCREEN / 2.0;
	YCTR = YMAXSCREEN / 2.0;
	ASPECT_RATIO = XMAXSCREEN / YMAXSCREEN;
	
	glViewport(0, 0, (GLsizei)(w), (GLsizei)(h));
	
	perspective_projection();
}

/*---------------------------*/
/* Act on a mousebutton push */
/*---------------------------*/

void mousebutton(int button, int state, int x, int y) {
	switch(button) {
	case GLUT_LEFT_BUTTON:
		/* action to take goes here */	
		break;
	case GLUT_RIGHT_BUTTON:
		/* action to take goes here */	
		break;
	default:
		break;
	}
} 

/*--------------------------------*/
/* Act on a keyboard button press */
/*--------------------------------*/

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27:					/* Exit on ESCAPE */
		exit(0);
		break;
	case 't':					/* Toggle text display on screen */
		text_display = !text_display;
		break;
	case 'q':                          /* change the view */
		if (view == OUT_THE_WINDOW)
			view = CHASEPLANE;
		else if (view == CHASEPLANE)
			view = TOWER;
		else
			view = OUT_THE_WINDOW;
		break;
	case '[':                        /* zoom out */
		zoom /= 1.1;
		break;
	case ']':                        /* zoom in */
		zoom *= 1.1;
		break;
	default:
		break;
	}
}

/*-----------------------*/
/* Act on mouse movement */
/*-----------------------*/

void mousepos(int x, int y) {
	double mx, my;
	
	mx = (double)(x) - XCTR;
	my = (double)(y) - YCTR;
	
//	control[AILERON] = mx/MOUSERADIUS;
//	control[ELEVATOR_2] = -my/MOUSERADIUS;
	
//	control[AILERON] = limit(control[AILERON], -1.0, 1.0);
//	control[ELEVATOR_2] = limit(control[ELEVATOR_2], -1.0, 1.0);
}


/*------------------------*/
/* Act on arrow key press */
/*------------------------*/

void arrowkeys(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_LEFT:			/* Pan Right (yaw inc.) left arrow */
		pan += 5.0;
		break;
	case GLUT_KEY_RIGHT:		/* Pan Left (yaw dec.) right arrow*/
		pan -= 5.0;
		break;
	case GLUT_KEY_UP:			/* Sway up (pitch dec.) up arrow */
		sway -= 5.0;
		break;
	case GLUT_KEY_DOWN:			/* Sway down (pitch inc.) down arrow */
		sway += 5.0;
		break;
	default:
		break;
	}
}

/*-----------------*/
/* get system time */
/*-----------------*/

clock_t gettime(void) {
	return clock();
}

/*----------------------------------------*/
/* limit val to remain between low and hi */
/*----------------------------------------*/

double limit(double val, double low, double hi) {
	if (val > hi)
		return hi;
	else if (val < low)
		return low;
	else
		return val;
}


/*-------------------------------------------*/
/* setup initial conditions, setup graphics  */
/*-------------------------------------------*/


void setup_graphics(void) {
	float fogcolor[] = {FOGR, FOGG, FOGB};
	
	// Set up terrain colors and altitudes 
	setup_terrain();
	
  	glClearColor(FOGR, FOGG, FOGB, 0.0);
	glShadeModel(GL_FLAT);
	
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
    
    
	Image* image = loadBMP("ocean_v3_256.bmp");
	//	_textureId_dam = loadTexture(image);
	_textureId_dam = loadMipmappedTexture(image);
	delete image;
    
	
    /*
	image = loadBMP("f4_cockpit3.bmp");
	_textureId_dam2 = loadMipmappedTexture(image);
	delete image;
	*/
    
	image = loadBMP("mountain_tex.bmp");
	_textureId_dam3 = loadMipmappedTexture(image);
	delete image;
    
	// set up fog depthcuing
	glEnable(GL_FOG); 
	{
		glFogi(GL_FOG_MODE, GL_EXP);
		glFogfv(GL_FOG_COLOR, fogcolor);
		glFogf(GL_FOG_DENSITY, FOG_DENSITY);
		glHint(GL_FOG_HINT, GL_DONT_CARE);
		glFogf(GL_FOG_START, FOG_START);
		glFogf(GL_FOG_END, FOG_END);
	}
	
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	// turn on display of textual information on screen 
	text_display = TRUE;
}

/*--------------------------------------------------*/
/* define a perspective projection transform matrix */
/*--------------------------------------------------*/

void perspective_projection(void) {
	double FOV2;
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	FOV2 = FOV*zoom;
	gluPerspective(FOV, ASPECT_RATIO, 10.0, 500000.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

//------------------------------
// draw the out-the-window view 
//------------------------------
 

void draw_view(void) {
	
	glPushMatrix();
	
	if (view == OUT_THE_WINDOW) {    // out-the-window view 
		
		// Transformations to set up coords to draw ground 
		glRotatef(90.0, 0.0, 1.0, 0.0);
		glRotatef(90.0, 1.0, 0.0, 0.0);
		
		// xyz are now oriented with body axis 
		
		// undo Euler angles 
		glRotatef(-state[PHI] * 180.0 / PI, 1.0, 0.0, 0.0);
		glRotatef(-state[THETA] * 180.0 / PI, 0.0, 1.0, 0.0);
		glRotatef(-state[PSI] * 180.0 /PI, 0.0, 0.0, 1.0);
		
		// xyz are now oriented with NED axes, but still centered on aircraft 
		
		// undo vehicle's position 
		glTranslatef(-state[NORTH], -state[EAST], -state[DOWN]);
	}
	 
	else if (view == CHASEPLANE) {     // chase plane view 
		
		glRotatef(90.0, 0.0, 1.0, 0.0);
		glRotatef(90.0, 1.0, 0.0, 0.0);
		
		glRotatef(-sway, 0.0, 1.0, 0.0);
		glRotatef(-pan, 0.0, 0.0, 1.0);
		
		 
//		glTranslatef(-state[NORTH]+100*zoom*cos(pan*PI/180)*cos(-sway*PI/180),-state[EAST]+100*zoom*sin(pan*PI/180)*cos(-sway*PI/180),-state[DOWN]-100*zoom*sin(sway*PI/180));
		glTranslatef(-state[NORTH]+100*zoom*cos(pan*PI/180)*cos(-sway*PI/180),-state[EAST]+100*zoom*sin(pan*PI/180)*cos(-sway*PI/180),1020-100*zoom*sin(sway*PI/180));
		
		
		// enter viewing transforms here to create proper chaseplane view 
	}
	
	else {     // (view == TOWER) --> tower view 
		
		glRotatef(90.0, 0.0, 1.0, 0.0);
		glRotatef(90.0, 1.0, 0.0, 0.0);
		glRotatef(atan2((state[DOWN]-TOWER_D),sqrt((state[NORTH]-TOWER_N)*(state[NORTH]-TOWER_N)+(state[EAST]-TOWER_E)*(state[EAST]-TOWER_E)))*180/PI,0.0,1.0,0.0);
		glRotatef(-atan2((state[EAST]-TOWER_E),(state[NORTH]-TOWER_N))*180/PI, 0.0, 0.0, 1.0);
		glTranslatef(-TOWER_N, -TOWER_E, -TOWER_D + 10);
		// enter viewing transforms here to create proper tower view 
	}
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GLfloat ambientLight[] = {0.8f, 0.8f, 0.8f, 1.0f};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);
	
	// draw a large ground plane to form horizon 
draw_ground_plane();
	
	// draw runways and control tower 
	draw_airport();
	
	// draw another aircraft in the environment 
//	draw_plane();
	
	// draw terrain field 
	draw_ground();
	
	// draw own aircraft if view from chaseplane or tower 
	if (view != OUT_THE_WINDOW){
		draw_canards(-1);
		draw_canards(1);
		draw_wings(1);		//	Right wing
		draw_wings(-1);		//  Left wing 
		//draw_fuselage();
		draw_propeller();
		draw_elevator(-1);
		draw_elevator(1);
		draw_new_plane();
		
	}
	 
	
	glPopMatrix();
	
	// draw the cockpit (currently only center indicator + text) 
	draw_cockpit();
}


/*-----------------------*/
/* draw the own aircraft */
/*-----------------------*/

void draw_propeller() {
	glPushMatrix();
	glTranslatef(state[NORTH], state[EAST], state[DOWN]);
	glRotatef(state[PSI]*180/PI, 0.0, 0.0, 1.0);
	glRotatef(state[THETA]*180/PI, 0.0, 1.0, 0.0);
	glRotatef(state[PHI]*180/PI, 1.0, 0.0, 0.0);
	
	glTranslatef(-7.5+cg_loc_x,0,-1.5+cg_loc_z);
	
	glRotatef(1.2*clock(), 1.0, 0.0, 0.0 );
	// two-bladed prop 
	glColor3f(.8, .8, .8);
	glBegin(GL_POLYGON);
	{
		glVertex3f(-10.3, 2.8, -0.25);
		glVertex3f(-10.3, -2.8, -0.25);
		glVertex3f(-10.3, -2.8 , 0.25);
		glVertex3f(-10.3, 2.8, 0.25);
	}
	glEnd();
	
	glPopMatrix();
	
}





void draw_wings(int right_or_left) {
	double root_chord;
	double tip_chord;
	double root_thickness;
	double tip_thickness;
	double aileron_start;
	double aileron_end;
	double aileron_back;
	double aileron_back2;
	double sweep_back;
	double dihedral;
	double half_span;
	
	double winglet_height;
	double winglet_tip_chord;
	double winglet_sweep_back;
	
	glPushMatrix();
	glTranslatef(state[NORTH], state[EAST], state[DOWN]);
	glRotatef(state[PSI]*180/PI, 0.0, 0.0, 1.0);
	glRotatef(state[THETA]*180/PI, 0.0, 1.0, 0.0);
	glRotatef(state[PHI]*180/PI, 1.0, 0.0, 0.0);

	glTranslatef(-12.0+cg_loc_x, 0.0, -2.0+cg_loc_z);
	glScalef(1.5, 1.5, 1.5);
	
	root_chord = -1.5;
	tip_chord = -1.5;
	root_thickness = .5;
	tip_thickness = .5;
	sweep_back = 0;
	dihedral = 0;
	half_span = 8.5*right_or_left;
	aileron_start = 6.5*right_or_left;
	aileron_end = 8.5*right_or_left;
	aileron_back = root_chord + ((tip_chord+sweep_back)-root_chord)*aileron_start/half_span;
	aileron_back2 = root_chord + ((tip_chord+sweep_back)-root_chord)*aileron_end/half_span;
	
	winglet_height = 6.25;
	winglet_tip_chord = -2;
	winglet_sweep_back = -1.4;
	
	// front panels 
	glColor3f(.3, .3, .3);
	glBegin(GL_POLYGON);
	{
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f((root_chord/4), 0.0, root_thickness/2);
		glVertex3f((sweep_back + tip_chord/4),half_span , tip_thickness/2 - dihedral);
		glVertex3f(sweep_back, half_span, -dihedral);
	}
	glEnd();
	glColor3f(.3, .3, .3);
	glBegin(GL_POLYGON);
	{
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f((root_chord/4), 0.0, -root_thickness/2);
		glVertex3f((sweep_back + tip_chord/4),half_span , -tip_thickness/2 - dihedral);
		glVertex3f(sweep_back, half_span, -dihedral);
	}
	glEnd();
	
	
	// Mid panels 
	glColor3f(.4, .4, .4);
	glBegin(GL_POLYGON);
	{
		glVertex3f((root_chord/4), 0.0, root_thickness/2);
		glVertex3f((sweep_back + tip_chord/4),half_span , tip_thickness/2 - dihedral);
		glVertex3f((sweep_back + 3*tip_chord/4),half_span , tip_thickness/6 - dihedral);
		glVertex3f(3*root_chord/4, 0.0, root_thickness/6);
	}
	glEnd();
	glColor3f(.4, .4, .4);
	glBegin(GL_POLYGON);
	{
		glVertex3f((root_chord/4), 0.0, -root_thickness/2);
		glVertex3f((sweep_back + tip_chord/4),half_span , -tip_thickness/2 - dihedral);
		glVertex3f((sweep_back + 3*tip_chord/4),half_span , -tip_thickness/6 - dihedral);
		glVertex3f(3*root_chord/4, 0.0, -root_thickness/6);
	}
	glEnd();
	
	
	// rear panels 
	glColor3f(.5, .5, .5);
	glBegin(GL_POLYGON);
	{
		glVertex3f(3*root_chord/4, 0.0, root_thickness/6);
		glVertex3f(root_chord, 0.0 , 0.0);
		glVertex3f(aileron_back,aileron_start,-dihedral*aileron_start/half_span);
	glVertex3f(0.75*root_chord - (.75*root_chord-(.75*tip_chord+sweep_back))*aileron_start/half_span,aileron_start,(-dihedral+root_thickness/6)*aileron_start/half_span);	}
	glEnd();
	glColor3f(.5, .5, .5);
	glBegin(GL_POLYGON);
	{
		glVertex3f(3*root_chord/4, 0.0, -root_thickness/6);
		glVertex3f(root_chord, 0.0 , 0.0);
		glVertex3f(aileron_back,aileron_start,-dihedral*aileron_start/half_span);
	glVertex3f(0.75*root_chord - (.75*root_chord-(.75*tip_chord+sweep_back))*aileron_start/half_span,aileron_start,(-dihedral-root_thickness/6)*aileron_start/half_span);	}
	glEnd();
	
	
	glColor3f(.5, .5, .5);
	glBegin(GL_POLYGON);
	{
		glVertex3f(0.75*root_chord - (.75*root_chord-(.75*tip_chord+sweep_back))*aileron_end/half_span,aileron_end,(-dihedral+root_thickness/6)*aileron_end/half_span);
		glVertex3f(aileron_back2,aileron_end,-dihedral*aileron_end/half_span);
		glVertex3f(tip_chord+sweep_back, half_span , -dihedral);
		glVertex3f(3*tip_chord/4+sweep_back, half_span, tip_thickness/6-dihedral);
	}
	glEnd();
	glColor3f(.5, .5, .5);
	glBegin(GL_POLYGON);
	{
		glVertex3f(0.75*root_chord - (.75*root_chord-(.75*tip_chord+sweep_back))*aileron_end/half_span,aileron_end,(-dihedral-root_thickness/6)*aileron_end/half_span);
		glVertex3f(aileron_back2,aileron_end,-dihedral*aileron_end/half_span);
		glVertex3f(tip_chord+sweep_back, half_span , -dihedral);
		glVertex3f(3*tip_chord/4+sweep_back, half_span, -tip_thickness/6-dihedral);
	}
	glEnd();
	
	// Aileron 
	glPushMatrix();
	glRotatef((control[AILERON]*180/PI*-right_or_left)/5, 0.0, 1.0, 0.0);
	glTranslatef(0.0,0.0,(.75*root_chord - (.75*root_chord-(.75*tip_chord+sweep_back))*aileron_start/half_span)*sin((control[AILERON]*-right_or_left)/5));
	
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_POLYGON);
	{
		glVertex3f(aileron_back,aileron_start,-dihedral*aileron_start/half_span);
		glVertex3f(.75*root_chord - (.75*root_chord-(.75*tip_chord+sweep_back))*aileron_start/half_span,aileron_start,(-dihedral+root_thickness/6)*aileron_start/half_span);
		glVertex3f(.75*root_chord - (.75*root_chord-(.75*tip_chord+sweep_back))*aileron_end/half_span,aileron_end,(-dihedral+root_thickness/6)*aileron_end/half_span);
		glVertex3f(aileron_back2,aileron_end,-dihedral*aileron_end/half_span);
		
	}
	glEnd();
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_POLYGON);
	{
		glVertex3f(aileron_back,aileron_start,-dihedral*aileron_start/half_span);
		glVertex3f(0.75*root_chord - (.75*root_chord-(.75*tip_chord+sweep_back))*aileron_start/half_span,aileron_start,(-dihedral-root_thickness/6)*aileron_start/half_span);
		glVertex3f(0.75*root_chord - (.75*root_chord-(.75*tip_chord+sweep_back))*aileron_end/half_span,aileron_end,(-dihedral-root_thickness/6)*aileron_end/half_span);
		glVertex3f(aileron_back2,aileron_end,-dihedral*aileron_end/half_span);
		
	}
	glEnd();
	glPopMatrix();
	
	glPopMatrix();
	
} 

void draw_canards(int right_or_left) {
	double canard_root_chord;
	double canard_tip_chord;
	double canard_root_thickness;
	double canard_tip_thickness;
	double elevator_start;
	double elevator_end;
	double elevator_back;
	double elevator_back2;
	double canard_sweep_back;
	double canard_dihedral;
	double canard_half_span;
	double LE_to_LE;
	
	glPushMatrix();
	glTranslatef(state[NORTH], state[EAST], state[DOWN]);
	glRotatef(state[PSI]*180/PI, 0.0, 0.0, 1.0);
	glRotatef(state[THETA]*180/PI, 0.0, 1.0, 0.0);
	glRotatef(state[PHI]*180/PI, 1.0, 0.0, 0.0);
	
	
	canard_root_chord = -1.2;
	canard_tip_chord = -1.2;
	canard_root_thickness = .5;
	canard_tip_thickness = 0.25;
	canard_sweep_back = 0.0;//-4.2;
	canard_dihedral = 0.0;//.5;
	canard_half_span = 6*right_or_left;
	elevator_start = 0*right_or_left;
	elevator_end = 5*right_or_left;
	elevator_back = canard_root_chord + ((canard_tip_chord+canard_sweep_back)-canard_root_chord)*elevator_start/canard_half_span;
	elevator_back2 = canard_root_chord + ((canard_tip_chord+canard_sweep_back)-canard_root_chord)*elevator_end/canard_half_span;
	
	LE_to_LE = -24;
	glTranslatef(LE_to_LE+cg_loc_x,0.0,-2.0+cg_loc_z);
	// front panels 
	glColor3f(.3, .3, .3);
	glBegin(GL_POLYGON);
	{
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f((canard_root_chord/4), 0.0, canard_root_thickness/2);
		glVertex3f((canard_sweep_back + canard_tip_chord/4),canard_half_span , canard_tip_thickness/2 - canard_dihedral);
		glVertex3f(canard_sweep_back, canard_half_span, -canard_dihedral);
	}
	glEnd();
	glColor3f(.3, .3, .3);
	glBegin(GL_POLYGON);
	{
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f((canard_root_chord/4), 0.0, -canard_root_thickness/2);
		glVertex3f((canard_sweep_back + canard_tip_chord/4),canard_half_span , -canard_tip_thickness/2 - canard_dihedral);
		glVertex3f(canard_sweep_back, canard_half_span, -canard_dihedral);
	}
	glEnd();
	
	
	// Mid panels 
	glColor3f(.5, .5, .5);
	glBegin(GL_POLYGON);
	{
		glVertex3f((canard_root_chord/4), 0.0, canard_root_thickness/2);
		glVertex3f((canard_sweep_back + canard_tip_chord/4),canard_half_span , canard_tip_thickness/2 - canard_dihedral);
		glVertex3f((canard_sweep_back + 3*canard_tip_chord/4),canard_half_span , canard_tip_thickness/6 - canard_dihedral);
		glVertex3f(3*canard_root_chord/4, 0.0, canard_root_thickness/6);
	}
	glEnd();
	glColor3f(.5, .5, .5);
	glBegin(GL_POLYGON);
	{
		glVertex3f((canard_root_chord/4), 0.0, -canard_root_thickness/2);
		glVertex3f((canard_sweep_back + canard_tip_chord/4),canard_half_span , -canard_tip_thickness/2 - canard_dihedral);
		glVertex3f((canard_sweep_back + 3*canard_tip_chord/4),canard_half_span , -canard_tip_thickness/6 - canard_dihedral);
		glVertex3f(3*canard_root_chord/4, 0.0, -canard_root_thickness/6);
	}
	glEnd();
	
	
	// rear panels 
	glColor3f(.7, .7, .7);
	glBegin(GL_POLYGON);
	{
		glVertex3f(3*canard_root_chord/4, 0.0, canard_root_thickness/6);
		glVertex3f(canard_root_chord, 0.0 , 0.0);
		glVertex3f(elevator_back,elevator_start,-canard_dihedral*elevator_start/canard_half_span);
		glVertex3f(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_start/canard_half_span,elevator_start,(-canard_dihedral+canard_root_thickness/6)*elevator_start/canard_half_span);
	}
	glEnd();
	glColor3f(.7, .7, .7);
	glBegin(GL_POLYGON);
	{
		glVertex3f(3*canard_root_chord/4, 0.0, -canard_root_thickness/6);
		glVertex3f(canard_root_chord, 0.0 , 0.0);
		glVertex3f(elevator_back,elevator_start,-canard_dihedral*elevator_start/canard_half_span);
		glVertex3f(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_start/canard_half_span,elevator_start,(-canard_dihedral-canard_root_thickness/6)*elevator_start/canard_half_span);
	}
	glEnd();
	
	
	glColor3f(.7, .7, .7);
	glBegin(GL_POLYGON);
	{
		glVertex3f(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_end/canard_half_span,elevator_end,(-canard_dihedral+canard_root_thickness/6)*elevator_end/canard_half_span);
		glVertex3f(elevator_back2,elevator_end,-canard_dihedral*elevator_end/canard_half_span);
		glVertex3f(canard_tip_chord+canard_sweep_back, canard_half_span , -canard_dihedral);
		glVertex3f(3*canard_tip_chord/4+canard_sweep_back, canard_half_span, canard_tip_thickness/6-canard_dihedral);
	}
	glEnd();
	glColor3f(.7, .7, .7);
	glBegin(GL_POLYGON);
	{
		glVertex3f(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_end/canard_half_span,elevator_end,(-canard_dihedral-canard_root_thickness/6)*elevator_end/canard_half_span);
		glVertex3f(elevator_back2,elevator_end,-canard_dihedral*elevator_end/canard_half_span);
		glVertex3f(canard_tip_chord+canard_sweep_back, canard_half_span , -canard_dihedral);
		glVertex3f(3*canard_tip_chord/4+canard_sweep_back, canard_half_span, -canard_tip_thickness/6-canard_dihedral);
	}
	glEnd();
	

	
	glPopMatrix();
	
} 

void draw_elevator(int right_or_left) {
	double canard_root_chord;
	double canard_tip_chord;
	double canard_root_thickness;
	double canard_tip_thickness;
	double elevator_start;
	double elevator_end;
	double elevator_back;
	double elevator_back2;
	double canard_sweep_back;
	double canard_dihedral;
	double canard_half_span;
	double LE_to_LE;
	
	glPushMatrix();


	glTranslatef(state[NORTH], state[EAST], state[DOWN]);
	glRotatef(state[PSI]*180/PI, 0.0, 0.0, 1.0);
	glRotatef(state[THETA]*180/PI, 0.0, 1.0, 0.0);
	glRotatef(state[PHI]*180/PI, 1.0, 0.0, 0.0);
	
	
	canard_root_chord = -1.2;
	canard_tip_chord = -1.2;
	canard_root_thickness = .5;
	canard_tip_thickness = 0.25;
	canard_sweep_back = 0.0;//-4.2;
	canard_dihedral = 0.0;//.5;
	canard_half_span = 6*right_or_left;
	elevator_start = 0*right_or_left;
	elevator_end = 5*right_or_left;
	elevator_back = canard_root_chord + ((canard_tip_chord+canard_sweep_back)-canard_root_chord)*elevator_start/canard_half_span;
	elevator_back2 = canard_root_chord + ((canard_tip_chord+canard_sweep_back)-canard_root_chord)*elevator_end/canard_half_span;
	
	LE_to_LE = -24;
	glTranslatef(LE_to_LE+cg_loc_x,0.0,-2.0+cg_loc_z);

	
	// Elevator 
	glPushMatrix();
	
	glTranslatef((0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_end/canard_half_span),0.0,0.0);
	glTranslatef(0.0,0.0,-(-canard_dihedral+canard_root_thickness/6)*elevator_end/canard_half_span);

	glRotatef(-control[ELEVATOR_2]*180/PI, 0.0, 1.0, 0.0);
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_POLYGON);
	{
		glVertex3f(0.0,elevator_end,0.0);
		glVertex3f(elevator_back2-(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_end/canard_half_span),elevator_end,-canard_dihedral*elevator_end/canard_half_span - (-canard_dihedral+canard_root_thickness/6)*elevator_end/canard_half_span);
		glVertex3f(elevator_back-(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_end/canard_half_span),elevator_start,-canard_dihedral*elevator_start/canard_half_span - (-canard_dihedral+canard_root_thickness/6)*elevator_end/canard_half_span);
		glVertex3f(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_start/canard_half_span-(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_end/canard_half_span),elevator_start,(-canard_dihedral+canard_root_thickness/6)*elevator_start/canard_half_span - (-canard_dihedral+canard_root_thickness/6)*elevator_end/canard_half_span);
	}
	glEnd();
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_POLYGON);
	{
		glVertex3f(0.0,elevator_end,(-canard_dihedral-canard_root_thickness/6)*elevator_end/canard_half_span);
		glVertex3f(elevator_back2-(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_end/canard_half_span),elevator_end,-canard_dihedral*elevator_end/canard_half_span);
		glVertex3f(elevator_back-(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_end/canard_half_span),elevator_start,-canard_dihedral*elevator_start/canard_half_span);
		glVertex3f(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_start/canard_half_span-(0.75*canard_root_chord - (.75*canard_root_chord-(.75*canard_tip_chord+canard_sweep_back))*elevator_end/canard_half_span),elevator_start,(-canard_dihedral-canard_root_thickness/6)*elevator_start/canard_half_span);
	}
	glEnd();
	
	glPopMatrix();
	
	
	glPopMatrix();
	
} 



 void draw_ground_plane(void) {
	
	// turn off the zbuffer while we draw the ground plane
	glDisable(GL_DEPTH_TEST);
	
	glColor3f(ccr[0][0], ccg[0][0], ccb[0][0]);
	glPushMatrix();
	glTranslatef(0.0, 0.0, 20.0);
	glRectf(-500000.0, -500000.0, 500000.0, 500000.0);
	glPopMatrix();
	
	glEnable(GL_DEPTH_TEST);
}

/*-------------------------------------------------------------------*/
/* draw two runways and a control tower near the center of the world */
/*-------------------------------------------------------------------*/


void draw_airport(void) {
	
	// Save the viewing transformation matrix 
	glPushMatrix();
	
	// draw 1st runway
	glTranslatef(TMAXX*FPB/2.0, TMAXY*FPB/2.0, 0.0);
	glRotatef(-90.0, 0.0, 0.0, 1.0);
	draw_runway();
	
	// restore and re-save transformation matrix 
	glPopMatrix();
	glPushMatrix();
	
	// draw 2nd runway 
	glTranslatef(TMAXX*FPB/2.0+400.0, TMAXY*FPB/2.0-3000.0, 2.0);
	glRotatef(-20.0, 0.0, 0.0, 1.0);
	draw_runway();
	
	// restore and re-save transformation matrix 
	glPopMatrix();
	glPushMatrix();
	
	glTranslatef(TMAXX*FPB/2.0-100.0, TMAXY*FPB/2.0-1000.0, 2.0);
	glRotatef(-90.0, 0.0, 0.0, 1.0);
	draw_runway();
	
	
	// restore and re-save transformation matrix 
	glPopMatrix();
	
	// draw the control tower 
	draw_tower();
	
	draw_waypoint();

}

/*--------------------------------------------*/
/* draw a runway rectangle, including stripes */
/*--------------------------------------------*/

void draw_runway(void) {
	int i, numstripes = 20;
	
	// draw runway rectangle 
	glPushMatrix();
	glColor3f(0.27, 0.27, 0.27);
	
	// turn on antialiasing? 
	glRectf(-150.0, 0.0, 150.0, 8000.0);
	
	glColor3f(1.0, 1.0, 1.0);
	
	// lift stripes of the runway so they are visible due to z-buffer 
	glTranslatef(0.0, 0.0, -5.0);
	
	// Draw Stripe 
	for (i = 0; i < numstripes; i++) {
		glRectf(-10.0, i * 8000.0 / numstripes,
				10.0, ( i + 0.5) * 8000.0 / numstripes);
		// antialias? 
	}
	
	glPopMatrix();
}


void draw_tower(void) {
	
	glPushMatrix();
	glColor3f(.8, .8, .8);
	glTranslatef(TOWER_N, TOWER_E, 0.0);
	
	glPushMatrix();
	glRotatef(90.0, 1.0, 0.0, 0.0);
	glTranslatef(0.0, 0.0, 20.0);
	glRectf(-20.0, 0.0, 20.0, TOWER_D);
	glTranslatef(0.0, 0.0, -40.0);
	glRectf(-20.0, 0.0, 20.0, TOWER_D);
	glPopMatrix();
	
	glPushMatrix();
	glRotatef(90.0, 0.0, 1.0, 0.0);
	glTranslatef(0.0, 0.0, 20.0);
	glRectf(0.0, -20.0, -TOWER_D, 20.0);
	glTranslatef(0.0, 0.0, -40.0);
	glRectf(0.0, -20.0, -TOWER_D, 20.0);
	glPopMatrix();
	glLineWidth(1);
	
	glPopMatrix();
}

void draw_waypoint(void) {
	
	glPushMatrix();
	glColor3f(1, 0, 0);
	glTranslatef(TARGET_N, TARGET_E, TARGET_D);
	
	glPushMatrix();
	glRotatef(90.0, 1.0, 0.0, 0.0);
	glTranslatef(0.0, 0.0, 20.0);
	glRectf(-20.0, 0.0, 20.0, 40);
	glTranslatef(0.0, 0.0, -40.0);
	glRectf(-20.0, 0.0, 20.0, 40);
	glPopMatrix();
	
	glPushMatrix();
	glRotatef(90.0, 0.0, 1.0, 0.0);
	glTranslatef(0.0, 0.0, 20.0);
	glRectf(0.0, -20.0, -40, 20.0);
	glTranslatef(0.0, 0.0, -40.0);
	glRectf(0.0, -20.0, -40, 20.0);
	glPopMatrix();
	glLineWidth(1);
	
	glPopMatrix();
}

/*---------------------*/
/* setup terrain model */
/*---------------------*/

void setup_terrain() {
	int i, j;
	
	//go thru all terrain blocks. Find altitude & color for each  
	// altitude stored in gnd[i][j], and color in cc[i][j]         
	// ground altitude in Down coords, so <0 means higher altitude 
	
	for(i=0; i<TMAXX; i++)
		for (j=0; j<TMAXY; j++) {
			if ((i<TMAXX/2-15 || i>TMAXX/2+15) || (j<TMAXY/2-15 || j>TMAXY/2+15))
				gnd[i][j] = getgnd(i, j);
			else
				gnd[i][j] = 10.0;
			getcolr(i, j);
		}
}

/*----------------------------------------------*/
/* calculate ground altitude given location i,j */
/*----------------------------------------------*/
 

 float getgnd(int i,int j) {
	//float getgnd(i, j) {
	return(XAMP*sin(i*XFREQ)*sin(j*YFREQ)+XAMP);
}

/*------------------------------------*/
/* calculate polygon color given i,j  */
/*------------------------------------*/

void getcolr(int i,int j) {
	float i1, i2, j1, j2;
	float cval, max;
	
	// use 2-way sum of sines in i and j 
	
	i1 = 0.2*sin(0.8*i);
	i2 = 0.3*sin(0.5*i);
	j1 = 0.2*sin(0.8*j);
	j2 = 0.3*sin(1.0*j);
	
	max = 0.2*0.2+0.3*0.3;
	
	cval = (i1*j1 + i2*j2 + max + 0.5);
	
	ccr[i][j] = GNDR * cval;
	ccg[i][j] = GNDG * cval;
	ccb[i][j] = GNDB * cval;
}

/*--------------------------------------*/
/* draw the ground polygons and a river */
/*--------------------------------------*/

void draw_ground(void) {
	int i, j;
	
	float xc, yc;
	
	xc = TMAXX/2.0*FPB;
	yc = TMAXY/2.0*FPB;
	
	// draw the river 
	glPushMatrix();
	glTranslatef(0.0, 0.0, -10.0);
//	glTranslatef(TOWER_N, TOWER_E, 0.0);
	glColor3f(0.14, 0.32, 0.73);

	
	GLfloat directedLight[] = {0.7f, 0.7f, 0.7f, 1.0f};
	GLfloat directedLightPos[] = {-10.0f, 15.0f, 20.0f, 0.0f};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, directedLight);
	glLightfv(GL_LIGHT0, GL_POSITION, directedLightPos);
	
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, _textureId_dam);
	glColor3f(1.0f, 1.0f, 1.0f);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,
					GL_TEXTURE_MIN_FILTER,
					GL_LINEAR_MIPMAP_LINEAR);
	

	glBegin(GL_POLYGON);
	{
		glTexCoord2f(1.0f, 0.0f);
		glVertex2f(xc-2000.0, yc-40000.0);
		glTexCoord2f(1.0f, 3.0f);
		glVertex2f(xc-000.0, yc+8000.0);
		glTexCoord2f(0.0f, 3.0f);
		glVertex2f(xc-4700.0, yc+8000.0);
		glTexCoord2f(0.0f, 0.0f);
		glVertex2f(xc-10200.0, yc-40000.0);
	}
	glEnd();
	

	
	glBegin(GL_POLYGON);
	{
		glTexCoord2f(1.0f, 3.0f);
		glVertex2f(xc-4700.0, yc+8000.0);
		glTexCoord2f(0.0f, 3.0f);
		glVertex2f(xc-000.0, yc+8000.0);
		glTexCoord2f(1.0f, 0.0f);
		glVertex2f(xc-6200.0, yc+22000.0);
		glTexCoord2f(0.0f, 0.0f);
		glVertex2f(xc-7500.0, yc+22000.0);
	}
	glEnd();
	
	glBegin(GL_POLYGON);
	{
		glTexCoord2f(1.0f, 0.0f);
		glVertex2f(xc-7500.0, yc+22000.0);
		glTexCoord2f(1.0f, 1.0f);
		glVertex2f(xc-6200.0, yc+22000.0);
		glTexCoord2f(0.0f, 1.0f);
		glVertex2f(xc-4300.0, yc+40000.0);
		glTexCoord2f(0.0f, 0.0f);
		glVertex2f(xc-4800.0, yc+40000.0);
	}
	glEnd();
 
	glDisable(GL_TEXTURE_2D);
 
	glPopMatrix();

	
	
    glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, _textureId_dam3);
	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,
					GL_TEXTURE_MIN_FILTER,
					GL_LINEAR_MIPMAP_LINEAR);
	
	
	// Go thru all blocks of terrain and draw them 
	for (j = 0; j < TMAXY; j++)
		for ( i = 0; i < TMAXX; i++)
			draw_rect(i,j);
} 

void  draw_rect(int i, int j) {
	float vv[3];
	
	vv[0] = i * FPB;
	vv[1] = j * FPB;
	vv[2] = gnd[i][j];
	
	glColor3f(ccr[i][j], ccg[i][j], ccb[i][j]);
	
	glBegin(GL_POLYGON);
	{
		glTexCoord2f(1.0f, 0.0f);
		glVertex3fv(vv);
		vv[0] = (i+1) * FPB;
		vv[2] = gnd[i+1][j];
		glTexCoord2f(1.0f, 1.0f);
		glVertex3fv(vv);
		vv[1] = (j+1) * FPB;
		vv[2] = gnd[i+1][j+1];
		glTexCoord2f(0.0f, 1.0f);
		glVertex3fv(vv);
		vv[0] = i * FPB;
		vv[2] = gnd[i][j+1];
		glTexCoord2f(0.0f, 0.0f);
		glVertex3fv(vv);
	}
	glEnd();
}

/*--------------------------------------------*/
/* draw another airplane circling the airport */
/*--------------------------------------------*/




void draw_new_plane(void) {
	
	glPushMatrix();

	glTranslatef(state[NORTH], state[EAST], state[DOWN]);
	glRotatef(state[PSI]*180/PI, 0.0, 0.0, 1.0);
	glRotatef(state[THETA]*180/PI, 0.0, 1.0, 0.0);
	glRotatef(state[PHI]*180/PI, 1.0, 0.0, 0.0);
	
	glTranslatef(cg_loc_x,0.0,0.0+cg_loc_z);

	// make it 2x bigger than modeled 
	glScalef(3, 3, 3);
	glColor3f(0.5, 0.5, .5);
	
	glCallList(display_list);
	
	glPopMatrix();
	
}



/*--------------------------------------------------------------*/
/* draw cockpit instrumentation                                 */
/* currently this is only a center indicator, rudder indicator, */
/* and some textual data for debugging                          */
/*--------------------------------------------------------------*/


void draw_cockpit(void) {
	int i;
	
	// cockpit instruments drawn in a 2-D orthographic projection 
	// first, we save away the perspective projection             
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	
	glLoadIdentity();
	gluOrtho2D(0.0, XMAXSCREEN, 0.0, YMAXSCREEN);
	glMatrixMode(GL_MODELVIEW);
	
	
	// draw rudder position indicator 
	glColor3f(.80, 0.0, 0.0);

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);   //DMT
	
	
	// draw text on the screen (toggle with the 't' key) 
	if (text_display)
		draw_text_readout();
	
	// return to the saved perspective projection 
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}


/*------------------------------*/
/* print the text on the screen */
/*------------------------------*/


void draw_text_readout(void) {
 double vel;
 char buffer[100];
 
 vel = sqrt(state[U]*state[U]+state[V]*state[V]+state[W]*state[W]);
 
 glColor3f(0.0, 0.0, 0.0);
 sprintf(buffer, "SPEED:%.1f fps  ALT:%.0f ft  V/S:%.0f fpm  ", 
 airspeed_read, -state[DOWN],
 airspeed_read*sin(state[THETA]-alpha_read)*60.0);
 glRasterPos2f(100.0,100.0);
 print_string1(buffer, H18);
 
 sprintf(buffer, "PITCH:%.2f deg ROLL:%.2f deg HEADING:%.2f deg  ",
 state[THETA]*180.0/PI, state[PHI]*180.0/PI, state[PSI]*180.0/PI);
 glRasterPos2f(100.0,80.0);
 print_string1(buffer, H18);
 
 sprintf(buffer,
 "ALPHA:%.1f deg  BETA:%.1f deg  ",
 alpha_read*180.0/PI,
 beta_read*180.0/PI);
 glRasterPos2f(100.0,60.0);
 print_string2(buffer, H18);
 
 sprintf(buffer, "Ail:%.2f deg Elev:%.2f deg  Rudd:%.2f deg     ",
 control[AILERON], control[ELEVATOR_2], control[RUDDER]);
 glRasterPos2f(100.0,40.0);
 print_string1(buffer, H18);
 
 sprintf(buffer, "throttle:%.0f   ",
 control[THROTTLE_L]);
 glRasterPos2f(100.0,20.0);
 print_string4(buffer, H18);
 }
 

/*--------------------------------------------*/
/* print a given string of text on the screen */
/*--------------------------------------------*/


void print_string1(char *s, void *fontname) {
	int i, len;
	int strlen(char *ss);
	
	len = 47; //strlen(s);
	for (i=0; i<len; i++)
		glutBitmapCharacter(fontname, s[i]);
}
void print_string2(char *s, void *fontname) {
	int i, len;
	int strlen(char *ss);
	
	len = 27; //strlen(s);
	for (i=0; i<len; i++)
		glutBitmapCharacter(fontname, s[i]);
}
void print_string3(char *s, void *fontname) {
	int i, len;
	int strlen(char *ss);
	
	len = 27; //strlen(s);
	for (i=0; i<len; i++)
		glutBitmapCharacter(fontname, s[i]);
}
void print_string4(char *s, void *fontname) {
	int i, len;
	int strlen(char *ss);
	
	len = 13; //strlen(s);
	for (i=0; i<len; i++)
		glutBitmapCharacter(fontname, s[i]);
}





void read_state(void) {
    double num1;
    double num2;
    double num3;
    double num4;
    double num5;
    double num6;
    double num7;
    double num8;
    double num9;
    double num10;
    double num11;
    double num12;
    double num13;
    double num14;
    double num15;
    double num16;
	 
	indata >> num1 >> num2 >> num3 >> num4 >> num5 >> num6 >> num7 >> num8 >> num9 >> num10 >> num11 >> num12 >> num13 >> num14 >> num15 >> num16;
	if(indata.eof() ){
		indata.close();
		cout << "End-of-file reached.." << endl;
		 
//		indata.open("animation_state.txt"); // opens the file
 
//		indata.open("animation_phugoid_2.txt"); // opens the file
//		indata.open("animation_phugoid_10.txt"); // opens the file
//		indata.open("animation_sp.txt"); // opens the file
//		indata.open("animation_dutch_roll.txt"); // opens the file
 //       indata.open("animation_dutch_roll2.txt"); // opens the file
//		indata.open("animation_airspeed_hold.txt"); // opens the file
//		indata.open("animation_airspeed_alt_hold.txt"); // opens the file
//		indata.open("animation_airspeed_hold_lowdamp.txt"); // opens the file
//		indata.open("animation_airspeed_hold_low_stab.txt"); // opens the file
//		indata.open("animation_airspeed_hold_low_stab_rate_feed.txt"); // opens the file
//		indata.open("animation_airspeed_alt_bank.txt"); // opens the file
//		indata.open("animation_airspeed_alt_bank_turn_coord.txt"); // opens the file
//		indata.open("animation_waypoint.txt"); // opens the file
//		indata.open("animation_waypoint_50.txt"); // opens the file
        indata.open("simRecording_20230225_1436.txt"); // opens the file

        
		if(!indata) { // file couldn't be opened
			cerr << "Error: file could not be opened" << endl;
			exit(0);
		}
		
		
//		exit(0);
	}	
	
	//cout << "The next number is " << state[NORTH] << state[EAST] << endl;
	
	state[PHI]=num1;
	state[THETA]=num2;
	state[PSI]=num3 + 0.0*PI/180.0;
	state[NORTH]=num4+240000;
	state[EAST]=num5+250000;
	state[DOWN]=-num6;
	TARGET_E = num7 + 250000;
	TARGET_N = num8 + 240000;
    TARGET_D = -num16;

	control[AILERON] = num9*180/PI;
	control[RUDDER] = num10*180/PI;
	control[ELEVATOR_2] = num11*180/PI;
	control[THROTTLE_L] = num12;
	control[THROTTLE_R] = num12;
	
	airspeed_read = num13;
	alpha_read = num14;
	beta_read = num15;
	 
	
}





/*----------------------------------*/
/* return # of seconds in last loop */
/*----------------------------------*/

double get_time_step(void) {
	double dt;
	clock_t nowtime;
	
	nowtime = gettime();
	dt = (nowtime - lasttime);
	lasttime = nowtime;
	
	return dt;
}







