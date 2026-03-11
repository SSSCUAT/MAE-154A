
  
// ============================================================
// flight_mae_154A.cpp  —  MAE 154A Flight Visualizer
// UCLA campus scene, Extra300 STL model, Simulink playback
// ============================================================

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <cstdlib>
#include <stdint.h>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "ac3d.h"
#include "imageloader.h"
#include <time.h>
#include <math.h>

#include "flight.h"
#include "graphics.h"
#include "vehicle.h"

// ============================================================
// STL MODEL LOADER  (ASCII + Binary)
// ============================================================

#define MAX_TRIS 200000

struct StlModel {
    float verts[MAX_TRIS][3][3];
    int   tri_count;
};

// ============================================================
// FORWARD DECLARATIONS
// ============================================================

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

float getgnd(int i, int j);
void setup_terrain(void);
void getcolr(int i, int j);
void draw_rect(int i, int j);
void draw_runway(void);
void draw_tower(void);
void draw_waypoint(void);
void draw_new_plane(void);
void draw_explosion(void);
void draw_stl_wireframe(StlModel& model);
void draw_ground(void);
void draw_airport(void);
void draw_ground_plane(void);
void draw_cockpit(void);
void draw_text_readout(void);
void print_string1(char *s, void *font);
void print_string2(char *s, void *font);
void print_string3(char *s, void *font);
void print_string4(char *s, void *font);
void draw_box(float w, float d, float h, float r, float g, float b);
void draw_dome(float radius, float r, float g, float b);

double get_time_step(void);

StlModel uav_model;

void load_stl(const char* filename, StlModel& model) {
    model.tri_count = 0;
    FILE* f = fopen(filename, "rb");
    if (!f) {
        printf("WARNING: Could not open %s\n", filename);
        return;
    }
    char header[6] = {0};
    fread(header, 1, 5, f);
    rewind(f);

    if (strncmp(header, "solid", 5) == 0) {
        // ASCII STL
        char line[256];
        float x, y, z;
        int vert_in_tri = 0;
        while (fgets(line, 256, f)) {
            char* l = line;
            while (*l == ' ' || *l == '\t') l++;
            if (strncmp(l, "vertex", 6) == 0) {
                sscanf(l + 7, "%f %f %f", &x, &y, &z);
                model.verts[model.tri_count][vert_in_tri][0] = x;
                model.verts[model.tri_count][vert_in_tri][1] = y;
                model.verts[model.tri_count][vert_in_tri][2] = z;
                vert_in_tri++;
                if (vert_in_tri == 3) { model.tri_count++; vert_in_tri = 0; }
            }
        }
    } else {
        // Binary STL
        fseek(f, 80, SEEK_SET);
        uint32_t num_tris;
        fread(&num_tris, 4, 1, f);
        for (uint32_t i = 0; i < num_tris && i < MAX_TRIS; i++) {
            float normal[3];
            fread(normal, 4, 3, f);
            for (int v = 0; v < 3; v++) {
                float vert[3];
                fread(vert, 4, 3, f);
                model.verts[i][v][0] = vert[0];
                model.verts[i][v][1] = vert[1];
                model.verts[i][v][2] = vert[2];
            }
            uint16_t attr;
            fread(&attr, 2, 1, f);
            model.tri_count++;
        }
    }
    fclose(f);
    printf("Loaded %s: %d triangles\n", filename, model.tri_count);
}

void draw_stl(StlModel& model) {
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < model.tri_count; i++) {
        glVertex3fv(model.verts[i][0]);
        glVertex3fv(model.verts[i][1]);
        glVertex3fv(model.verts[i][2]);
    }
    glEnd();
}

void draw_stl_wireframe(StlModel& model) {
    // Draw edges of the STL as wireframe
    glBegin(GL_LINES);
    for (int i = 0; i < model.tri_count; i++) {
        // Draw three edges of each triangle
        glVertex3fv(model.verts[i][0]);
        glVertex3fv(model.verts[i][1]);
        glVertex3fv(model.verts[i][1]);
        glVertex3fv(model.verts[i][2]);
        glVertex3fv(model.verts[i][2]);
        glVertex3fv(model.verts[i][0]);
    }
    glEnd();
}

// ============================================================
// GLOBALS
// ============================================================

ifstream indata;
double airspeed_read;
double alpha_read;
double beta_read;
double TARGET_N;
double TARGET_E;
double TARGET_D;
double cg_loc_x;
double cg_loc_z;

int display_list;

GLuint _textureId_dam;
GLuint _textureId_dam3;

// Plane visualization modes
int wireframe_mode = 0;  // 0=solid with edges, 1=wireframe only, 2=solid only

// Explosion state
int explosion_active = 0;
double explosion_time = 0.0;
double explosion_start_time = 0.0;
double explosion_x, explosion_y, explosion_z;
#define EXPLOSION_DURATION 5.5  // seconds

// Bombing state
int was_at_low_altitude = 0;
int bomb_dropped = 0;
double bomb_location_n, bomb_location_e, bomb_location_d;  // Where the bomb was dropped
#define LOW_ALTITUDE_THRESHOLD 400.0
#define ESCAPE_ALTITUDE 600.0

GLuint loadMipmappedTexture(Image *image) {
    GLuint id;
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_2D, id);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB,
                      image->width, image->height,
                      GL_RGB, GL_UNSIGNED_BYTE, image->pixels);
    return id;
}

// ============================================================
// MAIN
// ============================================================

int main(int argc, char** argv) {
    ACObject *ob;
    static char acFileName[] = "pioneer_body_tex.ac";

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(200, 200);
    glutInitWindowPosition(10, 10);
    glutCreateWindow("Flight");

    initial_states();
    glutFullScreen();
    setup_graphics();

    lasttime = gettime();

    glutDisplayFunc(display_and_dynamics);
    glutReshapeFunc(reshape);
    glutMouseFunc(mousebutton);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(arrowkeys);
    glutPassiveMotionFunc(mousepos);

    ob = ac_load_ac3d(acFileName);
    if (ob != NULL)
        display_list = ac_display_list_render_object(ob);
    else
        printf("WARNING: Could not load %s - using STL only\n", acFileName);

    glutMainLoop();
    return 0;
}

// ============================================================
// INITIAL STATES
// ============================================================

void initial_states(void) {
    zoom     = 1.0;
    pan      = 0.0;
    sway     = 0.0;
    view     = OUT_THE_WINDOW;
    cg_loc_x = 13.0;
    cg_loc_z = 1.0;

    load_stl("Prototype3part.stl", uav_model);

    indata.open("simRecording_20260309_2221.txt");
    if (!indata) {
        cerr << "Error: could not open simRecording_20260309_2221.txt" << endl;
        exit(0);
    }
}

// ============================================================
// MAIN LOOP
// ============================================================

void display_and_dynamics(void) {
    draw_view();
    read_state();
    glutSwapBuffers();
    Sleep(1);
    glutPostRedisplay();
}

// ============================================================
// RESHAPE
// ============================================================

void reshape(int w, int h) {
    XMAXSCREEN   = (double)(w);
    YMAXSCREEN   = (double)(h);
    XCTR         = XMAXSCREEN / 2.0;
    YCTR         = YMAXSCREEN / 2.0;
    ASPECT_RATIO = XMAXSCREEN / YMAXSCREEN;
    glViewport(0, 0, (GLsizei)(w), (GLsizei)(h));
    perspective_projection();
}

// ============================================================
// INPUT HANDLERS
// ============================================================

void mousebutton(int button, int state, int x, int y) {
    (void)button; (void)state; (void)x; (void)y;
}

void keyboard(unsigned char key, int x, int y) {
    (void)x; (void)y;
    switch (key) {
    case 27:  exit(0); break;
    case 't': text_display = !text_display; break;
    case 'q':
        if      (view == OUT_THE_WINDOW) view = CHASEPLANE;
        else if (view == CHASEPLANE)     view = TOWER;
        else                             view = OUT_THE_WINDOW;
        break;
    case '[': zoom /= 1.1; break;
    case ']': zoom *= 1.1; break;
    case 'w':
        wireframe_mode = (wireframe_mode + 1) % 3;
        printf("Wireframe mode: %d (0=solid+edges, 1=solid+wireframe, 2=wireframe)\n", wireframe_mode);
        break;
    default: break;
    }
}

void arrowkeys(int key, int x, int y) {
    (void)x; (void)y;
    switch (key) {
    case GLUT_KEY_LEFT:  pan  += 5.0; break;
    case GLUT_KEY_RIGHT: pan  -= 5.0; break;
    case GLUT_KEY_UP:    sway -= 5.0; break;
    case GLUT_KEY_DOWN:  sway += 5.0; break;
    default: break;
    }
}

void mousepos(int x, int y) { (void)x; (void)y; }

// ============================================================
// UTILITIES
// ============================================================

clock_t gettime(void) { return clock(); }

double limit(double val, double low, double hi) {
    if (val > hi)  return hi;
    if (val < low) return low;
    return val;
}

double get_time_step(void) {
    clock_t nowtime = gettime();
    double dt = (double)(nowtime - lasttime);
    lasttime = nowtime;
    return dt;
}

// ============================================================
// SETUP GRAPHICS
// ============================================================

void setup_graphics(void) {
    float fogcolor[] = {FOGR, FOGG, FOGB};

    setup_terrain();

    glClearColor(FOGR, FOGG, FOGB, 0.0f);
    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    Image* image = loadBMP("ocean_v3_256.bmp");
    _textureId_dam = loadMipmappedTexture(image);
    delete image;

    image = loadBMP("mountain_tex.bmp");
    _textureId_dam3 = loadMipmappedTexture(image);
    delete image;

    glEnable(GL_FOG);
    glFogi(GL_FOG_MODE,    GL_EXP);
    glFogfv(GL_FOG_COLOR,  fogcolor);
    glFogf(GL_FOG_DENSITY, FOG_DENSITY);
    glHint(GL_FOG_HINT,    GL_DONT_CARE);
    glFogf(GL_FOG_START,   FOG_START);
    glFogf(GL_FOG_END,     FOG_END);

    text_display = TRUE;
}

// ============================================================
// PERSPECTIVE PROJECTION
// ============================================================

void perspective_projection(void) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(FOV, ASPECT_RATIO, 10.0, 500000.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

// ============================================================
// DRAW VIEW  —  set camera, clear, draw world, draw plane, draw HUD
// ============================================================

void draw_view(void) {

    glPushMatrix();

    if (view == OUT_THE_WINDOW) {
        glRotatef( 90.0f, 0.0f, 1.0f, 0.0f);
        glRotatef( 90.0f, 1.0f, 0.0f, 0.0f);
        glRotatef((float)(-state[PHI]   * 180.0 / PI), 1.0f, 0.0f, 0.0f);
        glRotatef((float)(-state[THETA] * 180.0 / PI), 0.0f, 1.0f, 0.0f);
        glRotatef((float)(-state[PSI]   * 180.0 / PI), 0.0f, 0.0f, 1.0f);
        glTranslatef((float)(-state[NORTH]), (float)(-state[EAST]), (float)(-state[DOWN]));
    }
    else if (view == CHASEPLANE) {
        float dist = 200.0f * (float)zoom;
        float camN = (float)state[NORTH] - dist * (float)cos(state[PSI]);
        float camE = (float)state[EAST]  - dist * (float)sin(state[PSI]);
        float camD = (float)state[DOWN]  - 50.0f;
        gluLookAt(camN, camE, camD,
                  state[NORTH], state[EAST], state[DOWN],
                  0.0, 0.0, -1.0);
    }
    else {   // TOWER view
        glRotatef( 90.0f, 0.0f, 1.0f, 0.0f);
        glRotatef( 90.0f, 1.0f, 0.0f, 0.0f);
        float dN = (float)(state[NORTH] - TOWER_N);
        float dE = (float)(state[EAST]  - TOWER_E);
        float dD = (float)(state[DOWN]  - TOWER_D);
        glRotatef( atan2(dD, sqrtf(dN*dN + dE*dE)) * 180.0f / (float)PI,
                   0.0f, 1.0f, 0.0f);
        glRotatef(-atan2(dE, dN) * 180.0f / (float)PI,
                   0.0f, 0.0f, 1.0f);
        glTranslatef(-(float)TOWER_N, -(float)TOWER_E, -(float)TOWER_D + 10.0f);
    }

    // clear after camera is set so sky color is correct
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    GLfloat ambientLight[] = {0.8f, 0.8f, 0.8f, 1.0f};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);

    draw_ground_plane();   // sky blue backdrop
    draw_airport();        // UCLA campus buildings + target person
    draw_ground();         // colored terrain tiles

    if (view != OUT_THE_WINDOW)
        draw_new_plane();

    draw_explosion();      // explosion effect if active

    glPopMatrix();

    draw_cockpit();        // 2D HUD text — outside push/pop
}

// ============================================================
// SKY BLUE BACKDROP
// ============================================================

void draw_ground_plane(void) {
    glDisable(GL_DEPTH_TEST);
    glColor3f(0.53f, 0.81f, 0.98f);
    glPushMatrix();
    glTranslatef(0.0f, 0.0f, 20.0f);
    glRectf(-500000.0f, -500000.0f, 500000.0f, 500000.0f);
    glPopMatrix();
    glEnable(GL_DEPTH_TEST);
}

// ============================================================
// BUILDING HELPERS
// ============================================================

// Solid box. Origin at ground center. Top face at z = -h (up in NED).
void draw_box(float w, float d, float h, float r, float g, float b) {
    float hw = w * 0.5f, hd = d * 0.5f;
    // top face
    glColor3f(r, g, b);
    glBegin(GL_QUADS);
        glVertex3f(-hw, -hd, -h);
        glVertex3f( hw, -hd, -h);
        glVertex3f( hw,  hd, -h);
        glVertex3f(-hw,  hd, -h);
    glEnd();
    // north wall
    glColor3f(r*0.82f, g*0.82f, b*0.82f);
    glBegin(GL_QUADS);
        glVertex3f( hw, -hd,  0.0f);
        glVertex3f(-hw, -hd,  0.0f);
        glVertex3f(-hw, -hd,    -h);
        glVertex3f( hw, -hd,    -h);
    glEnd();
    // south wall
    glBegin(GL_QUADS);
        glVertex3f(-hw,  hd,  0.0f);
        glVertex3f( hw,  hd,  0.0f);
        glVertex3f( hw,  hd,    -h);
        glVertex3f(-hw,  hd,    -h);
    glEnd();
    // east wall (darker)
    glColor3f(r*0.68f, g*0.68f, b*0.68f);
    glBegin(GL_QUADS);
        glVertex3f(hw, -hd,  0.0f);
        glVertex3f(hw,  hd,  0.0f);
        glVertex3f(hw,  hd,    -h);
        glVertex3f(hw, -hd,    -h);
    glEnd();
    // west wall
    glBegin(GL_QUADS);
        glVertex3f(-hw,  hd,  0.0f);
        glVertex3f(-hw, -hd,  0.0f);
        glVertex3f(-hw, -hd,    -h);
        glVertex3f(-hw,  hd,    -h);
    glEnd();
}

// Hemisphere dome. Origin at base. Apex at z = -radius.
void draw_dome(float radius, float r, float g, float b) {
    int slices = 16, stacks = 8;
    glColor3f(r, g, b);
    for (int ii = 0; ii < stacks; ii++) {
        float phi1 = (float)ii     / (float)stacks * (float)PI * 0.5f;
        float phi2 = (float)(ii+1) / (float)stacks * (float)PI * 0.5f;
        glBegin(GL_QUAD_STRIP);
        for (int jj = 0; jj <= slices; jj++) {
            float theta = (float)jj / (float)slices * 2.0f * (float)PI;
            glVertex3f(radius*cosf(phi1)*cosf(theta),
                       radius*cosf(phi1)*sinf(theta),
                      -radius*sinf(phi1));
            glVertex3f(radius*cosf(phi2)*cosf(theta),
                       radius*cosf(phi2)*sinf(theta),
                      -radius*sinf(phi2));
        }
        glEnd();
    }
}

// ============================================================
// UCLA CAMPUS
// cx = campus center NORTH in world coords = TMAXX/2 * FPB
// cy = campus center EAST  in world coords = TMAXY/2 * FPB
// Offsets in feet from that center.
// ============================================================

void draw_airport(void) {
    float cx = (float)TMAXX * 0.5f * (float)FPB;
    float cy = (float)TMAXY * 0.5f * (float)FPB;
    int k;

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);

    // ========== MAJOR STREETS GRID ==========
    // North-South major streets (vertical)
    glColor3f(0.28f, 0.28f, 0.26f);
    for (float n = -45000.0f; n <= 45000.0f; n += 2000.0f) {
        glPushMatrix();
        glTranslatef(cx + n, cy, 1.0f);
        glBegin(GL_QUADS);
            glVertex3f(-50.0f, -35000.0f, 0.0f);
            glVertex3f( 50.0f, -35000.0f, 0.0f);
            glVertex3f( 50.0f,  35000.0f, 0.0f);
            glVertex3f(-50.0f,  35000.0f, 0.0f);
        glEnd();
        glPopMatrix();
    }
    // East-West major streets (horizontal)
    for (float e = -35000.0f; e <= 35000.0f; e += 2000.0f) {
        glPushMatrix();
        glTranslatef(cx, cy + e, 1.0f);
        glBegin(GL_QUADS);
            glVertex3f(-35000.0f, -50.0f, 0.0f);
            glVertex3f( 35000.0f, -50.0f, 0.0f);
            glVertex3f( 35000.0f,  50.0f, 0.0f);
            glVertex3f(-35000.0f,  50.0f, 0.0f);
        glEnd();
        glPopMatrix();
    }

    // ========== FILL ENTIRE MAP WITH DENSE RESIDENTIAL/COMMERCIAL ==========
    // Northern areas (Valley and Foothills)
    for (int i = -22; i < 22; i++) {
        for (int j = -22; j < -8; j++) {
            float bx = cx + i * 1500.0f;
            float by = cy + j * 1600.0f;
            float h = 30.0f + ((i + j) % 5) * 15.0f;
            float r = 0.50f + ((i * j) % 10) * 0.01f;
            glPushMatrix();
            glTranslatef(bx, by, 0.0f);
            draw_box(300.0f + (i % 3) * 30, 280.0f + (j % 3) * 25, h, r, r - 0.04f, r - 0.08f);
            glPopMatrix();
        }
    }

    // Eastern areas (Downtown & East LA)
    for (int i = 8; i < 30; i++) {
        for (int j = -18; j < 22; j++) {
            float bx = cx + i * 1400.0f;
            float by = cy + j * 1500.0f;
            float h = 50.0f + ((i + j) % 6) * 20.0f;
            float r = 0.48f + ((i * j) % 15) * 0.01f;
            glPushMatrix();
            glTranslatef(bx, by, 0.0f);
            draw_box(330.0f + (i % 4) * 20, 300.0f + (j % 4) * 20, h, r, r - 0.05f, r - 0.09f);
            glPopMatrix();
        }
    }

    // Western areas (Westwood & Santa Monica)
    for (int i = -28; i < -8; i++) {
        for (int j = -5; j < 20; j++) {
            float bx = cx + i * 1400.0f;
            float by = cy + j * 1500.0f;
            float h = 40.0f + ((i + j) % 5) * 15.0f;
            float r = 0.52f + ((i * j) % 12) * 0.01f;
            glPushMatrix();
            glTranslatef(bx, by, 0.0f);
            draw_box(310.0f + (i % 3) * 25, 290.0f + (j % 3) * 22, h, r, r - 0.04f, r - 0.07f);
            glPopMatrix();
        }
    }

    // Central areas (around UCLA)
    for (int i = -8; i < 8; i++) {
        for (int j = -5; j < 10; j++) {
            // Skip UCLA core area - will be drawn with custom buildings
            if (i >= -6 && i <= -2 && j >= 1 && j <= 7) continue;
            
            float bx = cx + i * 1500.0f;
            float by = cy + j * 1500.0f;
            float h = 35.0f + ((i + j) % 4) * 12.0f;
            float r = 0.51f + ((i * j) % 10) * 0.01f;
            glPushMatrix();
            glTranslatef(bx, by, 0.0f);
            draw_box(320.0f, 300.0f, h, r, r - 0.04f, r - 0.08f);
            glPopMatrix();
        }
    }

    // ========== UCLA CAMPUS (ORIGINAL ENHANCED) ==========
    // Main green quad
    glColor3f(0.22f, 0.50f, 0.18f);
    glPushMatrix();
    glTranslatef(cx - 4500.0f, cy + 4500.0f, 0.0f);
    glBegin(GL_QUADS);
        glVertex3f(-2000.0f, -1800.0f, 0.0f);
        glVertex3f( 2000.0f, -1800.0f, 0.0f);
        glVertex3f( 2000.0f,  1800.0f, 0.0f);
        glVertex3f(-2000.0f,  1800.0f, 0.0f);
    glEnd();
    glPopMatrix();

    // ------ ROYCE HALL — brick red, twin bell towers ------
    // main body
    glPushMatrix();
    glTranslatef(cx - 3500.0f, cy + 3000.0f, 0.0f);
    draw_box(600.0f, 400.0f, 120.0f, 0.72f, 0.42f, 0.18f);
    glPopMatrix();
    // west tower
    glPushMatrix();
    glTranslatef(cx - 2500.0f, cy + 2200.0f, 0.0f);
    draw_box(250.0f, 250.0f, 260.0f, 0.65f, 0.36f, 0.14f);
    glPopMatrix();
    // east tower
    glPushMatrix();
    glTranslatef(cx - 2500.0f, cy + 3800.0f, 0.0f);
    draw_box(250.0f, 250.0f, 260.0f, 0.65f, 0.36f, 0.14f);
    glPopMatrix();

    // ------ POWELL LIBRARY — brick + copper dome ------
    glPushMatrix();
    glTranslatef(cx - 4500.0f, cy + 8000.0f, 0.0f);
    draw_box(500.0f, 500.0f, 110.0f, 0.72f, 0.42f, 0.18f);
    // octagonal drum on top
    glTranslatef(0.0f, 0.0f, -110.0f);
    draw_box(320.0f, 320.0f, 70.0f, 0.68f, 0.38f, 0.16f);
    // copper green dome
    glTranslatef(0.0f, 0.0f, -70.0f);
    draw_dome(150.0f, 0.28f, 0.50f, 0.30f);
    glPopMatrix();

    // ------ KERCKHOFF HALL ------
    glPushMatrix();
    glTranslatef(cx - 2500.0f, cy + 5500.0f, 0.0f);
    draw_box(250.0f, 200.0f, 130.0f, 0.60f, 0.38f, 0.18f);
    glPopMatrix();

    // ------ MURPHY HALL — tall admin ------
    glPushMatrix();
    glTranslatef(cx - 3500.0f, cy + 8500.0f, 0.0f);
    draw_box(400.0f, 220.0f, 200.0f, 0.50f, 0.50f, 0.46f);
    glPopMatrix();

    // ------ BOELTER HALL (engineering) ------
    glPushMatrix();
    glTranslatef(cx - 2000.0f, cy + 6000.0f, 0.0f);
    draw_box(450.0f, 350.0f, 160.0f, 0.50f, 0.48f, 0.44f);
    glPopMatrix();

    // ------ BUNCHE HALL ------
    glPushMatrix();
    glTranslatef(cx - 1000.0f, cy + 5000.0f, 0.0f);
    draw_box(300.0f, 250.0f, 180.0f, 0.46f, 0.46f, 0.42f);
    glPopMatrix();

    // ------ MATH SCIENCES ------
    glPushMatrix();
    glTranslatef(cx - 3000.0f, cy + 6500.0f, 0.0f);
    draw_box(350.0f, 280.0f, 140.0f, 0.48f, 0.48f, 0.44f);
    glPopMatrix();

    // ------ ACKERMAN UNION ------
    glPushMatrix();
    glTranslatef(cx - 5500.0f, cy + 5000.0f, 0.0f);
    draw_box(500.0f, 600.0f, 110.0f, 0.54f, 0.52f, 0.46f);
    glPopMatrix();

    // ------ LUSKIN CONFERENCE CENTER ------
    glPushMatrix();
    glTranslatef(cx - 2000.0f, cy + 8000.0f, 0.0f);
    draw_box(350.0f, 300.0f, 150.0f, 0.38f, 0.50f, 0.58f);
    glPopMatrix();

    // ------ BRUIN BEAR STATUE ------
    glPushMatrix();
    glTranslatef(cx - 4000.0f, cy + 4500.0f, 0.0f);
    draw_box(80.0f, 80.0f, 50.0f, 0.52f, 0.42f, 0.20f);   // pedestal
    glTranslatef(0.0f, 0.0f, -50.0f);
    draw_box(100.0f, 60.0f, 75.0f, 0.38f, 0.28f, 0.10f);   // body
    glTranslatef(-20.0f, 0.0f, -70.0f);
    draw_box(50.0f, 45.0f, 45.0f, 0.34f, 0.24f, 0.08f);   // head
    glPopMatrix();

    // ------ JANSS STEPS — stone staircase ------
    glColor3f(0.78f, 0.72f, 0.60f);
    for (k = 0; k < 15; k++) {
        glPushMatrix();
        glTranslatef(cx - 2500.0f - k * 50.0f, cy + 3500.0f, 0.0f);
        glBegin(GL_QUADS);
            glVertex3f(-30.0f, -300.0f, (float)k * 8.0f);
            glVertex3f( 30.0f, -300.0f, (float)k * 8.0f);
            glVertex3f( 30.0f,  300.0f, (float)k * 8.0f);
            glVertex3f(-30.0f,  300.0f, (float)k * 8.0f);
        glEnd();
        glPopMatrix();
    }

    // ------ INVERTED FOUNTAIN ------
    glPushMatrix();
    glTranslatef(cx - 5000.0f, cy + 4500.0f, 0.0f);
    glColor3f(0.18f, 0.44f, 0.76f);
    glBegin(GL_TRIANGLE_FAN);
        glVertex3f(0.0f, 0.0f, 0.0f);
        for (k = 0; k <= 32; k++) {
            float a = (float)k * 2.0f * (float)PI / 32.0f;
            glVertex3f(cosf(a)*180.0f, sinf(a)*180.0f, 0.0f);
        }
    glEnd();
    glColor3f(0.38f, 0.62f, 0.88f);
    glBegin(GL_QUAD_STRIP);
    for (k = 0; k <= 32; k++) {
        float a = (float)k * 2.0f * (float)PI / 32.0f;
        glVertex3f(cosf(a)*200.0f, sinf(a)*200.0f, 0.0f);
        glVertex3f(cosf(a)*180.0f, sinf(a)*180.0f, 0.0f);
    }
    glEnd();
    glPopMatrix();

    // ========== WATER FEATURES ==========
    // Santa Monica Bay (south boundary)
    glColor3f(0.10f, 0.26f, 0.60f);
    glPushMatrix();
    glTranslatef(cx, cy + 32000.0f, 0.0f);
    glBegin(GL_QUADS);
        glVertex3f(-50000.0f, -3000.0f, 0.0f);
        glVertex3f( 50000.0f, -3000.0f, 0.0f);
        glVertex3f( 50000.0f,  3000.0f, 0.0f);
        glVertex3f(-50000.0f,  3000.0f, 0.0f);
    glEnd();
    glPopMatrix();

    // Los Angeles River
    glColor3f(0.20f, 0.40f, 0.70f);
    glPushMatrix();
    glTranslatef(cx + 15000.0f, cy - 5000.0f, 0.0f);
    glBegin(GL_QUADS);
        glVertex3f(-150.0f, -20000.0f, 0.0f);
        glVertex3f( 150.0f, -20000.0f, 0.0f);
        glVertex3f( 150.0f,  20000.0f, 0.0f);
        glVertex3f(-150.0f,  20000.0f, 0.0f);
    glEnd();
    glPopMatrix();

    // ========== MAJOR HIGHWAYS ==========
    // I-405 (north-south)
    glColor3f(0.32f, 0.32f, 0.30f);
    glPushMatrix();
    glTranslatef(cx - 8000.0f, cy, 1.0f);
    glBegin(GL_QUADS);
        glVertex3f(-180.0f, -35000.0f, 0.0f);
        glVertex3f( 180.0f, -35000.0f, 0.0f);
        glVertex3f( 180.0f,  35000.0f, 0.0f);
        glVertex3f(-180.0f,  35000.0f, 0.0f);
    glEnd();
    glPopMatrix();

    // I-10 (east-west)
    glPushMatrix();
    glTranslatef(cx + 12000.0f, cy - 10000.0f, 1.0f);
    glBegin(GL_QUADS);
        glVertex3f(-35000.0f, -180.0f, 0.0f);
        glVertex3f( 35000.0f, -180.0f, 0.0f);
        glVertex3f( 35000.0f,  180.0f, 0.0f);
        glVertex3f(-35000.0f,  180.0f, 0.0f);
    glEnd();
    glPopMatrix();

    glEnable(GL_LIGHTING);

    draw_waypoint();
}

void draw_runway(void) { /* replaced by UCLA campus */ }
void draw_tower(void)  { /* replaced by UCLA campus */ }

// ============================================================
// TARGET STICK FIGURE
// ============================================================

void draw_waypoint(void) {
    glPushMatrix();
    glTranslatef((float)TARGET_N, (float)TARGET_E, (float)TARGET_D);

    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    // shadow circle
    glColor3f(0.10f, 0.10f, 0.10f);
    glBegin(GL_TRIANGLE_FAN);
        glVertex3f(0.0f, 0.0f, 1.0f);
        for (int k = 0; k <= 16; k++) {
            float a = (float)k * 2.0f * (float)PI / 16.0f;
            glVertex3f(cosf(a)*22.0f, sinf(a)*22.0f, 1.0f);
        }
    glEnd();

    glEnable(GL_DEPTH_TEST);

    glColor3f(1.0f, 0.12f, 0.12f);
    glLineWidth(5.0f);

    // head
    glPushMatrix();
    glTranslatef(0.0f, 0.0f, -38.0f);
    glBegin(GL_QUADS);
        glVertex3f(-9.0f, -9.0f, 0.0f);
        glVertex3f( 9.0f, -9.0f, 0.0f);
        glVertex3f( 9.0f,  9.0f, 0.0f);
        glVertex3f(-9.0f,  9.0f, 0.0f);
    glEnd();
    glPopMatrix();

    // body
    glBegin(GL_LINES);
        glVertex3f(  0.0f, 0.0f, -27.0f);
        glVertex3f(  0.0f, 0.0f,  -7.0f);   // torso
        glVertex3f(-15.0f, 0.0f, -22.0f);
        glVertex3f( 15.0f, 0.0f, -22.0f);   // arms
        glVertex3f(  0.0f, 0.0f,  -7.0f);
        glVertex3f(-10.0f, 0.0f,   0.0f);   // left leg
        glVertex3f(  0.0f, 0.0f,  -7.0f);
        glVertex3f( 10.0f, 0.0f,   0.0f);   // right leg
    glEnd();

    glLineWidth(1.0f);
    glEnable(GL_LIGHTING);
    glPopMatrix();
}

// ============================================================
// EXPLOSION EFFECT
// ============================================================

void draw_explosion(void) {
    if (!explosion_active) return;

    glPushMatrix();
    glTranslatef((float)explosion_x, (float)explosion_y, (float)explosion_z);

    // Calculate explosion progress (0 to 1)
    double progress = explosion_time / EXPLOSION_DURATION;
    if (progress > 1.0) progress = 1.0;

    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);  // Ensure visibility

    // Expanding shockwave sphere - outer bright yellow flash
    double radius = 2000.0f + progress * 8000.0f;  // Much larger
    double alpha = 1.0f - progress;  // fade out

    glColor4f(1.0f, 1.0f, 0.3f, (float)(alpha * 0.7f));  // Bright yellow
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    int slices = 20, stacks = 12;  // More segments for detail
    for (int ii = 0; ii < stacks; ii++) {
        float phi1 = (float)ii     / (float)stacks * (float)PI * 0.5f;
        float phi2 = (float)(ii+1) / (float)stacks * (float)PI * 0.5f;
        glBegin(GL_QUAD_STRIP);
        for (int jj = 0; jj <= slices; jj++) {
            float theta = (float)jj / (float)slices * 2.0f * (float)PI;
            glVertex3f(radius*cosf(phi1)*cosf(theta),
                       radius*cosf(phi1)*sinf(theta),
                      -radius*sinf(phi1));
            glVertex3f(radius*cosf(phi2)*cosf(theta),
                       radius*cosf(phi2)*sinf(theta),
                      -radius*sinf(phi2));
        }
        glEnd();
    }

    // Middle layer - blue-purple for color variety
    double mid_radius = 1600.0f + progress * 6000.0f;
    glColor4f(0.4f, 0.2f, 1.0f, (float)(alpha * 0.4f));  // Purple-blue
    for (int ii = 0; ii < stacks; ii++) {
        float phi1 = (float)ii     / (float)stacks * (float)PI * 0.5f;
        float phi2 = (float)(ii+1) / (float)stacks * (float)PI * 0.5f;
        glBegin(GL_QUAD_STRIP);
        for (int jj = 0; jj <= slices; jj++) {
            float theta = (float)jj / (float)slices * 2.0f * (float)PI;
            glVertex3f(mid_radius*cosf(phi1)*cosf(theta),
                       mid_radius*cosf(phi1)*sinf(theta),
                      -mid_radius*sinf(phi1));
            glVertex3f(mid_radius*cosf(phi2)*cosf(theta),
                       mid_radius*cosf(phi2)*sinf(theta),
                      -mid_radius*sinf(phi2));
        }
        glEnd();
    }

    // Inner red-orange core
    double core_radius = 1200.0f + progress * 4500.0f;  // Much larger core
    glColor4f(1.0f, 0.2f, 0.0f, (float)(alpha * 0.6f));  // Red-orange
    for (int ii = 0; ii < stacks; ii++) {
        float phi1 = (float)ii     / (float)stacks * (float)PI * 0.5f;
        float phi2 = (float)(ii+1) / (float)stacks * (float)PI * 0.5f;
        glBegin(GL_QUAD_STRIP);
        for (int jj = 0; jj <= slices; jj++) {
            float theta = (float)jj / (float)slices * 2.0f * (float)PI;
            glVertex3f(core_radius*cosf(phi1)*cosf(theta),
                       core_radius*cosf(phi1)*sinf(theta),
                      -core_radius*sinf(phi1));
            glVertex3f(core_radius*cosf(phi2)*cosf(theta),
                       core_radius*cosf(phi2)*sinf(theta),
                      -core_radius*sinf(phi2));
        }
        glEnd();
    }

    // Burst particles - multi-colored (yellow to red)
    double particle_radius = 1000.0f + progress * 5000.0f;  // Much larger particle burst
    glBegin(GL_LINES);
    glLineWidth(2.0f);
    for (int i = 0; i < 16; i++) {
        float angle = (float)i * 2.0f * (float)PI / 16.0f;
        float distance = particle_radius;
        float color_blend = (float)i / 16.0f;
        // Alternate between yellow and red particles
        if (i % 2 == 0) {
            glColor4f(1.0f, 0.8f - color_blend * 0.3f, 0.0f, (float)(alpha * 0.5f));  // Yellow-red
        } else {
            glColor4f(1.0f, 0.2f + color_blend * 0.6f, 0.5f, (float)(alpha * 0.4f));  // Red-pink
        }
        // Radial burst lines
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(distance * cosf(angle), distance * sinf(angle), -distance * 0.5f);
    }
    glLineWidth(1.0f);
    glEnd();

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glPopMatrix();
}

// ============================================================
// TERRAIN
// ============================================================

void setup_terrain(void) {
    for (int i = 0; i < TMAXX; i++)
        for (int j = 0; j < TMAXY; j++) {
            gnd[i][j] = getgnd(i, j);
            getcolr(i, j);
        }
}

float getgnd(int i, int j) {
    float xc   = (float)TMAXX * 0.5f;
    float yc   = (float)TMAXY * 0.5f;
    float fi   = (float)i;
    float fj   = (float)j;
    float dist = sqrtf((fi-xc)*(fi-xc) + (fj-yc)*(fj-yc));

    if (dist < 14.0f) return 2.0f;   // flat campus core

    // Santa Monica Mountains — north = low i in NED terrain
    float hill = 0.0f;
    if (fi < xc - 12.0f) {
        float t = (xc - 12.0f - fi) / (xc - 12.0f);
        hill = (float)XAMP * 0.75f * t * t;
    }

    float urban = (float)XAMP * 0.15f
                * sinf(fi * (float)XFREQ * 0.6f)
                * sinf(fj * (float)YFREQ * 0.6f);

    return urban + hill;
}

void getcolr(int i, int j) {
    float xc   = (float)TMAXX * 0.5f;
    float yc   = (float)TMAXY * 0.5f;
    float fi   = (float)i;
    float fj   = (float)j;

    // UCLA campus core - green
    if (fi > xc - 6.0f && fi < xc - 2.0f && fj > yc + 2.0f && fj < yc + 8.0f) {
        ccr[i][j] = 0.22f; ccg[i][j] = 0.50f; ccb[i][j] = 0.18f;
    }
    // Santa Monica Bay (south)
    else if (fj > yc + 24.0f) {
        ccr[i][j] = 0.10f; ccg[i][j] = 0.26f; ccb[i][j] = 0.60f;
    }
    // Dense urban everywhere else - grey
    else {
        float urbanDensity = 0.50f + ((fi + fj) * 0.001f);
        urbanDensity = (urbanDensity > 0.54f) ? 0.54f : urbanDensity;
        ccr[i][j] = urbanDensity;
        ccg[i][j] = urbanDensity - 0.04f;
        ccb[i][j] = urbanDensity - 0.08f;
    }
}

// ============================================================
// DRAW GROUND TILES
// ============================================================

void draw_ground(void) {
    GLfloat directedLight[]    = {0.7f, 0.7f, 0.7f, 1.0f};
    GLfloat directedLightPos[] = {-10.0f, 15.0f, 20.0f, 0.0f};
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  directedLight);
    glLightfv(GL_LIGHT0, GL_POSITION, directedLightPos);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, _textureId_dam3);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

    for (int j = 0; j < TMAXY; j++)
        for (int i = 0; i < TMAXX; i++)
            draw_rect(i, j);

    glDisable(GL_TEXTURE_2D);
}

void draw_rect(int i, int j) {
    float fpb = (float)FPB;
    float vv[3];
    glColor3f(ccr[i][j], ccg[i][j], ccb[i][j]);
    glBegin(GL_POLYGON);
        vv[0] = (float)i     * fpb;  vv[1] = (float)j     * fpb;  vv[2] = gnd[i][j];
        glTexCoord2f(1.0f, 0.0f);  glVertex3fv(vv);
        vv[0] = (float)(i+1) * fpb;                                vv[2] = gnd[i+1][j];
        glTexCoord2f(1.0f, 1.0f);  glVertex3fv(vv);
        vv[1] = (float)(j+1) * fpb;                                vv[2] = gnd[i+1][j+1];
        glTexCoord2f(0.0f, 1.0f);  glVertex3fv(vv);
        vv[0] = (float)i     * fpb;                                vv[2] = gnd[i][j+1];
        glTexCoord2f(0.0f, 0.0f);  glVertex3fv(vv);
    glEnd();
}

// ============================================================
// DRAW AIRCRAFT
// ============================================================

void draw_new_plane(void) {
    glPushMatrix();

    glTranslatef((float)state[NORTH], (float)state[EAST], (float)state[DOWN]);
    glRotatef((float)(state[PSI]   * 180.0 / PI), 0.0f, 0.0f, 1.0f);
    glRotatef((float)(state[THETA] * 180.0 / PI), 0.0f, 1.0f, 0.0f);
    glRotatef((float)(state[PHI]   * 180.0 / PI), 1.0f, 0.0f, 0.0f);

    glTranslatef((float)cg_loc_x, 0.0f, (float)cg_loc_z);

    // orientation fix for Extra300170.stl (found by trial & error)
    glRotatef(180.0f, 0.0f, 1.0f, 0.0f);
    glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
    glRotatef(-90.0f, 0.0f, 1.0f, 0.0f);
    glRotatef(180.0f, 0.0f, 0.0f, 1.0f);
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);  // Fix 90° rotation around Y-axis

    // STL is in mm, sim is in feet  (~1 mm = 0.00328 ft, scale 0.05 gives smaller size)
    glScalef(0.05f, 0.05f, 0.05f);

    if (wireframe_mode == 2) {
        // Wireframe only mode
        glDisable(GL_LIGHTING);
        glColor3f(0.9f, 0.5f, 0.2f);  // Orange edges
        glLineWidth(1.5f);
        if (uav_model.tri_count > 0)
            draw_stl_wireframe(uav_model);
        glLineWidth(1.0f);
    } else if (wireframe_mode == 1) {
        // Solid wireframe mode (filled + edges)
        // First draw filled model
        glColor3f(0.35f, 0.65f, 0.85f);  // Nice light blue
        GLfloat specular[] = { 0.8f, 0.8f, 0.8f, 1.0f };
        GLfloat shininess[] = { 32.0f };
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
        
        if (uav_model.tri_count > 0)
            draw_stl(uav_model);
        else
            glCallList(display_list);
        
        // Then overlay wireframe edges
        glDisable(GL_LIGHTING);
        glColor3f(0.1f, 0.1f, 0.1f);  // Dark gray edges
        glLineWidth(0.8f);
        glPolygonOffset(1.0f, 1.0f);
        glEnable(GL_POLYGON_OFFSET_LINE);
        if (uav_model.tri_count > 0)
            draw_stl_wireframe(uav_model);
        glDisable(GL_POLYGON_OFFSET_LINE);
        glLineWidth(1.0f);
        glEnable(GL_LIGHTING);
    } else {
        // Mode 0: Solid with subtle edges (default)
        glColor3f(0.35f, 0.65f, 0.85f);  // Nice light blue
        GLfloat specular[] = { 0.8f, 0.8f, 0.8f, 1.0f };
        GLfloat shininess[] = { 32.0f };
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
        
        if (uav_model.tri_count > 0)
            draw_stl(uav_model);
        else
            glCallList(display_list);
    }

    glPopMatrix();
}

// ============================================================
// HUD  (2D overlay — drawn after glPopMatrix so it stays on screen)
// ============================================================

void draw_cockpit(void) {
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0, XMAXSCREEN, 0.0, YMAXSCREEN);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);

    if (text_display)
        draw_text_readout();

    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

void draw_text_readout(void) {
    char buffer[100];
    double current_altitude = -state[DOWN];
    
    glColor3f(0.0f, 0.0f, 0.0f);

    sprintf(buffer, "SPEED:%.1f fps  ALT:%.0f ft  V/S:%.0f fpm  ",
            airspeed_read,
            current_altitude,
            airspeed_read * sin(state[THETA] - alpha_read) * 60.0);
    glRasterPos2f(100.0f, 100.0f);
    print_string1(buffer, H18);

    sprintf(buffer, "PITCH:%.2f deg ROLL:%.2f deg HEADING:%.2f deg  ",
            state[THETA]*180.0/PI, state[PHI]*180.0/PI, state[PSI]*180.0/PI);
    glRasterPos2f(100.0f, 80.0f);
    print_string1(buffer, H18);

    sprintf(buffer, "ALPHA:%.1f deg  BETA:%.1f deg  ",
            alpha_read*180.0/PI, beta_read*180.0/PI);
    glRasterPos2f(100.0f, 60.0f);
    print_string2(buffer, H18);

    sprintf(buffer, "Ail:%.2f deg Elev:%.2f deg  Rudd:%.2f deg     ",
            control[AILERON], control[ELEVATOR_2], control[RUDDER]);
    glRasterPos2f(100.0f, 40.0f);
    print_string1(buffer, H18);

    sprintf(buffer, "throttle:%.0f   ", control[THROTTLE_L]);
    glRasterPos2f(100.0f, 20.0f);
    print_string4(buffer, H18);

    // ========== BOMBING STATUS ==========
    if (current_altitude < LOW_ALTITUDE_THRESHOLD && !bomb_dropped) {
        // In bombing zone - show red warning
        glColor3f(1.0f, 0.0f, 0.0f);
        sprintf(buffer, "*** LOW ALTITUDE - PREPARE TO DROP BOMB ***");
        glRasterPos2f(100.0f, 120.0f);
        print_string1(buffer, H18);
    } else if (bomb_dropped && current_altitude < ESCAPE_ALTITUDE) {
        // Bomb dropped, climbing out - show status in yellow
        glColor3f(1.0f, 1.0f, 0.0f);
        sprintf(buffer, "BOMB DROPPED! CLIMB TO %.0f FT FOR DETONATION (%.0f ft now)", 
                ESCAPE_ALTITUDE, current_altitude);
        glRasterPos2f(100.0f, 120.0f);
        print_string1(buffer, H18);
    }
}

void print_string1(char *s, void *fontname) {
    for (int i = 0; i < 47; i++) glutBitmapCharacter(fontname, s[i]);
}
void print_string2(char *s, void *fontname) {
    for (int i = 0; i < 27; i++) glutBitmapCharacter(fontname, s[i]);
}
void print_string3(char *s, void *fontname) {
    for (int i = 0; i < 27; i++) glutBitmapCharacter(fontname, s[i]);
}
void print_string4(char *s, void *fontname) {
    for (int i = 0; i < 13; i++) glutBitmapCharacter(fontname, s[i]);
}

// ============================================================
// READ STATE FROM SIMULINK FILE
// ============================================================

void read_state(void) {
    double num1,  num2,  num3,  num4,
           num5,  num6,  num7,  num8,
           num9,  num10, num11, num12,
           num13, num14, num15, num16;

    indata >> num1  >> num2  >> num3  >> num4
           >> num5  >> num6  >> num7  >> num8
           >> num9  >> num10 >> num11 >> num12
           >> num13 >> num14 >> num15 >> num16;

    if (indata.eof()) {
        indata.close();
        cout << "End-of-file reached.." << endl;
        indata.open("simRecording_20260309_2221.txt");
        if (!indata) {
            cerr << "Error: file could not be opened" << endl;
            exit(0);
        }
    }

    // num1 = time column (Simulink), skip it
    state[PHI]   = num2;
    state[THETA] = num3;
    state[PSI]   = num4;

    state[NORTH] = num5  + 240000.0;
    state[EAST]  = num6  + 250000.0;
    state[DOWN]  = -num7;

    TARGET_N = num8  + 240000.0;
    TARGET_E = num9  + 250000.0;
    TARGET_D = state[DOWN];   // same altitude as aircraft for chase

    control[AILERON]    = num10 * 180.0 / PI;
    control[RUDDER]     = num11 * 180.0 / PI;
    control[ELEVATOR_2] = num12 * 180.0 / PI;
    control[THROTTLE_L] = num13;
    control[THROTTLE_R] = num13;

    airspeed_read = num14;
    alpha_read    = num15;
    beta_read     = num16;

    // ========== BOMBING RUN LOGIC ==========
    double current_altitude = -state[DOWN];  // altitude is negative DOWN in NED coordinates

    // Check if we're at low altitude (below 400ft)
    if (current_altitude < LOW_ALTITUDE_THRESHOLD && !bomb_dropped) {
        was_at_low_altitude = 1;
        bomb_dropped = 1;  // Mark bomb as dropped
        // Store the location where bomb was dropped (plane's current position)
        bomb_location_n = state[NORTH];
        bomb_location_e = state[EAST];
        bomb_location_d = state[DOWN];
        printf("*** BOMB DROPPED at altitude %.0f ft! Location: N:%.0f E:%.0f ***\n", 
               current_altitude, bomb_location_n, bomb_location_e);
    }

    // Check if we've climbed back out (above 600ft) after dropping bomb
    if (bomb_dropped && current_altitude > ESCAPE_ALTITUDE && !explosion_active) {
        // Trigger explosion at the location where bomb was dropped
        explosion_active = 1;
        explosion_start_time = 0.0;
        explosion_x = bomb_location_n;
        explosion_y = bomb_location_e;
        explosion_z = bomb_location_d;
        bomb_dropped = 0;  // Reset for next run
        printf("*** BOMB HIT! Explosion at drop location! Altitude now: %.0f ft ***\n", current_altitude);
    }

    // Update explosion timer
    if (explosion_active) {
        explosion_time += 0.016;  // ~60 FPS
        if (explosion_time > EXPLOSION_DURATION) {
            explosion_active = 0;
            explosion_time = 0.0;
        }
    }
}
