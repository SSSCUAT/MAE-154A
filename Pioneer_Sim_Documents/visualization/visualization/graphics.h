/*                             graphics.h
                         Graphical Definitions 
             Flight Simulation 16.431  Spring, 2003  J. Kuchar
                   Massachusetts Institute of Technology
*/

/* Perspective projection field of view (vertical, degrees) */
#define FOV 40.0

/* screen parameters */
double ASPECT_RATIO, XCTR, YCTR, XMAXSCREEN, YMAXSCREEN;

/* Fog colors */
#define FOGR         0.49
#define FOGG         0.65
#define FOGB         0.79
//#define FOGR         0.05
//#define FOGG         0.05
//#define FOGB         0.1
#define FOG_DENSITY  0.00002
#define FOG_START    100.0     /* no fog before this distance           */
#define FOG_END      400000.0   /* completely fogged after this distance */

/* ground colors */
#define GNDR 0.90
#define GNDG 0.90
#define GNDB 0.90

/* ground colors */
//#define GNDRHor 0.90
//#define GNDGHor 0.90
//#define GNDBHor 1.0

/* ground colors */
#define GNDR 0.72
#define GNDG 0.94
#define GNDB 0.41

/* water colors */
#define WATR 0.41
#define WATG 0.72
#define WATB 0.94

/* The following are used to define the terrain in the environment       */
/* gnd[][] holds terrain altitude, ccr[][] holds red terrain color, etc. */
float gnd[500][1000], ccr[500][500], ccg[500][500], ccb[500][500];

/* terrain is divided into an array of blocks */
/* Max # of ground blocks in x & y directions */
#define TMAXX 50
#define TMAXY 50

/* Size of each block in terms of feet per block */
#define FPB 10000.0

/* Ground altitude amplitude and frequencies                 */
/* (altitude calculated using sum of sines to appear random) */
#define XAMP -0000.0
#define XFREQ 0.5
#define YFREQ 0.6

/* control tower position */
#define TOWER_N (TMAXX*FPB/2.0+700.0)
#define TOWER_E (TMAXY*FPB/2.0+700.0)
#define TOWER_D (-200.0)

/* view settings */
#define OUT_THE_WINDOW 0
#define CHASEPLANE 1
#define TOWER 2

/* variables used to control the view */
double pan, zoom, sway;
int view, text_display;

/* the following are used to keep track of another aircraft */
/* in the environment that flies in circles                 */
double oeast, onorth, odown, opsi, otheta, ophi;
float ang, val;

/* Other aircraft flight radius, origin */
#define ORAD           1000.0
#define OALT_RANGE     200.0
#define ONORTH_ORIGIN  (TMAXX*FPB/2.0+3000.0)
#define OEAST_ORIGIN   (TMAXY*FPB/2.0)
#define OALT_ORIGIN    -750.0

/* font names */
#define H10       GLUT_BITMAP_HELVETICA_10
#define H12       GLUT_BITMAP_HELVETICA_12
#define H18       GLUT_BITMAP_HELVETICA_18
#define C15       GLUT_BITMAP_9_BY_15      /* Courier */
#define C13       GLUT_BITMAP_8_BY_13




