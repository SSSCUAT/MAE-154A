/*                             flight.h
                          Misc. Definitions 
             Flight Simulation 16.431  Spring, 2003  J. Kuchar
                   Massachusetts Institute of Technology
*/

#define PI         3.14159
#define G          32.2              /* gravitational accel   (ft/sec^2)         */
#define RHO        0.002378          /* air density           (slug/ft^3)        */
#define DEG_TO_FT  364800.0          /* conversion from degrees latitude to feet */

/* State variables                                   */
/* States are placed in a vector of length NUMSTATES */
/*      state     vector index                       */
/*      -------   ------------                       */
#define NORTH     0
#define EAST      1
#define DOWN      2
#define U         3
#define V         4
#define W         5
#define PSI       6
#define THETA     7
#define PHI       8
#define P         9
#define Q         10
#define R         11
#define E0        12
#define E1        13
#define E2        14
#define E3        15
#define NUMSTATES 16   /* number of elements in state vector          */
/* if you add or remove states, be sure to change NUMSTATES as needed */

/* define the state vector */
double state[NUMSTATES];

/* Control variables */
#define ELEVATOR    0
#define AILERON     1
#define RUDDER      2
#define THROTTLE_L  3
#define THROTTLE_R  4
#define ELEVATOR_2  5
#define NUMCONTROLS 6   /* number of elements in control vector           */
/* if you add or remove controls, be sure to change NUMCONTROLS as needed */

/* define the control vector */
double control[NUMCONTROLS];

#define K 0.1              /* constant for quaternion correction */

#define MOUSERADIUS          100.0  /* range of allowable mouse motion (in pixels) */ 
#define JOYSTICK_SENSITIVITY 0.001  /* gain of joystick control */

clock_t lasttime;          /* stores time at last loop, to calc dt per loop */

//#define TRUE  1
//#define FALSE 0

