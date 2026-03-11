/*                             vehicle.h
                      Vehicle Parameter Definitions 
             Flight Simulation 16.431  Spring, 2003  J. Kuchar
                   Massachusetts Institute of Technology
*/

/*------------------------*/
/* KAWASAKI C-1 TRANSPORT */
/* (estimated parameters) */
/*------------------------*/

#define WEIGHT          99210.0                  /* gross weight            (lb)    */
#define MASS            (WEIGHT/G)               /* gross mass              (slug)  */

#define WING_SPAN       100.4                    /* wingspan                (ft)    */
#define WING_CHORD      13.3                     /* mean wing chord         (ft)    */
#define WING_AREA       (WING_SPAN*WING_CHORD)   /* wing planform area      (ft^2)  */
#define WING_ASPECT     (WING_SPAN/WING_CHORD)   /* wing aspect ratio       ()      */
#define WING_XCP        -4.2                     /* x-axis wing CP          (ft)    */
#define WING_YCP        18.7                     /* y-axis wing CP          (ft)    */
#define WING_ZCP        -10.0                    /* z-axis wing CP          (ft)    */
#define WING_SWEEP      (15.0*PI/180.0)          /* wing sweep angle        (rad)   */

#define HTAIL_SPAN      37.3                     /* horiz tail span         (ft)    */
#define HTAIL_CHORD     7.47                     /* horiz tail chord        (ft)    */
#define HTAIL_AREA      (HTAIL_SPAN*HTAIL_CHORD) /* horiz tail area         (ft^2)  */
#define HTAIL_ASPECT    (HTAIL_SPAN/HTAIL_CHORD) /* horiz tail aspect ratio ()      */
#define HTAIL_XCP       -51.4                    /* x-axis h-tail CP        (ft)    */
#define HTAIL_ZCP       -25.7                    /* z-axis h-tail CP        (ft)    */
#define HTAIL_INCIDENCE -4.7*PI/180.0            /* h-tail incidence angle  (rad)   */

#define VTAIL_SPAN      14.1                     /* vert tail span          (ft)    */
#define VTAIL_CHORD     14.1                     /* vert tail chord         (ft)    */
#define VTAIL_AREA      (VTAIL_SPAN*VTAIL_CHORD) /* vert tail area          (ft^2)  */
#define VTAIL_ASPECT    (VTAIL_SPAN/VTAIL_CHORD) /* vert tail aspect ratio  ()      */
#define VTAIL_XCP       -51.4                    /* x-axis v-tail CP        (ft)    */
#define VTAIL_ZCP       -18.3                    /* z-axis v-tail CP        (ft)    */

#define THRUST_MAX      14500.0                  /* max thrust per engine   (lb)    */

#define CL_ALPHA        (2*PI)                   /* base lift curve slope   (1/rad) */
#define K_A             0.2                      /* aileron gain            ()      */
#define K_E             0.8                      /* elevator gain           ()      */
#define K_R             1.2                      /* rudder gain             ()      */
#define CD0             0.015                    /* base coeff. of drag     ()      */
#define FUSE_SCD        50.0                     /* fuselage S*Cd           (ft^2)  */

/* Moments of Inertia                                                 */
/* (scaled-down from published values for the C5A Galaxy)             */
#define IX       255000.0 //2550000.0       /* x-axis                (slug-ft^2) */
#define IY       417000.0//4170000.0       /* y-axis                (slug-ft^2) */
#define IZ       627000.0//6270000.0       /* z-axis                (slug-ft^2) */
#define IXZ      33000.0//330000.0        /* cross-moment          (slug-ft^2) */

/* moment of inertia constants for computing pdot, qdot, rdot */
/* (using the notation in Stevens and Lewis)                  */
#define IGAMMA (IX*IZ-IXZ*IXZ)
#define C1     ((IY*IZ-IZ*IZ-IXZ*IXZ)/IGAMMA)
#define C2     ((IX-IY+IZ)*IXZ/IGAMMA)
#define C3     (IZ/IGAMMA)
#define C4     (IXZ/IGAMMA)
#define C5     ((IZ-IX)/IY)
#define C6     (IXZ/IY)
#define C7     (1.0/IY)
#define C8     ((IX*IX-IX*IY+IXZ*IXZ)/IGAMMA)
#define C9     (IX/IGAMMA)

