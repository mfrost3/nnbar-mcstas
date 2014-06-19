/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: /home/matt/Desktop/WORK/full_1b_moderator/source_test.instr (test)
 * Date:       Wed Jun 18 16:30:00 2014
 * File:       /home/matt/Desktop/WORK/full_1b_moderator/source_test.c
 */


#define MCCODE_STRING "McStas 2.0 - Dec. 21, 2012"
#define FLAVOR "mcstas"
#define FLAVOR_UPPER "MCSTAS"
#define MC_USE_DEFAULT_MAIN
#define MC_TRACE_ENABLED
#define MC_EMBEDDED_RUNTIME

#line 1 "mccode-r.h"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas 2.0
* Version: $Revision: 1.112 $
*
* Runtime system header for McStas.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas/McXtrace version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id: mcstas-r.h,v 1.112 2009/08/13 14:52:01 farhi Exp $
*
*******************************************************************************/

#ifndef MCCODE_R_H
#define MCCODE_R_H "$Revision: 1.0 $"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <float.h>

/* If the runtime is embedded in the simulation program, some definitions can
   be made static. */

#ifdef MC_EMBEDDED_RUNTIME
#define mcstatic static
#else
#define mcstatic
#endif

#ifdef __dest_os
#if (__dest_os == __mac_os)
#define MAC
#endif
#endif

#ifdef __FreeBSD__
#define NEED_STAT_H
#endif

#if defined(__APPLE__) && defined(__GNUC__)
#define NEED_STAT_H
#endif

#ifdef NEED_STAT_H
#include <sys/stat.h>
#endif

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#define random rand
#define srandom srand
#else  /* !WIN32 */
#ifdef MAC
#define MC_PATHSEP_C ':'
#define MC_PATHSEP_S ":"
#else  /* !MAC */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* MC_PATHSEP_C */

/* the version string is replaced when building distribution with mkdist */
#ifndef MCCODE_STRING
#define MCCODE_STRING "McStas 2.0 - Dec. 21, 2012"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "Dec. 21, 2012"
#endif

#ifndef MCCODE_VERSION
#define MCCODE_VERSION "2.0"
#endif

#ifndef MCCODE_NAME
#define MCCODE_NAME "McStas"
#endif

#ifndef MCCODE_PARTICULE
#define MCCODE_PARTICULE "neutron"
#endif

#ifndef MCCODE_LIBENV
#define MCCODE_LIBENV "MCSTAS"
#endif

#ifndef FLAVOR_UPPER
#define FLAVOR_UPPER MCCODE_NAME
#endif

#ifdef MC_PORTABLE
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#ifdef MAC
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (USE_MPI == 0)
#undef USE_MPI
#endif

#ifdef USE_MPI  /* default is to disable signals with MPI, as MPICH uses them to communicate */
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (NOSIGNALS == 0)
#undef NOSIGNALS
#endif

#if (USE_NEXUS == 0)
#undef USE_NEXUS
#endif

/* I/O section part ========================================================= */

/* Note: the enum instr_formal_types definition MUST be kept
   synchronized with the one in mcstas.h and with the
   instr_formal_type_names array in cogen.c. */
enum instr_formal_types
  {
    instr_type_double, instr_type_int, instr_type_string
  };
struct mcinputtable_struct { /* defines instrument parameters */
  char *name; /* name of parameter */
  void *par;  /* pointer to instrument parameter (variable) */
  enum instr_formal_types type;
  char *val;  /* default value */
};

typedef double MCNUM;
typedef struct {MCNUM x, y, z;} Coords;
typedef MCNUM Rotation[3][3];

/* the following variables are defined in the McStas generated C code
   but should be defined externally in case of independent library usage */
#ifndef DANSE
extern struct mcinputtable_struct mcinputtable[]; /* list of instrument parameters */
extern int    mcnumipar;                          /* number of instrument parameters */
extern char   mcinstrument_name[], mcinstrument_source[]; /* instrument name and filename */
extern MCNUM  mccomp_storein[]; /* 11 coords * number of components in instrument */
extern MCNUM  mcAbsorbProp[];
extern MCNUM  mcScattered;      /* number of SCATTER calls in current component */
#ifndef MC_ANCIENT_COMPATIBILITY
extern int mctraceenabled, mcdefaultmain;
#endif
#endif

/* file I/O definitions and function prototypes */

struct mcformats_struct {
  char *Name;       /* format name: may also specify: append, partial(hidden), binary */
  char *Extension;  /* file name extension set if not specified */
  char *Header;     /* to be written as header(begin) in output file */
  char *Footer;     /* to be written as foorter(end) in output file */
  char *BeginSection; /* defines how to start a new section */
  char *EndSection;   /* defines how to end a new section */
  char *AssignTag;    /* defines how to write e.g. 'par=value' */
  char *BeginData;    /* defines how to start a data block */
  char *EndData;      /* defines how to end a data block */
  char *BeginErrors;  /* defines how to start an error(statistical) block */
  char *EndErrors;    /* defines how to end an error(statistical) block */
  char *BeginNcount;  /* defines how to start an event (number of) block */
  char *EndNcount;    /* defines how to end an event (number of) block */
  };

#ifndef MC_EMBEDDED_RUNTIME /* the mcstatic variables (from mcstas-r.c) */
extern FILE * mcsiminfo_file;     /* handle to the output siminfo file */
extern int    mcgravitation;      /* flag to enable gravitation */
extern int    mcdotrace;          /* flag to print MCDISPLAY messages */
extern struct mcformats_struct mcformats[]; /* array of all possible formats */
extern struct mcformats_struct mcformat;    /* selected format to be used */
extern struct mcformats_struct mcformat_data; /* same as above, but may be different for data blocks */
#else
mcstatic FILE *mcsiminfo_file        = NULL;
#endif

/* Useful macros ============================================================ */

/* MPI stuff */

#ifdef USE_MPI
#include "mpi.h"

#ifdef OMPI_MPI_H  /* openmpi does not use signals: we may install our sighandler */
#undef NOSIGNALS
#endif

/*
 * MPI_MASTER(i):
 * execution of i only on master node
 */
#define MPI_MASTER(statement) { \
  if(mpi_node_rank == mpi_node_root)\
  { statement; } \
}

#ifndef MPI_REDUCE_BLOCKSIZE
#define MPI_REDUCE_BLOCKSIZE 1000
#endif

int mc_MPI_Sum(double* buf, long count);
int mc_MPI_Send(void *sbuf, long count, MPI_Datatype dtype, int dest);
int mc_MPI_Recv(void *rbuf, long count, MPI_Datatype dtype, int source);

/* MPI_Finalize exits gracefully and should be preferred to MPI_Abort */
#define exit(code) do {                                   \
    MPI_Finalize();                                       \
    exit(code);                                           \
  } while(0)

#else /* !USE_MPI */
#define MPI_MASTER(instr) instr
#endif /* USE_MPI */

#ifdef USE_MPI
static int mpi_node_count;
#endif

#ifdef USE_THREADS  /* user want threads */
#error Threading (USE_THREADS) support has been removed for very poor efficiency. Use MPI/SSH grid instead.
#endif

/* I/O function prototypes ================================================== */

#ifndef CHAR_BUF_LENGTH
#define CHAR_BUF_LENGTH 1024
#endif

/* main DETECTOR structure which stores most information to write to data files */
struct mcdetector_struct {
  char   prefix[CHAR_BUF_LENGTH];     /* prefix to prepend to format specifications */
  char   filename[CHAR_BUF_LENGTH];   /* file name of monitor */
  char   position[CHAR_BUF_LENGTH];   /* position of detector component */
  char   component[CHAR_BUF_LENGTH];  /* component instance name */
  char   instrument[CHAR_BUF_LENGTH]; /* instrument name */
  char   simulation[CHAR_BUF_LENGTH]; /* simulation name, e.g. directory */
  char   type[CHAR_BUF_LENGTH];       /* data type, e.g. 0d, 1d, 2d, 3d */
  char   user[CHAR_BUF_LENGTH];       /* user name, e.g. HOME */
  char   date[CHAR_BUF_LENGTH];       /* date of simulation end/write time */
  char   title[CHAR_BUF_LENGTH];      /* title of detector */
  char   xlabel[CHAR_BUF_LENGTH];     /* X axis label */
  char   ylabel[CHAR_BUF_LENGTH];     /* Y axis label */
  char   zlabel[CHAR_BUF_LENGTH];     /* Z axis label */
  char   xvar[CHAR_BUF_LENGTH];       /* X variable name */
  char   yvar[CHAR_BUF_LENGTH];       /* Y variable name */
  char   zvar[CHAR_BUF_LENGTH];       /* Z variable name */
  char   ncount[CHAR_BUF_LENGTH];     /* number of events initially generated */
  char   limits[CHAR_BUF_LENGTH];     /* X Y Z limits, e.g. [xmin xmax ymin ymax zmin zmax] */
  char   variables[CHAR_BUF_LENGTH];  /* variables written into data block */
  char   statistics[CHAR_BUF_LENGTH]; /* center, mean and half width along axis */
  char   signal[CHAR_BUF_LENGTH];     /* min max and mean of signal (data block) */
  char   values[CHAR_BUF_LENGTH];     /* integrated values e.g. [I I_err N] */
  double xmin,xmax;                   /* min max of axes */
  double ymin,ymax;
  double zmin,zmax;
  double intensity;                   /* integrated values for data block */
  double error;
  double events;
  double min;                         /* statistics for data block */
  double max;
  double mean;
  double centerX;                     /* statistics for axes */
  double halfwidthX;
  double centerY;
  double halfwidthY;
  int    rank;                        /* dimensionaly of monitor, e.g. 0 1 2 3 */
  char   istransposed;                /* flag to transpose matrix for some formats */

  long   m,n,p;                       /* dimensions of data block and along axes */
  long   date_l;                      /* same as date, but in sec since 1970 */

  double *p0, *p1, *p2;               /* pointers to saved data, NULL when freed */
  double *p0_orig,*p1_orig,*p2_orig;  /* initial pointer in case some post processing is done
                                         with e.g. lists and MPI handling */

  FILE   *file_handle;                /* file handle for detector file */
  struct mcformats_struct format;     /* format used for that detector */
};

typedef struct mcdetector_struct MCDETECTOR;

void   mcset_ncount(unsigned long long count);    /* wrapper to get mcncount */
unsigned long long int mcget_ncount(void);            /* wrapper to set mcncount */
unsigned long long mcget_run_num(void);           /* wrapper to get mcrun_num=0:mcncount */
char *mcfull_file(char *name, char *ext); /* builds "name+ext" */

/* Needed function signatures */
MCDETECTOR mcdetector_write_data(MCDETECTOR detector);

/* output functions */
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2, char *c, Coords pos);
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
                  char *xvar, double x1, double x2, int n,
                  double *p0, double *p1, double *p2, char *f, char *c, Coords pos);
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2, int m,
                  int n, double *p0, double *p1, double *p2, char *f,
                  char *c, Coords pos);
MCDETECTOR mcdetector_out_3D(char *t, char *xl, char *yl, char *zl,
      char *xvar, char *yvar, char *zvar,
                  double x1, double x2, double y1, double y2, double z1, double z2, int m,
                  int n, int p, double *p0, double *p1, double *p2, char *f,
                  char *c, Coords pos);
/* wrappers to output functions, that automatically set NAME and POSITION */
#define DETECTOR_OUT(p0,p1,p2) mcdetector_out_0D(NAME_CURRENT_COMP,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_0D(t,p0,p1,p2) mcdetector_out_0D(t,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f) \
     mcdetector_out_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f) \
     mcdetector_out_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_3D(t,xl,yl,zl,xv,yv,zv,x1,x2,y1,y2,z1,z2,m,n,p,p0,p1,p2,f) \
     mcdetector_out_3D(t,xl,yl,zl,xv,yv,zv,x1,x2,y1,y2,z1,z2,m,n,p,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_CUSTOM_HEADER(t)  if (t && strlen(t)) { \
     mcDetectorCustomHeader=malloc(strlen(t)); \
     if (mcDetectorCustomHeader) strcpy(mcDetectorCustomHeader, t); }

/* Following part is only embedded when not redundent with mcstas.h ========= */

#ifndef MCCODE_H

#ifndef NOSIGNALS
#include <signal.h>
#define SIG_MESSAGE(msg) strcpy(mcsig_message, msg);
#else
#define SIG_MESSAGE(msg)
#endif /* !NOSIGNALS */



/* Useful macros and constants ============================================== */

#ifndef FLT_MAX
#define FLT_MAX         3.40282347E+38F /* max decimal value of a "float" */
#endif

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef SQR
#define SQR(x) ( (x) * (x) )
#endif
#ifndef SIGN
#define SIGN(x) (((x)>0.0)?(1):(-1))
#endif

#define RAD2MIN  ((180*60)/PI)
#define MIN2RAD  (PI/(180*60))
#define DEG2RAD  (PI/180)
#define RAD2DEG  (180/PI)
#define FWHM2RMS 0.424660900144    /* Convert between full-width-half-max and */
#define RMS2FWHM 2.35482004503     /* root-mean-square (standard deviation) */
#define HBAR     1.05457168e-34    /* [Js] h bar Planck constant CODATA 2002 */
#define MNEUTRON 1.67492728e-27    /* [kg] mass of neutron CODATA 2002 */
#define GRAVITY  9.81              /* [m/s^2] gravitational acceleration */
#define NA       6.02214179e23     /* [#atoms/g .mole] Avogadro's number*/

#ifndef PI
# ifdef M_PI
#  define PI M_PI
# else
#  define PI 3.14159265358979323846
# endif
#endif

/* wrapper to get absolute and relative position of comp */
/* mccomp_posa and mccomp_posr are defined in McStas generated C code */
#define POS_A_COMP_INDEX(index) \
    (mccomp_posa[index])
#define POS_R_COMP_INDEX(index) \
    (mccomp_posr[index])
/* number of SCATTER calls in current comp: mcScattered defined in McStas generated C code */
#define SCATTERED mcScattered

/* Retrieve component information from the kernel */
/* Name, position and orientation (both absolute and relative)  */
/* Any component: For "redundancy", see comment by KN */
#define tmp_name_comp(comp) #comp
#define NAME_COMP(comp) tmp_name_comp(comp)
#define tmp_pos_a_comp(comp) (mcposa ## comp)
#define POS_A_COMP(comp) tmp_pos_a_comp(comp)
#define tmp_pos_r_comp(comp) (mcposr ## comp)
#define POS_R_COMP(comp) tmp_pos_r_comp(comp)
#define tmp_rot_a_comp(comp) (mcrota ## comp)
#define ROT_A_COMP(comp) tmp_rot_a_comp(comp)
#define tmp_rot_r_comp(comp) (mcrotr ## comp)
#define ROT_R_COMP(comp) tmp_rot_r_comp(comp)

/* Current component name, index, position and orientation */
#define NAME_CURRENT_COMP  NAME_COMP(mccompcurname)
#define INDEX_CURRENT_COMP mccompcurindex
#define POS_A_CURRENT_COMP POS_A_COMP(mccompcurname)
#define POS_R_CURRENT_COMP POS_R_COMP(mccompcurname)
#define ROT_A_CURRENT_COMP ROT_A_COMP(mccompcurname)
#define ROT_R_CURRENT_COMP ROT_R_COMP(mccompcurname)

/* Note: The two-stage approach to MC_GETPAR is NOT redundant; without it,
* after #define C sample, MC_GETPAR(C,x) would refer to component C, not to
* component sample. Such are the joys of ANSI C.

* Anyway the usage of MCGETPAR requires that we use sometimes bare names...
*/
#define MC_GETPAR2(comp, par) (mcc ## comp ## _ ## par)
#define MC_GETPAR(comp, par) MC_GETPAR2(comp,par)

/* MCDISPLAY/trace and debugging message sent to stdout */
#ifdef MC_TRACE_ENABLED
#define DEBUG
#endif

#ifdef DEBUG
#define mcDEBUG_INSTR() if(!mcdotrace); else { printf("INSTRUMENT:\n"); printf("Instrument '%s' (%s)\n", mcinstrument_name, mcinstrument_source); }
#define mcDEBUG_COMPONENT(name,c,t) if(!mcdotrace); else {\
  printf("COMPONENT: \"%s\"\n" \
         "POS: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         name, c.x, c.y, c.z, t[0][0], t[0][1], t[0][2], \
         t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]); \
  printf("Component %30s AT (%g,%g,%g)\n", name, c.x, c.y, c.z); \
  }
#define mcDEBUG_INSTR_END() if(!mcdotrace); else printf("INSTRUMENT END:\n");
#define mcDEBUG_ENTER() if(!mcdotrace); else printf("ENTER:\n");
#define mcDEBUG_COMP(c) if(!mcdotrace); else printf("COMP: \"%s\"\n", c);
#define mcDEBUG_LEAVE() if(!mcdotrace); else printf("LEAVE:\n");
#define mcDEBUG_ABSORB() if(!mcdotrace); else printf("ABSORB:\n");
#else
#define mcDEBUG_INSTR()
#define mcDEBUG_COMPONENT(name,c,t)
#define mcDEBUG_INSTR_END()
#define mcDEBUG_ENTER()
#define mcDEBUG_COMP(c)
#define mcDEBUG_LEAVE()
#define mcDEBUG_ABSORB()
#endif

// mcDEBUG_STATE and mcDEBUG_SCATTER are defined by mcstas-r.h and mcxtrace-r.h



#ifdef TEST
#define test_printf printf
#else
#define test_printf while(0) printf
#endif

/* send MCDISPLAY message to stdout to show gemoetry */
void mcdis_magnify(char *what);
void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2);
void mcdis_dashed_linemcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n);
void mcdis_multiline(int count, ...);
void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height);
void mcdis_box(double x, double y, double z,
	       double width, double height, double length);
void mcdis_circle(char *plane, double x, double y, double z, double r);

/* selection of random number generator. default is MT */
#ifndef MC_RAND_ALG
#define MC_RAND_ALG 1
#endif

#if MC_RAND_ALG == 0
   /* Use system random() (not recommended). */
#  define MC_RAND_MAX RAND_MAX
#elif MC_RAND_ALG == 1
   /* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
#  define MC_RAND_MAX ((unsigned long)0xffffffff)
#  define random mt_random
#  define srandom mt_srandom
#elif MC_RAND_ALG == 2
   /* Algorithm used in McStas CVS-080208 and earlier (not recommended). */
#  define MC_RAND_MAX 0x7fffffff
#  define random mc_random
#  define srandom mc_srandom
#else
#  error "Bad value for random number generator choice."
#endif

typedef int mc_int32_t;
mc_int32_t mc_random(void);
void mc_srandom (unsigned int x);
unsigned long mt_random(void);
void mt_srandom (unsigned long x);

double rand01();
double randpm1();
double rand0max(double max);
double randminmax(double min, double max);

double randnorm(void);
double randtriangle(void);

#ifndef DANSE
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);
#endif

/* simple vector algebra ==================================================== */
#define vec_prod(x, y, z, x1, y1, z1, x2, y2, z2) \
	vec_prod_func(&x, &y, &z, x1, y1, z1, x2, y2, z2)
mcstatic inline void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1, double x2, double y2, double z2);

mcstatic inline double scalar_prod(
		double x1, double y1, double z1, double x2, double y2, double z2);

#define NORM(x,y,z) \
	norm_func(&x, &y, &z)
mcstatic inline void norm_func(double *x, double *y, double *z) {
	double temp = (*x * *x) + (*y * *y) + (*z * *z);
	if (temp != 0) {
		temp = sqrt(temp);
		*x /= temp;
		*y /= temp;
		*z /= temp;
	}
}

void normal_vec(double *nx, double *ny, double *nz,
    double x, double y, double z);

/**
 * Rotate the given xyz in a certain way. Used in pol-lib
 */
#define rotate(x, y, z, vx, vy, vz, phi, ax, ay, az) \
  do { \
    double mcrt_tmpx = (ax), mcrt_tmpy = (ay), mcrt_tmpz = (az); \
    double mcrt_vp, mcrt_vpx, mcrt_vpy, mcrt_vpz; \
    double mcrt_vnx, mcrt_vny, mcrt_vnz, mcrt_vn1x, mcrt_vn1y, mcrt_vn1z; \
    double mcrt_bx, mcrt_by, mcrt_bz; \
    double mcrt_cos, mcrt_sin; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vp = scalar_prod((vx), (vy), (vz), mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vpx = mcrt_vp*mcrt_tmpx; \
    mcrt_vpy = mcrt_vp*mcrt_tmpy; \
    mcrt_vpz = mcrt_vp*mcrt_tmpz; \
    mcrt_vnx = (vx) - mcrt_vpx; \
    mcrt_vny = (vy) - mcrt_vpy; \
    mcrt_vnz = (vz) - mcrt_vpz; \
    vec_prod(mcrt_bx, mcrt_by, mcrt_bz, \
             mcrt_tmpx, mcrt_tmpy, mcrt_tmpz, mcrt_vnx, mcrt_vny, mcrt_vnz); \
    mcrt_cos = cos((phi)); mcrt_sin = sin((phi)); \
    mcrt_vn1x = mcrt_vnx*mcrt_cos + mcrt_bx*mcrt_sin; \
    mcrt_vn1y = mcrt_vny*mcrt_cos + mcrt_by*mcrt_sin; \
    mcrt_vn1z = mcrt_vnz*mcrt_cos + mcrt_bz*mcrt_sin; \
    (x) = mcrt_vpx + mcrt_vn1x; \
    (y) = mcrt_vpy + mcrt_vn1y; \
    (z) = mcrt_vpz + mcrt_vn1z; \
  } while(0)

/**
 * Mirror (xyz) in the plane given by the point (rx,ry,rz) and normal (nx,ny,nz)
 *
 * TODO: This define is seemingly never used...
 */
#define mirror(x,y,z,rx,ry,rz,nx,ny,nz) \
  do { \
    double mcrt_tmpx= (nx), mcrt_tmpy = (ny), mcrt_tmpz = (nz); \
    double mcrt_tmpt; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_tmpt=scalar_prod((rx),(ry),(rz),mcrt_tmpx,mcrt_tmpy,mcrt_tmpz); \
    (x) = rx -2 * mcrt_tmpt*mcrt_rmpx; \
    (y) = ry -2 * mcrt_tmpt*mcrt_rmpy; \
    (z) = rz -2 * mcrt_tmpt*mcrt_rmpz; \
  } while (0)

Coords coords_set(MCNUM x, MCNUM y, MCNUM z);
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z);
Coords coords_add(Coords a, Coords b);
Coords coords_sub(Coords a, Coords b);
Coords coords_neg(Coords a);
Coords coords_scale(Coords b, double scale);
double coords_sp(Coords a, Coords b);
Coords coords_xp(Coords b, Coords c);
void   coords_print(Coords a);
mcstatic inline void coords_norm(Coords* c);

void rot_set_rotation(Rotation t, double phx, double phy, double phz);
int  rot_test_identity(Rotation t);
void rot_mul(Rotation t1, Rotation t2, Rotation t3);
void rot_copy(Rotation dest, Rotation src);
void rot_transpose(Rotation src, Rotation dst);
Coords rot_apply(Rotation t, Coords a);

void mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
    double *vx, double *vy, double *vz, double *sx, double *sy, double *sz);
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz);

double mcestimate_error(double N, double p1, double p2);
void mcreadparams(void);

/* this is now in mcstas-r.h and mcxtrace-r.h as the number of state parameters is no longer equal*/
/* void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);
*/
void mcgenstate(void);

/* trajectory/shape intersection routines */
int inside_rectangle(double, double, double, double);
int box_intersect(double *dt_in, double *dt_out, double x, double y, double z,
    double vx, double vy, double vz, double dx, double dy, double dz);
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
    double vx, double vy, double vz, double r, double h);
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r);
/* second order equation roots */
int solve_2nd_order(double *t1, double *t2,
    double A,  double B,  double C);

/* random vector generation to shape */
void randvec_target_circle(double *xo, double *yo, double *zo,
    double *solid_angle, double xi, double yi, double zi, double radius);
#define randvec_target_sphere randvec_target_circle
void randvec_target_rect_angular(double *xo, double *yo, double *zo,
    double *solid_angle,
               double xi, double yi, double zi, double height, double width, Rotation A);
#define randvec_target_rect(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9)  randvec_target_rect_real(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,0,0,0,1)
void randvec_target_rect_real(double *xo, double *yo, double *zo,
    double *solid_angle,
	       double xi, double yi, double zi, double height, double width, Rotation A,
			 double lx, double ly, double lz, int order);

/* this is the main() */
int mccode_main(int argc, char *argv[]);


#endif /* !MCCODE_H */

#endif /* MCCODE_R_H */
/* End of file "mccode-r.h". */

#line 694 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

#line 1 "mcstas-r.h"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision: 1.112 $
*
* Runtime system header for McStas.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id: mcstas-r.h,v 1.112 2009/08/13 14:52:01 farhi Exp $
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#define MCSTAS_R_H "$Revision: 1.112 $"

/* Following part is only embedded when not redundent with mcstas.h ========= */

#ifndef MCCODE_H

#define AA2MS    629.622368        /* Convert k[1/AA] to v[m/s] */
#define MS2AA    1.58825361e-3     /* Convert v[m/s] to k[1/AA] */
#define K2V      AA2MS
#define V2K      MS2AA
#define Q2V      AA2MS
#define V2Q      MS2AA
#define SE2V     437.393377        /* Convert sqrt(E)[meV] to v[m/s] */
#define VS2E     5.22703725e-6     /* Convert (v[m/s])**2 to E[meV] */

#define SCATTER do {mcDEBUG_SCATTER(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcScattered++;} while(0)
#define ABSORB do {mcDEBUG_STATE(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcDEBUG_ABSORB(); MAGNET_OFF; goto mcabsorb;} while(0)

#define STORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcstore_neutron(mccomp_storein,index, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
#define RESTORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcrestore_neutron(mccomp_storein,index, &x, &y, &z, &vx, &vy, &vz, &t, &sx, &sy, &sz, &p);

#define MAGNET_ON \
  do { \
    mcMagnet = 1; \
  } while(0)

#define MAGNET_OFF \
  do { \
    mcMagnet = 0; \
  } while(0)

#define ALLOW_BACKPROP \
  do { \
    mcallowbackprop = 1; \
  } while(0)

#define DISALLOW_BACKPROP \
  do { \
    mcallowbackprop = 0; \
  } while(0)

#define PROP_MAGNET(dt) \
  do { \
  }while (0)
    /* change coordinates from local system to magnet system */
/*    Rotation rotLM, rotTemp; \
      Coords   posLM = coords_sub(POS_A_CURRENT_COMP, mcMagnetPos); \
      rot_transpose(ROT_A_CURRENT_COMP, rotTemp); \
      rot_mul(rotTemp, mcMagnetRot, rotLM); \
      mcMagnetPrecession(mcnlx, mcnly, mcnlz, mcnlt, mcnlvx, mcnlvy, mcnlvz, \
               &mcnlsx, &mcnlsy, &mcnlsz, dt, posLM, rotLM); \
      } while(0)
*/

#define mcPROP_DT(dt) \
  do { \
    if (mcMagnet && dt > 0) PROP_MAGNET(dt);\
    mcnlx += mcnlvx*(dt); \
    mcnly += mcnlvy*(dt); \
    mcnlz += mcnlvz*(dt); \
    mcnlt += (dt); \
    if (isnan(p) || isinf(p)) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
  } while(0)

/* ADD: E. Farhi, Aug 6th, 2001 PROP_GRAV_DT propagation with acceleration */
#define PROP_GRAV_DT(dt, Ax, Ay, Az) \
  do { \
    if(dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
    if (mcMagnet) printf("Spin precession gravity\n"); \
    mcnlx  += mcnlvx*(dt) + (Ax)*(dt)*(dt)/2; \
    mcnly  += mcnlvy*(dt) + (Ay)*(dt)*(dt)/2; \
    mcnlz  += mcnlvz*(dt) + (Az)*(dt)*(dt)/2; \
    mcnlvx += (Ax)*(dt); \
    mcnlvy += (Ay)*(dt); \
    mcnlvz += (Az)*(dt); \
    mcnlt  += (dt); \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_DT(dt) \
  do { \
    if(dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    if (mcgravitation) { Coords mcLocG; double mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    PROP_GRAV_DT(dt, mc_gx, mc_gy, mc_gz); } \
    else mcPROP_DT(dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_Z0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gz/2, -mcnlvz, -mcnlz); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Z0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_Z0 \
  do { \
    double mc_dt; \
    if(mcnlvz == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlz/mcnlvz; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlz = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_X0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gx/2, -mcnlvx, -mcnlx); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_X0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_X0 \
  do { \
    double mc_dt; \
    if(mcnlvx == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlx/mcnlvx; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlx = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_Y0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gy/2, -mcnlvy, -mcnly); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Y0; \
    DISALLOW_BACKPROP;\
  } while(0)


#define mcPROP_Y0 \
  do { \
    double mc_dt; \
    if(mcnlvy == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnly/mcnlvy; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnly = 0; \
    DISALLOW_BACKPROP; \
  } while(0)

/*moved from mccode-r.h*/
void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);

#ifdef DEBUG

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("STATE: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("SCATTER: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);

#else

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p)
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p)

#endif

#endif /* !MCCODE_H */

#endif /* MCSTAS_R_H */
/* End of file "mcstas-r.h". */

#line 926 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

#line 1 "nexus-lib.h"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2008, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/nexus-lib.h
*
* %Identification
* Written by: EF
* Date:    Jan 17, 2007
* Release: McStas CVS-080208
* Version: $Revision: 1.8 $
*
* NeXus Runtime system header for McStas.
* Overrides default mcstas runtime functions.
* Embedded within instrument in runtime mode.
*
* Usage: Automatically embbeded in the c code whenever required.
*
*******************************************************************************/

#ifdef USE_NEXUS

#include "napi.h"
#include <sys/stat.h>

/* NeXus variables to be used in functions */
NXhandle mcnxHandle;
char    *mcnxFilename=NULL;
char     mcnxversion[128];       /* init in cogen_init: 4,5 xml and compress */

/* NeXus output functions that replace calls to pfprintf in mcstas-r */
int mcnxfile_init(char *name, char *ext, char *mode, NXhandle *nxhandle);
int mcnxfile_close(NXhandle *nxHandle);

/* header/footer. f=mcsiminfo_file, datafile */
/* creates Entry=valid_parent+file+timestamp */
int mcnxinfo_header(NXhandle nxhandle, char *part,
    char *pre,                  /* %1$s  PRE  */
    char *instrname,            /* %2$s  SRC  */
    char *file,                 /* %3$s  FIL  */
    char *format_name,          /* %4$s  FMT  */
    char *date,                 /* %5$s  DAT  */
    char *user,                 /* %6$s  USR  */
    char *valid_parent,         /* %7$s  PAR = file */
    long  date_l);               /* %8$li DATL */

/* tag=value */
int mcnxinfo_tag(NXhandle nxhandle,
    char *pre,          /* %1$s PRE */
    char *valid_section,/* %2$s SEC */
    char *name,         /* %3$s NAM */
    char *value);        /* %4$s VAL */

/* output instrument/description */
int mcnxfile_instrcode(NXhandle nxhandle, 
    char *name,
    char *parent);

/* begin/end section */
int mcnxfile_section(NXhandle nxhandle, char *part,
    char *pre,          /* %1$s  PRE  */
    char *type,         /* %2$s  TYP  */
    char *name,         /* %3$s  NAM  */
    char *valid_name,   /* %4$s  VNA  */
    char *parent,       /* %5$s  PAR  */
    char *valid_parent, /* %6$s  VPA  */
    int   level);        /* %7$i  LVL */

/* data block begin/end */
int mcnxfile_data(NXhandle nxhandle, MCDETECTOR detector, char *part,
  char *valid_parent, char *valid_xlabel, char *valid_ylabel, char *valid_zlabel);

#endif
/* End of file "nexus-lib.h". */

#line 1007 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

#line 1 "nexus-lib.c"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/nexus-lib.c
*
* %Identification
* Written by: KN
* Date:    Jan 17, 2007
* Release: McStas 1.10
* Version: $Revision: 1.14 $
*
* NeXus Runtime output functions for McStas.
* Overrides default mcstas runtime functions.
* Embedded within instrument in runtime mode.
*
* Usage: Automatically embbeded in the C code whenever required.
*
*******************************************************************************/

#ifdef USE_NEXUS

/*******************************************************************************
* mcnxfile_init: Initialize NeXus file (open it). handles NeXus 4/5 and compression
* Returns: NX_ERROR or NX_OK
*******************************************************************************/
int mcnxfile_init(char *name, char *ext, char *mode, NXhandle *nxhandle)
{
  int mcnxMode=NXACC_CREATE5;
  char mcnxExt[CHAR_BUF_LENGTH];
  strcpy(mcnxExt, ext);
  char nxversion[CHAR_BUF_LENGTH];
  int i;
  if (!mcnxversion || !strlen(mcnxversion)) strcpy(nxversion, "5 zip");
  else for (i=0; i< strlen(mcnxversion) && i < 128; nxversion[i]=tolower(mcnxversion[i++]));

  if    (strstr(nxversion,"xml")) { mcnxMode =NXACC_CREATEXML; strcpy(mcnxExt, "xml"); }
  else if (strstr(nxversion,"4")) { mcnxMode =NXACC_CREATE;    strcpy(mcnxExt, "h4"); }
  else if (strstr(nxversion,"5")) { mcnxMode =NXACC_CREATE5;   strcpy(mcnxExt, "h5"); }

  if (!strcmp(mode, "a"))    mcnxMode |= NXACC_RDWR;
  mcnxFilename = mcfull_file(name, mcnxExt);
  if (NXopen(mcnxFilename, mcnxMode, nxhandle) == NX_ERROR) {
    fprintf(stderr, "Warning: NeXus: could not open file %s\n", mcnxFilename);
    mcsiminfo_file = NULL;
  } else { mcsiminfo_file=(FILE*)mcnxFilename; }
  return(mcsiminfo_file != NULL);
}

/*******************************************************************************
* mcnxfile_close: Close NeXus file
* Returns: NX_ERROR or NX_OK
*******************************************************************************/
int mcnxfile_close(NXhandle *nxHandle)
{
  return(NXclose(nxHandle));
}

/*******************************************************************************
* mcnxinfo_header: Write general header/footer information tags (header/footer)
*                  into e.g. mcsiminfo_file, datafile
*                  Group and DataSet (information) must have been opened
* Returns: NX_ERROR or NX_OK
*******************************************************************************/
int mcnxinfo_header(NXhandle nxhandle, char *part,
    char *pre,                  /* %1$s  PRE  */
    char *instrname,            /* %2$s  SRC  */
    char *file,                 /* %3$s  FIL  */
    char *format_name,          /* %4$s  FMT  */
    char *date,                 /* %5$s  DAT  */
    char *user,                 /* %6$s  USR  */
    char *valid_parent,         /* %7$s  PAR = file */
    long  date_l)               /* %8$li DATL */
{
  if (!strcmp(part, "header")) {
    if (NXputattr(nxhandle, "user_name", user, strlen(user), NX_CHAR) == NX_ERROR) {
      fprintf(stderr, "Warning: NeXus: could not write header information in /%s/%s/%s\n", 
        file, instrname, valid_parent);
      return(NX_ERROR);
    }
    char creator[CHAR_BUF_LENGTH];
    sprintf(creator, "%s gerenated with " MCCODE_STRING " (" MCCODE_DATE ") [www.mccode.org]", instrname);
    NXputattr(nxhandle, "creator", creator, strlen(creator), NX_CHAR);
    NXputattr(nxhandle, "simulation_begin", date, strlen(date), NX_CHAR);
    char *url="http://www.nexusformat.org/";
    NXputattr(nxhandle, "URL", url, strlen(url), NX_CHAR);
    char *browser="hdfview or NXbrowse or HDFExplorer";
    NXputattr(nxhandle, "Browser", browser, strlen(browser), NX_CHAR);
#if defined (USE_MPI) 
    NXputattr (nxhandle, "number_of_nodes", &mpi_node_count, 1, NX_INT32);
#endif
    return(NXputattr(nxhandle, "Format", format_name, strlen(format_name), NX_CHAR));
  } else
    return(NXputattr(nxhandle, "simulation_end", date, strlen(date), NX_CHAR));
  
} /* mcnxinfo_header */

/*******************************************************************************
* mcnxinfo_tag: write a single tag
*               Group and DataSet (information) must have been opened
* Returns: NX_ERROR or NX_OK
*******************************************************************************/
int mcnxinfo_tag(NXhandle nxhandle,
    char *pre,          /* %1$s PRE */
    char *valid_section,/* %2$s SEC */
    char *name,         /* %3$s NAM */
    char *value)        /* %4$s VAL */
{
  return(NXputattr(nxhandle, name, value, strlen(value), NX_CHAR));
} /* mcnxinfo_tag */

/*******************************************************************************
* mcnxfile_instrcode: writes the instrument description file
*                   open/close a new 'description' Data Set in the current opened Group
* Returns: NX_ERROR or NX_OK
*******************************************************************************/
int mcnxfile_instrcode(NXhandle nxhandle, 
    char *name,
    char *parent)
{
  FILE *f;
  char *instr_code=NULL;
  char nxname[CHAR_BUF_LENGTH];
  int length;
  
  struct stat stfile;
  if (stat(name,&stfile) != 0) {
    instr_code = (char*)malloc(CHAR_BUF_LENGTH);
    if (instr_code) 
      sprintf(instr_code, "File %s not found (instrument description %s is missing)", 
        name, parent);
  } else {
    long filesize = stfile.st_size;
    f=fopen(name, "r");
    instr_code = (char*)malloc(filesize);
    if (instr_code && f) fread(instr_code, 1, filesize, f);
    if (f) fclose(f);
  }
  length = instr_code ? strlen(instr_code) : 0;
  if (length) {
    time_t t;
    NXmakedata(nxhandle, "description", NX_CHAR, 1, &length);
    if (NXopendata(nxhandle, "description") == NX_ERROR) {
      fprintf(stderr, "Warning: NeXus: could not write instrument description %s (%s)\n", 
        name, parent);
      free(instr_code);
      return(NX_ERROR);
    }
    NXputdata (nxhandle, instr_code);
    free(instr_code);
    NXputattr (nxhandle, "file_name", name, strlen(name), NX_CHAR);
    NXputattr (nxhandle, "file_size", &length, 1, NX_INT32);
    t=stfile.st_mtime; strncpy(nxname, ctime(&t), CHAR_BUF_LENGTH);
    NXputattr (nxhandle, "file_date", nxname, strlen(nxname), NX_CHAR);
    NXputattr (nxhandle, "MCCODE_STRING", MCCODE_STRING, strlen(MCCODE_STRING), NX_CHAR);
    NXputattr (nxhandle, "name", parent, strlen(parent), NX_CHAR);
    
    return(NXclosedata(nxhandle));
  } else
  return(NX_ERROR);
} /* mcnxfile_instrcode */

/*******************************************************************************
* mcnxfile_section: begin/end a section
*                   open a new 'section' Group and an 'information' Data Set
*                   close the current Group when part="end"
*                   the Data Set 'information' is not closed here
* Returns: NX_ERROR or NX_OK
*******************************************************************************/
int mcnxfile_section(NXhandle nxhandle, char *part,
    char *pre,          /* %1$s  PRE  */
    char *type,         /* %2$s  TYP  */
    char *name,         /* %3$s  NAM  */
    char *valid_name,   /* %4$s  VNA  */
    char *parent,       /* %5$s  PAR  */
    char *valid_parent, /* %6$s  VPA  */
    int   level)        /* %7$i  LVL */
{
  if (!strcmp(part, "end")) {
    int ret;
    ret = NXclosegroup(nxhandle);
    if (ret == NX_ERROR)
      fprintf(stderr, "Warning: NeXus: could not close group %s (%s)\n", 
          valid_name, type);
    return(ret);
  }
  
  if (!strcmp(part, "begin")) {
    int length;
    char nxtype[CHAR_BUF_LENGTH];
    char nxname[CHAR_BUF_LENGTH];
    
    if (!strcmp(type, "instrument"))      strcpy(nxname, "instrument");
    else if (!strcmp(type, "simulation")) strcpy(nxname, "simulation");
    else strcpy(nxname, valid_name);
    snprintf(nxtype, CHAR_BUF_LENGTH, "NX%s", type);
    
    NXMDisableErrorReporting(); /* unactivate NeXus error messages */
    NXmakegroup(nxhandle, nxname, nxtype);
    NXMEnableErrorReporting();  /* enable NeXus error messages */
    
    if (NXopengroup(nxhandle, nxname, nxtype) == NX_ERROR) {
      fprintf(stderr, "Warning: NeXus: could not open group %s (%s) to store 'information'\n", 
          nxname, nxtype);
      return(NX_ERROR);
    }
    /* open a SDS to store attributes */
    snprintf(nxname, CHAR_BUF_LENGTH, "Information about %s of type %s is stored in attributes", name, nxtype);
    length = strlen(nxname);
    
    NXMDisableErrorReporting(); /* unactivate NeXus error messages */
    NXmakedata(nxhandle, "information", NX_CHAR, 1, &length);
    NXMEnableErrorReporting();  /* enable NeXus error messages */
    if (NXopendata(nxhandle, "information") == NX_ERROR) {
      fprintf(stderr, "Warning: NeXus: could not open 'information' for %s (%s)\n", 
          name, nxtype);
      return(NX_ERROR);
    }
    NXputdata (nxhandle, nxname);
    NXputattr(nxhandle, "name", name, strlen(name), NX_CHAR);
    NXputattr(nxhandle, "parent", parent, strlen(parent), NX_CHAR);

  }
  return(NX_OK);
} /* mcnxfile_section */

/*******************************************************************************
* mcnxfile_section: begin/end a data block (data/errors/events)
*                   open/close a 'part' Data Set
*                   open/close Axes (except for lists)
*                   handles compressed Data Set
* Returns: NX_ERROR or NX_OK
*******************************************************************************/
/* mcnxfile_datablock: data block begin/end. Returns: NX_ERROR or NX_OK */
int mcnxfile_data(NXhandle nxhandle, MCDETECTOR detector, char *part,
  char *valid_parent, char *valid_xlabel, char *valid_ylabel, char *valid_zlabel)
{
  /* write axes, only for data */
  if (strstr(part, "data")) {
    int i;
    if (!strstr(detector.format.Name, "list")) {
    
      if (detector.m > 1) {       /* X axis */
        double axis[detector.m];
        int dim=(int)detector.m;
        for(i = 0; i < detector.m; i++)
          axis[i] = detector.xmin+(detector.xmax-detector.xmin)*(i+0.5)/detector.m;
        if (strstr(mcnxversion,"compress") || strstr(mcnxversion,"zip"))
          NXcompmakedata(nxhandle, valid_xlabel, NX_FLOAT64, 1, &dim, NX_COMP_LZW, &dim);
        else
          NXmakedata(nxhandle, valid_xlabel, NX_FLOAT64, 1, &dim);

        if (NXopendata(nxhandle, valid_xlabel) == NX_ERROR) {
          fprintf(stderr, "Warning: could not open X axis %s in %s\n",
            valid_xlabel, detector.filename);
          return(NX_ERROR);
        }
        NXputdata (nxhandle, axis);
        NXputattr (nxhandle, "long_name", detector.xlabel, strlen(detector.xlabel), NX_CHAR);
        NXputattr (nxhandle, "short_name", detector.xvar, strlen(detector.xvar), NX_CHAR);
        int naxis=1;
        NXputattr (nxhandle, "axis", &naxis, 1, NX_INT32);
        NXputattr (nxhandle, "units", detector.xvar, strlen(detector.xvar), NX_CHAR);
        int nprimary=1;
        NXputattr (nxhandle, "primary", &nprimary, 1, NX_INT32);
        NXclosedata(nxhandle);
      }
      
      if (detector.n >= 1) {      /* Y axis */
        double axis[detector.n];
        int dim=(int)detector.n;
        for(i = 0; i < detector.n; i++)
          axis[i] = detector.ymin+(detector.ymax-detector.ymin)*(i+0.5)/detector.n;
        if (strstr(mcnxversion,"compress") || strstr(mcnxversion,"zip"))
          NXcompmakedata(nxhandle, valid_ylabel, NX_FLOAT64, 1, &dim, NX_COMP_LZW, &dim);
        else
          NXmakedata(nxhandle, valid_ylabel, NX_FLOAT64, 1, &dim);

        if (NXopendata(nxhandle, valid_ylabel) == NX_ERROR) {
          fprintf(stderr, "Warning: could not open Y axis %s in %s\n",
            valid_ylabel, detector.filename);
          return(NX_ERROR);
        }
        NXputdata (nxhandle, axis);
        NXputattr (nxhandle, "long_name", detector.ylabel, strlen(detector.ylabel), NX_CHAR);
        NXputattr (nxhandle, "short_name", detector.yvar, strlen(detector.yvar), NX_CHAR);
        int naxis=2;
        NXputattr (nxhandle, "axis", &naxis, 1, NX_INT32);
        NXputattr (nxhandle, "units", detector.yvar, strlen(detector.yvar), NX_CHAR);
        int nprimary=1;
        NXputattr (nxhandle, "primary", &nprimary, 1, NX_INT32);
        NXclosedata(nxhandle);
      }
      
      if (detector.p > 1) {     /* Z axis */
        double axis[detector.p];
        int dim=(int)detector.p;
        for(i = 0; i < detector.p; i++)
          axis[i] = detector.zmin+(detector.zmax-detector.zmin)*(i+0.5)/detector.p;
        if (strstr(mcnxversion,"compress") || strstr(mcnxversion,"zip"))
          NXcompmakedata(nxhandle, valid_zlabel, NX_FLOAT64, 1, &dim, NX_COMP_LZW, &dim);
        else
          NXmakedata(nxhandle, valid_zlabel, NX_FLOAT64, 1, &dim);

        if (NXopendata(nxhandle, valid_zlabel) == NX_ERROR) {
          fprintf(stderr, "Warning: could not open Z axis %s in %s\n",
            valid_ylabel, detector.filename);
          return(NX_ERROR);
        }
        NXputdata (nxhandle, axis);
        NXputattr (nxhandle, "long_name", detector.zlabel, strlen(detector.zlabel), NX_CHAR);
        NXputattr (nxhandle, "short_name", detector.zvar, strlen(detector.zvar), NX_CHAR);
        int naxis=3;
        NXputattr (nxhandle, "axis", &naxis, 1, NX_INT32);
        NXputattr (nxhandle, "units", detector.zvar, strlen(detector.zvar), NX_CHAR);
        int nprimary=1;
        NXputattr (nxhandle, "primary", &nprimary, 1, NX_INT32);
        NXclosedata(nxhandle);
      }
    } /* !list */
  } /* end format != list for data */
  
  /* write data */
  int dims[3]={detector.m,detector.n,detector.p};  /* number of elements to write */
  char *nxname=part;
  double *data;
  if (strstr(part,"data"))         { data=detector.p1; }
  else if (strstr(part,"errors"))  { data=detector.p2; }
  else if (strstr(part,"ncount"))  { data=detector.p0; }
  /* ignore errors for making/opening data (in case this has already been done */
  if (strstr(mcnxversion,"compress") || strstr(mcnxversion,"zip"))
    NXmakedata(nxhandle, nxname, NX_FLOAT64, detector.rank, dims);
  else
    NXcompmakedata(nxhandle, nxname, NX_FLOAT64, detector.rank, dims, NX_COMP_LZW, dims);

  if (NXopendata(nxhandle, nxname) == NX_ERROR) {
        fprintf(stderr, "Warning: could not open DataSet '%s' in %s\n",
          nxname, detector.filename);
        return(NX_ERROR);
      }
  NXputdata (nxhandle, data);
  NXputattr(nxhandle, "parent", valid_parent, strlen(valid_parent), NX_CHAR);
  int signal=1;
  if (strstr(part,"data")) {
    NXputattr(nxhandle, "signal", &signal, 1, NX_INT32);
    NXputattr(nxhandle, "short_name", detector.filename, strlen(detector.filename), NX_CHAR);
  }
  char nxtitle[CHAR_BUF_LENGTH];
  sprintf(nxtitle, "%s '%s'", nxname, detector.title);
  NXputattr(nxhandle, "long_name", nxtitle, strlen(nxtitle), NX_CHAR);
  /* first write attributes */
  char creator[CHAR_BUF_LENGTH];
  sprintf(creator, "%s/%s", mcinstrument_name, valid_parent);
  NXputattr(nxhandle, "creator", creator, strlen(creator), NX_CHAR);
  return(NXclosedata(nxhandle));
} /* mcnxfile_datablock */


#endif /* USE_NEXUS */
/* End of file "nexus-lib.c". */

#line 1374 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

#line 1 "mccode-r.c"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y/McXtrace X.Y
* Version: $Revision: 1.229 $
*
* Runtime system for McStas and McXtrace.
* Embedded within instrument in runtime mode.
* Contains SECTIONS:
*   MPI handling (sum, send, recv)
*   format definitions
*   I/O
*   mcdisplay support
*   random numbers
*   coordinates handling
*   vectors math (solve 2nd order, normals, randvec...)
*   parameter handling
*   signal and main handlers
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id: mcstas-r.c,v 1.229 2009/08/13 14:52:01 farhi Exp $
*
*******************************************************************************/

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/


/** Include header files to avoid implicit declarations (not allowed on LLVM) */
#include <ctype.h>
#include <sys/types.h>

// UNIX specific headers (non-Windows)
#if defined(__unix__) || defined(mac)
#include <unistd.h>
#include <sys/stat.h>
#endif


#ifndef DANSE
#ifdef MC_ANCIENT_COMPATIBILITY
int mctraceenabled = 0;
int mcdefaultmain  = 0;
#endif
/* else defined directly in the McStas generated C code */

static   long mcseed                 = 0; /* seed for random generator */
static   int  mcascii_only           = 0; /* flag for no header */
static   int  mcsingle_file          = 0; /* flag for storing all detectors into a single file */
static   long mcstartdate            = 0; /* start simulation time */
static   int  mcdisable_output_files = 0; /* --no-output-files */
mcstatic int  mcgravitation          = 0; /* use gravir=tation flag, for PROP macros */
int      mcMagnet                    = 0; /* megnet stack flag */
mcstatic int  mcdotrace              = 0; /* flag for --trace and messages for DISPLAY */
/* mcstatic FILE *mcsiminfo_file        = NULL; */
static   char *mcdirname             = NULL;      /* name of output directory */
static   char *mcsiminfo_name        = "mcstas";  /* default output sim file name */
int      mcallowbackprop             = 0;         /* flag to enable negative/backprop */
char*    mcDetectorCustomHeader      = NULL;      /* additional user output Tag in data files */
char*    mcopenedfiles               = "";        /* list of opened files (for append) */
long     mcopenedfiles_size          = 0;         /* size of that opened files list */
MCDETECTOR* mcDetectorArray          = NULL;      /* array of all opened detectors */
long     mcDetectorArray_size        = 0;         /* allocated detector array size */
long     mcDetectorArray_index       = 0;         /* current detector length (number of detectors so far) */

/* Number of particule histories to simulate. */
#ifdef NEUTRONICS
mcstatic unsigned long long int mcncount             = 1;
mcstatic unsigned long long int mcrun_num            = 0;
#else
mcstatic unsigned long long int mcncount             = 1e6;
mcstatic unsigned long long int mcrun_num            = 0;
#endif /* NEUTRONICS */

/* I/O structures */
mcstatic struct mcformats_struct mcformat;
mcstatic struct mcformats_struct mcformat_data;
#else
#include "mcstas-globals.h"
#endif /* !DANSE */

#ifndef MCSTAS_FORMAT
#define MCSTAS_FORMAT "McStas"  /* default format */
#endif

/* SECTION: parameters handling ============================================= */

/* Instrument input parameter type handling. */
/*******************************************************************************
* mcparm_double: extract double value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_double(char *s, void *vptr)
{
  char *p;
  double *v = (double *)vptr;

  if (!s) { *v = 0; return(1); }
  *v = strtod(s, &p);
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_double: display parameter type double
*******************************************************************************/
static char *
mcparminfo_double(char *parmname)
{
  return "double";
}

/*******************************************************************************
* mcparmerror_double: display error message when failed extract double
*******************************************************************************/
static void
mcparmerror_double(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for floating point parameter %s (mcparmerror_double)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_double: convert double to string
*******************************************************************************/
static void
mcparmprinter_double(char *f, void *vptr)
{
  double *v = (double *)vptr;
  sprintf(f, "%g", *v);
}

/*******************************************************************************
* mcparm_int: extract int value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_int(char *s, void *vptr)
{
  char *p;
  int *v = (int *)vptr;
  long x;

  if (!s) { *v = 0; return(1); }
  *v = 0;
  x = strtol(s, &p, 10);
  if(x < INT_MIN || x > INT_MAX)
    return 0;                        /* Under/overflow */
  *v = x;
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_int: display parameter type int
*******************************************************************************/
static char *
mcparminfo_int(char *parmname)
{
  return "int";
}

/*******************************************************************************
* mcparmerror_int: display error message when failed extract int
*******************************************************************************/
static void
mcparmerror_int(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for integer parameter %s (mcparmerror_int)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_int: convert int to string
*******************************************************************************/
static void
mcparmprinter_int(char *f, void *vptr)
{
  int *v = (int *)vptr;
  sprintf(f, "%d", *v);
}

/*******************************************************************************
* mcparm_string: extract char* value from 's' into 'vptr' (copy)
*******************************************************************************/
static int
mcparm_string(char *s, void *vptr)
{
  char **v = (char **)vptr;
  if (!s) { *v = NULL; return(1); }
  *v = (char *)malloc(strlen(s) + 1);
  if(*v == NULL)
  {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcparm_string).\n", (long)strlen(s) + 1));
  }
  strcpy(*v, s);
  return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_string: display parameter type string
*******************************************************************************/
static char *
mcparminfo_string(char *parmname)
{
  return "string";
}

/*******************************************************************************
* mcparmerror_string: display error message when failed extract string
*******************************************************************************/
static void
mcparmerror_string(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for string parameter %s (mcparmerror_string)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_string: convert string to string (including esc chars)
*******************************************************************************/
static void
mcparmprinter_string(char *f, void *vptr)
{
  char **v = (char **)vptr;
  char *p;

  if (!*v) { *f='\0'; return; }
  strcpy(f, "");
  for(p = *v; *p != '\0'; p++)
  {
    switch(*p)
    {
      case '\n':
        strcat(f, "\\n");
        break;
      case '\r':
        strcat(f, "\\r");
        break;
      case '"':
        strcat(f, "\\\"");
        break;
      case '\\':
        strcat(f, "\\\\");
        break;
      default:
        strncat(f, p, 1);
    }
  }
  /* strcat(f, "\""); */
} /* mcparmprinter_string */

/* now we may define the parameter structure, using previous functions */
static struct
  {
    int (*getparm)(char *, void *);
    char * (*parminfo)(char *);
    void (*error)(char *, char *);
    void (*printer)(char *, void *);
} mcinputtypes[] = {
  {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }, {
    mcparm_int, mcparminfo_int, mcparmerror_int,
    mcparmprinter_int
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }
};

/*******************************************************************************
* mcestimate_error: compute sigma from N,p,p2 in Gaussian large numbers approx
*******************************************************************************/
double mcestimate_error(double N, double p1, double p2)
{
  double pmean, n1;
  if(N <= 1)
    return p1;
  pmean = p1 / N;
  n1 = N - 1;
  /* Note: underflow may cause p2 to become zero; the fabs() below guards
     against this. */
  return sqrt((N/n1)*fabs(p2 - pmean*pmean));
}

double (*mcestimate_error_p)
  (double V2, double psum, double p2sum)=mcestimate_error;

/* SECTION: MPI handling ==================================================== */

#ifdef USE_MPI
/* MPI rank */
static int mpi_node_rank;
static int mpi_node_root = 0;


/*******************************************************************************
* mc_MPI_Reduce: Gathers arrays from MPI nodes using Reduce function.
*******************************************************************************/
int mc_MPI_Sum(double *sbuf, long count)
{
  if (!sbuf || count <= 0) return(MPI_ERR_COUNT);
  else {
    /* we must cut the buffer into blocks not exceeding the MPI max buffer size of 32000 */
    long   offset=0;
    double *rbuf=NULL;
    int    length=MPI_REDUCE_BLOCKSIZE; /* defined in mcstas.h */
    int    i=0;
    rbuf = calloc(count, sizeof(double));
    if (!rbuf)
      exit(-fprintf(stderr, "Error: Out of memory %li (mc_MPI_Sum)\n", count*sizeof(double)));
    while (offset < count) {
      if (!length || offset+length > count-1) length=count-offset;
      else length=MPI_REDUCE_BLOCKSIZE;
      MPI_Allreduce((double*)(sbuf+offset), (double*)(rbuf+offset),
              length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      offset += length;
    }

    for (i=0; i<count; i++) sbuf[i] = rbuf[i];
    free(rbuf);
  }
  return MPI_SUCCESS;
} /* mc_MPI_Sum */

/*******************************************************************************
* mc_MPI_Send: Send array to MPI node by blocks to avoid buffer limit
*******************************************************************************/
int mc_MPI_Send(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int dest)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mcstas.h */

  if (!sbuf || count <= 0) return(MPI_ERR_COUNT);
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    MPI_Send((void*)(sbuf+offset*dsize), length, dtype, dest, tag++, MPI_COMM_WORLD);
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Send */

/*******************************************************************************
* mc_MPI_Recv: Receives arrays from MPI nodes by blocks to avoid buffer limit
*             the buffer must have been allocated previously.
*******************************************************************************/
int mc_MPI_Recv(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int source)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mcstas.h */

  if (!sbuf || count <= 0) return(MPI_ERR_COUNT);
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    MPI_Recv((void*)(sbuf+offset*dsize), length, dtype, source, tag++,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Recv */

#endif /* USE_MPI */

/* SECTION: file i/o handling ================================================ */

/* Multiple output format support. ========================================== */
#ifdef USE_NEXUS
#define mcNUMFORMATS 10
#else
#define mcNUMFORMATS 9
#endif

/*******************************************************************************
* Definition of output formats. structure defined in mcstas-r.h
* Name aliases are defined in mcuse_format_* functions (below)
*******************************************************************************/

/* {Name, Extension, Header, Footer,
    BeginSection, EndSection, AssignTag;
    BeginData, EndData;
    BeginErrors, EndErrors;
    BeginNcount, EndNcount};
 */

mcstatic struct mcformats_struct mcformats[mcNUMFORMATS] = {
  { "McStas", "sim",
    "%PREFormat: %FMT file. Use mcplot/PGPLOT to view.\n"
      "%PREURL:    http://www.mccode.org/\n"
      "%PREEditor: %USR\n"
      "%PRECreator:%SRC simulation (" MCCODE_STRING ")\n"
      "%PREDate:   Simulation started (%DATL) %DAT\n"
      "%PREFile:   %FIL\n",
    "%PREEndDate:%DAT\n",
    "%PREbegin %TYP\n",
    "%PREend %TYP\n",
    "%PRE%NAM: %VAL\n",
    "%PREData [%PAR/%FIL] I: \n", "",
    "%PREErrors [%PAR/%FIL] E: \n", "",
    "%PREEvents [%PAR/%FIL] N: \n", "" },
  { "Scilab", "sci",
    "function mc_%VPA = get_%VPA(p)\n"
      "// %FMT function generated by McStas on %DAT\n"
      "// McStas simulation %SRC: %FIL %FMT\n"
      "// Import data using scilab> exec('%VPA.sci',-1); s=get_%VPA(); and s=get_%VPA('plot'); to plot\n\n"
      "mode(-1); //silent execution\n"
      "if argn(2) > 0, p=1; else p=0; end\n"
      "mc_%VPA = struct();\n"
      "mc_%VPA.Format ='%FMT';\n"
      "mc_%VPA.URL    ='http://www.mccode.org';\n"
      "mc_%VPA.Editor ='%USR';\n"
      "mc_%VPA.Creator='%SRC " MCCODE_STRING " simulation';\n"
      "mc_%VPA.Date   =%DATL; // for getdate\n"
      "mc_%VPA.File   ='%FIL';\n",
    "mc_%VPA.EndDate=%DATL; // for getdate\nendfunction\n"
    "function d=mcload_inline(d)\n"
      "// local inline func to load data\n"
      "execstr(['S=['+part(d.type,10:(length(d.type)-1))+'];']);\n"
      "if ~length(d.data)\n"
      " if ~length(strindex(d.format, 'binary'))\n"
      "  exec(d.filename,-1);p=d.parent;\n"
      "  if ~execstr('d2='+d.func+'();','errcatch'),d=d2; d.parent=p;end\n"
      " else\n"
      "  if length(strindex(d.format, 'float')), t='f';\n"
      "  elseif length(strindex(d.format, 'double')), t='d';\n"
      "  else return; end\n"
      "  fid=mopen(d.filename, 'rb');\n"
      "  pS = prod(S);\n"
      "  x = mget(3*pS, t, fid);\n"
      "  d.data  =matrix(x(1:pS), S);\n"
      "  if length(x) >= 3*pS,\n"
      "  d.errors=matrix(x((pS+1):(2*pS)), S);\n"
      "  d.events=matrix(x((2*pS+1):(3*pS)), S);end\n"
      "  mclose(fid);\n"
      "  return\n"
      " end\n"
      "end\n"
      "endfunction\n"
      "function d=mcplot_inline(d,p)\n"
      "// local inline func to plot data\n"
      "if ~length(strindex(d.type,'0d')), d=mcload_inline(d); end\n"
      "if ~p, return; end;\n"
      "execstr(['l=[',d.xylimits,'];']); S=size(d.data);\n"
      "t1=['['+d.parent+'] '+d.filename+': '+d.title];t = [t1;['  '+d.variables+'=['+d.values+']'];['  '+d.signal];['  '+d.statistics]];\n"
      "mprintf('%s\\n',t(:));\n"
      "if length(strindex(d.type,'0d')),return; end\n"
      "w=winsid();if length(w),w=w($)+1; else w=0; end\n"
      "xbasr(w); xset('window',w);\n"
      "if length(strindex(d.type,'2d'))\n"
      " if S(2) > 1, d.stepx=abs(l(1)-l(2))/(S(2)-1); else d.stepx=0; end\n"
      " if S(1) > 1, d.stepy=abs(l(3)-l(4))/(S(1)-1); else d.stepy=0; end\n"
      " d.x=linspace(l(1)+d.stepx/2,l(2)-d.stepx/2,S(2));\n"
      " d.y=linspace(l(3)+d.stepy/2,l(4)-d.stepy/2,S(1)); z=d.data;\n"
      " xlab=d.xlabel; ylab=d.ylabel; x=d.x; y=d.y;\n"
      " fz=max(abs(z));fx=max(abs(d.x));fy=max(abs(d.y));\n"
      " if fx>0,fx=round(log10(fx)); x=x/10^fx; xlab=xlab+' [*10^'+string(fx)+']'; end\n"
      " if fy>0,fy=round(log10(fy)); y=y/10^fy; ylab=ylab+' [*10^'+string(fy)+']'; end\n"
      " if fz>0,fz=round(log10(fz)); z=z/10^fz; t1=t1+' [*10^'+string(fz)+']'; end\n"
      " xset('colormap',hotcolormap(64));\n"
      " plot3d1(x,y,z',90,0,xlab+'@'+ylab+'@'+d.zlabel,[-1,2,4]); xtitle(t);\n"
      "else\n"
      " if max(S) > 1, d.stepx=abs(l(1)-l(2))/(max(S)-1); else d.stepx=0; end\n"
      " d.x=linspace(l(1)+d.stepx/2,l(2)-d.stepx/2,max(S));\n"
      " plot2d(d.x,d.data);xtitle(t,d.xlabel,d.ylabel);\n"
      "end\n"
      "xname(t1);\nendfunction\n"
    "mc_%VPA=get_%VPA();\n",
    "// Section %TYP [%NAM] (level %LVL)\n"
      "%PREt=[]; execstr('t=mc_%VNA.class','errcatch'); if ~length(t), mc_%VNA = struct(); end; mc_%VNA.class = '%TYP';",
    "%PREmc_%VPA.mc_%VNA = 0; mc_%VPA.mc_%VNA = mc_%VNA;\n",
    "%PREmc_%SEC.%NAM = '%VAL';\n",
    "%PREmc_%VPA.func='get_%VPA';\n%PREmc_%VPA.data = [ \n",
    " ]; // end of data\n%PREif length(mc_%VPA.data) == 0, single_file=0; else single_file=1; end\n%PREmc_%VPA=mcplot_inline(mc_%VPA,p);\n",
    "%PREerrors = [ \n",
    " ]; // end of errors\n%PREif single_file == 1, mc_%VPA.errors=errors; end\n",
    "%PREevents = [ \n",
    " ]; // end of events\n%PREif single_file == 1, mc_%VPA.events=events; end\n"},
  { "Matlab", "m",
    "function mc_%VPA = get_%VPA(p)\n"
      "%% %FMT function generated by McStas on %DAT\n"
      "%% McStas simulation %SRC: %FIL\n"
      "%% Import data using matlab> s=%VPA; and s=%VPA('plot'); to plot\n\n"
      "if nargout == 0 | nargin > 0, p=1; else p=0; end\n"
      "mc_%VPA.Format ='%FMT';\n"
      "mc_%VPA.URL    ='http://www.mccode.org';\n"
      "mc_%VPA.Editor ='%USR';\n"
      "mc_%VPA.Creator='%SRC " MCCODE_STRING " simulation';\n"
      "mc_%VPA.Date   =%DATL; %% for datestr\n"
      "mc_%VPA.File   ='%FIL';\n",
    "mc_%VPA.EndDate=%DATL; %% for datestr\n"
      "function d=mcload_inline(d)\n"
      "%% local inline function to load data\n"
      "S=d.type; eval(['S=[ ' S(10:(length(S)-1)) ' ];']);\n"
      "if isempty(d.data)\n"
      " if ~length(findstr(d.format, 'binary'))\n"
      "  if ~strcmp(d.filename,[d.func,'.m']) copyfile(d.filename,[d.func,'.m']); end\n"
      "  p=d.parent;path(path);\n"
      "  eval(['d=',d.func,';']);d.parent=p;\n"
      "  if ~strcmp(d.filename,[d.func,'.m']) delete([d.func,'.m']); end\n"
      " else\n"
      "  if length(findstr(d.format, 'float')), t='single';\n"
      "  elseif length(findstr(d.format, 'double')), t='double';\n"
      "  else return; end\n"
      "  if length(S) == 1, S=[S 1]; end\n"
      "  fid=fopen(d.filename, 'r');\n"
      "  pS = prod(S);\n"
      "  x = fread(fid, 3*pS, t);\n"
      "  d.data  =reshape(x(1:pS), S);\n"
      "  if prod(size(x)) >= 3*pS,\n"
      "  d.errors=reshape(x((pS+1):(2*pS)), S);\n"
      "  d.events=reshape(x((2*pS+1):(3*pS)), S);end\n"
      "  fclose(fid);\n"
      "  return\n"
      " end\n"
      "end\n"
      "return;\n"
      "function d=mcplot_inline(d,p)\n"
      "%% local inline function to plot data\n"
      "if isempty(findstr(d.type,'0d')), d=mcload_inline(d); end\nif ~p, return; end;\n"
      "eval(['l=[',d.xylimits,'];']); S=size(d.data);\n"
      "t1=['[',d.parent,'] ',d.filename,': ',d.title];t = strvcat(t1,['  ',d.variables,'=[',d.values,']'],['  ',d.signal],['  ',d.statistics]);\n"
      "disp(t);\n"
      "if ~isempty(findstr(d.type,'0d')), return; end\n"
      "figure; if ~isempty(findstr(d.type,'2d'))\n"
      " if S(2) > 1, d.stepx=abs(l(1)-l(2))/(S(2)-1); else d.stepx=0; end\n"
      " if S(1) > 1, d.stepy=abs(l(3)-l(4))/(S(1)-1); else d.stepy=0; end\n"
      " d.x=linspace(l(1)+d.stepx/2,l(2)-d.stepx/2,S(2));\n"
      " d.y=linspace(l(3)+d.stepy/2,l(4)-d.stepy/2,S(1));\n"
      " surface(d.x,d.y,d.data); xlim([l(1) l(2)]); ylim([l(3) l(4)]); shading flat;\n"
      "else\n"
      " if max(S) > 1, d.stepx=abs(l(1)-l(2))/(max(S)-1); else d.stepx=0; end\n"
      " d.x=linspace(l(1)+d.stepx/2,l(2)-d.stepx/2,max(S));\n"
      " plot(d.x,d.data); xlim([l(1) l(2)]);\n"
      "end\n"
      "xlabel(d.xlabel); ylabel(d.ylabel); title(t); \n"
      "set(gca,'position',[.18,.18,.7,.65]); set(gcf,'name',t1);grid on;\n"
      "if ~isempty(findstr(d.type,'2d')), colorbar; end\n",
    "%% Section %TYP [%NAM] (level %LVL)\n"
      "%PREmc_%VNA.class = '%TYP';",
    "mc_%VPA.mc_%VNA = mc_%VNA;\n",
    "%PREmc_%SEC.%NAM = '%VAL';\n",
    "%PREmc_%VPA.func='%VPA';\n%PREmc_%VPA.data = [ \n",
    " ]; %% end of data\nif length(mc_%VPA.data) == 0, single_file=0; else single_file=1; end\nmc_%VPA=mcplot_inline(mc_%VPA,p);\n",
    "%PREerrors = [ \n",
    " ]; %% end of errors\nif single_file, mc_%VPA.errors=errors; end\n",
    "%PREevents = [ \n",
    " ]; %% end of events\nif single_file, mc_%VPA.events=events; end\n"},
  { "IDL", "pro",
    "; %FMT function generated by McStas on %DAT\n"
      "; McStas simulation %SRC: %FIL\n"
      "; import using idl> s=%VPA() and s=%VPA(/plot) to plot\n\n"
      "function mcload_inline,d\n"
      "; local inline function to load external data\n"
      "S=d.type & a=execute('S=long(['+strmid(S,9,strlen(S)-10)+'])')\n"
      "if strpos(d.format, 'binary') lt 0 then begin\n"
      " p=d.parent\n"
      " x=read_binary(d.filename)\n"
      " get_lun, lun\n"
      " openw,lun,d.func+'.pro'\n"
      " writeu, lun,x\n"
      " free_lun,lun\n"
      " resolve_routine, d.func, /is_func, /no\n"
      " d=call_function(d.func)\n"
      "endif else begin\n"
      " if strpos(d.format, 'float') ge 0 then t=4 $\n"
      " else if strpos(d.format, 'double') ge 0 then t=5 $\n"
      " else return,d\n"
      " x=read_binary(d.filename, data_type=t)\n"
      " pS=n_elements(S)\nif pS eq 1 then pS=long(S) $\n"
      " else if pS eq 2 then pS=long(S(0)*S(1)) $\n"
      " else pS=long(S(0)*S(1)*S(2))\n"
      " pS=pS(0)\nstv,d,'data',reform(x(0:(pS-1)),S)\n"
      " d.data=transpose(d.data)\n"
      " if n_elements(x) ge long(3*pS) then begin\n"
      "  stv,d,'errors',reform(x(pS:(2*pS-1)),S)\n"
      "  stv,d,'events',reform(x((2*pS):(3*pS-1)),S)\n"
      "  d.errors=transpose(d.errors)\n"
      "  d.events=transpose(d.events)\n"
      " endif\n"
      "endelse\n"
      "return,d\nend ; FUN load\n"
    "function mcplot_inline,d,p\n"
      "; local inline function to plot data\n"
      "if size(d.data,/typ) eq 7 and strpos(d.type,'0d') lt 0 then d=mcload_inline(d)\n"
      "if p eq 0 or strpos(d.type,'0d') ge 0 then return, d\n"
      "S=d.type & a=execute('S=long(['+strmid(S,9,strlen(S)-10)+'])')\n"
      "stv,d,'data',reform(d.data,S,/over)\n"
      "if total(strpos(tag_names(d),'ERRORS')+1) gt 0 then begin\n"
      " stv,d,'errors',reform(d.errors,S,/over)\n"
      " stv,d,'events',reform(d.events,S,/over)\n"
      "endif\n"
      "d.xylimits=strjoin(strsplit(d.xylimits,' ',/extract),',') & a=execute('l=['+d.xylimits+']')\n"
      "t1='['+d.parent+'] '+d.filename+': '+d.title\n"
      "t=[t1,'  '+d.variables+'=['+d.values+']','  '+d.signal,'  '+d.statistics]\n"
      "print,t\n"
      "if strpos(d.type,'0d') ge 0 then return,d\n"
      "d.xlabel=strjoin(strsplit(d.xlabel,'`!\"^&*()-+=|\\,.<>/?@''~#{[}]',/extract),'_')\n"
      "d.ylabel=strjoin(strsplit(d.ylabel,'`!\"^&*()-+=|\\,.<>/?@''~#{[}]',/extract),'_')\n"
      "stv,d,'x',l(0)+indgen(S(0))*(l(1)-l(0))/S(0)\n"
      "if strpos(d.type,'2d') ge 0 then begin\n"
      "  name={DATA:d.func,IX:d.xlabel,IY:d.ylabel}\n"
      "  stv,d,'y',l(2)+indgen(S(1))*(l(3)-l(2))/S(1)\n"
      "  live_surface,d.data,xindependent=d.x,yindependent=d.y,name=name,reference_out=Win\n"
      "endif else begin\n"
      "  name={DATA:d.func,I:d.xlabel}\n"
      "  live_plot,d.data,independent=d.x,name=name,reference_out=Win\n"
      "endelse\n"
      "live_text,t,Window_In=Win.Win,location=[0.3,0.9]\n"
      "return,d\nend ; FUN plot\n"
    "pro stv,S,T,V\n"
      "; procedure set-tag-value that does S.T=V\n"
      "sv=size(V)\n"
      "T=strupcase(T)\n"
      "TL=strupcase(tag_names(S))\n"
      "id=where(TL eq T)\n"
      "sz=[0,0,0]\n"
      "vd=n_elements(sv)-2\n"
      "type=sv[vd]\n"
      "if id(0) ge 0 then d=execute('sz=SIZE(S.'+T+')')\n"
      "if (sz(sz(0)+1) ne sv(sv(0)+1)) or (sz(0) ne sv(0)) $\n"
      "  or (sz(sz(0)+2) ne sv(sv(0)+2)) $\n"
      "  or type eq 8 then begin\n"
      " ES = ''\n"
      " for k=0,n_elements(TL)-1 do begin\n"
      "  case TL(k) of\n"
      "   T:\n"
      "   else: ES=ES+','+TL(k)+':S.'+TL(k)\n"
      "  endcase\n"
      " endfor\n"
      " d=execute('S={'+T+':V'+ES+'}')\n"
      "endif else d=execute('S.'+T+'=V')\n"
      "end ; PRO stv\n"
    "function %VPA,plot=plot\n"
      "; %FMT function generated by McStas on %DAT\n"
      "; McStas simulation %SRC: %FIL\n"
      "; import using s=%VPA() and s=%VPA(/plot) to plot\n"
      "if keyword_set(plot) then p=1 else p=0\n"
      "%7$s={Format:'%FMT',URL:'http://www.mccode.org',"
      "Editor:'%USR',$\n"
      "Creator:'%SRC " MCCODE_STRING " simulation',$\n"
      "Date:%DATL,"
      "File:'%FIL'}\n",
    "stv,%VPA,'EndDate',%DATL ; for systime\nreturn, %VPA\nend\n",
    "; Section %TYP [%NAM] (level %LVL)\n"
      "%PRE%VNA={class:'%TYP'}\n",
    "%PREstv,%VPA,'%VNA',%VNA\n",
    "%PREstv,%SEC,'%NAM','%VAL'\n",
    "%PREstv,%VPA,'func','%VPA' & data=[ $\n",
    " ]\n%PREif size(data,/type) eq 7 then single_file=0 else single_file=1\n"
    "%PREstv,%VPA,'data',data & data=0 & %VPA=mcplot_inline(%VPA,p)\n",
    "%PREif single_file ne 0 then begin errors=[ ",
    " ]\n%PREstv,%VPA,'errors',reform(errors,%MDIM,%NDIM,/over) & errors=0\n%PREendif\n",
    "%PREif single_file ne 0 then begin events=[ ",
    " ]\n%PREstv,%VPA,'events',reform(events,%MDIM,%NDIM,/over) & events=0\n%PREendif\n\n"},
  { "XML", "xml",
    "<?xml version=\"1.0\" ?>\n<!--\n"
      "URL:    http://www.nexusformat.org/\n"
      "Editor: %USR\n"
      "Creator:%SRC " MCCODE_STRING " [www.mccode.org].\n"
      "Date:   Simulation started (%DATL) %DAT\n"
      "File:   %FIL\n"
      "View with Mozilla, InternetExplorer, gxmlviewer, kxmleditor\n-->\n"
      "<NX%PAR file_name=\"%FIL\" file_time=\"%DAT\" user=\"%USR\">\n"
        "<NXentry name=\"" MCCODE_STRING "\"><start_time>%DAT</start_time>\n",
    "<end_time>%DAT</end_time></NXentry></NX%PAR>\n<!-- EndDate:%DAT -->\n",
    "%PRE<NX%TYP name=\"%NAM\">\n",
    "%PRE</NX%TYP>\n",
    "%PRE<%NAM>%VAL</%NAM>\n",
    "%PRE<%XVL long_name=\"%XLA\" axis=\"1\" primary=\"1\" min=\"%XMIN\""
        " max=\"%XMAX\" dims=\"%MDIM\" range=\"1\"></%XVL>\n"
      "%PRE<%YVL long_name=\"%YLA\" axis=\"2\" primary=\"1\" min=\"%YMIN\""
        " max=\"%YMAX\" dims=\"%NDIM\" range=\"1\"></%YVL>\n"
      "%PRE<%ZVL long_name=\"%ZLA\" axis=\"3\" primary=\"1\" min=\"%ZMIN\""
        " max=\"%ZMAX\" dims=\"%PDIM\" range=\"1\"></%ZVL>\n"
      "%PRE<data long_name=\"%TITL\" signal=\"1\"  axis=\"[%XVL,%YVL,%ZVL]\" file_name=\"%FIL\">\n",
    "%PRE</data>\n",
    "%PRE<errors>\n", "%PRE</errors>\n",
    "%PRE<monitor>\n", "%PRE</monitor>\n"},
  { "HTML", "html",
    "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD %DAT//EN\"\n"
      "\"http://www.w3.org/TR/html4/strict.dtd\">\n"
      "<HTML><HEAD><META name=\"Author\" content=\"%PAR\">\n"
      "<META name=\"Creator\" content=\"%PAR (%SRC) " MCCODE_STRING " [www.mccode.org] simulation\">\n"
      "<META name=\"Date\" content=\"%DAT\">\n"
      "<TITLE>[McStas %PAR (%SRC)]%FIL</TITLE></HEAD>\n"
      "<BODY><center><h1><a name=\"%PAR\">"
        "McStas simulation %SRC (%SRC): Result file %FIL.html</a></h1></center><br>\n"
        "This simulation report was automatically created by"
        " <a href=\"http://www.mccode.org/\"><i>" MCCODE_STRING "</i></a><br>\n"
        "<pre>User:   %USR<br>\n"
        "%PRECreator: <a href=\"%SRC\">%SRC</a> %PAR McStas simulation<br>\n"
        "%PREFormat:  %FMT<br>\n"
        "%PREDate:    (%DATL) %DAT<br></pre>\n"
        "VRML viewers may be obtained at <a href=\"http://cic.nist.gov/vrml/vbdetect.html\">http://cic.nist.gov/vrml/vbdetect.html</a>\n",
    "<b>EndDate: </b>(%DATL) %DAT<br></BODY></HTML>\n",
    "%PRE<h%LVL><a name=\"%NAM\">%TYP %NAM</a></h%LVL> "
      "[child of <a href=\"#%PAR\">%PAR</a>]<br>\n",
    "[end of <a href=\"#%NAM\">%TYP %NAM</a>]<br>\n",
    "%PRE<b>%NAM: </b>%VAL<br>\n",
    "%PRE<b>DATA</b><br><center><embed src=\"%FIL\" type=\"model/vrml\" width=\"75%%\" height=\"50%%\"></embed><br>File <a href=\"%FIL\">%FIL [VRML format]</a></center><br>\n", "%PREEnd of DATA<br>\n",
    "%PRE<b>ERRORS</b><br>\n","%PREEnd of ERRORS<br>\n",
    "%PRE<b>EVENTS</b><br>\n", "%PREEnd of EVENTS<br>\n"},
  { "Octave", "m",
    "function mc_%VPA = get_%VPA(p)\n"
      "%% %FMT function generated by McStas on %DAT\n"
      "%% McStas simulation %SRC: %FIL\n"
      "%% Import data using octave> s=%VPA(); and plot with s=%VPA('plot');\n"
      "if nargin > 0, p=1; else p=0; end\n"
      "mc_%VPA.Format ='%FMT';\n"
      "mc_%VPA.URL    ='http://www.mccode.org';\n"
      "mc_%VPA.Editor ='%USR';\n"
      "mc_%VPA.Creator='%SRC " MCCODE_STRING " simulation';\n"
      "mc_%VPA.Date   =%DATL; %% for datestr\n"
      "mc_%VPA.File   ='%FIL';\n",
    "mc_%VPA.EndDate=%DATL; %% for datestr\nendfunction\n"
      "if exist('mcload_inline'), return; end\n"
      "function d=mcload_inline(d)\n"
      "%% local inline function to load data\n"
      "S=d.type; eval(['S=[ ' S(10:(length(S)-1)) ' ];']);\n"
      "if isempty(d.data)\n"
      " if ~length(findstr(d.format, 'binary'))\n"
      "  source(d.filename);p=d.parent;\n"
      "  eval(['d=get_',d.func,';']);d.parent=p;\n"
      " else\n"
      "  if length(findstr(d.format, 'float')), t='float';\n"
      "  elseif length(findstr(d.format, 'double')), t='double';\n"
      "  else return; end\n"
      "  if length(S) == 1, S=[S 1]; end\n"
      "  fid=fopen(d.filename, 'r');\n"
      "  pS = prod(S);\n"
      "  x = fread(fid, 3*pS, t);\n"
      "  d.data  =reshape(x(1:pS), S);\n"
      "  if prod(size(x)) >= 3*pS,\n"
      "  d.errors=reshape(x((pS+1):(2*pS)), S);\n"
      "  d.events=reshape(x((2*pS+1):(3*pS)), S);end\n"
      "  fclose(fid);\n"
      "  return\n"
      " end\n"
      "end\n"
      "return;\nendfunction\n\n"
      "function d=mcplot_inline(d,p)\n"
      "%% local inline function to plot data\n"
      "if isempty(findstr(d.type,'0d')), d=mcload_inline(d); end\nif ~p, return; end;\n"
      "eval(['l=[',d.xylimits,'];']); S=size(d.data);\n"
      "t1=['[',d.parent,'] ',d.filename,': ',d.title];t = strcat(t1,['  ',d.variables,'=[',d.values,']'],['  ',d.signal],['  ',d.statistics]);\n"
      "disp(t);\n"
      "if ~isempty(findstr(d.type,'0d')), return; end\n"
      "xlabel(d.xlabel); ylabel(d.ylabel); title(t);"
      "figure; if ~isempty(findstr(d.type,'2d'))\n"
      " if S(2) > 1, d.stepx=abs(l(1)-l(2))/(S(2)-1); else d.stepx=0; end\n"
      " if S(1) > 1, d.stepy=abs(l(3)-l(4))/(S(1)-1); else d.stepy=0; end\n"
      " d.x=linspace(l(1)+d.stepx/2,l(2)-d.stepx/2,S(2));\n"
      " d.y=linspace(l(3)+d.stepy/2,l(4)-d.stepy/2,S(1));\n"
      " mesh(d.x,d.y,d.data);\n"
      "else\n"
      " if max(S) > 1, d.stepx=abs(l(1)-l(2))/(max(S)-1); else d.stepx=0; end\n"
      " d.x=linspace(l(1)+d.stepx/2,l(2)-d.stepx/2,max(S));\n"
      " plot(d.x,d.data);\n"
      "end\nendfunction\n",
    "%% Section %TYP [%NAM] (level %LVL)\n"
      "mc_%VNA.class = '%TYP';",
    "mc_%VPA.mc_%VNA = mc_%VNA;\n",
    "mc_%SEC.%NAM = '%VAL';\n",
    "mc_%VPA.func='%VPA';\n%PREmc_%VPA.data = [ \n",
    " ]; %% end of data\nif length(mc_%VPA.data) == 0, single_file=0; else single_file=1; end\nmc_%VPA=mcplot_inline(mc_%VPA,p);\n",
    "errors = [ \n",
    " ]; %% end of errors\nif single_file, mc_%VPA.errors=errors; end\n",
    "events = [ \n",
    " ]; %% end of events\nif single_file, mc_%VPA.events=events; end\n"},
  { "VRML", "wrl",
    "#VRML V2.0 utf8\n# Format: %FMT file\n"
      "use freeWRL, openvrml, vrmlview, CosmoPlayer, Cortona, Octaga... to view file\n"
      "WorldInfo {\n"
      "title \"%SRC/%FIL simulation Data\"\n"
      "info [ \"URL:    http://www.mccode.org/\"\n"
      "       \"Editor: %USR\"\n"
      "       \"Creator:%SRC simulation (" MCCODE_STRING ")\"\n"
      "       \"Date:   Simulation started (%DATL) %DAT\"\n"
      "       \"File:   %FIL\" ]\n}\n"
      "Background { skyAngle [ 1.57 1.57] skyColor [0 0 1, 1 1 1, 0.1 0 0] }\n",
    "# EndDate:%DAT\n",
    "# begin %TYP %PAR\n",
    "# end %TYP %PAR\n",
    "%PRE%SEC.%NAM= '%VAL'\n",
    "# The Proto that contains data values and objects to plot these\n"
      "PROTO I_ERR_N_%VPA [\n"
      "# the PROTO parameters\n"
      "  field MFFloat Data [ ]\n"
      "  field MFFloat Errors [ ]\n"
      "  field MFFloat Ncounts [ ]\n"
      "] { # The plotting objects/methods in the Proto\n"
      "  # draw a small sphere at the origin\n"
      "  DEF Data_%VPA Group {\n"
      "  children [\n"
      "    DEF CoordinateOrigin Group {\n"
      "      children [\n"
      "        Transform { translation  0 0 0 }\n"
      "        Shape { \n"
      "          appearance Appearance { \n"
      "            material Material {\n"
      "              diffuseColor 1.0 1.0 0.0\n"
      "              transparency 0.5 } }\n"
      "          geometry Sphere { radius .01 } \n"
      "    } ] }\n"
      "    # defintion of the arrow allong Y axis\n"
      "    DEF ARROW Group {\n"
      "      children [\n"
      "        Transform {\n"
      "          translation 0 0.5 0\n"
      "          children [\n"
      "            Shape {\n"
      "              appearance DEF ARROW_APPEARANCE Appearance {\n"
      "                material Material {\n"
      "                  diffuseColor .3 .3 1\n"
      "                  emissiveColor .1 .1 .33\n"
      "                }\n"
      "              }\n"
      "              geometry Cylinder {\n"
      "                bottom FALSE\n"
      "                radius .005\n"
      "                height 1\n"
      "                top FALSE\n"
      "        } } ] }\n"
      "        Transform {\n"
      "          translation 0 1 0\n"
      "          children [\n"
      "            DEF ARROW_POINTER Shape {\n"
      "              geometry Cone {\n"
      "                bottomRadius .05\n"
      "                height .1\n"
      "              }\n"
      "              appearance USE ARROW_APPEARANCE\n"
      "    } ] } ] }\n"
      "    # the arrow along X axis\n"
      "    Transform {\n"
      "      translation 0 0 0\n"
      "      rotation 1 0 0 1.57\n"
      "      children [\n"
      "        Group {\n"
      "          children [ \n"
      "            USE ARROW\n"
      "    ] } ] }\n"
      "    # the arrow along Z axis\n"
      "    Transform {\n"
      "      translation 0 0 0\n"
      "      rotation 0 0 1 -1.57\n"
      "      children [\n"
      "        Group {\n"
      "          children [ \n"
      "            USE ARROW\n"
      "    ] } ] }\n"
      "    # the Y label (which is vertical)\n"
      "    DEF Y_Label Group {\n"
      "      children [\n"
      "        Transform {\n"
      "          translation 0 1 0\n"
      "          children [\n"
      "            Billboard {\n"
      "              children [\n"
      "                Shape {\n"
      "                  appearance DEF LABEL_APPEARANCE Appearance {\n"
      "                    material Material {\n"
      "                      diffuseColor 1 1 .3\n"
      "                      emissiveColor .33 .33 .1\n"
      "                    } }\n"
      "                  geometry Text { \n"
      "                    string [ \"%ZVAR: %ZLA\", \"%ZMIN:%ZMAX - %PDIM points\" ] \n"
      "                    fontStyle FontStyle {  size .2 }\n"
      "    } } ] } ] } ] }\n"
      "    # the X label\n"
      "    DEF X_Label Group {\n"
      "      children [\n"
      "        Transform {\n"
      "          translation 1 0 0\n"
      "          children [\n"
      "            Billboard {\n"
      "              children [\n"
      "                Shape {\n"
      "                  appearance DEF LABEL_APPEARANCE Appearance {\n"
      "                    material Material {\n"
      "                      diffuseColor 1 1 .3\n"
      "                      emissiveColor .33 .33 .1\n"
      "                    } }\n"
      "                  geometry Text { \n"
      "                    string [ \"%XVAR: %XLA\", \"%XMIN:%XMAX - %MDIM points\" ] \n"
      "                    fontStyle FontStyle {  size .2 }\n"
      "    } } ] } ] } ] }\n"
      "    # the Z label\n"
      "    DEF Z_Label Group {\n"
      "      children [\n"
      "        Transform {\n"
      "          translation 0 0 1\n"
      "          children [\n"
      "            Billboard {\n"
      "              children [\n"
      "                Shape {\n"
      "                  appearance DEF LABEL_APPEARANCE Appearance {\n"
      "                    material Material {\n"
      "                      diffuseColor 1 1 .3\n"
      "                      emissiveColor .33 .33 .1\n"
      "                    } }\n"
      "                  geometry Text { \n"
      "                    string [ \"%YVAR: %YLA\", \"%YMIN:%YMAX - %NDIM points\" ] \n"
      "                    fontStyle FontStyle {  size .2 }\n"
      "    } } ] } ] } ] }\n"
      "    # The text information (header data )\n"
      "    DEF Header Group {\n"
      "      children [\n"
      "        Transform {\n"
      "          translation 0 2 0\n"
      "          children [\n"
      "            Billboard {\n"
      "              children [\n"
      "                Shape {\n"
      "                  appearance Appearance {\n"
      "                    material Material { \n"
      "                      diffuseColor .9 0 0\n"
      "                      emissiveColor .9 0 0 }\n"
      "                  }\n"
      "                  geometry Text {\n"
      "                    string [ \"%PAR/%FIL\",\"%TITL\" ]\n"
      "                    fontStyle FontStyle {\n"
      "                        style \"BOLD\"\n"
      "                        size .2\n"
      "    } } } ] } ] } ] }\n"
      "    # The Data plot\n"
      "    DEF MonitorData Group {\n"
      "      children [\n"
      "        DEF TransformData Transform {\n"
      "          children [\n"
      "            Shape {\n"
      "              appearance Appearance {\n"
      "                material Material { emissiveColor 0 0.2 0 }\n"
      "              }\n"
      "              geometry ElevationGrid {\n"
      "                xDimension  %MDIM\n"
      "                zDimension  %NDIM\n"
      "                xSpacing    1\n"
      "                zSpacing    1\n"
      "                solid       FALSE\n"
      "                height IS Data\n"
      "    } } ] } ] }\n"
      "    # The VRMLScript that redimension x and z axis within 0:1\n"
      "    # and re-scale data within 0:1\n"
      "    DEF GetScale Script {\n"
      "      eventOut SFVec3f scale_vect\n"
      "      url \"javascript: \n"
      "        function initialize( ) {\n"
      "          scale_vect = new SFVec3f(1.0/%MDIM, 1.0/Math.abs(%ZMAX-%ZMIN), 1.0/%NDIM); }\" }\n"
      "  ] }\n"
      "ROUTE GetScale.scale_vect TO TransformData.scale\n} # end of PROTO\n"
      "# now we call the proto with Data values\n"
      "I_ERR_N_%VPA {\nData [\n",
    "] # End of Data\n",
    "Errors [\n",
    "] # End of Errors\n",
    "Ncounts [\n",
    "] # End of Ncounts\n}" },
    { "Python", "py",
    "from numpy import array\ndef %VPA(p):\n"
      "  '''\n  %PAR (%FIL) procedure generated from McStas on %DAT\n\n"
      "  McStas simulation %SRC: %FIL\n"
      "  import data using python> s=%VPA(p); with p=0/1 for no plot/plot\n  '''\n"
      "  %VPA={'Format':'%FMT','URL':'http://www.mccode.org',\\\n"
      "  'Editor':'%USR',\\\n"
      "  'Creator':'%SRC " MCCODE_STRING " [www.mccode.org]',\\\n"
      "  'Date':%DATL,\\\n"
      "  'File':'%FIL'}\n",
    "  %VPA['EndDate']=%DATL\nreturn %VPA\n",
    "  # Section Section %TYP [%NAM] (level %LVL)\n"
      "  %VNA['class'] = '%TYP'\n",
    "  %VPA['%VNA'] = %VNA\n  del %VNA\n",
    "  %SEC['%NAM'] = '%VAL'\n",
    "  %VPA['func'] ='%VPA'\n  %VPA['data'] = [ ",
    "  ] # end of data\n",
    "  %VPA['errors'] = [ ",
    "  ] # end of errors\n",
    "  %VPA['ncount'] = [ ",
    "  ] # end of ncount\n",
    }
#ifdef USE_NEXUS
    ,
    { "NeXus", "nxs",
    "%PREFormat: %FMT file. Use hdfview to view.\n"
      "%PREURL:    http://www.mccode.org/\n"
      "%PREEditor: %USR\n"
      "%PRECreator:%SRC simulation (" MCCODE_STRING ")\n"
      "%PREDate:   Simulation started (%DATL) %DAT\n"
      "%PREFile:   %FIL\n",
    "%PREEndDate:%DAT\n",
    "%PREbegin %TYP\n",
    "%PREend %TYP\n",
    "%PRE%NAM: %VAL\n",
    "", "",
    "%PREErrors [%PAR/%FIL]: \n", "",
    "%PREEvents [%PAR/%FIL]: \n", "" }
#endif
};

/*******************************************************************************
* mcfull_file: allocates a full file name=mcdirname+file. Catenate extension if missing.
*******************************************************************************/
char *mcfull_file(char *name, char *ext)
{
  int   dirlen=0;
  char *mem   =NULL;

  dirlen = mcdirname ? strlen(mcdirname) : 0;
  mem = malloc(dirlen + strlen(name) + CHAR_BUF_LENGTH);
  if(!mem) {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcfull_file)\n", (long)(dirlen + strlen(name) + 256)));
  }
  strcpy(mem, "");

  /* prepend directory name to path if name does not contain a path */
  if (dirlen > 0 && !strchr(name, MC_PATHSEP_C)) {
    strcat(mem, mcdirname);
    strcat(mem, MC_PATHSEP_S);
  } /* dirlen */

  strcat(mem, name);
  if (!strchr(name, '.') && ext && strlen(ext))
  { /* add extension if not in file name already */
    strcat(mem, ".");
    strcat(mem, ext);
  }
  return(mem);
} /* mcfull_file */

/*******************************************************************************
* mcnew_file: opens a new file within mcdirname if non NULL
*             if mode is non 0, then mode is used, else mode is 'w'
*******************************************************************************/
FILE *mcnew_file(char *name, char *ext, char *mode)
{
  char *mem;
  FILE *file;

  if (!name || strlen(name) == 0 || mcdisable_output_files) return(NULL);

  mem = mcfull_file(name, ext);
  file = fopen(mem, (mode ? mode : "w"));
  if(!file)
    fprintf(stderr, "Warning: could not open output file '%s' in mode '%s' (mcnew_file)\n", mem, mode);
  else {
    if (!mcopenedfiles ||
        (mcopenedfiles && mcopenedfiles_size <= strlen(mcopenedfiles)+strlen(mem))) {
      mcopenedfiles_size+=CHAR_BUF_LENGTH;
      if (!mcopenedfiles || !strlen(mcopenedfiles))
        mcopenedfiles = calloc(1, mcopenedfiles_size);
      else
        mcopenedfiles = realloc(mcopenedfiles, mcopenedfiles_size);
    }
    strcat(mcopenedfiles, " ");
    strcat(mcopenedfiles, mem);
  }
  free(mem);

  return file;
} /* mcnew_file */

/*******************************************************************************
* str_rep: Replaces a token by an other (of SAME length) in a string
* This function modifies 'string'
*******************************************************************************/
char *str_rep(char *string, char *from, char *to)
{
  char *p;

  if (!string || !strlen(string)) return(string);
  if (strlen(from) != strlen(to)) return(string);

  p   = string;

  while (( p = strstr(p, from) ) != NULL) {
    long index;
    for (index=0; index<strlen(to); index++) p[index]=to[index];
  }
  return(string);
} /* str_rep */

#define VALID_NAME_LENGTH 64
/*******************************************************************************
* mcvalid_name: makes a valid string for variable names.
*   copy 'original' into 'valid', replacing invalid characters by '_'
*   char arrays must be pre-allocated. n can be 0, or the maximum number of
*   chars to be copied/checked
*******************************************************************************/
static char *mcvalid_name(char *valid, char *original, int n)
{
  long i;


  if (original == NULL || strlen(original) == 0)
  { strcpy(valid, "noname"); return(valid); }
  if (n <= 0) n = strlen(valid);

  if (n > strlen(original)) n = strlen(original);
  else original += strlen(original)-n;
  strncpy(valid, original, n);

  for (i=0; i < n; i++)
  {
    if ( (valid[i] > 122)
      || (valid[i] < 32)
      || (strchr("!\"#$%&'()*+,-.:;<=>?@[\\]^`/ \n\r\t", valid[i]) != NULL) )
    {
      if (i) valid[i] = '_'; else valid[i] = 'm';
    }
  }
  valid[i] = '\0';

  return(valid);
} /* mcvalid_name */

/*******************************************************************************
* pfprintf: just as fprintf with positional arguments %N$t,
*   but with (char *)fmt_args being the list of arg type 't'.
*   Needed as the vfprintf is not correctly handled on some platforms.
*   1- look for the maximum %d$ field in fmt
*   2- look for all %d$ fields up to max in fmt and set their type (next alpha)
*   3- retrieve va_arg up to max, and save pointer to arg in local arg array
*   4- use strchr to split around '%' chars, until all pieces are written
*   returns number of arguments written.
* Warning: this function is restricted to only handles types t=s,g,i,li
*          without additional field formating, e.g. %N$t
*******************************************************************************/
static int pfprintf(FILE *f, char *fmt, char *fmt_args, ...)
{
  #define MyNL_ARGMAX 50
  char  *fmt_pos;

  char *arg_char[MyNL_ARGMAX];
  int   arg_int[MyNL_ARGMAX];
  long  arg_long[MyNL_ARGMAX];
  double arg_double[MyNL_ARGMAX];

  char *arg_posB[MyNL_ARGMAX];  /* position of '%' */
  char *arg_posE[MyNL_ARGMAX];  /* position of '$' */
  char *arg_posT[MyNL_ARGMAX];  /* position of type */

  int   arg_num[MyNL_ARGMAX];   /* number of argument (between % and $) */
  int   this_arg=0;
  int   arg_max=0;
  va_list ap;

  if (!f || !fmt_args || !fmt) return(-1);
  for (this_arg=0; this_arg<MyNL_ARGMAX;  arg_num[this_arg++] =0); this_arg = 0;
  fmt_pos = fmt;
  while(1)  /* analyse the format string 'fmt' */
  {
    char *tmp;

    arg_posB[this_arg] = (char *)strchr(fmt_pos, '%');
    tmp = arg_posB[this_arg];
    if (tmp)	/* found a percent */
    {
      char  printf_formats[]="dliouxXeEfgGcs\0";
      arg_posE[this_arg] = (char *)strchr(tmp, '$');
      if (arg_posE[this_arg] && isdigit(tmp[1]))
      { /* found a dollar following a percent  and a digit after percent */
        char  this_arg_chr[10];


        /* extract positional argument index %*$ in fmt */
        strncpy(this_arg_chr, arg_posB[this_arg]+1, arg_posE[this_arg]-arg_posB[this_arg]-1 < 10 ? arg_posE[this_arg]-arg_posB[this_arg]-1 : 9);
        this_arg_chr[arg_posE[this_arg]-arg_posB[this_arg]-1] = '\0';
        arg_num[this_arg] = atoi(this_arg_chr);
        if (arg_num[this_arg] <=0 || arg_num[this_arg] >= MyNL_ARGMAX)
          return(-fprintf(stderr,"pfprintf: invalid positional argument number (%i is <=0 or >=%i) from '%s'.\n", arg_num[this_arg], MyNL_ARGMAX, this_arg_chr));
        /* get type of positional argument: follows '%' -> arg_posE[this_arg]+1 */
        fmt_pos = arg_posE[this_arg]+1;
        fmt_pos[0] = tolower(fmt_pos[0]);
        if (!strchr(printf_formats, fmt_pos[0]))
          return(-fprintf(stderr,"pfprintf: invalid positional argument type (%c != expected %c).\n", fmt_pos[0], fmt_args[arg_num[this_arg]-1]));
        if (fmt_pos[0] == 'l' && (fmt_pos[1] == 'i' || fmt_pos[1] == 'd')) fmt_pos++;
        arg_posT[this_arg] = fmt_pos;
        /* get next argument... */
        this_arg++;
      }
      else
      { /* no dollar or no digit */
        if  (tmp[1] == '%') {
          fmt_pos = arg_posB[this_arg]+2;  /* found %% */
        } else if (strchr(printf_formats,tmp[1])) {
          fmt_pos = arg_posB[this_arg]+1;  /* found %s */
        } else {
          return(-fprintf(stderr,"pfprintf: must use only positional arguments (%s).\n", arg_posB[this_arg]));
        }
      }
    } else
      break;  /* no more % argument */
  }
  arg_max = this_arg;
  /* get arguments from va_arg list, according to their type */
  va_start(ap, fmt_args);
  for (this_arg=0; this_arg<strlen(fmt_args); this_arg++)
  {

    switch(tolower(fmt_args[this_arg]))
    {
      case 's':                       /* string */
              arg_char[this_arg] = va_arg(ap, char *);
              break;
      case 'd':
      case 'i':
      case 'c':                      /* int */
              arg_int[this_arg] = va_arg(ap, int);
              break;
      case 'l':                       /* long int */
              arg_long[this_arg] = va_arg(ap, long int);
              break;
      case 'f':
      case 'g':
      case 'e':                      /* double */
              arg_double[this_arg] = va_arg(ap, double);
              break;
      default: fprintf(stderr,"pfprintf: argument type is not implemented (arg %%%i$ type %c).\n", this_arg+1, fmt_args[this_arg]);
    }
  }
  va_end(ap);
  /* split fmt string into bits containing only 1 argument */
  fmt_pos = fmt;
  for (this_arg=0; this_arg<arg_max; this_arg++)
  {
    char *fmt_bit;
    int   arg_n;

    if (arg_posB[this_arg]-fmt_pos>0)
    {
      fmt_bit = (char*)malloc(arg_posB[this_arg]-fmt_pos+10);
      if (!fmt_bit) return(-fprintf(stderr,"pfprintf: not enough memory.\n"));
      strncpy(fmt_bit, fmt_pos, arg_posB[this_arg]-fmt_pos);
      fmt_bit[arg_posB[this_arg]-fmt_pos] = '\0';
      fprintf(f, "%s", fmt_bit); /* fmt part without argument */
    } else
    {
      fmt_bit = (char*)malloc(10);
      if (!fmt_bit) return(-fprintf(stderr,"pfprintf: not enough memory.\n"));
    }
    arg_n = arg_num[this_arg]-1; /* must be >= 0 */
    strcpy(fmt_bit, "%");
    strncat(fmt_bit, arg_posE[this_arg]+1, arg_posT[this_arg]-arg_posE[this_arg]);
    fmt_bit[arg_posT[this_arg]-arg_posE[this_arg]+1] = '\0';

    switch(tolower(fmt_args[arg_n]))
    {
      case 's': fprintf(f, fmt_bit, arg_char[arg_n]);
                break;
      case 'd':
      case 'i':
      case 'c':                      /* int */
              fprintf(f, fmt_bit, arg_int[arg_n]);
              break;
      case 'l':                       /* long */
              fprintf(f, fmt_bit, arg_long[arg_n]);
              break;
      case 'f':
      case 'g':
      case 'e':                       /* double */
              fprintf(f, fmt_bit, arg_double[arg_n]);
              break;
    }
    fmt_pos = arg_posT[this_arg]+1;
    if (this_arg == arg_max-1)
    { /* add eventual leading characters for last parameter */
      if (fmt_pos < fmt+strlen(fmt))
        fprintf(f, "%s", fmt_pos);
    }
    if (fmt_bit) free(fmt_bit);

  }
  return(this_arg);
} /* pfprintf */

/*******************************************************************************
* mcinfo_header: output header/footer tags/info using specific file format.
*   the 'part' may be 'header' or 'footer' depending on part to write.
* RETURN: negative==error, positive=all went fine
* Used by: mcsiminfo_init, mcsiminfo_close
*******************************************************************************/
static int mcinfo_header(MCDETECTOR detector, char *part)
{
  char*HeadFoot;              /* format string */
  long date_l;                /* date as a long number */
  char date[CHAR_BUF_LENGTH]; /* date as a string */
  char valid_parent[VALID_NAME_LENGTH];     /* who generates that: the simulation or mcstas */

  if (mcdisable_output_files) return(-1);
  if (!detector.file_handle && !strstr(detector.format.Name,"NeXus")) return(-1);
  if (strcmp(part,"header") && strstr(detector.format.Name, "no header")) return (-2);
  if (strcmp(part,"footer") && strstr(detector.format.Name, "no footer")) return (-3);

  /* initiate date and format string ======================================== */

  date_l = detector.date_l;   /* use current write time (from import) by default */
  strncpy(date, detector.date, CHAR_BUF_LENGTH);

  if (part && !strcmp(part,"footer"))
    HeadFoot = detector.format.Footer; /* footer has always detector (current) write time */
  else {
    HeadFoot = detector.format.Header; /* SIM header has simulation start time */
    if (detector.file_handle == mcsiminfo_file
        && !strcmp(detector.filename, mcsiminfo_name)) {
      date_l = mcstartdate;
      time_t t = (time_t)date_l;
      strncpy(date, ctime(&t), CHAR_BUF_LENGTH);
      if (strlen(date)) date[strlen(date)-1] = '\0';
    }
  }

  mcvalid_name(valid_parent,
    strlen(detector.filename) && detector.file_handle!=mcsiminfo_file ?
      detector.filename : detector.simulation, VALID_NAME_LENGTH);

  /* output header ========================================================== */

#ifdef USE_NEXUS
  if (strstr(detector.format.Name, "NeXus")) {
    if(mcnxinfo_header(mcnxHandle, part,
      detector.prefix,                          /* %1$s  PRE  */
      detector.instrument,                      /* %2$s  SRC  */
      strlen(detector.filename) ?
        detector.filename: detector.simulation, /* %3$s  FIL  */
      detector.format.Name,                     /* %4$s  FMT  */
      detector.date,                            /* %5$s  DAT  */
      detector.user,                            /* %6$s  USR  */
      valid_parent,                             /* %7$s  PAR  */
      detector.date_l) == NX_ERROR) {           /* %8$li DATL */
        fprintf(stderr, "Error: writing NeXus header file %s (mcinfo_header)\n",
          strlen(detector.filename) ?
            detector.filename: detector.simulation);
        /* close information data set opened by mcnxfile_section */
      return(-1);
    }
    else return(1);
  }
  else
#endif
  return(pfprintf(detector.file_handle, HeadFoot, "sssssssl",
    detector.prefix,                          /* %1$s  PRE  */
    detector.instrument,                      /* %2$s  SRC  */
    strlen(detector.filename) ?
      detector.filename: detector.simulation, /* %3$s  FIL  */
    detector.format.Name,                     /* %4$s  FMT  */
    detector.date,                            /* %5$s  DAT  */
    detector.user,                            /* %6$s  USR  */
    valid_parent,                             /* %7$s  PAR  */
    detector.date_l));                        /* %8$li DATL */
} /* mcinfo_header */

/*******************************************************************************
* mcinfo_tag: output tag/value using specific file format.
*   outputs a tag/value pair e.g. section.name=value
* RETURN: negative==error, positive=all went fine
* Used by: mcfile_section, mcinfo_instrument, mcinfo_simulation, mcinfo_data
*******************************************************************************/
static int mcinfo_tag(MCDETECTOR detector, char *section, char *name, char *value)
{
  char valid_section[VALID_NAME_LENGTH];
  char valid_name[VALID_NAME_LENGTH];
  int i;

  if (!detector.file_handle && !strstr(detector.format.Name,"NeXus")) return(-1);
  if (!strlen(detector.format.AssignTag)
   || !name || !strlen(name)
   || mcdisable_output_files) return(-1);

  mcvalid_name(valid_section, section, VALID_NAME_LENGTH);
  mcvalid_name(valid_name,    name,    VALID_NAME_LENGTH);

  /* remove quote chars in values =========================================== */
  if (strstr(detector.format.Name, "Scilab")
   || strstr(detector.format.Name, "Matlab")
   || strstr(detector.format.Name, "IDL"))
    for(i = 0; i < strlen(value); i++) {
      if (value[i] == '"' || value[i] == '\'')   value[i] = ' ';
      if (value[i] == '\n'  || value[i] == '\r') value[i] = ';';
    }

  /* output tag ============================================================= */
#ifdef USE_NEXUS
  if (strstr(detector.format.Name, "NeXus")) {
    if(mcnxinfo_tag(mcnxHandle, detector.prefix, valid_section, name, value) == NX_ERROR) {
      fprintf(stderr, "Error: writing NeXus tag file %s/%s=%s (mcinfo_tag)\n", section, name, value);
      return(-1);
    } else return(1); }
  else
#endif
  return(pfprintf(detector.file_handle, detector.format.AssignTag, "ssss",
    detector.prefix,  /* %1$s PRE */
    valid_section,    /* %2$s SEC */
    valid_name,       /* %3$s NAM */
    value));          /* %4$s VAL */
} /* mcinfo_tag */

/*******************************************************************************
* mcfile_section: output section start/end using specific file format.
*   'part' may be 'begin' or 'end' depending on section part to write
*   'type' may be e.g. 'instrument','simulation','component','data'
*   the prefix is automatically indented/un-indented
* RETURN: detector structure with updated prefix (indentation)
* Used by: mcsiminfo_init, mcsiminfo_close, mcinfo_simulation, mcdetector_write_sim
*******************************************************************************/
MCDETECTOR mcfile_section(MCDETECTOR detector, char *part, char *parent, char *section, char *type, int level)
{
  char *Section;
  char valid_section[VALID_NAME_LENGTH];
  char valid_parent[VALID_NAME_LENGTH];

  if((!detector.file_handle && !strstr(detector.format.Name,"NeXus"))
   || !section || !strlen(section)
   || mcdisable_output_files) return(detector);
  if (strcmp(part,"begin") && strstr(detector.format.Name, "no header")) return (detector);
  if (strcmp(part,"end")   && strstr(detector.format.Name, "no footer")) return (detector);

  /* initiate format string and handle prefix-- ============================= */

  if (part && !strcmp(part,"end")) Section = detector.format.EndSection;
  else                             Section = detector.format.BeginSection;

  if (!strlen(Section)) return (detector);

  mcvalid_name(valid_section, section, VALID_NAME_LENGTH);
  if (parent && strlen(parent)) mcvalid_name(valid_parent, parent, VALID_NAME_LENGTH);
  else                                strcpy(valid_parent, "root");

  if (!strcmp(part,"end") && strlen(detector.prefix) >= 2)
    detector.prefix[strlen(detector.prefix)-2]='\0'; /* un-indent prefix */

  /* output section ========================================================= */

#ifdef USE_NEXUS
  if (strstr(detector.format.Name, "NeXus")) {
    if (mcnxfile_section(mcnxHandle,part,
      detector.prefix, type, section, valid_section, parent, valid_parent, level) == NX_ERROR) {
      fprintf(stderr, "Error: writing NeXus section %s/%s=NX%s (mcfile_section)\n", parent, section, type);
      return(detector);
    }
  }
  else
#endif
  pfprintf(detector.file_handle, Section, "ssssssi",
    detector.prefix,/* %1$s  PRE  */
    type,           /* %2$s  TYP  */
    section,        /* %3$s  NAM  */
    valid_section,  /* %4$s  VNA  */
    parent,         /* %5$s  PAR  */
    valid_parent,   /* %6$s  VPA  */
    level);         /* %7$i  LVL */

  /* handle prefix++ and write section start ID ============================= */

  if (!strcmp(part,"begin"))
  {
    if (strlen(detector.prefix) < CHAR_BUF_LENGTH-2) strcat(detector.prefix,"  ");  /* indent prefix */
    if (section && strlen(section))
      mcinfo_tag(detector, section, "name",   section);
    if (parent && strlen(parent))
      mcinfo_tag(detector, section, "parent", parent);
  }

  return(detector);
} /* mcfile_section */

/*******************************************************************************
* mcinfo_instrument: output instrument tags/info. Only written to SIM file
* RETURN: negative==error, positive=all went fine
* Used by: mcsiminfo_init, mcdetector_write_sim
*******************************************************************************/
static int mcinfo_instrument(MCDETECTOR detector, char *name)
{
  char Parameters[CHAR_BUF_LENGTH] = "";
  int  i;

  if (!detector.file_handle && !strstr(detector.format.Name,"NeXus")) return(-1);
  if (!name || !strlen(name)
   || mcdisable_output_files) return(-1);

  /* create parameter string ================================================ */

  for(i = 0; i < mcnumipar; i++)
  {
    char ThisParam[VALID_NAME_LENGTH];
    if (strlen(mcinputtable[i].name) > VALID_NAME_LENGTH) break;
    snprintf(ThisParam, VALID_NAME_LENGTH, " %s(%s)", mcinputtable[i].name,
            (*mcinputtypes[mcinputtable[i].type].parminfo)
                (mcinputtable[i].name));
    strcat(Parameters, ThisParam);
    if (strlen(Parameters) >= CHAR_BUF_LENGTH-VALID_NAME_LENGTH) break;
  }

  /* output data ============================================================ */
  mcinfo_tag(detector, name, "Parameters",    Parameters);
  mcinfo_tag(detector, name, "Source",        mcinstrument_source);
  mcinfo_tag(detector, name, "Trace_enabled", mctraceenabled ? "yes" : "no");
  mcinfo_tag(detector, name, "Default_main",  mcdefaultmain ? "yes" : "no");
  mcinfo_tag(detector, name, "Embedded_runtime",
#ifdef MC_EMBEDDED_RUNTIME
         "yes"
#else
         "no"
#endif
         );

  fflush(NULL);

  return(1);
} /* mcinfo_instrument */

/*******************************************************************************
* mcinfo_simulation: output simulation tags/info
* RETURN: detector structure with updated prefix (indentation)
* Used by: mcsiminfo_init, mcdetector_write_sim
*******************************************************************************/
MCDETECTOR mcinfo_simulation(MCDETECTOR detector, char *instr)
{
  int i;
  char Parameters[CHAR_BUF_LENGTH];

  if (!detector.file_handle && !strstr(detector.format.Name,"NeXus")) return(detector);
  if (!instr || !strlen(instr)
   || mcdisable_output_files) return(detector);

  mcinfo_tag(detector, instr, "Ncount",      detector.ncount);
  mcinfo_tag(detector, instr, "Trace",       mcdotrace ? "yes" : "no");
  mcinfo_tag(detector, instr, "Gravitation", mcgravitation ? "yes" : "no");
  snprintf(Parameters, CHAR_BUF_LENGTH, "%ld", mcseed);
  mcinfo_tag(detector, instr, "Seed", Parameters);

  /* output parameter string ================================================ */
  if (!strstr(detector.format.Name, "McStas")) {
#ifdef USE_NEXUS
    /* close simulation/information */
    if (strstr(detector.format.Name, "NeXus"))
    if (NXclosedata(mcnxHandle) == NX_ERROR)
      fprintf(stderr, "Warning: NeXus: could not close data set simulation/information in %s\n", detector.filename);
#endif
    detector = mcfile_section(detector, "begin", instr, "parameters", "parameters", 3);
  }

  for(i = 0; i < mcnumipar; i++) {
    if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
      if (mcinputtable[i].par == NULL)
        strncpy(Parameters, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
      else
        (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
      if (strstr(detector.format.Name, "McStas"))
        fprintf(detector.file_handle, "%sParam: %s=%s\n", detector.prefix, mcinputtable[i].name, Parameters);
      else
        mcinfo_tag(detector, "parameters", mcinputtable[i].name, Parameters);
    }
  }
#ifdef USE_NEXUS
  if (strstr(mcformat.Name, "NeXus")) {
    /* close information SDS opened with mcfile_section(detector, "parameters") */
    if (NXclosedata(mcnxHandle) == NX_ERROR)
      fprintf(stderr, "Warning: NeXus: could not close data set simulation/parameters/information in %s\n", detector.filename);
    /* in the parameter group, create one SDS per parameter */
    mcnxfile_parameters(mcnxHandle);

    /* this can not be included in nexus-lib as it requires mcinputtypes which
      is defined only in mccode-r.c (points to function that are defined here) */


    /* re-open the simulation/information dataset */
    if (NXopendata(mcnxHandle, "information") == NX_ERROR) {
      fprintf(stderr, "Warning: NeXus: could not re-open 'information' for simulation in %s\n", detector.filename);
    }
  }
#endif
  if (!strstr(detector.format.Name, "McStas"))
    detector = mcfile_section(detector, "end", instr, "parameters", "parameters", 3);

  fflush(NULL);

  return(detector);
} /* mcinfo_simulation */

/*******************************************************************************
* mcinfo_data: output data tags/info
*              mcinfo_data(detector, NULL)        writes data info to SIM file
*              mcinfo_data(detector, "data file") writes info to data file
* Used by: mcdetector_write_sim, mcdetector_write_data
*******************************************************************************/
void mcinfo_data(MCDETECTOR detector, char *filename)
{
  char parent[CHAR_BUF_LENGTH];

  if (!detector.m || mcdisable_output_files) return;

  /* access either SIM or Data file */
  if (!filename) {
    /* use SIM file which has already be opened with instrument/simulation info */
    detector.file_handle = mcsiminfo_file;
    if (strstr(detector.format.Name, "McStas") && strlen(detector.prefix)>1) detector.prefix[0]=' ';
  } else {
    if (!detector.file_handle) return;
  }

  /* parent must be a valid name */
  mcvalid_name(parent, detector.filename, VALID_NAME_LENGTH);

  /* output data ============================================================ */
  mcinfo_tag(detector, parent, "type",       detector.type);
  mcinfo_tag(detector, parent, "Source",     mcinstrument_source);
  mcinfo_tag(detector, parent, (strstr(detector.format.Name,"McStas") ?
        "component" : "parent"),             detector.component);
  mcinfo_tag(detector, parent, "position",   detector.position);

  mcinfo_tag(detector, parent, "title",      detector.title);
  mcinfo_tag(detector, parent, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
      "Ncount" : "ratio",      detector.ncount);

  if (strlen(detector.filename)) {
    mcinfo_tag(detector, parent, "filename", detector.filename);
    mcinfo_tag(detector, parent, "format",   detector.format.Name);
  }

  mcinfo_tag(detector, parent, "statistics", detector.statistics);
  mcinfo_tag(detector, parent, strstr(detector.format.Name, "NeXus") ?
       "Signal" : "signal",                  detector.signal);
  mcinfo_tag(detector, parent, "values",     detector.values);

  if (detector.rank >= 1 || detector.rank == -1)
  {
    mcinfo_tag(detector, parent, (strstr(detector.format.Name," scan ") ?
          "xvars" : "xvar"),                 detector.xvar);
    mcinfo_tag(detector, parent, (strstr(detector.format.Name," scan ") ?
          "yvars" : "yvar"),                 detector.yvar);
    mcinfo_tag(detector, parent, "xlabel",   detector.xlabel);
    mcinfo_tag(detector, parent, "ylabel",   detector.ylabel);
    if (detector.rank > 1 || detector.rank == -1) {
      mcinfo_tag(detector, parent, "zvar",   detector.zvar);
      mcinfo_tag(detector, parent, "zlabel", detector.zlabel);
    }
  }

  mcinfo_tag(detector, parent, abs(detector.rank)==1 && strstr(detector.format.Name, "McStas") ?
                    "xlimits" : "xylimits", detector.limits);
  mcinfo_tag(detector, parent, "variables", detector.rank == -1 ?
                                            detector.zvar : detector.variables);

   if (mcDetectorCustomHeader && strlen(mcDetectorCustomHeader)) {
     mcinfo_tag(detector, parent, "custom",  mcDetectorCustomHeader);
   }

   fflush(NULL);
} /* mcinfo_data */

/*******************************************************************************
* mcdetector_import: build detector structure, merge non-lists from MPI
*                    compute basic stat, write "Detector:" line
*                    mcdetector_import(... filename=NULL ...) sets siminfo data
* RETURN:            detector structure. Invalid data if detector.p1 == NULL
*                    Invalid detector sets m=0 and filename=""
*                    Simulation data  sets m=0 and filename=mcsiminfo_name
* Used by: mcdetector_out_nD public functions, mcdetector_import_sim
*******************************************************************************/
MCDETECTOR mcdetector_import(struct mcformats_struct format,
  char *component, char *title,
  long m, long n,  long p,
  char *xlabel, char *ylabel, char *zlabel,
  char *xvar, char *yvar, char *zvar,
  double x1, double x2, double y1, double y2, double z1, double z2,
  char *filename,
  double *p0, double *p1, double *p2,
  Coords position)
{
  time_t t;       /* for detector.date */
  long   date_l;  /* date as a long number */
  char   istransposed=0;
  char   c[CHAR_BUF_LENGTH]; /* temp var for signal label */

  MCDETECTOR detector;

  /* build MCDETECTOR structure ============================================= */
  /* make sure we do not have NULL for char fields */

  /* these also apply to simfile */
  strncpy (detector.filename,  filename ? filename : "",    CHAR_BUF_LENGTH);
  /* add extension if missing */
  if (strlen(detector.filename) && !strchr(detector.filename, '.') && format.Extension)
  { /* add extension if not in file name already */
    strcat(detector.filename, ".");
    strcat(detector.filename, format.Extension);
  }
  strncpy (detector.component, component ? component : "McStas component", CHAR_BUF_LENGTH);
  snprintf(detector.simulation, CHAR_BUF_LENGTH, "%s%s%s",
      mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);
  detector.format=format;

  /* set prefix: # for McStas data files, VRML and Python */
  strcpy(detector.prefix,
      (strstr(format.Name, "McStas") && (strlen(detector.filename)) ?
      "# " : ""));

  snprintf(detector.instrument, CHAR_BUF_LENGTH, "%s (%s)", mcinstrument_name, mcinstrument_source);
  snprintf(detector.user, CHAR_BUF_LENGTH,      "%s on %s",
        getenv("USER") ? getenv("USER") : "mcstas",
        getenv("HOST") ? getenv("HOST") : "localhost");
  time(&t);         /* get current write time */
  date_l = (long)t; /* same but as a long */
  snprintf(detector.date, CHAR_BUF_LENGTH, "%s", ctime(&t));
  if (strlen(detector.date))   detector.date[strlen(detector.date)-1] = '\0'; /* remove last \n in date */
  detector.date_l = date_l;

  if (!mcget_run_num() || mcget_run_num() >= mcget_ncount())
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%g", (double)mcget_ncount());
  else
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%g/%g", (double)mcget_run_num(), (double)mcget_ncount());

  detector.p0         = p0; detector.p0_orig=p0;
  detector.p1         = p1; detector.p1_orig=p1;
  detector.p2         = p2; detector.p2_orig=p2;
  detector.file_handle= filename ? NULL : mcsiminfo_file;

  /* handle transposition */
  if (!strstr(format.Name, "NeXus")) {
    if (m<0 || n<0 || p<0)                istransposed = !istransposed;
    if (strstr(format.Name, "binary"))    istransposed = !istransposed;
    if (strstr(format.Name, "transpose")) istransposed = !istransposed;
    if (istransposed) { /* do the swap once for all */
      long i=m; m=abs(n); n=abs(i); p=abs(p);
    }
  }
  m=abs(m); n=abs(n); p=abs(p); /* make sure dimensions are positive */
  detector.istransposed = istransposed;

  /* determine detector rank (dimensionality) */
  if (!m || !n || !p || !p1) detector.rank = 4; /* invalid: exit with m=0 filename="" */
  else if (strstr(format.Name," scan ")) detector.rank=-1;  /* 1D scan: multiarray */
  else if (m*n*p == 1)       detector.rank = 0; /* 0D */
  else if (n == 1 || m == 1) detector.rank = 1; /* 1D */
  else if (p == 1)           detector.rank = 2; /* 2D */
  else                       detector.rank = 3; /* 3D */

  /* from rank, set type */
  switch (detector.rank) {
    case -1: sprintf(detector.type, "multiarray_1d(%ld)", n); p=1; break;
    case 0:  strcpy(detector.type,  "array_0d"); m=n=p=1; break;
    case 1:  sprintf(detector.type, "array_1d(%ld)", m*n*p); m *= n*p; n=p=1; break;
    case 2:  sprintf(detector.type, "array_2d(%ld, %ld)", m, n*p); n *= p; p=1; break;
    case 3:  sprintf(detector.type, "array_3d(%ld, %ld, %ld)", m, n, p); break;
    default: m=0; strcpy(detector.type, ""); strcpy(detector.filename, "");/* invalid */
  }

  detector.m    = m;
  detector.n    = n;
  detector.p    = p;

/* init default values for statistics */
  detector.intensity  = 0;
  detector.error      = 0;
  detector.events     = 0;
  detector.min        = 0;
  detector.max        = 0;
  detector.mean       = 0;
  detector.centerX    = 0;
  detector.halfwidthX = 0;
  detector.centerY    = 0;
  detector.halfwidthY = 0;

  /* these only apply to detector files ===================================== */

  snprintf(detector.position, CHAR_BUF_LENGTH, "%g %g %g", position.x, position.y, position.z);
  /* may also store actual detector orientation in the future */

  strncpy(detector.title,      title && strlen(title) ? title : component,       CHAR_BUF_LENGTH);
  strncpy(detector.xlabel,     xlabel && strlen(xlabel) ? xlabel : "X", CHAR_BUF_LENGTH); /* axis labels */
  strncpy(detector.ylabel,     ylabel && strlen(ylabel) ? ylabel : "Y", CHAR_BUF_LENGTH);
  strncpy(detector.zlabel,     zlabel && strlen(zlabel) ? zlabel : "Z", CHAR_BUF_LENGTH);
  strncpy(detector.xvar,       xvar && strlen(xvar) ? xvar :       "x", CHAR_BUF_LENGTH); /* axis variables */
  strncpy(detector.yvar,       yvar && strlen(yvar) ? yvar :       detector.xvar, CHAR_BUF_LENGTH);
  strncpy(detector.zvar,       zvar && strlen(zvar) ? zvar :       detector.yvar, CHAR_BUF_LENGTH);

  /* set "variables" as e.g. "I I_err N" */
  strcpy(c, "I ");
  if (strlen(detector.zvar))      strncpy(c, detector.zvar,32);
  else if (strlen(detector.yvar)) strncpy(c, detector.yvar,32);
  else if (strlen(detector.xvar)) strncpy(c, detector.xvar,32);

  if (detector.rank == 1)
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s %s_err N", detector.xvar, c, c);
  else
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s_err N", c, c);

  /* limits */
  detector.xmin = x1;
  detector.xmax = x2;
  detector.ymin = y1;
  detector.ymax = y2;
  detector.zmin = z1;
  detector.zmax = z2;
  if (abs(detector.rank) == 1 && strstr(format.Name, "McStas"))
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g", x1, x2);
  else if (detector.rank == 2)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g", x1, x2, y1, y2);
  else
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g %g %g", x1, x2, y1, y2, z1, z2);

  if (!m || !n || !p) {
    return(detector);
  }

  /* if MPI and nodes_nb > 1: reduce data sets when using MPI =============== */
  /* not for scan steps/multiarray (only processed by root */
#ifdef USE_MPI
  if (!strstr(format.Name," list ") && mpi_node_count > 1 && detector.rank != -1) {
    /* we save additive data: reduce everything into mpi_node_root */
    if (p0) mc_MPI_Sum(p0, m*n*p);
    if (p1) mc_MPI_Sum(p1, m*n*p);
    if (p2) mc_MPI_Sum(p2, m*n*p);
    if (!p0) {  /* additive signal must be then divided by the number of nodes */
      int i;
      for (i=0; i<m*n*p; i++) {
        p1[i] /= mpi_node_count;
        if (p2) p2[i] /= mpi_node_count;
      }
    }
  }
#endif /* USE_MPI */

  /* compute statistics and update MCDETECTOR structure ===================== */
  double sum_z   = 0;
  double min_z   = 0;
  double max_z   = 0;
  double fmon_x=0, smon_x=0, fmon_y=0, smon_y=0, mean_z=0;
  double Nsum=0;
  double P2sum=0;

  double sum_xz  = 0;
  double sum_yz  = 0;
  double sum_x   = 0;
  double sum_y   = 0;
  double sum_x2z = 0;
  double sum_y2z = 0;
  int    i,j;
  char   hasnan=0;
  char   hasinf=0;
  char   israw = (strstr(detector.format.Name," raw") != NULL);
  double *this_p1=NULL; /* new 1D McStas array [x I E N]. Freed after writing data */

  /* if McStas/PGPLOT and rank==1 we create a new m*4 data block=[x I E N] */
  if (detector.rank == 1 && strstr(detector.format.Name,"McStas")) {
    this_p1 = (double *)calloc(detector.m*detector.n*detector.p*4, sizeof(double));
    if (!this_p1)
      exit(-fprintf(stderr, "Error: Out of memory creating %li 1D McStas data set for file '%s' (mcdetector_import)\n", detector.m*detector.n*detector.p*4*sizeof(double*), detector.filename));
  }

  max_z = min_z = p1[0];

  for(j = 0; j < n*p; j++)
  {
    for(i = 0; i < m; i++)
    {
      double x,y,z;
      double N, E;
      long   index= !detector.istransposed ? i*n*p + j : i+j*m;

      if (m) x = x1 + (i + 0.5)/m*(x2 - x1); else x = 0;
      if (n && p) y = y1 + (j + 0.5)/n/p*(y2 - y1); else y = 0;
      z = p1[index];
      N = p0 ? p0[index] : 1;
      E = p2 ? p2[index] : 0;
      if (p2 && !israw) p2[index] = (*mcestimate_error_p)(p0[i],p1[i],p2[i]); /* set sigma */

      /* compute stats integrals */
      sum_xz += x*z;
      sum_yz += y*z;
      sum_x += x;
      sum_y += y;
      sum_z += z;
      sum_x2z += x*x*z;
      sum_y2z += y*y*z;
      if (z > max_z) max_z = z;
      if (z < min_z) min_z = z;

      Nsum += N;
      P2sum += E;

      if (isnan(z) || isnan(E) || isnan(N)) hasnan=1;
      if (isinf(z) || isinf(E) || isinf(N)) hasinf=1;

      if (detector.rank == 1 && this_p1 && strstr(detector.format.Name,"McStas")) {
        /* fill-in 1D McStas array [x I E N] */
        this_p1[index*4]   = x;
        this_p1[index*4+1] = z;
        this_p1[index*4+2] = p2 ? p2[index] : 0;
        this_p1[index*4+3] = N;
      }
    }
  } /* for j */

  /* compute 1st and 2nd moments */
  if (sum_z && n*m*p)
  {
    fmon_x = sum_xz/sum_z;
    fmon_y = sum_yz/sum_z;
    smon_x = sum_x2z/sum_z-fmon_x*fmon_x; smon_x = smon_x > 0 ? sqrt(smon_x) : 0;
    smon_y = sum_y2z/sum_z-fmon_y*fmon_y; smon_y = smon_y > 0 ? sqrt(smon_y) : 0;
    mean_z = sum_z/n/m/p;
  }
  /* store statistics into detector */
  detector.intensity = sum_z;
  detector.error     = (*mcestimate_error_p)(Nsum, sum_z, P2sum);
  detector.events    = Nsum;
  detector.min       = min_z;
  detector.max       = max_z;
  detector.mean      = mean_z;
  detector.centerX   = fmon_x;
  detector.halfwidthX= smon_x;
  detector.centerY   = fmon_y;
  detector.halfwidthY= smon_y;

  /* if McStas/PGPLOT and rank==1 replace p1 with new m*4 1D McStas and clear others */
  if (detector.rank == 1 && this_p1 && strstr(detector.format.Name,"McStas")) {
    detector.p1 = this_p1;
    detector.n=detector.m; detector.m  = 4;
    detector.p0 = NULL;
    detector.p2 = NULL;
    detector.istransposed = 1;
  }

  if (n*m*p > 1)
    sprintf(detector.signal, "Min=%g; Max=%g; Mean=%g;", detector.min, detector.max, detector.mean);
  else
    strcpy(detector.signal, "");
  sprintf(detector.values, "%g %g %g", detector.intensity, detector.error, detector.events);

  switch (detector.rank) {
    case 0:  strcpy(detector.statistics, ""); break;
    case 1:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g;",
      detector.centerX, detector.halfwidthX); break;
    case 2:
    case 3:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g; Y0=%g; dY=%g;",
      detector.centerX, detector.halfwidthX, detector.centerY, detector.halfwidthY);
      break;
    default: strcpy(detector.statistics, "");
  }

#ifdef USE_MPI
  /* slaves are done */
  if(mpi_node_rank != mpi_node_root) {
    return detector;
  }
#endif
  /* output "Detector:" line ================================================ */
  /* when this is a detector written by a component (not the SAVE from instrument),
     not a multiarray nor event lists */
  if (detector.rank == -1 || strstr(format.Name," list ")) return(detector);

  if (!strcmp(detector.component, mcinstrument_name)) {
    if (strlen(detector.filename))  /* we name it from its filename, or from its title */
      mcvalid_name(c, detector.filename, CHAR_BUF_LENGTH);
    else
      mcvalid_name(c, detector.title, CHAR_BUF_LENGTH);
  } else
    strncpy(c, detector.component, CHAR_BUF_LENGTH);  /* usual detectors written by components */

  printf("Detector: %s_I=%g %s_ERR=%g %s_N=%g",
         c, detector.intensity,
         c, detector.error,
         c, detector.events);
  printf(" \"%s\"\n", strlen(detector.filename) ? detector.filename : detector.component);

  if (hasnan)
    printf("WARNING: Nan detected in component %s\n", detector.component);
  if (hasinf)
    printf("WARNING: Inf detected in component %s\n", detector.component);

  /* add warning in case of low statistics or large number of bins in text format mode */
  if (detector.error > detector.intensity/4)
    printf("WARNING: file '%s': Low Statistics\n",    detector.filename);
  else if (strlen(detector.filename)) {
    if (m*n*p > 1000 && Nsum < m*n*p && Nsum)
       printf(
        "WARNING: file '%s': Low Statistics (%g events in %ldx%ldx%ld bins).\n",
        detector.filename, Nsum, m,n,p);
    if ( !strstr(format.Name, "binary")
      && (strstr(format.Name, "Scilab") || strstr(format.Name, "Matlab"))
      && (n*m*p > 10000 || m > 1000) ) printf(
        "WARNING: file '%s' (%s format)\n"
        "         Large matrices (%ldx%ldx%ld) in text mode may be\n"
        "         slow or fail at import. Prefer binary mode e.g. --format='Matlab_binary'.\n",
        detector.filename, format.Name, m,n,p);
  }

  if (mcDetectorCustomHeader && strlen(mcDetectorCustomHeader)) {
   if (strstr(detector.format.Name, "Octave") || strstr(detector.format.Name, "Matlab"))
     str_rep(mcDetectorCustomHeader, "%PRE", "%   ");
   else if (strstr(detector.format.Name, "IDL"))    str_rep(mcDetectorCustomHeader, "%PRE", ";   ");
   else if (strstr(detector.format.Name, "Scilab")) str_rep(mcDetectorCustomHeader, "%PRE", "//  ");
   else if (strstr(detector.format.Name, "McStas")) str_rep(mcDetectorCustomHeader, "%PRE", "#   ");
   else str_rep(mcDetectorCustomHeader, "%PRE", "    ");
  }

  return(detector);
} /* mcdetector_import */

/*******************************************************************************
* mcdetector_register: adds detector to the detector array, for scans and sim abstract
* RETURN:            detector array pointer
* Used by: mcdetector_out_xD
*******************************************************************************/
MCDETECTOR* mcdetector_register(MCDETECTOR detector)
{
#ifdef USE_MPI
  /* only for Master */
  if(mpi_node_rank != mpi_node_root)                      return(NULL);
#endif

  if (detector.m) {
    /* add detector entry in the mcDetectorArray */
    if (!mcDetectorArray ||
        (mcDetectorArray && mcDetectorArray_size <= mcDetectorArray_index)) {
      mcDetectorArray_size+=CHAR_BUF_LENGTH;
      if (!mcDetectorArray || !mcDetectorArray_index)
        mcDetectorArray = (MCDETECTOR*)calloc(mcDetectorArray_size, sizeof(MCDETECTOR));
      else
        mcDetectorArray = (MCDETECTOR*)realloc(mcDetectorArray, mcDetectorArray_size*sizeof(MCDETECTOR));
    }
    mcDetectorArray[mcDetectorArray_index++]=detector;
  }

  return(mcDetectorArray);

} /* mcdetector_register */

MCDETECTOR mcdetector_init(void) {
  Coords zero={0.0,0.0,0.0};
  MCDETECTOR detector=mcdetector_import(mcformat, "mcstas", NULL,
    0,0,0,            /* mnp */
    NULL, NULL, NULL, /* labels */
    NULL, NULL, NULL, /* vars */
    0,0,0,0,0,0,      /* limits */
    NULL,             /* filename */
    NULL, NULL, NULL, /* p012 */
    zero);
  strncpy(detector.filename, mcsiminfo_name, CHAR_BUF_LENGTH);
  return(detector);
}

/*******************************************************************************
* mcdetector_import_sim: build detector structure as SIM data
* RETURN:            detector structure for SIM (m=0 and filename=mcsiminfo_name).
* Used by: mcsiminfo_init, mcsiminfo_close, mcdetector_write_sim
*******************************************************************************/
MCDETECTOR mcdetector_import_sim(void) {
  MCDETECTOR detector  = mcdetector_init();
  detector.file_handle = mcsiminfo_file;
  return(detector);
}

/*******************************************************************************
* mcsiminfo_init: writes simulation structured description file (mcstas.sim)
*                 f may be NULL or e.g. stdout
*                 opens NX file if this is NeXus format
* Used by: mcsave/mcfinally from code generation (cogen), mcinfo, mcstas_main
*******************************************************************************/
static int mcsiminfo_init(FILE *f)
{

#ifdef USE_MPI
  /* only for Master */
  if(mpi_node_rank != mpi_node_root)                      return(-1);
#endif
  if (mcdisable_output_files)                             return(-2);
  if (!f && (!mcsiminfo_name || !strlen(mcsiminfo_name))) return(-3);

  /* clear list of opened files to start new save session */
  if (mcopenedfiles && strlen(mcopenedfiles) > 0) strcpy(mcopenedfiles, "");

  /* open sim file ========================================================== */

#ifdef USE_NEXUS
  /* NeXus sim info is the NeXus root file. */
  if(strstr(mcformat.Name, "NeXus")) {
    mcsiminfo_file = NULL;
    if (mcnxfile_init(mcsiminfo_name, mcformat.Extension,
      strstr(mcformat.Name, "append") || strstr(mcformat.Name, "catenate") ? "a":"w",
      &mcnxHandle) == NX_ERROR) {
      return(-fprintf(stderr, "Error: opening NeXus file %s (mcsiminfo_init)\n", mcsiminfo_name));
    } else {
      mcsiminfo_file = (FILE*)mcnxHandle; /* make it non NULL */
    }
  } else /* not NX - includes next if+else */
#endif
  {
    if (!f || (!mcsiminfo_file && f != stdout))
      mcsiminfo_file = mcnew_file(mcsiminfo_name, mcformat.Extension, "w");
    else mcsiminfo_file = f;
  }

  if(!mcsiminfo_file)
    return(-fprintf(stderr,
            "Warning: could not open simulation description file '%s' (mcsiminfo_init)\n",
            mcsiminfo_name));

  /* initialize sim file information, sets detector.file_handle=mcsiminfo_file */
  MCDETECTOR mcsiminfo = mcdetector_import_sim();

  /* flag true for McStas or NX data format */
  int  ismcstas_nx= (strstr(mcformat.Name, "McStas") || strstr(mcformat.Name, "NeXus"));

  /* start to write meta data =============================================== */

#ifdef USE_NEXUS
  if (strstr(mcformat.Name, "NeXus")) {
    /* NXentry class */
    char file_time[CHAR_BUF_LENGTH];
    snprintf(file_time, CHAR_BUF_LENGTH, "%s_%li", mcinstrument_name, mcstartdate);
    mcsiminfo = mcfile_section(mcsiminfo, "begin", "mcstas", file_time, "entry", 1);
  }
#endif
  mcinfo_header(mcsiminfo, "header");
#ifdef USE_NEXUS
  /* close information SDS opened with mcfile_section(siminfo, "mcstas (header)") */
  if (strstr(mcformat.Name, "NeXus"))
    if (NXclosedata(mcnxHandle) == NX_ERROR)
      fprintf(stderr, "Warning: NeXus: could not close data set mcstas/information (header)\n");
#endif

  /* begin-end instrument =================================================== */
  mcsiminfo = mcfile_section(mcsiminfo, "begin", mcsiminfo.simulation, mcinstrument_name, "instrument", 1);
  mcinfo_instrument(mcsiminfo, mcinstrument_name);

#ifdef USE_NEXUS
  if (strstr(mcformat.Name, "NeXus")) {
    /* close information SDS opened with mcfile_section(siminfo, "instrument") */
    if (NXclosedata(mcnxHandle) == NX_ERROR)
      fprintf(stderr, "Warning: NeXus: could not close data set instrument/information\n");
    /* insert instrument source code */
    mcnxfile_instrcode(mcnxHandle, mcinstrument_source, mcinstrument_name);
  }
#endif
  if (ismcstas_nx)
    mcsiminfo = mcfile_section(mcsiminfo, "end", mcsiminfo.simulation, mcinstrument_name, "instrument", 1);

  /* begin-end simulation =================================================== */
  mcsiminfo = mcfile_section(mcsiminfo, "begin", mcsiminfo.simulation, mcsiminfo_name, "simulation", 2);
  mcsiminfo = mcinfo_simulation(mcsiminfo, mcsiminfo_name);

#ifdef USE_NEXUS
  /* close information SDS opened with mcfile_section(siminfo,"simulation") */
  if (strstr(mcformat.Name, "NeXus"))
    if (NXclosedata(mcnxHandle) == NX_ERROR)
      fprintf(stderr, "Warning: NeXus: could not close data set simulation/information\n");
#endif
  if (ismcstas_nx)
    mcsiminfo = mcfile_section(mcsiminfo, "end", mcsiminfo.simulation, mcsiminfo_name, "simulation", 2);

  return(1);
} /* mcsiminfo_init */

/*******************************************************************************
* mcdetector_write_content: write mcstas.dat, which has integrated intensities for all monitors
* Used by: mcsiminfo_close
*******************************************************************************/
int mcdetector_write_content(MCDETECTOR *DetectorArray, long DetectorArray_index)
{
#ifdef USE_MPI
  /* only for Master */
  if(mpi_node_rank != mpi_node_root)                      return(-1);
#endif
  if (mcdisable_output_files)                             return(-2);
  if (!mcsiminfo_name || !strlen(mcsiminfo_name))         return(-3);
  if (!DetectorArray_index)                               return(-4);

  /* build p1 array from all detector integrated counts */
  double *this_p1 = (double *)calloc(DetectorArray_index*2, sizeof(double));
  int i;

  if (this_p1 && DetectorArray) {
    char   *labels=NULL;
    long    labels_size=0;
    long    index=0;
    for (i=0; i < DetectorArray_index; i++) {
        char mem[CHAR_BUF_LENGTH];
        char valid[CHAR_BUF_LENGTH];
        /* store I I_err */
        this_p1[index++] = DetectorArray[i].intensity;
        this_p1[index++] = DetectorArray[i].error;
        /* add corresponding label */
        mcvalid_name(valid, DetectorArray[i].filename, CHAR_BUF_LENGTH);
        sprintf(mem, "%s_I %s_Err ", valid, valid);
        if (!labels ||
          (labels && labels_size <= strlen(labels)+strlen(mem))) {
          labels_size+=CHAR_BUF_LENGTH;
          if (!labels || !strlen(labels))
            labels = calloc(1, labels_size);
          else
            labels = realloc(labels, labels_size);
        }
        strcat(labels, " ");
        strcat(labels, mem);
    } /* for */

    struct mcformats_struct format=mcformat;
    strcat(format.Name, " scan step");
    Coords zero={0.0,0.0,0.0};

    /* now create detector and write 'abstract' file */
    MCDETECTOR detector = mcdetector_import(format,
      mcinstrument_name, "Monitor integrated counts",
      index, 1, 1,
      "Monitors", "I", "Integrated Signal",
      "Index", labels, "I",
      0, index-1, 0, 0, 0, 0, "content",  /* use name from OpenOffice content.xml description file */
      NULL, this_p1, NULL, zero);

    mcdetector_write_data(detector);

    free(labels); labels=NULL; labels_size=0;

    /* free DETECTOR array */
    free(this_p1); this_p1=NULL;
    free(mcDetectorArray); mcDetectorArray=NULL;
  } /* if this_p1 */

} /* mcdetector_write_content */


/*******************************************************************************
* mcsiminfo_close: close simulation file (mcstas.sim) and write mcstas.dat
* Used by: mcsave/mcfinally from code generation (cogen), mcinfo, mcstas_main
*******************************************************************************/
void mcsiminfo_close(void)
{
  int i;
#ifdef USE_MPI
  if(mpi_node_rank != mpi_node_root) return;
#endif
  if (mcdisable_output_files || !mcsiminfo_file) return;

  int  ismcstas_nx  = (strstr(mcformat.Name, "McStas") || strstr(mcformat.Name, "NeXus"));

  mcdetector_write_content(mcDetectorArray, mcDetectorArray_index);

  /* initialize sim file information, sets detector.file_handle=mcsiminfo_file */
  MCDETECTOR mcsiminfo = mcdetector_import_sim();

  /* close those sections which were opened in mcsiminfo_init =============== */
  if (!ismcstas_nx) {
    mcfile_section(mcsiminfo, "end", mcsiminfo.simulation, mcsiminfo_name, "simulation", 2);
    mcfile_section(mcsiminfo, "end", mcsiminfo.simulation, mcinstrument_name, "instrument", 1);
  }
#ifdef USE_NEXUS
  /* re-open the simulation/information dataset */
  if (NXopendata(mcnxHandle, "information") == NX_ERROR)
    fprintf(stderr, "Warning: NeXus: could not re-open 'information' for mcstas NXentry\n");
#endif
  mcinfo_header(mcsiminfo, "footer");
#ifdef USE_NEXUS
  /* close information SDS opened with mcfile_section(siminfo, "mcstas(footer)" */
  if (strstr(mcformat.Name, "NeXus")) {
    if (NXclosedata(mcnxHandle) == NX_ERROR)
      fprintf(stderr, "Warning: NeXus: could not close data set mcstas/information (footer)\n");
    mcfile_section(mcsiminfo, "end", "mcstas", mcinstrument_name, "entry", 1);
  }
#endif

  /* close sim file ========================================================= */
#ifdef USE_NEXUS
  if (strstr(mcformat.Name, "NeXus")) mcnxfile_close(&mcnxHandle);
  else
#endif
  if (mcsiminfo_file != stdout && mcsiminfo_file) {
    fclose(mcsiminfo_file);
  }
  mcsiminfo_file = NULL;

} /* mcsiminfo_close */

/*******************************************************************************
* mcfile_data: output a single data block using specific file format, e.g. "part=[ array ]"
*   'part' can be 'data','errors','ncount'
* Used by: mcdetector_write_sim, mcdetector_write_data
*******************************************************************************/
int mcfile_data(MCDETECTOR detector, char *part)
{
  double *this_p1=NULL;
  char   *Begin  =NULL;
  char   *End    =NULL;
  char    valid_xlabel[VALID_NAME_LENGTH];
  char    valid_ylabel[VALID_NAME_LENGTH];
  char    valid_zlabel[VALID_NAME_LENGTH];
  char    valid_parent[VALID_NAME_LENGTH];

  if (strstr(part,"data"))
  { this_p1=detector.p1; Begin = detector.format.BeginData; End = detector.format.EndData; }
  if (strstr(part,"errors"))
  { this_p1=detector.p2;Begin = detector.format.BeginErrors; End = detector.format.EndErrors; }
  if (strstr(part,"ncount") || strstr(part,"events"))
  { this_p1=detector.p0;Begin = detector.format.BeginNcount; End = detector.format.EndNcount; }

  if (!this_p1 || mcdisable_output_files) return(-1);
  if (!detector.file_handle && !strstr(detector.format.Name,"NeXus")) return(-1);
  if (!detector.rank) return(-1);


  mcvalid_name(valid_xlabel, detector.xlabel, VALID_NAME_LENGTH);
  mcvalid_name(valid_ylabel, detector.ylabel, VALID_NAME_LENGTH);
  mcvalid_name(valid_zlabel, detector.zlabel, VALID_NAME_LENGTH);
  mcvalid_name(valid_parent, detector.filename, VALID_NAME_LENGTH);

  /* output Begin =========================================================== */
  if (!strstr(detector.format.Name,"NeXus")
   && (!strstr(detector.format.Name,"binary") || strstr(part,"empty array"))
   && !strstr(detector.format.Name,"no header"))
  pfprintf(detector.file_handle, Begin, "sssssssssssssiiigggggg",
      detector.prefix,       /* %1$s   PRE  */
      valid_parent,          /* %2$s   PAR  */
      detector.filename,     /* %3$s   FIL  */
      detector.xlabel,       /* %4$s   XLA  */
      valid_xlabel,          /* %5$s   XVL  */
      detector.ylabel,       /* %6$s   YLA  */
      valid_ylabel,          /* %7$s   YVL  */
      detector.zlabel,       /* %8$s   ZLA  */
      valid_zlabel,          /* %9$s   ZVL  */
      detector.title,        /* %10$s  TITL */
      detector.xvar,         /* %11$s  XVAR */
      detector.yvar,         /* %12$s  YVAR */
      detector.zvar,         /* %13$s  ZVAR */
      detector.m,            /* %14$i  MDIM */
      detector.n,            /* %15$i  NDIM */
      detector.p,            /* %16$i  PDIM */
      detector.xmin,         /* %17$g  XMIN */
      detector.xmax,         /* %18$g  XMAX */
      detector.ymin,         /* %19$g  YMIN */
      detector.ymax,         /* %20$g  YMAX */
      detector.zmin,         /* %21$g  ZMIN */
      detector.zmax);        /* %22$g  ZMAX */

  /* output array */
  if (strstr(part,"empty array")) {
    /* skip array output: set as empty */

    /* special case of IDL: can not have empty vectors. Init to 'external' */
    if (strstr(detector.format.Name, "IDL")) fprintf(detector.file_handle, "'external'");
  } else {
#ifdef USE_NEXUS
    /* array NeXus ========================================================== */
    if (strstr(detector.format.Name,"NeXus")) {
      if(mcnxfile_data(mcnxHandle, detector, part,
        valid_parent, valid_xlabel, valid_ylabel, valid_zlabel) == NX_ERROR) {
        fprintf(stderr, "Error: writing NeXus data %s/%s (mcfile_data)\n", valid_parent, detector.filename);
        return(-1);
      }
      return(1);
    } /* end if NeXus */
#endif
    /* array non binary (text) ============================================== */
    if (!strstr(detector.format.Name, "binary")
     && !strstr(detector.format.Name, "float") && !strstr(detector.format.Name, "double"))
    {
      char eol_char[3];
      int isIDL    = (strstr(detector.format.Name, "IDL") != NULL);
      int isPython = (strstr(detector.format.Name, "Python") != NULL);
      int i,j;
      if (isIDL) strcpy(eol_char,"$\n");
      else if (isPython) strcpy(eol_char,"\\\n");
      else strcpy(eol_char,"\n");

      for(j = 0; j < detector.n*detector.p; j++)  /* loop on rows (y) */
      {
        for(i = 0; i < detector.m; i++)  /* write all columns (x) */
        {
          long index   = !detector.istransposed ? i*detector.n*detector.p + j : i+j*detector.m;
          double value = this_p1[index];
          fprintf(detector.file_handle, "%.12g", value);
          fprintf(detector.file_handle, "%s",
            (isIDL || isPython) && ((i+1)*(j+1) < detector.m*detector.n*detector.p) ?
            "," : " ");
        }
        fprintf(detector.file_handle, "%s", eol_char);
      } /* end 2 loops if not Binary */
    } /* end if !Binary */

    /* array binary double =================================================== */
    if (strstr(detector.format.Name, "double"))
    {
      long count=0;
      count = fwrite(this_p1, sizeof(double), detector.m*detector.n*detector.p, detector.file_handle);
      if (count != detector.m*detector.n*detector.p)
        fprintf(stderr, "Error: writing double binary file '%s' (%li instead of %li, mcfile_data)\n", detector.filename,count, detector.m*detector.n*detector.p);
    } /* end if Binary double */

    /* array binary float =================================================== */
    if (strstr(detector.format.Name, "binary") || strstr(detector.format.Name, "float"))
    {
      float *s;
      s = (float*)malloc(detector.m*detector.n*detector.p*sizeof(float));
      if (s)
      {
        long    i, count=0;
        for (i=0; i<detector.m*detector.n*detector.p; i++) {
          s[i] = (float)this_p1[i]; }
        count = fwrite(s, sizeof(float), detector.m*detector.n*detector.p, detector.file_handle);
        if (count != detector.m*detector.n*detector.p)
          fprintf(stderr, "Error writing float binary file '%s' (%li instead of %li, mcfile_data)\n",
            detector.filename,count, detector.m*detector.n*detector.p);
        free(s);
      } else
      fprintf(stderr, "Error: Out of memory for writing %li float binary file '%s' (mcfile_data)\n",
        detector.m*detector.n*detector.p*sizeof(float), detector.filename);
    } /* end if Binary float */
  } /* end if not empty array */

  /* output End ============================================================= */
  if (!strstr(detector.format.Name,"NeXus")
   && (!strstr(detector.format.Name,"binary") || strstr(part,"empty array"))
   && !strstr(detector.format.Name,"no footer"))
  pfprintf(detector.file_handle, End, "sssssssssssssiiigggggg",
      detector.prefix,       /* %1$s   PRE  */
      valid_parent,          /* %2$s   PAR  */
      detector.filename,     /* %3$s   FIL  */
      detector.xlabel,       /* %4$s   XLA  */
      valid_xlabel,          /* %5$s   XVL  */
      detector.ylabel,       /* %6$s   YLA  */
      valid_ylabel,          /* %7$s   YVL  */
      detector.zlabel,       /* %8$s   ZLA  */
      valid_zlabel,          /* %9$s   ZVL  */
      detector.title,        /* %10$s  TITL */
      detector.xvar,         /* %11$s  XVAR */
      detector.yvar,         /* %12$s  YVAR */
      detector.zvar,         /* %13$s  ZVAR */
      detector.m,            /* %14$i  MDIM */
      detector.n,            /* %15$i  NDIM */
      detector.p,            /* %16$i  PDIM */
      detector.xmin,         /* %17$g  XMIN */
      detector.xmax,         /* %18$g  XMAX */
      detector.ymin,         /* %19$g  YMIN */
      detector.ymax,         /* %20$g  YMAX */
      detector.zmin,         /* %21$g  ZMIN */
      detector.zmax);        /* %22$g  ZMAX */

  return(2);
} /* mcfile_data */

/*******************************************************************************
* mcdetector_write_sim: write information to sim file
*******************************************************************************/
MCDETECTOR mcdetector_write_sim(MCDETECTOR detector)
{
  /* MPI: only write sim by Master ========================================== */
#ifdef USE_MPI
  if(mpi_node_rank != mpi_node_root) return(detector);
#endif
  /* skip invalid detectors */
  if (!detector.m || mcdisable_output_files) return(detector);

  /* sim file has been initialized when starting simulation and when calling
   * mcsave ; this defines mcsiminfo_file as the SIM file handle
   * and calls:
   *    mcinfo_header
   *    mcfile_section begin instrument
   *    mcinfo_instrument
   *    mcfile_section end instrument
   *    mcfile_section begin simulation
   *    mcinfo_simulation
   *    mcfile_section end simulation
   *    sets "mcopenedfiles" list to empty string
   * When using NeXus format, this information is stored into mcnxHandle NX
   *
   * The sim file (and mcnxHandle) is closed when ending mcsave
   */

  MCDETECTOR simfile = mcdetector_import_sim(); /* store reference structure to SIM file */

  /* add component and data section to SIM file */
  if (!strstr(detector.format.Name,"NeXus"))
    simfile = mcfile_section(simfile, "begin", detector.simulation, detector.component, "component", 3);
  simfile = mcfile_section(simfile, "begin", detector.component, detector.filename, "data", 4);

  /* handle indentation for simfile */
  char prefix[CHAR_BUF_LENGTH];
  strncpy(prefix, detector.prefix, CHAR_BUF_LENGTH);
  strncpy(detector.prefix, simfile.prefix, CHAR_BUF_LENGTH);

  /* WRITE detector information to sim file ================================= */
  mcinfo_data(detector, NULL);
#ifdef USE_NEXUS
  /* close information SDS opened with mcfile_section(siminfo,"data") */
  if (strstr(mcformat.Name, "NeXus"))
    if (NXclosedata(mcnxHandle) == NX_ERROR)
      fprintf(stderr, "Warning: NeXus: could not close data set %s/information\n",
        detector.filename);
#endif
  /* WRITE data link to SIM file */
  if (!strstr(detector.format.Name,"McStas") && !mcsingle_file) {
    FILE *file=detector.file_handle;
    detector.file_handle = mcsiminfo_file;
    mcfile_data(detector, "data (empty array)");
    if (detector.p2) mcfile_data(detector, "errors (empty array)");
    if (detector.p0) mcfile_data(detector, "events (empty array)");
    detector.file_handle = file;
  }

  /* close component and data sections in SIM file */
  simfile = mcfile_section(simfile, "end", detector.component, detector.filename, "data", 4);
  if (!strstr(detector.format.Name,"NeXus"))
    simfile = mcfile_section(simfile, "end", detector.simulation, detector.component, "component", 3);

  strncpy(detector.prefix, prefix, CHAR_BUF_LENGTH);

  return(detector);
} /* mcdetector_write_sim */

/*******************************************************************************
* mcdetector_write_data: write information to data file and catenate lists
*                        this is where we open/close data files
*                        also free detector.p1=this_p1 field which may hold
*                        1D McStas [x I Ierr N] array which is set in mcdetector_import
*******************************************************************************/
MCDETECTOR mcdetector_write_data(MCDETECTOR detector)
{
  /* skip if 0D or no filename or no data (only stored in sim file) */
  if (!detector.rank || !strlen(detector.filename)
   || !detector.m) return(detector);

#ifdef USE_MPI
  /* only by MASTER for non lists (MPI reduce has been done in detector_import) */
  if (!strstr(detector.format.Name," list ") && mpi_node_rank != mpi_node_root) {
    if (detector.rank == 1 && detector.p1 && strstr(detector.format.Name,"McStas")) {
      free(detector.p1);  /* 'this_p1' allocated in mcdetector_import for 1D McStas data sets [x I E N] */
      detector.p0 = detector.p0_orig;
      detector.p1 = detector.p1_orig;
      detector.p2 = detector.p2_orig;
      detector.m = detector.n; detector.n=1; detector.istransposed=0;
    }
    return(detector);
  }
#endif

  /* OPEN data file (possibly appending if already opened) ================== */
  if (mcformat_data.Name && strlen(mcformat_data.Name) && !mcsingle_file)
    detector.format = mcformat_data;

  if (!strstr(detector.format.Name, "NeXus")) {
    /* explicitely open data file if not NeXus format */
    char mode[10];

    /* open or re-open file for write/append */
    strcpy(mode,
           (strstr(detector.format.Name, "no header")
            || strstr(detector.format.Name, "append") || strstr(detector.format.Name, "catenate")
            || strstr(mcopenedfiles, detector.filename) ?
           "a" : "w"));
    if (strstr(detector.format.Name, "binary")) strcat(mode, "b");
    detector.file_handle = mcnew_file(!mcsingle_file ?
      detector.filename : mcsiminfo_name,
      detector.format.Extension, mode);
  } else {
    /* NeXus file mcnxHandle has been opened in mcsiminfo_init */
    /* but we must open the proper Data Set group and write information */
    mcfile_section(detector, "begin", detector.component, detector.filename, "data", 4);
  }

  /* WRITE data header (except when in 'no header') */
  if (!strstr(detector.format.Name, "no header")) {
    if (!strstr(detector.format.Name, "binary") && !mcascii_only)  {
      /* skip in data-only mode or binary */
      mcinfo_header(detector, "header");
      mcinfo_simulation(detector, mcsiminfo_name);
      mcinfo_data(detector, detector.filename);
      /* mcinfo_simulation(detector, detector.filename); */
    }
  }
#ifdef USE_NEXUS
  /* close information SDS opened with mcfile_section(detector,"data") */
  if (strstr(mcformat.Name, "NeXus"))
    mcinfo_header(detector, "footer");
    if (NXclosedata(mcnxHandle) == NX_ERROR)
      fprintf(stderr, "Warning: NeXus: could not close data set %s/information (footer)\n",
        detector.filename);
  /* group remains opened for data to be written */
#endif

  /* case: list of events =================================================== */
  if (strstr(detector.format.Name," list ")) {

#ifdef USE_MPI
    if (mpi_node_rank != mpi_node_root) {
      /* MPI slave: all slaves send their data to master */
      /* m, n, p must be sent too, since all slaves do not have the same number of events */
      int mnp[3]={detector.m,detector.n,detector.p};

      if (mc_MPI_Send(mnp, 3, MPI_INT, mpi_node_root)!= MPI_SUCCESS)
        fprintf(stderr, "Warning: proc %i to master: MPI_Send mnp list error (mcdetector_write_data)", mpi_node_rank);
      if (!detector.p1
       || mc_MPI_Send(detector.p1, abs(mnp[0]*mnp[1]*mnp[2]), MPI_DOUBLE, mpi_node_root) != MPI_SUCCESS)
        fprintf(stderr, "Warning: proc %i to master: MPI_Send p1 list error (mcdetector_write_data)", mpi_node_rank);
      /* slaves are done */
      return (detector);
    }
#endif

    /* master node or non-MPI list: store data block */
    int     master_mnp[3]={detector.m,detector.n,detector.p};

#ifdef USE_MPI
    int     node_i=0;
    /* MPI master: receive data from slaves */
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      double *this_p1=NULL;                               /* buffer to hold the list to save */
      int     mnp[3]={detector.m,detector.m,detector.p};  /* size of this buffer */
      if (node_i != mpi_node_root) { /* get data from slaves */
        if (mc_MPI_Recv(mnp, 3, MPI_INT, node_i) != MPI_SUCCESS)
          fprintf(stderr, "Warning: master from proc %i: "
            "MPI_Recv mnp list error (mcdetector_write_data)", node_i);
        this_p1 = (double *)calloc(abs(mnp[0]*mnp[1]*mnp[2]), sizeof(double));
        if (!this_p1 || mc_MPI_Recv(this_p1, abs(mnp[0]*mnp[1]*mnp[2]), MPI_DOUBLE, node_i)!= MPI_SUCCESS)
          fprintf(stderr, "Warning: master from proc %i: "
            "MPI_Recv p1 list error (mcdetector_write_data)", node_i);
        detector.p1 = this_p1;
        detector.m  = mnp[0]; detector.n  = mnp[1]; detector.p  = mnp[2];
      } else
#endif
      { /* MASTER/single node use its own detector structure */
        detector.p1 = detector.p1_orig;
        detector.m  = master_mnp[0]; detector.n  = master_mnp[1]; detector.p  = master_mnp[2];
      }

      /* WRITE list data (will be catenated in case of MPI as file is already opened) */
      mcfile_data(detector, "data");

#ifdef USE_MPI
      /* free temporary MPI block used for Recv from slaves */
      if (this_p1 && node_i != mpi_node_root)
        free(this_p1);
      detector.p0 = detector.p0_orig;
      detector.p1 = detector.p1_orig;
      detector.p2 = detector.p2_orig;
      detector.m  = master_mnp[0]; detector.n  = master_mnp[1]; detector.p  = master_mnp[2];

    } /* for node_i: end loop on nodes */
#endif

  } /* end list case */
  else
  {
  /* case: normal detector ================================================== */

    /* WRITE data */
    mcfile_data(detector, "data");

    /* write errors (not for lists) */
    if (detector.p2) mcfile_data(detector, "errors");

    /* write events (not for lists) */
    if (detector.p0) mcfile_data(detector, "events");

    if (detector.rank == 1 && detector.p1 && strstr(detector.format.Name,"McStas")) {
      free(detector.p1);  /* 'this_p1' allocated in mcdetector_write_sim for 1D McStas data sets [x I E N] */
      detector.p0 = detector.p0_orig;
      detector.p1 = detector.p1_orig;
      detector.p2 = detector.p2_orig;
      detector.m = detector.n; detector.n=1; detector.istransposed=0;
    }
  } /* end normal detector case */

  /* WRITE data footer (except when 'no footer') */
  if (!strstr(detector.format.Name, "no footer") && !strstr(mcformat.Name, "NeXus")) {
    /* skip in data-only mode or binary */
    if (!strstr(detector.format.Name, "binary") && !mcascii_only)
      mcinfo_header(detector, "footer");
  }

  /* close data file and free memory */
  if (mcDetectorCustomHeader && strlen(mcDetectorCustomHeader)) {
    free(mcDetectorCustomHeader); mcDetectorCustomHeader=NULL; }

  if (!strstr(mcformat.Name, "NeXus") && detector.file_handle != mcsiminfo_file
   && !mcdisable_output_files && detector.file_handle)
    fclose(detector.file_handle);
  else if (strstr(mcformat.Name, "NeXus")) {
    /* close data set group */
    mcfile_section(detector, "end", detector.component, detector.filename, "data", 4);
  }

  return(detector);
} /* mcdetector_write_data */

/*******************************************************************************
* mcdetector_out_0D: wrapper for 0D (single value).
*******************************************************************************/
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2,
                         char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI reduce except for lists) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : "McStas data"),
    1, 1, 1,
    "I", "", "",
    "I", "", "",
    0, 0, 0, 0, 0, 0, "",
    &p0, &p1, &p2, posa); /* write Detector: line */

  /* write detector to simulation file (incl custom header if any) */
  detector = mcdetector_write_sim(detector);
  mcdetector_register(detector);
  return(detector);
}

/*******************************************************************************
* mcdetector_out_1D: wrapper for 1D.
*******************************************************************************/
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
        char *xvar, double x1, double x2,
        int n,
        double *p0, double *p1, double *p2, char *f,
        char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI_Reduce except for lists) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : "McStas 1D data"),
    n, 1, 1,
    xl, yl, (n > 1 ? "Signal per bin" : " Signal"),
    xvar, "(I,I_err)", "I",
    x1, x2, 0, 0, 0, 0, f,
    p0, p1, p2, posa); /* write Detector: line */

  /* write detector to simulation and data file (incl custom header if any) */
  detector = mcdetector_write_sim(detector);
  detector = mcdetector_write_data(detector);  /* will also merge lists */
  mcdetector_register(detector);
  return(detector);
}

/*******************************************************************************
* mcdetector_out_2D: wrapper for 2D.
*******************************************************************************/
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2,
                  int m, int n,
                  double *p0, double *p1, double *p2, char *f,
                  char *c, Coords posa)
{
  char xvar[CHAR_BUF_LENGTH];
  char yvar[CHAR_BUF_LENGTH];

  if (xl && strlen(xl)) { strncpy(xvar, xl, CHAR_BUF_LENGTH); xvar[2]='\0'; }
  else strcpy(xvar, "x");
  if (yl && strlen(yl)) { strncpy(yvar, yl, CHAR_BUF_LENGTH); yvar[2]='\0'; }
  else strcpy(yvar, "y");

  /* is it in fact a 1D call ? */
  if (m == 1)      return(mcdetector_out_1D(
                    t, yl, "I", yvar, y1, y2, n, p0, p1, p2, f, c, posa));
  else if (n == 1) return(mcdetector_out_1D(
                    t, xl, "I", xvar, x1, x2, m, p0, p1, p2, f, c, posa));

  /* import and perform basic detector analysis (and handle MPI_Reduce except for lists) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : "McStas 2D data"),
    m, n, 1,
    xl, yl, "Signal per bin",
    xvar, yvar, "I",
    x1, x2, y1, y2, 0, 0, f,
    p0, p1, p2, posa); /* write Detector: line */

  /* write detector to simulation and data file (incl custom header if any) */
  detector = mcdetector_write_sim(detector);
  detector = mcdetector_write_data(detector); /* will also merge lists */
  mcdetector_register(detector);
  return(detector);
}

/*******************************************************************************
* mcdetector_out_3D: wrapper for 3D.
*   exported as a large 2D array, but the 3 dims are given in the header
*******************************************************************************/
MCDETECTOR mcdetector_out_3D(char *t, char *xl, char *yl, char *zl,
      char *xvar, char *yvar, char *zvar,
      double x1, double x2, double y1, double y2, double z1, double z2,
      int m, int n, int p,
      double *p0, double *p1, double *p2, char *f,
      char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI_Reduce except for lists) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : "McStas 3D data"),
    m, n, p,
    xl, yl, zl,
    xvar, yvar, zvar,
    x1, x2, y1, y2, z1, z2, f,
    p0, p1, p2, posa); /* write Detector: line */

  /* write detector to simulation and data file (incl custom header if any) */
  detector = mcdetector_write_sim(detector);
  detector = mcdetector_write_data(detector); /* will also merge lists */
  mcdetector_register(detector);
  return(detector);
}

/*******************************************************************************
* mcuse_format_header: Replaces aliases names in format fields (header part)
*******************************************************************************/
char *mcuse_format_header(char *format_const)
{ /* Header Footer */
  char *format=NULL;

  if (!format_const) return(NULL);
  format = (char *)malloc(strlen(format_const)+1);
  if (!format) exit(fprintf(stderr, "Error: insufficient memory (mcuse_format_header)\n"));
  strcpy(format, format_const);
  str_rep(format, "%PRE", "%1$s"); /* prefix */
  str_rep(format, "%SRC", "%2$s"); /* name of instrument source file */
  str_rep(format, "%FIL", "%3$s"); /* output file name (data) */
  str_rep(format, "%FMT", "%4$s"); /* format name */
  str_rep(format, "%DATL","%8$li");/* Time elapsed since Jan 1st 1970, in seconds */
  str_rep(format, "%DAT", "%5$s"); /* Date as a string */
  str_rep(format, "%USR", "%6$s"); /* User/machine name */
  str_rep(format, "%PAR", "%7$s"); /* Parent name (root/mcstas) valid_parent */
  str_rep(format, "%VPA", "%7$s");
  return(format);
}

/*******************************************************************************
* mcuse_format_tag: Replaces aliases names in format fields (tag part)
*******************************************************************************/
char *mcuse_format_tag(char *format_const)
{ /* AssignTag */
  char *format=NULL;

  if (!format_const) return(NULL);
  format = (char *)malloc(strlen(format_const)+1);
  if (!format) exit(fprintf(stderr, "Error: insufficient memory (mcuse_format_tag)\n"));
  strcpy(format, format_const);
  str_rep(format, "%PRE", "%1$s"); /* prefix */
  str_rep(format, "%SEC", "%2$s"); /* section/parent name valid_parent/section */
  str_rep(format, "%PAR", "%2$s");
  str_rep(format, "%VPA", "%2$s");
  str_rep(format, "%NAM", "%3$s"); /* name of field valid_name */
  str_rep(format, "%VNA", "%3$s");
  str_rep(format, "%VAL", "%4$s"); /* value of field */
  return(format);
}

/*******************************************************************************
* mcuse_format_section: Replaces aliases names in format fields (section part)
*******************************************************************************/
char *mcuse_format_section(char *format_const)
{ /* BeginSection EndSection */
  char *format=NULL;

  if (!format_const) return(NULL);
  format = (char *)malloc(strlen(format_const)+1);
  if (!format) exit(fprintf(stderr, "Error: insufficient memory (mcuse_format_section)\n"));
  strcpy(format, format_const);
  str_rep(format, "%PRE", "%1$s"); /* prefix */
  str_rep(format, "%TYP", "%2$s"); /* type/class of section */
  str_rep(format, "%NAM", "%3$s"); /* name of section */
  str_rep(format, "%SEC", "%3$s");
  str_rep(format, "%VNA", "%4$s"); /* valid name (letters/digits) of section */
  str_rep(format, "%PAR", "%5$s"); /* parent name (simulation) */
  str_rep(format, "%VPA", "%6$s"); /* valid parent name (letters/digits) */
  str_rep(format, "%LVL", "%7$i"); /* level index */
  return(format);
}

/*******************************************************************************
* mcuse_format_data: Replaces aliases names in format fields (data part)
*******************************************************************************/
char *mcuse_format_data(char *format_const)
{ /* BeginData EndData BeginErrors EndErrors BeginNcount EndNcount */
  char *format=NULL;

  if (!format_const) return(NULL);
  format = (char *)malloc(strlen(format_const)+1);
  if (!format) exit(fprintf(stderr, "Error: insufficient memory (mcuse_format_data)\n"));
  strcpy(format, format_const);
  str_rep(format, "%PRE", "%1$s"); /* prefix */
  str_rep(format, "%PAR", "%2$s"); /* parent name (component instance name) valid_parent */
  str_rep(format, "%VPA", "%2$s");
  str_rep(format, "%FIL", "%3$s"); /* output file name (data) */
  str_rep(format, "%XLA", "%4$s"); /* x axis label */
  str_rep(format, "%XVL", "%5$s"); /* x axis valid label (letters/digits) */
  str_rep(format, "%YLA", "%6$s"); /* y axis label */
  str_rep(format, "%YVL", "%7$s"); /* y axis valid label (letters/digits) */
  str_rep(format, "%ZLA", "%8$s"); /* z axis label */
  str_rep(format, "%ZVL", "%9$s"); /* z axis valid label (letters/digits) */
  str_rep(format, "%TITL", "%10$s"); /* data title */
  str_rep(format, "%XVAR", "%11$s"); /* x variables */
  str_rep(format, "%YVAR", "%12$s"); /* y variables */
  str_rep(format, "%ZVAR", "%13$s"); /* z variables */
  str_rep(format, "%MDIM", "%14$i"); /* m dimension of x axis */
  str_rep(format, "%NDIM", "%15$i"); /* n dimension of y axis */
  str_rep(format, "%PDIM", "%16$i"); /* p dimension of z axis */
  str_rep(format, "%XMIN", "%17$g"); /* x min axis value (m bins) */
  str_rep(format, "%XMAX", "%18$g"); /* x max axis value (m bins) */
  str_rep(format, "%YMIN", "%19$g"); /* y min axis value (n bins) */
  str_rep(format, "%YMAX", "%20$g"); /* y max axis value (n bins) */
  str_rep(format, "%ZMIN", "%21$g"); /* z min axis value (usually min of signal, p bins) */
  str_rep(format, "%ZMAX", "%22$g"); /* z max axis value (usually max of signal, p bins) */
  return(format);
}

/*******************************************************************************
* mcuse_format: selects an output format for sim and data files. returns format
*******************************************************************************/
struct mcformats_struct mcuse_format(char *format)
{
  int i,j;
  int i_format=-1;
  char *tmp;
  char low_format[256];
  struct mcformats_struct usedformat;

  usedformat.Name = NULL;
  /* get the format to lower case */
  if (!format) exit(fprintf(stderr, "Error: Invalid NULL format. Exiting (mcuse_format)\n"));
  strcpy(low_format, format);
  for (i=0; i<strlen(low_format); i++) low_format[i]=tolower(format[i]);
  /* handle format aliases */
  if (!strcmp(low_format, "pgplot")) strcpy(low_format, "mcstas");
  if (!strcmp(low_format, "hdf"))    strcpy(low_format, "nexus");
#ifndef USE_NEXUS
  if (!strcmp(low_format, "nexus"))
    fprintf(stderr, "WARNING: to enable NeXus format you must have the NeXus libs installed.\n"
                    "         You should then re-compile with the -DUSE_NEXUS define.\n");
#endif

  tmp = (char *)malloc(256);
  if(!tmp) exit(fprintf(stderr, "Error: insufficient memory (mcuse_format)\n"));

  /* look for a specific format in mcformats.Name table */
  for (i=0; i < mcNUMFORMATS; i++)
  {
    strcpy(tmp, mcformats[i].Name);
    for (j=0; j<strlen(tmp); j++) tmp[j] = tolower(tmp[j]);
    if (strstr(low_format, tmp))  i_format = i;
  }
  if (i_format < 0)
  {
    i_format = 0; /* default format is #0 McStas */
    fprintf(stderr, "Warning: unknown output format '%s'. Using default (%s, mcuse_format).\n", format, mcformats[i_format].Name);
  }

  usedformat = mcformats[i_format];
  strcpy(tmp, usedformat.Name);
  usedformat.Name = tmp;
  if (strstr(low_format,"raw")) strcat(usedformat.Name," raw");
  if (strstr(low_format,"binary"))
  {
    if (strstr(low_format,"double")) strcat(usedformat.Name," binary double data");
    else strcat(usedformat.Name," binary float data");
    mcascii_only = 1;
  }

  /* Replaces vfprintf parameter name aliases */
  /* Header Footer */
  usedformat.Header       = mcuse_format_header(usedformat.Header);
  usedformat.Footer       = mcuse_format_header(usedformat.Footer);
  /* AssignTag */
  usedformat.AssignTag    = mcuse_format_tag(usedformat.AssignTag);
  /* BeginSection EndSection */
  usedformat.BeginSection = mcuse_format_section(usedformat.BeginSection);
  usedformat.EndSection   = mcuse_format_section(usedformat.EndSection);
  /*  BeginData  EndData  BeginErrors  EndErrors  BeginNcount  EndNcount  */
  usedformat.BeginData    = mcuse_format_data(usedformat.BeginData  );
  usedformat.EndData      = mcuse_format_data(usedformat.EndData    );
  usedformat.BeginErrors  = mcuse_format_data(usedformat.BeginErrors);
  usedformat.EndErrors    = mcuse_format_data(usedformat.EndErrors  );
  usedformat.BeginNcount  = mcuse_format_data(usedformat.BeginNcount);
  usedformat.EndNcount    = mcuse_format_data(usedformat.EndNcount  );

  return(usedformat);
} /* mcuse_format */

/*******************************************************************************
* mcclear_format: free format structure
*******************************************************************************/
void mcclear_format(struct mcformats_struct usedformat)
{
/* free format specification strings */
  if (usedformat.Name        ) free(usedformat.Name        );
  else return;
  if (usedformat.Header      ) free(usedformat.Header      );
  if (usedformat.Footer      ) free(usedformat.Footer      );
  if (usedformat.AssignTag   ) free(usedformat.AssignTag   );
  if (usedformat.BeginSection) free(usedformat.BeginSection);
  if (usedformat.EndSection  ) free(usedformat.EndSection  );
  if (usedformat.BeginData   ) free(usedformat.BeginData   );
  if (usedformat.EndData     ) free(usedformat.EndData     );
  if (usedformat.BeginErrors ) free(usedformat.BeginErrors );
  if (usedformat.EndErrors   ) free(usedformat.EndErrors   );
  if (usedformat.BeginNcount ) free(usedformat.BeginNcount );
  if (usedformat.EndNcount   ) free(usedformat.EndNcount   );
} /* mcclear_format */

/*******************************************************************************
* mcuse_file: will save data/sim files
*******************************************************************************/
static void mcuse_file(char *file)
{
  if (file && strcmp(file, "NULL"))
    mcsiminfo_name = file;
  else {
    char *filename=(char*)malloc(CHAR_BUF_LENGTH);
    sprintf(filename, "%s_%li", mcinstrument_name, mcstartdate);
    mcsiminfo_name = filename;
  }
  mcsingle_file  = 1;
}

/*******************************************************************************
* mcset_ncount: set total number of rays to generate
*******************************************************************************/
void mcset_ncount(unsigned long long int count)
{
  mcncount = count;
}

/* mcget_ncount: get total number of rays to generate */
unsigned long long int mcget_ncount(void)
{
  return mcncount;
}

/* mcget_run_num: get curent number of rays in TRACE */
unsigned long long int mcget_run_num(void)
{
  return mcrun_num;
}

/* mcsetn_arg: get ncount from a string argument */
static void
mcsetn_arg(char *arg)
{
  mcset_ncount((long long int) strtod(arg, NULL));
}

/* mcsetseed: set the random generator seed from a string argument */
static void
mcsetseed(char *arg)
{
  mcseed = atol(arg);
  if(mcseed) {
#ifdef USE_MPI
    if (mpi_node_count > 1) srandom(mcseed+mpi_node_rank);
    else
#endif
    srandom(mcseed);
  } else {
    fprintf(stderr, "Error: seed must not be zero (mcsetseed)\n");
    exit(1);
  }
}

/* Following part is only embedded when not redundent with mcstas.h ========= */

#ifndef MCCODE_H

/* SECTION: MCDISPLAY support. =============================================== */

/*******************************************************************************
* Just output MCDISPLAY keywords to be caught by an external plotter client.
*******************************************************************************/

void mcdis_magnify(char *what){
  printf("MCDISPLAY: magnify('%s')\n", what);
}

void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2){
  printf("MCDISPLAY: multiline(2,%g,%g,%g,%g,%g,%g)\n",
         x1,y1,z1,x2,y2,z2);
}

void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n){
  int i;
  const double dx = (x2-x1)/(2*n+1);
  const double dy = (y2-y1)/(2*n+1);
  const double dz = (z2-z1)/(2*n+1);

  for(i = 0; i < n+1; i++)
    mcdis_line(x1 + 2*i*dx,     y1 + 2*i*dy,     z1 + 2*i*dz,
	       x1 + (2*i+1)*dx, y1 + (2*i+1)*dy, z1 + (2*i+1)*dz);
}

void mcdis_multiline(int count, ...){
  va_list ap;
  double x,y,z;

  printf("MCDISPLAY: multiline(%d", count);
  va_start(ap, count);
  while(count--)
    {
    x = va_arg(ap, double);
    y = va_arg(ap, double);
    z = va_arg(ap, double);
    printf(",%g,%g,%g", x, y, z);
    }
  va_end(ap);
  printf(")\n");
}

void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height){
  /* draws a rectangle in the plane           */
  /* x is ALWAYS width and y is ALWAYS height */
  if (strcmp("xy", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y - height/2, z,
		    x + width/2, y - height/2, z,
		    x + width/2, y + height/2, z,
		    x - width/2, y + height/2, z,
		    x - width/2, y - height/2, z);
  } else if (strcmp("xz", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y, z - height/2,
		    x + width/2, y, z - height/2,
		    x + width/2, y, z + height/2,
		    x - width/2, y, z + height/2,
		    x - width/2, y, z - height/2);
  } else if (strcmp("yz", plane)==0) {
    mcdis_multiline(5,
		    x, y - height/2, z - width/2,
		    x, y - height/2, z + width/2,
		    x, y + height/2, z + width/2,
		    x, y + height/2, z - width/2,
		    x, y - height/2, z - width/2);
  } else {

    fprintf(stderr, "Error: Definition of plane %s unknown\n", plane);
    exit(1);
  }
}

/*  draws a box with center at (x, y, z) and
    width (deltax), height (deltay), length (deltaz) */
void mcdis_box(double x, double y, double z,
	       double width, double height, double length){

  mcdis_rectangle("xy", x, y, z-length/2, width, height);
  mcdis_rectangle("xy", x, y, z+length/2, width, height);
  mcdis_line(x-width/2, y-height/2, z-length/2,
	     x-width/2, y-height/2, z+length/2);
  mcdis_line(x-width/2, y+height/2, z-length/2,
	     x-width/2, y+height/2, z+length/2);
  mcdis_line(x+width/2, y-height/2, z-length/2,
	     x+width/2, y-height/2, z+length/2);
  mcdis_line(x+width/2, y+height/2, z-length/2,
	     x+width/2, y+height/2, z+length/2);
}

void mcdis_circle(char *plane, double x, double y, double z, double r){
  printf("MCDISPLAY: circle('%s',%g,%g,%g,%g)\n", plane, x, y, z, r);
}

/* SECTION: coordinates handling ============================================ */

/*******************************************************************************
* Since we use a lot of geometric calculations using Cartesian coordinates,
* we collect some useful routines here. However, it is also permissible to
* work directly on the underlying struct coords whenever that is most
* convenient (that is, the type Coords is not abstract).
*
* Coordinates are also used to store rotation angles around x/y/z axis.
*
* Since coordinates are used much like a basic type (such as double), the
* structure itself is passed and returned, rather than a pointer.
*
* At compile-time, the values of the coordinates may be unknown (for example
* a motor position). Hence coordinates are general expressions and not simple
* numbers. For this we used the type Coords_exp which has three CExp
* fields. For runtime (or calculations possible at compile time), we use
* Coords which contains three double fields.
*******************************************************************************/

/* coords_set: Assign coordinates. */
Coords
coords_set(MCNUM x, MCNUM y, MCNUM z)
{
  Coords a;

  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

/* coords_get: get coordinates. Required when 'x','y','z' are #defined as ray pars */
Coords
coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z)
{
  *x = a.x;
  *y = a.y;
  *z = a.z;
  return a;
}

/* coords_add: Add two coordinates. */
Coords
coords_add(Coords a, Coords b)
{
  Coords c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_sub: Subtract two coordinates. */
Coords
coords_sub(Coords a, Coords b)
{
  Coords c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_neg: Negate coordinates. */
Coords
coords_neg(Coords a)
{
  Coords b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;
  return b;
}

/* coords_scale: Scale a vector. */
Coords coords_scale(Coords b, double scale) {
  Coords a;

  a.x = b.x*scale;
  a.y = b.y*scale;
  a.z = b.z*scale;
  return a;
}

/* coords_sp: Scalar product: a . b */
double coords_sp(Coords a, Coords b) {
  double value;

  value = a.x*b.x + a.y*b.y + a.z*b.z;
  return value;
}

/* coords_xp: Cross product: a = b x c. */
Coords coords_xp(Coords b, Coords c) {
  Coords a;

  a.x = b.y*c.z - c.y*b.z;
  a.y = b.z*c.x - c.z*b.x;
  a.z = b.x*c.y - c.x*b.y;
  return a;
}

/* coords_mirror: Mirror a in plane (through the origin) defined by normal n*/
Coords coords_mirror(Coords a, Coords n) {
  double t = scalar_prod(n.x, n.y, n.z, n.x, n.y, n.z);
  Coords b;
  if (t!=1) {
    t = sqrt(t);
    n.x /= t;
    n.y /= t;
    n.z /= t;
  }
  t=scalar_prod(a.x, a.y, a.z, n.x, n.y, n.z);
  b.x = a.x-2*t*n.x;
  b.y = a.y-2*t*n.y;
  b.z = a.z-2*t*n.z;
  return b;
}

/* coords_print: Print out vector values. */
void coords_print(Coords a) {

  fprintf(stdout, "(%f, %f, %f)\n", a.x, a.y, a.z);
  return;
}

mcstatic inline void coords_norm(Coords* c) {
	double temp = coords_sp(*c,*c);

	// Skip if we will end dividing by zero
	if (temp == 0) return;

	temp = sqrt(temp);

	c->x /= temp;
	c->y /= temp;
	c->z /= temp;
}

/*******************************************************************************
* The Rotation type implements a rotation transformation of a coordinate
* system in the form of a double[3][3] matrix.
*
* Contrary to the Coords type in coords.c, rotations are passed by
* reference. Functions that yield new rotations do so by writing to an
* explicit result parameter; rotations are not returned from functions. The
* reason for this is that arrays cannot by returned from functions (though
* structures can; thus an alternative would have been to wrap the
* double[3][3] array up in a struct). Such are the ways of C programming.
*
* A rotation represents the tranformation of the coordinates of a vector when
* changing between coordinate systems that are rotated with respect to each
* other. For example, suppose that coordinate system Q is rotated 45 degrees
* around the Z axis with respect to coordinate system P. Let T be the
* rotation transformation representing a 45 degree rotation around Z. Then to
* get the coordinates of a vector r in system Q, apply T to the coordinates
* of r in P. If r=(1,0,0) in P, it will be (sqrt(1/2),-sqrt(1/2),0) in
* Q. Thus we should be careful when interpreting the sign of rotation angles:
* they represent the rotation of the coordinate systems, not of the
* coordinates (which has opposite sign).
*******************************************************************************/

/*******************************************************************************
* rot_set_rotation: Get transformation for rotation first phx around x axis,
* then phy around y, then phz around z.
*******************************************************************************/
void
rot_set_rotation(Rotation t, double phx, double phy, double phz)
{
  if ((phx == 0) && (phy == 0) && (phz == 0)) {
    t[0][0] = 1.0;
    t[0][1] = 0.0;
    t[0][2] = 0.0;
    t[1][0] = 0.0;
    t[1][1] = 1.0;
    t[1][2] = 0.0;
    t[2][0] = 0.0;
    t[2][1] = 0.0;
    t[2][2] = 1.0;
  } else {
    double cx = cos(phx);
    double sx = sin(phx);
    double cy = cos(phy);
    double sy = sin(phy);
    double cz = cos(phz);
    double sz = sin(phz);

    t[0][0] = cy*cz;
    t[0][1] = sx*sy*cz + cx*sz;
    t[0][2] = sx*sz - cx*sy*cz;
    t[1][0] = -cy*sz;
    t[1][1] = cx*cz - sx*sy*sz;
    t[1][2] = sx*cz + cx*sy*sz;
    t[2][0] = sy;
    t[2][1] = -sx*cy;
    t[2][2] = cx*cy;
  }
}

/*******************************************************************************
* rot_test_identity: Test if rotation is identity
*******************************************************************************/
int
rot_test_identity(Rotation t)
{
  return (t[0][0] + t[1][1] + t[2][2] == 3);
}

/*******************************************************************************
* rot_mul: Matrix multiplication of transformations (this corresponds to
* combining transformations). After rot_mul(T1, T2, T3), doing T3 is
* equal to doing first T2, then T1.
* Note that T3 must not alias (use the same array as) T1 or T2.
*******************************************************************************/
void
rot_mul(Rotation t1, Rotation t2, Rotation t3)
{
  if (rot_test_identity(t1)) {
    rot_copy(t3, t2);
  } else if (rot_test_identity(t2)) {
    rot_copy(t3, t1);
  } else {
    int i,j;
    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
	t3[i][j] = t1[i][0]*t2[0][j] + t1[i][1]*t2[1][j] + t1[i][2]*t2[2][j];
  }
}

/*******************************************************************************
* rot_copy: Copy a rotation transformation (arrays cannot be assigned in C).
*******************************************************************************/
void
rot_copy(Rotation dest, Rotation src)
{
  int i,j;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      dest[i][j] = src[i][j];
}

/*******************************************************************************
* rot_transpose: Matrix transposition, which is inversion for Rotation matrices
*******************************************************************************/
void
rot_transpose(Rotation src, Rotation dst)
{
  dst[0][0] = src[0][0];
  dst[0][1] = src[1][0];
  dst[0][2] = src[2][0];
  dst[1][0] = src[0][1];
  dst[1][1] = src[1][1];
  dst[1][2] = src[2][1];
  dst[2][0] = src[0][2];
  dst[2][1] = src[1][2];
  dst[2][2] = src[2][2];
}

/*******************************************************************************
* rot_apply: returns t*a
*******************************************************************************/
Coords
rot_apply(Rotation t, Coords a)
{
  Coords b;
  if (rot_test_identity(t)) {
    return a;
  } else {
    b.x = t[0][0]*a.x + t[0][1]*a.y + t[0][2]*a.z;
    b.y = t[1][0]*a.x + t[1][1]*a.y + t[1][2]*a.z;
    b.z = t[2][0]*a.x + t[2][1]*a.y + t[2][2]*a.z;
    return b;
  }
}

/**
 * Pretty-printing of rotation matrices.
 */
void rot_print(Rotation rot) {
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[0][0], rot[0][1], rot[0][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[1][0], rot[1][1], rot[1][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n\n",
			rot[2][0], rot[2][1], rot[2][2]);
}

/**
 * Vector product: used by vec_prod (mccode-r.h). Use coords_xp for Coords.
 */
mcstatic inline void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
    *x = (y1)*(z2) - (y2)*(z1);
    *y = (z1)*(x2) - (z2)*(x1);
    *z = (x1)*(y2) - (x2)*(y1);
}

/**
 * Scalar product: use coords_sp for Coords.
 */
mcstatic inline double scalar_prod(
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return ((x1 * x2) + (y1 * y2) + (z1 * z2));
}

/*******************************************************************************
* mccoordschange: applies rotation to (x y z) and (vx vy vz) and Spin (sx,sy,sz)
*******************************************************************************/
void
mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *sx, double *sy, double *sz)
{
  Coords b, c;

  if (x != NULL && y != NULL && z != NULL) {
    b.x = *x;
    b.y = *y;
    b.z = *z;
    c = rot_apply(t, b);
    b = coords_add(c, a);
    *x = b.x;
    *y = b.y;
    *z = b.z;
  }

  if (vx != NULL && vy != NULL && vz != NULL) mccoordschange_polarisation(t, vx, vy, vz);

  if (sx != NULL && sy != NULL && sz != NULL) mccoordschange_polarisation(t, sx, sy, sz);

}

/*******************************************************************************
* mccoordschange_polarisation: applies rotation to vector (sx sy sz)
*******************************************************************************/
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *sx;
  b.y = *sy;
  b.z = *sz;
  c = rot_apply(t, b);
  *sx = c.x;
  *sy = c.y;
  *sz = c.z;
}

/* SECTION: vector math  ==================================================== */

/* normal_vec: Compute normal vector to (x,y,z). */
void normal_vec(double *nx, double *ny, double *nz,
                double x, double y, double z)
{
  double ax = fabs(x);
  double ay = fabs(y);
  double az = fabs(z);
  double l;
  if(x == 0 && y == 0 && z == 0)
  {
    *nx = 0;
    *ny = 0;
    *nz = 0;
    return;
  }
  if(ax < ay)
  {
    if(ax < az)
    {                           /* Use X axis */
      l = sqrt(z*z + y*y);
      *nx = 0;
      *ny = z/l;
      *nz = -y/l;
      return;
    }
  }
  else
  {
    if(ay < az)
    {                           /* Use Y axis */
      l = sqrt(z*z + x*x);
      *nx = z/l;
      *ny = 0;
      *nz = -x/l;
      return;
    }
  }
  /* Use Z axis */
  l = sqrt(y*y + x*x);
  *nx = y/l;
  *ny = -x/l;
  *nz = 0;
} /* normal_vec */

/*******************************************************************************
 * solve_2nd_order: second order equation solve: A*t^2 + B*t + C = 0
 * solve_2nd_order(&t1, NULL, A,B,C)
 *   returns 0 if no solution was found, or set 't1' to the smallest positive
 *   solution.
 * solve_2nd_order(&t1, &t2, A,B,C)
 *   same as with &t2=NULL, but also returns the second solution.
 * EXAMPLE usage for intersection of a trajectory with a plane in gravitation
 * field (gx,gy,gz):
 * The neutron starts at point r=(x,y,z) with velocityv=(vx vy vz). The plane
 * has a normal vector n=(nx,ny,nz) and contains the point W=(wx,wy,wz).
 * The problem consists in solving the 2nd order equation:
 *      1/2.n.g.t^2 + n.v.t + n.(r-W) = 0
 * so that A = 0.5 n.g; B = n.v; C = n.(r-W);
 * Without acceleration, t=-n.(r-W)/n.v
 ******************************************************************************/
int solve_2nd_order(double *t1, double *t2,
                  double A,  double B,  double C)
{
  int ret=0;

  if (!t1) return 0;
  *t1 = 0;
  if (t2) *t2=0;

  if (fabs(A) < 1E-10) /* approximate to linear equation: A ~ 0 */
  {
    if (B) {  *t1 = -C/B; ret=1; if (t2) *t2=*t1; }
    /* else no intersection: A=B=0 ret=0 */
  }
  else
  {
    double D;
    D = B*B - 4*A*C;
    if (D >= 0) /* Delta > 0: two solutions */
    {
      double sD, dt1, dt2;
      sD = sqrt(D);
      dt1 = (-B + sD)/2/A;
      dt2 = (-B - sD)/2/A;
      /* we identify very small values with zero */
      if (fabs(dt1) < 1e-10) dt1=0.0;
      if (fabs(dt2) < 1e-10) dt2=0.0;

      /* now we choose the smallest positive solution */
      if      (dt1<=0.0 && dt2>0.0) ret=2; /* dt2 positive */
      else if (dt2<=0.0 && dt1>0.0) ret=1; /* dt1 positive */
      else if (dt1> 0.0 && dt2>0.0)
      {  if (dt1 < dt2) ret=1; else ret=2; } /* all positive: min(dt1,dt2) */
      /* else two solutions are negative. ret=-1 */
      if (ret==1) { *t1 = dt1;  if (t2) *t2=dt2; }
      else        { *t1 = dt2;  if (t2) *t2=dt1; }
      ret=2;  /* found 2 solutions and t1 is the positive one */
    } /* else Delta <0: no intersection. ret=0 */
  }
  return(ret);
} /* solve_2nd_order */

/*******************************************************************************
 * randvec_target_circle: Choose random direction towards target at (x,y,z)
 * with given radius.
 * If radius is zero, choose random direction in full 4PI, no target.
 ******************************************************************************/
void
randvec_target_circle(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double radius)
{
  double l2, phi, theta, nx, ny, nz, xt, yt, zt, xu, yu, zu;

  if(radius == 0.0)
  {
    /* No target, choose uniformly a direction in full 4PI solid angle. */
    theta = acos (1 - rand0max(2));
    phi = rand0max(2 * PI);
    if(solid_angle)
      *solid_angle = 4*PI;
    nx = 1;
    ny = 0;
    nz = 0;
    yi = sqrt(xi*xi+yi*yi+zi*zi);
    zi = 0;
    xi = 0;
  }
  else
  {
    double costheta0;
    l2 = xi*xi + yi*yi + zi*zi; /* sqr Distance to target. */
    costheta0 = sqrt(l2/(radius*radius+l2));
    if (radius < 0) costheta0 *= -1;
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
        *solid_angle = 2*PI*(1 - costheta0);
    }

    /* Now choose point uniformly on circle surface within angle theta0 */
    theta = acos (1 - rand0max(1 - costheta0)); /* radius on circle */
    phi = rand0max(2 * PI); /* rotation on circle at given radius */
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around a
       perpendicular axis u=i x n and then angle phi around i. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, xu, yu, zu);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xi, yi, zi);
} /* randvec_target_circle */

/*******************************************************************************
 * randvec_target_rect_angular: Choose random direction towards target at
 * (xi,yi,zi) with given ANGULAR dimension height x width. height=phi_x,
 * width=phi_y (radians)
 * If height or width is zero, choose random direction in full 4PI, no target.
 *******************************************************************************/
void
randvec_target_rect_angular(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double width, double height, Rotation A)
{
  double theta, phi, nx, ny, nz, xt, yt, zt, xu, yu, zu;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
      *solid_angle = 2*fabs(width*sin(height/2));
    }

    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Now choose point uniformly on quadrant within angle theta0/phi0 */
    theta = width*randpm1()/2.0;
    phi   = height*randpm1()/2.0;
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle phi around
       n, and then theta around u. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, phi, nx, ny, nz);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, theta, xu,  yu,  zu);

  /* Go back to local coordinate system */
  tmp = coords_set(*xo, *yo, *zo);
  tmp = rot_apply(A, tmp);
  coords_get(tmp, &*xo, &*yo, &*zo);

} /* randvec_target_rect_angular */

/*******************************************************************************
 * randvec_target_rect_real: Choose random direction towards target at (xi,yi,zi)
 * with given dimension height x width (in meters !).
 *
 * Local emission coordinate is taken into account and corrected for 'order' times.
 * (See remarks posted to mcstas-users by George Apostolopoulus <gapost@ipta.demokritos.gr>)
 *
 * If height or width is zero, choose random direction in full 4PI, no target.
 *
 * Traditionally, this routine had the name randvec_target_rect - this is now a
 * a define (see mcstas-r.h) pointing here. If you use the old rouine, you are NOT
 * taking the local emmission coordinate into account.
*******************************************************************************/

void
randvec_target_rect_real(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi,
               double width, double height, Rotation A,
               double lx, double ly, double lz, int order)
{
  double dx, dy, dist, dist_p, nx, ny, nz, mx, my, mz, n_norm, m_norm;
  double cos_theta;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {

    /* Now choose point uniformly on rectangle within width x height */
    dx = width*randpm1()/2.0;
    dy = height*randpm1()/2.0;

    /* Determine distance to target plane*/
    dist = sqrt(xi*xi + yi*yi + zi*zi);
    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Determine vector normal to trajectory axis (z) and gravity [0 1 0] */
    vec_prod(nx, ny, nz, xi, yi, zi, 0, 1, 0);

    /* This now defines the x-axis, normalize: */
    n_norm=sqrt(nx*nx + ny*ny + nz*nz);
    nx = nx/n_norm;
    ny = ny/n_norm;
    nz = nz/n_norm;

    /* Now, determine our y-axis (vertical in many cases...) */
    vec_prod(mx, my, mz, xi, yi, zi, nx, ny, nz);
    m_norm=sqrt(mx*mx + my*my + mz*mz);
    mx = mx/m_norm;
    my = my/m_norm;
    mz = mz/m_norm;

    /* Our output, random vector can now be defined by linear combination: */

    *xo = xi + dx * nx + dy * mx;
    *yo = yi + dx * ny + dy * my;
    *zo = zi + dx * nz + dy * mz;

    /* Go back to local coordinate system */
    tmp = coords_set(*xo, *yo, *zo);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &*xo, &*yo, &*zo);

    /* Go back to local coordinate system */
    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    if (solid_angle) {
      /* Calculate vector from local point to remote random point */
      lx = *xo - lx;
      ly = *yo - ly;
      lz = *zo - lz;
      dist_p = sqrt(lx*lx + ly*ly + lz*lz);

      /* Adjust the 'solid angle' */
      /* 1/r^2 to the chosen point times cos(\theta) between the normal */
      /* vector of the target rectangle and direction vector of the chosen point. */
      cos_theta = (xi * lx + yi * ly + zi * lz) / (dist * dist_p);
      *solid_angle = width * height / (dist_p * dist_p);
      int counter;
      for (counter = 0; counter < order; counter++) {
	*solid_angle = *solid_angle * cos_theta;
      }
    }
  }
} /* randvec_target_rect_real */

/* SECTION: random numbers ================================================== */

/*
 * Copyright (c) 1983 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the University of California, Berkeley.  The name of the
 * University may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 * This is derived from the Berkeley source:
 *        @(#)random.c        5.5 (Berkeley) 7/6/88
 * It was reworked for the GNU C Library by Roland McGrath.
 * Rewritten to use reentrant functions by Ulrich Drepper, 1995.
 */

/*******************************************************************************
* Modified for McStas from glibc 2.0.7pre1 stdlib/random.c and
* stdlib/random_r.c.
*
* This way random() is more than four times faster compared to calling
* standard glibc random() on ix86 Linux, probably due to multithread support,
* ELF shared library overhead, etc. It also makes McStas generated
* simulations more portable (more likely to behave identically across
* platforms, important for parrallel computations).
*******************************************************************************/


#define        TYPE_3                3
#define        BREAK_3                128
#define        DEG_3                31
#define        SEP_3                3

static mc_int32_t randtbl[DEG_3 + 1] =
  {
    TYPE_3,

    -1726662223, 379960547, 1735697613, 1040273694, 1313901226,
    1627687941, -179304937, -2073333483, 1780058412, -1989503057,
    -615974602, 344556628, 939512070, -1249116260, 1507946756,
    -812545463, 154635395, 1388815473, -1926676823, 525320961,
    -1009028674, 968117788, -123449607, 1284210865, 435012392,
    -2017506339, -911064859, -370259173, 1132637927, 1398500161,
    -205601318,
  };

static mc_int32_t *fptr = &randtbl[SEP_3 + 1];
static mc_int32_t *rptr = &randtbl[1];
static mc_int32_t *state = &randtbl[1];
#define rand_deg DEG_3
#define rand_sep SEP_3
static mc_int32_t *end_ptr = &randtbl[sizeof (randtbl) / sizeof (randtbl[0])];

mc_int32_t
mc_random (void)
{
  mc_int32_t result;

  *fptr += *rptr;
  /* Chucking least random bit.  */
  result = (*fptr >> 1) & 0x7fffffff;
  ++fptr;
  if (fptr >= end_ptr)
  {
    fptr = state;
    ++rptr;
  }
  else
  {
    ++rptr;
    if (rptr >= end_ptr)
      rptr = state;
  }
  return result;
}

void
mc_srandom (unsigned int x)
{
  /* We must make sure the seed is not 0.  Take arbitrarily 1 in this case.  */
  state[0] = x ? x : 1;
  {
    long int i;
    for (i = 1; i < rand_deg; ++i)
    {
      /* This does:
         state[i] = (16807 * state[i - 1]) % 2147483647;
         but avoids overflowing 31 bits.  */
      long int hi = state[i - 1] / 127773;
      long int lo = state[i - 1] % 127773;
      long int test = 16807 * lo - 2836 * hi;
      state[i] = test + (test < 0 ? 2147483647 : 0);
    }
    fptr = &state[rand_sep];
    rptr = &state[0];
    for (i = 0; i < 10 * rand_deg; ++i)
      random ();
  }
}

/* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
/* See http://www.math.keio.ac.jp/~matumoto/emt.html for original source. */


/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using mt_srandom(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

#include <stdio.h>

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void mt_srandom(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(init_key, key_length)
unsigned long init_key[], key_length;
{
    int i, j, k;
    mt_srandom(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt_random(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if mt_srandom() has not been called, */
            mt_srandom(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK

/* End of "Mersenne Twister". */

/* End of McStas random number routine. */

/* randnorm: generate a random number from normal law */
double
randnorm(void)
{
  static double v1, v2, s;
  static int phase = 0;
  double X, u1, u2;

  if(phase == 0)
  {
    do
    {
      u1 = rand01();
      u2 = rand01();
      v1 = 2*u1 - 1;
      v2 = 2*u2 - 1;
      s = v1*v1 + v2*v2;
    } while(s >= 1 || s == 0);

    X = v1*sqrt(-2*log(s)/s);
  }
  else
  {
    X = v2*sqrt(-2*log(s)/s);
  }

  phase = 1 - phase;
  return X;
}

/**
 * Generate a random number from -1 to 1 with triangle distribution
 */
double randtriangle(void) {
	double randnum = rand01();
	if (randnum>0.5) return(1-sqrt(2*(randnum-0.5)));
	else return(sqrt(2*randnum)-1);
}

/**
 * Random number between 0.0 and 1.0 (including?)
 */
double rand01() {
	double rand;
	rand = (double) random();
	rand /= (double) MC_RAND_MAX + 1;
	return rand;
}

/**
 * Return a random number between 1 and -1
 */
double randpm1() {
	double rand;
	rand = (double) random();
	rand /= ((double) MC_RAND_MAX + 1) / 2;
	rand -= 1;
	return rand;
}

/**
 * Return a random number between 0 and max.
 */
double rand0max(double max) {
	double rand;
	rand = (double) random();
	rand /= ((double) MC_RAND_MAX + 1) / max;
	return rand;
}

/**
 * Return a random number between min and max.
 */
double randminmax(double min, double max) {
	return rand0max(max - min) + max;
}

/* SECTION: main and signal handlers ======================================== */

/*******************************************************************************
* mchelp: displays instrument executable help with possible options
*******************************************************************************/
static void
mchelp(char *pgmname)
{
  int i;

  fprintf(stderr, "%s (%s) instrument simulation, generated with " MCCODE_STRING " (" MCCODE_DATE ")\n", mcinstrument_name, mcinstrument_source);
  fprintf(stderr, "Usage: %s [options] [parm=value ...]\n", pgmname);
  fprintf(stderr,
"Options are:\n"
"  -s SEED   --seed=SEED      Set random seed (must be != 0)\n"
"  -n COUNT  --ncount=COUNT   Set number of @MCCODE_PARTICULE@s to simulate.\n"
"  -d DIR    --dir=DIR        Put all data files in directory DIR.\n"
"  -f FILE   --file=FILE      Put all data in a single file.\n"
"  -t        --trace          Enable trace of @MCCODE_PARTICULE@s through instrument.\n"
"  -g        --gravitation    Enable gravitation for all trajectories.\n"
"  -a        --data-only      Do not put any headers in the data files.\n"
"  --no-output-files          Do not write any data files.\n"
"  -h        --help           Show this help message.\n"
"  -i        --info           Detailed instrument information.\n"
"  --format=FORMAT            Output data files using format FORMAT\n"
"                             (use option +a to include text header in files\n"
#ifdef USE_MPI
"This instrument has been compiled with MPI support. Use 'mpirun'.\n"
#endif
"\n"
);
  if(mcnumipar > 0)
  {
    fprintf(stderr, "Instrument parameters are:\n");
    for(i = 0; i < mcnumipar; i++)
      if (mcinputtable[i].val && strlen(mcinputtable[i].val))
        fprintf(stderr, "  %-16s(%s) [default='%s']\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name),
        mcinputtable[i].val);
      else
        fprintf(stderr, "  %-16s(%s)\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name));
  }
  fprintf(stderr, "Available output formats are (default is %s):\n  ", mcformat.Name);
  for (i=0; i < mcNUMFORMATS; fprintf(stderr,"\"%s\" " , mcformats[i++].Name) );
  fprintf(stderr, "\n  Format modifiers: FORMAT may be followed by 'binary float' or \n");
  fprintf(stderr, "  'binary double' to save data blocks as binary. This removes text headers.\n");
  fprintf(stderr, "  The " FLAVOR_UPPER "_FORMAT environment variable may set the default FORMAT to use.\n");
#ifndef NOSIGNALS
  fprintf(stderr, "Known signals are: "
#ifdef SIGUSR1
  "USR1 (status) "
#endif
#ifdef SIGUSR2
  "USR2 (save) "
#endif
#ifdef SIGBREAK
  "BREAK (save) "
#endif
#ifdef SIGTERM
  "TERM (save and exit)"
#endif
  "\n");
#endif /* !NOSIGNALS */
} /* mchelp */


/* mcshowhelp: show help and exit with 0 */
static void
mcshowhelp(char *pgmname)
{
  mchelp(pgmname);
  exit(0);
}

/* mcusage: display usage when error in input arguments and exit with 1 */
static void
mcusage(char *pgmname)
{
  fprintf(stderr, "Error: incorrect command line arguments\n");
  mchelp(pgmname);
  exit(1);
}

/* mcenabletrace: enable trace/mcdisplay or error if requires recompile */
static void
mcenabletrace(void)
{
 if(mctraceenabled)
  mcdotrace = 1;
 else
 {
   fprintf(stderr,
           "Error: trace not enabled (mcenabletrace)\n"
           "Please re-run the McStas compiler "
                   "with the --trace option, or rerun the\n"
           "C compiler with the MC_TRACE_ENABLED macro defined.\n");
   exit(1);
 }
}

/*******************************************************************************
 * mcuse_dir: set data/sim storage directory and create it,
 * or exit with error if exists
 ******************************************************************************/
static void
mcuse_dir(char *dir)
{
  if (!dir || !strlen(dir)) return;
#ifdef MC_PORTABLE
  fprintf(stderr, "Error: "
          "Directory output cannot be used with portable simulation (mcuse_dir)\n");
  exit(1);
#else  /* !MC_PORTABLE */
#ifdef USE_MPI
  if(mpi_node_rank == mpi_node_root)
  {
#endif
    if(mkdir(dir, 0777)) {
#ifndef DANSE
      fprintf(stderr, "Error: unable to create directory '%s' (mcuse_dir)\n", dir);
      fprintf(stderr, "(Maybe the directory already exists?)\n");
      exit(1);
#endif
    }
#ifdef USE_MPI
  }
#endif
  /* handle file://directory URL type */
  if (strncmp(dir, "file:/", strlen("file:/")))
    mcdirname = dir;
  else
    mcdirname = dir+strlen("file:/");

  /* remove trailing PATHSEP (if any) */
  while (strlen(mcdirname) && mcdirname[strlen(mcdirname) - 1] == MC_PATHSEP_C)
    mcdirname[strlen(mcdirname) - 1]='\0';
#endif /* !MC_PORTABLE */
} /* mcuse_dir */

/*******************************************************************************
* mcinfo: display instrument simulation info to stdout and exit
*******************************************************************************/
static void
mcinfo(void)
{
  if(strstr(mcformat.Name,"NeXus"))
    exit(fprintf(stderr,"Error: Can not display instrument information in NeXus binary format\n"));
  mcsiminfo_init(stdout);
  mcsiminfo_close();
  exit(0); /* includes MPI_Finalize in MPI mode */
} /* mcinfo */

/*******************************************************************************
* mcreadparams: request parameters from the prompt (or use default)
*******************************************************************************/
void
mcreadparams(void)
{
  int i,j,status;
  static char buf[CHAR_BUF_LENGTH];
  char *p;
  int len;

  MPI_MASTER(printf("Instrument parameters for %s (%s)\n",
                    mcinstrument_name, mcinstrument_source));

  for(i = 0; mcinputtable[i].name != 0; i++)
  {
    do
    {
      MPI_MASTER(
                 if (mcinputtable[i].val && strlen(mcinputtable[i].val))
                   printf("Set value of instrument parameter %s (%s) [default='%s']:\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name), mcinputtable[i].val);
                 else
                   printf("Set value of instrument parameter %s (%s):\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name));
                 fflush(stdout);
                 );
#ifdef USE_MPI
      if(mpi_node_rank == mpi_node_root)
        {
          p = fgets(buf, CHAR_BUF_LENGTH, stdin);
          if(p == NULL)
            {
              fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
              exit(1);
            }
        }
      else
        p = buf;
      MPI_Bcast(buf, CHAR_BUF_LENGTH, MPI_CHAR, mpi_node_root, MPI_COMM_WORLD);
#else /* !USE_MPI */
      p = fgets(buf, CHAR_BUF_LENGTH, stdin);
      if(p == NULL)
        {
          fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
          exit(1);
        }
#endif /* USE_MPI */
      len = strlen(buf);
      if (!len || (len == 1 && (buf[0] == '\n' || buf[0] == '\r')))
      {
        if (mcinputtable[i].val && strlen(mcinputtable[i].val)) {
          strncpy(buf, mcinputtable[i].val, CHAR_BUF_LENGTH);  /* use default value */
          len = strlen(buf);
        }
      }
      for(j = 0; j < 2; j++)
      {
        if(len > 0 && (buf[len - 1] == '\n' || buf[len - 1] == '\r'))
        {
          len--;
          buf[len] = '\0';
        }
      }

      status = (*mcinputtypes[mcinputtable[i].type].getparm)
                   (buf, mcinputtable[i].par);
      if(!status)
      {
        (*mcinputtypes[mcinputtable[i].type].error)(mcinputtable[i].name, buf);
        if (!mcinputtable[i].val || strlen(mcinputtable[i].val)) {
          fprintf(stderr, "       Change %s default value in instrument definition.\n", mcinputtable[i].name);
          exit(1);
        }
      }
    } while(!status);
  }
} /* mcreadparams */

/*******************************************************************************
* mcparseoptions: parse command line arguments (options, parameters)
*******************************************************************************/
void
mcparseoptions(int argc, char *argv[])
{
  int i, j;
  char *p;
  int paramset = 0, *paramsetarray;
  char *usedir=NULL;

  /* Add one to mcnumipar to avoid allocating zero size memory block. */
  paramsetarray = malloc((mcnumipar + 1)*sizeof(*paramsetarray));
  if(paramsetarray == NULL)
  {
    fprintf(stderr, "Error: insufficient memory (mcparseoptions)\n");
    exit(1);
  }
  for(j = 0; j < mcnumipar; j++)
    {
      paramsetarray[j] = 0;
      if (mcinputtable[j].val != NULL && strlen(mcinputtable[j].val))
      {
        int  status;
        char buf[CHAR_BUF_LENGTH];
        strncpy(buf, mcinputtable[j].val, CHAR_BUF_LENGTH);
        status = (*mcinputtypes[mcinputtable[j].type].getparm)
                   (buf, mcinputtable[j].par);
        if(!status) fprintf(stderr, "Invalid '%s' default value %s in instrument definition (mcparseoptions)\n", mcinputtable[j].name, buf);
        else paramsetarray[j] = 1;
      } else {
        (*mcinputtypes[mcinputtable[j].type].getparm)
          (NULL, mcinputtable[j].par);
        paramsetarray[j] = 0;
      }
    }
  for(i = 1; i < argc; i++)
  {
    if(!strcmp("-s", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("-s", argv[i], 2))
      mcsetseed(&argv[i][2]);
    else if(!strcmp("--seed", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("--seed=", argv[i], 7))
      mcsetseed(&argv[i][7]);
    else if(!strcmp("-n", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("-n", argv[i], 2))
      mcsetn_arg(&argv[i][2]);
    else if(!strcmp("--ncount", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("--ncount=", argv[i], 9))
      mcsetn_arg(&argv[i][9]);
    else if(!strcmp("-d", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];  /* will create directory after parsing all arguments (end of this function) */
    else if(!strncmp("-d", argv[i], 2))
      usedir=&argv[i][2];
    else if(!strcmp("--dir", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];
    else if(!strncmp("--dir=", argv[i], 6))
      usedir=&argv[i][6];
    else if(!strcmp("-f", argv[i]) && (i + 1) < argc)
      mcuse_file(argv[++i]);
    else if(!strncmp("-f", argv[i], 2))
      mcuse_file(&argv[i][2]);
    else if(!strcmp("--file", argv[i]) && (i + 1) < argc)
      mcuse_file(argv[++i]);
    else if(!strncmp("--file=", argv[i], 7))
      mcuse_file(&argv[i][7]);
    else if(!strcmp("-h", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("--help", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("-i", argv[i])) {
      mcformat=mcuse_format(MCSTAS_FORMAT);
      mcinfo();
    }
    else if(!strcmp("--info", argv[i]))
      mcinfo();
    else if(!strcmp("-t", argv[i]))
      mcenabletrace();
    else if(!strcmp("--trace", argv[i]))
      mcenabletrace();
    else if(!strcmp("-a", argv[i]))
      mcascii_only = 1;
    else if(!strcmp("+a", argv[i]))
      mcascii_only = 0;
    else if(!strcmp("--data-only", argv[i]))
      mcascii_only = 1;
    else if(!strcmp("--gravitation", argv[i]))
      mcgravitation = 1;
    else if(!strcmp("-g", argv[i]))
      mcgravitation = 1;
    else if(!strncmp("--format=", argv[i], 9)) {
      mcascii_only = 0;
      mcformat=mcuse_format(&argv[i][9]);
    }
    else if(!strcmp("--format", argv[i]) && (i + 1) < argc) {
      mcascii_only = 0;
      mcformat=mcuse_format(argv[++i]);
    }
    else if(!strncmp("--format_data=", argv[i], 14) || !strncmp("--format-data=", argv[i], 14)) {
      mcascii_only = 0;
      mcformat_data=mcuse_format(&argv[i][14]);
    }
    else if((!strcmp("--format_data", argv[i]) || !strcmp("--format-data", argv[i])) && (i + 1) < argc) {
      mcascii_only = 0;
      mcformat_data=mcuse_format(argv[++i]);
    }
    else if(!strcmp("--no-output-files", argv[i]))
      mcdisable_output_files = 1;
    else if(argv[i][0] != '-' && (p = strchr(argv[i], '=')) != NULL)
    {
      *p++ = '\0';

      for(j = 0; j < mcnumipar; j++)
        if(!strcmp(mcinputtable[j].name, argv[i]))
        {
          int status;
          status = (*mcinputtypes[mcinputtable[j].type].getparm)(p,
                        mcinputtable[j].par);
          if(!status || !strlen(p))
          {
            (*mcinputtypes[mcinputtable[j].type].error)
              (mcinputtable[j].name, p);
            exit(1);
          }
          paramsetarray[j] = 1;
          paramset = 1;
          break;
        }
      if(j == mcnumipar)
      {                                /* Unrecognized parameter name */
        fprintf(stderr, "Error: unrecognized parameter %s (mcparseoptions)\n", argv[i]);
        exit(1);
      }
    }
    else if(argv[i][0] == '-') {
      fprintf(stderr, "Error: unrecognized option argument %s (mcparseoptions). Ignored.\n", argv[i++]);
    }
    else {
      fprintf(stderr, "Error: unrecognized argument %s (mcparseoptions). Aborting.\n", argv[i]);
      mcusage(argv[0]);
    }
  }
  if (!mcascii_only) {
    if (strstr(mcformat.Name,"binary")) fprintf(stderr, "Warning: %s files will contain text headers.\n         Use -a option to clean up.\n", mcformat.Name);
    strcat(mcformat.Name, " with text headers");
  }
  if(!paramset)
    mcreadparams();                /* Prompt for parameters if not specified. */
  else
  {
    for(j = 0; j < mcnumipar; j++)
      if(!paramsetarray[j])
      {
        fprintf(stderr, "Error: Instrument parameter %s left unset (mcparseoptions)\n",
                mcinputtable[j].name);
        exit(1);
      }
  }
  free(paramsetarray);
#ifdef USE_MPI
  if (mcdotrace) mpi_node_count=1; /* disable threading when in trace mode */
#endif
  if (usedir && strlen(usedir)) mcuse_dir(usedir);
} /* mcparseoptions */

#ifndef NOSIGNALS
mcstatic char  mcsig_message[256];


/*******************************************************************************
* sighandler: signal handler that makes simulation stop, and save results
*******************************************************************************/
void sighandler(int sig)
{
  /* MOD: E. Farhi, Sep 20th 2001: give more info */
  time_t t1, t0;
#define SIG_SAVE 0
#define SIG_TERM 1
#define SIG_STAT 2
#define SIG_ABRT 3

  printf("\n# McStas: [pid %i] Signal %i detected", getpid(), sig);
#ifdef USE_MPI
  printf("[proc %i]", mpi_node_rank);
#endif
#if defined(SIGUSR1) && defined(SIGUSR2) && defined(SIGKILL)
  if (!strcmp(mcsig_message, "sighandler") && (sig != SIGUSR1) && (sig != SIGUSR2))
  {
    printf("\n# Fatal : unrecoverable loop ! Suicide (naughty boy).\n");
    kill(0, SIGKILL); /* kill myself if error occurs within sighandler: loops */
  }
#endif
  switch (sig) {
#ifdef SIGINT
    case SIGINT : printf(" SIGINT (interrupt from terminal, Ctrl-C)"); sig = SIG_TERM; break;
#endif
#ifdef SIGILL
    case SIGILL  : printf(" SIGILL (Illegal instruction)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGFPE
    case SIGFPE  : printf(" SIGFPE (Math Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGSEGV
    case SIGSEGV : printf(" SIGSEGV (Mem Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGTERM
    case SIGTERM : printf(" SIGTERM (Termination)"); sig = SIG_TERM; break;
#endif
#ifdef SIGABRT
    case SIGABRT : printf(" SIGABRT (Abort)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGQUIT
    case SIGQUIT : printf(" SIGQUIT (Quit from terminal)"); sig = SIG_TERM; break;
#endif
#ifdef SIGTRAP
    case SIGTRAP : printf(" SIGTRAP (Trace trap)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGPIPE
    case SIGPIPE : printf(" SIGPIPE (Broken pipe)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGUSR1
    case SIGUSR1 : printf(" SIGUSR1 (Display info)"); sig = SIG_STAT; break;
#endif
#ifdef SIGUSR2
    case SIGUSR2 : printf(" SIGUSR2 (Save simulation)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGHUP
    case SIGHUP  : printf(" SIGHUP (Hangup/update)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGBUS
    case SIGBUS  : printf(" SIGBUS (Bus error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGURG
    case SIGURG  : printf(" SIGURG (Urgent socket condition)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGBREAK
    case SIGBREAK: printf(" SIGBREAK (Break signal, Ctrl-Break)"); sig = SIG_SAVE; break;
#endif
    default : printf(" (look at signal list for signification)"); sig = SIG_ABRT; break;
  }
  printf("\n");
  printf("# Simulation: %s (%s) \n", mcinstrument_name, mcinstrument_source);
  printf("# Breakpoint: %s ", mcsig_message);
  if (strstr(mcsig_message, "Save") && (sig == SIG_SAVE))
    sig = SIG_STAT;
  SIG_MESSAGE("sighandler");
  if (mcget_ncount() == 0)
    printf("(0 %%)\n" );
  else
  {
    printf("%.2f %% (%10.1f/%10.1f)\n", 100.0*mcget_run_num()/mcget_ncount(), 1.0*mcget_run_num(), 1.0*mcget_ncount());
  }
  t0 = (time_t)mcstartdate;
  t1 = time(NULL);
  printf("# Date:      %s", ctime(&t1));
  printf("# Started:   %s", ctime(&t0));

  if (sig == SIG_STAT)
  {
    printf("# McStas: Resuming simulation (continue)\n");
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_SAVE)
  {
    printf("# McStas: Saving data and resume simulation (continue)\n");
    mcsave(NULL);
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_TERM)
  {
    printf("# McStas: Finishing simulation (save results and exit)\n");
    mcfinally();
    exit(0);
  }
  else
  {
    fflush(stdout);
    perror("# Last I/O Error");
    printf("# McStas: Simulation stop (abort)\n");
    exit(-1);
  }
#undef SIG_SAVE
#undef SIG_TERM
#undef SIG_STAT
#undef SIG_ABRT

} /* sighandler */
#endif /* !NOSIGNALS */

/*******************************************************************************
* mccode_main: McCode main() function.
*******************************************************************************/
int mccode_main(int argc, char *argv[])
{
/*  double run_num = 0; */
  time_t t;
#ifdef USE_MPI
  char mpi_node_name[MPI_MAX_PROCESSOR_NAME];
  int  mpi_node_name_len;
#endif /* USE_MPI */

#ifdef MAC
  argc = ccommand(&argv);
#endif

#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_node_count); /* get number of nodes */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_node_rank);
  MPI_Get_processor_name(mpi_node_name, &mpi_node_name_len);
#endif /* USE_MPI */

#ifdef USE_MPI
/* *** print number of nodes *********************************************** */
  if (mpi_node_count > 1) {
    MPI_MASTER(
    printf("Simulation '%s' (%s): running on %i nodes (master is '%s', MPI version %i.%i).\n",
      mcinstrument_name, mcinstrument_source, mpi_node_count, mpi_node_name, MPI_VERSION, MPI_SUBVERSION);
    );
    /* adapt random seed for each node */
    mcseed=(long)(time(&t) + mpi_node_rank);
    srandom(mcseed);
    t += mpi_node_rank;
  } else
#else /* !USE_MPI */
  mcseed=(long)time(&t);
  srandom(mcseed);
#endif /* USE_MPI */
  mcstartdate = (long)t;  /* set start date before parsing options and creating sim file */

/* *** parse options ******************************************************* */
  SIG_MESSAGE("main (Start)");
  mcformat=mcuse_format(getenv(FLAVOR_UPPER "_FORMAT") ?
                        getenv(FLAVOR_UPPER "_FORMAT") : MCSTAS_FORMAT);
  /* default is to output as McStas format */
  mcformat_data.Name=NULL;
  /* read simulation parameters and options */
  mcparseoptions(argc, argv); /* sets output dir and format */
  if (strstr(mcformat.Name, "NeXus")) {
    if (mcformat_data.Name) mcclear_format(mcformat_data);
    mcformat_data.Name=NULL;
  }
  if (!mcformat_data.Name && strstr(mcformat.Name, "HTML"))
    mcformat_data = mcuse_format("VRML");

/* *** install sig handler, but only once !! after parameters parsing ******* */
#ifndef NOSIGNALS
#ifdef SIGQUIT
  if (signal( SIGQUIT ,sighandler) == SIG_IGN)
    signal( SIGQUIT,SIG_IGN);   /* quit (ASCII FS) */
#endif
#ifdef SIGABRT
  if (signal( SIGABRT ,sighandler) == SIG_IGN)
    signal( SIGABRT,SIG_IGN);   /* used by abort, replace SIGIOT in the future */
#endif
#ifdef SIGTERM
  if (signal( SIGTERM ,sighandler) == SIG_IGN)
    signal( SIGTERM,SIG_IGN);   /* software termination signal from kill */
#endif
#ifdef SIGUSR1
  if (signal( SIGUSR1 ,sighandler) == SIG_IGN)
    signal( SIGUSR1,SIG_IGN);   /* display simulation status */
#endif
#ifdef SIGUSR2
  if (signal( SIGUSR2 ,sighandler) == SIG_IGN)
    signal( SIGUSR2,SIG_IGN);
#endif
#ifdef SIGHUP
  if (signal( SIGHUP ,sighandler) == SIG_IGN)
    signal( SIGHUP,SIG_IGN);
#endif
#ifdef SIGILL
  if (signal( SIGILL ,sighandler) == SIG_IGN)
    signal( SIGILL,SIG_IGN);    /* illegal instruction (not reset when caught) */
#endif
#ifdef SIGFPE
  if (signal( SIGFPE ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* floating point exception */
#endif
#ifdef SIGBUS
  if (signal( SIGBUS ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* bus error */
#endif
#ifdef SIGSEGV
  if (signal( SIGSEGV ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);   /* segmentation violation */
#endif
#endif /* !NOSIGNALS */
  if (!strstr(mcformat.Name,"NeXus")) {
    mcsiminfo_init(NULL); mcsiminfo_close();  /* makes sure we can do that */
  }
  SIG_MESSAGE("main (Init)");
  mcinit();
#ifndef NOSIGNALS
#ifdef SIGINT
  if (signal( SIGINT ,sighandler) == SIG_IGN)
    signal( SIGINT,SIG_IGN);    /* interrupt (rubout) only after INIT */
#endif
#endif /* !NOSIGNALS */

/* ================ main particule generation/propagation loop ================ */
#if defined (USE_MPI)
  /* sliced Ncount on each MPI node */
  mcncount = mpi_node_count > 1 ?
    floor(mcncount / mpi_node_count) :
    mcncount; /* number of rays per node */
#endif

/* main particule event loop */
while(mcrun_num < mcncount || mcrun_num < mcget_ncount())
  {
#ifndef NEUTRONICS
    mcgenstate();
#endif
    /* old init: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
    mcraytrace();
    mcrun_num++;
  }

#ifdef USE_MPI
 /* merge run_num from MPI nodes */
  if (mpi_node_count > 1) {
  double mcrun_num_double = (double)mcrun_num;
  mc_MPI_Sum(&mcrun_num_double, 1);
  mcrun_num = (unsigned long long)mcrun_num_double;
  }
#endif

/* save/finally executed by master node/thread */
  mcfinally();
  mcclear_format(mcformat);
  if (mcformat_data.Name) mcclear_format(mcformat_data);

#ifdef USE_MPI
  MPI_Finalize();
#endif /* USE_MPI */

  return 0;
} /* mccode_main */

#ifdef NEUTRONICS
/*Main neutronics function steers the McStas calls, initializes parameters etc */
/* Only called in case NEUTRONICS = TRUE */
void neutronics_main_(float *inx, float *iny, float *inz, float *invx, float *invy, float *invz, float *intime, float *insx, float *insy, float *insz, float *inw, float *outx, float *outy, float *outz, float *outvx, float *outvy, float *outvz, float *outtime, float *outsx, float *outsy, float *outsz, float *outwgt)
{

  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  /* External code governs iteration - McStas is iterated once per call to neutronics_main. I.e. below counter must be initiancated for each call to neutronics_main*/
  mcrun_num=0;

  time_t t;
  t = (time_t)mcstartdate;
  mcstartdate = t;  /* set start date before parsing options and creating sim file */
  mcinit();

  /* *** parse options *** */
  SIG_MESSAGE("main (Start)");
  mcformat=mcuse_format(getenv(FLAVOR_UPPER "_FORMAT") ?
                        getenv(FLAVOR_UPPER "_FORMAT") : MCSTAS_FORMAT);
  /* default is to output as McStas format */
  mcformat_data.Name=NULL;
  /* read simulation parameters and options */
  if (strstr(mcformat.Name, "NeXus")) {
    if (mcformat_data.Name) mcclear_format(mcformat_data);
    mcformat_data.Name=NULL;
  }
  if (!mcformat_data.Name && strstr(mcformat.Name, "HTML"))
    mcformat_data = mcuse_format("VRML");

  /* Set neutron state based on input from neutronics code */
  mcsetstate(*inx,*iny,*inz,*invx,*invy,*invz,*intime,*insx,*insy,*insz,*inw);

  /* main neutron event loop - runs only one iteration */

  //mcstas_raytrace(&mcncount); /* prior to McStas 1.12 */

  mcallowbackprop = 1; //avoid absorbtion from negative dt
  int argc=1;
  char *argv[0];
  int dummy = mccode_main(argc, argv);

  /*  save/finally executed by master node/thread. Needed for McStas 1.12 and older */
  /*  mcfinally(); */
  /*  mcclear_format(mcformat); */

  if (mcformat_data.Name) mcclear_format(mcformat_data);

  *outx =  mcnx;
  *outy =  mcny;
  *outz =  mcnz;
  *outvx =  mcnvx;
  *outvy =  mcnvy;
  *outvz =  mcnvz;
  *outtime =  mcnt;
  *outsx =  mcnsx;
  *outsy =  mcnsy;
  *outsz =  mcnsz;
  *outwgt =  mcnp;

  return;
} /* neutronics_main */

#endif /*NEUTRONICS*/

#ifdef USE_NEXUS
/*******************************************************************************
* mcnxfile_parameters: writes the simulation parameters
*                   open/close a new Data Set per parameter in the current simulation Group
* NOTE: this function can not be included in nexus-lib as it depends on mcinputtypes
*       and mcinputtable are defined at compile time in here.
* Returns: NX_OK
*******************************************************************************/
int mcnxfile_parameters(NXhandle nxhandle)
{
  int i;
  char Parameters[CHAR_BUF_LENGTH];
  /* in the parameter group, create one SDS per parameter */
  for(i = 0; i < mcnumipar; i++) {
    if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
      char nxname[CHAR_BUF_LENGTH];
      int length;
      if (mcinputtable[i].par == NULL)
        strncpy(Parameters, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
      else
        (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
      sprintf(nxname, "%s", Parameters);
      length = strlen(nxname);
      NXmakedata(mcnxHandle, mcinputtable[i].name, NX_CHAR, 1, &length);
      NXopendata(mcnxHandle, mcinputtable[i].name);
      NXputdata (mcnxHandle, nxname);
      strcpy(nxname, (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name));
      NXputattr (mcnxHandle, "type", nxname, strlen(nxname), NX_CHAR);
      strcpy(nxname, mcinputtable[i].val ? mcinputtable[i].val : "");
      NXputattr (mcnxHandle, "default_value", nxname, strlen(nxname), NX_CHAR);
      NXputattr (mcnxHandle, "name", mcinputtable[i].name, strlen(mcinputtable[i].name), NX_CHAR);
      NXclosedata(mcnxHandle);
    }
  }
  return(NX_OK);
} /* mcnxfile_parameters */
#endif /* USE_NEXUS */

#endif /* !MCCODE_H */
/* End of file "mccode-r.c". */
/* End of file "mccode-r.c". */

#line 6385 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

#line 1 "mcstas-r.c"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision: 1.229 $
*
* Runtime system for McStas.
* Embedded within instrument in runtime mode.
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id: mcstas-r.c,v 1.229 2009/08/13 14:52:01 farhi Exp $
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#include "mcstas-r.h"
#endif
#ifdef DANSE
#include "mcstas-globals.h"
#endif

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/

/*the magnet stack*/
#ifdef MC_POL_COMPAT
void (*mcMagnetPrecession) (double, double, double, double, double, double,
    double, double*, double*, double*, double, Coords, Rotation)=NULL;
Coords   mcMagnetPos;
Rotation mcMagnetRot;
double*  mcMagnetData                = NULL;
/* mcMagneticField(x, y, z, t, Bx, By, Bz) */
int (*mcMagneticField) (double, double, double, double,
    double*, double*, double*, void *) = NULL;
#endif

#ifndef MCSTAS_H

/*******************************************************************************
* mcstore_neutron: stores neutron coodinates into global array (per component)
*******************************************************************************/
void
mcstore_neutron(MCNUM *s, int index, double x, double y, double z,
               double vx, double vy, double vz, double t,
               double sx, double sy, double sz, double p)
{
    double *dptr = &s[11*index];
    *dptr++  = x;
    *dptr++  = y ;
    *dptr++  = z ;
    *dptr++  = vx;
    *dptr++  = vy;
    *dptr++  = vz;
    *dptr++  = t ;
    *dptr++  = sx;
    *dptr++  = sy;
    *dptr++  = sz;
    *dptr    = p ;
} /* mcstore_neutron */

/*******************************************************************************
* mcrestore_neutron: restores neutron coodinates from global array
*******************************************************************************/
void
mcrestore_neutron(MCNUM *s, int index, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *t,
               double *sx, double *sy, double *sz, double *p)
{
    double *dptr = &s[11*index];
    *x  =  *dptr++;
    *y  =  *dptr++;
    *z  =  *dptr++;
    *vx =  *dptr++;
    *vy =  *dptr++;
    *vz =  *dptr++;
    *t  =  *dptr++;
    *sx =  *dptr++;
    *sy =  *dptr++;
    *sz =  *dptr++;
    *p  =  *dptr;
} /* mcrestore_neutron */

/*******************************************************************************
* mcsetstate: transfer parameters into global McStas variables 
*******************************************************************************/
void
mcsetstate(double x, double y, double z, double vx, double vy, double vz,
           double t, double sx, double sy, double sz, double p)
{
  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  mcnx = x;
  mcny = y;
  mcnz = z;
  mcnvx = vx;
  mcnvy = vy;
  mcnvz = vz;
  mcnt = t;
  mcnsx = sx;
  mcnsy = sy;
  mcnsz = sz;
  mcnp = p;
} /* mcsetstate */

/*******************************************************************************
* mcgenstate: set default neutron parameters 
*******************************************************************************/
void
mcgenstate(void)
{
  mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
  /* old initialisation: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
}

/* intersection routines ==================================================== */

/*******************************************************************************
* inside_rectangle: Check if (x,y) is inside rectangle (xwidth, yheight) 
* return 0 if outside and 1 if inside 
*******************************************************************************/
int inside_rectangle(double x, double y, double xwidth, double yheight)
{
  if (x>-xwidth/2 && x<xwidth/2 && y>-yheight/2 && y<yheight/2)
    return 1;
  else
    return 0;
}

/*******************************************************************************
 * box_intersect: compute time intersection with a box
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times dt_in and dt_out
 * This function written by Stine Nyborg, 1999. 
 *******************************************************************************/
int box_intersect(double *dt_in, double *dt_out,
                  double x, double y, double z,
                  double vx, double vy, double vz,
                  double dx, double dy, double dz)
{
  double x_in, y_in, z_in, tt, t[6], a, b;
  int i, count, s;

      /* Calculate intersection time for each of the six box surface planes
       *  If the box surface plane is not hit, the result is zero.*/

  if(vx != 0)
   {
    tt = -(dx/2 + x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[0] = tt;
    else
      t[0] = 0;

    tt = (dx/2 - x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[1] = tt;
    else
      t[1] = 0;
   }
  else
    t[0] = t[1] = 0;

  if(vy != 0)
   {
    tt = -(dy/2 + y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[2] = tt;
    else
      t[2] = 0;

    tt = (dy/2 - y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[3] = tt;
    else
      t[3] = 0;
   }
  else
    t[2] = t[3] = 0;

  if(vz != 0)
   {
    tt = -(dz/2 + z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[4] = tt;
    else
      t[4] = 0;

    tt = (dz/2 - z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[5] = tt;
    else
      t[5] = 0;
   }
  else
    t[4] = t[5] = 0;

  /* The intersection is evaluated and *dt_in and *dt_out are assigned */

  a = b = s = 0;
  count = 0;

  for( i = 0; i < 6; i = i + 1 )
    if( t[i] == 0 )
      s = s+1;
    else if( count == 0 )
    {
      a = t[i];
      count = 1;
    }
    else
    {
      b = t[i];
      count = 2;
    }

  if ( a == 0 && b == 0 )
    return 0;
  else if( a < b )
  {
    *dt_in = a;
    *dt_out = b;
    return 1;
  }
  else
  {
    *dt_in = b;
    *dt_out = a;
    return 1;
  }

} /* box_intersect */

/*******************************************************************************
 * cylinder_intersect: compute intersection with a cylinder
 * returns 0 when no intersection is found
 *      or 2/4/8/16 bits depending on intersection,
 *     and resulting times t0 and t1
 * Written by: EM,NB,ABA 4.2.98 
  *******************************************************************************/
int
cylinder_intersect(double *t0, double *t1, double x, double y, double z,
                   double vx, double vy, double vz, double r, double h)
{
  double D, t_in, t_out, y_in, y_out;
  int ret=1;

  D = (2*vx*x + 2*vz*z)*(2*vx*x + 2*vz*z)
    - 4*(vx*vx + vz*vz)*(x*x + z*z - r*r);

  if (D>=0)
  {
    if (vz*vz + vx*vx) {
      t_in  = (-(2*vz*z + 2*vx*x) - sqrt(D))/(2*(vz*vz + vx*vx));
      t_out = (-(2*vz*z + 2*vx*x) + sqrt(D))/(2*(vz*vz + vx*vx));
    } else if (vy) { /* trajectory parallel to cylinder axis */
      t_in = (-h/2-y)/vy;
      t_out = (h/2-y)/vy;
      if (t_in>t_out){
        double tmp=t_in;
        t_in=t_out;t_out=tmp;
      }
    } else return 0;
    y_in = vy*t_in + y;
    y_out =vy*t_out + y;

    if ( (y_in > h/2 && y_out > h/2) || (y_in < -h/2 && y_out < -h/2) )
      return 0;
    else
    {
      if (y_in > h/2)
        { t_in = ((h/2)-y)/vy; ret += 2; }
      else if (y_in < -h/2)
        { t_in = ((-h/2)-y)/vy; ret += 4; }
      if (y_out > h/2)
        { t_out = ((h/2)-y)/vy; ret += 8; }
      else if (y_out < -h/2)
        { t_out = ((-h/2)-y)/vy; ret += 16; }
    }
    *t0 = t_in;
    *t1 = t_out;
    return ret;
  }
  else
  {
    *t0 = *t1 = 0;
    return 0;
  }
} /* cylinder_intersect */


/*******************************************************************************
 * sphere_intersect: Calculate intersection between a line and a sphere.
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times t0 and t1 
 *******************************************************************************/
int
sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r)
{
  double A, B, C, D, v;

  v = sqrt(vx*vx + vy*vy + vz*vz);
  A = v*v;
  B = 2*(x*vx + y*vy + z*vz);
  C = x*x + y*y + z*z - r*r;
  D = B*B - 4*A*C;
  if(D < 0)
    return 0;
  D = sqrt(D);
  *t0 = (-B - D) / (2*A);
  *t1 = (-B + D) / (2*A);
  return 1;
} /* sphere_intersect */

/*******************************************************************************
 * plane_intersect: Calculate intersection between a plane and a line.
 * returns 0 when no intersection is found (i.e. line is parallel to the plane)
 * returns 1 or -1 when intersection time is positive and negative respectively
 *******************************************************************************/
int
plane_intersect(double *t, double x, double y, double z,
                 double vx, double vy, double vz, double nx, double ny, double nz, double wx, double wy, double wz)
{
  double s;
  if (fabs(s=scalar_prod(nx,ny,nz,vx,vy,vz))<FLT_EPSILON) return 0;
  *t = - scalar_prod(nx,ny,nz,x-wx,y-wy,z-wz)/s;
  if (*t<0) return -1;
  else return 1;
} /* plane_intersect */

#endif /* !MCSTAS_H */
/* End of file "mcstas-r.c". */

#line 6745 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/usr/local/lib/mcstas-2.0/"
int mcdefaultmain = 1;
char mcinstrument_name[] = "test";
char mcinstrument_source[] = "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr";
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Virtual_input'. */
#line 82 "/usr/local/lib/mcstas-2.0/sources/Virtual_input.comp"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision: 1.21 $
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions.
*
* This library may be used directly as an external library. It has no dependency
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#define READ_TABLE_LIB_H "$Revision: 1.21 $"

#define READ_TABLE_STEPTOL  0.04 /* tolerancy for constant step approx */

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#ifdef MAC
#define MC_PATHSEP_C ':'
#define MC_PATHSEP_S ":"
#else  /* !MAC */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MC_PATHSEP_C */

#ifndef MCSTAS
#ifdef WIN32
#define MCSTAS "C:\\mcstas\\lib"
#else  /* !WIN32 */
#ifdef MAC
#define MCSTAS ":mcstas:lib" /* ToDo: What to put here? */
#else  /* !MAC */
#define MCSTAS "/usr/local/lib/mcstas"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MCSTAS */

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

  typedef struct struct_table
  {
    char    filename[256];
    long    filesize;
    char   *header;  /* text header, e.g. comments */
    double *data;    /* vector { x[0], y[0], ... x[n-1], y[n-1]... } */
    double  min_x;   /* min value of first column */
    double  max_x;   /* max value of first column */
    double  step_x;  /* minimal step value of first column */
    long    rows;    /* number of rows in matrix block */
    long    columns; /* number of columns in matrix block */

    long    begin;   /* start fseek index of block */
    long    end;     /* stop  fseek index of block */
    long    block_number;  /* block index. 0 is catenation of all */
    long    array_length;  /* number of elements in the t_Table array */
    char    monotonic;     /* true when 1st column/vector data is monotonic */
    char    constantstep;  /* true when 1st column/vector data has constant step */
    char    method[32];    /* interpolation method: nearest, linear */
  } t_Table;

/* read_table-lib function prototypes */
/* ========================================================================= */

/* 'public' functions */
long     Table_Read              (t_Table *Table, char *File, long block_number);
long     Table_Read_Offset       (t_Table *Table, char *File, long block_number,
                                  long *offset, long max_lines);
long     Table_Read_Offset_Binary(t_Table *Table, char *File, char *Type,
                                  long *Offset, long Rows, long Columns);
long     Table_Rebin(t_Table *Table); /* rebin table with regular 1st column and interpolate all columns 2:end */
long     Table_Info (t_Table Table);
double   Table_Index(t_Table Table,   long i, long j); /* get indexed value */
double   Table_Value(t_Table Table, double X, long j); /* search X in 1st column and return interpolated value in j-column */
t_Table *Table_Read_Array(char *File, long *blocks);
void     Table_Free_Array(t_Table *Table);
long     Table_Info_Array(t_Table *Table);
int      Table_SetElement(t_Table *Table, long i, long j, double value);
long     Table_Init(t_Table *Table, long rows, long columns); /* create a Table */
double   Table_Value2d(t_Table Table, double X, double Y);    /* same as Table_Index with non-integer indices and 2d interpolation */
MCDETECTOR Table_Write(t_Table Table, char*file, char*xl, char*yl, 
           double x1, double x2, double y1, double y2); /* write Table to disk */

#define Table_ParseHeader(header, ...) \
  Table_ParseHeader_backend(header,__VA_ARGS__,NULL);

char **Table_ParseHeader_backend(char *header, ...);

/* private functions */
void Table_Free(t_Table *Table);
long Table_Read_Handle(t_Table *Table, FILE *fid, long block_number, long max_lines, char *name);
static void Table_Stat(t_Table *Table);
double Table_Interp1d(double x, double x1, double y1, double x2, double y2);
double Table_Interp1d_nearest(double x, double x1, double y1, double x2, double y2);
double Table_Interp2d(double x, double y, double x1, double y1, double x2, double y2,
  double z11, double z12, double z21, double z22);

#endif

/* end of read_table-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas CVS_090504
* Version: $Revision: 1.42 $
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#endif

/*******************************************************************************
* long Read_Table(t_Table *Table, char *name, int block_number)
*   ACTION: read a single Table from a text file
*   input   Table: pointer to a t_Table structure
*           name:  file name from which table should be extracted
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* File is opened, read and closed
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebined with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read(t_Table *Table, char *File, long block_number)
  { /* reads all or a single data block from 'file' and returns a Table structure  */
    return(Table_Read_Offset(Table, File, block_number, NULL, 0));
  } /* end Table_Read */

/*******************************************************************************
* long Table_Read_Offset(t_Table *Table, char *name, int block_number, long *offset
*                        long max_rows)
*   ACTION: read a single Table from a text file, starting at offset
*     Same as Table_Read(..) except:
*   input   offset:    pointer to an offset (*offset should be 0 at start)
*           max_rows: max number of data rows to read from file (0 means all)
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset(t_Table *Table, char *File,
                         long block_number, long *offset,
                         long max_rows)
  { /* reads all/a data block in 'file' and returns a Table structure  */
    FILE *hfile;
    long  nelements;
    long  begin;
    long  filesize=0;
    char  name[256];
    char  path[1024];
    struct stat stfile;

    if (!Table) return(-1);
    Table_Init(Table, 0, 0);
    if (!File || File[0]=='\0')  return(-1);
    if (!strcmp(File,"NULL") || !strcmp(File,"0"))  return(-1);
    strncpy(path, File, 1024);
    hfile = fopen(path, "r");
    if(!hfile)
    {
      char dir[1024];

      if (!hfile) /* search in instrument location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_source, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_source;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_source, path_length);
            sprintf(path, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, "r");
          }
        }
      }
      if (!hfile) /* search in HOME */
      {
        strcpy(dir, getenv("HOME") ? getenv("HOME") : ".");
        sprintf(path, "%s%c%s", dir, MC_PATHSEP_C, File);
        hfile = fopen(path, "r");
      }
      if (!hfile) /* search in MCSTAS data */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        sprintf(path, "%s%c%s%c%s", dir, MC_PATHSEP_C, "data", MC_PATHSEP_C, File);
        hfile = fopen(path, "r");
      }
      if (!hfile) /* search in MVCSTAS/contrib */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        sprintf(path, "%s%c%s%c%s", dir, MC_PATHSEP_C, "contrib", MC_PATHSEP_C, File);
        hfile = fopen(path, "r");
      }
      if(!hfile)
      {
        fprintf(stderr, "Error: Could not open input file '%s' (Table_Read_Offset_Binary)\n", File);
        return (-1);
      } else if (!offset || (offset && !*offset))
        printf("Opening input file '%s' (Table_Read)\n", path);
    }
    stat(path,&stfile); filesize = stfile.st_size;
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    if (offset && *offset) sprintf(name, "%s@%li", File, *offset);
    else                   strncpy(name, File, 128);
    nelements = Table_Read_Handle(Table, hfile, block_number, max_rows, name);
    Table->begin = begin;
    Table->end   = ftell(hfile);
    Table->filesize = (filesize>0 ? filesize : 0);
    Table_Stat(Table);
    if (offset) *offset=Table->end;
    fclose(hfile);
    return(nelements);

  } /* end Table_Read_Offset */

/*******************************************************************************
* long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
*                               long *offset, long rows, long columns)
*   ACTION: read a single Table from a binary file, starting at offset
*     Same as Table_Read_Offset(..) except that it handles binary files.
*   input   type: may be "float"/NULL or "double"
*           offset: pointer to an offset (*offset should be 0 at start)
*           rows   : number of rows (0 means read all)
*           columns: number of columns
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
                                long *offset, long rows, long columns)
  { /* reads all/a data block in binary 'file' and returns a Table structure  */
    long    nelements, sizeofelement;
    long    filesize;
    FILE   *hfile;
    struct stat stfile;
    double *data;
    long    i;
    long    begin;

    if (!Table) return(-1);

    Table_Init(Table, 0, 0);
    if (!File || File[0]=='\0')  return(-1);
    if (!strcmp(File,"NULL") || !strcmp(File,"0"))  return(-1);

    hfile = fopen(File, "r");
    if(!hfile)
    {
      char path[1024];
      char dir[1024];

      if (!hfile) /* search in instrument location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_source, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_source;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_source, path_length);
            sprintf(path, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, "r");
          }
        }
      }
      if (!hfile) /* search in HOME */
      {
        strcpy(dir, getenv("HOME") ? getenv("HOME") : ".");
        sprintf(path, "%s%c%s", dir, MC_PATHSEP_C, File);
        hfile = fopen(path, "r");
      }
      if (!hfile) /* search in MCSTAS data */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        sprintf(path, "%s%c%s%c%s", dir, MC_PATHSEP_C, "data", MC_PATHSEP_C, File);
        hfile = fopen(path, "r");
      }
      if (!hfile) /* search in MVCSTAS/contrib */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        sprintf(path, "%s%c%s%c%s", dir, MC_PATHSEP_C, "contrib", MC_PATHSEP_C, File);
        hfile = fopen(path, "r");
      }
      if(!hfile)
      {
        fprintf(stderr, "Error: Could not open input file '%s' (Table_Read_Offset_Binary)\n", File);
        return (-1);
      } else
        printf("Opening input file '%s' (Table_Read)\n", path);
    }
    stat(File,&stfile);
    filesize = stfile.st_size;
    Table->filesize=filesize;
    if (type && !strcmp(type,"double")) sizeofelement = sizeof(double);
    else  sizeofelement = sizeof(float);
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    if (rows && filesize > sizeofelement*columns*rows)
      nelements = columns*rows;
    else nelements = (long)(filesize/sizeofelement);
    if (!nelements || filesize <= *offset) return(0);
    data    = (double*)malloc(nelements*sizeofelement);
    if (!data) {
      fprintf(stderr,"Error: allocating %ld elements for %s file '%s'. Too big (Table_Read_Offset_Binary).\n", nelements, type, File);
      exit(-1);
    }
    nelements = fread(data, sizeofelement, nelements, hfile);

    if (!data || !nelements)
    {
      fprintf(stderr,"Error: reading %ld elements from %s file '%s' (Table_Read_Offset_Binary)\n", nelements, type, File);
      exit(-1);
    }
    Table->begin   = begin;
    Table->end     = ftell(hfile);
    if (offset) *offset=Table->end;
    fclose(hfile);
    data = (double*)realloc(data, (double)nelements*sizeofelement);
    /* copy file data into Table */
    if (type && !strcmp(type,"double")) Table->data = data;
    else {
      float  *s;
      double *dataf;
      s     = (float*)data;
      dataf = (double*)malloc(sizeof(double)*nelements);
      for (i=0; i<nelements; i++)
        dataf[i]=s[i];
      free(data);
      Table->data = dataf;
    }
    strcpy(Table->filename, File);
    Table->rows    = nelements/columns;
    Table->columns = columns;
    Table->array_length = 1;
    Table->block_number = 1;

    Table_Stat(Table);

    return(nelements);
  } /* end Table_Read_Offset_Binary */

/*******************************************************************************
* long Table_Read_Handle(t_Table *Table, FILE *fid, int block_number, long max_rows, char *name)
*   ACTION: read a single Table from a text file handle (private)
*   input   Table:pointer to a t_Table structure
*           fid:  pointer to FILE handle
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*           max_rows: if non 0, only reads that number of lines
*   return  initialized single Table t_Table structure containing data, header, ...
*           modified Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebined with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read_Handle(t_Table *Table, FILE *hfile,
                         long block_number, long max_rows, char *name)
  { /* reads all/a data block from 'file' handle and returns a Table structure  */
    double *Data;
    char *Header              = NULL;
    long  malloc_size         = CHAR_BUF_LENGTH;
    long  malloc_size_h       = 4096;
    long  Rows = 0,   Columns = 0;
    long  count_in_array      = 0;
    long  count_in_header     = 0;
    long  block_Current_index = 0;
    char  flag_End_row_loop   = 0;

    if (!Table) return(-1);
    Table_Init(Table, 0, 0);
    if (name && name[0]!='\0') strcpy(Table->filename, name);

    if(!hfile) {
       fprintf(stderr, "Error: File handle is NULL (Table_Read_Handle).\n");
       return (-1);
    }
    Header = (char*)  calloc(malloc_size_h, sizeof(char));
    Data   = (double*)calloc(malloc_size,   sizeof(double));
    if ((Header == NULL) || (Data == NULL)) {
       fprintf(stderr, "Error: Could not allocate Table and Header (Table_Read_Handle).\n");
       return (-1);
    }

    int flag_In_array = 0;
    do { /* while (!flag_End_row_loop) */
      char  line[1024*CHAR_BUF_LENGTH];
      long  back_pos=0;   /* ftell start of line */

      back_pos = ftell(hfile);
      if (fgets(line, 1024*CHAR_BUF_LENGTH, hfile) != NULL) { /* analyse line */
        /* first skip blank and tabulation characters */
        int i = strspn(line, " \t");

        /* handle comments: stored in header */
        if (NULL != strchr("#%;/", line[i]))
        { /* line is a comment */
          count_in_header += strlen(line);
          if (count_in_header >= malloc_size_h) {
            /* if succeed and in array : add (and realloc if necessary) */
            malloc_size_h = count_in_header+4096;
            Header        = (char*)realloc(Header, malloc_size_h*sizeof(char));
          }
          strncat(Header, line, 4096);
          flag_In_array=0;
          /* exit line and file if passed desired block */
          if (block_number > 0 && block_number == block_Current_index) {
            flag_End_row_loop = 1;
          }

          /* Continue with next line */
          continue;
        }

        /* get the number of columns splitting line with strtok */
        char  *lexeme;
        char  flag_End_Line = 0;
        long  block_Num_Columns = 0;
        const char seps[] = " ,;\t\n\r";

        lexeme = strtok(line, seps);
        while (!flag_End_Line) {
          if ((lexeme != NULL) && (lexeme[0] != '\0')) {
            /* reading line: the token is not empty */
            double X;
            if (sscanf(lexeme,"%lg",&X) == 1) {
              /* reading line: the token is a number in the line */
              if (!flag_In_array) {
                /* reading num: not already in a block: starts a new data block */
                block_Current_index++;
                flag_In_array    = 1;
                block_Num_Columns= 0;
                if (block_number > 0) {
                  /* initialise a new data block */
                  Rows = 0;
                  count_in_array = 0;
                } /* else append */
              }
              /* reading num: all blocks or selected block */
              if (flag_In_array && (block_number == 0 ||
                  block_number == block_Current_index)) {
                /* starting block: already the desired number of rows ? */
                if (block_Num_Columns == 0 &&
                    max_rows > 0 && Rows >= max_rows) {
                  flag_End_Line      = 1;
                  flag_End_row_loop  = 1;
                  flag_In_array      = 0;
                  /* reposition to begining of line (ignore line) */
                  fseek(hfile, back_pos, SEEK_SET);
                } else { /* store into data array */
                  if (count_in_array >= malloc_size) {
                    /* realloc data buffer if necessary */
                    malloc_size = count_in_array+CHAR_BUF_LENGTH;
                    Data = (double*) realloc(Data, malloc_size*sizeof(double));
                    if (Data == NULL) {
                      fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Handle).\n",
                              malloc_size*sizeof(double));
                      return (-1);
                    }
                  }
                  if (0 == block_Num_Columns) Rows++;
                  Data[count_in_array] = X;
                  count_in_array++;
                  block_Num_Columns++;
                }
              } /* reading num: end if flag_In_array */
            } /* end reading num: end if sscanf lexeme -> numerical */
            else {
              /* reading line: the token is not numerical in that line. end block */
              if (block_Current_index == block_number) {
                flag_End_Line = 1;
                flag_End_row_loop = 1;
              } else {
                flag_In_array = 0;
                flag_End_Line = 1;
              }
            }
          }
          else {
            /* no more tokens in line */
            flag_End_Line = 1;
            if (block_Num_Columns > 0) Columns = block_Num_Columns;
          }

          // parse next token
          lexeme = strtok(NULL, seps);

        } /* while (!flag_End_Line) */
      } /* end: if fgets */
      else flag_End_row_loop = 1; /* else fgets : end of file */

    } while (!flag_End_row_loop); /* end while flag_End_row_loop */

    Table->block_number = block_number;
    Table->array_length = 1;

    // shrink header to actual size (plus terminating 0-byte)
    if (count_in_header) {
      Header = (char*)realloc(Header, count_in_header*sizeof(char) + 1);
    }
    Table->header = Header;

    if (count_in_array*Rows*Columns == 0)
    {
      Table->rows         = 0;
      Table->columns      = 0;
      free(Data);
      return (0);
    }
    if (Rows * Columns != count_in_array)
    {
      fprintf(stderr, "Warning: Read_Table :%s %s Data has %li values that should be %li x %li\n",
        (Table->filename ? Table->filename : ""),
        (!block_number ? " catenated" : ""),
        count_in_array, Rows, Columns);
      Columns = count_in_array; Rows = 1;
    }
    Data     = (double*)realloc(Data, count_in_array*sizeof(double));
    Table->data         = Data;
    Table->rows         = Rows;
    Table->columns      = Columns;

    return (count_in_array);

  } /* end Table_Read_Handle */

/*******************************************************************************
* long Table_Rebin(t_Table *Table)
*   ACTION: rebin a single Table, sorting 1st column in ascending order
*   input   Table: single table containing data.
*                  The data block is reallocated in this process
*   return  updated Table with increasing, evenly spaced first column (index 0)
*           number of data elements (-1: error, 0:empty data)
*******************************************************************************/
  long Table_Rebin(t_Table *Table)
  {
    double new_step=0;
    long   i;
    /* performs linear interpolation on X axis (0-th column) */

    if (!Table) return(-1);
    if (!Table->data
    || Table->rows*Table->columns == 0 || !Table->step_x)
      return(0);
    Table_Stat(Table); /* recompute statitstics and minimal step */
    new_step = Table->step_x; /* minimal step in 1st column */

    if (!(Table->constantstep)) /* not already evenly spaced */
    {
      long Length_Table;
      double *New_Table;

      Length_Table = ceil(fabs(Table->max_x - Table->min_x)/new_step);
      New_Table    = (double*)malloc(Length_Table*Table->columns*sizeof(double));

      for (i=0; i < Length_Table; i++)
      {
        long   j;
        double X;
        X = Table->min_x + i*new_step;
        New_Table[i*Table->columns] = X;
        for (j=1; j < Table->columns; j++)
          New_Table[i*Table->columns+j]
                = Table_Value(*Table, X, j);
      } /* end for i */

      Table->rows = Length_Table;
      Table->step_x = new_step;
      free(Table->data);
      Table->data = New_Table;
    } /* end else (!constantstep) */
    return (Table->rows*Table->columns);
  } /* end Table_Rebin */

/*******************************************************************************
* double Table_Index(t_Table Table, long i, long j)
*   ACTION: read an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*   return  Value = data[i][j]
* Returns Value from the i-th row, j-th column of Table
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

double Table_Index(t_Table Table, long i, long j)
{
  long AbsIndex;

  if (Table.rows == 1 || Table.columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table.columns*Table.rows - 1);
    i = 0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table.rows - 1);
    j = MIN(MAX(0, j), Table.columns - 1);
  }

  /* handle vectors specifically */
  AbsIndex = i*(Table.columns)+j;

  if (Table.data != NULL)
    return (Table.data[AbsIndex]);
  else
    return 0;
} /* end Table_Index */

/*******************************************************************************
* void Table_SetElement(t_Table *Table, long i, long j, double value)
*   ACTION: set an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*           value = data[i][j]
* Returns 0 in case of error
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/
int Table_SetElement(t_Table *Table, long i, long j,
                     double value)
{
  long AbsIndex;

  if (Table->rows == 1 || Table->columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table->columns*Table->rows - 1); i=0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table->rows - 1);
    j = MIN(MAX(0, j), Table->columns - 1);
  }

  AbsIndex = i*(Table->columns)+j;
  if (Table->data != NULL) {
    Table->data[AbsIndex] = value;
    return 1;
  }

  return 0;
} /* end Table_SetElement */

/*******************************************************************************
* double Table_Value(t_Table Table, double X, long j)
*   ACTION: read column [j] of a single Table at row which 1st column is X
*   input   Table: table containing data.
*           X : data value in the first column (index 0)
*           j : index of column from which is extracted the Value (0:Columns-1)
*   return  Value = data[index for X][j] with linear interpolation
* Returns Value from the j-th column of Table corresponding to the
* X value for the 1st column (index 0)
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
double Table_Value(t_Table Table, double X, long j)
{
  long   Index = -1;
  double X1=0, Y1=0, X2=0, Y2=0;
  double ret=0;

  if (X > Table.max_x) return Table_Index(Table,Table.rows-1  ,j);
  if (X < Table.min_x) return Table_Index(Table,0  ,j);

  // Use constant-time lookup when possible
  if(Table.constantstep) {
    Index = (long)floor(
              (X - Table.min_x) / (Table.max_x - Table.min_x) * Table.rows);
    X1 = Table_Index(Table,Index  ,0);
    X2 = Table_Index(Table,Index+1,0);
  }
  // Use binary search on large, monotonic tables
  else if(Table.monotonic && Table.rows > 100) {
    long left = Table.min_x;
    long right = Table.max_x;

    while (!((X1 <= X) && (X < X2)) && (right - left > 1)) {
      Index = (left + right) / 2;

      X1 = Table_Index(Table, Index-1, 0);
      X2 = Table_Index(Table, Index,   0);

      if (X < X1) {
        right = Index;
      } else {
        left  = Index;
      }
    }
  }

  // Fall back to linear search, if no-one else has set X1, X2 correctly
  if (!((X1 <= X) && (X < X2))) {
    /* look for index surrounding X in the table -> Index */
    for (Index=1; Index < Table.rows-1; Index++) {
        X1 = Table_Index(Table, Index-1,0);
        X2 = Table_Index(Table, Index  ,0);
        if ((X1 <= X) && (X < X2)) break;
      } /* end for Index */
  }

  Y1 = Table_Index(Table,Index-1,j);
  Y2 = Table_Index(Table,Index  ,j);

  if (!strcmp(Table.method,"linear")) {
    ret = Table_Interp1d(X, X1,Y1, X2,Y2);
  }
  else if (!strcmp(Table.method,"nearest")) {
    ret = Table_Interp1d_nearest(X, X1,Y1, X2,Y2);
  }

  return ret;
} /* end Table_Value */

/*******************************************************************************
* double Table_Value2d(t_Table Table, double X, double Y)
*   ACTION: read element [X,Y] of a matrix Table
*   input   Table: table containing data.
*           X : row index, may be non integer
*           Y : column index, may be non integer
*   return  Value = data[index X][index Y] with bi-linear interpolation
* Returns Value for the indices [X,Y]
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
  double Table_Value2d(t_Table Table, double X, double Y)
  {
    long   x1,x2,y1,y2;
    double z11,z12,z21,z22;
    double ret=0;

    x1 = (long)floor(X);
    y1 = (long)floor(Y);

    if (x1 > Table.rows-1 || x1 < 0) {
      x2 = x1;
    } else {
      x2 = x1 + 1;
    }

    if (y1 > Table.columns-1 || y1 < 0) {
      y2 = y1;
    } else {
      y2 = y1 + 1;
    }

    z11 = Table_Index(Table, x1, y1);

    if (y2 != y1) z12=Table_Index(Table, x1, y2); else z12 = z11;
    if (x2 != x1) z21=Table_Index(Table, x2, y1); else z21 = z11;
    if (y2 != y1) z22=Table_Index(Table, x2, y2); else z22 = z21;

    if (!strcmp(Table.method,"linear"))
      ret = Table_Interp2d(X,Y, x1,y1,x2,y2, z11,z12,z21,z22);
    else {
      if (fabs(X-x1) < fabs(X-x2)) {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z11; else ret = z12;
      } else {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z21; else ret = z22;
      }
    }
    return ret;
  } /* end Table_Value2d */


/*******************************************************************************
* void Table_Free(t_Table *Table)
*   ACTION: free a single Table
*   return: empty Table
*******************************************************************************/
  void Table_Free(t_Table *Table)
  {
    if (!Table) return;
    if (Table->data   != NULL) free(Table->data);
    if (Table->header != NULL) free(Table->header);
    Table->data   = NULL;
    Table->header = NULL;
  } /* end Table_Free */

/******************************************************************************
* void Table_Info(t_Table Table)
*    ACTION: print informations about a single Table
*******************************************************************************/
  long Table_Info(t_Table Table)
  {
    char buffer[256];
    long ret=0;

    if (!Table.block_number) strcpy(buffer, "catenated");
    else sprintf(buffer, "block %li", Table.block_number);
    printf("Table from file '%s' (%s)",
      Table.filename && strlen(Table.filename) ? Table.filename : "", buffer);
    if ((Table.data != NULL) && (Table.rows*Table.columns))
    {
      printf(" is %li x %li ", Table.rows, Table.columns);
      if (Table.rows*Table.columns > 1)
        printf("(x=%g:%g)", Table.min_x, Table.max_x);
      else printf("(x=%g) ", Table.min_x);
      ret = Table.rows*Table.columns;
      if (Table.monotonic)    printf(", monotonic");
      if (Table.constantstep) printf(", constant step");
      printf(". interpolation: %s\n", Table.method);
    }
    else printf(" is empty.\n");

    if (Table.header && strlen(Table.header)) {
      char *header;
      int  i;
      header = malloc(80);
      if (!header) return(ret);
      for (i=0; i<80; header[i++]=0);
      strncpy(header, Table.header, 75);
      if (strlen(Table.header) > 75) {
        strcat( header, " ...");
      }
      for (i=0; i<strlen(header); i++)
        if (header[i] == '\n' || header[i] == '\r') header[i] = ';';
      printf("  '%s'\n", header);
      free(header);
    }
    return(ret);
  } /* end Table_Info */

/******************************************************************************
* long Table_Init(t_Table *Table, m, n)
*   ACTION: initialise a Table to empty m by n table
*   return: empty Table
******************************************************************************/
long Table_Init(t_Table *Table, long rows, long columns)
{
  double *data=NULL;
  long   i;

  if (!Table) return(0);

  Table->header  = NULL;
  Table->filename[0]= '\0';
  Table->filesize= 0;
  Table->min_x   = 0;
  Table->max_x   = 0;
  Table->step_x  = 0;
  Table->block_number = 0;
  Table->array_length = 0;
  Table->monotonic    = 0;
  Table->constantstep = 0;
  Table->begin   = 0;
  Table->end     = 0;
  strcpy(Table->method,"linear");

  if (rows*columns >= 1) {
    data    = (double*)malloc(rows*columns*sizeof(double));
    if (data) for (i=0; i < rows*columns; data[i++]=0);
    else {
      fprintf(stderr,"Error: allocating %ld double elements."
                     "Too big (Table_Init).\n", rows*columns);
      rows = columns = 0;
    }
  }
  Table->rows    = (rows >= 1 ? rows : 0);
  Table->columns = (columns >= 1 ? columns : 0);
  Table->data    = data;
  return(Table->rows*Table->columns);
} /* end Table_Init */

/******************************************************************************
* long Table_Write(t_Table Table, char *file, x1,x2, y1,y2)
*   ACTION: write a Table to disk (ascii).
*     when x1=x2=0 or y1=y2=0, the table default limits are used.
*   return: 0=all is fine, non-0: error
*******************************************************************************/
MCDETECTOR Table_Write(t_Table Table, char *file, char *xl, char *yl, 
  double x1, double x2, double y1, double y2)
{
  long    i =0;
  MCDETECTOR detector;

  if ((Table.data == NULL) && (Table.rows*Table.columns)) {
    detector.m = 0;
    return(detector); /* Table is empty - nothing to do */
  }
  if (!x1 && !x2) {
    x1 = Table.min_x;
    x2 = Table.max_x;
  }
  if (!y1 && !y2) {
    y1 = 1;
    y2 = Table.columns;
  }

  /* transfer content of the Table into a 2D detector */
  Coords coords = { 0, 0, 0};

  if (Table.rows == 1 || Table.columns == 1) {
    detector = mcdetector_out_1D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      "x", x1, x2,
                      Table.rows * Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  } else {
    detector = mcdetector_out_2D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      x1, x2, y1, y2,
                      Table.rows, Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  }
  return(detector);
}

/******************************************************************************
* void Table_Stat(t_Table *Table)
*   ACTION: computes min/max/mean step of 1st column for a single table (private)
*   return: updated Table
*******************************************************************************/
  static void Table_Stat(t_Table *Table)
  {
    long   i;
    double max_x, min_x;
    double row=1;
    char   monotonic=1;
    char   constantstep=1;
    double step=0;
    double n;

    if (!Table) return;
    if (!Table->rows || !Table->columns) return;
    if (Table->rows == 1) row=0;
    max_x = -FLT_MAX;
    min_x =  FLT_MAX;
    n     = (row ? Table->rows : Table->columns);
    /* get min and max of first column/vector */
    for (i=0; i < n; i++)
    {
      double X;
      X = (row ? Table_Index(*Table,i  ,0)
                               : Table_Index(*Table,0, i));
      if (X < min_x) min_x = X;
      if (X > max_x) max_x = X;
    } /* for */
    /* test for monotonicity and constant step if the table is an XY or single vector */
    if (n > 1) {
      /* mean step */
      step = (max_x - min_x)/(n-1);
      /* now test if table is monotonic on first column, and get minimal step size */
      for (i=0; i < n-1; i++) {
        double X, diff;;
        X    = (row ? Table_Index(*Table,i  ,0)
                    : Table_Index(*Table,0,  i));
        diff = (row ? Table_Index(*Table,i+1,0)
                    : Table_Index(*Table,0,  i+1)) - X;
        if (fabs(diff) < fabs(step)) step = diff;
        /* change sign ? */
        if ((max_x - min_x)*diff < 0 && monotonic)
          monotonic = 0;
      } /* end for */
      /* now test if steps are constant within READ_TABLE_STEPTOL */
      if(!step){
        /*means there's a disconitnuity -> not constantstep*/
        constantstep=0;
      }else if (monotonic) {
        for (i=0; i < n-1; i++) {
          double X, diff;
          X    = (row ? Table_Index(*Table,i  ,0)
              : Table_Index(*Table,0,  i));
          diff = (row ? Table_Index(*Table,i+1,0)
              : Table_Index(*Table,0,  i+1)) - X;
          if ( fabs(step)*(1+READ_TABLE_STEPTOL) < fabs(diff) ||
                fabs(diff) < fabs(step)*(1-READ_TABLE_STEPTOL) )
          { constantstep = 0; break; }
        }
      }
    }
    Table->step_x= step;
    Table->max_x = max_x;
    Table->min_x = min_x;
    Table->monotonic = monotonic;
    Table->constantstep = constantstep;
  } /* end Table_Stat */

/******************************************************************************
* t_Table *Table_Read_Array(char *File, long *blocks)
*   ACTION: read as many data blocks as available, iteratively from file
*   return: initialized t_Table array, last element is an empty Table.
*           the number of extracted blocks in non NULL pointer *blocks
*******************************************************************************/
  t_Table *Table_Read_Array(char *File, long *blocks)
  {
    t_Table *Table_Array=NULL;
    long offset=0;
    long block_number=0;
    long allocated=256;
    long nelements=1;

    /* fisrt allocate an initial empty t_Table array */
    Table_Array = (t_Table *)malloc(allocated*sizeof(t_Table));
    if (!Table_Array) {
      fprintf(stderr, "Error: Can not allocate memory %li (Table_Read_Array).\n",
         allocated*sizeof(t_Table));
      *blocks = 0;
      return (NULL);
    }

    while (nelements > 0)
    {
      t_Table Table;

      /* access file at offset and get following block */
      nelements = Table_Read_Offset(&Table, File, 1, &offset,0);
      /* if ok, set t_Table block number else exit loop */
      block_number++;
      Table.block_number = block_number;
      /* if t_Table array is not long enough, expand and realocate */
      if (block_number >= allocated-1) {
        allocated += 256;
        Table_Array = (t_Table *)realloc(Table_Array,
           allocated*sizeof(t_Table));
        if (!Table_Array) {
          fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Array).\n",
              allocated*sizeof(t_Table));
          *blocks = 0;
          return (NULL);
        }
      }
      /* store it into t_Table array */
      sprintf(Table.filename, "%s#%li", File, block_number-1);
      Table_Array[block_number-1] = Table;
      /* continues until we find an empty block */
    }
    /* send back number of extracted blocks */
    if (blocks) *blocks = block_number-1;

    /* now store total number of elements in Table array */
    for (offset=0; offset < block_number;
      Table_Array[offset++].array_length = block_number-1);

    return(Table_Array);
  } /* end Table_Read_Array */
/*******************************************************************************
* void Table_Free_Array(t_Table *Table)
*   ACTION: free a Table array
*******************************************************************************/
  void Table_Free_Array(t_Table *Table)
  {
    long index=0;
    if (!Table) return;
    do {
        if (Table[index].data || Table[index].header)
          Table_Free(&Table[index]);
        else index=-1;
    } while (index>= 0);
    free(Table);
  } /* end Table_Free_Array */

/******************************************************************************
* long Table_Info_Array(t_Table *Table)
*    ACTION: print informations about a Table array
*    return: number of elements in the Table array
*******************************************************************************/
  long Table_Info_Array(t_Table *Table)
  {
    long index=0;

    if (!Table) return(-1);
    while (index < Table[index].array_length
       && (Table[index].data || Table[index].header)
       && (Table[index].rows*Table[index].columns) ) {
      Table_Info(Table[index]);
      index++;
    }
    printf("This Table array contains %li elements\n", index);
    return(index);
  } /* end Table_Info_Array */

/******************************************************************************
* char **Table_ParseHeader(char *header, symbol1, symbol2, ..., NULL)
*    ACTION: search for char* symbols in header and return their value or NULL
*            Last argument MUST be NULL
*    return: array of char* with line following each symbol, or NULL if not found
*******************************************************************************/
#ifndef MyNL_ARGMAX
#define MyNL_ARGMAX 50
#endif

char **Table_ParseHeader_backend(char *header, ...){
  va_list ap;
  char exit_flag=0;
  int counter   =0;
  char **ret    =NULL;
  if (!header || header[0]=='\0') return(NULL);

  ret = (char**)calloc(MyNL_ARGMAX, sizeof(char*));
  if (!ret) {
    printf("Table_ParseHeader: Cannot allocate %i values array for Parser (Table_ParseHeader).\n",
      MyNL_ARGMAX);
    return(NULL);
  }
  for (counter=0; counter < MyNL_ARGMAX; ret[counter++] = NULL);
  counter=0;

  va_start(ap, header);
  while(!exit_flag && counter < MyNL_ARGMAX-1)
  {
    char *arg_char=NULL;
    char *pos     =NULL;
    /* get variable argument value as a char */
    arg_char = va_arg(ap, char *);
    if (!arg_char || arg_char[0]=='\0'){
      exit_flag = 1; break;
    }
    /* search for the symbol in the header */
    pos = strstr(header, arg_char);
    if (pos) {
      char *eol_pos;
      eol_pos = strchr(pos+strlen(arg_char), '\n');
      if (!eol_pos)
        eol_pos = strchr(pos+strlen(arg_char), '\r');
      if (!eol_pos)
        eol_pos = pos+strlen(pos)-1;
      ret[counter] = (char*)malloc(eol_pos - pos);
      if (!ret[counter]) {
        printf("Table_ParseHeader: Cannot allocate value[%i] array for Parser searching for %s (Table_ParseHeader).\n",
          counter, arg_char);
        exit_flag = 1; break;
      }
      strncpy(ret[counter], pos+strlen(arg_char), eol_pos - pos - strlen(arg_char));
      ret[counter][eol_pos - pos - strlen(arg_char)]='\0';
    }
    counter++;
  }
  va_end(ap);
  return(ret);
} /* Table_ParseHeader */

/******************************************************************************
* double Table_Interp1d(x, x1, y1, x2, y2)
*    ACTION: interpolates linearly at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d(double x,
  double x1, double y1,
  double x2, double y2)
{
  double slope;
  if (x2 == x1) return (y1+y2)/2;
  if (y1 == y2) return  y1;
  slope = (y2 - y1)/(x2 - x1);
  return y1+slope*(x - x1);
} /* Table_Interp1d */

/******************************************************************************
* double Table_Interp1d_nearest(x, x1, y1, x2, y2)
*    ACTION: table lookup with nearest method at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d_nearest(double x,
  double x1, double y1,
  double x2, double y2)
{
  if (fabs(x-x1) < fabs(x-x2)) return (y1);
  else return(y2);
} /* Table_Interp1d_nearest */

/******************************************************************************
* double Table_Interp2d(x,y, x1,y1, x2,y2, z11,z12,z21,z22)
*    ACTION: interpolates bi-linearly at (x,y) between z1=f(x1,y1) and z2=f(x2,y2)
*    return: z=f(x,y) value
*    x,y |   x1   x2
*    ----------------
*     y1 |   z11  z21
*     y2 |   z12  z22
*******************************************************************************/
double Table_Interp2d(double x, double y,
  double x1, double y1,
  double x2, double y2,
  double z11, double z12, double z21, double z22)
{
  double ratio_x, ratio_y;
  if (x2 == x1) return Table_Interp1d(y, y1,z11, y2,z12);
  if (y1 == y2) return Table_Interp1d(x, x1,z11, x2,z21);

  ratio_y = (y - y1)/(y2 - y1);
  ratio_x = (x - x1)/(x2 - x1);
  return (1-ratio_x)*(1-ratio_y)*z11 + ratio_x*(1-ratio_y)*z21
    + ratio_x*ratio_y*z22         + (1-ratio_x)*ratio_y*z12;
} /* Table_Interp2d */

/* end of read_table-lib.c */


long Virtual_input_Read_Input(char *aFile, char *aType, t_Table *aTable, long *aOffset)
  {
    long max_lines = 50000;
    long length=0;
    char bType[32];

    if (!aFile) return (0);
    if (aType) strcpy(bType, aType);
    else strcpy(bType, "???");

    Table_Free(aTable);

    /* Try to Open neutron input text filename. */
    if((aFile && aType == NULL) || !strcmp(bType,"text")) {
      Table_Read_Offset(aTable, aFile, 0, aOffset, max_lines);  /* read data from filename into rTable */
      strcpy(bType, "text");
    }
    if (!aTable->data && aType && aType[0] != 't')
      Table_Read_Offset_Binary(aTable, aFile, aType, aOffset, max_lines, 11);

    return(aTable->rows);
  }
#line 8026 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

/* Shared user declarations for all components 'Monitor_nD'. */
#line 224 "/usr/local/lib/mcstas-2.0/monitors/Monitor_nD.comp"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/monitor_nd-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision: 1.19 $
*
* This file is to be imported by the monitor_nd related components
* It handles some shared functions.
*
* Usage: within SHARE
* %include "monitor_nd-lib"
*
*******************************************************************************/

#ifndef MONITOR_ND_LIB_H

#define MONITOR_ND_LIB_H "$Revision: 1.19 $"
#define MONnD_COORD_NMAX  30  /* max number of variables to record */

  typedef struct MonitornD_Defines
  {
    int COORD_NONE  ;
    int COORD_X     ;
    int COORD_Y     ;
    int COORD_Z     ;
    int COORD_RADIUS; 
    int COORD_VX    ;
    int COORD_VY    ;
    int COORD_VZ    ;
    int COORD_V     ;
    int COORD_T     ;
    int COORD_P     ;
    int COORD_SX    ;
    int COORD_SY    ;
    int COORD_SZ    ;
    int COORD_KX    ;
    int COORD_KY    ;
    int COORD_KZ    ;
    int COORD_K     ;
    int COORD_ENERGY;
    int COORD_LAMBDA;
    int COORD_KXY   ;
    int COORD_KYZ   ;
    int COORD_KXZ   ;
    int COORD_VXY   ;
    int COORD_VYZ   ;
    int COORD_VXZ   ;
    int COORD_HDIV  ;
    int COORD_VDIV  ;
    int COORD_ANGLE ;
    int COORD_NCOUNT;
    int COORD_THETA ;
    int COORD_PHI   ;
    int COORD_USER1 ;
    int COORD_USER2 ;
    int COORD_USER3 ;
    int COORD_XY    ;
    int COORD_XZ    ;
    int COORD_YZ    ;

    /* token modifiers */
    int COORD_VAR   ; /* next token should be a variable or normal option */
    int COORD_MIN   ; /* next token is a min value */
    int COORD_MAX   ; /* next token is a max value */
    int COORD_DIM   ; /* next token is a bin value */
    int COORD_FIL   ; /* next token is a filename */
    int COORD_EVNT  ; /* next token is a buffer size value */
    int COORD_3HE   ; /* next token is a 3He pressure value */
    int COORD_LOG   ; /* next variable will be in log scale */
    int COORD_ABS   ; /* next variable will be in abs scale */
    int COORD_SIGNAL; /* next variable will be the signal var */
    int COORD_AUTO  ; /* set auto limits */

    char TOKEN_DEL[32]; /* token separators */

    char SHAPE_SQUARE; /* shape of the monitor */
    char SHAPE_DISK  ;
    char SHAPE_SPHERE;
    char SHAPE_CYLIND;
    char SHAPE_BANANA; /* cylinder without top/bottom, on restricted angular area */
    char SHAPE_BOX   ;
    char SHAPE_PREVIOUS;

  } MonitornD_Defines_type;

  typedef struct MonitornD_Variables
  {
    double area;
    double Sphere_Radius     ;
    double Cylinder_Height   ;
    char   Flag_With_Borders ;   /* 2 means xy borders too */
    char   Flag_List         ;   /* 1 store 1 buffer, 2 is list all, 3 list all+append */
    char   Flag_Multiple     ;   /* 1 when n1D, 0 for 2D */
    char   Flag_Verbose      ;
    int    Flag_Shape        ;
    char   Flag_Auto_Limits  ;   /* get limits from first Buffer */
    char   Flag_Absorb       ;   /* monitor is also a slit */
    char   Flag_Exclusive    ;   /* absorb neutrons out of monitor limits */
    char   Flag_per_cm2      ;   /* flux is per cm2 */
    char   Flag_log          ;   /* log10 of the flux */
    char   Flag_parallel     ;   /* set neutron state back after detection (parallel components) */
    char   Flag_Binary_List  ;
    char   Flag_capture      ;   /* lambda monitor with lambda/lambda(2200m/s = 1.7985 Angs) weightening */
    int    Flag_signal       ;   /* 0:monitor p, else monitor a mean value */

    long   Coord_Number      ;   /* total number of variables to monitor, plus intensity (0) */
    long   Buffer_Block      ;   /* Buffer size for list or auto limits */
    long   Neutron_Counter   ;   /* event counter, simulation total counts is mcget_ncount() */
    long   Buffer_Counter    ;   /* index in Buffer size (for realloc) */
    long   Buffer_Size       ;
    int    Coord_Type[MONnD_COORD_NMAX];    /* type of variable */
    char   Coord_Label[MONnD_COORD_NMAX][30];       /* label of variable */
    char   Coord_Var[MONnD_COORD_NMAX][30]; /* short id of variable */
    long   Coord_Bin[MONnD_COORD_NMAX];             /* bins of variable array */
    double Coord_Min[MONnD_COORD_NMAX];
    double Coord_Max[MONnD_COORD_NMAX];
    char   Monitor_Label[MONnD_COORD_NMAX*30];      /* Label for monitor */
    char   Mon_File[128];    /* output file name */

    double cx,cy,cz;
    double cvx, cvy, cvz;
    double ckx, cky, ckz;
    double csx, csy, csz;
    double cEx, cEy, cEz;
    double cs1, cs2, ct, cphi, cp;
    double He3_pressure;
    char   Flag_UsePreMonitor    ;   /* use a previously stored neutron parameter set */
    char   UserName1[128];
    char   UserName2[128];
    char   UserName3[128];
    double UserVariable1;
    double UserVariable2;
    double UserVariable3;
    char   option[CHAR_BUF_LENGTH];

    double Nsum;
    double psum, p2sum;
    double **Mon2D_N;
    double **Mon2D_p;
    double **Mon2D_p2;
    double *Mon2D_Buffer;

    double mxmin,mxmax,mymin,mymax,mzmin,mzmax;
    double mean_dx, mean_dy, min_x, min_y, max_x, max_y, mean_p;

    char   compcurname[128];
    Coords compcurpos;

  } MonitornD_Variables_type;

/* monitor_nd-lib function prototypes */
/* ========================================================================= */

void Monitor_nD_Init(MonitornD_Defines_type *, MonitornD_Variables_type *, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM);
double Monitor_nD_Trace(MonitornD_Defines_type *, MonitornD_Variables_type *);
MCDETECTOR Monitor_nD_Save(MonitornD_Defines_type *, MonitornD_Variables_type *);
void Monitor_nD_Finally(MonitornD_Defines_type *, MonitornD_Variables_type *);
void Monitor_nD_McDisplay(MonitornD_Defines_type *,
 MonitornD_Variables_type *);

#define MONND_DECLARE(monname) \
  struct MonitornD_Variables *mcmonnd ## monname;
#define MONND_USER_TITLE(monname, num, title) \
  { mcmonnd ## monname = &(MC_GETPAR(monname, Vars)); \
    strcpy(mcmonnd ## monname->UserName ## num, title); }
#define MONND_USER_VALUE(monname, num, value) \
  { mcmonnd ## monname = &(MC_GETPAR(monname, Vars)); \
    mcmonnd ## monname->UserVariable ## num = (value); }

#endif

/* end of monitor_nd-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/monitor_nd-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision: 1.46 $
*
* This file is to be imported by the monitor_nd related components
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "monitor_nd-lib"
*
*******************************************************************************/

#ifndef MONITOR_ND_LIB_H
#error McStas : please import this library with %include "monitor_nd-lib"
#endif

/* ========================================================================= */
/* Monitor_nD_Init: this routine is used to parse options                    */
/* ========================================================================= */

void Monitor_nD_Init(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars,
  MCNUM xwidth,
  MCNUM yheight,
  MCNUM zdepth,
  MCNUM xmin,
  MCNUM xmax,
  MCNUM ymin,
  MCNUM ymax,
  MCNUM zmin,
  MCNUM zmax)
  {
    long carg = 1;
    char *option_copy, *token;
    char Flag_New_token = 1;
    char Flag_End       = 1;
    char Flag_All       = 0;
    char Flag_No        = 0;
    char Flag_abs       = 0;
    int  Flag_auto      = 0;  /* -1: all, 1: the current variable */
    int  Set_Vars_Coord_Type;
    char Set_Vars_Coord_Label[64];
    char Set_Vars_Coord_Var[64];
    char Short_Label[MONnD_COORD_NMAX][64];
    int  Set_Coord_Mode;
    long i=0, j=0;
    double lmin, lmax, XY=0;
    long t;


    t = (long)time(NULL);

/* initialize DEFS */
/* Variables to monitor */
    DEFS->COORD_NONE   =0;
    DEFS->COORD_X      =1;
    DEFS->COORD_Y      =2;
    DEFS->COORD_Z      =3;
    DEFS->COORD_RADIUS =19;
    DEFS->COORD_VX     =4;
    DEFS->COORD_VY     =5;
    DEFS->COORD_VZ     =6;
    DEFS->COORD_V      =16;
    DEFS->COORD_T      =7;
    DEFS->COORD_P      =8;
    DEFS->COORD_SX     =9;
    DEFS->COORD_SY     =10;
    DEFS->COORD_SZ     =11;
    DEFS->COORD_KX     =12;
    DEFS->COORD_KY     =13;
    DEFS->COORD_KZ     =14;
    DEFS->COORD_K      =15;
    DEFS->COORD_ENERGY =17;
    DEFS->COORD_LAMBDA =18;
    DEFS->COORD_HDIV   =20;
    DEFS->COORD_VDIV   =21;
    DEFS->COORD_ANGLE  =22;
    DEFS->COORD_NCOUNT =23;
    DEFS->COORD_THETA  =24;
    DEFS->COORD_PHI    =25;
    DEFS->COORD_USER1  =26;
    DEFS->COORD_USER2  =27;
    DEFS->COORD_USER3  =28;
    DEFS->COORD_XY     =37;
    DEFS->COORD_YZ     =31;
    DEFS->COORD_XZ     =32;
    DEFS->COORD_VXY    =30;
    DEFS->COORD_VYZ    =34;
    DEFS->COORD_VXZ    =36;
    DEFS->COORD_KXY    =29;
    DEFS->COORD_KYZ    =33;
    DEFS->COORD_KXZ    =35;



/* token modifiers */
    DEFS->COORD_VAR    =0;    /* next token should be a variable or normal option */
    DEFS->COORD_MIN    =1;    /* next token is a min value */
    DEFS->COORD_MAX    =2;    /* next token is a max value */
    DEFS->COORD_DIM    =3;    /* next token is a bin value */
    DEFS->COORD_FIL    =4;    /* next token is a filename */
    DEFS->COORD_EVNT   =5;    /* next token is a buffer size value */
    DEFS->COORD_3HE    =6;    /* next token is a 3He pressure value */
    DEFS->COORD_LOG    =64;   /* next variable will be in log scale */
    DEFS->COORD_ABS    =128;  /* next variable will be in abs scale */
    DEFS->COORD_SIGNAL =256;  /* next variable will be the signal var */
    DEFS->COORD_AUTO   =512;  /* set auto limits */

    strcpy(DEFS->TOKEN_DEL, " =,;[](){}:");  /* token separators */

    DEFS->SHAPE_SQUARE =0;    /* shape of the monitor */
    DEFS->SHAPE_DISK   =1;
    DEFS->SHAPE_SPHERE =2;
    DEFS->SHAPE_CYLIND =3;
    DEFS->SHAPE_BANANA =4;
    DEFS->SHAPE_BOX    =5;
    DEFS->SHAPE_PREVIOUS=6;

    Vars->Sphere_Radius     = 0;
    Vars->Cylinder_Height   = 0;
    Vars->Flag_With_Borders = 0;   /* 2 means xy borders too */
    Vars->Flag_List         = 0;   /* 1=store 1 buffer, 2=list all, 3=re-use buffer */
    Vars->Flag_Multiple     = 0;   /* 1 when n1D, 0 for 2D */
    Vars->Flag_Verbose      = 0;
    Vars->Flag_Shape        = DEFS->SHAPE_SQUARE;
    Vars->Flag_Auto_Limits  = 0;   /* get limits from first Buffer */
    Vars->Flag_Absorb       = 0;   /* monitor is also a slit */
    Vars->Flag_Exclusive    = 0;   /* absorb neutrons out of monitor limits */
    Vars->Flag_per_cm2      = 0;   /* flux is per cm2 */
    Vars->Flag_log          = 0;   /* log10 of the flux */
    Vars->Flag_parallel     = 0;   /* set neutron state back after detection (parallel components) */
    Vars->Flag_Binary_List  = 0;   /* save list as a binary file (smaller) */
    Vars->Coord_Number      = 0;   /* total number of variables to monitor, plus intensity (0) */
    Vars->Buffer_Block      = 10000;     /* Buffer size for list or auto limits */
    Vars->Neutron_Counter   = 0;   /* event counter, simulation total counts is mcget_ncount() */
    Vars->Buffer_Counter    = 0;   /* index in Buffer size (for realloc) */
    Vars->Buffer_Size       = 0;
    Vars->UserVariable1     = 0;
    Vars->UserVariable2     = 0;
    Vars->He3_pressure      = 0;
    Vars->Flag_capture      = 0;
    Vars->Flag_signal       = DEFS->COORD_P;
    Vars->mean_dx=Vars->mean_dy=0;
    Vars->min_x = Vars->max_x  =0;
    Vars->min_y = Vars->max_y  =0;

    Set_Vars_Coord_Type = DEFS->COORD_NONE;
    Set_Coord_Mode = DEFS->COORD_VAR;

    /* handle size parameters */
    /* normal use is with xwidth, yheight, zdepth */
    /* if xmin,xmax,ymin,ymax,zmin,zmax are non 0, use them */
    if (fabs(xmin-xmax) == 0)
      { Vars->mxmin = -fabs(xwidth)/2; Vars->mxmax = fabs(xwidth)/2; }
    else
      { if (xmin < xmax) {Vars->mxmin = xmin; Vars->mxmax = xmax;}
        else {Vars->mxmin = xmax; Vars->mxmax = xmin;}
      }
    if (fabs(ymin-ymax) == 0)
      { Vars->mymin = -fabs(yheight)/2; Vars->mymax = fabs(yheight)/2; }
    else
      { if (ymin < ymax) {Vars->mymin = ymin; Vars->mymax = ymax;}
        else {Vars->mymin = ymax; Vars->mymax = ymin;}
      }
    if (fabs(zmin-zmax) == 0)
      { Vars->mzmin = -fabs(zdepth)/2; Vars->mzmax = fabs(zdepth)/2; }
    else
      { if (zmin < zmax) {Vars->mzmin = zmin; Vars->mzmax = zmax; }
        else {Vars->mzmin = zmax; Vars->mzmax = zmin; }
      }

    if (fabs(Vars->mzmax-Vars->mzmin) == 0)
      Vars->Flag_Shape        = DEFS->SHAPE_SQUARE;
    else
      Vars->Flag_Shape        = DEFS->SHAPE_BOX;

    /* parse option string */

    option_copy = (char*)malloc(strlen(Vars->option)+1);
    if (option_copy == NULL)
    {
      fprintf(stderr,"Monitor_nD: %s cannot allocate 'options' copy (%li). Fatal.\n", Vars->compcurname, (long)strlen(Vars->option));
      exit(-1);
    }

    if (strlen(Vars->option))
    {
      Flag_End = 0;
      strcpy(option_copy, Vars->option);
    }

    if (strstr(Vars->option, "cm2") || strstr(Vars->option, "cm^2")) Vars->Flag_per_cm2 = 1;

    if (strstr(Vars->option, "binary") || strstr(Vars->option, "float"))
      Vars->Flag_Binary_List  = 1;
    if (strstr(Vars->option, "double"))
      Vars->Flag_Binary_List  = 2;

    strcpy(Vars->Coord_Label[0],"Intensity");
    strncpy(Vars->Coord_Var[0],"p",30);
    Vars->Coord_Type[0] = DEFS->COORD_P;
    Vars->Coord_Bin[0] = 1;
    Vars->Coord_Min[0] = 0;
    Vars->Coord_Max[0] = FLT_MAX;

    /* default file name is comp_name+dateID */
    sprintf(Vars->Mon_File, "%s_%li", Vars->compcurname, t);

    carg = 1;
    while((Flag_End == 0) && (carg < 128))
    {

      if (Flag_New_token) /* retain previous token or get a new one */
      {
        if (carg == 1) token=(char *)strtok(option_copy,DEFS->TOKEN_DEL);
        else token=(char *)strtok(NULL,DEFS->TOKEN_DEL);
        if (token == NULL) Flag_End=1;
      }
      Flag_New_token = 1;
      if ((token != NULL) && (strlen(token) != 0))
      {
        char iskeyword=0;
        int  old_Mode;
        /* change token to lower case */
        for (i=0; i<strlen(token); i++) token[i]=tolower(token[i]);
        /* first handle option values from preceeding keyword token detected */
        old_Mode = Set_Coord_Mode;
        if (Set_Coord_Mode == DEFS->COORD_MAX)  /* max=%i */
        {
          if (!Flag_All)
            Vars->Coord_Max[Vars->Coord_Number] = atof(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Max[i++] = atof(token));
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_MIN)  /* min=%i */
        {
          if (!Flag_All)
            Vars->Coord_Min[Vars->Coord_Number] = atof(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Min[i++] = atof(token));
          Set_Coord_Mode = DEFS->COORD_MAX;
        }
        if (Set_Coord_Mode == DEFS->COORD_DIM)  /* bins=%i */
        {
          if (!Flag_All)
            Vars->Coord_Bin[Vars->Coord_Number] = atoi(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Bin[i++] = atoi(token));
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_FIL)  /* file=%s */
        {
          if (!Flag_No) strncpy(Vars->Mon_File,token,128);
          else { strcpy(Vars->Mon_File,""); Vars->Coord_Number = 0; Flag_End = 1;}
          Set_Coord_Mode = DEFS->COORD_VAR;
        }
        if (Set_Coord_Mode == DEFS->COORD_EVNT) /* list=%i */
        {
          if (!strcmp(token, "all") || Flag_All) Vars->Flag_List = 2;
          else { i = (long)ceil(atof(token)); if (i) Vars->Buffer_Block = i;
            Vars->Flag_List = 1; }
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_3HE)  /* pressure=%g */
        {
            Vars->He3_pressure = atof(token);
            Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }

        /* now look for general option keywords */
        if (!strcmp(token, "borders"))  {Vars->Flag_With_Borders = 1; iskeyword=1; }
        if (!strcmp(token, "verbose"))  {Vars->Flag_Verbose      = 1; iskeyword=1; }
        if (!strcmp(token, "log"))      {Vars->Flag_log          = 1; iskeyword=1; }
        if (!strcmp(token, "abs"))      {Flag_abs                = 1; iskeyword=1; }
        if (!strcmp(token, "multiple")) {Vars->Flag_Multiple     = 1; iskeyword=1; }
        if (!strcmp(token, "exclusive")){Vars->Flag_Exclusive    = 1; iskeyword=1; }
        if (!strcmp(token, "list") || !strcmp(token, "events")) {
          Vars->Flag_List = 1; Set_Coord_Mode = DEFS->COORD_EVNT;  }
        if (!strcmp(token, "limits") || !strcmp(token, "min"))
          Set_Coord_Mode = DEFS->COORD_MIN;
        if (!strcmp(token, "slit") || !strcmp(token, "absorb")) {
          Vars->Flag_Absorb = 1;  iskeyword=1; }
        if (!strcmp(token, "max"))  Set_Coord_Mode = DEFS->COORD_MAX;
        if (!strcmp(token, "bins") || !strcmp(token, "dim")) Set_Coord_Mode = DEFS->COORD_DIM;
        if (!strcmp(token, "file") || !strcmp(token, "filename")) {
          Set_Coord_Mode = DEFS->COORD_FIL;
          if (Flag_No) { strcpy(Vars->Mon_File,""); Vars->Coord_Number = 0; Flag_End = 1; }
        }
        if (!strcmp(token, "unactivate")) {
          Flag_End = 1; Vars->Coord_Number = 0; iskeyword=1; }
        if (!strcmp(token, "all"))    { Flag_All = 1;  iskeyword=1; }
        if (!strcmp(token, "sphere")) { Vars->Flag_Shape = DEFS->SHAPE_SPHERE; iskeyword=1; }
        if (!strcmp(token, "cylinder")) { Vars->Flag_Shape = DEFS->SHAPE_CYLIND; iskeyword=1; }
        if (!strcmp(token, "banana")) { Vars->Flag_Shape = DEFS->SHAPE_BANANA; iskeyword=1; }
        if (!strcmp(token, "square")) { Vars->Flag_Shape = DEFS->SHAPE_SQUARE; iskeyword=1; }
        if (!strcmp(token, "disk"))   { Vars->Flag_Shape = DEFS->SHAPE_DISK; iskeyword=1; }
        if (!strcmp(token, "box"))     { Vars->Flag_Shape = DEFS->SHAPE_BOX; iskeyword=1; }
        if (!strcmp(token, "previous")) { Vars->Flag_Shape = DEFS->SHAPE_PREVIOUS; iskeyword=1; }
        if (!strcmp(token, "parallel")){ Vars->Flag_parallel = 1; iskeyword=1; }
        if (!strcmp(token, "capture")) { Vars->Flag_capture = 1; iskeyword=1; }
        if (!strcmp(token, "auto") && (Flag_auto != -1)) {
          Vars->Flag_Auto_Limits = 1;
          if (Flag_All) Flag_auto = -1;
          else Flag_auto = 1;
          iskeyword=1; }
        if (!strcmp(token, "premonitor")) {
          Vars->Flag_UsePreMonitor = 1; iskeyword=1; }
        if (!strcmp(token, "3He_pressure") || !strcmp(token, "pressure")) {
          Vars->He3_pressure = 3; iskeyword=1; }
        if (!strcmp(token, "no") || !strcmp(token, "not")) { Flag_No = 1;  iskeyword=1; }
        if (!strcmp(token, "signal")) Set_Coord_Mode = DEFS->COORD_SIGNAL;

        /* Mode has changed: this was a keyword or value  ? */
        if (Set_Coord_Mode != old_Mode) iskeyword=1;

        /* now look for variable names to monitor */
        Set_Vars_Coord_Type = DEFS->COORD_NONE; lmin = 0; lmax = 0;

        if (!strcmp(token, "x"))
          { Set_Vars_Coord_Type = DEFS->COORD_X; strcpy(Set_Vars_Coord_Label,"x [m]"); strcpy(Set_Vars_Coord_Var,"x");
          lmin = Vars->mxmin; lmax = Vars->mxmax;
          Vars->Coord_Min[Vars->Coord_Number+1] = Vars->mxmin;
          Vars->Coord_Max[Vars->Coord_Number+1] = Vars->mxmax;}
        if (!strcmp(token, "y"))
          { Set_Vars_Coord_Type = DEFS->COORD_Y; strcpy(Set_Vars_Coord_Label,"y [m]"); strcpy(Set_Vars_Coord_Var,"y");
          lmin = Vars->mymin; lmax = Vars->mymax;
          Vars->Coord_Min[Vars->Coord_Number+1] = Vars->mymin;
          Vars->Coord_Max[Vars->Coord_Number+1] = Vars->mymax;}
        if (!strcmp(token, "z"))
          { Set_Vars_Coord_Type = DEFS->COORD_Z; strcpy(Set_Vars_Coord_Label,"z [m]"); strcpy(Set_Vars_Coord_Var,"z"); lmin = Vars->mzmin; lmax = Vars->mzmax; }
        if (!strcmp(token, "k") || !strcmp(token, "wavevector"))
          { Set_Vars_Coord_Type = DEFS->COORD_K; strcpy(Set_Vars_Coord_Label,"|k| [Angs-1]"); strcpy(Set_Vars_Coord_Var,"k"); lmin = 0; lmax = 10; }
        if (!strcmp(token, "v"))
          { Set_Vars_Coord_Type = DEFS->COORD_V; strcpy(Set_Vars_Coord_Label,"Velocity [m/s]"); strcpy(Set_Vars_Coord_Var,"v"); lmin = 0; lmax = 10000; }
        if (!strcmp(token, "t") || !strcmp(token, "time"))
          { Set_Vars_Coord_Type = DEFS->COORD_T; strcpy(Set_Vars_Coord_Label,"TOF [s]"); strcpy(Set_Vars_Coord_Var,"t"); lmin = 0; lmax = .1; }
        if ((!strcmp(token, "p") || !strcmp(token, "i") || !strcmp(token, "intensity") || !strcmp(token, "flux")))
          { Set_Vars_Coord_Type = DEFS->COORD_P;
            strcpy(Set_Vars_Coord_Label,"Intensity");
            strncat(Set_Vars_Coord_Label, " [n/s", 30);
            if (Vars->Flag_per_cm2) strncat(Set_Vars_Coord_Label, "/cm2", 30);
            if (XY > 1 && Vars->Coord_Number)
              strncat(Set_Vars_Coord_Label, "/bin", 30);
            strncat(Set_Vars_Coord_Label, "]", 30);
            strcpy(Set_Vars_Coord_Var,"I");
            lmin = 0; lmax = FLT_MAX;
          }

        if (!strcmp(token, "vx"))
          { Set_Vars_Coord_Type = DEFS->COORD_VX; strcpy(Set_Vars_Coord_Label,"vx [m/s]"); strcpy(Set_Vars_Coord_Var,"vx"); lmin = -1000; lmax = 1000; }
        if (!strcmp(token, "vy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VY; strcpy(Set_Vars_Coord_Label,"vy [m/s]"); strcpy(Set_Vars_Coord_Var,"vy"); lmin = -1000; lmax = 1000; }
        if (!strcmp(token, "vz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VZ; strcpy(Set_Vars_Coord_Label,"vz [m/s]"); strcpy(Set_Vars_Coord_Var,"vz"); lmin = -10000; lmax = 10000; }
        if (!strcmp(token, "kx"))
          { Set_Vars_Coord_Type = DEFS->COORD_KX; strcpy(Set_Vars_Coord_Label,"kx [Angs-1]"); strcpy(Set_Vars_Coord_Var,"kx"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "ky"))
          { Set_Vars_Coord_Type = DEFS->COORD_KY; strcpy(Set_Vars_Coord_Label,"ky [Angs-1]"); strcpy(Set_Vars_Coord_Var,"ky"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "kz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KZ; strcpy(Set_Vars_Coord_Label,"kz [Angs-1]"); strcpy(Set_Vars_Coord_Var,"kz"); lmin = -10; lmax = 10; }
        if (!strcmp(token, "sx"))
          { Set_Vars_Coord_Type = DEFS->COORD_SX; strcpy(Set_Vars_Coord_Label,"sx [1]"); strcpy(Set_Vars_Coord_Var,"sx"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "sy"))
          { Set_Vars_Coord_Type = DEFS->COORD_SY; strcpy(Set_Vars_Coord_Label,"sy [1]"); strcpy(Set_Vars_Coord_Var,"sy"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "sz"))
          { Set_Vars_Coord_Type = DEFS->COORD_SZ; strcpy(Set_Vars_Coord_Label,"sz [1]"); strcpy(Set_Vars_Coord_Var,"sz"); lmin = -1; lmax = 1; }

        if (!strcmp(token, "energy") || !strcmp(token, "omega") || !strcmp(token, "e"))
          { Set_Vars_Coord_Type = DEFS->COORD_ENERGY; strcpy(Set_Vars_Coord_Label,"Energy [meV]"); strcpy(Set_Vars_Coord_Var,"E"); lmin = 0; lmax = 100; }
        if (!strcmp(token, "lambda") || !strcmp(token, "wavelength") || !strcmp(token, "l"))
          { Set_Vars_Coord_Type = DEFS->COORD_LAMBDA; strcpy(Set_Vars_Coord_Label,"Wavelength [Angs]"); strcpy(Set_Vars_Coord_Var,"L"); lmin = 0; lmax = 100; }
        if (!strcmp(token, "radius") || !strcmp(token, "r"))
          { Set_Vars_Coord_Type = DEFS->COORD_RADIUS; strcpy(Set_Vars_Coord_Label,"Radius [m]"); strcpy(Set_Vars_Coord_Var,"xy"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "xy"))
          { Set_Vars_Coord_Type = DEFS->COORD_XY; strcpy(Set_Vars_Coord_Label,"Radius (xy) [m]"); strcpy(Set_Vars_Coord_Var,"xy"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "yz"))
          { Set_Vars_Coord_Type = DEFS->COORD_YZ; strcpy(Set_Vars_Coord_Label,"Radius (yz) [m]"); strcpy(Set_Vars_Coord_Var,"yz"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "xz"))
          { Set_Vars_Coord_Type = DEFS->COORD_XZ; strcpy(Set_Vars_Coord_Label,"Radius (xz) [m]"); strcpy(Set_Vars_Coord_Var,"xz"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "vxy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VXY; strcpy(Set_Vars_Coord_Label,"Radial Velocity (xy) [m]"); strcpy(Set_Vars_Coord_Var,"Vxy"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kxy"))
          { Set_Vars_Coord_Type = DEFS->COORD_KXY; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (xy) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kxy"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "vyz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VYZ; strcpy(Set_Vars_Coord_Label,"Radial Velocity (yz) [m]"); strcpy(Set_Vars_Coord_Var,"Vyz"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kyz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KYZ; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (yz) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kyz"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "vxz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VXZ; strcpy(Set_Vars_Coord_Label,"Radial Velocity (xz) [m]"); strcpy(Set_Vars_Coord_Var,"Vxz"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kxz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KXZ; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (xz) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kxz"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "angle") || !strcmp(token, "a"))
          { Set_Vars_Coord_Type = DEFS->COORD_ANGLE; strcpy(Set_Vars_Coord_Label,"Angle [deg]"); strcpy(Set_Vars_Coord_Var,"A"); lmin = -50; lmax = 50; }
        if (!strcmp(token, "hdiv")|| !strcmp(token, "divergence") || !strcmp(token, "xdiv") || !strcmp(token, "hd") || !strcmp(token, "dx"))
          { Set_Vars_Coord_Type = DEFS->COORD_HDIV; strcpy(Set_Vars_Coord_Label,"Hor. Divergence [deg]"); strcpy(Set_Vars_Coord_Var,"hd"); lmin = -5; lmax = 5; }
        if (!strcmp(token, "vdiv") || !strcmp(token, "ydiv") || !strcmp(token, "vd") || !strcmp(token, "dy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VDIV; strcpy(Set_Vars_Coord_Label,"Vert. Divergence [deg]"); strcpy(Set_Vars_Coord_Var,"vd"); lmin = -5; lmax = 5; }
        if (!strcmp(token, "theta") || !strcmp(token, "longitude") || !strcmp(token, "th"))
          { Set_Vars_Coord_Type = DEFS->COORD_THETA; strcpy(Set_Vars_Coord_Label,"Longitude [deg]"); strcpy(Set_Vars_Coord_Var,"th"); lmin = -180; lmax = 180; }
        if (!strcmp(token, "phi") || !strcmp(token, "lattitude") || !strcmp(token, "ph"))
          { Set_Vars_Coord_Type = DEFS->COORD_PHI; strcpy(Set_Vars_Coord_Label,"Lattitude [deg]"); strcpy(Set_Vars_Coord_Var,"ph"); lmin = -180; lmax = 180; }
        if (!strcmp(token, "ncounts") || !strcmp(token, "n"))
          { Set_Vars_Coord_Type = DEFS->COORD_NCOUNT; strcpy(Set_Vars_Coord_Label,"Neutrons [1]"); strcpy(Set_Vars_Coord_Var,"n"); lmin = 0; lmax = 1e10; }
        if (!strcmp(token, "user") || !strcmp(token, "user1") || !strcmp(token, "u1"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER1; strncpy(Set_Vars_Coord_Label,Vars->UserName1,30); strcpy(Set_Vars_Coord_Var,"U1"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "user2") || !strcmp(token, "u2"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER2; strncpy(Set_Vars_Coord_Label,Vars->UserName2,30); strcpy(Set_Vars_Coord_Var,"U2"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "user3") || !strcmp(token, "u3"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER3; strncpy(Set_Vars_Coord_Label,Vars->UserName3,30); strcpy(Set_Vars_Coord_Var,"U3"); lmin = -1e10; lmax = 1e10; }

        /* now stores variable keywords detected, if any */
        if (Set_Vars_Coord_Type != DEFS->COORD_NONE)
        {
          int Coord_Number = Vars->Coord_Number;
          if (Vars->Flag_log) { Set_Vars_Coord_Type |= DEFS->COORD_LOG; Vars->Flag_log = 0; }
          if (Flag_abs) { Set_Vars_Coord_Type |= DEFS->COORD_ABS; Flag_abs = 0; }
          if (Flag_auto != 0) { Set_Vars_Coord_Type |= DEFS->COORD_AUTO; Flag_auto = 0; }
          if (Set_Coord_Mode == DEFS->COORD_SIGNAL)
          {
            Coord_Number = 0;
            Vars->Flag_signal = Set_Vars_Coord_Type;
          }
          else
          {
            if (Coord_Number < MONnD_COORD_NMAX)
            { Coord_Number++;
              Vars->Coord_Number = Coord_Number; }
            else if (Vars->Flag_Verbose) printf("Monitor_nD: %s reached max number of variables (%i).\n", Vars->compcurname, MONnD_COORD_NMAX);
          }
          Vars->Coord_Type[Coord_Number] = Set_Vars_Coord_Type;
          strncpy(Vars->Coord_Label[Coord_Number], Set_Vars_Coord_Label,30);
          strncpy(Vars->Coord_Var[Coord_Number], Set_Vars_Coord_Var,30);
          if (lmin > lmax) { XY = lmin; lmin=lmax; lmax = XY; }
          Vars->Coord_Min[Coord_Number] = lmin;
          Vars->Coord_Max[Coord_Number] = lmax;
          if (Set_Coord_Mode != DEFS->COORD_SIGNAL) Vars->Coord_Bin[Coord_Number] = 20;
          Set_Coord_Mode = DEFS->COORD_VAR;
          Flag_All = 0;
          Flag_No  = 0;
        } else {
          /* no variable name could be read from options */
          if (!iskeyword) {
            if (strcmp(token, "cm2") && strcmp(token, "incoming")
             && strcmp(token, "outgoing") && strcmp(token, "cm2")
             && strcmp(token, "cm^2") && strcmp(token, "float")
             && strcmp(token, "double") && strcmp(token, "binary")
             && strcmp(token, "steradian") && Vars->Flag_Verbose)
              printf("Monitor_nD: %s: unknown '%s' keyword in 'options'. Ignoring.\n", Vars->compcurname, token);
          }
        }
      carg++;
      } /* end if token */
    } /* end while carg */
    free(option_copy);
    if (carg == 128) printf("Monitor_nD: %s reached max number of tokens (%i). Skipping.\n", Vars->compcurname, 128);

    if ((Vars->Flag_Shape == DEFS->SHAPE_BOX) && (fabs(Vars->mzmax - Vars->mzmin) == 0)) Vars->Flag_Shape = DEFS->SHAPE_SQUARE;

    if (Vars->Flag_log == 1) Vars->Coord_Type[0] |= DEFS->COORD_LOG;
    if (Vars->Coord_Number == 0)
    { Vars->Flag_Auto_Limits=0; Vars->Flag_Multiple=0; Vars->Flag_List=0; }

    /* now setting Monitor Name from variable labels */
    strcpy(Vars->Monitor_Label,"");
    XY = 1; /* will contain total bin number */
    for (i = 0; i <= Vars->Coord_Number; i++)
    {
      if (Flag_auto != 0) Vars->Coord_Type[i] |= DEFS->COORD_AUTO;
      Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
      if ((Set_Vars_Coord_Type == DEFS->COORD_X)
       || (Set_Vars_Coord_Type == DEFS->COORD_Y)
       || (Set_Vars_Coord_Type == DEFS->COORD_Z))
       strcpy(Short_Label[i],"Position");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_THETA)
       || (Set_Vars_Coord_Type == DEFS->COORD_PHI)
       || (Set_Vars_Coord_Type == DEFS->COORD_ANGLE))
       strcpy(Short_Label[i],"Angle");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_XY)
       || (Set_Vars_Coord_Type == DEFS->COORD_XZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_YZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_RADIUS))
       strcpy(Short_Label[i],"Radius");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_VX)
       || (Set_Vars_Coord_Type == DEFS->COORD_VY)
       || (Set_Vars_Coord_Type == DEFS->COORD_VZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_V)
       || (Set_Vars_Coord_Type == DEFS->COORD_VXY)
       || (Set_Vars_Coord_Type == DEFS->COORD_VYZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_VXZ))
       strcpy(Short_Label[i],"Velocity");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_KX)
       || (Set_Vars_Coord_Type == DEFS->COORD_KY)
       || (Set_Vars_Coord_Type == DEFS->COORD_KZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_KXY)
       || (Set_Vars_Coord_Type == DEFS->COORD_KYZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_KXZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_K))
       strcpy(Short_Label[i],"Wavevector");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_SX)
       || (Set_Vars_Coord_Type == DEFS->COORD_SY)
       || (Set_Vars_Coord_Type == DEFS->COORD_SZ))
       strcpy(Short_Label[i],"Spin");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_HDIV)
       || (Set_Vars_Coord_Type == DEFS->COORD_VDIV))
       strcpy(Short_Label[i],"Divergence");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_ENERGY)
       strcpy(Short_Label[i],"Energy");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_LAMBDA)
       strcpy(Short_Label[i],"Wavelength");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_NCOUNT)
       strcpy(Short_Label[i],"Neutron counts");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_T)
          strcpy(Short_Label[i],"Time Of Flight");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_P)
          strcpy(Short_Label[i],"Intensity");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER1)
          strncpy(Short_Label[i],Vars->UserName1,30);
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER2)
          strncpy(Short_Label[i],Vars->UserName2,30);
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER3)
          strncpy(Short_Label[i],Vars->UserName3,30);
      else
          strcpy(Short_Label[i],"Unknown");

      if (Vars->Coord_Type[i] & DEFS->COORD_ABS)
      { strcat(Vars->Coord_Label[i]," (abs)"); }

      if (Vars->Coord_Type[i] & DEFS->COORD_LOG)
      { strcat(Vars->Coord_Label[i]," (log)"); }

      strcat(Vars->Monitor_Label, " ");
      strcat(Vars->Monitor_Label, Short_Label[i]);
      XY *= Vars->Coord_Bin[i];
    } /* end for Short_Label */

    if ((Vars->Coord_Type[0] & (DEFS->COORD_LOG-1)) == DEFS->COORD_P) {
      strncat(Vars->Coord_Label[0], " [n/s", 30);
      if (Vars->Flag_per_cm2) strncat(Vars->Coord_Label[0], "/cm2", 30);

      if (XY > 1 && Vars->Coord_Number)
        strncat(Vars->Coord_Label[0], "/bin", 30);
      strncat(Vars->Coord_Label[0], "]", 30);
    }

    /* update label 'signal per bin' if more than 1 bin */
    if (XY > 1 && Vars->Coord_Number) {
      if (Vars->Flag_capture)
        printf("Monitor_nD: %s: Using capture flux weightening on %ld bins.\n"
               "WARNING     Use binned data with caution, and prefer monitor integral value (I,Ierr).\n", Vars->compcurname, (long)XY);
    }

    strcat(Vars->Monitor_Label, " Monitor");
    if (Vars->Flag_Shape == DEFS->SHAPE_SQUARE) strcat(Vars->Monitor_Label, " (Square)");
    if (Vars->Flag_Shape == DEFS->SHAPE_DISK)   strcat(Vars->Monitor_Label, " (Disk)");
    if (Vars->Flag_Shape == DEFS->SHAPE_SPHERE) strcat(Vars->Monitor_Label, " (Sphere)");
    if (Vars->Flag_Shape == DEFS->SHAPE_CYLIND) strcat(Vars->Monitor_Label, " (Cylinder)");
    if (Vars->Flag_Shape == DEFS->SHAPE_BANANA) strcat(Vars->Monitor_Label, " (Banana)");
    if (Vars->Flag_Shape == DEFS->SHAPE_BOX)    strcat(Vars->Monitor_Label, " (Box)");
    if (Vars->Flag_Shape == DEFS->SHAPE_PREVIOUS) strcat(Vars->Monitor_Label, " (on PREVIOUS)");
    if ((Vars->Flag_Shape == DEFS->SHAPE_CYLIND) || (Vars->Flag_Shape == DEFS->SHAPE_BANANA) || (Vars->Flag_Shape == DEFS->SHAPE_SPHERE) || (Vars->Flag_Shape == DEFS->SHAPE_BOX))
    {
      if (strstr(Vars->option, "incoming"))
      {
        Vars->Flag_Shape = abs(Vars->Flag_Shape);
        strcat(Vars->Monitor_Label, " [in]");
      }
      else /* if strstr(Vars->option, "outgoing")) */
      {
        Vars->Flag_Shape = -abs(Vars->Flag_Shape);
        strcat(Vars->Monitor_Label, " [out]");
      }
    }
    if (Vars->Flag_UsePreMonitor == 1)
    {
        strcat(Vars->Monitor_Label, " at ");
        strncat(Vars->Monitor_Label, Vars->UserName1,30);
    }
    if (Vars->Flag_log == 1) strcat(Vars->Monitor_Label, " [log] ");

    /* now allocate memory to store variables in TRACE */

    /* Vars->Coord_Number  0   : intensity or signal
     * Vars->Coord_Number  1:n : detector variables */

    if ((Vars->Coord_Number != 2) && !Vars->Flag_Multiple && !Vars->Flag_List)
    { Vars->Flag_Multiple = 1; Vars->Flag_List = 0; } /* default is n1D */

    /* list and auto limits case : Vars->Flag_List or Vars->Flag_Auto_Limits
     * -> Buffer to flush and suppress after Vars->Flag_Auto_Limits
     */
    if ((Vars->Flag_Auto_Limits || Vars->Flag_List) && Vars->Coord_Number)
    { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
      Vars->Mon2D_Buffer = (double *)malloc((Vars->Coord_Number+1)*Vars->Buffer_Block*sizeof(double));
      if (Vars->Mon2D_Buffer == NULL)
      { printf("Monitor_nD: %s cannot allocate Vars->Mon2D_Buffer (%li). No list and auto limits.\n", Vars->compcurname, Vars->Buffer_Block*(Vars->Coord_Number+1)*sizeof(double)); Vars->Flag_List = 0; Vars->Flag_Auto_Limits = 0; }
      else
      {
        for (i=0; i < (Vars->Coord_Number+1)*Vars->Buffer_Block; Vars->Mon2D_Buffer[i++] = (double)0);
      }
      Vars->Buffer_Size = Vars->Buffer_Block;
    }

    /* 1D and n1D case : Vars->Flag_Multiple */
    if (Vars->Flag_Multiple && Vars->Coord_Number)
    { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors */
      Vars->Mon2D_N  = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      Vars->Mon2D_p  = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      Vars->Mon2D_p2 = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
      { fprintf(stderr,"Monitor_nD: %s n1D cannot allocate Vars->Mon2D_N/p/p2 (%li). Fatal.\n", Vars->compcurname, (Vars->Coord_Number)*sizeof(double *)); exit(-1); }
      for (i= 1; i <= Vars->Coord_Number; i++)
      {
        Vars->Mon2D_N[i-1]  = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        Vars->Mon2D_p[i-1]  = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        Vars->Mon2D_p2[i-1] = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
        { fprintf(stderr,"Monitor_nD: %s n1D cannot allocate %s Vars->Mon2D_N/p/p2[%li] (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[i], i, (Vars->Coord_Bin[i])*sizeof(double *)); exit(-1); }
        else
        {
          for (j=0; j < Vars->Coord_Bin[i]; j++ )
          { Vars->Mon2D_N[i-1][j] = (double)0; Vars->Mon2D_p[i-1][j] = (double)0; Vars->Mon2D_p2[i-1][j] = (double)0; }
        }
      }
    }
    else /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
    if ((Vars->Coord_Number == 2) && !Vars->Flag_Multiple)
    { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
      Vars->Mon2D_N  = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      Vars->Mon2D_p  = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      Vars->Mon2D_p2 = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
      { fprintf(stderr,"Monitor_nD: %s 2D cannot allocate %s Vars->Mon2D_N/p/p2 (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[1], (Vars->Coord_Bin[1])*sizeof(double *)); exit(-1); }
      for (i= 0; i < Vars->Coord_Bin[1]; i++)
      {
        Vars->Mon2D_N[i]  = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        Vars->Mon2D_p[i]  = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        Vars->Mon2D_p2[i] = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
        { fprintf(stderr,"Monitor_nD: %s 2D cannot allocate %s Vars->Mon2D_N/p/p2[%li] (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[1], i, (Vars->Coord_Bin[2])*sizeof(double *)); exit(-1); }
        else
        {
          for (j=0; j < Vars->Coord_Bin[2]; j++ )
          { Vars->Mon2D_N[i][j] = (double)0; Vars->Mon2D_p[i][j] = (double)0; Vars->Mon2D_p2[i][j] = (double)0; }
        }
      }
    }
      /* no Mon2D allocated for
       * (Vars->Coord_Number != 2) && !Vars->Flag_Multiple && Vars->Flag_List */

    Vars->psum  = 0;
    Vars->p2sum = 0;
    Vars->Nsum  = 0;

    Vars->area  = fabs(Vars->mxmax - Vars->mxmin)*fabs(Vars->mymax - Vars->mymin)*1E4; /* in cm**2 for square and box shapes */
    Vars->Sphere_Radius = fabs(Vars->mxmax - Vars->mxmin)/2;
    if ((abs(Vars->Flag_Shape) == DEFS->SHAPE_DISK) || (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE))
    {
      Vars->area = PI*Vars->Sphere_Radius*Vars->Sphere_Radius*1E4; /* disk shapes */
    }
    if (Vars->area == 0 && abs(Vars->Flag_Shape) != DEFS->SHAPE_PREVIOUS)
      Vars->Coord_Number = 0;
    if (Vars->Coord_Number == 0 && Vars->Flag_Verbose)
      printf("Monitor_nD: %s is unactivated (0D)\n", Vars->compcurname);
    Vars->Cylinder_Height = fabs(Vars->mymax - Vars->mymin);

    if (Vars->Flag_Verbose)
    {
      printf("Monitor_nD: %s is a %s.\n", Vars->compcurname, Vars->Monitor_Label);
      printf("Monitor_nD: version %s with options=%s\n", MONITOR_ND_LIB_H, Vars->option);
    }
  } /* end Monitor_nD_Init */

/* ========================================================================= */
/* Monitor_nD_Trace: this routine is used to monitor one propagating neutron */
/* ========================================================================= */

double Monitor_nD_Trace(MonitornD_Defines_type *DEFS, MonitornD_Variables_type *Vars)
{

  double  XY=0;
  long    i=0,j;
  double  pp=0;
  double  Coord[MONnD_COORD_NMAX];
  long    Coord_Index[MONnD_COORD_NMAX];
  char    While_End   =0;
  long    While_Buffer=0;
  char    Set_Vars_Coord_Type = DEFS->COORD_NONE;

  /* Vars->Flag_Auto_Limits: we read the Buffer, and determine min and max bounds */
  if ((Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 1) && (Vars->Coord_Number > 0))
  {
    /* auto limits case : get limits in Buffer for each variable */
          /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
    if (Vars->Flag_Verbose) printf("Monitor_nD: %s getting %li Auto Limits from List (%li) in TRACE.\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
    for (i = 1; i <= Vars->Coord_Number; i++)
    {
      if (Vars->Coord_Type[i] & DEFS->COORD_AUTO)
      {
        Vars->Coord_Min[i] =  FLT_MAX;
        Vars->Coord_Max[i] = -FLT_MAX;
        for (j = 0; j < Vars->Buffer_Block; j++)
        {
          XY = Vars->Mon2D_Buffer[i+j*(Vars->Coord_Number+1)];  /* scanning variables in Buffer */
          if (XY < Vars->Coord_Min[i]) Vars->Coord_Min[i] = XY;
          if (XY > Vars->Coord_Max[i]) Vars->Coord_Max[i] = XY;
        }
      }
    }
    Vars->Flag_Auto_Limits = 2;  /* pass to 2nd auto limits step (read Buffer and generate new events to store in histograms) */
  } /* end if Flag_Auto_Limits == 1 */

  /* manage realloc for 'list all' if Buffer size exceeded: flush Buffer to file */
  if ((Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_List >= 2))
  {
    if (Vars->Buffer_Size >= 20000 || Vars->Flag_List == 3)
    { /* save current (possibly append) and re-use Buffer */
      Monitor_nD_Save(DEFS, Vars);
      Vars->Flag_List = 3;
      Vars->Buffer_Block = Vars->Buffer_Size;
      Vars->Buffer_Counter  = 0;
      Vars->Neutron_Counter = 0;
    }
    else
    {
      Vars->Mon2D_Buffer  = (double *)realloc(Vars->Mon2D_Buffer, (Vars->Coord_Number+1)*(Vars->Neutron_Counter+Vars->Buffer_Block)*sizeof(double));
      if (Vars->Mon2D_Buffer == NULL)
            { printf("Monitor_nD: %s cannot reallocate Vars->Mon2D_Buffer[%li] (%li). Skipping.\n", Vars->compcurname, i, (Vars->Neutron_Counter+Vars->Buffer_Block)*sizeof(double)); Vars->Flag_List = 1; }
      else { Vars->Buffer_Counter = 0; Vars->Buffer_Size = Vars->Neutron_Counter+Vars->Buffer_Block; }
    }
  } /* end if Buffer realloc */

  while (!While_End)
  { /* we generate Coord[] and Coord_index[] from Buffer (auto limits) or passing neutron */
    if ((Vars->Flag_Auto_Limits == 2) && (Vars->Coord_Number > 0))
    { /* Vars->Flag_Auto_Limits == 2: read back from Buffer (Buffer is filled or auto limits have been computed) */
      if (While_Buffer < Vars->Buffer_Block)
      {
        /* first while loops (While_Buffer) */
        /* auto limits case : scan Buffer within limits and store in Mon2D */
        pp = Vars->Mon2D_Buffer[While_Buffer*(Vars->Coord_Number+1)];
        /* For some reason the Intel c compiler version 10.1 gives 0 counts with Monitor_nD!
           An ugly patch seems to be the following printf */
        printf("%s", "");
        Coord[0] = pp;

        for (i = 1; i <= Vars->Coord_Number; i++)
        {
          /* scanning variables in Buffer */
          XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);

          Coord[i] = Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)];
          if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
          else Coord_Index[i] = 0;
          if (Vars->Flag_With_Borders)
          {
            if (Coord_Index[i] < 0) Coord_Index[i] = 0;
            if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
          }
        } /* end for */
        While_Buffer++;
      } /* end if in Buffer */
      else /* (While_Buffer >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 2) */
      {
        Vars->Flag_Auto_Limits = 0;
        if (!Vars->Flag_List) /* free Buffer not needed anymore (no list to output) */
        { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, p2) */
          free(Vars->Mon2D_Buffer); Vars->Mon2D_Buffer = NULL;
        }
        if (Vars->Flag_Verbose) printf("Monitor_nD: %s flushed %li Auto Limits from List (%li) in TRACE.\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
      }
    }
    if (Vars->Flag_Auto_Limits != 2 || !Vars->Coord_Number) /* Vars->Flag_Auto_Limits == 0 (no auto limits/list) or 1 (store events into Buffer) */
    {
      /* automatically compute area and steradian solid angle when in AUTO mode */
      /* compute the steradian solid angle incoming on the monitor */
      double v;
      v=sqrt(Vars->cvx*Vars->cvx
            +Vars->cvy*Vars->cvy
            +Vars->cvz*Vars->cvz);
      if (Vars->min_x > Vars->cx) Vars->min_x = Vars->cx;
      if (Vars->max_x < Vars->cx) Vars->max_x = Vars->cx;
      if (Vars->min_y > Vars->cy) Vars->min_y = Vars->cy;
      if (Vars->max_y < Vars->cy) Vars->max_y = Vars->cy;
      Vars->mean_p  += Vars->cp;
      if (v) {
        Vars->mean_dx += Vars->cp*fabs(Vars->cvx/v);
        Vars->mean_dy += Vars->cp*fabs(Vars->cvy/v);
      }

      for (i = 0; i <= Vars->Coord_Number; i++)
      { /* handle current neutron : last while */
        XY = 0;
        Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
        /* get values for variables to monitor */
        if (Set_Vars_Coord_Type == DEFS->COORD_X) XY = Vars->cx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_Y) XY = Vars->cy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_Z) XY = Vars->cz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VX) XY = Vars->cvx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VY) XY = Vars->cvy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VZ) XY = Vars->cvz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KX) XY = V2K*Vars->cvx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KY) XY = V2K*Vars->cvy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KZ) XY = V2K*Vars->cvz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SX) XY = Vars->csx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SY) XY = Vars->csy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SZ) XY = Vars->csz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_T) XY = Vars->ct;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_P) XY = Vars->cp;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_HDIV) XY = RAD2DEG*atan2(Vars->cvx,Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VDIV) XY = RAD2DEG*atan2(Vars->cvy,Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_V) XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_RADIUS)
          XY = sqrt(Vars->cx*Vars->cx+Vars->cy*Vars->cy+Vars->cz*Vars->cz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_XY)
          XY = sqrt(Vars->cx*Vars->cx+Vars->cy*Vars->cy)*(Vars->cx > 0 ? 1 : -1);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_YZ) XY = sqrt(Vars->cy*Vars->cy+Vars->cz*Vars->cz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_XZ)
          XY = sqrt(Vars->cx*Vars->cx+Vars->cz*Vars->cz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VXY) XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VXZ) XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvz*Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VYZ) XY = sqrt(Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_K) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KXY) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KXZ) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvz*Vars->cvz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KYZ) { XY = sqrt(Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_ENERGY) { XY = Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz;  XY *= VS2E; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_LAMBDA) { XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy+Vars->cvz*Vars->cvz);  XY *= V2K; if (XY != 0) XY = 2*PI/XY; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_NCOUNT) XY = Coord[i]+1;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_ANGLE)
        {  XY = sqrt(Vars->cvx*Vars->cvx+Vars->cvy*Vars->cvy);
           if (Vars->cvz != 0)
                XY = RAD2DEG*atan2(XY,Vars->cvz)*(Vars->cx > 0 ? 1 : -1);
           else XY = 0;
        }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_THETA)  { if (Vars->cz != 0) XY = RAD2DEG*atan2(Vars->cx,Vars->cz); }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_PHI) { if (Vars->cz != 0) XY = RAD2DEG*asin(Vars->cy/Vars->cz); }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER1) XY = Vars->UserVariable1;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER2) XY = Vars->UserVariable2;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER3) XY = Vars->UserVariable3;
        else
        XY = 0;
        /* handle 'abs' and 'log' keywords */
        if (Vars->Coord_Type[i] & DEFS->COORD_ABS) XY=fabs(XY);

        if (i && (Vars->Coord_Type[i] & DEFS->COORD_LOG)) /* not for the flux */
        {  if (XY > 0) XY = log(XY)/log(10);
           else XY = -100; }

        Coord[i] = XY;
        if (i == 0) { pp = XY; Coord_Index[i] = 0; }
        else if (!Vars->Flag_Auto_Limits)
        { /* compute index in histograms for each variable to monitor */
          XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);
          if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
          else Coord_Index[i] = 0;
          if (Vars->Flag_With_Borders)
          {
            if (Coord_Index[i] < 0) Coord_Index[i] = 0;
            if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
          }
        } /* else will get Index later from Buffer when Flag_Auto_Limits == 2 */
      } /* end for i */
      While_End = 1;
    }/* end else if Vars->Flag_Auto_Limits == 2 */

    if (Vars->Flag_Auto_Limits != 2) /* not when reading auto limits Buffer */
    { /* now store Coord into Buffer (no index needed) if necessary (list or auto limits) */
      if ((Vars->Buffer_Counter < Vars->Buffer_Block) && ((Vars->Flag_List) || (Vars->Flag_Auto_Limits == 1)))
      {
        for (i = 0; i <= Vars->Coord_Number; i++)
        {
          Vars->Mon2D_Buffer[i + Vars->Neutron_Counter*(Vars->Coord_Number+1)] = Coord[i];
        }
        Vars->Buffer_Counter++;
        if (Vars->Flag_Verbose && (Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_List == 1)) printf("Monitor_nD: %s %li neutrons stored in List.\n", Vars->compcurname, Vars->Buffer_Counter);
      }
      Vars->Neutron_Counter++;
    } /* end (Vars->Flag_Auto_Limits != 2) */

    /* ====================================================================== */
    /* store n1d/2d neutron from Buffer (Auto_Limits == 2) or current neutron in while */
    if (Vars->Flag_Auto_Limits != 1) /* not when storing auto limits Buffer */
    {
      /* apply per cm2 or per st */
      if (Vars->Flag_per_cm2 && Vars->area      != 0)
        pp /= Vars->area;

      /* 1D and n1D case : Vars->Flag_Multiple */
      if (Vars->Flag_Multiple)
      { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors (intensity is not included) */
        /* check limits: monitors define a phase space to record */
        char within_limits=1;
        for (i= 1; i <= Vars->Coord_Number; i++)
        {
          j = Coord_Index[i];
          if (j < 0 || j >= Vars->Coord_Bin[i])
            within_limits=0;
        }
        if (within_limits)
        { for (i= 1; i <= Vars->Coord_Number; i++)
          {
            j = Coord_Index[i];
            if (j >= 0 && j < Vars->Coord_Bin[i])
            {
              Vars->Mon2D_N[i-1][j]++;
              Vars->Mon2D_p[i-1][j] += pp;
              Vars->Mon2D_p2[i-1][j] += pp*pp;
            }
          }
        }
        else if (Vars->Flag_Exclusive)
        { pp = 0.0;
        }
      }
      else /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
      if ((Vars->Coord_Number == 2) && !Vars->Flag_Multiple)
      { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
        i = Coord_Index[1];
        j = Coord_Index[2];
        if (i >= 0 && i < Vars->Coord_Bin[1] && j >= 0 && j < Vars->Coord_Bin[2])
        {
          Vars->Mon2D_N[i][j]++;
          Vars->Mon2D_p[i][j] += pp;
          Vars->Mon2D_p2[i][j] += pp*pp;
        }
        else if (Vars->Flag_Exclusive)
        { pp = 0.0;
        }
      }
    } /* end (Vars->Flag_Auto_Limits != 1) */
  } /* end while */
  return pp;
} /* end Monitor_nD_Trace */

/* ========================================================================= */
/* Monitor_nD_Save: this routine is used to save data files                  */
/* ========================================================================= */

MCDETECTOR Monitor_nD_Save(MonitornD_Defines_type *DEFS, MonitornD_Variables_type *Vars)
  {
    char   *fname;
    long    i,j;
    double *p0m = NULL;
    double *p1m = NULL;
    double *p2m = NULL;
    char    Coord_X_Label[CHAR_BUF_LENGTH];
    double  min1d, max1d;
    double  min2d, max2d;
    long    bin1d, bin2d;
    char    While_End = 0;
    long    While_Buffer = 0;
    double  XY=0, pp;
    double  Coord[MONnD_COORD_NMAX];
    long    Coord_Index[MONnD_COORD_NMAX];
    char    label[CHAR_BUF_LENGTH];
    double  ratio;

    MCDETECTOR detector;

    ratio = 100.0*mcget_run_num()/mcget_ncount();
    if (Vars->Flag_Verbose && Vars->Flag_per_cm2) {
      printf("Monitor_nD: %s: active flat detector area is %g [cm^2], total area is %g [cm^2]\n",
        Vars->compcurname, (Vars->max_x-Vars->min_x)
                          *(Vars->max_y-Vars->min_y)*1E4, Vars->area);
      printf("Monitor_nD: %s: beam solid angle is %g [st] (%g x %g [deg^2])\n",
        Vars->compcurname,
        2*fabs(2*atan(Vars->mean_dx/Vars->mean_p)
         *sin(2*atan(Vars->mean_dy/Vars->mean_p)/2)),
        atan(Vars->mean_dx/Vars->mean_p)*RAD2DEG,
        atan(Vars->mean_dy/Vars->mean_p)*RAD2DEG);
    }

    /* check Buffer flush when end of simulation reached */
    if ((Vars->Buffer_Counter <= Vars->Buffer_Block) && Vars->Flag_Auto_Limits && Vars->Mon2D_Buffer && Vars->Buffer_Counter)
    {
      /* Get Auto Limits */
      if (Vars->Flag_Verbose) printf("Monitor_nD: %s getting %li Auto Limits from List (%li).\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);

      for (i = 1; i <= Vars->Coord_Number; i++)
      {
        if (Vars->Coord_Type[i] & DEFS->COORD_AUTO)
        {
          Vars->Coord_Min[i] = FLT_MAX;
          Vars->Coord_Max[i] = -FLT_MAX;
          for (j = 0; j < Vars->Buffer_Block; j++)
          {
            XY = Vars->Mon2D_Buffer[i+j*(Vars->Coord_Number+1)];  /* scanning variables in Buffer */
            if (XY < Vars->Coord_Min[i]) Vars->Coord_Min[i] = XY;
            if (XY > Vars->Coord_Max[i]) Vars->Coord_Max[i] = XY;
          }
        }
      }
      Vars->Flag_Auto_Limits = 2;  /* pass to 2nd auto limits step */
      Vars->Buffer_Block = Vars->Buffer_Counter;

      while (!While_End)
      { /* we generate Coord[] and Coord_index[] from Buffer (auto limits) */
        /* simulation ended before Buffer was filled. Limits have to be computed, and stored events must be sent into histograms */
        if (While_Buffer < Vars->Buffer_Block)
        {
          /* first while loops (While_Buffer) */
          Coord[0] = Vars->Mon2D_Buffer[While_Buffer*(Vars->Coord_Number+1)];

          /* auto limits case : scan Buffer within limits and store in Mon2D */
          for (i = 1; i <= Vars->Coord_Number; i++)
          {
            /* scanning variables in Buffer */
            XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);
            Coord[i] = Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)];
            if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
            else Coord_Index[i] = 0;
            if (Vars->Flag_With_Borders)
            {
              if (Coord_Index[i] < 0) Coord_Index[i] = 0;
              if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
            }
          } /* end for */
          While_Buffer++;
        } /* end if in Buffer */
        else /* (While_Buffer >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 2) */
        {
          Vars->Flag_Auto_Limits = 0;
          While_End = 1;
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s flushed %li Auto Limits from List (%li).\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
        }

        /* store n1d/2d section from Buffer */

        pp = Coord[0];
        /* apply per cm2 or per st */
        if (Vars->Flag_per_cm2 && Vars->area      != 0)
          pp /= Vars->area;

        /* 1D and n1D case : Vars->Flag_Multiple */
        if (Vars->Flag_Multiple)
        { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors (intensity is not included) */
          for (i= 0; i < Vars->Coord_Number; i++)
          {
            j = Coord_Index[i+1];
            if (j >= 0 && j < Vars->Coord_Bin[i+1])
            {
              Vars->Mon2D_N[i][j]++;
              Vars->Mon2D_p[i][j] += pp;
              Vars->Mon2D_p2[i][j] += pp*pp;
            }
          }
        }
        else /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
        if ((Vars->Coord_Number == 2) && !Vars->Flag_Multiple)
        { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
          i = Coord_Index[1];
          j = Coord_Index[2];
          if (i >= 0 && i < Vars->Coord_Bin[1] && j >= 0 && j < Vars->Coord_Bin[2])
          {
            Vars->Mon2D_N[i][j]++;
            Vars->Mon2D_p[i][j] += pp;
            Vars->Mon2D_p2[i][j] += pp*pp;
          }
        } /* end store 2D/1D */
      } /* end while */
    } /* end Force Get Limits */

    /* write output files (sent to file as p[i*n + j] vectors) */
    if (Vars->Coord_Number == 0)
    {
      double Nsum;
      double psum, p2sum;
      Nsum = Vars->Nsum;
      psum = Vars->psum;
      p2sum= Vars->p2sum;
      if (Vars->Flag_signal != DEFS->COORD_P && Nsum > 0)
      { psum /=Nsum; p2sum /= Nsum*Nsum; }
      /* DETECTOR_OUT_0D(Vars->Monitor_Label, Vars->Nsum, Vars->psum, Vars->p2sum); */
      detector = mcdetector_out_0D(Vars->Monitor_Label, Nsum, psum, p2sum, Vars->compcurname, Vars->compcurpos);
    }
    else
    if (strlen(Vars->Mon_File) > 0)
    {
      fname = (char*)malloc(strlen(Vars->Mon_File)+10*Vars->Coord_Number);
      if (Vars->Flag_List && Vars->Mon2D_Buffer) /* List: DETECTOR_OUT_2D */
      {
        int  ascii_only_orig;
        char formatName[256];
        char *formatName_orig;

        if (Vars->Flag_List >= 2) Vars->Buffer_Size = Vars->Neutron_Counter;
        if (Vars->Buffer_Size >= Vars->Neutron_Counter)
          Vars->Buffer_Size = Vars->Neutron_Counter;
        strcpy(fname,Vars->Mon_File);
        if (strchr(Vars->Mon_File,'.') == NULL) strcat(fname, "_list");

        min1d = 1; max1d = Vars->Coord_Number+1;
        min2d = 0; max2d = Vars->Buffer_Size;
        bin1d = Vars->Coord_Number+1; bin2d = Vars->Buffer_Size;
        strcpy(Coord_X_Label,"");
        for (i= 0; i <= Vars->Coord_Number; i++)
        {
          if (min2d < Vars->Coord_Min[i]) min2d = Vars->Coord_Min[i];
          if (max2d < Vars->Coord_Max[i]) max2d = Vars->Coord_Max[i];
          strcat(Coord_X_Label, Vars->Coord_Var[i]);
          strcat(Coord_X_Label, " ");
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[i]); }
        }
        if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s List (%lix%li).\n", Vars->compcurname, fname,bin2d,bin1d);

        /* handle the type of list output */
        ascii_only_orig = mcascii_only;
        formatName_orig = mcformat.Name;  /* copy the pointer position */
        strcpy(formatName, mcformat.Name);
        if (Vars->Flag_List >= 1)
        { /* Flag_List mode:
               1=store 1 buffer
               2=list all, triggers 3 when 1st buffer reallocated
               3=re-used buffer (file already opened)
             Format modifiers for Flag_List
               1= normal monitor file (no modifier, export in one go)
               2= write data+header, and footer if buffer not full (Vars->Buffer_Counter < Vars->Buffer_Block)
               3= write data, and footer if buffer not full (final)
           */
          strcat(formatName, " list ");
          if (Vars->Flag_List == 3) strcat(formatName, " no header ");
          if (Vars->Flag_List >= 2 && Vars->Buffer_Counter >= Vars->Buffer_Block)
            strcat(formatName, " no footer ");

          if (Vars->Flag_Binary_List) mcascii_only = 1;
          if (Vars->Flag_Binary_List == 1)
            strcat(formatName, " binary float ");
          else if (Vars->Flag_Binary_List == 2)
            strcat(formatName, " binary double ");
        }
        if (min2d == max2d) max2d = min2d+1e-6;
        if (min1d == max1d) max1d = min1d+1e-6;
        strcpy(label, Vars->Monitor_Label);
        if (!Vars->Flag_Binary_List)
        { bin2d=-bin2d; }
        mcformat.Name = formatName;
        detector = mcdetector_out_2D(
              label,
              "List of neutron events",
              Coord_X_Label,
              min2d, max2d,
              min1d, max1d,
              bin2d,
              bin1d,
            NULL,Vars->Mon2D_Buffer,NULL,
            fname, Vars->compcurname, Vars->compcurpos);

        /* reset the original type of output */
        mcascii_only = ascii_only_orig;
        mcformat.Name= formatName_orig;
      }
      if (Vars->Flag_Multiple) /* n1D: DETECTOR_OUT_1D */
      {
        for (i= 0; i < Vars->Coord_Number; i++)
        {

          strcpy(fname,Vars->Mon_File);
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[i+1]); }
          sprintf(Coord_X_Label, "%s monitor", Vars->Coord_Label[i+1]);
          strcpy(label, Coord_X_Label);
          if (Vars->Coord_Bin[i+1] > 0) { /* 1D monitor */
            if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s 1D (%li).\n", Vars->compcurname, fname, Vars->Coord_Bin[i+1]);
            min1d = Vars->Coord_Min[i+1];
            max1d = Vars->Coord_Max[i+1];
            if (min1d == max1d) max1d = min1d+1e-6;
            p1m = (double *)malloc(Vars->Coord_Bin[i+1]*sizeof(double));
            p2m = (double *)malloc(Vars->Coord_Bin[i+1]*sizeof(double));
            if (p2m == NULL) /* use Raw Buffer line output */
            {
              if (Vars->Flag_Verbose) printf("Monitor_nD: %s cannot allocate memory for output. Using raw data.\n", Vars->compcurname);
              if (p1m != NULL) free(p1m);
              detector = mcdetector_out_1D(
              label,
              Vars->Coord_Label[i+1],
              Vars->Coord_Label[0],
              Vars->Coord_Var[i+1],
              min1d, max1d,
              Vars->Coord_Bin[i+1],
              Vars->Mon2D_N[i],Vars->Mon2D_p[i],Vars->Mon2D_p2[i],
              fname, Vars->compcurname, Vars->compcurpos);
            } /* if (p2m == NULL) */
            else
            {
              if (Vars->Flag_log != 0)
              {
                XY = FLT_MAX;
                for (j=0; j < Vars->Coord_Bin[i+1]; j++) /* search min of signal */
                  if ((XY > Vars->Mon2D_p[i][j]) && (Vars->Mon2D_p[i][j] > 0)) XY = Vars->Mon2D_p[i][j];
                if (XY <= 0) XY = -log(FLT_MAX)/log(10); else XY = log(XY)/log(10)-1;
              } /* if */

              for (j=0; j < Vars->Coord_Bin[i+1]; j++)
              {
                p1m[j] = Vars->Mon2D_p[i][j];
                p2m[j] = Vars->Mon2D_p2[i][j];
                if (Vars->Flag_signal != DEFS->COORD_P && Vars->Mon2D_N[i][j] > 0)
                { /* normalize mean signal to the number of events */
                  p1m[j] /= Vars->Mon2D_N[i][j];
                  p2m[j] /= Vars->Mon2D_N[i][j]*Vars->Mon2D_N[i][j];
                }
                if (Vars->Flag_log != 0)
                {
                  if ((p1m[j] > 0) && (p2m[j] > 0))
                  {
                    p2m[j] /= p1m[j]*p1m[j];
                    p1m[j] = log(p1m[j])/log(10);
                  }
                  else
                  {
                    p1m[j] = XY;
                    p2m[j] = 0;
                  }
                }
              } /* for */
              detector = mcdetector_out_1D(
                label,
                Vars->Coord_Label[i+1],
                Vars->Coord_Label[0],
                Vars->Coord_Var[i+1],
                min1d, max1d,
                Vars->Coord_Bin[i+1],
                Vars->Mon2D_N[i],p1m,p2m,
                fname, Vars->compcurname, Vars->compcurpos);

            } /* else */
            /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
            if (p1m != NULL) free(p1m); p1m=NULL;
            if (p2m != NULL) free(p2m); p2m=NULL;
            */
          } else { /* 0d monitor */
            detector = mcdetector_out_0D(label, Vars->Mon2D_p[i][0], Vars->Mon2D_p2[i][0], Vars->Mon2D_N[i][0], Vars->compcurname, Vars->compcurpos);
          }


        } /* for */
      } /* if 1D */
      else
      if (Vars->Coord_Number == 2)  /* 2D: DETECTOR_OUT_2D */
      {
        strcpy(fname,Vars->Mon_File);

        p0m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        p1m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        p2m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        if (p2m == NULL)
        {
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s cannot allocate memory for 2D array (%li). Skipping.\n", Vars->compcurname, 3*Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
          /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
          if (p0m != NULL) free(p0m);
          if (p1m != NULL) free(p1m);
          */
        }
        else
        {
          if (Vars->Flag_log != 0)
          {
            XY = FLT_MAX;
            for (i= 0; i < Vars->Coord_Bin[1]; i++)
              for (j= 0; j < Vars->Coord_Bin[2]; j++) /* search min of signal */
                if ((XY > Vars->Mon2D_p[i][j]) && (Vars->Mon2D_p[i][j]>0)) XY = Vars->Mon2D_p[i][j];
            if (XY <= 0) XY = -log(FLT_MAX)/log(10); else XY = log(XY)/log(10)-1;
          }
          for (i= 0; i < Vars->Coord_Bin[1]; i++)
          {
            for (j= 0; j < Vars->Coord_Bin[2]; j++)
            {
              long index;
              index = j + i*Vars->Coord_Bin[2];
              p0m[index] = Vars->Mon2D_N[i][j];
              p1m[index] = Vars->Mon2D_p[i][j];
              p2m[index] = Vars->Mon2D_p2[i][j];
              if (Vars->Flag_signal != DEFS->COORD_P && p0m[index] > 0)
              {
                  p1m[index] /= p0m[index];
                  p2m[index] /= p0m[index]*p0m[index];
              }

              if (Vars->Flag_log != 0)
              {
                if ((p1m[index] > 0) && (p2m[index] > 0))
                {
                  p2m[index] /= (p1m[index]*p1m[index]);
                  p1m[index] = log(p1m[index])/log(10);

                }
                else
                {
                  p1m[index] = XY;
                  p2m[index] = 0;
                }
              }
            }
          }
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[1]);
              strcat(fname, "_"); strcat(fname, Vars->Coord_Var[2]); }
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s 2D (%lix%li).\n", Vars->compcurname, fname, Vars->Coord_Bin[1], Vars->Coord_Bin[2]);

          min1d = Vars->Coord_Min[1];
          max1d = Vars->Coord_Max[1];
          if (min1d == max1d) max1d = min1d+1e-6;
          min2d = Vars->Coord_Min[2];
          max2d = Vars->Coord_Max[2];
          if (min2d == max2d) max2d = min2d+1e-6;
          strcpy(label, Vars->Monitor_Label);
          if (Vars->Coord_Bin[1]*Vars->Coord_Bin[2] > 1
           && Vars->Flag_signal == DEFS->COORD_P)
            strcat(label, " per bin");

          detector = mcdetector_out_2D(
            label,
            Vars->Coord_Label[1],
            Vars->Coord_Label[2],
            min1d, max1d,
            min2d, max2d,
            Vars->Coord_Bin[1],
            Vars->Coord_Bin[2],
            p0m,p1m,p2m,
            fname, Vars->compcurname, Vars->compcurpos);

          /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
          if (p0m != NULL) free(p0m);
          if (p1m != NULL) free(p1m);
          if (p2m != NULL) free(p2m);
          */
        }
      }
      free(fname);
    }
    return(detector);
  } /* end Monitor_nD_Save */

/* ========================================================================= */
/* Monitor_nD_Finally: this routine is used to free memory                   */
/* ========================================================================= */

void Monitor_nD_Finally(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars)
  {
    int i;

    /* Now Free memory Mon2D.. */
    if ((Vars->Flag_Auto_Limits || Vars->Flag_List) && Vars->Coord_Number)
    { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
      if (Vars->Mon2D_Buffer != NULL) free(Vars->Mon2D_Buffer);
    }

    /* 1D and n1D case : Vars->Flag_Multiple */
    if (Vars->Flag_Multiple && Vars->Coord_Number)
    { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors */
      for (i= 0; i < Vars->Coord_Number; i++)
      {
        free(Vars->Mon2D_N[i]);
        free(Vars->Mon2D_p[i]);
        free(Vars->Mon2D_p2[i]);
      }
      free(Vars->Mon2D_N);
      free(Vars->Mon2D_p);
      free(Vars->Mon2D_p2);
    }


    /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
    if ((Vars->Coord_Number == 2) && !Vars->Flag_Multiple)
    { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
      for (i= 0; i < Vars->Coord_Bin[1]; i++)
      {
        free(Vars->Mon2D_N[i]);
        free(Vars->Mon2D_p[i]);
        free(Vars->Mon2D_p2[i]);
      }
      free(Vars->Mon2D_N);
      free(Vars->Mon2D_p);
      free(Vars->Mon2D_p2);
    }
  } /* end Monitor_nD_Finally */

/* ========================================================================= */
/* Monitor_nD_McDisplay: this routine is used to display component           */
/* ========================================================================= */

void Monitor_nD_McDisplay(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars)
  {
    double radius, h;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    int    i;
    double hdiv_min=-180, hdiv_max=180, vdiv_min=-180, vdiv_max=180;
    char   restricted = 0;

    radius = Vars->Sphere_Radius;
    h = Vars->Cylinder_Height;
    xmin = Vars->mxmin;
    xmax = Vars->mxmax;
    ymin = Vars->mymin;
    ymax = Vars->mymax;
    zmin = Vars->mzmin;
    zmax = Vars->mzmax;

    /* determine if there are angular limits set at start (no auto) in coord_types
     * cylinder/banana: look for hdiv
     * sphere: look for angle, radius (->atan2(val,radius)), hdiv, vdiv
     * this activates a 'restricted' flag, to draw a region as blades on cylinder/sphere
     */
    for (i= 0; i <= Vars->Coord_Number; i++)
    {
      int Set_Vars_Coord_Type;
      Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
      if (Set_Vars_Coord_Type == DEFS->COORD_HDIV || Set_Vars_Coord_Type == DEFS->COORD_THETA)
      { hdiv_min = Vars->Coord_Min[i]; hdiv_max = Vars->Coord_Max[i]; restricted = 1; }
      else if (Set_Vars_Coord_Type == DEFS->COORD_VDIV || Set_Vars_Coord_Type == DEFS->COORD_PHI)
      { vdiv_min = Vars->Coord_Min[i]; vdiv_max = Vars->Coord_Max[i];restricted = 1;  }
      else if (Set_Vars_Coord_Type == DEFS->COORD_ANGLE)
      { hdiv_min = vdiv_min = Vars->Coord_Min[i];
        hdiv_max = vdiv_max = Vars->Coord_Max[i];
        restricted = 1; }
      else if (Set_Vars_Coord_Type == DEFS->COORD_RADIUS)
      { double angle;
        angle = RAD2DEG*atan2(Vars->Coord_Max[i], radius);
        hdiv_min = vdiv_min = angle;
        hdiv_max = vdiv_max = angle;
        restricted = 1; }
    }
    /* full sphere */
    if ((!restricted && (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE))
    || abs(Vars->Flag_Shape) == DEFS->SHAPE_PREVIOUS)
    {
      mcdis_magnify("");
      mcdis_circle("xy",0,0,0,radius);
      mcdis_circle("xz",0,0,0,radius);
      mcdis_circle("yz",0,0,0,radius);
    }
    /* banana/cylinder/sphere portion */
    else
    if (restricted && ((abs(Vars->Flag_Shape) == DEFS->SHAPE_CYLIND)
                    || (abs(Vars->Flag_Shape) == DEFS->SHAPE_BANANA)
                    || (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE)))
    {
      int NH=24, NV=24;
      int ih, iv;
      double width, height;
      int issphere;
      issphere = (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE);
      width = (hdiv_max-hdiv_min)/NH;
      if (!issphere) NV=1;
      else height= (vdiv_max-vdiv_min)/NV;
      mcdis_magnify("xyz");
      for(ih = 0; ih < NH; ih++)
        for(iv = 0; iv < NV; iv++)
        {
          double theta0, phi0, theta1, phi1;
          double x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3;
          phi0 = (hdiv_min+ width*ih)*DEG2RAD; /* in xz plane */
          phi1 = (hdiv_min+ width*(ih+1))*DEG2RAD;
          if (issphere)
          {
            theta0= (90-vdiv_min+height* iv)   *DEG2RAD;
            theta1= (90-vdiv_min+height*(iv+1))*DEG2RAD;
            y0    = radius*cos(theta0);
            y1    = radius*cos(theta1);
          } else {
            y0 = ymin;
            y1 = ymax;
            theta0=theta1=90*DEG2RAD;
          }

          z0 = radius*sin(theta0)*cos(phi0);
          x0 = radius*sin(theta0)*sin(phi0);
          z1 = radius*sin(theta1)*cos(phi0);
          x1 = radius*sin(theta1)*sin(phi0);
          z2 = radius*sin(theta1)*cos(phi1);
          x2 = radius*sin(theta1)*sin(phi1);
          y2 = y1;
          z3 = radius*sin(theta0)*cos(phi1);
          x3 = radius*sin(theta0)*sin(phi1);
          y3 = y0;
          mcdis_multiline(5,
            x0,y0,z0,
            x1,y1,z1,
            x2,y2,z2,
            x3,y3,z3,
            x0,y0,z0);
        }
    }
    /* disk (circle) */
    else
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_DISK)
    {
      mcdis_magnify("");
      mcdis_circle("xy",0,0,0,radius);
    }
    /* rectangle (square) */
    else
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_SQUARE)
    {
      mcdis_magnify("xy");
      mcdis_multiline(5, (double)xmin, (double)ymin, 0.0,
             (double)xmax, (double)ymin, 0.0,
             (double)xmax, (double)ymax, 0.0,
             (double)xmin, (double)ymax, 0.0,
             (double)xmin, (double)ymin, 0.0);
    }
    /* full cylinder/banana */
    else
    if (!restricted && ((abs(Vars->Flag_Shape) == DEFS->SHAPE_CYLIND) || (abs(Vars->Flag_Shape) == DEFS->SHAPE_BANANA)))
    {
      mcdis_magnify("xyz");
      mcdis_circle("xz", 0,  h/2.0, 0, radius);
      mcdis_circle("xz", 0, -h/2.0, 0, radius);
      mcdis_line(-radius, -h/2.0, 0, -radius, +h/2.0, 0);
      mcdis_line(+radius, -h/2.0, 0, +radius, +h/2.0, 0);
      mcdis_line(0, -h/2.0, -radius, 0, +h/2.0, -radius);
      mcdis_line(0, -h/2.0, +radius, 0, +h/2.0, +radius);
    }
    else
    /* box */
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_BOX)
    {
      mcdis_magnify("xyz");
      mcdis_multiline(5, xmin, ymin, zmin,
                   xmax, ymin, zmin,
                   xmax, ymax, zmin,
                   xmin, ymax, zmin,
                   xmin, ymin, zmin);
      mcdis_multiline(5, xmin, ymin, zmax,
                   xmax, ymin, zmax,
                   xmax, ymax, zmax,
                   xmin, ymax, zmax,
                   xmin, ymin, zmax);
      mcdis_line(xmin, ymin, zmin, xmin, ymin, zmax);
      mcdis_line(xmax, ymin, zmin, xmax, ymin, zmax);
      mcdis_line(xmin, ymax, zmin, xmin, ymax, zmax);
      mcdis_line(xmax, ymax, zmin, xmax, ymax, zmax);
    }
  } /* end Monitor_nD_McDisplay */

/* end of monitor_nd-lib.c */


/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2008, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/interoff.h
*
* %Identification
* Written by: Reynald Arnerin
* Date:    Jun 12, 2008
* Release: 
* Version: 
*
* Object File Format intersection header for McStas. Requires the qsort function.
*
* Such files may be obtained with e.g.
*   qhull < points.xyz Qx Qv Tv o > points.off
* where points.xyz has format:
*   3
*   <nb_points>
*   <x> <y> <z>
*   ...
* The resulting file should have its first line being changed from '3' into 'OFF'.
* It can then be displayed with geomview.
* A similar, but somewhat older solution is to use 'powercrust' with e.g.
*   powercrust -i points.xyz
* which will generate a 'pc.off' file to be renamed as suited.
*
*******************************************************************************/

#ifndef INTEROFF_LIB_H
#define INTEROFF_LIB_H "$Revision: 1.4 $"

#ifndef EPSILON
#define EPSILON 1e-13
#endif

//#include <float.h>

#define N_VERTEX_DISPLAYED    2000

typedef struct intersection {
	MCNUM time;  	  //time of the intersection
	Coords v;	      //intersection point
	Coords normal;  //normal vector of the surface intersected
	short in_out;	  //1 if the ray enters the volume, -1 otherwise
	short edge;	    //1 if the intersection is on the boundary of the polygon, and error is possible
} intersection;

typedef struct polygon {
  MCNUM* p;       //vertices of the polygon in adjacent order, this way : x1 | y1 | z1 | x2 | y2 | z2 ...
  int npol;       //number of vertices
  Coords normal;
} polygon;

typedef struct off_struct {
    long vtxSize;
    long polySize;
    long faceSize;
    Coords* vtxArray;
    Coords* normalArray;
    unsigned long* faceArray;
} off_struct;

/*******************************************************************************
* long off_init(  char *offfile, double xwidth, double yheight, double zdepth, off_struct* data)
* ACTION: read an OFF file, optionally center object and rescale, initialize OFF data structure
* INPUT: 'offfile' OFF file to read
*        'xwidth,yheight,zdepth' if given as non-zero, apply bounding box. 
*           Specifying only one of these will also use the same ratio on all axes
*        'notcenter' center the object to the (0,0,0) position in local frame when set to zero
* RETURN: number of polyhedra and 'data' OFF structure 
*******************************************************************************/
long off_init(  char *offfile, double xwidth, double yheight, double zdepth, 
                int notcenter, off_struct* data);

/*******************************************************************************
* int off_intersect(double* t0, double* t3, 
     Coords *n0, Coords *n3,
     double x, double y, double z, 
     double vx, double vy, double vz, 
     off_struct data )
* ACTION: computes intersection of neutron trajectory with an object. 
* INPUT:  x,y,z and vx,vy,vz are the position and velocity of the neutron
*         data points to the OFF data structure
* RETURN: the number of polyhedra which trajectory intersects
*         t0 and t3 are the smallest incoming and outgoing intersection times
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
int off_intersect(double* t0, double* t3, 
     Coords *n0, Coords *n3,
     double x, double y, double z, 
     double vx, double vy, double vz, 
     off_struct data );

/*****************************************************************************
* int off_intersectx(double* l0, double* l3, 
     Coords *n0, Coords *n3,
     double x, double y, double z, 
     double kx, double ky, double kz, 
     off_struct data )
* ACTION: computes intersection of an xray trajectory with an object.
* INPUT:  x,y,z and kx,ky,kz, are spatial coordinates and wavevector of the x-ray
*         respectively. data points to the OFF data structure.
* RETURN: the number of polyhedra the trajectory intersects
*         l0 and l3 are the smallest incoming and outgoing intersection lengths
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
int off_x_intersect(double *l0,double *l3,
     Coords *n0, Coords *n3,
     double x,  double y,  double z, 
     double kx, double ky, double kz, 
     off_struct data );

/*******************************************************************************
* void off_display(off_struct data)
* ACTION: display up to N_VERTEX_DISPLAYED points from the object
*******************************************************************************/
void off_display(off_struct);

#endif

/* end of interoff-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2008, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/interoff-lib.c
*
* %Identification
* Written by: Reynald Arnerin
* Date:    Jun 12, 2008
* Origin: ILL
* Release: $Revision: 1.8 $
* Version: McStas X.Y
*
* Object File Format intersection library for McStas. Requires the qsort function.
*
* Such files may be obtained with e.g.
*   qhull < points.xyz Qx Qv Tv o > points.off
* where points.xyz has format (it supports comments):
*   3
*   <nb_points>
*   <x> <y> <z>
*   ...
* The resulting file should have its first line being changed from '3' into 'OFF'.
* It can then be displayed with geomview.
* A similar, but somewhat older solution is to use 'powercrust' with e.g.
*   powercrust -i points.xyz
* which will generate a 'pc.off' file to be renamed as suited.
*
*******************************************************************************/

#ifndef INTEROFF_LIB_H
#include "interoff-lib.h"
#endif

double off_F(double x, double y,double z,double A,double B,double C,double D) {
  return ( A*x + B*y + C*z + D );
}

char off_sign(double a) {
  if (a<0)       return(-1);
  else if (a==0) return(0);
  else           return(1);
}

// off_normal ******************************************************************
//gives the normal vector of p
void off_normal(Coords* n, polygon p)
{
  //using Newell method
  int i=0,j=0;
  n->x=0;n->y=0;n->z=0;
  for (i = 0, j = p.npol-1; i < p.npol; j = i++)
  {
    MCNUM x1=p.p[3*i],
          y1=p.p[3*i+1],
          z1=p.p[3*i+2];
    MCNUM x2=p.p[3*j],
          y2=p.p[3*j+1],
          z2=p.p[3*j+2];
    // n is the cross product of v1*v2
    n->x += (y1 - y2) * (z1 + z2);
    n->y += (z1 - z2) * (x1 + x2);
    n->z += (x1 - x2) * (y1 + y2);
  }
} /* off_normal */

// off_pnpoly ******************************************************************
//based on http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
//return 0 if the vertex is out
//    1 if it is in
//   -1 if on the boundary
int off_pnpoly(polygon p, Coords v)
{
  int i=0, c = 0;
  MCNUM minx=FLT_MAX,maxx=-FLT_MAX,miny=FLT_MAX,maxy=-FLT_MAX,minz=FLT_MAX,maxz=-FLT_MAX;
  MCNUM rangex=0,rangey=0,rangez=0;

  int pol2dx=0,pol2dy=1;          //2d restriction of the poly
  MCNUM x=v.x,y=v.y;


  //take the most relevant 2D projection (prevent from instability)
  for (i=0; i<p.npol; ++i)
  {
    if (p.p[3*i]<minx)   minx=p.p[3*i];
    if (p.p[3*i]>maxx)   maxx=p.p[3*i];
    if (p.p[3*i+1]<miny) miny=p.p[3*i+1];
    if (p.p[3*i+1]>maxy) maxy=p.p[3*i+1];
    if (p.p[3*i+2]<minz) minz=p.p[3*i+2];
    if (p.p[3*i+2]>maxz) maxz=p.p[3*i+2];
  }
  rangex=maxx-minx;
  rangey=maxy-miny;
  rangez=maxz-minz;

  if (rangex<rangez)
  {
    if (rangex<rangey) {
      pol2dx=2;
      x=v.z;
    } else {
      pol2dy=2;
      y=v.z;
    }
  }
  else if (rangey<rangez) {
    pol2dy=2;
    y=v.z;
  }

  //trace rays and test number of intersection
  int j;
  for (i = 0, j = p.npol-1; i < p.npol; j = i++) {
    if (((((p.p[3*i+pol2dy])<=y) && (y<(p.p[3*j+pol2dy]))) ||
         (((p.p[3*j+pol2dy])<=y) && (y<(p.p[3*i+pol2dy])))) &&
        (x < ( (p.p[3*j+pol2dx] - p.p[3*i+pol2dx]) * (y - p.p[3*i+pol2dy])
             / (p.p[3*j+pol2dy] - p.p[3*i+pol2dy]) + p.p[3*i+pol2dx]) ))
      c = !c;

    if (((fabs(p.p[3*i+pol2dy]-y)<=EPSILON) || ((fabs(p.p[3*j+pol2dy]-y)<=EPSILON))) &&
        fabs(x -((p.p[3*j+pol2dx] - p.p[3*i+pol2dx]) * (y - p.p[3*i+pol2dy])
          / (p.p[3*j+pol2dy] - p.p[3*i+pol2dy]) + p.p[3*i+pol2dx])) < EPSILON)
    {
      //the point lies on the edge
      c=-1;
      break;
    }
  }

  return c;
} /* off_pnpoly */

// off_intersectPoly ***********************************************************
//gives the intersection vertex between ray [a,b) and polygon p and its parametric value on (a b)
//based on http://geometryalgorithms.com/Archive/algorithm_0105/algorithm_0105.htm
int off_intersectPoly(intersection *inter, Coords a, Coords b, polygon p)
{
  //direction vector of [a,b]
  Coords dir = {b.x-a.x, b.y-a.y, b.z-a.z};

  //the normal vector to the polygon
  Coords normale=p.normal;
  //off_normal(&normale, p); done at the init stage

  //direction vector from a to a vertex of the polygon
  Coords w0 = {a.x-p.p[0], a.y-p.p[1], a.z-p.p[2]};

  //scalar product
  MCNUM nw0  =-scalar_prod(normale.x,normale.y,normale.z,w0.x,w0.y,w0.z);
  MCNUM ndir = scalar_prod(normale.x,normale.y,normale.z,dir.x,dir.y,dir.z);
  inter->time = inter->edge = inter->in_out=0;
  inter->v = inter->normal = coords_set(0,0,1);

  if (fabs(ndir) < EPSILON)    // ray is parallel to polygon plane
  {
    if (nw0 == 0)              // ray lies in polygon plane (infinite number of solution)
      return 0;
    else return 0;             // ray disjoint from plane (no solution)
  }

  // get intersect point of ray with polygon plane
  inter->time = nw0 / ndir;            //parametric value the point on line (a,b)

  inter->v = coords_set(a.x + inter->time * dir.x,// intersect point of ray and plane
    a.y + inter->time * dir.y,
    a.z + inter->time * dir.z);

  int res=off_pnpoly(p,inter->v);

  inter->edge=(res==-1);
  if (ndir<0)
    inter->in_out=1;  //the negative dot product means we enter the surface
  else
    inter->in_out=-1;

  inter->normal=p.normal;

  return res;         //true if the intersection point lies inside the poly
} /* off_intersectPoly */


// off_getBlocksIndex **********************************************************
/*reads the indexes at the beginning of the off file as this :
line 1  OFF
line 2  nbVertex nbFaces nbEdges
*/
long off_getBlocksIndex(char* filename, long* vtxIndex, long* vtxSize, long* faceIndex, long* polySize )
{
  if (!filename)             return(0);
  if (strlen(filename) == 0) return (0);
  if (!strcmp(filename,"NULL") || !strcmp(filename,"0"))  return(0);
  FILE* f = fopen(filename,"r");
  if (!f) {
    char mc_rt_path[CHAR_BUF_LENGTH];
    char mc_rt_dir[CHAR_BUF_LENGTH];

    if (!f)
    {
      strcpy(mc_rt_dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
      sprintf(mc_rt_path, "%s%c%s%c%s", mc_rt_dir, MC_PATHSEP_C, "data", MC_PATHSEP_C, filename);
      f = fopen(mc_rt_path, "r");
    }
    if (!f)
    {
      strcpy(mc_rt_dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
      sprintf(mc_rt_path, "%s%c%s%c%s", mc_rt_dir, MC_PATHSEP_C, "contrib", MC_PATHSEP_C, filename);
      f = fopen(mc_rt_path, "r");
    }
    if (!f)
    {
      fprintf(stderr, "Error: Could not open input file '%s' (interoff/off_getBlocksIndex)\n", filename);
      return (0);
    }
  }
  MPI_MASTER(
  printf("Loading geometry file (OFF/PLY): %s\n", filename);
  );
  
  char line[CHAR_BUF_LENGTH];
  char *ret=0;
  *vtxIndex = *vtxSize = *faceIndex = *polySize = 0;

  /* **************** start to read the file header */
  /* OFF file:
     'OFF' or '3'
   */

  ret=fgets(line,CHAR_BUF_LENGTH , f);// line 1 = "OFF"
  if (ret == NULL)
  {
    fprintf(stderr, "Error: Can not read 1st line in file %s (interoff/off_getBlocksIndex)\n", filename);
    exit(1);
  }

  if (strncmp(line,"OFF",3) && strncmp(line,"3",1) && strncmp(line,"ply",1))
  {
    fprintf(stderr, "Error: %s is probably not an OFF, NOFF or PLY file (interoff/off_getBlocksIndex).\n"
                    "       Requires first line to be 'OFF', '3' or 'ply'.\n",filename);
    fclose(f);
    return(0);
  }

  *vtxIndex+= strlen(line);
  if (!strncmp(line,"OFF",3) || !strncmp(line,"3",1)) {
    do  /* OFF file: skip # comments which may be there */
    {
      ret=fgets(line,CHAR_BUF_LENGTH , f);
      if (ret == NULL)
      {
        fprintf(stderr, "Error: Can not read line in file %s (interoff/off_getBlocksIndex)\n", filename);
        exit(1);
      }
      *vtxIndex+= strlen(line);
    } while (line[0]=='#');
    //line = nblines of vertex,faces and edges arrays
    sscanf(line,"%lu %lu",vtxSize,polySize);
  } else {
    do  /* PLY file: read all lines until find 'end_header'
           and locate 'element faces' and 'element vertex' */
    {
      ret=fgets(line,CHAR_BUF_LENGTH , f);
      if (ret == NULL)
      {
        fprintf(stderr, "Error: Can not read line in file %s (interoff/off_getBlocksIndex)\n", filename);
        exit(1);
      }
      if (!strncmp(line,"element face",12))
        sscanf(line,"element face %lu",polySize);
      else if (!strncmp(line,"element vertex",14))
        sscanf(line,"element vertex %lu",vtxSize);
      else if (!strncmp(line,"format binary",13))
        exit(fprintf(stderr,
          "Error: Can not read binary PLY file %s, only 'format ascii' (interoff/off_getBlocksIndex)\n%s\n",
          filename, line));
      *vtxIndex+= strlen(line);
    } while (strncmp(line,"end_header",10));
  }

  /* **************** read the vertices and polygons */
  *faceIndex=*vtxIndex;
  int i=0;
  for (i=0; i<*vtxSize; )
  {
    ret=fgets(line,CHAR_BUF_LENGTH,f);
    if (ret == NULL)
    {
      fprintf(stderr,
        "Error: Can not read vertex %i in file %s (interoff/off_getBlocksIndex)\n",
        i, filename);
      exit(1);
    }
    *faceIndex+=strlen(line);
    if (line[0]!='#' && strncmp(line,"comment",7)) i++;   /* do not count comments */
  }

  fclose(f);
  return(*vtxIndex);
} /* off_getBlocksIndex */

// off_init_planes *************************************************************
//gives the equations of 2 perpandicular planes of [ab]
void off_init_planes(Coords a, Coords b,
  MCNUM* A1, MCNUM* C1, MCNUM* D1, MCNUM *A2, MCNUM* B2, MCNUM* C2, MCNUM* D2)
{
  //direction vector of [a b]
  Coords dir={b.x-a.x, b.y-a.y, b.z-a.z};

  //the plane parallel to the 'y' is computed with the normal vector of the projection of [ab] on plane 'xz'
  *A1= dir.z;
  *C1=-dir.x;
  if(*A1!=0 || *C1!=0)
    *D1=-(a.x)*(*A1)-(a.z)*(*C1);
  else
  {
    //the plane does not support the vector, take the one parallel to 'z''
    *A1=1;
    //B1=dir.x=0
    *D1=-(a.x);
  }
  //the plane parallel to the 'x' is computed with the normal vector of the projection of [ab] on plane 'yz'
  *B2= dir.z;
  *C2=-dir.y;
  *A2= 0;
  if (*B2==0 && *C2==0)
  {
    //the plane does not support the vector, take the one parallel to 'z'
    *B2=1;
    //B1=dir.x=0
    *D2=-(a.y);
  }
  else {
    if (dir.z==0)
    {
      //the planes are the same, take the one parallel to 'z'
      *A2= dir.y;
      *B2=-dir.x;
      *D2=-(a.x)*(*A2)-(a.y)*(*B2);
    }
    else
      *D2=-(a.y)**B2-(a.z)**C2;
  }
} /* off_init_planes */

// off_clip_3D_mod *************************************************************
int off_clip_3D_mod(intersection* t, Coords a, Coords b,
  Coords* vtxArray, unsigned long vtxSize, unsigned long* faceArray,
  unsigned long faceSize, Coords* normalArray)
{
  MCNUM A1=0, C1=0, D1=0, A2=0, B2=0, C2=0, D2=0;      //perpendicular plane equations to [a,b]
  off_init_planes(a, b, &A1, &C1, &D1, &A2, &B2, &C2, &D2);

  int t_size=0;
  //unsigned long vtxSize=vtxTable.rows, faceSize=faceTable.columns;  //Size of the corresponding tables
  char sg[vtxSize];  //array telling if vertex is left or right of the plane
  MCNUM popol[3*CHAR_BUF_LENGTH];
  unsigned long i=0,indPoly=0;
  for (i=0; i < vtxSize; ++i)
  {
    sg[i]=off_sign(off_F(vtxArray[i].x,vtxArray[i].y,vtxArray[i].z,A1,0,C1,D1));
  }

  //exploring the polygons :
  i=indPoly=0;
  while (i<faceSize)
  {
    polygon pol;
    pol.npol  = faceArray[i];                //nb vertex of polygon
    pol.p     = popol;
    pol.normal= coords_set(0,0,1);
    unsigned long indVertP1=faceArray[++i];  //polygon's first vertex index in vtxTable
    int j=1;
    while (j<pol.npol)
    {
      //polygon's j-th vertex index in vtxTable
      if (sg[indVertP1]!=sg[faceArray[i+j]]) //if the plane intersect the polygon
        break;

      ++j;
    }

    if (j<pol.npol)          //ok, let's test with the second plane
    {
      char sg1=off_sign(off_F(vtxArray[indVertP1].x,vtxArray[indVertP1].y,vtxArray[indVertP1].z,A2,B2,C2,D2));//tells if vertex is left or right of the plane

      j=1;
      while (j<pol.npol)
      {
        //unsigned long indVertPi=faceArray[i+j];  //polyg's j-th vertex index in vtxTable
        Coords vertPi=vtxArray[faceArray[i+j]];
        if (sg1!=off_sign(off_F(vertPi.x,vertPi.y,vertPi.z,A2,B2,C2,D2)))//if the plane intersect the polygon
          break;
        ++j;
      }
      if (j<pol.npol)
      {
        if (t_size>CHAR_BUF_LENGTH)
        {
          fprintf(stderr, "Warning: number of intersection exceeded (%d) (interoff-lib/off_clip_3D_mod)\n", CHAR_BUF_LENGTH);
            return (t_size);
        }
        //both planes intersect the polygon, let's find the intersection point
        //our polygon :
        int k;
        for (k=0; k<pol.npol; ++k)
        {
          Coords vertPk=vtxArray[faceArray[i+k]];
          pol.p[3*k]  =vertPk.x;
          pol.p[3*k+1]=vertPk.y;
          pol.p[3*k+2]=vertPk.z;
        }
        pol.normal=normalArray[indPoly];
        intersection x;
        if (off_intersectPoly(&x, a, b, pol))
        {
          t[t_size++]=x;
        }
      } /* if (j<pol.npol) */
    } /* if (j<pol.npol) */
    i += pol.npol;
    indPoly++;
  } /* while i<faceSize */
  return t_size;
} /* off_clip_3D_mod */


// off_compare *****************************************************************
int off_compare (void const *a, void const *b)
{
   intersection const *pa = a;
   intersection const *pb = b;

   return off_sign(pa->time - pb->time);
} /* off_compare */

// off_cleanDouble *************************************************************
//given an array of intersections throw those which appear several times
//returns 1 if there is a possibility of error
int off_cleanDouble(intersection* t, int* t_size)
{
  int i=1;
  intersection prev=t[0];
  while (i<*t_size)
  {
    int j=i;
    //for each intersection with the same time
    while (j<*t_size && fabs(prev.time-t[j].time)<EPSILON)
    {
      //if the intersection is the exact same erase it
      if (prev.in_out==t[j].in_out)
      {
        int k;
        for (k=j+1; k<*t_size; ++k)
        {
          t[k-1]=t[k];
        }
        *t_size-=1;
      }
      else
        ++j;
    }
    prev=t[i];
    ++i;

  }
  return 1;
} /* off_cleanDouble */

// off_cleanInOut **************************************************************
//given an array of intesections throw those which enter and exit in the same time
//Meaning the ray passes very close to the volume
//returns 1 if there is a possibility of error
int off_cleanInOut(intersection* t, int* t_size)
{
  int i=1;
  intersection prev=t[0];
  while (i<*t_size)
  {
    //if two intersection have the same time but one enters and the other exits erase both
    //(such intersections must be adjacent in the array : run off_cleanDouble before)
    if (fabs(prev.time-t[i].time)<EPSILON && prev.in_out!=t[i].in_out)
    {
      int j=0;
      for (j=i+1; j<*t_size; ++j)
      {
        t[j-2]=t[j];
      }
      *t_size-=2;
      prev=t[i-1];
    }
    else
    {
      prev=t[i];
      ++i;
    }
  }
  return (*t_size);
} /* off_cleanInOut */

/* PUBLIC functions ******************************************************** */

/*******************************************************************************
* long off_init(  char *offfile, double xwidth, double yheight, double zdepth, off_struct* data)
* ACTION: read an OFF file, optionally center object and rescale, initialize OFF data structure
* INPUT: 'offfile' OFF file to read
*        'xwidth,yheight,zdepth' if given as non-zero, apply bounding box.
*           Specifying only one of these will also use the same ratio on all axes
*        'notcenter' center the object to the (0,0,0) position in local frame when set to zero
* RETURN: number of polyhedra and 'data' OFF structure
*******************************************************************************/
long off_init(  char *offfile, double xwidth, double yheight, double zdepth,
                int notcenter, off_struct* data)
{
  //data to be initialized
  long faceSize=0, vtxSize =0, polySize=0, vtxIndex=0, faceIndex=0;
  Coords* vtxArray        =NULL;
  Coords* normalArray     =NULL;
  unsigned long* faceArray=NULL;
  t_Table vtxTable, faceTable;

  double minx=FLT_MAX,maxx=-FLT_MAX,miny=FLT_MAX,maxy=-FLT_MAX,minz=FLT_MAX,maxz=-FLT_MAX;

  // get the indexes
  if (!data || off_getBlocksIndex(offfile,&vtxIndex,&vtxSize,&faceIndex, &polySize) <=0)
    return(0);

  //read vertex table = [x y z | x y z | ...]
  Table_Read_Offset(&vtxTable, offfile, 0, &vtxIndex, vtxSize);

  //read face table = [nbvertex v1 v2 vn nbvertex v1 v2 vn ...]
  Table_Read_Offset(&faceTable, offfile, 0, &faceIndex, 0);

  //initialize Arrays
  faceSize   = faceTable.columns;
  MPI_MASTER(
  printf("  Number of polygons: %ld\n", polySize);
  printf("  Number of vertices: %ld\n", vtxSize);
  );

  vtxArray   = malloc(vtxSize*sizeof(Coords));
  normalArray= malloc(polySize*sizeof(Coords));
  if (!vtxArray || !normalArray) return(0);

  long i=0,j=0;
  for(i=0; i<vtxSize; ++i)
  {
    vtxArray[i].x=Table_Index(vtxTable, i,0);
    vtxArray[i].y=Table_Index(vtxTable, i,1);
    vtxArray[i].z=Table_Index(vtxTable, i,2);

    //bounding box
    if (vtxArray[i].x<minx) minx=vtxArray[i].x;
    if (vtxArray[i].x>maxx) maxx=vtxArray[i].x;
    if (vtxArray[i].y<miny) miny=vtxArray[i].y;
    if (vtxArray[i].y>maxy) maxy=vtxArray[i].y;
    if (vtxArray[i].z<minz) minz=vtxArray[i].z;
    if (vtxArray[i].z>maxz) maxz=vtxArray[i].z;
  }

  //resizing and repositioning params
  double centerx=0, centery=0, centerz=0;
  if (!notcenter) {
    centerx=(minx+maxx)*0.5;
    centery=(miny+maxy)*0.5;
    centerz=(minz+maxz)*0.5;
  }

  double rangex=-minx+maxx,
         rangey=-miny+maxy,
         rangez=-minz+maxz;

  double ratiox=1,ratioy=1,ratioz=1;

  if (xwidth && rangex)
  {
    ratiox=xwidth/rangex;
    ratioy=ratiox;
    ratioz=ratiox;
  }

  if (yheight && rangey)
  {
    ratioy=yheight/rangey;
    if(!xwidth)  ratiox=ratioy;
    ratioz=ratioy;
  }

  if (zdepth && rangez)
  {
    ratioz=zdepth/rangez;
    if(!xwidth)  ratiox=ratioz;
    if(!yheight) ratioy=ratioz;
  }

  rangex *= ratiox;
  rangey *= ratioy;
  rangez *= ratioz;

  //center and resize the object
  for (i=0; i<vtxSize; ++i)
  {
    vtxArray[i].x=(vtxArray[i].x-centerx)*ratiox+(!notcenter ? 0 : centerx);
    vtxArray[i].y=(vtxArray[i].y-centery)*ratioy+(!notcenter ? 0 : centery);
    vtxArray[i].z=(vtxArray[i].z-centerz)*ratioz+(!notcenter ? 0 : centerz);
  }

  //table_read create a table on one line if the number of columns is not constant, so there are 2 cases :
  if (faceTable.rows==1)
  {
    //copy the table in a 1-row array
    faceArray=malloc(faceSize*sizeof(unsigned long));
    if (!faceArray) return(0);
    for (i=0; i<faceSize; ++i)
    {
      faceArray[i]=Table_Index(faceTable, 0, i);
    }
  }
  else
  {
    //read each row of the table and concatenate in a 1-row array
    faceArray=malloc(polySize*(faceSize)*sizeof(unsigned long));
    if (!faceArray) return(0);
    for(i=0; i<polySize; ++i)
    {
      for(j=0; j<faceSize; ++j)
        faceArray[i*(faceSize)+j]=Table_Index(faceTable, i, j);
    }
    faceSize*=polySize;
  }

  //precomputes normals
  long indNormal=0;//index in polyArray
  i=0;    //index in faceArray
  while (i<faceSize)
  {
    int    nbVertex=faceArray[i];//nb of vertices of this polygon
    double vertices[3*nbVertex];
    int j;

    for (j=0; j<nbVertex; ++j)
    {
      unsigned long indVertPj=faceArray[i+j+1];
      vertices[3*j]  =vtxArray[indVertPj].x;
      vertices[3*j+1]=vtxArray[indVertPj].y;
      vertices[3*j+2]=vtxArray[indVertPj].z;
    }

    polygon p;
    p.p   =vertices;
    p.npol=nbVertex;
    off_normal(&(p.normal),p);

    normalArray[indNormal]=p.normal;

    i += nbVertex+1;
    indNormal++;

  }
  MPI_MASTER(
  if (ratiox!=ratioy || ratiox!=ratioz || ratioy!=ratioz)
    printf("Warning: Aspect ratio of the sample was modified.\n"
           "         If you want to keep the original proportions, specifiy only one of the dimensions.\n");
  
  printf("  Bounding box dimensions:\n");
  printf("    Length=%f (%.3f%%)\n", rangex, ratiox*100);
  printf("    Width= %f (%.3f%%)\n", rangey, ratioy*100);
  printf("    Depth= %f (%.3f%%)\n", rangez, ratioz*100);
  );

  data->vtxArray   = vtxArray;
  data->normalArray= normalArray;
  data->faceArray  = faceArray;
  data->vtxSize    = vtxSize;
  data->polySize   = polySize;
  data->faceSize   = faceSize;
  return(polySize);
} /* off_init */

/*******************************************************************************
* int off_intersect(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double vx, double vy, double vz,
     off_struct data )
* ACTION: computes intersection of neutron trajectory with an object.
* INPUT:  x,y,z and vx,vy,vz are the position and velocity of the neutron
*         data points to the OFF data structure
* RETURN: the number of polyhedra which trajectory intersects
*         t0 and t3 are the smallest incoming and outgoing intersection times
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
int off_intersect(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x,  double y,  double z,
     double vx, double vy, double vz,
     off_struct data )
{
    intersection t[CHAR_BUF_LENGTH];
    Coords A={x, y, z};
    Coords B={x+vx, y+vy, z+vz};
    int t_size=off_clip_3D_mod(t, A, B,
      data.vtxArray, data.vtxSize, data.faceArray, data.faceSize, data.normalArray );
    qsort(t, t_size, sizeof(intersection),  off_compare);
    off_cleanDouble(t, &t_size);
    off_cleanInOut(t,  &t_size);

    if(t_size>0)
    {
      int i=0;
      if (t0) *t0 = t[0].time;
      if (n0) *n0 = t[0].normal;;
      if(t_size>1) {
        for (i=1; i < t_size; i++)
          if (t[i].time > 0 && t[i].time > t[0].time) break;
        if (t3) *t3 = t[i].time;
        if (n3) *n3 = t[i].normal;
      }
      return t_size;
    }
    return 0;
} /* off_intersect */


/*****************************************************************************
* int off_intersectx(double* l0, double* l3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double kx, double ky, double kz,
     off_struct data )
* ACTION: computes intersection of an xray trajectory with an object.
* INPUT:  x,y,z and kx,ky,kz, are spatial coordinates and wavevector of the x-ray
*         respectively. data points to the OFF data structure.
* RETURN: the number of polyhedra the trajectory intersects
*         l0 and l3 are the smallest incoming and outgoing intersection lengths
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
int off_x_intersect(double *l0,double *l3,
     Coords *n0, Coords *n3,
     double x,  double y,  double z,
     double kx, double ky, double kz,
     off_struct data )
{
  /*This function simply reformats and calls off_intersect (as for neutrons)
   *by normalizing the wavevector - this will yield the intersection lengths
   *in m*/
  double jx,jy,jz,invk;
  int n;
  invk=1/sqrt(scalar_prod(kx,ky,kz,kx,ky,kz));
  jx=kx*invk;jy=ky*invk;jz=kz*invk;
  n=off_intersect(l0,l3,n0,n3,x,y,z,jx,jy,jz,data);
  return n;
}


/*******************************************************************************
* void off_display(off_struct data)
* ACTION: display up to N_VERTEX_DISPLAYED polygons from the object
*******************************************************************************/
void off_display(off_struct data)
{
  unsigned int i;
  double ratio=(double)(N_VERTEX_DISPLAYED)/(double)data.faceSize;
  for (i=0; i<data.faceSize-1; i++) {
    int j;
    int nbVertex = data.faceArray[i];
    double x0,y0,z0;
    x0 = data.vtxArray[data.faceArray[i+1]].x;
    y0 = data.vtxArray[data.faceArray[i+1]].y;
    z0 = data.vtxArray[data.faceArray[i+1]].z;
    if (ratio > 1 || rand01() < ratio) {
      double x1=x0,y1=y0,z1=z0;
      for (j=2; j<=nbVertex; j++) {
        double x2,y2,z2;
        x2 = data.vtxArray[data.faceArray[i+j]].x;
        y2 = data.vtxArray[data.faceArray[i+j]].y;
        z2 = data.vtxArray[data.faceArray[i+j]].z;
        mcdis_line(x1,y1,z1,x2,y2,z2);
        x1 = x2; y1 = y2; z1 = z2;
      }
      mcdis_line(x1,y1,z1,x0,y0,z0);
    }
    i += nbVertex;
  }
} /* off_display */

/* end of interoff-lib.c */

#line 10738 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

/* Instrument parameters. */
MCNUM mcipdistance;

#define mcNUMIPAR 1
int mcnumipar = 1;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "distance", &mcipdistance, instr_type_double, "1.0", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  test
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposatest coords_set(0,0,0)
#define distance mcipdistance
#undef distance
#undef mcposatest
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*7];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[7];
Coords mccomp_posr[7];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[7];
MCNUM  mcPCounter[7];
MCNUM  mcP2Counter[7];
#define mcNUMCOMP 6 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[7];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Declarations of component definition and setting parameters. */

/* Setting parameters for component 'Origin' [1]. */
char mccOrigin_profile[16384];
MCNUM mccOrigin_percent;
MCNUM mccOrigin_flag_save;
MCNUM mccOrigin_minutes;

/* Setting parameters for component 'ESS_1b_source' [2]. */
char mccESS_1b_source_filename[16384];
char mccESS_1b_source_type[16384];
MCNUM mccESS_1b_source_verbose;
MCNUM mccESS_1b_source_repeat_count;
MCNUM mccESS_1b_source_smooth;

/* Definition parameters for component 'v_mon' [3]. */
#define mccv_mon_options "v limits=[0,10000], bins=1000" /* declared as a string. May produce warnings at compile */
#define mccv_mon_filename "v_monitor" /* declared as a string. May produce warnings at compile */
#define mccv_mon_geometry 0 /* declared as a string. May produce warnings at compile */
#define mccv_mon_user1 FLT_MAX
#define mccv_mon_user2 FLT_MAX
#define mccv_mon_user3 FLT_MAX
#define mccv_mon_username1 0
#define mccv_mon_username2 0
#define mccv_mon_username3 0
/* Setting parameters for component 'v_mon' [3]. */
MCNUM mccv_mon_xwidth;
MCNUM mccv_mon_yheight;
MCNUM mccv_mon_zdepth;
MCNUM mccv_mon_xmin;
MCNUM mccv_mon_xmax;
MCNUM mccv_mon_ymin;
MCNUM mccv_mon_ymax;
MCNUM mccv_mon_zmin;
MCNUM mccv_mon_zmax;
MCNUM mccv_mon_bins;
MCNUM mccv_mon_min;
MCNUM mccv_mon_max;
MCNUM mccv_mon_restore_neutron;
MCNUM mccv_mon_radius;

/* Definition parameters for component 'image' [4]. */
#define mccimage_nx 90
#define mccimage_ny 90
#define mccimage_restore_neutron 1
/* Setting parameters for component 'image' [4]. */
char mccimage_filename[16384];
MCNUM mccimage_xmin;
MCNUM mccimage_xmax;
MCNUM mccimage_ymin;
MCNUM mccimage_ymax;
MCNUM mccimage_xwidth;
MCNUM mccimage_yheight;

/* Definition parameters for component 'lamb' [5]. */
#define mcclamb_nL 1000
/* Setting parameters for component 'lamb' [5]. */
char mcclamb_filename[16384];
MCNUM mcclamb_xmin;
MCNUM mcclamb_xmax;
MCNUM mcclamb_ymin;
MCNUM mcclamb_ymax;
MCNUM mcclamb_xwidth;
MCNUM mcclamb_yheight;
MCNUM mcclamb_Lmin;
MCNUM mcclamb_Lmax;
MCNUM mcclamb_restore_neutron;

/* User component declarations. */

/* User declarations for component 'Origin' [1]. */
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define profile mccOrigin_profile
#define percent mccOrigin_percent
#define flag_save mccOrigin_flag_save
#define minutes mccOrigin_minutes
#line 47 "/usr/local/lib/mcstas-2.0/misc/Progress_bar.comp"
#ifndef PROGRESS_BAR
#define PROGRESS_BAR
#else
#error Only one Progress_bar component may be used in an instrument definition.
#endif

  double IntermediateCnts=0;
  time_t StartTime       =0;
  time_t EndTime         =0;
  time_t CurrentTime     =0;
#line 10869 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
#undef minutes
#undef flag_save
#undef percent
#undef profile
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'ESS_1b_source' [2]. */
#define mccompcurname  ESS_1b_source
#define mccompcurtype  Virtual_input
#define mccompcurindex 2
#define read_block mccESS_1b_source_read_block
#define pos mccESS_1b_source_pos
#define nrows mccESS_1b_source_nrows
#define Offset mccESS_1b_source_Offset
#define rTable mccESS_1b_source_rTable
#define repeat_number mccESS_1b_source_repeat_number
#define filename_ncount mccESS_1b_source_filename_ncount
#define mean_vx mccESS_1b_source_mean_vx
#define mean_vy mccESS_1b_source_mean_vy
#define mean_vz mccESS_1b_source_mean_vz
#define mean_dx mccESS_1b_source_mean_dx
#define mean_dy mccESS_1b_source_mean_dy
#define mean_dz mccESS_1b_source_mean_dz
#define n_neutrons mccESS_1b_source_n_neutrons
#define min_x mccESS_1b_source_min_x
#define min_y mccESS_1b_source_min_y
#define min_z mccESS_1b_source_min_z
#define max_x mccESS_1b_source_max_x
#define max_y mccESS_1b_source_max_y
#define max_z mccESS_1b_source_max_z
#define min_vx mccESS_1b_source_min_vx
#define min_vy mccESS_1b_source_min_vy
#define min_vz mccESS_1b_source_min_vz
#define max_vx mccESS_1b_source_max_vx
#define max_vy mccESS_1b_source_max_vy
#define max_vz mccESS_1b_source_max_vz
#define first_block mccESS_1b_source_first_block
#define mean_x mccESS_1b_source_mean_x
#define mean_y mccESS_1b_source_mean_y
#define mean_z mccESS_1b_source_mean_z
#define end_reading mccESS_1b_source_end_reading
#define n_count_extrapolated mccESS_1b_source_n_count_extrapolated
#define repeat_cnt mccESS_1b_source_repeat_cnt
#define filename mccESS_1b_source_filename
#define type mccESS_1b_source_type
#define verbose mccESS_1b_source_verbose
#define repeat_count mccESS_1b_source_repeat_count
#define smooth mccESS_1b_source_smooth
#line 110 "/usr/local/lib/mcstas-2.0/sources/Virtual_input.comp"
  int     repeat_number=1; /* Neutron repeat of the filename */
  int     repeat_cnt;     /* Repeat count, MPI taken into account */
  long    pos=0;        /* current pos in block */
  long    nrows=0;      /* total nrows in block */
  long    Offset=0;     /* offset in filename */
  double  filename_ncount=0;  /* total number of neutrons in filename */
  double  n_neutrons=0;
  char    read_block=1; /* flag to start by reading block */
  char    first_block=1;
  char    end_reading=0;
  t_Table rTable;

  /* statistics on first block */
  double mean_x =0, mean_y =0, mean_z =0;
  double mean_vx=0, mean_vy=0, mean_vz=0;
  double mean_dx=0, mean_dy=0, mean_dz=0;
  double min_x = FLT_MAX, min_y = FLT_MAX, min_z = FLT_MAX;
  double max_x =-FLT_MAX, max_y =-FLT_MAX, max_z =-FLT_MAX;
  double min_vx= FLT_MAX, min_vy= FLT_MAX, min_vz= FLT_MAX;
  double max_vx=-FLT_MAX, max_vy=-FLT_MAX, max_vz=-FLT_MAX;
  double n_count_extrapolated=0;
#line 10945 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
#undef smooth
#undef repeat_count
#undef verbose
#undef type
#undef filename
#undef repeat_cnt
#undef n_count_extrapolated
#undef end_reading
#undef mean_z
#undef mean_y
#undef mean_x
#undef first_block
#undef max_vz
#undef max_vy
#undef max_vx
#undef min_vz
#undef min_vy
#undef min_vx
#undef max_z
#undef max_y
#undef max_x
#undef min_z
#undef min_y
#undef min_x
#undef n_neutrons
#undef mean_dz
#undef mean_dy
#undef mean_dx
#undef mean_vz
#undef mean_vy
#undef mean_vx
#undef filename_ncount
#undef repeat_number
#undef rTable
#undef Offset
#undef nrows
#undef pos
#undef read_block
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'v_mon' [3]. */
#define mccompcurname  v_mon
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define options mccv_mon_options
#define filename mccv_mon_filename
#define geometry mccv_mon_geometry
#define user1 mccv_mon_user1
#define user2 mccv_mon_user2
#define user3 mccv_mon_user3
#define username1 mccv_mon_username1
#define username2 mccv_mon_username2
#define username3 mccv_mon_username3
#define DEFS mccv_mon_DEFS
#define Vars mccv_mon_Vars
#define detector mccv_mon_detector
#define xwidth mccv_mon_xwidth
#define yheight mccv_mon_yheight
#define zdepth mccv_mon_zdepth
#define xmin mccv_mon_xmin
#define xmax mccv_mon_xmax
#define ymin mccv_mon_ymin
#define ymax mccv_mon_ymax
#define zmin mccv_mon_zmin
#define zmax mccv_mon_zmax
#define bins mccv_mon_bins
#define min mccv_mon_min
#define max mccv_mon_max
#define restore_neutron mccv_mon_restore_neutron
#define radius mccv_mon_radius
#line 232 "/usr/local/lib/mcstas-2.0/monitors/Monitor_nD.comp"
  MonitornD_Defines_type DEFS;
  MonitornD_Variables_type Vars;
  MCDETECTOR detector;
  off_struct offdata;
#line 11023 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
#undef radius
#undef restore_neutron
#undef max
#undef min
#undef bins
#undef zmax
#undef zmin
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef zdepth
#undef yheight
#undef xwidth
#undef detector
#undef Vars
#undef DEFS
#undef username3
#undef username2
#undef username1
#undef user3
#undef user2
#undef user1
#undef geometry
#undef filename
#undef options
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'image' [4]. */
#define mccompcurname  image
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
#define nx mccimage_nx
#define ny mccimage_ny
#define restore_neutron mccimage_restore_neutron
#define PSD_N mccimage_PSD_N
#define PSD_p mccimage_PSD_p
#define PSD_p2 mccimage_PSD_p2
#define filename mccimage_filename
#define xmin mccimage_xmin
#define xmax mccimage_xmax
#define ymin mccimage_ymin
#define ymax mccimage_ymax
#define xwidth mccimage_xwidth
#define yheight mccimage_yheight
#line 58 "/usr/local/lib/mcstas-2.0/monitors/PSD_monitor.comp"
    double PSD_N[nx][ny];
    double PSD_p[nx][ny];
    double PSD_p2[nx][ny];
#line 11075 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_neutron
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'lamb' [5]. */
#define mccompcurname  lamb
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mcclamb_nL
#define L_N mcclamb_L_N
#define L_p mcclamb_L_p
#define L_p2 mcclamb_L_p2
#define filename mcclamb_filename
#define xmin mcclamb_xmin
#define xmax mcclamb_xmax
#define ymin mcclamb_ymin
#define ymax mcclamb_ymax
#define xwidth mcclamb_xwidth
#define yheight mcclamb_yheight
#define Lmin mcclamb_Lmin
#define Lmax mcclamb_Lmax
#define restore_neutron mcclamb_restore_neutron
#line 58 "/usr/local/lib/mcstas-2.0/monitors/L_monitor.comp"
    double L_N[nL];
    double L_p[nL], L_p2[nL];
#line 11114 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
#undef restore_neutron
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

Coords mcposaOrigin, mcposrOrigin;
Rotation mcrotaOrigin, mcrotrOrigin;
Coords mcposaESS_1b_source, mcposrESS_1b_source;
Rotation mcrotaESS_1b_source, mcrotrESS_1b_source;
Coords mcposav_mon, mcposrv_mon;
Rotation mcrotav_mon, mcrotrv_mon;
Coords mcposaimage, mcposrimage;
Rotation mcrotaimage, mcrotrimage;
Coords mcposalamb, mcposrlamb;
Rotation mcrotalamb, mcrotrlamb;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  test
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposatest coords_set(0,0,0)
#define distance mcipdistance
#undef distance
#undef mcposatest
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname
  /* Computation of coordinate transformations. */
  {
    Coords mctc1, mctc2;
    Rotation mctr1;

    mcDEBUG_INSTR()
  /* Component initializations. */
    /* Component Origin. */
    SIG_MESSAGE("Origin (Init:Place/Rotate)");
    rot_set_rotation(mcrotaOrigin,
      (0.0)*DEG2RAD,
      (0.0)*DEG2RAD,
      (0.0)*DEG2RAD);
#line 11172 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
    rot_copy(mcrotrOrigin, mcrotaOrigin);
    mcposaOrigin = coords_set(
#line 58 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0,
#line 58 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0,
#line 58 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0);
#line 11181 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
    mctc1 = coords_neg(mcposaOrigin);
    mcposrOrigin = rot_apply(mcrotaOrigin, mctc1);
    mcDEBUG_COMPONENT("Origin", mcposaOrigin, mcrotaOrigin)
    mccomp_posa[1] = mcposaOrigin;
    mccomp_posr[1] = mcposrOrigin;
    mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
    mcAbsorbProp[1]= 0;
  /* Setting parameters for component Origin. */
  SIG_MESSAGE("Origin (Init:SetPar)");
#line 42 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  if(0) strncpy(mccOrigin_profile,0, 16384); else mccOrigin_profile[0]='\0';
#line 42 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccOrigin_percent = 10;
#line 42 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccOrigin_flag_save = 0;
#line 42 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccOrigin_minutes = 0;
#line 11199 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

    /* Component ESS_1b_source. */
    SIG_MESSAGE("ESS_1b_source (Init:Place/Rotate)");
    rot_set_rotation(mctr1,
#line 63 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD,
#line 63 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD,
#line 63 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD);
#line 11210 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
    rot_mul(mctr1, mcrotaOrigin, mcrotaESS_1b_source);
    rot_transpose(mcrotaOrigin, mctr1);
    rot_mul(mcrotaESS_1b_source, mctr1, mcrotrESS_1b_source);
    mctc1 = coords_set(
#line 62 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0,
#line 62 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0,
#line 62 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0);
#line 11221 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
    rot_transpose(mcrotaOrigin, mctr1);
    mctc2 = rot_apply(mctr1, mctc1);
    mcposaESS_1b_source = coords_add(mcposaOrigin, mctc2);
    mctc1 = coords_sub(mcposaOrigin, mcposaESS_1b_source);
    mcposrESS_1b_source = rot_apply(mcrotaESS_1b_source, mctc1);
    mcDEBUG_COMPONENT("ESS_1b_source", mcposaESS_1b_source, mcrotaESS_1b_source)
    mccomp_posa[2] = mcposaESS_1b_source;
    mccomp_posr[2] = mcposrESS_1b_source;
    mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
    mcAbsorbProp[2]= 0;
  /* Setting parameters for component ESS_1b_source. */
  SIG_MESSAGE("ESS_1b_source (Init:SetPar)");
#line 61 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  if("1b_source_ess") strncpy(mccESS_1b_source_filename,"1b_source_ess", 16384); else mccESS_1b_source_filename[0]='\0';
#line 61 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  if("float") strncpy(mccESS_1b_source_type,"float", 16384); else mccESS_1b_source_type[0]='\0';
#line 61 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccESS_1b_source_verbose = 1;
#line 72 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccESS_1b_source_repeat_count = 1;
#line 61 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccESS_1b_source_smooth = 0;
#line 11244 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

    /* Component v_mon. */
    SIG_MESSAGE("v_mon (Init:Place/Rotate)");
    rot_set_rotation(mctr1,
#line 69 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD,
#line 69 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD,
#line 69 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD);
#line 11255 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
    rot_mul(mctr1, mcrotaOrigin, mcrotav_mon);
    rot_transpose(mcrotaESS_1b_source, mctr1);
    rot_mul(mcrotav_mon, mctr1, mcrotrv_mon);
    mctc1 = coords_set(
#line 68 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0,
#line 68 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0,
#line 68 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      mcipdistance);
#line 11266 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
    rot_transpose(mcrotaOrigin, mctr1);
    mctc2 = rot_apply(mctr1, mctc1);
    mcposav_mon = coords_add(mcposaOrigin, mctc2);
    mctc1 = coords_sub(mcposaESS_1b_source, mcposav_mon);
    mcposrv_mon = rot_apply(mcrotav_mon, mctc1);
    mcDEBUG_COMPONENT("v_mon", mcposav_mon, mcrotav_mon)
    mccomp_posa[3] = mcposav_mon;
    mccomp_posr[3] = mcposrv_mon;
    mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
    mcAbsorbProp[3]= 0;
  /* Setting parameters for component v_mon. */
  SIG_MESSAGE("v_mon (Init:SetPar)");
#line 67 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_xwidth = 5;
#line 67 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_yheight = 5;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_zdepth = 0;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_xmin = 0;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_xmax = 0;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_ymin = 0;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_ymax = 0;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_zmin = 0;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_zmax = 0;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_bins = 0;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_min = -1e40;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_max = 1e40;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_restore_neutron = 0;
#line 217 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccv_mon_radius = 0;
#line 11307 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

    /* Component image. */
    SIG_MESSAGE("image (Init:Place/Rotate)");
    rot_set_rotation(mctr1,
#line 74 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD,
#line 74 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD,
#line 74 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD);
#line 11318 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
    rot_mul(mctr1, mcrotaOrigin, mcrotaimage);
    rot_transpose(mcrotav_mon, mctr1);
    rot_mul(mcrotaimage, mctr1, mcrotrimage);
    mctc1 = coords_set(
#line 73 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0,
#line 73 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0,
#line 73 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      mcipdistance);
#line 11329 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
    rot_transpose(mcrotaOrigin, mctr1);
    mctc2 = rot_apply(mctr1, mctc1);
    mcposaimage = coords_add(mcposaOrigin, mctc2);
    mctc1 = coords_sub(mcposav_mon, mcposaimage);
    mcposrimage = rot_apply(mcrotaimage, mctc1);
    mcDEBUG_COMPONENT("image", mcposaimage, mcrotaimage)
    mccomp_posa[4] = mcposaimage;
    mccomp_posr[4] = mcposrimage;
    mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
    mcAbsorbProp[4]= 0;
  /* Setting parameters for component image. */
  SIG_MESSAGE("image (Init:SetPar)");
#line 72 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  if("imager") strncpy(mccimage_filename,"imager", 16384); else mccimage_filename[0]='\0';
#line 52 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccimage_xmin = -0.05;
#line 52 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccimage_xmax = 0.05;
#line 52 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccimage_ymin = -0.05;
#line 52 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccimage_ymax = 0.05;
#line 72 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccimage_xwidth = 5.0;
#line 72 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mccimage_yheight = 5.0;
#line 11356 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

    /* Component lamb. */
    SIG_MESSAGE("lamb (Init:Place/Rotate)");
    rot_set_rotation(mctr1,
#line 80 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD,
#line 80 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD,
#line 80 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      (0)*DEG2RAD);
#line 11367 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
    rot_mul(mctr1, mcrotaOrigin, mcrotalamb);
    rot_transpose(mcrotaimage, mctr1);
    rot_mul(mcrotalamb, mctr1, mcrotrlamb);
    mctc1 = coords_set(
#line 79 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0,
#line 79 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      0,
#line 79 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
      mcipdistance);
#line 11378 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
    rot_transpose(mcrotaOrigin, mctr1);
    mctc2 = rot_apply(mctr1, mctc1);
    mcposalamb = coords_add(mcposaOrigin, mctc2);
    mctc1 = coords_sub(mcposaimage, mcposalamb);
    mcposrlamb = rot_apply(mcrotalamb, mctc1);
    mcDEBUG_COMPONENT("lamb", mcposalamb, mcrotalamb)
    mccomp_posa[5] = mcposalamb;
    mccomp_posr[5] = mcposrlamb;
    mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
    mcAbsorbProp[5]= 0;
  /* Setting parameters for component lamb. */
  SIG_MESSAGE("lamb (Init:SetPar)");
#line 77 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  if("lambda") strncpy(mcclamb_filename,"lambda", 16384); else mcclamb_filename[0]='\0';
#line 52 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mcclamb_xmin = -0.05;
#line 52 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mcclamb_xmax = 0.05;
#line 52 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mcclamb_ymin = -0.05;
#line 52 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mcclamb_ymax = 0.05;
#line 77 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mcclamb_xwidth = 5.0;
#line 77 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mcclamb_yheight = 5.0;
#line 78 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mcclamb_Lmin = 0.0;
#line 78 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mcclamb_Lmax = 10.0;
#line 53 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.instr"
  mcclamb_restore_neutron = 0;
#line 11411 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"

  /* Component initializations. */
  /* Initializations for component Origin. */
  SIG_MESSAGE("Origin (Init)");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define profile mccOrigin_profile
#define percent mccOrigin_percent
#define flag_save mccOrigin_flag_save
#define minutes mccOrigin_minutes
#line 60 "/usr/local/lib/mcstas-2.0/misc/Progress_bar.comp"
{
  fprintf(stdout, "[%s] Initialize\n", mcinstrument_name);
  if (percent*mcget_ncount()/100 < 1e5) {
    percent=1e5*100.0/mcget_ncount();
  }
}
#line 11433 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
#undef minutes
#undef flag_save
#undef percent
#undef profile
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component ESS_1b_source. */
  SIG_MESSAGE("ESS_1b_source (Init)");
#define mccompcurname  ESS_1b_source
#define mccompcurtype  Virtual_input
#define mccompcurindex 2
#define read_block mccESS_1b_source_read_block
#define pos mccESS_1b_source_pos
#define nrows mccESS_1b_source_nrows
#define Offset mccESS_1b_source_Offset
#define rTable mccESS_1b_source_rTable
#define repeat_number mccESS_1b_source_repeat_number
#define filename_ncount mccESS_1b_source_filename_ncount
#define mean_vx mccESS_1b_source_mean_vx
#define mean_vy mccESS_1b_source_mean_vy
#define mean_vz mccESS_1b_source_mean_vz
#define mean_dx mccESS_1b_source_mean_dx
#define mean_dy mccESS_1b_source_mean_dy
#define mean_dz mccESS_1b_source_mean_dz
#define n_neutrons mccESS_1b_source_n_neutrons
#define min_x mccESS_1b_source_min_x
#define min_y mccESS_1b_source_min_y
#define min_z mccESS_1b_source_min_z
#define max_x mccESS_1b_source_max_x
#define max_y mccESS_1b_source_max_y
#define max_z mccESS_1b_source_max_z
#define min_vx mccESS_1b_source_min_vx
#define min_vy mccESS_1b_source_min_vy
#define min_vz mccESS_1b_source_min_vz
#define max_vx mccESS_1b_source_max_vx
#define max_vy mccESS_1b_source_max_vy
#define max_vz mccESS_1b_source_max_vz
#define first_block mccESS_1b_source_first_block
#define mean_x mccESS_1b_source_mean_x
#define mean_y mccESS_1b_source_mean_y
#define mean_z mccESS_1b_source_mean_z
#define end_reading mccESS_1b_source_end_reading
#define n_count_extrapolated mccESS_1b_source_n_count_extrapolated
#define repeat_cnt mccESS_1b_source_repeat_cnt
#define filename mccESS_1b_source_filename
#define type mccESS_1b_source_type
#define verbose mccESS_1b_source_verbose
#define repeat_count mccESS_1b_source_repeat_count
#define smooth mccESS_1b_source_smooth
#line 134 "/usr/local/lib/mcstas-2.0/sources/Virtual_input.comp"
{
  Table_Init(&rTable, 0, 0);

  if (!filename || !repeat_count)
  {
    fprintf(stderr,"Virtual_input: %s: please give me a filename name (filename) to read (repeat_count>0).\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (filename && strlen(filename) && repeat_count>0) {
    if (type && strstr(type, "Vitess"))
    { fprintf(stderr, "Virtual_input: %s: Vitess files may be read using the Vitess_input component\n", NAME_CURRENT_COMP); exit(-1); }

    if (verbose)
      printf("Virtual_input: %s: Reading neutron events from filename '%s'. Repeat %g time(s)\n", NAME_CURRENT_COMP, filename, repeat_count);

#if defined (USE_MPI)
    if (!smooth && mpi_node_count > 1) {
      if (verbose)
      printf("Virtual_input: %s: smoothing (smooth=1) is recommended when running MPI execution\n", NAME_CURRENT_COMP);
    }
#endif

    double min_dv=fabs(max_vx-min_vx);
    if (min_dv > fabs(max_vy-min_vy)) min_dv = fabs(max_vy-min_vy);
    if (min_dv > fabs(max_vz-min_vz)) min_dv = fabs(max_vz-min_vz);
    min_vx = min_dv;

    if (verbose && smooth)
        printf("* Beam will be smoothed\n");
  } else if (!filename)
    exit(fprintf(stderr,"Virtual_input: %s: please give me a filename name (filename).\n", NAME_CURRENT_COMP));
  
  repeat_cnt = repeat_count;
#if defined (USE_MPI)
  repeat_cnt = ceil(1.0*repeat_cnt/mpi_node_count);
#endif
}
#line 11527 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
#undef smooth
#undef repeat_count
#undef verbose
#undef type
#undef filename
#undef repeat_cnt
#undef n_count_extrapolated
#undef end_reading
#undef mean_z
#undef mean_y
#undef mean_x
#undef first_block
#undef max_vz
#undef max_vy
#undef max_vx
#undef min_vz
#undef min_vy
#undef min_vx
#undef max_z
#undef max_y
#undef max_x
#undef min_z
#undef min_y
#undef min_x
#undef n_neutrons
#undef mean_dz
#undef mean_dy
#undef mean_dx
#undef mean_vz
#undef mean_vy
#undef mean_vx
#undef filename_ncount
#undef repeat_number
#undef rTable
#undef Offset
#undef nrows
#undef pos
#undef read_block
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component v_mon. */
  SIG_MESSAGE("v_mon (Init)");
#define mccompcurname  v_mon
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define options mccv_mon_options
#define filename mccv_mon_filename
#define geometry mccv_mon_geometry
#define user1 mccv_mon_user1
#define user2 mccv_mon_user2
#define user3 mccv_mon_user3
#define username1 mccv_mon_username1
#define username2 mccv_mon_username2
#define username3 mccv_mon_username3
#define DEFS mccv_mon_DEFS
#define Vars mccv_mon_Vars
#define detector mccv_mon_detector
#define xwidth mccv_mon_xwidth
#define yheight mccv_mon_yheight
#define zdepth mccv_mon_zdepth
#define xmin mccv_mon_xmin
#define xmax mccv_mon_xmax
#define ymin mccv_mon_ymin
#define ymax mccv_mon_ymax
#define zmin mccv_mon_zmin
#define zmax mccv_mon_zmax
#define bins mccv_mon_bins
#define min mccv_mon_min
#define max mccv_mon_max
#define restore_neutron mccv_mon_restore_neutron
#define radius mccv_mon_radius
#line 239 "/usr/local/lib/mcstas-2.0/monitors/Monitor_nD.comp"
{
  char tmp[CHAR_BUF_LENGTH];
  strcpy(Vars.compcurname, NAME_CURRENT_COMP);
  if (options != NULL)
    strncpy(Vars.option, options, CHAR_BUF_LENGTH);
  else {
    strcpy(Vars.option, "x y");
    printf("Monitor_nD: %s has no option specified. Setting to PSD ('x y') monitor.\n", NAME_CURRENT_COMP);
  }
  Vars.compcurpos = POS_A_CURRENT_COMP;

  if (strstr(Vars.option, "source"))
    strcat(Vars.option, " list, x y z vx vy vz t sx sy sz ");

  if (bins) { sprintf(tmp, " all bins=%ld ", (long)bins); strcat(Vars.option, tmp); }
  if (min > -FLT_MAX && max < FLT_MAX) { sprintf(tmp, " all limits=[%g %g]", min, max); strcat(Vars.option, tmp); }
  else if (min > -FLT_MAX) { sprintf(tmp, " all min=%g", min); strcat(Vars.option, tmp); }
  else if (max <  FLT_MAX) { sprintf(tmp, " all max=%g", max); strcat(Vars.option, tmp); }

  strncpy(Vars.UserName1, username1 && strlen(username1) ? username1 : "", 32);
  strncpy(Vars.UserName2, username2 && strlen(username2) ? username2 : "", 32);
  strncpy(Vars.UserName3, username3 && strlen(username3) ? username3 : "", 32);
  if (radius) { 
    xwidth = zdepth = 2*radius;
    if (yheight && !strstr(Vars.option, "cylinder") && !strstr(Vars.option, "banana")) 
      strcat(Vars.option, " banana");
    else if (!yheight && !strstr(Vars.option ,"sphere")) {
      strcat(Vars.option, " sphere");
      yheight=2*radius;
    }
  }
  
  if (geometry && strlen(geometry))
    if (!off_init(  geometry, xwidth, yheight, zdepth, 0, &offdata )) {
      printf("Monitor_nD: %s could not initiate the OFF geometry. \n"
             "            Defaulting to normal Monitor dimensions.\n", NAME_CURRENT_COMP);
      strcpy(geometry, "");
    }
  
  if (!radius && !xwidth && !yheight && !zdepth && !xmin && !xmax && !ymin && !ymax && !strstr(Vars.option, "previous") && (!geometry || !strlen(geometry)))
    exit(printf("Monitor_nD: %s has no dimension specified. Aborting (radius, xwidth, yheight, zdepth, previous, geometry).\n", NAME_CURRENT_COMP));

  Monitor_nD_Init(&DEFS, &Vars, xwidth, yheight, zdepth, xmin,xmax,ymin,ymax,zmin,zmax);

  if (filename && strcmp(filename,"NULL") && strcmp(filename,"0") && strlen(filename))
    strncpy(Vars.Mon_File, filename, 128);

  /* check if user given filename with ext will be used more than once */
  if ( ((Vars.Flag_Multiple && Vars.Coord_Number > 1) || Vars.Flag_List) && strchr(Vars.Mon_File,'.') )
  { char *XY; XY = strrchr(Vars.Mon_File,'.'); *XY='_'; }
  
  if (restore_neutron) Vars.Flag_parallel=1;
  detector.m = 0;
  
#ifdef USE_MPI
  if (strstr(Vars.option, "auto") && mpi_node_count > 1)
    printf("Monitor_nD: %s is using automatic limits option 'auto' together with MPI.\n"
           "WARNING     this may create incorrect distributions (but integrated flux will be right).\n", NAME_CURRENT_COMP);
#endif
}
#line 11662 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
#undef radius
#undef restore_neutron
#undef max
#undef min
#undef bins
#undef zmax
#undef zmin
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef zdepth
#undef yheight
#undef xwidth
#undef detector
#undef Vars
#undef DEFS
#undef username3
#undef username2
#undef username1
#undef user3
#undef user2
#undef user1
#undef geometry
#undef filename
#undef options
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component image. */
  SIG_MESSAGE("image (Init)");
#define mccompcurname  image
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
#define nx mccimage_nx
#define ny mccimage_ny
#define restore_neutron mccimage_restore_neutron
#define PSD_N mccimage_PSD_N
#define PSD_p mccimage_PSD_p
#define PSD_p2 mccimage_PSD_p2
#define filename mccimage_filename
#define xmin mccimage_xmin
#define xmax mccimage_xmax
#define ymin mccimage_ymin
#define ymax mccimage_ymax
#define xwidth mccimage_xwidth
#define yheight mccimage_yheight
#line 63 "/usr/local/lib/mcstas-2.0/monitors/PSD_monitor.comp"
{
    int i,j;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("PSD_monitor: %s: Null detection area !\n"
                   "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nx; i++)
     for (j=0; j<ny; j++)
     {
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
     }
}
#line 11733 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_neutron
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component lamb. */
  SIG_MESSAGE("lamb (Init)");
#define mccompcurname  lamb
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mcclamb_nL
#define L_N mcclamb_L_N
#define L_p mcclamb_L_p
#define L_p2 mcclamb_L_p2
#define filename mcclamb_filename
#define xmin mcclamb_xmin
#define xmax mcclamb_xmax
#define ymin mcclamb_ymin
#define ymax mcclamb_ymax
#define xwidth mcclamb_xwidth
#define yheight mcclamb_yheight
#define Lmin mcclamb_Lmin
#define Lmax mcclamb_Lmax
#define restore_neutron mcclamb_restore_neutron
#line 62 "/usr/local/lib/mcstas-2.0/monitors/L_monitor.comp"
{
    int i;

    if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("L_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(0);
    }

    for (i=0; i<nL; i++)
    {
      L_N[i] = 0;
      L_p[i] = 0;
      L_p2[i] = 0;
    }
}
#line 11791 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
#undef restore_neutron
#undef Lmax
#undef Lmin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if(mcdotrace) mcdisplay();
    mcDEBUG_INSTR_END()
  }

/* NeXus support */

#ifdef USE_NEXUS

strncmp(mcnxversion,"5 zip",128);

#endif

} /* end init */

void mcraytrace(void) {
  /* Neutronics-specific defines */
#ifdef NEUTRONICS
extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;
#endif
  /* End of Neutronics-specific defines */
  /* Copy neutron state to local variables. */
  MCNUM mcnlx = mcnx;
  MCNUM mcnly = mcny;
  MCNUM mcnlz = mcnz;
  MCNUM mcnlvx = mcnvx;
  MCNUM mcnlvy = mcnvy;
  MCNUM mcnlvz = mcnvz;
  MCNUM mcnlt = mcnt;
  MCNUM mcnlsx = mcnsx;
  MCNUM mcnlsy = mcnsy;
  MCNUM mcnlsz = mcnsz;
  MCNUM mcnlp = mcnp;

  mcDEBUG_ENTER()
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define mcabsorb mcabsorbAll
  /* TRACE Component Origin [1] */
  mccoordschange(mcposrOrigin, mcrotrOrigin,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Origin (without coords transformations) */
  mcJumpTrace_Origin:
  SIG_MESSAGE("Origin (Trace)");
  mcDEBUG_COMP("Origin")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp
  STORE_NEUTRON(1,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcNCounter[1]++;
  mcPCounter[1] += p;
  mcP2Counter[1] += p*p;
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
{   /* Declarations of Origin=Progress_bar() SETTING parameters. */
char* profile = mccOrigin_profile;
MCNUM percent = mccOrigin_percent;
MCNUM flag_save = mccOrigin_flag_save;
MCNUM minutes = mccOrigin_minutes;
#line 68 "/usr/local/lib/mcstas-2.0/misc/Progress_bar.comp"
{
  double ncount;
  ncount = mcget_run_num();
  if (!StartTime) {
    time(&StartTime); /* compute starting time */
    IntermediateCnts = 1e3;
  }
  time_t NowTime;
  time(&NowTime);
  /* compute initial estimate of computation duration */
  if (!EndTime && ncount >= IntermediateCnts) {
    CurrentTime = NowTime;
    if (difftime(NowTime,StartTime) > 10) { /* wait 10 sec before writing ETA */
      EndTime = StartTime + (time_t)(difftime(NowTime,StartTime)
				     *(double)mcget_ncount()/ncount);
      IntermediateCnts = 0;
      fprintf(stdout, "\nTrace ETA ");
      if (difftime(EndTime,StartTime) < 60.0)
        fprintf(stdout, "%g [s] %% ", difftime(EndTime,StartTime));
      else if (difftime(EndTime,StartTime) > 3600.0)
        fprintf(stdout, "%g [h] %% ", difftime(EndTime,StartTime)/3600.0);
      else
        fprintf(stdout, "%g [min] %% ", difftime(EndTime,StartTime)/60.0);
    } else IntermediateCnts += 1e3;
    fflush(stdout);
  }
  
  /* display percentage when percent or minutes have reached step */
  if (EndTime &&
    (    (minutes && difftime(NowTime,CurrentTime) > minutes*60)
      || (percent && !minutes && ncount >= IntermediateCnts))   )
  {
    fprintf(stdout, "%d ", (int)(ncount*100.0/mcget_ncount())); fflush(stdout);
    CurrentTime = NowTime;
    
    IntermediateCnts = ncount + percent*mcget_ncount()/100;
    /* check that next intermediate ncount check is a multiple of the desired percentage */
    IntermediateCnts = floor(IntermediateCnts*100/percent/mcget_ncount())*percent*mcget_ncount()/100;
    /* raise flag to indicate that we did something */
    SCATTER;
    if (flag_save) mcsave(NULL);
  }
}
#line 11963 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component ESS_1b_source [2] */
  mccoordschange(mcposrESS_1b_source, mcrotrESS_1b_source,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component ESS_1b_source (without coords transformations) */
  mcJumpTrace_ESS_1b_source:
  SIG_MESSAGE("ESS_1b_source (Trace)");
  mcDEBUG_COMP("ESS_1b_source")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp
  STORE_NEUTRON(2,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcNCounter[2]++;
  mcPCounter[2] += p;
  mcP2Counter[2] += p*p;
#define mccompcurname  ESS_1b_source
#define mccompcurtype  Virtual_input
#define mccompcurindex 2
#define read_block mccESS_1b_source_read_block
#define pos mccESS_1b_source_pos
#define nrows mccESS_1b_source_nrows
#define Offset mccESS_1b_source_Offset
#define rTable mccESS_1b_source_rTable
#define repeat_number mccESS_1b_source_repeat_number
#define filename_ncount mccESS_1b_source_filename_ncount
#define mean_vx mccESS_1b_source_mean_vx
#define mean_vy mccESS_1b_source_mean_vy
#define mean_vz mccESS_1b_source_mean_vz
#define mean_dx mccESS_1b_source_mean_dx
#define mean_dy mccESS_1b_source_mean_dy
#define mean_dz mccESS_1b_source_mean_dz
#define n_neutrons mccESS_1b_source_n_neutrons
#define min_x mccESS_1b_source_min_x
#define min_y mccESS_1b_source_min_y
#define min_z mccESS_1b_source_min_z
#define max_x mccESS_1b_source_max_x
#define max_y mccESS_1b_source_max_y
#define max_z mccESS_1b_source_max_z
#define min_vx mccESS_1b_source_min_vx
#define min_vy mccESS_1b_source_min_vy
#define min_vz mccESS_1b_source_min_vz
#define max_vx mccESS_1b_source_max_vx
#define max_vy mccESS_1b_source_max_vy
#define max_vz mccESS_1b_source_max_vz
#define first_block mccESS_1b_source_first_block
#define mean_x mccESS_1b_source_mean_x
#define mean_y mccESS_1b_source_mean_y
#define mean_z mccESS_1b_source_mean_z
#define end_reading mccESS_1b_source_end_reading
#define n_count_extrapolated mccESS_1b_source_n_count_extrapolated
#define repeat_cnt mccESS_1b_source_repeat_cnt
{   /* Declarations of ESS_1b_source=Virtual_input() SETTING parameters. */
char* filename = mccESS_1b_source_filename;
char* type = mccESS_1b_source_type;
MCNUM verbose = mccESS_1b_source_verbose;
MCNUM repeat_count = mccESS_1b_source_repeat_count;
MCNUM smooth = mccESS_1b_source_smooth;
#line 174 "/usr/local/lib/mcstas-2.0/sources/Virtual_input.comp"
{
  if (filename && strlen(filename)) {
    while (read_block && !end_reading) {
      /* read block and increase Offset for next reading */
      nrows = Virtual_input_Read_Input(filename, type, &rTable, &Offset);

      if (!nrows) { /* nrows is 0 if end of filename/no filename */
        if (!filename_ncount) {
          filename_ncount = (double)mcget_run_num();  /* ncount in filename */
          if (verbose)
            printf("Virtual_input: %s: filename '%s' contains %g events\n", NAME_CURRENT_COMP, filename, filename_ncount);
          /* set ncount from filename length */
          mcset_ncount(filename_ncount*repeat_cnt);
        }
        Offset = 0;       /* reposition to begining of filename */
        repeat_number++;  /* we start a new repeat_cnt loop */

        /* end of simulation if ... */
        if (repeat_number > repeat_cnt) {
          if (verbose)
            printf("Virtual_input: %s: Ending after %g events (%li repeat count)\n", NAME_CURRENT_COMP, (double)mcget_run_num(), (long)repeat_cnt);
          read_block=0; mcset_ncount(mcget_run_num()); pos=0;
          end_reading = 1;
        }
        /* else continue reading blocks */

      } else { /* block at Offset could be read */
        pos = 0;  /* position at begining of new block */
        read_block = 0;
      }
    }

      /* &p, &x, &y, &z, &vx, &vy, &vz, &t, &sx, &sy, &sz */
    if (!end_reading) {
      /* BEWARE! This is a non-standard way of using mcrestore_neutron, order of parameters
	 MUST be  &p, &x, &y, &z, &vx, &vy, &vz, &t, &sx, &sy,  &sz 
         to match the data format from the Virtual_input files */
      mcrestore_neutron(rTable.data,pos, &p, &x, &y, &z, &vx, &vy, &vz, &t, &sx, &sy,  &sz);

      if (first_block) {
        double v;
        mean_x  += p*x;  mean_y  += p*y;  mean_z  += p*z;
        mean_vx += p*vx; mean_vy += p*vy; mean_vz += p*vz;
        v = sqrt(vx*vx+vy*vy+vz*vz);
        if (v)
          { mean_dx += p*fabs(vx/v); mean_dy += p*fabs(vy/v); mean_dz += p*fabs(vz/v); }
        if (x  < min_x)  min_x  = x;
        if (y  < min_y)  min_y  = y;
        if (z  < min_z)  min_z  = z;
        if (vx < min_vx) min_vx = vx;
        if (vy < min_vy) min_vy = vy;
        if (vz < min_vz) min_vz = z;
        if (x  > max_x)  max_x  = x;
        if (y  > max_y)  max_y  = y;
        if (z  > max_z)  max_z  = z;
        if (vx > max_vx) max_vx = vx;
        if (vy > max_vy) max_vy = vy;
        if (vz > max_vz) max_vz = z;
        n_neutrons += p;
      }

      pos++;
      p /= repeat_cnt;
#ifdef USE_MPI
      /* We always repeat by the number of nodes in an MPI run */
      p /= mpi_node_count;
#endif

      SCATTER;

      if (pos >= nrows) { /* reached end of block */
        read_block = 1;
        if (first_block) {
          double mean_v;
          /* display statitics for 1st block */
          mean_x  /= n_neutrons;
          mean_y  /= n_neutrons;
          mean_z  /= n_neutrons;
          mean_vx /= n_neutrons;
          mean_vy /= n_neutrons;
          mean_vz /= n_neutrons;
          mean_dx /= n_neutrons;
          mean_dy /= n_neutrons;
          mean_dz /= n_neutrons;
          /* now estimates total ncount */
          mean_v = sqrt(mean_vx*mean_vx+mean_vy*mean_vy+mean_vz*mean_vz);
          n_count_extrapolated = (double)nrows*rTable.filesize/Offset;
          if (verbose) {
            double mean_k, mean_w=0, mean_L=0;

            mean_k = V2K*mean_v;
            if (mean_k) mean_L = 2*PI/mean_k;
            mean_w = VS2E*mean_v*mean_v;
            printf("McStas Virtual Source filename %s\nContains about %g events, intensity=%g\n", filename, n_count_extrapolated, n_neutrons*rTable.filesize/Offset);

            printf("  Source size (full width in [m]):      ");
            printf("    dX=%g dY=%g dZ=%g\n", max_x-min_x, max_y-min_y, max_z-min_z);
            printf("  Source center (in [m]):               ");
            printf("    X0=%g Y0=%g Z0=%g\n", mean_x, mean_y, mean_z);
            printf("  Beam divergence (full width in [deg]):");
            printf("    dVx=%g dVy=%g dVz=%g\n",
              atan(mean_dx)*RAD2DEG,
              atan(mean_dy)*RAD2DEG,
              atan(mean_dz)*RAD2DEG);
            printf("  Beam speed (in [m/s]):                ");
            printf("    Vx=%g Vy=%g Vz=%g\n", mean_vx, mean_vy, mean_vz);
            printf("  Beam mean energy:\n");
            printf("    speed=%g [m/s] energy=%g [meV]\n    wavelength=%g [Angs] wavevector=%g [Angs-1]\n", mean_v, mean_w, mean_L, mean_k);
          }
          /* set ncount from estimate with security margin 10 % */
          mcset_ncount(n_count_extrapolated*repeat_cnt*1.1);
        }
        first_block= 0;
      }
  #if defined (USE_MPI)
      if (smooth && n_count_extrapolated)
  #else
      if (smooth && repeat_number>1 && n_count_extrapolated)
  #endif
      {
        /* apply smmothing */
        x += randnorm()*(max_x-min_x)/n_count_extrapolated/2;
        y += randnorm()*(max_y-min_y)/n_count_extrapolated/2;
        z += randnorm()*(max_z-min_z)/n_count_extrapolated/2;
        vx += randnorm()*min_vx/n_count_extrapolated/2;
        vy += randnorm()*min_vx/n_count_extrapolated/2;
        vx += randnorm()*min_vx/n_count_extrapolated/2;
      }
    } else { ABSORB; }
  }
}
#line 12223 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of ESS_1b_source=Virtual_input() SETTING parameter declarations. */
#undef repeat_cnt
#undef n_count_extrapolated
#undef end_reading
#undef mean_z
#undef mean_y
#undef mean_x
#undef first_block
#undef max_vz
#undef max_vy
#undef max_vx
#undef min_vz
#undef min_vy
#undef min_vx
#undef max_z
#undef max_y
#undef max_x
#undef min_z
#undef min_y
#undef min_x
#undef n_neutrons
#undef mean_dz
#undef mean_dy
#undef mean_dx
#undef mean_vz
#undef mean_vy
#undef mean_vx
#undef filename_ncount
#undef repeat_number
#undef rTable
#undef Offset
#undef nrows
#undef pos
#undef read_block
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component v_mon [3] */
  mccoordschange(mcposrv_mon, mcrotrv_mon,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component v_mon (without coords transformations) */
  mcJumpTrace_v_mon:
  SIG_MESSAGE("v_mon (Trace)");
  mcDEBUG_COMP("v_mon")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp
  STORE_NEUTRON(3,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcNCounter[3]++;
  mcPCounter[3] += p;
  mcP2Counter[3] += p*p;
#define mccompcurname  v_mon
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define options mccv_mon_options
#define filename mccv_mon_filename
#define geometry mccv_mon_geometry
#define user1 mccv_mon_user1
#define user2 mccv_mon_user2
#define user3 mccv_mon_user3
#define username1 mccv_mon_username1
#define username2 mccv_mon_username2
#define username3 mccv_mon_username3
#define DEFS mccv_mon_DEFS
#define Vars mccv_mon_Vars
#define detector mccv_mon_detector
{   /* Declarations of v_mon=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccv_mon_xwidth;
MCNUM yheight = mccv_mon_yheight;
MCNUM zdepth = mccv_mon_zdepth;
MCNUM xmin = mccv_mon_xmin;
MCNUM xmax = mccv_mon_xmax;
MCNUM ymin = mccv_mon_ymin;
MCNUM ymax = mccv_mon_ymax;
MCNUM zmin = mccv_mon_zmin;
MCNUM zmax = mccv_mon_zmax;
MCNUM bins = mccv_mon_bins;
MCNUM min = mccv_mon_min;
MCNUM max = mccv_mon_max;
MCNUM restore_neutron = mccv_mon_restore_neutron;
MCNUM radius = mccv_mon_radius;
#line 301 "/usr/local/lib/mcstas-2.0/monitors/Monitor_nD.comp"
{
  double  XY=0;
  double  t0 = 0;
  double  t1 = 0;
  double  pp;
  int     intersect   = 0;
  char    Flag_Restore = 0;

  if (user1 != FLT_MAX) Vars.UserVariable1 = user1;
  if (user2 != FLT_MAX) Vars.UserVariable2 = user2;
  if (user3 != FLT_MAX) Vars.UserVariable3 = user3;

  /* this is done automatically
    STORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  */
  
  if (geometry && strlen(geometry))
  {
    /* determine intersections with object */
    intersect = off_intersect(&t0, &t1, NULL, NULL,
       x,y,z, vx, vy, vz, offdata );
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_SQUARE) /* square xy */
  {
    PROP_Z0;
    intersect = (x>=Vars.mxmin && x<=Vars.mxmax && y>=Vars.mymin && y<=Vars.mymax);
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_DISK)   /* disk xy */
  {
    PROP_Z0;
    intersect = ((x*x + y*y) <= Vars.Sphere_Radius*Vars.Sphere_Radius);
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) /* sphere */
  {
    intersect = sphere_intersect(&t0, &t1, x, y, z, vx, vy, vz, Vars.Sphere_Radius);
  /*      intersect = (intersect && t0 > 0); */
  }
  else if ((abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA)) /* cylinder */
  {
    intersect = cylinder_intersect(&t0, &t1, x, y, z, vx, vy, vz, Vars.Sphere_Radius, Vars.Cylinder_Height);
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX) /* box */
  {
    intersect = box_intersect(&t0, &t1, x, y, z, vx, vy, vz, fabs(Vars.mxmax-Vars.mxmin), fabs(Vars.mymax-Vars.mymin), fabs(Vars.mzmax-Vars.mzmin));
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_PREVIOUS) /* previous comp */
  { intersect = 1; }

  if (intersect)
  {
    if ((abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND) 
     || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA)
     || (geometry && strlen(geometry)) )
    {
      /* check if we have to remove the top/bottom with BANANA shape */
      if ((abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA) && (intersect != 1)) {
        double y0,y1;
        /* propagate to intersection point as temporary variable to check top/bottom */
        y0 = y+t0*vy; 
        y1 = y+t1*vy;
        if (fabs(y0) >= Vars.Cylinder_Height/2*0.99) t0 = t1;
        if (fabs(y1) >= Vars.Cylinder_Height/2*0.99) t1 = t0;
      }
      if (t0 < 0 && t1 > 0)
        t0 = t;  /* neutron was already inside ! */
      if (t1 < 0 && t0 > 0) /* neutron exit before entering !! */
        t1 = t;
      /* t0 is now time of incoming intersection with the detection area */
      if ((Vars.Flag_Shape < 0) && (t1 > 0))
        PROP_DT(t1); /* t1 outgoing beam */
      else
        PROP_DT(t0); /* t0 incoming beam */
      /* Final test if we are on lid / bottom of banana */
      if (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA) {
	if (fabs(y) >= Vars.Cylinder_Height/2*0.99) ABSORB;
      }
    }

    /* Now get the data to monitor: current or keep from PreMonitor */
    if (Vars.Flag_UsePreMonitor != 1)
    {
      Vars.cp  = p;
      Vars.cx  = x;
      Vars.cvx = vx;
      Vars.csx = sx;
      Vars.cy  = y;
      Vars.cvy = vy;
      Vars.csy = sy;
      Vars.cz  = z;
      Vars.cvz = vz;
      Vars.csz = sz;
      Vars.ct  = t;
    }

    if ((Vars.He3_pressure > 0) && (t1 != t0) && ((abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX)))
    {
      XY = exp(-7.417*Vars.He3_pressure*fabs(t1-t0)*2*PI*K2V);
      /* will monitor the absorbed part */
      Vars.cp *= 1-XY;
      /* and modify the neutron weight after monitor, only remains 1-p_detect */
      p *= XY;
    }

    if (Vars.Flag_capture)
    {
      XY = sqrt(Vars.cvx*Vars.cvx+Vars.cvy*Vars.cvy+Vars.cvz*Vars.cvz);
      XY *= V2K;
      if (XY != 0) XY = 2*PI/XY; /* lambda. lambda(2200 m/2) = 1.7985 Angs  */
      Vars.cp *= XY/1.7985;
    }

    pp = Monitor_nD_Trace(&DEFS, &Vars);
    if (pp==0.0)
    { ABSORB;
    }
    else
    { Vars.Nsum++;
      Vars.psum += pp;
      Vars.p2sum += pp*pp;
      SCATTER;
    }

    if (Vars.Flag_parallel) /* back to neutron state before detection */
      Flag_Restore = 1;
  } /* end if intersection */
  else {
    if (Vars.Flag_Absorb && !Vars.Flag_parallel) ABSORB;
    else Flag_Restore = 1;  /* no intersection, back to previous state */
  }

  if (Flag_Restore)
  {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 12505 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of v_mon=Monitor_nD() SETTING parameter declarations. */
#undef detector
#undef Vars
#undef DEFS
#undef username3
#undef username2
#undef username1
#undef user3
#undef user2
#undef user1
#undef geometry
#undef filename
#undef options
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component image [4] */
  mccoordschange(mcposrimage, mcrotrimage,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component image (without coords transformations) */
  mcJumpTrace_image:
  SIG_MESSAGE("image (Trace)");
  mcDEBUG_COMP("image")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp
  STORE_NEUTRON(4,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcNCounter[4]++;
  mcPCounter[4] += p;
  mcP2Counter[4] += p*p;
#define mccompcurname  image
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
#define nx mccimage_nx
#define ny mccimage_ny
#define restore_neutron mccimage_restore_neutron
#define PSD_N mccimage_PSD_N
#define PSD_p mccimage_PSD_p
#define PSD_p2 mccimage_PSD_p2
{   /* Declarations of image=PSD_monitor() SETTING parameters. */
char* filename = mccimage_filename;
MCNUM xmin = mccimage_xmin;
MCNUM xmax = mccimage_xmax;
MCNUM ymin = mccimage_ymin;
MCNUM ymax = mccimage_ymax;
MCNUM xwidth = mccimage_xwidth;
MCNUM yheight = mccimage_yheight;
#line 85 "/usr/local/lib/mcstas-2.0/monitors/PSD_monitor.comp"
{
    int i,j;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      i = floor((x - xmin)*nx/(xmax - xmin));
      j = floor((y - ymin)*ny/(ymax - ymin));
      PSD_N[i][j]++;
      PSD_p[i][j] += p;
      PSD_p2[i][j] += p*p;
      SCATTER;
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 12635 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of image=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_neutron
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component lamb [5] */
  mccoordschange(mcposrlamb, mcrotrlamb,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component lamb (without coords transformations) */
  mcJumpTrace_lamb:
  SIG_MESSAGE("lamb (Trace)");
  mcDEBUG_COMP("lamb")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp
  STORE_NEUTRON(5,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcNCounter[5]++;
  mcPCounter[5] += p;
  mcP2Counter[5] += p*p;
#define mccompcurname  lamb
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mcclamb_nL
#define L_N mcclamb_L_N
#define L_p mcclamb_L_p
#define L_p2 mcclamb_L_p2
{   /* Declarations of lamb=L_monitor() SETTING parameters. */
char* filename = mcclamb_filename;
MCNUM xmin = mcclamb_xmin;
MCNUM xmax = mcclamb_xmax;
MCNUM ymin = mcclamb_ymin;
MCNUM ymax = mcclamb_ymax;
MCNUM xwidth = mcclamb_xwidth;
MCNUM yheight = mcclamb_yheight;
MCNUM Lmin = mcclamb_Lmin;
MCNUM Lmax = mcclamb_Lmax;
MCNUM restore_neutron = mcclamb_restore_neutron;
#line 83 "/usr/local/lib/mcstas-2.0/monitors/L_monitor.comp"
{
    int i;
    double L;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      L = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
      i = floor((L-Lmin)*nL/(Lmax-Lmin));
      if(i >= 0 && i < nL)
      {
        L_N[i]++;
        L_p[i] += p;
        L_p2[i] += p*p;
        SCATTER;
      }
    } 
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 12764 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of lamb=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  mcabsorbAll:
  mcDEBUG_LEAVE()
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)
  /* Copy neutron state to global variables. */
  mcnx = mcnlx;
  mcny = mcnly;
  mcnz = mcnlz;
  mcnvx = mcnlvx;
  mcnvy = mcnlvy;
  mcnvz = mcnlvz;
  mcnt = mcnlt;
  mcnsx = mcnlsx;
  mcnsy = mcnlsy;
  mcnsz = mcnlsz;
  mcnp = mcnlp;

} /* end trace */

void mcsave(FILE *handle) {
  if (!handle) mcsiminfo_init(NULL);
  /* User component SAVE code. */

  /* User SAVE code for component 'Origin'. */
  SIG_MESSAGE("Origin (Save)");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
{   /* Declarations of Origin=Progress_bar() SETTING parameters. */
char* profile = mccOrigin_profile;
MCNUM percent = mccOrigin_percent;
MCNUM flag_save = mccOrigin_flag_save;
MCNUM minutes = mccOrigin_minutes;
#line 113 "/usr/local/lib/mcstas-2.0/misc/Progress_bar.comp"
{
  MPI_MASTER(fprintf(stdout, "\nSave [%s]\n", mcinstrument_name););
  if (profile && strlen(profile)) {
    char filename[256];
    if (!strlen(profile)) strcpy(filename, mcinstrument_name);
    else strcpy(filename, profile);
    DETECTOR_OUT_1D(
        "Intensity profiler",
        "Component index [1]",
        "Intensity",
        "prof", 1, mcNUMCOMP, mcNUMCOMP-1,
        &mcNCounter[1],&mcPCounter[1],&mcP2Counter[1],
        filename);

  }
}
#line 12860 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'v_mon'. */
  SIG_MESSAGE("v_mon (Save)");
#define mccompcurname  v_mon
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define options mccv_mon_options
#define filename mccv_mon_filename
#define geometry mccv_mon_geometry
#define user1 mccv_mon_user1
#define user2 mccv_mon_user2
#define user3 mccv_mon_user3
#define username1 mccv_mon_username1
#define username2 mccv_mon_username2
#define username3 mccv_mon_username3
#define DEFS mccv_mon_DEFS
#define Vars mccv_mon_Vars
#define detector mccv_mon_detector
{   /* Declarations of v_mon=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccv_mon_xwidth;
MCNUM yheight = mccv_mon_yheight;
MCNUM zdepth = mccv_mon_zdepth;
MCNUM xmin = mccv_mon_xmin;
MCNUM xmax = mccv_mon_xmax;
MCNUM ymin = mccv_mon_ymin;
MCNUM ymax = mccv_mon_ymax;
MCNUM zmin = mccv_mon_zmin;
MCNUM zmax = mccv_mon_zmax;
MCNUM bins = mccv_mon_bins;
MCNUM min = mccv_mon_min;
MCNUM max = mccv_mon_max;
MCNUM restore_neutron = mccv_mon_restore_neutron;
MCNUM radius = mccv_mon_radius;
#line 438 "/usr/local/lib/mcstas-2.0/monitors/Monitor_nD.comp"
{
  /* save results, but do not free pointers */
  detector = Monitor_nD_Save(&DEFS, &Vars);
}
#line 12906 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of v_mon=Monitor_nD() SETTING parameter declarations. */
#undef detector
#undef Vars
#undef DEFS
#undef username3
#undef username2
#undef username1
#undef user3
#undef user2
#undef user1
#undef geometry
#undef filename
#undef options
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'image'. */
  SIG_MESSAGE("image (Save)");
#define mccompcurname  image
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
#define nx mccimage_nx
#define ny mccimage_ny
#define restore_neutron mccimage_restore_neutron
#define PSD_N mccimage_PSD_N
#define PSD_p mccimage_PSD_p
#define PSD_p2 mccimage_PSD_p2
{   /* Declarations of image=PSD_monitor() SETTING parameters. */
char* filename = mccimage_filename;
MCNUM xmin = mccimage_xmin;
MCNUM xmax = mccimage_xmax;
MCNUM ymin = mccimage_ymin;
MCNUM ymax = mccimage_ymax;
MCNUM xwidth = mccimage_xwidth;
MCNUM yheight = mccimage_yheight;
#line 103 "/usr/local/lib/mcstas-2.0/monitors/PSD_monitor.comp"
{
    DETECTOR_OUT_2D(
        "PSD monitor",
        "X position [cm]",
        "Y position [cm]",
        xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
        nx, ny,
        &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
        filename);
}
#line 12954 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of image=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_neutron
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'lamb'. */
  SIG_MESSAGE("lamb (Save)");
#define mccompcurname  lamb
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mcclamb_nL
#define L_N mcclamb_L_N
#define L_p mcclamb_L_p
#define L_p2 mcclamb_L_p2
{   /* Declarations of lamb=L_monitor() SETTING parameters. */
char* filename = mcclamb_filename;
MCNUM xmin = mcclamb_xmin;
MCNUM xmax = mcclamb_xmax;
MCNUM ymin = mcclamb_ymin;
MCNUM ymax = mcclamb_ymax;
MCNUM xwidth = mcclamb_xwidth;
MCNUM yheight = mcclamb_yheight;
MCNUM Lmin = mcclamb_Lmin;
MCNUM Lmax = mcclamb_Lmax;
MCNUM restore_neutron = mcclamb_restore_neutron;
#line 105 "/usr/local/lib/mcstas-2.0/monitors/L_monitor.comp"
{
    DETECTOR_OUT_1D(
        "Wavelength monitor",
        "Wavelength [AA]",
        "Intensity",
        "L", Lmin, Lmax, nL,
        &L_N[0],&L_p[0],&L_p2[0],
        filename);
}
#line 12996 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of lamb=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  if (!handle) mcsiminfo_close(); 
} /* end save */
void mcfinally(void) {
  /* User component FINALLY code. */
  mcsiminfo_init(NULL);
  mcsave(mcsiminfo_file); /* save data when simulation ends */

  /* User FINALLY code for component 'Origin'. */
  SIG_MESSAGE("Origin (Finally)");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
{   /* Declarations of Origin=Progress_bar() SETTING parameters. */
char* profile = mccOrigin_profile;
MCNUM percent = mccOrigin_percent;
MCNUM flag_save = mccOrigin_flag_save;
MCNUM minutes = mccOrigin_minutes;
#line 131 "/usr/local/lib/mcstas-2.0/misc/Progress_bar.comp"
{
  time_t NowTime;
  time(&NowTime);
  fprintf(stdout, "\nFinally [%s/%s]. Time: ", mcinstrument_name, mcdirname ? mcdirname : ".");
  if (difftime(NowTime,StartTime) < 60.0)
    fprintf(stdout, "%g [s] ", difftime(NowTime,StartTime));
  else if (difftime(NowTime,StartTime) > 3600.0)
    fprintf(stdout, "%g [h] ", difftime(NowTime,StartTime)/3660.0);
  else
    fprintf(stdout, "%g [min] ", difftime(NowTime,StartTime)/60.0);
  fprintf(stdout, "\n");
}
#line 13039 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] Origin\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] Origin=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
  /* User FINALLY code for component 'ESS_1b_source'. */
  SIG_MESSAGE("ESS_1b_source (Finally)");
#define mccompcurname  ESS_1b_source
#define mccompcurtype  Virtual_input
#define mccompcurindex 2
#define read_block mccESS_1b_source_read_block
#define pos mccESS_1b_source_pos
#define nrows mccESS_1b_source_nrows
#define Offset mccESS_1b_source_Offset
#define rTable mccESS_1b_source_rTable
#define repeat_number mccESS_1b_source_repeat_number
#define filename_ncount mccESS_1b_source_filename_ncount
#define mean_vx mccESS_1b_source_mean_vx
#define mean_vy mccESS_1b_source_mean_vy
#define mean_vz mccESS_1b_source_mean_vz
#define mean_dx mccESS_1b_source_mean_dx
#define mean_dy mccESS_1b_source_mean_dy
#define mean_dz mccESS_1b_source_mean_dz
#define n_neutrons mccESS_1b_source_n_neutrons
#define min_x mccESS_1b_source_min_x
#define min_y mccESS_1b_source_min_y
#define min_z mccESS_1b_source_min_z
#define max_x mccESS_1b_source_max_x
#define max_y mccESS_1b_source_max_y
#define max_z mccESS_1b_source_max_z
#define min_vx mccESS_1b_source_min_vx
#define min_vy mccESS_1b_source_min_vy
#define min_vz mccESS_1b_source_min_vz
#define max_vx mccESS_1b_source_max_vx
#define max_vy mccESS_1b_source_max_vy
#define max_vz mccESS_1b_source_max_vz
#define first_block mccESS_1b_source_first_block
#define mean_x mccESS_1b_source_mean_x
#define mean_y mccESS_1b_source_mean_y
#define mean_z mccESS_1b_source_mean_z
#define end_reading mccESS_1b_source_end_reading
#define n_count_extrapolated mccESS_1b_source_n_count_extrapolated
#define repeat_cnt mccESS_1b_source_repeat_cnt
{   /* Declarations of ESS_1b_source=Virtual_input() SETTING parameters. */
char* filename = mccESS_1b_source_filename;
char* type = mccESS_1b_source_type;
MCNUM verbose = mccESS_1b_source_verbose;
MCNUM repeat_count = mccESS_1b_source_repeat_count;
MCNUM smooth = mccESS_1b_source_smooth;
#line 307 "/usr/local/lib/mcstas-2.0/sources/Virtual_input.comp"
{
  Table_Free(&rTable);
  if (!filename_ncount) {
    printf("Warning: Virtual_input: %s: filename '%s' was not used entirely.\n"
           "               Intensities may be wrong. Increase ncount value\n",
           NAME_CURRENT_COMP, filename);
  } else {
    double tmp;
    tmp = (double)mcget_ncount()/filename_ncount;
    if (fabs(rint(tmp)/tmp-1) > 0.02)
      printf("Warning: Virtual_input: %s: simulation finished in the middle of filename '%s'\n"
             "               ncount=%g but filename contains %g events\n"
             "               Intensities may be wrong.\n"
             "               Increase ncount value to %g or higher.\n",
	     NAME_CURRENT_COMP, filename, (double)mcget_ncount(), filename_ncount, filename_ncount*repeat_cnt);
    if (mcget_ncount() < filename_ncount*repeat_cnt)
      printf("Warning: Virtual_input: %s: not all source %s repetitions were generated.\n"
             "               Intensities may be wrong.\n"
             "               Increase ncount value to %g or higher.\n",
             NAME_CURRENT_COMP, filename, filename_ncount*repeat_cnt);
  }
}
#line 13117 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of ESS_1b_source=Virtual_input() SETTING parameter declarations. */
#undef repeat_cnt
#undef n_count_extrapolated
#undef end_reading
#undef mean_z
#undef mean_y
#undef mean_x
#undef first_block
#undef max_vz
#undef max_vy
#undef max_vx
#undef min_vz
#undef min_vy
#undef min_vx
#undef max_z
#undef max_y
#undef max_x
#undef min_z
#undef min_y
#undef min_x
#undef n_neutrons
#undef mean_dz
#undef mean_dy
#undef mean_dx
#undef mean_vz
#undef mean_vy
#undef mean_vx
#undef filename_ncount
#undef repeat_number
#undef rTable
#undef Offset
#undef nrows
#undef pos
#undef read_block
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] ESS_1b_source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] ESS_1b_source=Virtual_input()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
  /* User FINALLY code for component 'v_mon'. */
  SIG_MESSAGE("v_mon (Finally)");
#define mccompcurname  v_mon
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define options mccv_mon_options
#define filename mccv_mon_filename
#define geometry mccv_mon_geometry
#define user1 mccv_mon_user1
#define user2 mccv_mon_user2
#define user3 mccv_mon_user3
#define username1 mccv_mon_username1
#define username2 mccv_mon_username2
#define username3 mccv_mon_username3
#define DEFS mccv_mon_DEFS
#define Vars mccv_mon_Vars
#define detector mccv_mon_detector
{   /* Declarations of v_mon=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccv_mon_xwidth;
MCNUM yheight = mccv_mon_yheight;
MCNUM zdepth = mccv_mon_zdepth;
MCNUM xmin = mccv_mon_xmin;
MCNUM xmax = mccv_mon_xmax;
MCNUM ymin = mccv_mon_ymin;
MCNUM ymax = mccv_mon_ymax;
MCNUM zmin = mccv_mon_zmin;
MCNUM zmax = mccv_mon_zmax;
MCNUM bins = mccv_mon_bins;
MCNUM min = mccv_mon_min;
MCNUM max = mccv_mon_max;
MCNUM restore_neutron = mccv_mon_restore_neutron;
MCNUM radius = mccv_mon_radius;
#line 444 "/usr/local/lib/mcstas-2.0/monitors/Monitor_nD.comp"
{
  /* free pointers */
  Monitor_nD_Finally(&DEFS, &Vars);
}
#line 13195 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of v_mon=Monitor_nD() SETTING parameter declarations. */
#undef detector
#undef Vars
#undef DEFS
#undef username3
#undef username2
#undef username1
#undef user3
#undef user2
#undef user1
#undef geometry
#undef filename
#undef options
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] v_mon\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] v_mon=Monitor_nD()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] image\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] image=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] lamb\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] lamb=L_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
  mcsiminfo_close(); 
} /* end finally */
#define magnify mcdis_magnify
#define line mcdis_line
#define dashed_line mcdis_dashed_line
#define multiline mcdis_multiline
#define rectangle mcdis_rectangle
#define box mcdis_box
#define circle mcdis_circle
void mcdisplay(void) {
  printf("MCDISPLAY: start\n");
  /* Components MCDISPLAY code. */

  /* MCDISPLAY code for component 'ESS_1b_source'. */
  SIG_MESSAGE("ESS_1b_source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "ESS_1b_source");
#define mccompcurname  ESS_1b_source
#define mccompcurtype  Virtual_input
#define mccompcurindex 2
#define read_block mccESS_1b_source_read_block
#define pos mccESS_1b_source_pos
#define nrows mccESS_1b_source_nrows
#define Offset mccESS_1b_source_Offset
#define rTable mccESS_1b_source_rTable
#define repeat_number mccESS_1b_source_repeat_number
#define filename_ncount mccESS_1b_source_filename_ncount
#define mean_vx mccESS_1b_source_mean_vx
#define mean_vy mccESS_1b_source_mean_vy
#define mean_vz mccESS_1b_source_mean_vz
#define mean_dx mccESS_1b_source_mean_dx
#define mean_dy mccESS_1b_source_mean_dy
#define mean_dz mccESS_1b_source_mean_dz
#define n_neutrons mccESS_1b_source_n_neutrons
#define min_x mccESS_1b_source_min_x
#define min_y mccESS_1b_source_min_y
#define min_z mccESS_1b_source_min_z
#define max_x mccESS_1b_source_max_x
#define max_y mccESS_1b_source_max_y
#define max_z mccESS_1b_source_max_z
#define min_vx mccESS_1b_source_min_vx
#define min_vy mccESS_1b_source_min_vy
#define min_vz mccESS_1b_source_min_vz
#define max_vx mccESS_1b_source_max_vx
#define max_vy mccESS_1b_source_max_vy
#define max_vz mccESS_1b_source_max_vz
#define first_block mccESS_1b_source_first_block
#define mean_x mccESS_1b_source_mean_x
#define mean_y mccESS_1b_source_mean_y
#define mean_z mccESS_1b_source_mean_z
#define end_reading mccESS_1b_source_end_reading
#define n_count_extrapolated mccESS_1b_source_n_count_extrapolated
#define repeat_cnt mccESS_1b_source_repeat_cnt
{   /* Declarations of ESS_1b_source=Virtual_input() SETTING parameters. */
char* filename = mccESS_1b_source_filename;
char* type = mccESS_1b_source_type;
MCNUM verbose = mccESS_1b_source_verbose;
MCNUM repeat_count = mccESS_1b_source_repeat_count;
MCNUM smooth = mccESS_1b_source_smooth;
#line 331 "/usr/local/lib/mcstas-2.0/sources/Virtual_input.comp"
{
  mcdis_box(mean_x, mean_y, mean_z, max_x-min_x, max_y-min_y, max_z-min_z);
  /* a line in the main beam direction */ 
  line(mean_x, mean_y, mean_z,mean_dx,mean_dy,mean_dz);
}
#line 13283 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of ESS_1b_source=Virtual_input() SETTING parameter declarations. */
#undef repeat_cnt
#undef n_count_extrapolated
#undef end_reading
#undef mean_z
#undef mean_y
#undef mean_x
#undef first_block
#undef max_vz
#undef max_vy
#undef max_vx
#undef min_vz
#undef min_vy
#undef min_vx
#undef max_z
#undef max_y
#undef max_x
#undef min_z
#undef min_y
#undef min_x
#undef n_neutrons
#undef mean_dz
#undef mean_dy
#undef mean_dx
#undef mean_vz
#undef mean_vy
#undef mean_vx
#undef filename_ncount
#undef repeat_number
#undef rTable
#undef Offset
#undef nrows
#undef pos
#undef read_block
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'v_mon'. */
  SIG_MESSAGE("v_mon (McDisplay)");
  printf("MCDISPLAY: component %s\n", "v_mon");
#define mccompcurname  v_mon
#define mccompcurtype  Monitor_nD
#define mccompcurindex 3
#define options mccv_mon_options
#define filename mccv_mon_filename
#define geometry mccv_mon_geometry
#define user1 mccv_mon_user1
#define user2 mccv_mon_user2
#define user3 mccv_mon_user3
#define username1 mccv_mon_username1
#define username2 mccv_mon_username2
#define username3 mccv_mon_username3
#define DEFS mccv_mon_DEFS
#define Vars mccv_mon_Vars
#define detector mccv_mon_detector
{   /* Declarations of v_mon=Monitor_nD() SETTING parameters. */
MCNUM xwidth = mccv_mon_xwidth;
MCNUM yheight = mccv_mon_yheight;
MCNUM zdepth = mccv_mon_zdepth;
MCNUM xmin = mccv_mon_xmin;
MCNUM xmax = mccv_mon_xmax;
MCNUM ymin = mccv_mon_ymin;
MCNUM ymax = mccv_mon_ymax;
MCNUM zmin = mccv_mon_zmin;
MCNUM zmax = mccv_mon_zmax;
MCNUM bins = mccv_mon_bins;
MCNUM min = mccv_mon_min;
MCNUM max = mccv_mon_max;
MCNUM restore_neutron = mccv_mon_restore_neutron;
MCNUM radius = mccv_mon_radius;
#line 450 "/usr/local/lib/mcstas-2.0/monitors/Monitor_nD.comp"
{
  Monitor_nD_McDisplay(&DEFS, &Vars);
}
#line 13359 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of v_mon=Monitor_nD() SETTING parameter declarations. */
#undef detector
#undef Vars
#undef DEFS
#undef username3
#undef username2
#undef username1
#undef user3
#undef user2
#undef user1
#undef geometry
#undef filename
#undef options
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'image'. */
  SIG_MESSAGE("image (McDisplay)");
  printf("MCDISPLAY: component %s\n", "image");
#define mccompcurname  image
#define mccompcurtype  PSD_monitor
#define mccompcurindex 4
#define nx mccimage_nx
#define ny mccimage_ny
#define restore_neutron mccimage_restore_neutron
#define PSD_N mccimage_PSD_N
#define PSD_p mccimage_PSD_p
#define PSD_p2 mccimage_PSD_p2
{   /* Declarations of image=PSD_monitor() SETTING parameters. */
char* filename = mccimage_filename;
MCNUM xmin = mccimage_xmin;
MCNUM xmax = mccimage_xmax;
MCNUM ymin = mccimage_ymin;
MCNUM ymax = mccimage_ymax;
MCNUM xwidth = mccimage_xwidth;
MCNUM yheight = mccimage_yheight;
#line 115 "/usr/local/lib/mcstas-2.0/monitors/PSD_monitor.comp"
{
  magnify("xy");
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 13406 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of image=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef restore_neutron
#undef ny
#undef nx
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'lamb'. */
  SIG_MESSAGE("lamb (McDisplay)");
  printf("MCDISPLAY: component %s\n", "lamb");
#define mccompcurname  lamb
#define mccompcurtype  L_monitor
#define mccompcurindex 5
#define nL mcclamb_nL
#define L_N mcclamb_L_N
#define L_p mcclamb_L_p
#define L_p2 mcclamb_L_p2
{   /* Declarations of lamb=L_monitor() SETTING parameters. */
char* filename = mcclamb_filename;
MCNUM xmin = mcclamb_xmin;
MCNUM xmax = mcclamb_xmax;
MCNUM ymin = mcclamb_ymin;
MCNUM ymax = mcclamb_ymax;
MCNUM xwidth = mcclamb_xwidth;
MCNUM yheight = mcclamb_yheight;
MCNUM Lmin = mcclamb_Lmin;
MCNUM Lmax = mcclamb_Lmax;
MCNUM restore_neutron = mcclamb_restore_neutron;
#line 116 "/usr/local/lib/mcstas-2.0/monitors/L_monitor.comp"
{
  magnify("xy");
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 13448 "/home/matt/Desktop/WORK/full_1b_moderator/source_test.c"
}   /* End of lamb=L_monitor() SETTING parameter declarations. */
#undef L_p2
#undef L_p
#undef L_N
#undef nL
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  printf("MCDISPLAY: end\n");
} /* end display */
#undef magnify
#undef line
#undef dashed_line
#undef multiline
#undef rectangle
#undef box
#undef circle
/* end of generated C code /home/matt/Desktop/WORK/full_1b_moderator/source_test.c */
