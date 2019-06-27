/*
 * parameter.h
 *
 *  Created on: 2019/03/13
 *      Author: yamasita
 */

#ifndef PARAMETER_H_
#define PARAMETER_H_

/* Work base directory */
#define WORKBASE         "../images/"

/* Default parameters: switch whether making template or matching, file names. They can be changed by command line parameters */
#define MKTEMPLATE     0 /* 0: matching, 1: make template */
#define TEMPLATEDNAME "img1.template"  /* Template data file name. It used for both input (matching) and output (making template) */
#define IMAGENAME     "img2"           /* Input image file name */

/* Floading point type (float or double) */
#define FPTYPE float
//#define FPTYPE double
/* Selection of image database */

#define GRAF /* GRAF, BOAT, BARK, UBC, MNIST, or LENNA*/

#define BLURPRE "/" /* Input image is blurred prefix "" or "blur3/" or "blur6" */
//#define BLURPRE "blur3/"
#define INPRE  ""
#define OUTPRE "OutMatch/"
#define TMPLPRE "OutMatch/"

/* To make template */
#define NTEMPLATEROT       4       /* Number of types of rotation for template image 0 or 4 */
//#define NTEMPLATEROT      128
#define NTEMPLATEENLONG    2 /* Number of types of anisotropic elongation for template image 0 or 2 */
#define TEMPLATEENCOEF  0.25       /* The coefficient for anisotropic elongation */
#define TEMPLSCALERATIO 1.41421356 /* the next scale = the present scale * SCALERATIO */

/* Image information */
#ifdef GRAF
#define IMGPRE           "Graffiti/"
#define TEMPLATEX       430
#define TEMPLATEY       310
#define TEMPLSCALEINIT  150.0     /* The initial scale */
#define TEMPLSCALEMAX   250.0     /* The maximum scale to be handled */
#define NX              800      /* Horizontal # of pixel in image */
#define NY              640      /* Vertical # of pixel in image */
#define SCALEINIT       80.0     /* The initial scale */
#define SCALEMAX        161.0    /* The maximum scale to be handled */
#endif

#ifdef BOAT
#define IMGPRE           "Boat/"
#define TEMPLATEX       425
#define TEMPLATEY       340
#define TEMPLSCALEINIT  200.0     /* The initial scale */
#define TEMPLSCALEMAX   300.0     /* The maximum scale to be handled */
#define NX              850      /* Horizontal # of pixel in image */
#define NY              680      /* Vertical # of pixel in image */
#define SCALEINIT       75.0     /* The initial scale */
#define SCALEMAX        301.0    /* The maximum scale to be handled */
#endif

#ifdef BARK
#define IMGPRE          "Bark/"
#define TEMPLATEX       360
#define TEMPLATEY       255
#define TEMPLSCALEINIT  100.0     /* The initial scale for Boat */
#define TEMPLSCALEMAX   150.0     /* The maximum scale to be handled for Boat */
#define NX              765      /* Horizontal # of pixel in image */
#define NY              512      /* Vertical # of pixel in image */
#define SCALEINIT       25.0     /* The initial scale for BOAT */
#define SCALEMAX        101.0    /* The maximum scale to be handled BOAT */
#endif

#ifdef UBC
#define IMGPRE           "UBC/"
#define TEMPLATEX       410
#define TEMPLATEY       240
#define TEMPLSCALEINIT  100.0     /* The initial scale for Boat */
#define TEMPLSCALEMAX   150.0     /* The maximum scale to be handled for Boat */
#define NX              800      /* Horizontal # of pixel in image */
#define NY              640      /* Vertical # of pixel in image */
#define SCALEINIT       50.0     /* The initial scale for BOAT */
#define SCALEMAX        201.0    /* The maximum scale to be handled BOAT */
#endif


#ifdef LENNA  /* This file is used only for debug (Not matching) */
#define IMGPRE           "Lenna/"
#define TEMPLATEX       256
#define TEMPLATEY       256
#define TEMPLSCALEINIT   80.0     /* The initial scale */
#define TEMPLSCALEMAX    81.0     /* The maximum scale to be handled */
#define NX              512      /* Horizontal # of pixel in image */
#define NY              512      /* Vertical # of pixel in image */
#define SCALEINIT        80.0     /* The initial scale for */
#define SCALEMAX         81.0    /* The maximum scale to be handled */
#endif

#ifdef MNIST
#define IMGPRE           "mtmp/"
#define TEMPLATEX       360
#define TEMPLATEY       255
#define TEMPLSCALEINIT  8.0     /* The initial scale for Boat */
#define TEMPLSCALEMAX   14.0     /* The maximum scale to be handled for Boat */
#define NX              765      /* Horizontal # of pixel in image */
#define NY              512      /* Vertical # of pixel in image */
#define SCALEINIT       30.0     /* The initial scale for BOAT */
#define SCALEMAX        61.1    /* The maximum scale to be handled BOAT */
#endif


/* Parameters for calculating feature */
#define PGRAD      4        /* P to calculate of gradient and blur*/
#define PBLUR      2        /* P to calculate  of only blur */
#define PANGLE     4        /* P to express edge angle */
#define SIGMA2     1.1      /* sigma when P = 2 */
#define SIGMA4     0.80     /* sigma when P = 4 */
#define MAXK       1024     /* The maximum K (to allocate memory) */
#define SCALERATIO 2.0      /* The next scale = the present scale * SCALERATIO */
#define GRADRATIO  0.5      /* Ratio for blur of gradient comparing to simple blur */
#define STEPRATIO  0.2      /* The sample points is scale * STEPRATIO */
#define GRIDSIZE   4        /* The grid size for feature vector = GRIDSIZE x GRIDSIZE */
#define NMAXTHETA  16       /* Number for thetas to search the maximum point of histogram */
//#define NMAXTHETA  512
#define HISTREG              /* Use regularization for edge histogram */
#define EDGEPOWER  1.0e-12   /* Regularization parameter or minimum power of gradient */

/* Calculate large scale or approximate from small scale */
#define USELARGESCALE  /* Define Calculate large scale. Not define: Approximate from small scale*/
/* Parameter for large scale approximation */
#define NAPPROPOINT     (N0AP + N1AP + N2AP + N3AP)  /* 1 + 8 + 16 + 24 = 49*/
#define N0AP   1
#define N1AP   8
#define N2AP   16
#define N3AP   24
#define R0AP   0.0
#define R1AP   0.6
#define R2AP   1.2
#define R3AP   1.8
#define A0AP   0.73975
#define A1AP   0.59329
#define A2AP   0.35835
#define A3AP   0.14622

/* Fixed values */
#define WHITE          255
#define BLACK          0
#define MAX_BUFFERSIZE 256                /* Buffer size to read text from a file */
#define MAX_FILENAME   256                /* Filename length limit  */
#define PI             3.141592653589793
#define ROUNDFRAC      0.49999            /* For round */
#define EPS            1.0e-10

#endif /* PARAMETER_H_ */
