// ===========================================================================
// 
// warpingSIMDSmallTest.C --
// simplified test code for WarpingCompound class
// 
// This source code file is part of the following software:
// 
//    - the low-level C++ template SIMD library
//    - the SIMD implementation of the MinWarping and the 2D-Warping methods 
//      for local visual homing.
// 
// The software is provided based on the accompanying license agreement
// in the file LICENSE or LICENSE.doc. The software is provided "as is"
// without any warranty by the licensor and without any liability of the
// licensor, and the software may not be distributed by the licensee; see
// the license agreement for details.
// 
// (C) Ralf Möller
//     Computer Engineering
//     Faculty of Technology
//     Bielefeld University
//     www.ti.uni-bielefeld.de
// 
// ===========================================================================

// This simple example demonstrates the initialization of MinWarping
// data structures and the computation (in a class). It is restricted
// to MinWarping with NSAD2 distance measure, compass acceleration,
// and double search; for all alternatives, see warpingSIMDTest.C. No
// random azimuthal rotation is applied to the input images.

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <cmath>
#include <malloc.h>
#include <iostream>
#include <vector>
#include "SIMDImageFunctions.H"
#include "WarpingCompound.H"
#include "WarpingSPSComputationCollection.H"
#include "MinWarpingSearchCollection.H"
#include "WarpingSearchCollection.H"

using namespace ns_simd;

// ===========================================================================
// types used for Warping processing
// ===========================================================================

// types
#define IMGTYPE SIMDByte
#define PROCTYPE SIMDFloat
#define MEASTYPE SIMDShort
#define SPSTYPE SIMDByte
#define MATCHTYPE SIMDWord
// COMPASSTYPE needs to be signed (uses hadds which only exists for signed)
// for MATCHTYPE == SIMDWord, we can't use SIMDShort since the pixelScale is
// to large, therefore we use COMPASSTYPE = SIMDInt
#define COMPASSTYPE SIMDInt

#ifdef _SIMD_VEC_64_FULL_AVAIL_
# define SW 64
# define SA 64
#else
# ifdef _SIMD_VEC_32_FULL_AVAIL_
#  define SW 32
#  define SA 32
# else
#  define SW 16
#  define SA 16
# endif
#endif

// ===========================================================================
// complete warping initialization and computation
// ===========================================================================

template <typename ImgType, typename ProcType,
	  typename MeasType, typename SPSType,
	  typename MatchType, typename CompassType,
	  int SIMD_WIDTH, int SIMD_ALIGN>
class WarpingExample
{
protected:

  WarpingSPSComputation<ImgType,ProcType,MeasType,SPSType,
			SIMD_WIDTH,SIMD_ALIGN> *spsComp;
  MinWarpingPartial<SPSType,MatchType,
		    SIMD_WIDTH,SIMD_ALIGN> *partialSearcherMinWarping;
  WarpingCompound<SPSType,MatchType,CompassType,
		  SIMD_WIDTH,SIMD_ALIGN> *wsc;

public:

  // number of scale planes
  int nScalePlanes;
  // maximal scale factor
  double maxScaleFactor;
  // maximal threshold
  double maxThresholdMinWarping;
  // restriction of min-warping template to certain values of rho;
  // (0.0, 100.0) means that there is no restriction
  double rhoMinMinWarping;
  double rhoMaxMinWarping;
  double rhoMinWarping;
  double rhoMaxWarping;
  int nRhoWarping;
  // switch for double search (second run with exchanged snapshot/current view)
  // affects both partial and full search
  int doubleSearch;
  // interpolation: interpolation used in image magnification 
  // 0: nearest neighbor, 1: linear
  int interpolation;
  // parameter for compass acceleration
  double psiFraction;

  WarpingExample(int w, int nAlpha, int nPsi) :
    // 31. Jul 17 (rm): BUGFIX: wsc() -> wsc(0)
    spsComp(0), partialSearcherMinWarping(0), wsc(0),
    nScalePlanes(9), maxScaleFactor(2.0), maxThresholdMinWarping(2.5),
    rhoMinMinWarping(0.0), rhoMaxMinWarping(100.0),
    rhoMinWarping(0.0), rhoMaxWarping(1.0),
    nRhoWarping(20),
    doubleSearch(1),
    interpolation(0),
    psiFraction(0.3)
  {
    // spsComp and partialSearcherMinWarping don't depend on parameters
    // first phase: NSAD2
    spsComp = 
      new WarpingSPSComputationEdgeNSAD2<ImgType,ProcType,MeasType,SPSType,
					 SIMD_WIDTH,SIMD_ALIGN>();
    // second phase: MinWarping partial searcher
    partialSearcherMinWarping = 
      new MinWarpingPartial_XPAY<SPSType,MatchType,
				 SIMD_WIDTH,SIMD_ALIGN>();
    // parameter-dependent initialization
    init(w, nAlpha, nPsi);
  }

  // parameter members can be changed and a new compound can then be
  // created by calling init()
  void
  init(int w, int nAlpha, int nPsi)
  {
    // compound
    if (wsc) delete wsc;
    wsc =
      new WarpingCompound<SPSType,MatchType,CompassType,SIMD_WIDTH,SIMD_ALIGN>
      (w, nAlpha, nPsi,
       nScalePlanes, maxScaleFactor, maxThresholdMinWarping,
       rhoMinMinWarping, rhoMaxMinWarping,
       rhoMinWarping, rhoMaxWarping, nRhoWarping);
  }
	
  ~WarpingExample()
  {
    delete spsComp;
    delete partialSearcherMinWarping;
    delete wsc;
  }
  
  void run(SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &ss,
	   SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &cv,
	   const std::vector<double> &pixelScale, double postScale,
	   double &alphaMin, double &psiMin, double &dMin)
  {
    // first phase: NSAD2
    // (why can't we specify <ProcType,MeasType> here as in warpingSIMDTest.C?)
    // 27. Feb 18 (rm): vert.res. and hor. taken from images
    wsc->computeSPS(ss, cv, 
		    interpolation, pixelScale, postScale,
		    *spsComp);
    // second phase: MinWarping with compass acceleration
    wsc->minWarping.compassAcceleration
      (*partialSearcherMinWarping, doubleSearch, psiFraction);
    // find best match (without fine search)
    int iAlphaMin, iPsiMin;	// ignored
    wsc->minWarping.bestMatchFull(iAlphaMin, alphaMin, iPsiMin, psiMin, dMin);
  }
 
};

// ===========================================================================
// main
// ===========================================================================

// database
#define DATAPATH "DATA2"
#define BASE "living3"
#define SS "day"
#define CV "day"
#define SUFFIX "Hh288sh"
#define BW "0.10"
// grid coordinates of snapshot and current view
#define SSX 5
#define SSY 4
#define CVX 8
#define CVY 2
// number of search steps in alpha and phi direction
// note that there are restrictions on the choice depending on image width
#define NALPHA 96
#define NPSI 96
// pixelScale: scaling of pixel values when input images are converted
// into into internal representation, 16 for w = 288, same for larger
// height with same bww (for edge measures)
#define PIXELSCALE 16
// postScale: scaling of computed image distance values, use
// idealPostScale function to determine optimal value, ASC 100, NSAD
// 200 for w = 288, adjust for other w!
#define POSTSCALE 200

// should produce the same results as runWarpingSIMDTestVis for
// day-day combination, with
// setenv rotateImages 0
// setenv compassAcceleration 1
// setenv fineSearch 0

int
main()
{
  char dataPath[256], dbFile[275];
  char ssPath[274], cvPath[274], ssFile[292], cvFile[292];
  int xmin, xmax, ymin, ymax, w;
  double vertRes, horizon, alphaMin, psiMin, dMin;
  std::vector<double> pixelScale(1, PIXELSCALE);
  
  // path names differ whether integrated into PROG system or not
  char *warpingSIMDSmallTest_loc = getenv("warpingSIMDSmallTest_loc");
  if (warpingSIMDSmallTest_loc)
    sprintf(dataPath, "%s/data/%s", warpingSIMDSmallTest_loc, DATAPATH);
  else
    sprintf(dataPath, "../%s", DATAPATH);
  // build directory names and get bw
  sprintf(dbFile, "%s/%s_%s.db", dataPath, BASE, SUFFIX);
  sprintf(ssPath, "%s/%s%s%s", dataPath, BASE, SS, SUFFIX);
  sprintf(cvPath, "%s/%s%s%s", dataPath, BASE, CV, SUFFIX);
  // read db file: width, grid coordinates etc.
  FILE *f;
  assert((f = fopen(dbFile, "r")));
  assert((fscanf(f, "%d %d %d %d %d %lf %lf", 
		 &w, &xmin, &xmax, &ymin, &ymax, &vertRes, &horizon) == 7));
  fprintf(stderr, 
	  "w = %d, x: %d %d, y: %d %d, vertRes = %g, horizon = %g\n", 
	  w, xmin, xmax, ymin, ymax, vertRes, horizon);
  fclose(f);
  // read images (without random rotation)
  SIMDImage<IMGTYPE,SW,SA,Panorama> ss, cv;
  sprintf(ssFile, "%s/cv_%d_%d_bw%s.pgm", ssPath, SSX, SSY, BW);
  sprintf(cvFile, "%s/cv_%d_%d_bw%s.pgm", cvPath, CVX, CVY, BW);
  fprintf(stderr, "ssFile = %s\n", ssFile);
  fprintf(stderr, "cvFile = %s\n", cvFile);
  assert(loadPGM(ssFile, ss, Panorama(vertRes, horizon)));
  assert(loadPGM(cvFile, cv, Panorama(vertRes, horizon)));
  // save images
  assert(savePGM("ss.pgm", ss));
  assert(savePGM("cv.pgm", cv));
  // create class
  WarpingExample<IMGTYPE,PROCTYPE,MEASTYPE,SPSTYPE,MATCHTYPE,COMPASSTYPE,
		 SW,SA> we(w, NALPHA, NPSI);
  // run warping
  // 27. Feb 18 (rm): vert.res. and hor. taken from images
  we.run(ss, cv, 
	 pixelScale, POSTSCALE, 
	 alphaMin, psiMin, dMin);
  // home vector angle (relative to grid x axis!)
  double betaMin = -alphaMin + psiMin + M_PI;
  // true home vector (in grid coordinates!)
  double xTrueHome = SSX - CVX, yTrueHome = SSY - CVY;
  // true home vector angle (relative to grid x axis!)
  double betaTrue = atan2(yTrueHome, xTrueHome);
  // true psi (0.0 since images are not rotated)
  double psiTrue = 0.0;
  // print result
  printf("alphaMin = %g, psiMin = %g, dMin = %g\n",
	 alphaMin, psiMin, dMin);
  printf("betaMin = %g, betaTrue = %g, betaAE = %g\n",
	 betaMin, betaTrue, acos(cos(betaMin - betaTrue)));
  printf("psiMin = %g, psiTrue = %g, psiAE = %g\n",
	 psiMin, psiTrue, acos(cos(psiMin - psiTrue)));
  return 0;
}
