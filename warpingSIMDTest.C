// ===========================================================================
// 
// warpingSIMDTest.C --
// test code for WarpingCompound class
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

// test

// un-comment these to activate debug output
// #define SHARED_SIMD_PTR_DEBUG_ON
// #define SIMD_IMAGE_DEBUG_ON
// #define SHARED_SIMD_PTR_HEAP_COUNT

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
#include "WarpingFactories.H"
#include "WarpingBundle.H"
#include "SIMDExperimental.H"

// ===========================================================================
// simple array class from std::vector
// ===========================================================================

template <class TYPE>
class Array : public std::vector< std::vector<TYPE> >
{ 
public:
  Array(int nx, int ny): 
    std::vector< std::vector<TYPE> >(nx, std::vector<TYPE>(ny)) {}
};

using namespace ns_simd;

// ===========================================================================
// random image
// ===========================================================================

template <class TYPE>
void randomImage(TYPE *data, size_t size)
{
    for (unsigned int i = 0; i < size; ++i)
      data[i] = (TYPE) (drand48() * 255);
}

// ===========================================================================
// max function on std::vector<double>
// ===========================================================================

std::vector<double>
maxVec(const std::vector<double> &x,
       const std::vector<double> &y)
{
  assert(x.size() == y.size());
  std::vector<double> z(x.size());
  for (size_t i = 0; i < z.size(); i++)
    z[i] = std::max(x[i], y[i]);
  return z;
}

// ===========================================================================
// heap checks
// ===========================================================================

#ifdef SHARED_SIMD_PTR_HEAP_COUNT
SHARED_SIMD_PTR_HEAP_COUNT_INIT
void heapInfo(FILE *f)
{
  fprintf(f, "heap %d desc %lu mem %lu\n", 
	  mallinfo().uordblks, 
	  SharedSIMDPtrDesc::totalDescNo,
	  SharedSIMDPtrBase::totalSize);
}
#define HEAP_INFO(TXT) fprintf(stderr, TXT " "); heapInfo(stderr);
#else
#define HEAP_INFO(TXT)
#endif

// ===========================================================================
// types used for Warping processing
// ===========================================================================

// types
//
// restrictions explored by Benedikt Volkmer:
// IMGTYPE: arbitrary
// PROCTYPE:  only explored for SIMDFloat (req'd by all Harris versions) 
// MEASTYPE:  arbitrary signed type (SIMDSignedByte,SIMDShort,SIMDInt,SIMDFloat)
// SPSTYPE:   requires sizeof(SPSTYPE) <= sizeof(MEASTYPE)
// COMPASSTYPE: SIMDShort,SIMDInt,SIMDFloat
// MATCHTYPE: depends on SPSTYPE:
//       SPSTYPE            MATCHTYPE
// ------------------- -----------------------
//         Byte         Word,Short,Int,Float
//         Word            Word,Int,Float
//   SignedByte/Short     Short,Int,Float
//       Int,Float          Int,Float
//
// TODO: explore reason for relation between MEASTYPE and SPSTYPE:
// TODO: otherwise images of MEASTYPE can't be swizzled to to match wSPS
//
// TODO: PROCTYPE can presently only be SIMDFloat:
// TODO: req'd by Harris versions in WarpingFactories (doesn't compile)
// TODO: problem: PROCTYPE is template parameter of WarpingSPSComputation
// TODO: (base class), has to be set
// TODO: one could rewrite the Harris versions such that internally SIMDFloat
// TODO: is used and convert from PROCTYPE to SIMDFloat and from SIMDFloat
// TODO: to PROCTYPE, but then the problem is, that the values could get
// TODO: too large in these versions, so one would need a "preprocScale"
// TODO: (comment by Benedikt Volkmer)
//
#define IMGTYPE SIMDByte
#define PROCTYPE SIMDFloat
#define MEASTYPE SIMDShort
#define SPSTYPE SIMDByte
#define MATCHTYPE SIMDWord
// COMPASSTYPE needs to be signed (uses hadds which only exists for signed)
// for MATCHTYPE == SIMDWord, we can't use SIMDShort since the pixelScale is
// too large, therefore we use COMPASSTYPE = SIMDInt
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

// ----- informative output, not results (these go to stdout) -----
FILE *finfo = stderr;

// ----- parameter handling -----
// parameters can be changed by setenv (tcsh) or export (bash)

// should parameters be printed?
#define PRINT_NO_PARAMETERS 0
#define PRINT_OVERWRITTEN_PARAMETERS 1
#define PRINT_ALL_PARAMETERS 2
int printParameters = PRINT_ALL_PARAMETERS;

int 
getEnvInt(const char *name, const int defaultVal)
{
  char *valstr = getenv(name);
  int val = valstr ? atoi(valstr) : defaultVal;
  if (((printParameters == PRINT_OVERWRITTEN_PARAMETERS) && valstr) ||
      (printParameters == PRINT_ALL_PARAMETERS))
    fprintf(finfo, "PARAMETER int %s = %d (%d)\n", 
	    name, val, defaultVal);
  return val;
}

double
getEnvDouble(const char *name, const double defaultVal)
{
  char *valstr = getenv(name);
  double val = valstr ? atof(valstr) : defaultVal;
  if (((printParameters == PRINT_OVERWRITTEN_PARAMETERS) && valstr) ||
      (printParameters == PRINT_ALL_PARAMETERS))
    fprintf(finfo, "PARAMETER double %s = %g (%g)\n", 
	    name, val, defaultVal);
  return val;
}

char*
getEnvString(const char *name, const char *defaultVal)
{
  char *valstr = getenv(name);
  char *str = new char[256];
  if (valstr) strcpy(str, valstr); else strcpy(str, defaultVal);
  if (((printParameters == PRINT_OVERWRITTEN_PARAMETERS) && valstr) ||
      (printParameters == PRINT_ALL_PARAMETERS))
    fprintf(finfo, "PARAMETER char* %s = \"%s\" (\"%s\")\n", 
	    name, str, defaultVal);
  return str;
}

std::vector<double>
getEnvDoubleVector(const char *name, const double defaultVal[],
                   const unsigned int defaultSize)
{
  char *valstr = getenv(name);
  std::vector<double> result;
  if (!valstr) {
    for (unsigned int i = 0; i < defaultSize; i++) {
      result.push_back(defaultVal[i]);
    }
  }
  else {
    int size = 1;
    for (char *p=valstr; *p; p++)
      if (*p == ',')
        size++;
    result = std::vector<double>(size);
    int pos = 0;
    char *ptr = strtok(valstr, ",");
    while (ptr) {
      result[pos++] = atof(ptr);
      ptr = strtok(NULL, ",");
    }
  }
  fprintf(finfo, "PARAMETER std::vector<double> %s = {", name);
  for (unsigned int i = 0; i <  result.size(); i++)
      fprintf(finfo, "%g,", result[i]);
  fputs("} ({", finfo);
  for (unsigned int i = 0; i < defaultSize; i++)
      fprintf(finfo, "%g,", defaultVal[i]);
  fputs("})\n", finfo);
  return result;
}

#define INT(NAME,DEFAULTVAL) int NAME = getEnvInt(#NAME, DEFAULTVAL)
#define DOUBLE(NAME,DEFAULTVAL) double NAME = getEnvDouble(#NAME, DEFAULTVAL)
#define STRING(NAME,DEFAULTVAL) char *NAME = getEnvString(#NAME, DEFAULTVAL)
#define DOUBLE_VEC(NAME,DEFAULTVAL, DEFAULTSIZE)	\
  std::vector<double> NAME =				\
    getEnvDoubleVector(#NAME, DEFAULTVAL, DEFAULTSIZE)

// ----- path names -----
// path to data files generated by this program
STRING(testPathName, "TEST");
// path to image databases used by this program
STRING(dataPathName, "DATA2");

// use some test code
// INT(useTestCode, 0);

STRING(expId, "");

// rotate images
INT(rotateImages, 1);

// write rotated images to output
INT(writeDatabaseImages, 0);

// ----- parameter for constructor -----
// number of search steps in alpha and phi direction
// note that there are restrictions on the choice depending on image width
INT(nAlpha, 96);
INT(nPsi, 96);
// number of scale planes
INT(nScalePlanes, 9);
// maximal scale factor
DOUBLE(maxScaleFactor, 2.0);
// maximal threshold
DOUBLE(maxThresholdMinWarping, 2.5);
// restriction of min-warping template to certain values of rho;
// (0.0, 100.0) means that there is no restriction
DOUBLE(rhoMinMinWarping, 0.0);
DOUBLE(rhoMaxMinWarping, 100.0);
DOUBLE(rhoMinWarping, 0.0);
DOUBLE(rhoMaxWarping, 1.0);
INT(nRhoWarping, 20);

// ----- parameter for first phase -----
INT(firstPhase, 
    WarpingSPSComputationSelector::firstPhaseNSAD2);

// ............................................................................
// parameter for distance measures
// interpolation:   interpolation used in image magnification
//                  0: nearest neighbor, 1: linear
// pixelScale:      scaling of pixel values when input images are converted
//                  into into internal representation,
//                  use idealPixelScale function to determine optimal value
// postScale:       scaling of computed image distance values,
//                  use idealPostScale function to determine optimal value
// ............................................................................
INT(interpolation, 0);
// 16 for w = 288, same for larger height with same bww (for edge measures)
const double pixelScaleDefault[] = {16};
DOUBLE_VEC(pixelScale, pixelScaleDefault, 1);
// DOUBLE(pixelScale, 16);
// ASC 100, NSAD 200 for w = 288, adjust for other w! 
DOUBLE(postScale, 200);
// Harris Parameters (contributed by Benedikt Volkmer)
DOUBLE(harrisK, 0.04);
INT(harrisBinomialFilterCount, 1);
// Sigmoid weights
DOUBLE_VEC(sigmoidW, NULL, 0);
DOUBLE_VEC(sigmoidW0, NULL, 0);

// ----- parameter for second phase -----
// const int minWarpingSearcher = 0;
// const int warpingSearcher = 1;
INT(searcher,
    WarpingSearcherSelector::minWarpingSearcher);
// XPAY: search loop order x, psi, alpha, y (typically fastest)
INT(searchMethodMinWarping,
    MinWarpingSearchSelector::searchMethodMinWarpingXPAY);
// XPAY: search loop order x, psi, alpha, y (typically fastest)
INT(searchMethodWarping, 
    WarpingSearchSelector::searchMethodWarpingPARX);
// switch for compass acceleration (1=sum, 2=quantil)
INT(compassAcceleration, 0);
// quantil fraction for quantil-based compass acceleration
DOUBLE(quantilCompassFraction, 0.2);
// parameter for compass acceleration
DOUBLE(psiFraction, 0.3);
// switch for partial search (otherwise full search)
INT(partialSearch, 0);
// partial searcher, but full search range (tests only)
INT(partialFullSearch,0);
// search radii for partial search
DOUBLE(alphaRad, M_PI / 4);
DOUBLE(psiRad, M_PI / 4);
// switch for double search (second run with exchanged snapshot/current view)
// affects both partial and full search
INT(doubleSearch, 1);

// ----- parameter for third phase -----
// fine search is executed
INT(fineSearch, 1);
// fine search uses double search
INT(doubleSearchFine, 1);
// fine search only uses immediate surroundings
INT(fineSearchSurround, 0);

// ===========================================================================
// main
// ===========================================================================

int
main(int argc, char *argv[])
{
  char id[256], dbFile[261], ssPath[257], cvPath[257], imageFile[290], 
    ssFile[547], cvFile[547], bw[16], action, dataPath[256], testPath[256];
  int xss, yss, xcv, ycv, xmin, xmax, ymin, ymax, w, nx, ny;
  double vertRes, horizon;
  FILE *f;

  printf("SW = %d\n", SW);
  printf("__cplusplus = %ld\n", (long int) __cplusplus);
  // initialize random number generator (fixed seed) 3513581345
  srand48(351358);
  // check if minimal number of arguments is available
  if (argc < 2) {
    fprintf(stderr, 
	    "warpingSIMDTest <action> [<base> <ss> <cv> <suffix> <bw>] "
	    "{<option>}*\n");
    exit(-1);
  }
  action = argv[1][0];
  // except of action x and Z
  if ((action != 'x') && (action != 'Z')) {
    // all actions need the remaining 5 parameters (<base> ... <bw>)
    if (argc < 7) {
      fprintf(stderr,
	      "warpingSIMDTest <action> <base> <ss> <cv> <suffix> <bw> "
	      "{<option>}*\n");
      exit(-1);
    }
    // path names differ whether integrated into PROG system or not
    char *warpingSIMDTest_loc = getenv("warpingSIMDTest_loc");
    sprintf(testPath, "./%s", testPathName);
    if (warpingSIMDTest_loc)
      sprintf(dataPath, "%s/data/%s", warpingSIMDTest_loc, dataPathName);
    else
      sprintf(dataPath, "../%s", dataPathName);
    fprintf(finfo, "using database at %s\n", dataPath);
    fprintf(finfo, "saving test data to %s\n", testPath);
    // build directory names and get bw
    sprintf(id, "%s_%s_%s_%s_%s%s", argv[2], argv[3], argv[4], argv[5], argv[6], expId);
    fprintf(finfo, "id = %s\n", id);
    sprintf(dbFile, "%s/%s_%s.db", dataPath, argv[2], argv[5]);
    fprintf(finfo, "dbFile = %s\n", dbFile);
    sprintf(ssPath, "%s/%s%s%s", dataPath, argv[2], argv[3], argv[5]);
    fprintf(finfo, "ssPath = %s\n", ssPath);
    sprintf(cvPath, "%s/%s%s%s", dataPath, argv[2], argv[4], argv[5]);
    fprintf(finfo, "cvPath = %s\n", cvPath);
    strcpy(bw, argv[6]);
    fprintf(finfo, "bww = %s\n", bw);
    // read db file: width, grid coordinates etc.
    assert((f = fopen(dbFile, "r")));
    assert((fscanf(f, "%d %d %d %d %d %lf %lf", 
		   &w, &xmin, &xmax, &ymin, &ymax, &vertRes, &horizon) == 7));
    fprintf(finfo, 
	    "DBINFO: w = %d, x: %d %d, y: %d %d, vertRes = %g, horizon = %g\n", 
	    w, xmin, xmax, ymin, ymax, vertRes, horizon);
    fclose(f);
    // create array of images and load images
    nx = (xmax - xmin) + 1;
    ny = (ymax - ymin) + 1;
    fprintf(finfo, "nx = %d, ny = %d\n", nx, ny);
    // SIMDImage<IMGTYPE,SW,SA,Panorama> ss[nx][ny], cv[nx][ny];
    // double rotAngle[nx][ny];
    HEAP_INFO("before ss, cv creation");
    Array< SIMDImage<IMGTYPE,SW,SA,Panorama> > ss(nx, ny), cv(nx, ny);
    HEAP_INFO("after ss, cv creation");
    // array contains rotations in mathematically negative angles
    Array< double > rotAngle(nx, ny);
    HEAP_INFO("after rotAngle creation");   
    SIMDImage<IMGTYPE,SW,SA,Panorama> ssOrig, cvOrig;
    HEAP_INFO("after ssOrig, cvOrig creation");   
    for (int x = 0; x < nx; x++) {
      // fprintf(finfo, "%2d: ", x); 
      for (int y = 0; y < ny; y++) {
	// fprintf(finfo, "%2d,", y);
	// read image file
	sprintf(imageFile, "cv_%d_%d_bw%s.pgm", x+xmin, y+ymin, bw);
	// fprintf(finfo, "reading %s\n", imageFile);
	sprintf(ssFile, "%s/%s", ssPath, imageFile);
	sprintf(cvFile, "%s/%s", cvPath, imageFile);
	// puts(ssFile);
	assert(loadPGM(ssFile, ssOrig, Panorama(vertRes, horizon)));
	// puts(cvFile);
	assert(loadPGM(cvFile, cvOrig, Panorama(vertRes, horizon)));
	if (rotateImages) {
	  // randomly rotate images to make task more demanding for min-warping
	  int rotIndex = (int) (drand48() * w);
	  rotAngle[x][y] = rotIndex * (2.0 * M_PI) / w;
	  // we use the same rotAngle field for both databases
	  rotateHor(ssOrig, rotIndex, ss[x][y]);
	  rotateHor(cvOrig, rotIndex, cv[x][y]);
	}
	else {
	  // don't rotate images
	  rotAngle[x][y] = 0.0;
	  ss[x][y] = ssOrig;
	  cv[x][y] = cvOrig;
	}
	// fprintf(finfo, "\n");
      }
    }
    HEAP_INFO("after image loading");   
    // write files (test only)
    if (writeDatabaseImages)
      for (int x = 0; x < nx; x++) {
	for (int y = 0; y < ny; y++) {
	  sprintf(imageFile, "%s/test_ss_%d_%d.pgm", testPath, x, y);
	  assert(savePGM(imageFile, ss[x][y]));
	  sprintf(imageFile, "%s/test_cv_%d_%d.pgm", testPath, x, y);
	  assert(savePGM(imageFile, cv[x][y]));
	}
      }
    
    // =========================================================================
    // create warping
    // =========================================================================
    
    fprintf(finfo, "creating warping bundle\n");
    WarpingBundle<IMGTYPE,PROCTYPE,
		  MEASTYPE,SPSTYPE,
		  MATCHTYPE,COMPASSTYPE,
		  SW,SA>
      warping(w, nAlpha, nPsi,
	      nScalePlanes, maxScaleFactor, maxThresholdMinWarping,
	      rhoMinMinWarping, rhoMaxMinWarping,
	      rhoMinWarping, rhoMaxWarping, nRhoWarping,
	      firstPhase,
	      searchMethodWarping,
	      searchMethodMinWarping,
	      HarrisParameter(harrisK, harrisBinomialFilterCount),
	      sigmoidW, sigmoidW0);
    
    // =========================================================================
    // run different commands
    // =========================================================================

    switch (action) {
      
    case 'p':
      // -----------------------------------------------------------------------
      // p: compute parameters 
      // -----------------------------------------------------------------------
      // example: ./warpingSIMDTest p living3 day night Hh288sh 0.10   
      {
	if (argc != 7) {
	  fprintf(stderr, "excess parameters given\n");
	  exit(-1);
	}
	int w = ss[0][0].w, h = ss[0][0].h;
	fprintf(finfo, "%d x %d\n", w, h);
	const int numPreprocs = warping.spsComp->numPreprocs();
	std::vector<double> maxDenom(numPreprocs, 0.0);
	for (int x = 0; x < nx; x++)
	  for (int y = 0; y < ny; y++) {
	    maxDenom = maxVec(maxDenom, warping.spsComp->maxDenom(ss[x][y]));
	    maxDenom = maxVec(maxDenom, warping.spsComp->maxDenom(cv[x][y]));
	  }
	fprintf(finfo, "maxDenom = ");
        for (int i = 0; i < numPreprocs-1; i++)
          fprintf(finfo, "%g,", maxDenom[i]);
        fprintf(finfo, "%g\n", maxDenom[numPreprocs-1]);
	std::vector<double> localPixelScale;
	localPixelScale = warping.spsComp->idealPixelScale(maxDenom);
	printf("SCALING: %s firstPhase %d pixelScale ",
	       id,
	       firstPhase);
        for (int i = 0; i < numPreprocs-1; i++)
          printf("%g,", localPixelScale[i]);
        printf("%g postScale %g\n",
               localPixelScale[numPreprocs-1],
	       warping.spsComp->idealPostScale
	       (double(SIMDTypeInfo<MATCHTYPE>::max()) / w));
      }
      break;
      
    case 's':
      // -----------------------------------------------------------------------
      // s: single snapshot and all current views
      // -----------------------------------------------------------------------
      // example: ./warpingSIMDTest s living3 day night Hh288sh 0.10 10 4  
      {
	char gnuplotFile[551], gnuplotTrueFile[555];
	if (argc != 9) {
	  fprintf(stderr, "<xss> <yss>\n");
	  exit(-1);
	}
        assert((int)pixelScale.size() == warping.spsComp->numPreprocs());
	int xssExt = atoi(argv[7]), yssExt = atoi(argv[8]) ; 
	// for the filename we use the external ss coordinates
	sprintf(gnuplotFile, "%s/homevectors_%s_%d_%d.gp", 
		testPath, id, xssExt, yssExt);
	sprintf(gnuplotTrueFile, "%s/truehomevectors_%s_%d_%d.gp", 
		testPath, id, xssExt, yssExt);
	// internal snapshot coordinates
	xss = xssExt - xmin; yss = yssExt - ymin;
	fprintf(finfo, "ss internal (%d, %d)\n", xss, yss);
	assert((xss >= 0) && (xss < nx) && (yss >= 0) && (yss < ny));
	FILE *fg, *fgt;
	assert(fg = fopen(gnuplotFile, "w"));
	assert(fgt = fopen(gnuplotTrueFile, "w"));
	fprintf(fg, "%d %d\n\n\n", xss, yss);
	fprintf(fgt, "%d %d\n\n\n", xss, yss);
	// signed long long t1, t2, t1sum = 0LL, t2sum = 0LL;
	double t1, t2, t1sum = 0, t2sum = 0;
	double alphaMin, psiMin, dMin, betaErrorSum = 0.0, psiErrorSum = 0.0;
	int iAlphaMin, iPsiMin, iAlphaMinInv, iPsiMinInv;
	int jAlphaOff, jPsiOff;
	SIMDImage<MATCHTYPE,1,1> fineMatch;
	int n = 0;
	for (int x = 0; x < nx; x++)
	  for (int y = 0; y < ny; y++)
	    if ((x != xss) || (y != yss)) {
	      // printf("%d %d\n", x, y);
	      // note that WarpingSIMD assumes mathematically positive
	      // angles, but the azimuth coordinate in the database images
	      // increases with mathematically negative azimuth angle, this
	      // explains the sign inversions in the pre- and post-angle
	      // computations
	      // true psi
	      double psiTrue = rotAngle[x][y] - rotAngle[xss][yss];
	      // true home vector (in grid coordinates!)
	      double xTrueHome = xss - x, yTrueHome = yss - y;
	      // true home vector angle (relative to grid x axis!)
	      double betaTrue = atan2(yTrueHome, xTrueHome);
	      // true alpha value
	      double alphaTrue = psiTrue - betaTrue + M_PI - rotAngle[x][y];
	      // run warping
	      // we use alphaTrue and psiTrue as estimates for partial
	      // 27. Feb 18 (rm): vert.res. and hor. taken from images
	      warping.run(ss[xss][yss], cv[x][y], 
			  interpolation, pixelScale, postScale,
			  alphaTrue, psiTrue,
			  searcher,
			  compassAcceleration, 
			  psiFraction, quantilCompassFraction,
			  partialSearch,
			  alphaRad, psiRad,
			  partialFullSearch,
			  doubleSearch,
			  fineSearch,
			  doubleSearchFine,
			  fineSearchSurround,
			  alphaMin, psiMin, dMin, 
			  iAlphaMin, iPsiMin, 
			  iAlphaMinInv, iPsiMinInv,
			  jAlphaOff, jPsiOff, fineMatch,
			  t1, t2);
	      // home vector angle (relative to grid x axis!)
	      double betaMin = -alphaMin + psiMin + M_PI - rotAngle[x][y];
	      // errors in [0,pi]
	      psiErrorSum += acos(cos(psiMin - psiTrue));
	      betaErrorSum += acos(cos(betaMin - betaTrue));
	      // data for gnuplot files
	      fprintf(fg, "%d %d %g %g\n", 
		      x, y, cos(betaMin), sin(betaMin));
	      fprintf(fgt, "%d %d %g %g\n",
		      x, y, cos(betaTrue), sin(betaTrue));
	      // sum up computation times
	      t1sum += t1;
	      t2sum += t2;
	      n++;
	    }
	fclose(fg);
	fclose(fgt);
	// print average values, times are in micro-secs
	printf("RESULTS: betaAE %g psiAE %g t1 %.0f t2 %.0f n %d\n", 
	       betaErrorSum / n, psiErrorSum / n, t1sum / n, t2sum / n, n);
      }
      break;
      
    case 'a':
    case 'C':
    case 'S':
      // -----------------------------------------------------------------------
      // a, C, S: all snapshots and all current views
      // -----------------------------------------------------------------------
      // a: ss[xss][yss] <-> cv[x][y]
      // C: ss[xss][yss] <-> cv[x][y] and
      //    cv[xss][yss] <-> ss[x][y]
      // S: ss[xss][yss] <-> ss[x][y] and
      //    cv[xss][yss] <-> cv[x][y]
      // (C and S only make sense if different ss/cv sets are provided)
      //
      // example: 
      // ./warpingSIMDTest a living3 day night Hh288sh 0.10 gridStep boxRad
      {
	int pairs = (action == 'a') ? 1 : 2;
	std::vector<double> betaErrorList, psiErrorList, alphaErrorList;
	// no further cmdline parameters necessary
	if (argc != 9) {
	  fprintf(stderr, "<gridStep> <boxRad>\n");
	  exit(-1);
	}
        assert((int)pixelScale.size() == warping.spsComp->numPreprocs());
	char allFileName[522];
	sprintf(allFileName, "%s/all_%c_%s.dat", testPath, action, id);
	FILE *fAll;
	assert(fAll = fopen(allFileName, "w"));
	fprintf(fAll,
		"# xss yss x y dist betaError psiError alphaError pairing\n");
	int gridStep = atoi(argv[7]);
	int boxRad = atoi(argv[8]);
	fprintf(finfo, "ALL: gridStep = %d, boxRad = %d\n", gridStep, boxRad);
	// signed long long t1, t2, t1sum = 0LL, t2sum = 0LL;
	double t1, t2, t1sum = 0, t2sum = 0;
	double alphaMin, psiMin, dMin, betaError, psiError, alphaError;
	double betaErrorSum = 0.0, psiErrorSum = 0.0, alphaErrorSum = 0.0;
	int iAlphaMin, iPsiMin, iAlphaMinInv, iPsiMinInv;
	int jAlphaOff, jPsiOff;
	SIMDImage<MATCHTYPE,1,1> fineMatch;
	int n = 0;
	for (int xss = 0; xss < nx; xss += gridStep) {
	  // box with radius boxRad around xss, but inside database
	  int xl = std::max(xss - boxRad, 0);
	  int xr = std::min(xss + boxRad, nx-1);
	  for (int yss = 0; yss < ny; yss += gridStep) {
	    // box with radius boxRad around yss, but inside database
	    int yt = std::max(yss - boxRad, 0);
	    int yb = std::min(yss + boxRad, ny-1);
	    for (int x = xl; x <= xr; x += gridStep)
	      for (int y = yt; y <= yb; y += gridStep)
		if ((x != xss) || (y != yss)) {
		  // note that WarpingSIMD assumes mathematically positive
		  // angles, but the azimuth coordinate in the database images
		  // increases with mathematically negative azimuth angle, this
		  // explains the sign inversions in the pre- and post-angle
		  // computations
		  // true psi
		  double psiTrue = rotAngle[x][y] - rotAngle[xss][yss];
		  // true home vector (in grid coordinates!)
		  double xTrueHome = xss - x, yTrueHome = yss - y;
		  // true home vector angle (relative to grid x axis!)
		  double betaTrue = atan2(yTrueHome, xTrueHome);
		  // true alpha value
		  double alphaTrue = psiTrue - betaTrue + M_PI - rotAngle[x][y];
		  // different pairings ('a' vs. 'A' mode)
		  SIMDImage<IMGTYPE,SW,SA,Panorama> _ss_, _cv_;
		  for (int p = 0; p < pairs; p++) {
		    switch (action) {
		    case 'a':
		      // only p == 0 possible (pairs == 1)
		      _ss_ = ss[xss][yss];
		      _cv_ = cv[x][y];
		      break;
		    case 'C':
		      // both cross-db tests
		      _ss_ = (p == 0) ? ss[xss][yss] : cv[xss][yss];
		      _cv_ = (p == 0) ? cv[x][y] :     ss[x][y];
		      break;
		    case 'S':
		      // both same-db tests
		      _ss_ = (p == 0) ? ss[xss][yss] : cv[xss][yss];
		      _cv_ = (p == 0) ? ss[x][y] :     cv[x][y];
		      break;
		    }
		    // run warping
		    // we use alphaTrue and psiTrue as estimates for partial
		    // 27. Feb 18 (rm): vert.res. and hor. taken from images
		    warping.run(_ss_, _cv_, 
				interpolation, pixelScale, postScale,
				alphaTrue, psiTrue,
				searcher,
				compassAcceleration,
				psiFraction, quantilCompassFraction,
				partialSearch,
				alphaRad, psiRad,
				partialFullSearch,
				doubleSearch,
				fineSearch,
				doubleSearchFine,
				fineSearchSurround,
				alphaMin, psiMin, dMin,
				iAlphaMin, iPsiMin, 
				iAlphaMinInv, iPsiMinInv,
				jAlphaOff, jPsiOff, fineMatch,
				t1, t2);
		    // home vector angle (relative to grid x axis!)
		    double betaMin = -alphaMin + psiMin + M_PI - rotAngle[x][y];
		    // errors in [0,pi]
		    psiError = acos(cos(psiMin - psiTrue));
		    betaError = acos(cos(betaMin - betaTrue));
		    alphaError = acos(cos(alphaMin - alphaTrue));
		    psiErrorSum += psiError;
		    betaErrorSum += betaError;
		    alphaErrorSum += alphaError;
#ifndef SHARED_SIMD_PTR_HEAP_COUNT
		    // we omit this when we check the heap size
		    psiErrorList.push_back(psiError);
		    betaErrorList.push_back(betaError);
		    alphaErrorList.push_back(alphaError);
#endif
		    // sum up computation times
		    t1sum += t1;
		    t2sum += t2;
		    n++;
		    // write to file
		    double d = hypot(xss - x, yss - y);
		    fprintf(fAll, "%d %d %d %d %g %g %g %g %d\n",
			    xss, yss, x, y, d,
			    betaError, psiError, alphaError, p);
		    HEAP_INFO("after warping run");
		  }
		}
	  }
	}
	fclose(fAll);
#ifndef SHARED_SIMD_PTR_HEAP_COUNT
	// we have to omit this since the lists are empty when checking heap
	int underMedianIndex = (n - 1) / 2;
	nth_element(psiErrorList.begin(),
		    psiErrorList.begin() + underMedianIndex,
		    psiErrorList.end());
	nth_element(betaErrorList.begin(),
		    betaErrorList.begin() + underMedianIndex,
		    betaErrorList.end());
	nth_element(alphaErrorList.begin(),
		    alphaErrorList.begin() + underMedianIndex,
		    alphaErrorList.end());
	// print average values, times are in micro-secs
	printf("RESULTS: betaAE %g %g psiAE %g %g alphaAE %g %g "
	       "t1 %.0f t2 %.0f n %d\n", 
	       betaErrorSum / n, betaErrorList[underMedianIndex],
	       psiErrorSum / n, psiErrorList[underMedianIndex],
	       alphaErrorSum / n, alphaErrorList[underMedianIndex],
	       t1sum / n, t2sum / n, n);
#endif
      }
      break;
     
    case 'v':
      // -----------------------------------------------------------------------
      // v: single snapshot, single current view, visualization
      // -----------------------------------------------------------------------
      // example: ./warpingSIMDTest v living3 day night Hh288sh 0.10 9 4 5 2  
      {
	char fnSuffix[304];
	
	if (argc != 11) {
	  fprintf(stderr, "<xss> <yss> <xcv> <ycv>\n");
	  exit(-1);
	}
        assert((int)pixelScale.size() == warping.spsComp->numPreprocs());
	int xssExt = atoi(argv[7]), yssExt = atoi(argv[8]);
	int xcvExt = atoi(argv[9]), ycvExt = atoi(argv[10]);
	xss = xssExt - xmin; yss = yssExt - ymin;
	xcv = xcvExt - xmin; ycv = ycvExt - ymin;
	sprintf(fnSuffix, "%s_%d_%d_%d_%d", 
		id, xssExt, yssExt, xcvExt, ycvExt);
	fprintf(finfo, "ss internal (%d, %d)\n", xss, yss);
	fprintf(finfo, "cv internal (%d, %d)\n", xcv, ycv);
	assert((xss >= 0) && (xss < nx) && (yss >= 0) && (yss < ny));
	assert((xcv >= 0) && (xcv < nx) && (ycv >= 0) && (ycv < ny));
	// print scale factors and thresholds,
	// scaleFactors are chosen inversion-symmetric and with the same
	// factor between neighboring scale factors,
	// thresholds are chosen in the middle between two scale factors
	printf("symScaleFac: ");
	for (unsigned i = 0;
	     i < warping.wsc->symmScaleFac.scaleFactors.size(); i++)
	  printf("%6.3f ", warping.wsc->symmScaleFac.scaleFactors[i]);
	printf("\n");
	printf("thresholds:   ");
	for (unsigned i = 0; 
	     i < warping.wsc->minWarping.searchTemplate->thresholds.size(); i++)
	  printf("%6.3f ", 
		 warping.wsc->minWarping.searchTemplate->thresholds[i]);
	printf("\n");
	// for single scale plane
	printf("singleScaleFac: ");
	for (unsigned i = 0;
	     i < warping.wsc->singleScaleFac.scaleFactors.size(); i++)
	  printf("%6.3f ", warping.wsc->singleScaleFac.scaleFactors[i]);
	printf("\n");
	// compute ground truth
	double alphaMin, psiMin, dMin;
	// true psi
	double psiTrue = rotAngle[xcv][ycv] - rotAngle[xss][yss];
	// true home vector (in grid coordinates!)
	double xTrueHome = xss - xcv, yTrueHome = yss - ycv;
	// true home vector angle (relative to grid x axis!)
	double betaTrue = atan2(yTrueHome, xTrueHome);
	// true alpha value
	double alphaTrue = psiTrue - betaTrue + M_PI - rotAngle[xcv][ycv];
	// run warping
	// we use alphaTrue and psiTrue as estimates for partial
	// signed long long t1, t2;
	double t1, t2;
	int iAlphaMin, iPsiMin, iAlphaMinInv, iPsiMinInv;
	int jAlphaOff, jPsiOff;
	SIMDImage<MATCHTYPE,1,1> fineMatch;
	// we use alphaTrue and psiTrue as estimates for partial
	// 27. Feb 18 (rm): vert.res. and hor. taken from images
	warping.run(ss[xss][yss], cv[xcv][ycv], 
		    interpolation, pixelScale, postScale,
		    alphaTrue, psiTrue,
		    searcher,
		    compassAcceleration,
		    psiFraction, quantilCompassFraction,
		    partialSearch,
		    alphaRad, psiRad,
		    partialFullSearch,
		    doubleSearch,
		    fineSearch,
		    doubleSearchFine,
		    fineSearchSurround,
		    alphaMin, psiMin, dMin,
		    iAlphaMin, iPsiMin, 
		    iAlphaMinInv, iPsiMinInv,
		    jAlphaOff, jPsiOff, fineMatch,
		    t1, t2);
#if 1
	// save min-template as readable data file
	char templateFile[577];
	sprintf(templateFile, "%s/mintemplate_%s.dat", testPath, fnSuffix);
	FILE *ftmpl;
	assert(ftmpl = fopen(templateFile, "w"));
	warping.wsc->minWarping.searchTemplate->save(ftmpl);
	fclose(ftmpl);
	// save warping-template as readable data file
	sprintf(templateFile, "%s/warptemplate_%s.dat", testPath, fnSuffix);
        assert(ftmpl = fopen(templateFile, "w"));
	warping.wsc->warping.searchTemplate->save(ftmpl);
	fclose(ftmpl);
	// save min-template as image
	char templateImageFile[576];
	sprintf(templateImageFile, "%s/mintemplate_%s.pgm", testPath, fnSuffix);
	SIMDImage<MATCHTYPE,SW,SA> minTemplate;
	warping.wsc->minWarping.searchTemplate->getImage
	  (minTemplate, 
	   MATCHTYPE(255), -20.0, 200.0);
	assert(savePGM(templateImageFile, minTemplate));
	// save input images (no scaling, requires SIMDByte format)
	sprintf(ssFile, "%s/ss_%d_%d.pgm", testPath, xss, yss);
	assert(savePGM(ssFile, ss[xss][yss]));
	sprintf(cvFile, "%s/cv_%d_%d.pgm", testPath, xcv, ycv);
	assert(savePGM(cvFile, cv[xcv][ycv]));
	// result
	printf("SOLUTION: alphaMin = %g, psiMin = %g, dMin = %g\n", 
	       alphaMin, psiMin, dMin);
	printf("INDICES: iAlphaMin = %d, iPsiMin = %d\n", 
	       iAlphaMin, iPsiMin);
	fprintf(finfo, "rotAngle = ss %g cv %g\n", 
		rotAngle[xss][yss], rotAngle[xcv][ycv]);
	// home vector angle (relative to grid x axis!)
	double betaMin = -alphaMin + psiMin + M_PI - rotAngle[xcv][ycv];
	// errors in [0,pi]
	double psiError = acos(cos(psiMin - psiTrue));
	double betaError = acos(cos(betaMin - betaTrue));
	double alphaError = acos(cos(alphaMin - alphaTrue));
	printf("EST: betaMin = %g, psiMin = %g, alphaMin = %g\n", 
	       betaMin, psiMin, alphaMin);
	printf("TRUE: betaTrue = %g, psiTrue = %g, alphaTrue = %g\n", 
	       betaTrue, psiTrue, alphaTrue);
	printf("ERROR: betaError = %g, psiError = %g, alphaError = %g\n", 
	       betaError, psiError, alphaError);
	// solution and inverted solution
	int iAlphaMinVec[2], iPsiMinVec[2];
	iAlphaMinVec[0] = iAlphaMin;
	iPsiMinVec[0] = iPsiMin;
	iAlphaMinVec[1] = iAlphaMinInv;
	iPsiMinVec[1] = iPsiMinInv;
	fprintf(finfo, "MIN: (iAlphaMin, iPsiMin) = (%d, %d) | (%d, %d)\n",
		iAlphaMinVec[0], iPsiMinVec[0],
		iAlphaMinVec[1], iPsiMinVec[1]);
	// spsArray of current searcher method
	WarpingSPS<SPSTYPE,SW,SA> **spsArray 
	  = warping.wsc->getSPSArray(searcher);
	// planes
	for (int spsNo = warping.wsc->SPS_ORIG; 
	     spsNo <= warping.wsc->SPS_INV; spsNo++) {
	  fprintf(finfo, "processing sps %d\n", spsNo);
	  // save SPS in binary format
	  char spsFile[570];
	  sprintf(spsFile, "%s/sps_%d_%s.dat",
		  testPath, spsNo, fnSuffix);
	  // fprintf(finfo, "writing entire sps %d to %s\n", spsNo, spsFile);
	  FILE *fsps;
	  assert(fsps = fopen(spsFile, "w"));
	  spsArray[spsNo]->saveAll(fsps);
	  fclose(fsps);
	  // get solution
	  fprintf(finfo, "getting solution\n");
	  std::vector<int> jDeltaVec, planeIndexVec, jYVec, jXVec, jThetaVec;
	  MATCHTYPE matchSum;
	  switch (searcher) {
	  case WarpingSearcherSelector::minWarpingSearcher:
	    {
	      std::vector<SPSTYPE> minVec; // currently unused
	      minWarpingCurve(*(warping.wsc->minWarping.searchTemplate), 
			      *(warping.wsc->minWarping.spsArray[spsNo]),
			      iAlphaMinVec[spsNo], iPsiMinVec[spsNo], 
			      jDeltaVec, planeIndexVec, minVec, jYVec,
			      jXVec, jThetaVec,
			      matchSum);
	    }
	    break;
	  case WarpingSearcherSelector::warpingSearcher:
	    {
	      int iRhoMin;
	      warpingCurve(*(warping.wsc->warping.searchTemplate), 
			   *(warping.wsc->warping.spsArray[spsNo]),
			   iAlphaMinVec[spsNo], iPsiMinVec[spsNo],
			   jDeltaVec, planeIndexVec, iRhoMin,
			   matchSum);
	      fprintf(finfo, "iRhoMin = %d\n", iRhoMin);
	    }
	    // BUGFIX: 28. Oct 18 (rm) added break
	    break;
	  case WarpingSearcherSelector::minWarpingSearcher1SP:
	    {
	      std::vector<SPSTYPE> minVec; // currently unused
	      minWarpingCurve(*(warping.wsc->minWarping1SP.searchTemplate), 
			      *(warping.wsc->minWarping1SP.spsArray[spsNo]),
			      iAlphaMinVec[spsNo], iPsiMinVec[spsNo], 
			      jDeltaVec, planeIndexVec, minVec, jYVec,
			      jXVec, jThetaVec,
			      matchSum);
	    }
	    break;
	  case WarpingSearcherSelector::warpingSearcher1SP:
	    {
	      int iRhoMin;
	      warpingCurve(*(warping.wsc->warping1SP.searchTemplate), 
			   *(warping.wsc->warping1SP.spsArray[spsNo]),
			   iAlphaMinVec[spsNo], iPsiMinVec[spsNo],
			   jDeltaVec, planeIndexVec, iRhoMin,
			   matchSum);
	      fprintf(finfo, "iRhoMin = %d\n", iRhoMin);
	    }
	    break;
	  default:
	    fprintf(stderr, "invalid searcher %d\n", searcher);
	    exit(-1);
	  }
	  fprintf(finfo, "warpCurve: match = %s\n",
		  SIMDDecimal<MATCHTYPE>(matchSum).str);
	  char solutionFile[575];
	  sprintf(solutionFile, "%s/solution_%d_%s.dat",
		  testPath, spsNo, fnSuffix);
	  // fprintf(finfo, "writing solution of sps %d to %s\n", 
	  //         spsNo, solutionFile);
	  FILE *fsol;
	  assert(fsol = fopen(solutionFile, "w"));
	  for (int jTheta = 0; jTheta < w; jTheta++)
	    fprintf(fsol, "%3d %3d %3d\n", 
		    jTheta, jDeltaVec[jTheta], planeIndexVec[jTheta]);
	  fclose(fsol);
	  SPSTYPE spsMin = min(spsArray[spsNo]->stack());
	  SPSTYPE spsMax = max(spsArray[spsNo]->stack());
	  fprintf(finfo, "sps %d: min %s max %s\n",
	  	  spsNo, 
	  	  SIMDDecimal<SPSTYPE>(spsMin).str, 
	  	  SIMDDecimal<SPSTYPE>(spsMax).str);
	  for (int i = 0; i < spsArray[spsNo]->stack.numPlanes; i++) {
	    SPSTYPE planeMin = min(spsArray[spsNo]->stack[i]);
	    SPSTYPE planeMax = max(spsArray[spsNo]->stack[i]);
	    fprintf(finfo, "sps %d plane %d: min %s max %s\n",
	    	    spsNo, i,
	    	    SIMDDecimal<SPSTYPE>(planeMin).str, 
	    	    SIMDDecimal<SPSTYPE>(planeMax).str);
	    char planeFile[592];
	    char minPlaneFile[584];
	    SIMDImage<SPSTYPE,SW,SA> plane;
	    SIMDImage<SPSTYPE,SW,SA> minPlane;
	    for (int spsFormat = 0; spsFormat <= 2; spsFormat++) {
	      sprintf(planeFile, "%s/plane_%d_%u_format%d_%s.pgm",
		      testPath, spsNo, i, spsFormat, fnSuffix);
	      // fprintf(finfo, "writing plane %d of sps %d to %s\n",
	      //         i, spsNo, planeFile);
	      spsArray[spsNo]->getPlane
		(i, spsFormat, plane, -spsMin, 255.0 / (spsMax - spsMin));
	      assert(savePGM(planeFile, plane));
	      // min-plane
	      sprintf(minPlaneFile, "%s/minplane_%d_format%d_%s.pgm",
		      testPath, spsNo, spsFormat, fnSuffix);
	      // fprintf(finfo, "writing minplane of sps %d to %s\n",
	      // 	 spsNo, minPlaneFile);
	      spsArray[spsNo]->getPlaneMinimum
		(spsFormat, minPlane, -spsMin, 255.0 / (spsMax - spsMin));
	      assert(savePGM(minPlaneFile, minPlane));
	    }	      
	    // visualize solution
	    // only unshuffled!:
	    sprintf(planeFile, "%s/plane_solution_%d_%u_%s.pgm",
		    testPath, spsNo, i, fnSuffix);
	    // fprintf(finfo, 
	    //         "writing plane %d of sps %d with solution to %s\n",
	    // 	    i, spsNo, planeFile);
	    spsArray[spsNo]->getPlaneUnshuffled
	      (i, plane, -spsMin, 255.0 / (spsMax - spsMin));
	    for (int jTheta = 0; jTheta < w; jTheta++)
	      if (planeIndexVec[jTheta] == (int) i)
		if (jDeltaVec[jTheta] >= 0)
		  plane[jDeltaVec[jTheta]][jTheta] 
		    = SIMDTypeInfo<SPSTYPE>::max();
	    assert(savePGM(planeFile, plane));
	    // only unshuffled!
	    sprintf(minPlaneFile, "%s/minplane_solution_%d_%s.pgm",
		    testPath, spsNo, fnSuffix);
	    // fprintf(finfo, "writing minplane of sps %d with solution to %s\n",
	    //		  spsNo, minPlaneFile);
	    spsArray[spsNo]->getPlaneMinimumUnshuffled
	      (minPlane, -spsMin, 255.0 / (spsMax - spsMin));
	    for (int jTheta = 0; jTheta < w; jTheta++)
	      if (jDeltaVec[jTheta] >= 0)
		minPlane[jDeltaVec[jTheta]][jTheta] 
		  = SIMDTypeInfo<SPSTYPE>::max();
	    assert(savePGM(minPlaneFile, minPlane));
	  }
	  fprintf(finfo, "sps %d done\n", spsNo);
	}
	// match
	for (int matchNo = 0; matchNo < warping.wsc->NUM_MATCH; matchNo++) {
	  char matchFile[572];
	  SIMDImage<MATCHTYPE,SW,SA> match;
	  sprintf(matchFile, "%s/match_%d_%s.pgm",
		  testPath, matchNo, fnSuffix);
	  // fprintf(finfo, "writing match %d to %s\n",
	  //	  matchNo, matchFile);
	  MATCHTYPE matchMin, matchMax;
	  warping.wsc->matchArray[matchNo]->extremaExceptInvalid
	    (matchMin, matchMax);
	  fprintf(finfo, "match%d: min %s max %s\n",
		  matchNo,
		  SIMDDecimal<MATCHTYPE>(matchMin).str, 
		  SIMDDecimal<MATCHTYPE>(matchMax).str);
	  warping.wsc->matchArray[matchNo]->getImage
	    (match, -matchMin, 255.0 / (matchMax - matchMin));
	  assert(savePGM(matchFile, match));
	}
	// compass
	SIMDImage<SIMDFloat,SW,SA> compass, quantilCompass;
	for (int spsNo = 0; spsNo < warping.wsc->NUM_SPS; spsNo++) {
	  char compassFile[574];
	  spsArray[spsNo]->getCompassImage<COMPASSTYPE>(compass);
	  spsArray[spsNo]->getQuantilCompassImage
	    (quantilCompass,
	     quantilCompassFraction);
	  sprintf(compassFile, "%s/compass_%d_%s.dat",
		  testPath, spsNo, fnSuffix);
	  // fprintf(finfo, "writing compass data of sps %d to %s\n",
	  //	  spsNo, compassFile);
	  assert(f = fopen(compassFile, "w"));
	  for (size_t i = 0; i < compass.size; i++)
	    fprintf(f, "%g %g\n", compass.data[i], quantilCompass.data[i]);
	  fclose(f);
	}
	// fine match
	if (fineSearch) {
	  char fineMatchFile[574];
	  sprintf(fineMatchFile, "%s/finematch_%s.dat",
		  testPath, fnSuffix);
	  FILE *ffm;
	  assert(ffm = fopen(fineMatchFile, "w"));
	  fprintf(ffm, "iAlphaMin %d jAlphaOff %d iPsiMin %d jPsiOff %d\n",
		  iAlphaMin, jAlphaOff, iPsiMin, jPsiOff);  
	  fprintf(ffm, "alphaMin %g psiMin %g dMin %g\n",
		  alphaMin, psiMin, dMin);
	  fprintf(ffm, 
		  "alpha[iAlphaMin] %g psi[iPsiMin] %g dTheta %g\n",
		  warping.wsc->param.alphaVec[iAlphaMin],
		  warping.wsc->param.psiVec[iPsiMin],
		  warping.wsc->param.dTheta);
	  for (int j = 0; j < fineMatch.h; j++) {
	    for (int i = 0; i < fineMatch.w; i++)
	      fprintf(ffm, "%g ", double(fineMatch[j][i]));
	    fprintf(ffm, "\n");
	  }
	  fclose(ffm);
	}
	// save SPS input images as one diagram
	// here we refer to the SPS compute for symmScaleFac
	MEASTYPE allMax, allMin, allRange;
	SIMDImage<SIMDByte,SW,SA> spsInputDiagram;
	// 27. Feb 18 (rm): vert.res. and hor. taken from images
	warping.spsComp->computeSPSInputDiagram
	  (ss[xss][yss], cv[xcv][ycv], 
	   warping.wsc->symmScaleFac.scaleFactors,
	   interpolation, pixelScale, postScale,
	   *(warping.wsc->spsArray[0]),
	   allMax, allMin, allRange, spsInputDiagram);
	fprintf(finfo,
		"SPS inputs: allMax = %s, allMin = %s, allRange = %s\n",
		SIMDDecimal<MEASTYPE>(allMax).str,
		SIMDDecimal<MEASTYPE>(allMin).str,
		SIMDDecimal<MEASTYPE>(allRange).str);
	char spsInputFile[573];
	sprintf(spsInputFile, "%s/spsinput_%s.pgm",
		testPath, fnSuffix);
	assert(savePGM(spsInputFile, spsInputDiagram));
#endif
      }
      break;

    case 'l':
      // -----------------------------------------------------------------------
      // l: localization
      // -----------------------------------------------------------------------
      // example: ./warpingSIMDTest l living3 day night Hh288sh 0.10 10 4
      // adapted from code by Michael Horst
      {	char gnuplotFile[551];
	if (argc != 9) {
	  fprintf(stderr, "<xss> <yss>\n");
	  exit(-1);
	}
        assert((int)pixelScale.size() == warping.spsComp->numPreprocs());
	int xssExt = atoi(argv[7]), yssExt = atoi(argv[8]); 
	FILE *fg;
	double t1, dMin;
	sprintf(gnuplotFile, "%s/dissim_%s_%d_%d.gp", 
		testPath, id, xssExt, yssExt);
        assert(fg = fopen(gnuplotFile, "w"));
	// internal snapshot coordinates
	xss = xssExt - xmin; yss = yssExt - ymin;
	fprintf(finfo, "ss internal (%d, %d)\n", xss, yss);
	assert((xss >= 0) && (xss < nx) && (yss >= 0) && (yss < ny));
	// print ss position as separate data set
	// 27. Feb 18 (rm): vert.res. and hor. taken from images
	warping.compassMatch(ss[xss][yss], cv[xss][yss], 
			     interpolation, 
			     pixelScale, postScale, 
			     dMin, t1);
	fprintf(fg, "%d %d %f\n\n\n", xss, yss, dMin);
	double t1sum = 0;
	int n = 0;
        for (xcv = 0; xcv < nx; xcv++) {
	  for (ycv = 0; ycv < ny; ycv++) {
	    // 27. Feb 18 (rm): vert.res. and hor. taken from images
	    warping.compassMatch(ss[xss][yss], cv[xcv][ycv], 
				 interpolation, 
				 pixelScale, postScale, 
				 dMin, t1);
	    fprintf(fg, "%d %d %f\n", xcv, ycv, dMin);
	    t1sum += t1;
	    n++;
	  }
	  fprintf(fg, "\n");
        }
        fclose(fg);
	printf("RESULTS: t1 %.0f n %d\n", 
	       t1sum / n, n);
      }
      break;

    case 'z':
      // -----------------------------------------------------------------------
       // z: do nothing except initialization stuff
      // -----------------------------------------------------------------------
      if (argc != 7) {
	fprintf(stderr, "excess parameters given\n");
	exit(-1);
      }
      break;

    case 'Z':
      // -----------------------------------------------------------------------
      // Z: do absolutely nothing
      // -----------------------------------------------------------------------
      if (argc != 2) {
	fprintf(stderr, "excess parameters given\n");
	exit(-1);
      }
      break;
	
    default:
      fprintf(stderr, "invalid action %c\n", action);
      exit(-1);
    }
  }
  else {
    // action == x
    // -----------------------------------------------------------------------
    // x: explore possible sizes
    // -----------------------------------------------------------------------
    // example: ./warpingSIMDTest x
    //
    // TODO: very strange: for
    // TODO: #define MEASTYPE SIMDFloat
    // TODO: #define MATCHTYPE SIMDShort
    // TODO: and SW==32 (but not for SW==16 and SW==64),
    // TODO: and for C++11 and C++17 (but not for C++98)
    // TODO: we get an assertion error for
    // TODO: w = 416, nAlpha_nPsi = 416, step =   1
    // TODO: warpingSIMDTest: ComplexSearch.H:386:
    // TODO: void ns_simd::ComplexSearch<SPSType, MatchType, CompassType,
    // TODO: SIMD_WIDTH, SIMD_ALIGN,
    // TODO: SearchTemplateClass>::fine(PartialSearcherClass<SPSType,
    // TODO: MatchType, SIMD_WIDTH, SIMD_ALIGN>&, bool, int, int, int&,
    // TODO: double&, int&, double&, double&, ns_simd::SIMDImage<MatchType,
    // TODO: 1, 1>&, bool) [with PartialSearcherClass =
    // TODO: ns_simd::MinWarpingPartial; SPSType = unsigned char; MatchType
    // TODO: = short int; CompassType = int; int SIMD_WIDTH = 32; int
    // TODO: SIMD_ALIGN = 32; SearchTemplateClass =
    // TODO: ns_simd::MinWarpingTemplate]: Assertion `_dMin !=
    // TODO: matchFineRes->invalid' failed.
    // 
    double alphaMin, psiMin, dMin;
    int iAlphaMin, iPsiMin, iAlphaMinInv, iPsiMinInv;
    int jAlphaOff, jPsiOff;
    SIMDImage<MATCHTYPE,1,1> fineMatch;
    // signed long long t1, t2;
    double t1, t2;
    int simd_img_elems = SW / sizeof(IMGTYPE);
    for (int w = simd_img_elems; w < 500; w += simd_img_elems) {
      for (int step = 1; step < 6; step++)
	if (w % step == 0) {
	  int nAlpha_nPsi = w / step;
	  bool valid = true;
	  try {
	    SIMDImage<IMGTYPE,SW,SA,Panorama>
	      ss(w, 10, Panorama(0.0, 9.0)), cv(w, 10, Panorama(0.0, 9.0));
	    // 21. Apr 20 (rm): B. Volkmer noticed that images were
	    // used uninitialized, leading to weird effects, now
	    // initialized randomly
	    randomImage(ss.data, ss.size);
	    randomImage(cv.data, cv.size);
	    WarpingBundle<IMGTYPE,PROCTYPE,
			  MEASTYPE,SPSTYPE,
			  MATCHTYPE,COMPASSTYPE,
			  SW,SA>
	      warping
	      (w, nAlpha_nPsi, nAlpha_nPsi,
	       nScalePlanes, maxScaleFactor, maxThresholdMinWarping,
	       rhoMinMinWarping, rhoMaxMinWarping,
	       rhoMinWarping, rhoMaxWarping, nRhoWarping,
	       firstPhase,
	       searchMethodWarping,
	       searchMethodMinWarping,
	       HarrisParameter(harrisK, harrisBinomialFilterCount),
	       sigmoidW, sigmoidW0);
	    // 29. Sep 19 (rm): localPixelScale: all elements 1.0
	    std::vector<double> localPixelScale(warping.spsComp->numPreprocs(),
						1.0);
	    // 27. Feb 18 (rm): vert.res. and hor. taken from images
	    warping.run(ss, cv, 
			/* interpolation */ 0,
			/* pixelScale */ localPixelScale,
			/* postScale */ 100,
			/* alphaEst */ 0.0, /* psiEst */ 0.0,
			searcher,
			compassAcceleration,
			psiFraction, quantilCompassFraction,
			partialSearch,
			alphaRad, psiRad,
			partialFullSearch,
			doubleSearch,
			fineSearch,
			doubleSearchFine,
			fineSearchSurround,
			alphaMin, psiMin, dMin, 
			iAlphaMin, iPsiMin, 
			iAlphaMinInv, iPsiMinInv,
			jAlphaOff, jPsiOff, fineMatch,
			t1, t2);
	  }
	  catch (SIMDException &e) {
	    valid = false;
#if 0
	    fprintf(finfo, 
		    "\tw = %3d, nAlpha_nPsi = %3d, step = %3d %s %s\n",
		    w, nAlpha_nPsi, step,
		    e.loc.c_str(),
		    e.err.c_str());
#endif
	  }
	  if (valid) {
	    printf("w = %3d, nAlpha_nPsi = %3d, step = %3d\n",
		   w, nAlpha_nPsi, step);
	  }
	}
    }	
  }
  HEAP_INFO("at end");
  return 0;
}
