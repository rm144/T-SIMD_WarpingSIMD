// ===========================================================================
// 
// WarpingSPSComputation.H --
// abstract classes for algorithms computing the scale-plane stack
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

#ifndef _WARPING_SPS_COMPUTATION_H_
#define _WARPING_SPS_COMPUTATION_H_

#include <vector>
#include <limits>
#include "SIMDException.H"
#include "WarpingSPS.H"
#include "SIMDVecAll.H"
#include "SIMDImage.H"

namespace ns_simd {

  // ===========================================================================
  // WarpingSPSComputation
  // ===========================================================================

  // this is the abstract interface for the computation of the
  // scale-plane stack; core function is computeSPS, 
  //
  // ImgType is the type used for the input images
  // ProcType is the type used for image preprocessing
  // MeasType is the type used for the components of the distance measure
  // SPSType is the type use to store the scale-plane stack
  // 
  // TODO: this is specific for normalized measures, and probably only works
  // TODO: for normalized measures with sum-of-abs in the denominator

  template <typename ImgType,
	    typename ProcType, typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputation
  {
    
  public:

    // ======================== ctor, dtor =========================

    WarpingSPSComputation() {}
    virtual ~WarpingSPSComputation() {}
    
    // ======================== virtual methods ============================

    // computation of scale-plane stack 
    // (not possible to restrict ImgType to this function since
    // templates may not be virtual)
    // 27. Feb 18 (rm): vert.res. and hor. taken from images
    virtual void 
    computeSPS
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &ss,
     const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &cv,
     const std::vector<double> &scaleFactors, 
     int interpolation,
     const std::vector<double> &pixelScale,
     double postScale,
     WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps) = 0;
  
    // given an example image, compute maximum denominator of the
    // normalized measure; a maximum over multiple such images should be
    // computed and passed to idealPixelScale
    // (not possible to restrict ImgType to this function since
    // templates may not be virtual)
    // TODO: this mechanism may not be general enough, what about
    // TODO: sumabs * sumabs in the denominator: sumabs may overflow, and
    // TODO: the product may also overflow
    virtual std::vector<double>
    maxDenom
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img) = 0;

    virtual int
    numPreprocs() const = 0;

    // if computed with real numbers, maximal value of the normalized measure
    // (used for WarpingSPSComputation::idealPostScale())
    virtual double 
    measureMax() const = 0;

    // whether computation considers invalid pixels (default: false)
    virtual bool
    considersInvalidPixels() const
    {
      return false;
    }
    
    // compute a diagram showing all intermediate input images for SPS
    // 27. Feb 18 (rm): vert.res. and hor. taken from images
    virtual void
    computeSPSInputDiagram
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &ss,
     const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &cv,
     const std::vector<double> &scaleFactors, 
     int interpolation,
     const std::vector<double> &pixelScale,
     double postScale,
     WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps,
     // outputs:
     MeasType &allMax, MeasType &allMin, MeasType &allRange,
     SIMDImage<SIMDByte,SIMD_WIDTH,SIMD_ALIGN> &spsInput) = 0;


    // ===================== non-virtual methods ===========================

    // given the maximum denominator of a normalized measure, compute
    // the ideal pixel scale (scaling of pixels in transformation to
    // MeasType representation) for each preprocessing,
    // for MeasType==SIMDFloat it may be useful to pass measTypeMax instead
    // of the default argument
    std::vector<double>
    idealPixelScale(const std::vector<double> &maxDenom,
		    double measTypeMax
		    = double(SIMDTypeInfo<MeasType>::max())) const
    {
      std::vector<double> pixelScale(maxDenom.size());
      for (size_t i = 0; i < maxDenom.size(); i++)
        pixelScale[i] = measTypeMax / maxDenom[i];
      // maximum denominator of measure multiplied be
      // pixelScale, should not exceed max of MeasType
      return pixelScale;
    }

    // given the maximal value an element in the SPS is allowed to have
    // in the search phase (where there is a summation over all rows),
    // compute the ideal post scale (scaling of distance measure in
    // transformation to integer representation)
    // uses virtual measureMax()
    // for SPSType==SIMDFloat it may be useful to pass spsTypeMax instead
    // of the default argument
    double 
    idealPostScale(double maxElemSearchPhase,
		   double spsTypeMax 
		   = double(SIMDTypeInfo<SPSType>::max())) const
    {
      // two conditions: normalized measure is in [0,this->measureMax()], 
      // therefore this->measureMax() * postScale should not exceed
      // (1) spsTypeMax, the maximum of the storage type of the SPS,
      // (2) maximum element value in summation in search phase 
      //     (loose coupling to search phase, we don't know MatchType here)
      return std::min(spsTypeMax, maxElemSearchPhase) / this->measureMax();
    }
 
  };
  
  // ===========================================================================
  // WarpingSPSComputationSingle
  // ===========================================================================

  // abstract interface defining preprocessing, column measure
  // (applied to individual images), joint measure (applied to a pair
  // of image), normalization method, for the case of a single
  // preprocessed image (typically vertically edge-filtered)

  template <typename ImgType,
	    typename ProcType, typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationSingle :
    public WarpingSPSComputation<ImgType,ProcType,MeasType,SPSType,
				 SIMD_WIDTH,SIMD_ALIGN>
  {
    
  protected:

    // computeSPS should store input images
    bool storeSPSInput;
    // input images in sps computation
    std::vector< SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> >
    ssInputs, cvInputs;

    // ===================== non-virtual methods ===========================

    // compute SPS and intermediate input images
    // 27. Feb 18 (rm): vert.res. and horizon taken from images
    void
    computeSPSInput
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &ss,
     const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &cv,
     const std::vector<double> &scaleFactors, 
     int interpolation,
     const std::vector<double> &pixelScale,
     double postScale,
     WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps)    
    {
      storeSPSInput = true;
      ssInputs.resize(scaleFactors.size());
      cvInputs.resize(scaleFactors.size());
      this->computeSPS(ss, cv, scaleFactors,
		       interpolation, pixelScale, postScale, sps);
      storeSPSInput = false;
    }

  public:

    // ======================== ctor, dtor =========================
    
    WarpingSPSComputationSingle() :
      WarpingSPSComputation<ImgType,ProcType,MeasType,SPSType,
			    SIMD_WIDTH,SIMD_ALIGN>(),
      storeSPSInput(false) {}

    virtual ~WarpingSPSComputationSingle() {}

    // ======================== virtual methods ============================

    virtual void
    convertImage
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &in,
     SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &out) = 0;

    // magnify and scale image
    // 27. Feb 18 (rm): verticalResolution and horizon taken from img
    virtual void
    magnifyImage
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     int interpolation,
     double magnifyScale,
     double multiplyScale,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outImg) = 0;
    
    // copy and scale image
    virtual void
    copyImage
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     double multiplyScale,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outImg) = 0;
    
    // re-arrange image (makeView instead of pure virtual)
    virtual void
    reArrange
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &imgR,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // TODO: this is fast but may be dangerous, true copy may be advisable
      imgR.view(img);
      // imgR = img;
    }

    // un-re-arrange image (copy instead of pure virtual)
    // typically used for visualization, may be slow
    virtual void
    unReArrange
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &imgR,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // this doesn't work, true copy is needed
      // img.view(imgR);
      img = imgR;
    }

    // preprocessing applied to the input images
    virtual void
    preprocessing
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &output,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par) = 0;
  
    // measure computed over all columns of one input image 
    virtual void 
    columnMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &output,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par) = 0;
  
    // measure computed over pairs of columns in two input images
    // also takes raw (intermediate) column measures as input
    // and produces final column measures as output
    // (it also takes raw columns measures colMeasRaw1/2 as input
    // and computes column measures colMeas1/2 and a combined measure
    // colMeas12 (e.g a sum))
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas12,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par) = 0;
  
    // normalization of jointMeasure and ssMeasure/cvMeasure
    virtual void
    normalization
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &ssMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &cvMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &combMeasure,
     double postScale,
     SIMDImage<SPSType,SIMD_WIDTH,SIMD_ALIGN> &scalePlane,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par) = 0;
    
    // check whether input images with given width could be processed
    // this definition should be useful for many cases:
    virtual bool
    validForWidth(int width) const
    {
      return 
	// ProcType->MeasType: uses magnifyImage, copyImage:
	((width % NUM_SIMDVECS_ELEMENTS(MeasType,ProcType,SIMD_WIDTH) == 0) &&
	 // MeasType->SPSType: uses computeScalePlaneNormalizedAdd:
	 (width % NUM_SIMDVECS_ELEMENTS(SPSType,MeasType,SIMD_WIDTH) == 0));
    }
    
    // compute a diagram showing all intermediate input images for SPS
    // 27. Feb 18 (rm): vert.res. and hor. taken from images
    // TODO: option to leave out an "invalid" value
    // TODO: needs "max/minExceptInvalid()" etc.
    void
    computeSPSInputDiagram
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &ss,
     const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &cv,
     const std::vector<double> &scaleFactors, 
     int interpolation,
     const std::vector<double> &pixelScale,
     double postScale,
     WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps,
     // outputs:
     MeasType &allMax, MeasType &allMin, MeasType &allRange,
     SIMDImage<SIMDByte,SIMD_WIDTH,SIMD_ALIGN> &spsInput)
    {
      computeSPSInput(ss, cv, scaleFactors,
		      interpolation, pixelScale, postScale, sps);
      unsigned nPlanes = scaleFactors.size();
      MeasType invalid = SIMDTypeInfo<MeasType>::max();
      // first: determine max and min over all images
      allMax = SIMDTypeInfo<MeasType>::min();
      allMin = SIMDTypeInfo<MeasType>::max();
      for (unsigned planeIdx = 0; planeIdx < nPlanes; planeIdx++) {
	MeasType _min, _max;
	_max =  this->considersInvalidPixels() ?
	  maxExcept(ssInputs[planeIdx], invalid) : max(ssInputs[planeIdx]);
	_min =  min(ssInputs[planeIdx]);
	allMax = std::max(allMax, _max);
	allMin = std::min(allMin, _min);
	_max =  this->considersInvalidPixels() ?
	  maxExcept(cvInputs[planeIdx], invalid) : max(cvInputs[planeIdx]);
	_min =  min(cvInputs[planeIdx]);
	allMax = std::max(allMax, _max);
	allMin = std::max(allMin, _min);
      }
      allRange = std::max(allMax, MeasType(-allMin));
      // write into same image
      SIMDImage<SIMDByte,SIMD_WIDTH,SIMD_ALIGN> tmp;
      int w = ssInputs[0].w, h = ssInputs[0].h;
      const int margin = SIMDVec<SIMDByte,SIMD_WIDTH>::elements;
      spsInput.resize(2 * (w + margin) + margin,
		      nPlanes * (h + margin) + margin);
      spsInput.fill(255);
      double shift = 128.0, scale = 255.0 / allRange;  
      SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> tmpSS, tmpCV;
      for (unsigned planeIdx = 0; planeIdx < nPlanes; planeIdx++) {
	if (this->considersInvalidPixels()) {
	  replace(ssInputs[planeIdx], invalid, allRange, tmpSS);
	  replace(cvInputs[planeIdx], invalid, allRange, tmpCV);
	}
	else {
	  tmpSS.view(ssInputs[planeIdx]);
	  tmpCV.view(cvInputs[planeIdx]);
	}
	copyScaleAndShift(tmpSS, scale, shift, tmp);
	insert(tmp, margin,
	       margin + planeIdx * (h + margin), spsInput);
	copyScaleAndShift(tmpCV, scale, shift, tmp);
	insert(tmp, w + 2 * margin, 
	       margin + planeIdx * (h + margin), spsInput);
      }
    }

    int
    numPreprocs() const {
      return 1;
    }

  };
  
  // ===========================================================================
  // WarpingSPSComputationStack
  // ===========================================================================

  // abstract interface defining preprocessing, column measure
  // (applied to individual image stacks), joint measure (applied to a pair
  // of image stacks), normalization method, for the case of a stack of
  // preprocessed image (typically vertically edge-filtered)

  template <typename ImgType,
	    typename ProcType, typename MeasType, typename SPSType,
            int NPREPROCS,
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationStack :
    public WarpingSPSComputation<ImgType,ProcType,MeasType,SPSType,
				 SIMD_WIDTH,SIMD_ALIGN>
  {
    
  protected:

    // computeSPS should store input images
    bool storeSPSInput;
    // input images in sps computation
    std::vector< SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> >
    ssInputs, cvInputs;

    // ===================== non-virtual methods ===========================

    // compute SPS and intermediate input images
    // 27. Feb 18 (rm): vert.res. and horizon taken from images
    void
    computeSPSInput
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &ss,
     const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &cv,
     const std::vector<double> &scaleFactors, 
     int interpolation,
     const std::vector<double> &pixelScale,
     double postScale,
     WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps)
    {
      storeSPSInput = true;
      ssInputs.resize(scaleFactors.size());
      cvInputs.resize(scaleFactors.size());
      this->computeSPS(ss, cv, scaleFactors,
		       interpolation, pixelScale, postScale, sps);
      storeSPSInput = false;
    }

  public:

    // ======================== ctor, dtor =========================
    
    WarpingSPSComputationStack() :
      WarpingSPSComputation<ImgType,ProcType,MeasType,SPSType,
			    SIMD_WIDTH,SIMD_ALIGN>(),
      storeSPSInput(false) {}

    virtual ~WarpingSPSComputationStack() {}

    // ======================== virtual methods ============================
    

    virtual void
    convertImage
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &in,
     SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &out) = 0;

    // copy and scale image
    virtual void
    copyImage
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     double multiplyScale,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outImg) = 0;

    // copy and scale stack while each plane has its own multiplyScale
    // multiplyScale[NPREPROCS] is assumed
    virtual void
    copyStack
    (const SIMDStack<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &stk,
     const std::vector<double> &multiplyScale,
     SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outStk) = 0;

    // re-arrange image (makeView instead of pure virtual)
    virtual void
    reArrange
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &imgR,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // TODO: this is fast but may be dangerous, true copy may be advisable
      imgR.view(img);
      // imgR = img;
    }

    // un-re-arrange image (copy instead of pure virtual)
    // typically used for visualization, may be slow
    virtual void
    unReArrange
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &imgR,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // this doesn't work, true copy is needed
      // img.view(imgR);
      img = imgR;
    }

    // re-arrange image (makeView instead of pure virtual)
    virtual void
    reArrange
    (const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &stk,
     SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &stkR,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // TODO: this is fast but may be dangerous, true copy may be advisable
      // TODO: bv: enable following when SIMDStack has view function
      // stkR.view(stk);
      stkR = stk;
    }

    // un-re-arrange image (copy instead of pure virtual)
    // typically used for visualization, may be slow
    virtual void
    unReArrange
    (const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &stkR,
     SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &stk,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // this doesn't work, true copy is needed
      // stk.view(stkR);
      stk = stkR;
    }

    // preprocessing applied to the input images
    virtual void
    preprocessing
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDStack<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &output,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par) = 0;
  
    // measure computed over all columns of one input image 
    virtual void 
    columnMeasure
    (const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &output,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par) = 0;
  
    // measure computed over pairs of columns in two input images
    // also takes raw (intermediate) column measures as input
    // and produces final column measures as output
    // (it also takes raw columns measures colMeasRaw1/2 as input
    // and computes column measures colMeas1/2 and a combined measure
    // colMeas12 (e.g a sum))
    virtual void 
    jointMeasure
    (const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par) = 0;
  
    // normalization of jointMeasure and ssMeasure/cvMeasure
    virtual void
    normalization
    (const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeasure,
     const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &ssMeasure,
     const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &cvMeasure,
     const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &combMeasure,
     double postScale,
     SIMDImage<SPSType,SIMD_WIDTH,SIMD_ALIGN> &scalePlane,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par) = 0;
    
    // check whether input images with given width could be processed
    // this definition should be useful for many cases:
    virtual bool
    validForWidth(int width) const
    {
      return 
	// ProcType->MeasType: uses magnifyImage, copyImage:
	((width % NUM_SIMDVECS_ELEMENTS(MeasType,ProcType,SIMD_WIDTH) == 0) &&
	 // MeasType->SPSType: uses computeScalePlaneNormalizedAdd:
	 (width % NUM_SIMDVECS_ELEMENTS(SPSType,MeasType,SIMD_WIDTH) == 0));
    }
    
    // compute a diagram showing all intermediate input images for SPS
    // 27. Feb 18 (rm): vert.res. and hor. taken from images
    void
    computeSPSInputDiagram
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &ss,
     const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &cv,
     const std::vector<double> &scaleFactors, 
     int interpolation,
     const std::vector<double> &pixelScale,
     double postScale,
     WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps,
     // outputs:
     MeasType &allMax, MeasType &allMin, MeasType &allRange,
     SIMDImage<SIMDByte,SIMD_WIDTH,SIMD_ALIGN> &spsInput)
    {
      // TODO: this always shows only the following preprocessing image
      computeSPSInput(ss, cv, scaleFactors,
		      interpolation, pixelScale, postScale, sps);
      unsigned nPlanes = scaleFactors.size();
      // first: determine max and min over all images
      allMax = SIMDTypeInfo<MeasType>::min();
      allMin = SIMDTypeInfo<MeasType>::max();
      for (unsigned planeIdx = 0; planeIdx < nPlanes; planeIdx++) {
        for (int d = 0; d < this->numPreprocs(); d++) {
          MeasType _min, _max;
          _max =  max(ssInputs[planeIdx][d]);
          _min =  min(ssInputs[planeIdx][d]);
          allMax = std::max(allMax, _max);
          allMin = std::min(allMin, _min);
          _max =  max(cvInputs[planeIdx][d]);
          _min =  min(cvInputs[planeIdx][d]);
          allMax = std::max(allMax, _max);
          allMin = std::max(allMin, _min);
        }
      }
      allRange = std::max(allMax, MeasType(-allMin));
      // write into same image
      SIMDImage<SIMDByte,SIMD_WIDTH,SIMD_ALIGN> tmp;
      // asserting all preprocs have same width and height
      int w = ssInputs[0][0].w, h = ssInputs[0][0].h;
      const int margin = SIMDVec<SIMDByte,SIMD_WIDTH>::elements;
      spsInput.resize(this->numPreprocs() * 2 * (w + margin) + margin,
          nPlanes * (h + margin) + margin);
      spsInput.fill(255);
      double shift = 128.0, scale = 255.0 / allRange;
      int verticalOffset = 0;
      for (int d = 0; d < this->numPreprocs(); d++) {
        for (unsigned planeIdx = 0; planeIdx < nPlanes; planeIdx++) {
          copyScaleAndShift(ssInputs[planeIdx][d], scale, shift, tmp);
          insert(tmp, verticalOffset + margin,
                 margin + planeIdx * (h + margin), spsInput);
          copyScaleAndShift(cvInputs[planeIdx][d], scale, shift, tmp);
          insert(tmp, verticalOffset + w + 2 * margin, 
                 margin + planeIdx * (h + margin), spsInput);
        }
        verticalOffset += 2 * w + 2* margin;
      }
    }

    int
    numPreprocs() const {
      return NPREPROCS;
    }

  };
  
}

#endif
