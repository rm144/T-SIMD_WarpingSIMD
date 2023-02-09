// ===========================================================================
// 
// WarpingSPSComputationCollection.H --
// collection of different classes derived from WarpingSPSComputation
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

#ifndef _WARPING_SPS_COMPUTATION_COLLECTION_H_
#define _WARPING_SPS_COMPUTATION_COLLECTION_H_

#include "WarpingSPSComputation.H"
// 18. Apr 18 (rm): added, thanks to Benedikt Volkmer
#include "SIMDImageFunctions.H"

namespace ns_simd {

  // TODO: there is some over-/underflow in the normalized measure if
  // TODO: pixelScale is too large, even though everything is in 
  // TODO: saturated arithmetics, check whether this over-/underflow can be
  // TODO: avoided
  // TODO: analyze effect of pixelScale/postScale on performance
  // TODO: idea: nonlinear pixel functions?

  // class hierarchy:
  //
  // WarpingSPSComputation
  // -> WarpingSPSComputationSingle
  //   -> WarpingSPSComputation1 (defines computeSPS)
  //    -> WarpingSPSComputation2 (def's convertImage, magnifyImage, copyImage)
  //      -> WarpingSPSComputationEdgeAbs (for edge-filtered input + abs. sum)
  //        -> WarpingSPSComputationEdgeAbs1 (element-wise rewrite)
  //          -> WarpingSPSComputationEdgeNSAD1 (for NSAD joint measure)
  //          -> WarpingSPSComputationASC1 (for ASC joint measure)
  //          -> more...
  //        -> WarpingSPSComputationEdgeAbs2 (with reArrange)
  //          -> WarpingSPSComputationEdgeNSAD2 (with reArrange)
  //          -> WarpingSPSComputationASC2 (with reArrange)
  //          -> more...
  //
  // contributed by Annika Hoffmann: (for images with undefined regions)
  //                                 (only for unswizzled SPS)
  //
  // WarpingSPSComputation
  // -> WarpingSPSComputationSingle
  //    -> WarpingSPSComputation1
  //       -> WarpingSPSComputationTilt (convert~,magnify~,copyImage)
  //          -> WarpingSPSComputationAbsTilt (preprocess.,col.Meas.,maxDen.)
  //             -> WarpingSPSComputationNSADTilt (jointMeas., norm., meas.Max.)
  //
  // contributed by Christoph Berganski: (images with undefined regions,
  //                                     with edge filter)
  //
  // WarpingSPSComputation
  // -> WarpingSPSComputationSingle
  //    -> WarpingSPSComputation1
  //       -> WarpingSPSComputationTilt (convert~,magnify~,copyImage)
  //          -> WarpingSPSComputationEdgeAbsTilt (C.B.)
  //             -> WarpingSPSComputationEdgeNSADTilt (C.B.)
  //
  // contributed by Benedikt Volkmer: (edge and corner detectors,
  //                                   some with logistic regression)
  //
  // WarpingSPSComputation
  // -> WarpingSPSComputationSingle
  //    -> WarpingSPSComputation1Alt (computeSPS)
  //       -> WarpingSPSComputation2Alt (convert~,magnify~,copyImage)
  //          -> WarpingSPSComputationCornerAbs (preProc.,col.Meas.,maxDen.)
  //             -> WarpingSPSComputationCornerAbs2 (reArrange,norm.)
  //                -> WarpingSPSComputationCornerNSAD2 (jointM.,meas.Max.)
  //          -> WarpingSPSCompuationEdgeAbsAlt (preProc.,col.Meas.,maxDen.)
  //             -> WarpingSPSComputationEdgeAbs2Alt (reArrange,norm.)
  //                -> WarpingSPSCompuationEdgeNSAD2Alt (jointMeas.,meas.Max.)
  // -> WarpingSPSComputationStack
  //    -> WarpingSPSComputationSNSAD
  //      -> WarpingSPSComputationSNSADHarris
  
  // ===========================================================================
  // WarpingSPSComputation1
  // ===========================================================================

  // only defines method computeSPS
  template <typename ImgType,
	    typename ProcType, typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputation1 
    : public WarpingSPSComputationSingle<ImgType, 
					 ProcType, MeasType, SPSType,
					 SIMD_WIDTH, SIMD_ALIGN>
  {    
    // 2. Mar 18 (rm): made images used in computeSPS() members of
    // WarpingSPSComputation1 (should accelerate allocation and avoids
    // deallocation), was slower, undone
    
  public:

    // 27. Feb 18 (rm): vert.res. and hor. now info. in images
    virtual void 
    computeSPS
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &ss,
     const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &cv,
     const std::vector<double> &scaleFactors, 
     int interpolation,
     const std::vector<double> &pixelScale,
     double postScale,
     WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps)
    {
      if (pixelScale.size() != 1)
	throw SIMDException("WarpingSPSComputation1::computeSPS",
			    "invalid size of pixelScale");
      const int w = sps.param.w, h = ss.h;
      if ((h != cv.h) || (ss.w != cv.w) || (w != ss.w)) {
	char errStr[256];
	sprintf(errStr,
		"mismatch in ss, cv, sps dimensions: " 
		"ss %d x %d, cv %d x %d, sps.param.w = %d",
		ss.w, ss.h, cv.w, cv.h, w);
	throw SIMDException("WarpingSPSComputation1::computeSPS",
			    errStr);
      }
      if (ss.addOn != cv.addOn) {
	char errStr[256];
	sprintf(errStr,
		"panorama parameters of ss (%g, %g) and cv (%g, %g) differ",
		ss.addOn.verticalResolution,
		ss.addOn.horizon,
		cv.addOn.verticalResolution,
		cv.addOn.horizon);
	throw SIMDException("WarpingSPSComputation1::computeSPS",
			    errStr);
      }
      // check validity for given image width
      if (!this->validForWidth(w))
	throw SIMDException("WarpingSPSComputation1::computeSPS",
			    "computation not possible for given image width");
      // get input images
        // convert from input type to preprocessing type
      SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> ssOrigF, cvOrigF;
      this->convertImage(ss, ssOrigF);
      this->convertImage(cv, cvOrigF);
      // preprocess input images
      SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> ssF, cvF;
      this->preprocessing(ssOrigF, ssF, sps.param);
      this->preprocessing(cvOrigF, cvF, sps.param);
      // we first process the unmagnified images since the results are needed
      // multiple times
      // 2. Mar 18 (rm): *SR images are the re-arranged versions,
      // earlier in-place re-arrange required heap operations (width differs)
      SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> ssS, cvS, magS,
	ssSR, cvSR, magSR, *cvSPtr, *ssSPtr;
      SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> ssMeasureS, cvMeasureS,
	cvMeasFinal, ssMeasFinal, combMeas,
	magMeasureS, jointMeas, *cvMeasureSPtr, *ssMeasureSPtr;
      this->copyImage(cvF, pixelScale[0], cvS);
      this->reArrange(cvS, cvSR, sps.param);
      this->columnMeasure(cvSR, cvMeasureS, sps.param);
      this->copyImage(ssF, pixelScale[0], ssS);
      this->reArrange(ssS, ssSR, sps.param);
      this->columnMeasure(ssSR, ssMeasureS, sps.param);
  
      // create all scale plane images
      for (unsigned int spIndex = 0; spIndex < scaleFactors.size(); spIndex++) {

	// edge images are scaled
	double scale = scaleFactors[spIndex];

	// printf("spIndex %u, scale %g\n", spIndex, scale);
	if (scale > 1.0) {
	
	  // ===== magnify cv by 1/scale and use unmodified ss =====
	  // 27. Feb 18 (rm): vert.res. and hor. taken from image
	  this->magnifyImage(cvF, 
			     interpolation, 1.0 / scale, 
			     pixelScale[0], magS);
	  this->reArrange(magS, magSR, sps.param);
	  cvSPtr = &magSR;
	  this->columnMeasure(magSR, magMeasureS, sps.param);
	  cvMeasureSPtr = &magMeasureS;
	  ssSPtr = &ssSR;
	  ssMeasureSPtr = &ssMeasureS;
	
	} else if (scale < 1.0) {
	
	  // ===== magnify ss by scale and use unmodified cv =====
	  cvSPtr = &cvSR;
	  cvMeasureSPtr = &cvMeasureS;
	  // 27. Feb 18 (rm): vert.res. and hor. taken from image
	  this->magnifyImage(ssF,
			     interpolation, scale, 
			     pixelScale[0], magS);
	  this->reArrange(magS, magSR, sps.param);
	  ssSPtr = &magSR;
	  this->columnMeasure(magSR, magMeasureS, sps.param);
	  ssMeasureSPtr = &magMeasureS;
	
	} else {
	
	  // ===== scale = 1.0: use unmodified ss and cv =====
	  cvSPtr = &cvSR;
	  cvMeasureSPtr = &cvMeasureS;
	  ssSPtr = &ssSR;
	  ssMeasureSPtr = &ssMeasureS;
	}
	// store input images if requested
	if (this->storeSPSInput) {
	  // 27. Jul 18 (rm)
	  this->unReArrange(*ssSPtr, this->ssInputs[spIndex], sps.param);
	  this->unReArrange(*cvSPtr, this->cvInputs[spIndex], sps.param);
	  // this->ssInputs[spIndex] = *ssSPtr;
	  // this->cvInputs[spIndex] = *cvSPtr;
	}
	// compute distance term
	this->jointMeasure(*cvSPtr, *ssSPtr, 
			   *cvMeasureSPtr, *ssMeasureSPtr,
			   jointMeas, cvMeasFinal, ssMeasFinal, combMeas,
			   sps.param);
	// compute normalized measure
	this->normalization(jointMeas, ssMeasFinal, cvMeasFinal, combMeas,
			    postScale,
			    sps.stack[spIndex], sps.param);

      }
      sps.complement();
    }
  };

  // ===========================================================================
  // WarpingSPSComputation2
  // ===========================================================================

  // defines convertImage, magnifyImage, copyImage
  template <typename ImgType,
	    typename ProcType, typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputation2 
    : public WarpingSPSComputation1<ImgType, 
				    ProcType, MeasType, SPSType,
				    SIMD_WIDTH, SIMD_ALIGN>
  {
  public:
    // convert image
    virtual void
    convertImage
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &in,
     SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &out)
    {
      // TODO: this assignment is based on type casts, is that always ok?
      // TODO: what is actually happening between signed/unsigned types and 
      // TODO: between integer/float types?
      out = in;
    }

    // magnify and scale image
    // 27. Feb 18 (rm): vert.res. and hor. taken from img
    virtual void
    magnifyImage
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     int interpolation,
     double magnifyScale,
     double multiplyScale,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outImg)
    {
      magnifyAndScale(img,
		      interpolation, magnifyScale, multiplyScale,
		      outImg);
    }
    
    // copy and scale image
    virtual void
    copyImage
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     double multiplyScale,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outImg)
    {
      copyAndScale(img, multiplyScale, outImg);
    }
  };

  // ===========================================================================
  // WarpingSPSComputationEdgeAbs
  // ===========================================================================

  // only defines methods preprocessing, columnMeasure,
  // maxDenom
  template <typename ImgType,
	    typename ProcType, typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeAbs
    : public WarpingSPSComputation2<ImgType, 
				    ProcType, MeasType, SPSType,
				    SIMD_WIDTH, SIMD_ALIGN>
  {
  public:
  
    virtual void 
    preprocessing
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &preprocessed,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      verticalEdgeFilter(input, preprocessed);
    }

    virtual void 
    columnMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      sumOfAbs(input, colMeas);
    }

    virtual std::vector<double>
    maxDenom
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img)
    {
      std::vector<double> md(1);
      // all processing is done in SIMDFloat to avoid overflow
      SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN,Panorama> inputImg, edgeImg;
      SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN> sumOfAbsImg;
      inputImg = img;
      verticalEdgeFilter(inputImg, edgeImg);
      sumOfAbs(edgeImg, sumOfAbsImg);
      // 2.0* : two sums are added in the denominator
      md[0] = (2.0 * max(sumOfAbsImg));
      return md;
    }

  };

  // ===========================================================================
  // measures with element-wise rewrite to SPS (1)
  // ===========================================================================

  template <typename ImgType, 
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN, int INVERT>
  class WarpingSPSComputationEdgeAbs1
    :  public WarpingSPSComputationEdgeAbs<ImgType, 
					   ProcType, MeasType, SPSType,
					   SIMD_WIDTH, SIMD_ALIGN>
  {
  public:

    virtual void
    normalization
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &ssMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &cvMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     double postScale,
     SIMDImage<SPSType,SIMD_WIDTH,SIMD_ALIGN> &scalePlane,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      // INVERT: false: without inversion (distance measure), true:
      // with inversion (correlation)
      computeScalePlaneNormalizedAdd1<INVERT>
	(jointMeasure, ssMeasure, cvMeasure, 
	 par.wSPS,
	 par.modulo2wpw, par.shuffle, 
	 postScale, scalePlane);
    }

  };

  // ===========================================================================
  // NSAD with edge filter (derived from WarpingSPSComputationEdgeAbs1)
  // ===========================================================================

  template <typename ImgType, 
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeNSAD1
    :  public WarpingSPSComputationEdgeAbs1<ImgType, 
					    ProcType, MeasType, SPSType,
					    SIMD_WIDTH, SIMD_ALIGN, false>
  {
  public:
  
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // puts("NSAD");
      normSumAbsDiffTerm1(input1, input2, jointMeas);
      // colMeas1 = colMeasRaw1;
      // colMeas2 = colMeasRaw2;
      colMeas1.view(colMeasRaw1);
      colMeas2.view(colMeasRaw2);
    }

    // NSAD measure is in [0,1]
    virtual double 
    measureMax() const { return 1.0; }

  };

  // ===========================================================================
  // ASC with edge filter (derived from WarpingSPSComputationEdgeAbs1)
  // ===========================================================================

  template <typename ImgType,
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeASC1 
    :  public WarpingSPSComputationEdgeAbs1<ImgType,
					    ProcType, MeasType, SPSType,
					    SIMD_WIDTH, SIMD_ALIGN, true>
  {
  public:
  
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // puts("ASC");
      asc2Term1(input1, input2, jointMeas);
      // colMeas1 = colMeasRaw1;
      // colMeas2 = colMeasRaw2;
      colMeas1.view(colMeasRaw1);
      colMeas2.view(colMeasRaw2);
    }
  
    // ASC+ measure is in [0,2]
    virtual double 
    measureMax() const { return 2.0; }

  };

  // ===========================================================================
  // NSADA with edge filter (derived from WarpingSPSComputationEdgeAbs1)
  // ===========================================================================

  template <typename ImgType, 
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeNSADA1
    :  public WarpingSPSComputationEdgeAbs1<ImgType, 
					    ProcType, MeasType, SPSType,
					    SIMD_WIDTH, SIMD_ALIGN, false>
  {
  public:
  
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // puts("NSADA");
      normSumAbsDiffAbsTerm1(input1, input2, jointMeas);
      // colMeas1 = colMeasRaw1;
      // colMeas2 = colMeasRaw2;
      colMeas1.view(colMeasRaw1);
      colMeas2.view(colMeasRaw2);
    }

    // NSAD measure is in [0,1]
    virtual double 
    measureMax() const { return 1.0; }

  };

  // ===========================================================================
  // NSADA with edge filter (derived from WarpingSPSComputationEdgeAbs1)
  // ===========================================================================

  template <typename ImgType, 
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeNSADNSADA1
    :  public WarpingSPSComputationEdgeAbs1<ImgType, 
					    ProcType, MeasType, SPSType,
					    SIMD_WIDTH, SIMD_ALIGN, false>
  {
  public:
  
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // puts("NSAD-NSADA");
      normSumAvgAbsDiffAbsDiffAbsTerm1(input1, input2, jointMeas);
      // colMeas1 = colMeasRaw1;
      // colMeas2 = colMeasRaw2;
      colMeas1.view(colMeasRaw1);
      colMeas2.view(colMeasRaw2);
    }

    // NSAD measure is in [0,1]
    virtual double 
    measureMax() const { return 1.0; }

  };

  // ===========================================================================
  // measures with reArrange and block-wise rewrite to SPS (2)
  // ===========================================================================

  // with reArrange
  template <typename ImgType, 
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN, int INVERT>
  class WarpingSPSComputationEdgeAbs2
    :  public WarpingSPSComputationEdgeAbs<ImgType, 
					   ProcType, MeasType, SPSType,
					   SIMD_WIDTH, SIMD_ALIGN>
  {
  public:

    virtual void
    reArrange
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &imgR,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      swizzleCyclic(img, par.stepAlpha,
		    NUM_INPUT_SIMDVECS(SPSType,MeasType),
		    imgR);
    }

    virtual void
    unReArrange
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &imgR,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      // serial, slow!
      unSwizzleCyclic(imgR, par.stepAlpha,
		      NUM_INPUT_SIMDVECS(SPSType,MeasType),
		      img);
    }
    
    virtual void
    normalization
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &ssMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &cvMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     double postScale,
     SIMDImage<SPSType,SIMD_WIDTH,SIMD_ALIGN> &scalePlane,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      // INVERT: false: without inversion (distance measure), true:
      // with inversion (correlation)
      computeScalePlaneNormalizedAdd2<INVERT>
	(jointMeasure, ssMeasure, cvMeasure, 
	 par.w, par.wSPS, par.stepAlpha, par.nAlpha, par.blkSize,
	 par.modNAlpha,
	 par.shuffleModulo, par.unshuffle,
	 postScale, scalePlane);
    }
    
  };

  // ===========================================================================
  // NSAD with edge filter (derived from WarpingSPSComputationEdgeAbs2)
  // ===========================================================================

  // with reArrange
  template <typename ImgType, 
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeNSAD2 
    :  public WarpingSPSComputationEdgeAbs2<ImgType, 
					    ProcType, MeasType, SPSType,
					    SIMD_WIDTH, SIMD_ALIGN, false>
  {
  public:
  
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      /*
      printf("input1 %d x %d (%g, %g)\n", input1.w, input1.h,
	     input1.addOn.verticalResolution, input1.addOn.horizon);
      printf("input2 %d x %d (%g, %g)\n", input2.w, input2.h,
	     input2.addOn.verticalResolution, input2.addOn.horizon);
      */
      
      // puts("NSAD");
      normSumAbsDiffTerm2
	(input1, input2, 
	 par.w, par.wSPS, par.stepAlpha, par.nAlpha, par.blkSize,
	 par.modNAlpha, jointMeas);
      // colMeas1 = colMeasRaw1;
      // colMeas2 = colMeasRaw2;
      colMeas1.view(colMeasRaw1);
      colMeas2.view(colMeasRaw2);
    }

    // NSAD measure is in [0,1]
    virtual double 
    measureMax() const { return 1.0; }

  };

  // ===========================================================================
  // ASC with edge filter (derived from WarpingSPSComputationEdgeAbs2)
  // ===========================================================================

  // with reArrange
  template <typename ImgType,
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeASC2 
    :  public WarpingSPSComputationEdgeAbs2<ImgType,
					    ProcType, MeasType, SPSType,
					    SIMD_WIDTH, SIMD_ALIGN, true>
  {
  public:
  
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      // puts("ASC");
      asc2Term2(input1, input2, 
		par.w, par.wSPS, par.stepAlpha, par.nAlpha, par.blkSize,
		par.modNAlpha, jointMeas);
      // colMeas1 = colMeasRaw1;
      // colMeas2 = colMeasRaw2;
      colMeas1.view(colMeasRaw1);
      colMeas2.view(colMeasRaw2);
    }
  
    // ASC+ measure is in [0,2]
    virtual double 
    measureMax() const { return 2.0; }
  
  };

  // ===========================================================================
  // ASCDirect (slow, tests only, derived from WarpingSPSComputationEdgeAbs2)
  // ===========================================================================

  // with reArrange
  template <typename ImgType, 
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeASCDirect2 
    :  public WarpingSPSComputationEdgeAbs2<ImgType, 
					    ProcType, MeasType, SPSType,
					    SIMD_WIDTH, SIMD_ALIGN, false>
  {
  public:
  
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      ascDirectTerm(input1, input2, 
		    par.w, par.wSPS, par.stepAlpha, par.nAlpha, par.blkSize,
		    par.modNAlpha, jointMeas);
      // colMeas1 = colMeasRaw1;
      // colMeas2 = colMeasRaw2;
      colMeas1.view(colMeasRaw1);
      colMeas2.view(colMeasRaw2);
    }

    // ASCDirect in [0,2]
    virtual double 
    measureMax() const { return 2.0; }

  };

  // ===========================================================================
  // NSADA with edge filter (derived from WarpingSPSComputationEdgeAbs2)
  // ===========================================================================

  // with reArrange
  template <typename ImgType, 
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeNSADA2 
    :  public WarpingSPSComputationEdgeAbs2<ImgType, 
					    ProcType, MeasType, SPSType,
					    SIMD_WIDTH, SIMD_ALIGN, false>
  {
  public:
  
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      // puts("NSAD-NSADA");
      normSumAbsDiffAbsTerm2
	(input1, input2, 
	 par.w, par.wSPS, par.stepAlpha, par.nAlpha, par.blkSize,
	 par.modNAlpha, jointMeas);
      // colMeas1 = colMeasRaw1;
      // colMeas2 = colMeasRaw2;
      colMeas1.view(colMeasRaw1);
      colMeas2.view(colMeasRaw2);
    }

    // NSADA measure is in [0,1]
    virtual double 
    measureMax() const { return 1.0; }

  };

  // ===========================================================================
  // NSAD-NSADA mix with edge filter (deriv. from
  // WarpingSPSComputationEdgeAbs2)
  // ===========================================================================

  // with reArrange
  template <typename ImgType, 
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeNSADNSADA2 
    :  public WarpingSPSComputationEdgeAbs2<ImgType, 
					    ProcType, MeasType, SPSType,
					    SIMD_WIDTH, SIMD_ALIGN, false>
  {
  public:
  
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      // puts("NSAD-NSADA");
      normSumAvgAbsDiffAbsDiffAbsTerm2
	(input1, input2, 
	 par.w, par.wSPS, par.stepAlpha, par.nAlpha, par.blkSize,
	 par.modNAlpha, jointMeas);
      // colMeas1 = colMeasRaw1;
      // colMeas2 = colMeasRaw2;
      colMeas1.view(colMeasRaw1);
      colMeas2.view(colMeasRaw2);
    }

    // NSADA measure is in [0,1]
    virtual double 
    measureMax() const { return 1.0; }

  };

  // ===========================================================================
  // WarpingSPSComputationTilt
  // contributed by Annika Hoffmann
  // ===========================================================================

  // adapted functions to compute converted, magnified and scaled
  // image with invalid values!
  template <typename ImgType,
	    typename ProcType, typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationTilt 
    : public WarpingSPSComputation1<ImgType, 
				    ProcType, MeasType, SPSType,
				    SIMD_WIDTH, SIMD_ALIGN>
  {
  public:
 
    bool
    considersInvalidPixels() const
    {
      return true;
    }

    // convert image: keep values from ImgType, but use invalid value
    // for proc type
    // if value is invalid (max(ImgType)) in input image, then this
    // value is set to max(ProcType)
      
    virtual void
    convertImage(const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &in,
		 SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &out)
    {
      // code assumes that ProcType is larger or equal in size then ImgType
      SIMD_STATIC_ASSERT(sizeof(ProcType) >= sizeof(ImgType));

      const int numElems = NUM_SIMDVECS_ELEMENTS(ImgType,ProcType,SIMD_WIDTH);
      const size_t size = in.size;
      
      const int numInVecs = NUM_INPUT_SIMDVECS(ProcType,ImgType);
      const int numOutVecs = NUM_OUTPUT_SIMDVECS(ProcType,ImgType);
      
      SIMDVec<ImgType,SIMD_WIDTH> inVecs[numInVecs];
      SIMDVec<ProcType, SIMD_WIDTH> convertedInVecs[numOutVecs];
      SIMDVec<ProcType,SIMD_WIDTH> outVecs[numOutVecs];
      // the invalid value is needed in ProcType for comparison
      SIMDVec<ProcType,SIMD_WIDTH> maxVecIn;
      SIMDVec<ProcType,SIMD_WIDTH> maxVecOut;
      SIMDVec<ProcType, SIMD_WIDTH> cmp;
      
      maxVecIn = set1<ProcType,SIMD_WIDTH>
	((ProcType)(SIMDTypeInfo<ImgType>::max()));
      maxVecOut = set1<ProcType,SIMD_WIDTH>
	(SIMDTypeInfo<ProcType>::max());
      out.resize(in.w, in.h, in.addOn);
      
      ProcType *dst = out.data;
      
      for (size_t i = 0; i < size; i += numElems) {
	inVecs[0] = load<SIMD_WIDTH>(in.data+i);
	convert(inVecs, convertedInVecs);
	for (int j = 0; j < numOutVecs; j++) {
            cmp = cmpeq(convertedInVecs[j], maxVecIn);
            outVecs[j] = ifelse(cmp, maxVecOut, convertedInVecs[j]);	    
	}
	store(dst, outVecs, numOutVecs);
	dst+=numElems;
      }
    }

    // magnify and scale image
    virtual void
    magnifyImage(const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
		 int interpolation,
		 double magnifyScale,
		 double multiplyScale,
		 SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outImg)
    {
      magnifyAndScaleTilt(img,
			  interpolation, magnifyScale, multiplyScale,
			  outImg);
    }
    
    // copy and scale image
    virtual void
    copyImage(const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
	      double multiplyScale,
	      SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outImg)
    {
      copyAndScaleTilt(img, multiplyScale, outImg);
    }
  };

  // ===========================================================================
  // WarpingSPSComputationAbsTilt
  // contributed by Annika Hoffmann
  // ===========================================================================

  template <typename ImgType,
	    typename ProcType, typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationAbsTilt
    : public WarpingSPSComputationTilt<ImgType, 
				    ProcType, MeasType, SPSType,
				    SIMD_WIDTH, SIMD_ALIGN>
  {
  public:
  
    virtual void 
    preprocessing
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &output,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      //  verticalEdgeFilter(input, output);
      // we use the original images here:
      // images are edge-filtered in advance!
      output = input;
    }

    virtual void 
    columnMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &output,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      sumOfAbsAccuInvalid(input, output);
    }

    virtual std::vector<double>
    maxDenom
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img)
    {
      std::vector<double> md(1);
      // all processing is done in SIMDFloat to avoid overflow
      SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN> inputImg, sumOfAbsImg;
      inputImg = img;
      // no edge filtering for tilted images, we use bg-edges which
      // are computed separately
      // verticalEdgeFilter(inputImg, edgeImg);
      sumOfAbsAccuInvalid(inputImg, sumOfAbsImg);
      md[0] = (2.0 * max(sumOfAbsImg));
      return md;
    }

  };

  // ===========================================================================
  // WarpingSPSComputationNSADTilt
  // contributed by Annika Hoffmann
  // ===========================================================================

  // NSAD with edge filter on tilted image, which include invalid pixels

    template <typename ImgType, 
	    typename ProcType, typename MeasType, typename SPSType, 
	    int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationNSADTilt
    :  public WarpingSPSComputationAbsTilt<ImgType, 
					   ProcType, MeasType, SPSType,
					   SIMD_WIDTH, SIMD_ALIGN>
  {
  public:
  
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &output,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas12,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // with NSAD match operator
      columnPairMatchingTilt1<MeasType,SIMD_WIDTH,SIMD_ALIGN,
			      MatchAbsDiff<MeasType,SIMD_WIDTH> >
	(input1, input2, output,
	 colMeasRaw1, colMeasRaw2,
	 colMeas1, colMeas2, colMeas12);
    }

    virtual void
    normalization
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &combinedMeasure,
     double postScale,
     SIMDImage<SPSType,SIMD_WIDTH,SIMD_ALIGN> &scalePlane,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      // here we only use the combined measure 
      computeScalePlaneNormalizedAdd1Tilt<false>
	(jointMeasure, combinedMeasure, 
	 par.wSPS,
	 par.modulo2wpw, par.shuffle, 
	 postScale, scalePlane);	
    }

    // NSAD measure is in [0,1]
    virtual double 
    measureMax() const { return 1.0; }
    
  };

  // ===========================================================================
  // HarrisParameter
  // contributed by Benedikt Volkmer
  // ===========================================================================

  class HarrisParameter
  {
  public:
    HarrisParameter(double k, int binomialFilterCount) :
      k(k), binomialFilterCount(binomialFilterCount) {}
    double k;
    int binomialFilterCount;
  };
  
  // ===========================================================================
  // WarpingSPSComputation1Alt
  // contributed by Benedikt Volkmer
  // ===========================================================================

  // only defines method computeSPS
  template <typename ImgType, typename ProcType,
	    typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN,
	    bool MAGNIFY_FIRST, bool SCALE_MAGNIFY>
  class WarpingSPSComputation1Alt :
    public WarpingSPSComputationSingle<ImgType, ProcType,
				       MeasType, SPSType,
				       SIMD_WIDTH, SIMD_ALIGN>
  {
  public:
    
    virtual void
    computeSPS
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &ss,
     const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &cv,
     const std::vector<double> &scaleFactors,
     int interpolation,
     const std::vector<double> &pixelScale,
     double postScale,
     WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps)
    {
      if (pixelScale.size() != 1)
	throw SIMDException("WarpingSPSComputation1Alt::computeSPS",
			    "invalid size of pixelScale");
      const int w = sps.param.w, h = ss.h;
      if ((h != cv.h) || (ss.w != cv.w) || (w != ss.w)) {
	char errStr[256];
	sprintf(errStr,
                "mismatch in ss, cv, sps dimensions: "
                "ss %d x %d, cv %d x %d, sps.param.w = %d",
                ss.w, ss.h, cv.w, cv.h, w);
	throw SIMDException("WarpingSPSComputation1Alt::computeSPS",
			    errStr);
      }
      if (ss.addOn != cv.addOn) {
	char errStr[256];
	sprintf(errStr,
                "panorama parameters of ss (%g, %g) and cv (%g, %g) differ",
                ss.addOn.verticalResolution,
                ss.addOn.horizon,
                cv.addOn.verticalResolution,
                cv.addOn.horizon);
	throw SIMDException("WarpingSPSComputation1Alt::computeSPS",
			    errStr);
      }
      // check validity for given image width
      if (!this->validForWidth(w))
	throw SIMDException("WarpingSPSComputation1Alt::computeSPS",
			    "computation not possible for given image width");
      // get input images
      // convert from input type to preprocessing type
      SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> ssOrigF, cvOrigF;
      this->convertImage(ss, ssOrigF);
      this->convertImage(cv, cvOrigF);
      // preprocess input images
      SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> ssF, cvF;
      this->preprocessing(ssOrigF, ssF, sps.param);
      this->preprocessing(cvOrigF, cvF, sps.param);
      // we first process the unmagnified images since the results are needed
      // multiple times
      // 2. Mar 18 (rm): *SR images are the re-arranged versions,
      // earlier in-place re-arrange required heap operations (width differs)
      SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> ssS, cvS, magS,
	ssSR, cvSR, magSR, *cvSPtr, *ssSPtr;
      SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> ssMeasureS, cvMeasureS,
	cvMeasFinal, ssMeasFinal, combMeas,
	magMeasureS, jointMeas, *cvMeasureSPtr, *ssMeasureSPtr;
      this->copyImage(cvF, pixelScale[0], cvS);
      this->reArrange(cvS, cvSR, sps.param);
      this->columnMeasure(cvSR, cvMeasureS, sps.param);
      this->copyImage(ssF, pixelScale[0], ssS);
      this->reArrange(ssS, ssSR, sps.param);
      this->columnMeasure(ssSR, ssMeasureS, sps.param);
      // create all scale plane images
      for (unsigned int spIndex = 0; spIndex < scaleFactors.size(); spIndex++) {
	
	// edge images are scaled
	double scale = scaleFactors[spIndex];
	
	// printf("spIndex %u, scale %g\n", spIndex, scale);
	if (scale > 1.0) {
	  
	  // ===== magnify cv by 1/scale and use unmodified ss =====
	  if (MAGNIFY_FIRST) {
	    SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> magF, magFPre;
	    if (SCALE_MAGNIFY) {
	      magnifyAndScale(cvOrigF, interpolation, 1.0 / scale, scale,
			      magF);
	    }
	    else {
	      magnifyAndScale(cvOrigF, interpolation, 1.0 / scale, 1.0,
			      magF);
	    }
	    this->preprocessing(magF, magFPre, sps.param);
	    this->copyImage(magFPre, pixelScale[0], magS);
	  }
	  else {
	    this->magnifyImage(cvF,
			       interpolation, 1.0 / scale,
			       pixelScale[0], magS);
	  }
	  this->reArrange(magS, magSR, sps.param);
	  cvSPtr = &magSR;
	  this->columnMeasure(magSR, magMeasureS, sps.param);
	  cvMeasureSPtr = &magMeasureS;
	  ssSPtr = &ssSR;
	  ssMeasureSPtr = &ssMeasureS;
	  
	} else if (scale < 1.0) {
	  
	  // ===== magnify ss by scale and use unmodified cv =====
	  cvSPtr = &cvSR;
	  cvMeasureSPtr = &cvMeasureS;
	  if (MAGNIFY_FIRST) {
	    SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> magF, magFPre;
	    if (SCALE_MAGNIFY) {
	      magnifyAndScale(ssOrigF, interpolation, scale, 1/scale, magF);
	    }
	    else {
	      magnifyAndScale(ssOrigF, interpolation, scale, 1, magF);
	    }
	    this->preprocessing(magF, magFPre, sps.param);
	    this->copyImage(magFPre, pixelScale[0], magS);
	  }
	  else {
	    this->magnifyImage(ssF,
			       interpolation, scale,
			       pixelScale[0], magS);
	  }
	  this->reArrange(magS, magSR, sps.param);
	  ssSPtr = &magSR;
	  this->columnMeasure(magSR, magMeasureS, sps.param);
	  ssMeasureSPtr = &magMeasureS;
	  
	} else {
	  
	  // ===== scale = 1.0: use unmodified ss and cv =====
	  cvSPtr = &cvSR;
	  cvMeasureSPtr = &cvMeasureS;
	  ssSPtr = &ssSR;
	  ssMeasureSPtr = &ssMeasureS;
	}
	// store input images if requested
	if (this->storeSPSInput) {
	  // 27. Jul 18 (rm)
	  this->unReArrange(*ssSPtr, this->ssInputs[spIndex], sps.param);
	  this->unReArrange(*cvSPtr, this->cvInputs[spIndex], sps.param);
	  // this->ssInputs[spIndex] = *ssSPtr;
	  // this->cvInputs[spIndex] = *cvSPtr;
	}
	// compute distance term
	this->jointMeasure(*cvSPtr, *ssSPtr,
			   *cvMeasureSPtr, *ssMeasureSPtr,
			   jointMeas, cvMeasFinal, ssMeasFinal, combMeas,
			   sps.param);
	// compute normalized measure
	this->normalization(jointMeas, ssMeasFinal, cvMeasFinal, combMeas,
			    postScale,
			    sps.stack[spIndex], sps.param);
	
      }
      sps.complement();
    }
  };
  
  // ===========================================================================
  // WarpingSPSComputation2Alt
  // contributed by Benedikt Volkmer
  // ===========================================================================

  // defines convertImage, magnifyImage, copyImage
  template <typename ImgType, typename ProcType,
	    typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN,
	    bool MAGNIFY_FIRST, bool SCALE_MAGNIFY>
  class WarpingSPSComputation2Alt
    : public WarpingSPSComputation1Alt<ImgType, ProcType,
				       MeasType, SPSType,
				       SIMD_WIDTH, SIMD_ALIGN,
				       MAGNIFY_FIRST, SCALE_MAGNIFY>
  {
    public:
    // convert image
    virtual void
    convertImage
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &in,
     SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &out)
    {
      // TODO: this assignment is based on type casts, is that always ok?
      // TODO: what is actually happening between signed/unsigned types and
      // TODO: between integer/float types?
      out = in;
    }

    // magnify and scale image
    virtual void
    magnifyImage
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     int interpolation,
     double magnifyScale,
     double multiplyScale,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outImg)
    {
      magnifyAndScale(img,
		      interpolation, magnifyScale, multiplyScale,
		      outImg);
    }

    // copy and scale image
    virtual void
    copyImage
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     double multiplyScale,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outImg)
    {
      copyAndScale(img, multiplyScale, outImg);
    }
  };

  // ===========================================================================
  // WarpingSPSComputationCornerAbs
  // contributed by Benedikt Volkmer
  // ===========================================================================

  // only defines methods preprocessing, columnMeasure,
  // maxDenom

  // ProcType==SIMDFloat (omitted)
  template <typename ImgType,
	    typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN,
	    bool MAGNIFY_FIRST, bool SCALE_MAGNIFY,
	    int NORMALIZE, bool VARIANT>
  class WarpingSPSComputationCornerAbs
    : public WarpingSPSComputation2Alt<ImgType, SIMDFloat,
				       MeasType, SPSType,
				       SIMD_WIDTH, SIMD_ALIGN,
				       MAGNIFY_FIRST, SCALE_MAGNIFY>
  {
    public:

    //ctor
    WarpingSPSComputationCornerAbs(const HarrisParameter &param):
      k(param.k), binomialFilterCount(param.binomialFilterCount) {}

    virtual void
    preprocessing
    (const SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN,Panorama> &jointMeas,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      harrisFilterCyclic<NORMALIZE, VARIANT>
	(input, jointMeas, k, binomialFilterCount);
    }
    
    virtual void
    columnMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      sumOfAbs(input, jointMeas);
    }
    
    virtual std::vector<double>
    maxDenom
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img)
    {
      std::vector<double> md(1);
      // all processing is done in SIMDFloat to avoid overflow
      SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN,Panorama> inputImg, cornerImg;
      SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN> sumOfAbsImg;
      inputImg = img;
      harrisFilterCyclic<NORMALIZE, VARIANT>
	(inputImg, cornerImg, k, binomialFilterCount);
      sumOfAbs(cornerImg, sumOfAbsImg);
      // 2.0* : two sums are added in the denominator
      md[0] = (2.0 * max(sumOfAbsImg));
      return md;
    }
    
  private:
    SIMDFloat k;
    int binomialFilterCount;
  };

  // ===========================================================================
  // measures with reArrange and block-wise rewrite to SPS (2)
  // contributed by Benedikt Volkmer
  // ===========================================================================

  // with reArrange
  // ProcType==SIMDFloat (omitted)
  template <typename ImgType,
	    typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN,
	    int INVERT,
	    bool MAGNIFY_FIRST, bool SCALE_MAGNIFY,
	    int NORMALIZE, bool VARIANT>
  class WarpingSPSComputationCornerAbs2
    :  public WarpingSPSComputationCornerAbs<ImgType,
					     MeasType, SPSType,
					     SIMD_WIDTH, SIMD_ALIGN,
					     MAGNIFY_FIRST, SCALE_MAGNIFY,
					     NORMALIZE, VARIANT>
  {
  public:
    
    //ctor
    WarpingSPSComputationCornerAbs2(const HarrisParameter &param)
      : WarpingSPSComputationCornerAbs
	<ImgType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN, MAGNIFY_FIRST,
	 SCALE_MAGNIFY, NORMALIZE, VARIANT>(param) {}
    
    virtual void
    reArrange
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &imgR,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      swizzleCyclic(img, par.stepAlpha,
		    NUM_INPUT_SIMDVECS(SPSType,MeasType),
		    imgR);
    }

    // 9. Aug 18 (rm): added unReArrange
    virtual void
    unReArrange
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &imgR,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      // serial, slow!
      unSwizzleCyclic(imgR, par.stepAlpha,
		      NUM_INPUT_SIMDVECS(SPSType,MeasType),
		      img);
    }
    
    virtual void
    normalization
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &ssMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &cvMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     double postScale,
     SIMDImage<SPSType,SIMD_WIDTH,SIMD_ALIGN> &scalePlane,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      // INVERT: false: without inversion (distance measure), true:
      // with inversion (correlation)
      computeScalePlaneNormalizedAdd2<INVERT>
	(jointMeasure, ssMeasure, cvMeasure,
	 par.w, par.wSPS, par.stepAlpha, par.nAlpha, par.blkSize,
	 par.modNAlpha,
	 par.shuffleModulo, par.unshuffle,
	 postScale, scalePlane);
    }
    
  };

  // ===========================================================================
  // NSAD with corner filter (derived from WarpingSPSComputationCornerAbs2)
  // contributed by Benedikt Volkmer
  // ===========================================================================

  // with reArrange
  // ProcType==SIMDFloat (omitted)
  template <typename ImgType,
	    typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN,
	    bool MAGNIFY_FIRST, bool SCALE_MAGNIFY,
	    int NORMALIZE, bool VARIANT>
  class WarpingSPSComputationCornerNSAD2
    :  public WarpingSPSComputationCornerAbs2<ImgType,
					      MeasType, SPSType,
					      SIMD_WIDTH, SIMD_ALIGN,
					      /*INVERT*/ false,
					      MAGNIFY_FIRST, SCALE_MAGNIFY,
					      NORMALIZE, VARIANT>
  {
    public:
      //ctor
      WarpingSPSComputationCornerNSAD2(const HarrisParameter &param)
        : WarpingSPSComputationCornerAbs2<ImgType, MeasType, SPSType,
					  SIMD_WIDTH, SIMD_ALIGN,
					  /*INVERT*/ false,
					  MAGNIFY_FIRST, SCALE_MAGNIFY,
					  NORMALIZE, VARIANT>(param) {}

      virtual void
      jointMeasure
      (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
       const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
       const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
       const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
       SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
       SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
       SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
       SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
       const WarpingParameter<SPSType,SIMD_WIDTH> &par)
      {
	/*
	  printf("input1 %d x %d (%g, %g)\n", input1.w, input1.h,
	  input1.addOn.verticalResolution, input1.addOn.horizon);
	  printf("input2 %d x %d (%g, %g)\n", input2.w, input2.h,
	  input2.addOn.verticalResolution, input2.addOn.horizon);
	*/
	
	// puts("NSAD");
	normSumAbsDiffTerm2
	  (input1, input2,
	   par.w, par.wSPS, par.stepAlpha, par.nAlpha, par.blkSize,
	   par.modNAlpha, jointMeas);
	// colMeas1 = colMeasRaw1;
	// colMeas2 = colMeasRaw2;
	colMeas1.view(colMeasRaw1);
	colMeas2.view(colMeasRaw2);
      }
    
    // NSAD measure is in [0,1]
    virtual double
    measureMax() const { return 1.0; }

  };

  // ===========================================================================
  // WarpingSPSComputationEdgeAbs
  // contributed by Benedikt Volkmer
  // ===========================================================================

  // only defines methods preprocessing, columnMeasure,
  // maxDenom
  template <typename ImgType, typename ProcType,
	    typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN,
	    bool MAGNIFY_FIRST, bool SCALE_MAGNIFY>
  class WarpingSPSComputationEdgeAbsAlt
    : public WarpingSPSComputation2Alt<ImgType, ProcType,
				       MeasType, SPSType,
				       SIMD_WIDTH, SIMD_ALIGN, MAGNIFY_FIRST,
				       SCALE_MAGNIFY>
  {
  public:
  
    virtual void 
    preprocessing
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &jointMeas,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      verticalEdgeFilter(input, jointMeas);
    }
    
    virtual void 
    columnMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      sumOfAbs(input, jointMeas);
    }
    
    virtual std::vector<double>
    maxDenom
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img)
    {
      std::vector<double> md(1);
      // all processing is done in SIMDFloat to avoid overflow
      SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN,Panorama> inputImg, edgeImg;
      SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN> sumOfAbsImg;
      inputImg = img;
      verticalEdgeFilter(inputImg, edgeImg);
      sumOfAbs(edgeImg, sumOfAbsImg);
      // 2.0* : two sums are added in the denominator
      md[0] = (2.0 * max(sumOfAbsImg));
      return md;
    }

  };

  // ===========================================================================
  // measures with reArrange and block-wise rewrite to SPS (2)
  // contributed by Benedikt Volkmer
  // ===========================================================================

  // with reArrange
  template <typename ImgType, typename ProcType,
	    typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN,
	    int INVERT, bool MAGNIFY_FIRST,
	    bool SCALE_MAGNIFY>
  class WarpingSPSComputationEdgeAbs2Alt
    :  public WarpingSPSComputationEdgeAbsAlt<ImgType, ProcType,
					      MeasType, SPSType,
					      SIMD_WIDTH, SIMD_ALIGN,
					      MAGNIFY_FIRST, SCALE_MAGNIFY>
  {
  public:
    
    virtual void
    reArrange
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &imgR,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      swizzleCyclic(img, par.stepAlpha,
		    NUM_INPUT_SIMDVECS(SPSType,MeasType),
		    imgR);
    }
  
    virtual void
    normalization
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &ssMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &cvMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     double postScale,
     SIMDImage<SPSType,SIMD_WIDTH,SIMD_ALIGN> &scalePlane,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      // INVERT: false: without inversion (distance measure), true:
      // with inversion (correlation)
      computeScalePlaneNormalizedAdd2<INVERT>
	(jointMeasure, ssMeasure, cvMeasure, 
	 par.w, par.wSPS, par.stepAlpha, par.nAlpha, par.blkSize,
	 par.modNAlpha,
	 par.shuffleModulo, par.unshuffle,
	 postScale, scalePlane);
    }
    
  };

  // ===========================================================================
  // NSAD with edge filter (derived from WarpingSPSComputationEdgeAbs2)
  // contributed by Benedikt Volkmer
  // ===========================================================================

  // with reArrange
  template <typename ImgType, typename ProcType,
	    typename MeasType, typename SPSType,
	    int SIMD_WIDTH, int SIMD_ALIGN,
	    bool MAGNIFY_FIRST, bool SCALE_MAGNIFY>
  class WarpingSPSComputationEdgeNSAD2Alt
    :  public WarpingSPSComputationEdgeAbs2Alt<ImgType, ProcType,
					       MeasType, SPSType,
					       SIMD_WIDTH, SIMD_ALIGN, false,
					       MAGNIFY_FIRST, SCALE_MAGNIFY>
  {
  public:
  
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      /*
      printf("input1 %d x %d (%g, %g)\n", input1.w, input1.h,
	     input1.addOn.verticalResolution, input1.addOn.horizon);
      printf("input2 %d x %d (%g, %g)\n", input2.w, input2.h,
	     input2.addOn.verticalResolution, input2.addOn.horizon);
      */
      
      // puts("NSAD");
      normSumAbsDiffTerm2
	(input1, input2, 
	 par.w, par.wSPS, par.stepAlpha, par.nAlpha, par.blkSize,
	 par.modNAlpha, jointMeas);
      // colMeas1 = colMeasRaw1;
      // colMeas2 = colMeasRaw2;
      colMeas1.view(colMeasRaw1);
      colMeas2.view(colMeasRaw2);
    }

    // NSAD measure is in [0,1]
    virtual double 
    measureMax() const { return 1.0; }

  };
  
  //===========================================================================
  // Multi-Dimensional Sigmoid applied to NSADs of multiple preprocessed images
  // - magnify first, preprocess second
  //===========================================================================

  template <typename ImgType, typename ProcType, typename MeasType,
            typename SPSType, int SIMD_WIDTH, int SIMD_ALIGN, int NPREPROCS>
  class WarpingSPSComputationSNSAD
    :   public WarpingSPSComputationStack<ImgType, ProcType, MeasType, SPSType,
                                          NPREPROCS,
                                          SIMD_WIDTH, SIMD_ALIGN>
  {
    public:
      //ctor
      WarpingSPSComputationSNSAD
      (const std::vector<double> &sigmoidW,
       const std::vector<double> &sigmoidW0)
      {
	if (sigmoidW.size() != (size_t) NPREPROCS)
	  throw SIMDException("WarpingSPSComputationSNSAD",
			      "invalid size of sigmoidW");
	if (sigmoidW0.size() != (size_t) NPREPROCS)
	  throw SIMDException("WarpingSPSComputationSNSAD",
			      "invalid size of sigmoidW0");
	this->sigmoidW = sigmoidW;
	this->sigmoidW0 = sigmoidW0;
      };

      // computation of scale-plane stack 
      // (not possible to restrict ImgType to this function since
      // templates may not be virtual)
      void
      computeSPS
      (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &ss,
       const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &cv,
       const std::vector<double> &scaleFactors,
       int interpolation,
       const std::vector<double> &pixelScale,
       double postScale,
       WarpingSPS<SPSType,SIMD_WIDTH,SIMD_ALIGN> &sps)
      {
	if (pixelScale.size() != (size_t) NPREPROCS)
          throw SIMDException("WarpingSPSComputationSNSAD::computeSPS",
			      "pixelScale array has wrong size");
        const int w = sps.param.w, h = ss.h;
        if ((h != cv.h) || (ss.w != cv.w) || (w != ss.w)) {
          char errStr[256];
          sprintf(errStr,
                  "mismatch in ss, cv, sps dimensions: "
                  "ss %d x %d, cv %d x %d, sps.param.w = %d",
                  ss.w, ss.h, cv.w, cv.h, w);
          throw SIMDException("WarpingSPSComputationSNSAD::computeSPS",
                              errStr);
        }
        if (ss.addOn != cv.addOn) {
          char errStr[256];
          sprintf(errStr,
                  "panorama parameters of ss (%g, %g) and cv (%g, %g) differ",
                  ss.addOn.verticalResolution,
                  ss.addOn.horizon,
                  cv.addOn.verticalResolution,
                  cv.addOn.horizon);
          throw SIMDException("WarpingSPSComputationSNSAD::computeSPS",
                              errStr);
        }
        // check validity for given image width
        if (!this->validForWidth(w))
          throw SIMDException("WarpingSPSComputationSNSAD::computeSPS",
                              "computation not possible for given image width");
        // get input images
        // convert from input type to preprocessing type
        SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> ssOrigF, cvOrigF;
        this->convertImage(ss, ssOrigF);
        this->convertImage(cv, cvOrigF);
        // preprocess input images
        SIMDStack<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> ssF, cvF;
        this->preprocessing(ssOrigF, ssF, sps.param);
        this->preprocessing(cvOrigF, cvF, sps.param);
        // we first process the unmagnified images since the results are needed
        // multiple times
        // 2. Mar 18 (rm): *SR images are the re-arranged versions,
        // earlier in-place re-arrange required heap operations (width differs)
        SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> ssS, cvS, magS,
          ssSR, cvSR, magSR, *cvSPtr, *ssSPtr;
        SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> ssMeasureS, cvMeasureS,
          cvMeasFinal, ssMeasFinal, combMeas,
          magMeasureS, jointMeas, *cvMeasureSPtr, *ssMeasureSPtr;
        this->copyStack(cvF, pixelScale, cvS);
        this->reArrange(cvS, cvSR, sps.param);
        this->columnMeasure(cvSR, cvMeasureS, sps.param);
        this->copyStack(ssF, pixelScale, ssS);
        this->reArrange(ssS, ssSR, sps.param);
        this->columnMeasure(ssSR, ssMeasureS, sps.param);


        // create all scale plane images
        for (unsigned int spIndex = 0; spIndex < scaleFactors.size(); spIndex++) {
          
          // edge images are scaled
          double scale = scaleFactors[spIndex];
          
          // printf("spIndex %u, scale %g\n", spIndex, scale);
          if (scale > 1.0) {
            // ===== magnify cv by 1/scale and use unmodified ss =====
            SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> magF;
            SIMDStack<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> magFPre;
            magnifyAndScale(cvOrigF, interpolation, 1.0 / scale, 1.0,
                            magF);
            this->preprocessing(magF, magFPre, sps.param);
            this->copyStack(magFPre, pixelScale, magS);
            this->reArrange(magS, magSR, sps.param);
            cvSPtr = &magSR;
            this->columnMeasure(magSR, magMeasureS, sps.param);
            cvMeasureSPtr = &magMeasureS;
            ssSPtr = &ssSR;
            ssMeasureSPtr = &ssMeasureS;
            
          } else if (scale < 1.0) {
            
            // ===== magnify ss by scale and use unmodified cv =====
            cvSPtr = &cvSR;
            cvMeasureSPtr = &cvMeasureS;
            SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> magF;
            SIMDStack<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> magFPre;
            magnifyAndScale(ssOrigF, interpolation, scale, 1, magF);
            this->preprocessing(magF, magFPre, sps.param);
            this->copyStack(magFPre, pixelScale, magS);
            this->reArrange(magS, magSR, sps.param);
            ssSPtr = &magSR;
            this->columnMeasure(magSR, magMeasureS, sps.param);
            ssMeasureSPtr = &magMeasureS;
            
          } else {
            
            // ===== scale = 1.0: use unmodified ss and cv =====
            cvSPtr = &cvSR;
            cvMeasureSPtr = &cvMeasureS;
            ssSPtr = &ssSR;
            ssMeasureSPtr = &ssMeasureS;
          }
          // store input images if requested
          if (this->storeSPSInput) {
            // 27. Jul 18 (rm)
            this->unReArrange(*ssSPtr, this->ssInputs[spIndex], sps.param);
            this->unReArrange(*cvSPtr, this->cvInputs[spIndex], sps.param);
            // this->ssInputs[spIndex] = *ssSPtr;
            // this->cvInputs[spIndex] = *cvSPtr;
          }
          // compute distance term
          this->jointMeasure(*cvSPtr, *ssSPtr, jointMeas, sps.param);
          // compute normalized measure
          this->normalization(jointMeas, *ssMeasureSPtr, *cvMeasureSPtr,
                              combMeas, postScale, sps.stack[spIndex],
                              sps.param);
          
        }
        sps.complement();
      }

      double
      measureMax() const { return 1.0; }

      void
      convertImage
      (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &in,
       SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &out)
      {
        //TODO: see WarpingSPSComputation2
        out = in;
      }

      void
      copyImage
      (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img,
       double multiplyScale,
       SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outImg)
      {
        copyAndScale(img, multiplyScale, outImg);
      }

      void
      copyStack
      (const SIMDStack<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &stk,
       const std::vector<double> &multiplyScale,
       SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &outStk)
      {
        copyAndScale(stk, multiplyScale, outStk);
      }

      void 
      columnMeasure
      (const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &in,
       SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeas,
       const WarpingParameter<SPSType,SIMD_WIDTH> &)
      {
        //jointMeas.h == 1 i.e. one value for each column
        jointMeas.resize(in.numPlanes, in.w, 1);
        // TODO: replace loop with own functions
        for (int i=0; i < in.numPlanes; i++){
          sumOfAbs(in[i], jointMeas[i]);
        }
      }

      void
      jointMeasure
      (const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &inputStk1,
       const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &inputStk2,
       SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeasStk,
       const WarpingParameter<SPSType,SIMD_WIDTH> &)
      {
        //jointMeasStk.h == input.w
        jointMeasStk.resize(inputStk1.numPlanes, inputStk1.w, inputStk1.w);
        normSumAbsDiffTerm1(inputStk1[0], inputStk2[0], jointMeasStk[0]);
        for (int i=0; i < inputStk1.numPlanes; i++)
        {
          // puts("NSAD");
          normSumAbsDiffTerm1(inputStk1[i], inputStk2[i], jointMeasStk[i]);
        }
      }

      // normalization of jointMeasure and ssMeasure/cvMeasure
      void
      normalization
      (const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeasure,
       const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &ssMeasure,
       const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &cvMeasure,
       const SIMDStack<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
       double postScale,
       SIMDImage<SPSType,SIMD_WIDTH,SIMD_ALIGN> &scalePlane,
       const WarpingParameter<SPSType,SIMD_WIDTH> &par)
      {
        // INVERT: false: without inversion (distance measure), true:
        // with inversion (correlation)
        computeScalePlaneSigmoidNormalizedAdd1<false, NPREPROCS>
          (jointMeasure, ssMeasure, cvMeasure, 
           par.wSPS,
           par.modulo2wpw, par.shuffle, 
           sigmoidW, sigmoidW0, postScale, scalePlane);
        //ns_simd::computeScalePlaneNormalizedAdd1<false>
        //  (jointMeasure[0], ssMeasure[0], cvMeasure[0], par.wSPS,
        //   par.modulo2wpw, par.shuffle, postScale, scalePlane);
      }

      static const unsigned int npreprocs = NPREPROCS;
    private:
      std::vector<double> sigmoidW, sigmoidW0;
  };

  template <typename ImgType, typename MeasType, typename SPSType,
            int SIMD_WIDTH, int SIMD_ALIGN, HarrisOutputSelector HOS,
            int OUTPUT_COUNT>
  class WarpingSPSComputationSNSADHarris
    :   public WarpingSPSComputationSNSAD<ImgType, SIMDFloat, MeasType, SPSType,
                                          SIMD_WIDTH, SIMD_ALIGN,
                                          OUTPUT_COUNT>
  {
    public:
      //ctor
      WarpingSPSComputationSNSADHarris
      (const HarrisParameter &param,
       const std::vector<double> &sigmoidW,
       const std::vector<double> &sigmoidW0):
	WarpingSPSComputationSNSAD<ImgType, SIMDFloat, MeasType, SPSType,
				   SIMD_WIDTH, SIMD_ALIGN,
				   OUTPUT_COUNT>
	(sigmoidW, sigmoidW0),
	k(param.k), binomialFilterCount(param.binomialFilterCount){};
    
      // preprocessing applied to the input images
      void
      preprocessing
      (const SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
       SIMDStack<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN,Panorama> &output,
       const WarpingParameter<SPSType,SIMD_WIDTH> &)
      {
        harrisFilterCyclic<2, 0, HOS>(input, output, k, binomialFilterCount);
      }

      std::vector<double>
      maxDenom
      (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img)
      {
	std::vector<double> md(this->npreprocs);
        // all processing is done in SIMDFloat to avoid overflow
        SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN,Panorama> inputImg;
        SIMDStack<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN,Panorama> preprocStk;
        SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN> sumOfAbsImg;
        inputImg = img;
        harrisFilterCyclic<2, 0, HOS>
                          (inputImg, preprocStk, k, binomialFilterCount);
        for (unsigned int pp=0; pp < this->npreprocs; pp++) {
          sumOfAbs(preprocStk[pp], sumOfAbsImg);
          // 2.0* : two sums are added in the denominator
          md[pp] = 2.0 * max(sumOfAbsImg);
        }
	return md;
      }

    private:
      double k;
      int binomialFilterCount;
  };

  //===========================================================================
  // version for tilted images but with vertical edge filter
  // code contributed by C. Berganski (based on code by A. Hoffmann)
  //===========================================================================

  template<typename ImgType, typename ProcType, typename MeasType,
	   typename SPSType, int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeAbsTilt 
    : public WarpingSPSComputationTilt<ImgType, ProcType, MeasType, SPSType,
				       SIMD_WIDTH, SIMD_ALIGN>
  {
  public:
    
    virtual void 
    preprocessing
    (const SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<ProcType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &output,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      verticalEdgeFilterTilt(input, output);
    }
    
    virtual void 
    columnMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &output,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      sumOfAbsAccuInvalid(input, output);
    }
    
    // 29. Sep 19 (rm): modified
    virtual std::vector<double>
    maxDenom
    (const SIMDImage<ImgType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &img)
    {
      std::vector<double> md(1);
      // all processing is done in SIMDFloat to avoid overflow
      SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN,Panorama> inputImg;
      SIMDImage<SIMDFloat,SIMD_WIDTH,SIMD_ALIGN> sumOfAbsImg, edgeImg;
      this->convertImage(img, inputImg);
      verticalEdgeFilter(inputImg, edgeImg);
      sumOfAbsAccuInvalid(edgeImg, sumOfAbsImg);
      md[0] = (2.0 * max(sumOfAbsImg));
      return md;
    }
    
  };

  template<typename ImgType, 
	   typename ProcType, typename MeasType, typename SPSType, 
	   int SIMD_WIDTH, int SIMD_ALIGN>
  class WarpingSPSComputationEdgeNSADTilt
    :  public WarpingSPSComputationEdgeAbsTilt<ImgType, ProcType, MeasType,
					       SPSType, SIMD_WIDTH,
					       SIMD_ALIGN>
  {
  public:
    
    virtual void 
    jointMeasure
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN,Panorama> &input2,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw1,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeasRaw2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &output,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas1,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas2,
     SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &colMeas12,
     const WarpingParameter<SPSType,SIMD_WIDTH> &)
    {
      // with NSAD match operator
      columnPairMatchingTilt1<MeasType,SIMD_WIDTH,SIMD_ALIGN,
			      MatchAbsDiff<MeasType,SIMD_WIDTH> >
        (input1, input2, output,
	 colMeasRaw1, colMeasRaw2,
	 colMeas1, colMeas2, colMeas12);
    }
    
    virtual void
    normalization
    (const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &jointMeasure,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &,
     const SIMDImage<MeasType,SIMD_WIDTH,SIMD_ALIGN> &combinedMeasure,
     double postScale,
     SIMDImage<SPSType,SIMD_WIDTH,SIMD_ALIGN> &scalePlane,
     const WarpingParameter<SPSType,SIMD_WIDTH> &par)
    {
      // here we only use the combined measure 
      computeScalePlaneNormalizedAdd1Tilt<false>
        (jointMeasure, combinedMeasure, 
	 par.wSPS,
	 par.modulo2wpw, par.shuffle, 
	 postScale, scalePlane);	
    }
    
    // NSAD measure is in [0,1]
    virtual double 
    measureMax() const { return 1.0; }
  };
  
}

#endif
