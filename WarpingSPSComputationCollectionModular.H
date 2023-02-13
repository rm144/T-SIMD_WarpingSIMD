// ===========================================================================
//
// WarpingSPSComputationCollectionModular.H
// novel policy-based scheme for WarpingSPSComputation
// Author: Benedikt Volkmer (Bielefeld University)
// Year of creation: 2020
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
// (C) Benedikt Volkmer, Ralf Möller
//     Computer Engineering
//     Faculty of Technology
//     Bielefeld University
//     www.ti.uni-bielefeld.de
//
// ===========================================================================

#ifndef _WARPING_SPS_COMPUTATION_COLLECTION_MODULAR_H_
#define _WARPING_SPS_COMPUTATION_COLLECTION_MODULAR_H_

#include "SIMDImageFunctions.H"
#include "WarpingSPSComputation.H"

namespace simd {
namespace modular_sps_computation {

template <bool B>
struct Rearranged
{};
template <bool B>
struct ConsidersInvalid
{};

// ============================================================================
// Single
// ============================================================================
// Single is the base class for single preprocessing output. It is based on
// former WarpingSPSComputationSingle and WarpingSPSComputation1/Alt. It
// implements all virtual functions of WarpingSPSComputation using static
// functions of policy classes given by template parameters:
//
// - ImageOpsPolicy: (convertImage, copyImage, magnifyImage) (for both single
//                                                            and stack)
//   - Simple (former WarpingSPSComputation2/Alt) (TODO: remove if Tilt is
//                                                 universal)
//   - Invalid (former WarpingSPSComputationTilt)
// - SinglePreprocessingPolicy: (preprocessing)
//   - Edge (former WarpingSPSComputationEdgeAbs/Alt)
//   - Harris (former WarpingSPSComputationHarrisAbs)
//   - NonAccuInvalid (former WarpingSPSComputationAbsTilt)
//   - EdgeAccuInvalid (former WarpingSPSComputationEdgeAbsTilt)
// - SingleColumnMeasurePolicy: (columnMeasure)
//   - SingleColumnMeasurePolicyAbs (former WarpingSPSComputation*Abs)
//   - SingleColumnMeasurePolicyAbsAccu (former WarpingSPSComputation*AbsTilt
// - NormalizePolicy: (normalization, reArrange, unReArrange)
//   - NonRearranged (virtual, included in former WarpingSPSComputationSingle,
//                    reArrange and unReArrange just copy,
//                    normalization virtual, for both stack and single)
//     - Single*Abs (former WarpingSPSComputationEdgeAbs1)
//   - Single*AbsRearranged (former WarpingSPSComputationEdgeAbs2/Alt or
//                    WarpingSPSComputationHarrisAbs2)
// - MeasurePolicy: (jointMeasure, measureMax) (for both single and stack)
//   - NSADA (former WarpingSPSComputationEdgeNSADA1)
//   - NSADNSADA (former WarpingSPSComputationEdgeNSADNSADA1)
//   - ASCDirectRearranged (former WarpingSPSComputationASCDirect2)
//   - NSADRearranged (former WarpingSPSComputationEdgeNSAD2/Alt or
//                     WarpingSPSComputationHarrisNSAD2)
//   - NSADARearranged (former WarpingSPSComputationEdgeNSADA2)
//   - NSADNSADARearranged (former WarpingSPSComputationEdgeNSADNSADA2)
//   - NSADInvalid (former WarpingSPSComputationNSADTilt or
//                               WarpingSPSComputationEdgeNSADTilt)
//
// Example usage: Preprocess first, simple image ops, edge filtered
//                preprocessing, rearranged abs normalization (non inverted),
//                NSAD measure
// ```
// using ns_simd::modular_sps_computation;
// WarpingSPSComputation* spsComp = new Single<
//        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, false,
//        ImageOpsPolicySimple<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
//        SinglePreprocessingPolicyEdge<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
//        SingleNormalizePolicyAbsRearranged<SPS_COMPUTATION_TEMPLATE_ARGUMENTS,
//                                           0>,
//        SingleMeasurePolicyNSADRearranged<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
//        >();
// ```
// TODO: example is outdated, please consult WarpingFactories.H
#define SPS_COMPUTATION_TEMPLATE_PARAMETER                                     \
  typename ImgType, typename ProcType, typename MeasType, typename SPSType,    \
    int SIMD_WIDTH, int SIMD_ALIGN
#define SPS_COMPUTATION_TEMPLATE_ARGUMENTS                                     \
  ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN

template <SPS_COMPUTATION_TEMPLATE_PARAMETER, bool MAGNIFY_FIRST,
          class ImageOpsPolicy, class SinglePreprocessingPolicy,
          class SingleNormalizePolicy, class SingleMeasurePolicy,
          class SingleColumnMeasurePolicy>
class Single : public WarpingSPSComputation<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
{
protected:
  // computeSPS should store input images
  bool storeSPSInput;
  // input images in sps computation
  std::vector<SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama>> ssInputs,
    cvInputs;

  // compute SPS and intermediate input images
  // 27. Feb 18 (rm): vert.res. and horizon taken from images
  void computeSPSInput(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
    const std::vector<double> &scaleFactors, int interpolation,
    const std::vector<double> &pixelScale, double postScale,
    WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps)
  {
    storeSPSInput = true;
    ssInputs.resize(scaleFactors.size());
    cvInputs.resize(scaleFactors.size());
    this->computeSPS(ss, cv, scaleFactors, interpolation, pixelScale, postScale,
                     sps);
    storeSPSInput = false;
  }

public:
  // ======================== ctor, dtor =========================

  Single(typename SinglePreprocessingPolicy::Parameter preproc_param)
    : WarpingSPSComputation<ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH,
                            SIMD_ALIGN>(),
      storeSPSInput(false),
      singlePreprocessingPolicy(SinglePreprocessingPolicy(preproc_param))
  {}

  virtual ~Single() {}

  virtual void computeSPS(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
    const std::vector<double> &scaleFactors, int interpolation,
    const std::vector<double> &pixelScale, double postScale,
    WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps)
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
      throw SIMDException("WarpingSPSComputation1Alt::computeSPS", errStr);
    }
    if (ss.addOn != cv.addOn) {
      char errStr[256];
      sprintf(errStr,
              "panorama parameters of ss (%g, %g) and cv (%g, %g) differ",
              ss.addOn.verticalResolution, ss.addOn.horizon,
              cv.addOn.verticalResolution, cv.addOn.horizon);
      throw SIMDException("WarpingSPSComputation1Alt::computeSPS", errStr);
    }
    // check validity for given image width
    if (!this->validForWidth(w))
      throw SIMDException("WarpingSPSComputation1Alt::computeSPS",
                          "computation not possible for given image width");
    // get input images
    // convert from input type to preprocessing type
    SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> ssOrigF, cvOrigF;
    ImageOpsPolicy::convertImage(ss, ssOrigF);
    ImageOpsPolicy::convertImage(cv, cvOrigF);
    // preprocess input images
    SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> ssF, cvF;
    singlePreprocessingPolicy.preprocessing(considersInvalid, ssOrigF, ssF);
    singlePreprocessingPolicy.preprocessing(considersInvalid, cvOrigF, cvF);
    // we first process the unmagnified images since the results are needed
    // multiple times
    // 2. Mar 18 (rm): *SR images are the re-arranged versions,
    // earlier in-place re-arrange required heap operations (width differs)
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> ssS, cvS, magS, ssSR,
      cvSR, magSR, *cvSPtr, *ssSPtr;
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> ssMeasureS, cvMeasureS,
      cvMeasFinal, ssMeasFinal, combMeas, magMeasureS, jointMeas,
      *cvMeasureSPtr, *ssMeasureSPtr;
    ImageOpsPolicy::copyImage(cvF, pixelScale[0], cvS);
    SingleNormalizePolicy::reArrange(cvS, cvSR, sps.param);
    SingleColumnMeasurePolicy::columnMeasure(considersInvalid, cvSR,
                                             cvMeasureS);
    ImageOpsPolicy::copyImage(ssF, pixelScale[0], ssS);
    SingleNormalizePolicy::reArrange(ssS, ssSR, sps.param);
    SingleColumnMeasurePolicy::columnMeasure(considersInvalid, ssSR,
                                             ssMeasureS);
    // create all scale plane images
    for (unsigned int spIndex = 0; spIndex < scaleFactors.size(); spIndex++) {
      // edge images are scaled
      double scale = scaleFactors[spIndex];

      // printf("spIndex %u, scale %g\n", spIndex, scale);
      if (scale > 1.0) {
        // ===== magnify cv by 1/scale and use unmodified ss =====
        if (MAGNIFY_FIRST) {
          SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> magF, magFPre;
          magnifyAndScale(cvOrigF, interpolation, 1.0 / scale, 1.0, magF);
          singlePreprocessingPolicy.preprocessing(considersInvalid, magF,
                                                  magFPre);
          ImageOpsPolicy::copyImage(magFPre, pixelScale[0], magS);
        } else {
          ImageOpsPolicy::magnifyImage(cvF, interpolation, 1.0 / scale,
                                       pixelScale[0], magS);
        }
        SingleNormalizePolicy::reArrange(magS, magSR, sps.param);
        cvSPtr = &magSR;
        SingleColumnMeasurePolicy::columnMeasure(considersInvalid, magSR,
                                                 magMeasureS);
        cvMeasureSPtr = &magMeasureS;
        ssSPtr        = &ssSR;
        ssMeasureSPtr = &ssMeasureS;

      } else if (scale < 1.0) {
        // ===== magnify ss by scale and use unmodified cv =====
        cvSPtr        = &cvSR;
        cvMeasureSPtr = &cvMeasureS;
        if (MAGNIFY_FIRST) {
          SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> magF, magFPre;
          magnifyAndScale(ssOrigF, interpolation, scale, 1 / scale, magF);
          singlePreprocessingPolicy.preprocessing(considersInvalid, magF,
                                                  magFPre);
          ImageOpsPolicy::copyImage(magFPre, pixelScale[0], magS);
        } else {
          ImageOpsPolicy::magnifyImage(ssF, interpolation, scale, pixelScale[0],
                                       magS);
        }
        SingleNormalizePolicy::reArrange(magS, magSR, sps.param);
        ssSPtr = &magSR;
        SingleColumnMeasurePolicy::columnMeasure(considersInvalid, magSR,
                                                 magMeasureS);
        ssMeasureSPtr = &magMeasureS;

      } else {
        // ===== scale = 1.0: use unmodified ss and cv =====
        cvSPtr        = &cvSR;
        cvMeasureSPtr = &cvMeasureS;
        ssSPtr        = &ssSR;
        ssMeasureSPtr = &ssMeasureS;
      }
      // store input images if requested
      if (this->storeSPSInput) {
        // 27. Jul 18 (rm)
        SingleNormalizePolicy::unReArrange(*ssSPtr, this->ssInputs[spIndex],
                                           sps.param);
        SingleNormalizePolicy::unReArrange(*cvSPtr, this->cvInputs[spIndex],
                                           sps.param);
        // this->ssInputs[spIndex] = *ssSPtr;
        // this->cvInputs[spIndex] = *cvSPtr;
      }
      // compute distance term
      SingleMeasurePolicy::jointMeasure(rearranged, considersInvalid, *cvSPtr,
                                        *ssSPtr, *cvMeasureSPtr, *ssMeasureSPtr,
                                        jointMeas, cvMeasFinal, ssMeasFinal,
                                        combMeas, sps.param);
      // compute normalized measure
      SingleNormalizePolicy::normalization(
        considersInvalid, jointMeas, ssMeasFinal, cvMeasFinal, combMeas,
        postScale, sps.stack[spIndex], sps.param);
    }
    sps.complement();
  }

  std::vector<double> maxDenom(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img)
  {
    std::vector<double> md(1);
    // all processing is done in SIMDFloat to avoid overflow
    SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> inputImg, preprocImg;
    SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN> measImg;
    ImageOpsPolicy::convertImage(img, inputImg);
    this->singlePreprocessingPolicy.preprocessing(considersInvalid, inputImg,
                                                  preprocImg);
    SingleColumnMeasurePolicy::columnMeasure(considersInvalid, preprocImg,
                                             measImg);
    md[0] = (2.0 * max(measImg));
    return md;
  };

  int numPreprocs() const { return 1; };

  // if computed with real numbers, maximal value of the normalized measure
  // (used for WarpingSPSComputation::idealPostScale())
  virtual double measureMax() const
  {
    return SingleMeasurePolicy::measureMax();
  };

  // compute a diagram showing all intermediate input images for SPS
  // 27. Feb 18 (rm): vert.res. and hor. taken from images
  virtual void computeSPSInputDiagram(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
    const std::vector<double> &scaleFactors, int interpolation,
    const std::vector<double> &pixelScale, double postScale,
    WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps,
    // outputs:
    MeasType &allMax, MeasType &allMin, MeasType &allRange,
    SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN> &spsInput)
  {
    computeSPSInput(ss, cv, scaleFactors, interpolation, pixelScale, postScale,
                    sps);
    unsigned nPlanes = scaleFactors.size();
    MeasType invalid = SIMDTypeInfo<MeasType>::max();
    // first: determine max and min over all images
    allMax = SIMDTypeInfo<MeasType>::min();
    allMin = SIMDTypeInfo<MeasType>::max();
    for (unsigned planeIdx = 0; planeIdx < nPlanes; planeIdx++) {
      MeasType _min, _max;
      _max   = ImageOpsPolicy::considersInvalidPixels ?
                 maxExcept(ssInputs[planeIdx], invalid) :
                 max(ssInputs[planeIdx]);
      _min   = min(ssInputs[planeIdx]);
      allMax = std::max(allMax, _max);
      allMin = std::min(allMin, _min);
      _max   = ImageOpsPolicy::considersInvalidPixels ?
                 maxExcept(cvInputs[planeIdx], invalid) :
                 max(cvInputs[planeIdx]);
      _min   = min(cvInputs[planeIdx]);
      allMax = std::max(allMax, _max);
      allMin = std::max(allMin, _min);
    }
    allRange = std::max(allMax, MeasType(-allMin));
    // write into same image
    SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN> tmp;
    int w = ssInputs[0].w, h = ssInputs[0].h;
    const int margin = SIMDVec<SIMDByte, SIMD_WIDTH>::elements;
    spsInput.resize(2 * (w + margin) + margin, nPlanes * (h + margin) + margin);
    spsInput.fill(255);
    double shift = 128.0, scale = 255.0 / allRange;
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> tmpSS, tmpCV;
    for (unsigned planeIdx = 0; planeIdx < nPlanes; planeIdx++) {
      if (ImageOpsPolicy::considersInvalidPixels) {
        replace(ssInputs[planeIdx], invalid, allRange, tmpSS);
        replace(cvInputs[planeIdx], invalid, allRange, tmpCV);
      } else {
        tmpSS.view(ssInputs[planeIdx]);
        tmpCV.view(cvInputs[planeIdx]);
      }
      copyScaleAndShift(tmpSS, scale, shift, tmp);
      insert(tmp, margin, margin + planeIdx * (h + margin), spsInput);
      copyScaleAndShift(tmpCV, scale, shift, tmp);
      insert(tmp, w + 2 * margin, margin + planeIdx * (h + margin), spsInput);
    }
  }

private:
  SinglePreprocessingPolicy singlePreprocessingPolicy;

  // check whether input images with given width could be processed
  // this definition should be useful for many cases:
  virtual bool validForWidth(int width) const
  {
    return
      // ProcType->MeasType: uses magnifyImage, copyImage:
      ((width % numSIMDVecsElements<MeasType, ProcType, SIMD_WIDTH>() == 0) &&
       // MeasType->SPSType: uses computeScalePlaneNormalizedAdd:
       (width % numSIMDVecsElements<SPSType, MeasType, SIMD_WIDTH>() == 0));
  }

  static const Rearranged<SingleNormalizePolicy::isRearranged> rearranged;
  static const ConsidersInvalid<ImageOpsPolicy::considersInvalidPixels>
    considersInvalid;
};

// To prevent clang++ warning for c++98 (no instantiation of variable). Can't be
// in class, because it is static (this is relaxed in c++17). This is meant to
// (most likely) be removed with compiler optimization (empty struct and unused
// variables)
template <SPS_COMPUTATION_TEMPLATE_PARAMETER, bool MAGNIFY_FIRST,
          class ImageOpsPolicy, class SinglePreprocessingPolicy,
          class SingleNormalizePolicy, class SingleMeasurePolicy,
          class SingleColumnMeasurePolicy>
const Rearranged<SingleNormalizePolicy::isRearranged>
  Single<SPS_COMPUTATION_TEMPLATE_ARGUMENTS, MAGNIFY_FIRST, ImageOpsPolicy,
         SinglePreprocessingPolicy, SingleNormalizePolicy, SingleMeasurePolicy,
         SingleColumnMeasurePolicy>::rearranged = {};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER, bool MAGNIFY_FIRST,
          class ImageOpsPolicy, class SinglePreprocessingPolicy,
          class SingleNormalizePolicy, class SingleMeasurePolicy,
          class SingleColumnMeasurePolicy>
const ConsidersInvalid<ImageOpsPolicy::considersInvalidPixels>
  Single<SPS_COMPUTATION_TEMPLATE_ARGUMENTS, MAGNIFY_FIRST, ImageOpsPolicy,
         SinglePreprocessingPolicy, SingleNormalizePolicy, SingleMeasurePolicy,
         SingleColumnMeasurePolicy>::considersInvalid = {};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class ImageOpsPolicySimple
{
public:
  static const bool considersInvalidPixels = false;

  static void convertImage(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &in,
    SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &out)
  {
    // TODO: this assignment is based on type casts, is that always ok?
    // TODO: what is actually happening between signed/unsigned types and
    // TODO: between integer/float types?
    out = in;
  }

  static void magnifyImage(
    const SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img,
    int interpolation, double magnifyScale, double multiplyScale,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &outImg)
  {
    magnifyAndScale(img, interpolation, magnifyScale, multiplyScale, outImg);
  }

  static void copyImage(
    const SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img,
    double multiplyScale,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &outImg)
  {
    copyAndScale(img, multiplyScale, outImg);
  }

  static void copyStack(
    const SIMDStack<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &stk,
    const std::vector<double> &multiplyScale,
    SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &outStk)
  {
    copyAndScale(stk, multiplyScale, outStk);
  }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class SingleImageOpsPolicyInvalid
{
public:
  static const bool considersInvalidPixels = true;

  // convert image: keep values from ImgType, but use invalid value
  // for proc type
  // if value is invalid (max(ImgType)) in input image, then this
  // value is set to max(ProcType)

  static void convertImage(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &in,
    SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &out)
  {
    // code assumes that ProcType is larger or equal in size then ImgType
    SIMD_STATIC_ASSERT(sizeof(ImgType) <= sizeof(ProcType));

    const int numElems = numSIMDVecsElements<ImgType, ProcType, SIMD_WIDTH>();
    const size_t size  = in.size;

    const int numInVecs  = numInputSIMDVecs<ProcType, ImgType>();
    const int numOutVecs = numOutputSIMDVecs<ProcType, ImgType>();

    SIMDVec<ImgType, SIMD_WIDTH> inVecs[numInVecs];
    SIMDVec<ProcType, SIMD_WIDTH> convertedInVecs[numOutVecs];
    SIMDVec<ProcType, SIMD_WIDTH> outVecs[numOutVecs];
    // the invalid value is needed in ProcType for comparison
    SIMDVec<ProcType, SIMD_WIDTH> maxVecIn;
    SIMDVec<ProcType, SIMD_WIDTH> maxVecOut;
    SIMDVec<ProcType, SIMD_WIDTH> cmp;

    maxVecIn =
      set1<ProcType, SIMD_WIDTH>((ProcType) (SIMDTypeInfo<ImgType>::max()));
    maxVecOut = set1<ProcType, SIMD_WIDTH>(SIMDTypeInfo<ProcType>::max());
    out.resize(in.w, in.h, in.addOn);

    ProcType *dst = out.data;

    for (size_t i = 0; i < size; i += numElems) {
      inVecs[0] = load<SIMD_WIDTH>(in.data + i);
      convert(inVecs, convertedInVecs);
      for (int j = 0; j < numOutVecs; j++) {
        cmp        = cmpeq(convertedInVecs[j], maxVecIn);
        outVecs[j] = ifelse(cmp, maxVecOut, convertedInVecs[j]);
      }
      store(dst, outVecs, numOutVecs);
      dst += numElems;
    }
  }

  // magnify and scale image
  static void magnifyImage(
    const SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img,
    int interpolation, double magnifyScale, double multiplyScale,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &outImg)
  {
    magnifyAndScaleTilt(img, interpolation, magnifyScale, multiplyScale,
                        outImg);
  }

  // copy and scale image
  static void copyImage(
    const SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img,
    double multiplyScale,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &outImg)
  {
    copyAndScaleTilt(img, multiplyScale, outImg);
  }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class SinglePreprocessingPolicyEdge
{
public:
  class Parameter
  {};

  // ctor
  SinglePreprocessingPolicyEdge(SinglePreprocessingPolicyEdge::Parameter) {};

  static void preprocessing(
    ConsidersInvalid<false>,
    const SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input,
    SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &preprocessed)
  {
    verticalEdgeFilter(input, preprocessed);
  }

  static void preprocessing(
    ConsidersInvalid<true>,
    const SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input,
    SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &output)
  {
    verticalEdgeFilterTilt(input, output);
  }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class SinglePreprocessingPolicyNon
{
public:
  class Parameter
  {};

  // ctor
  SinglePreprocessingPolicyNon(SinglePreprocessingPolicyNon::Parameter) {};

  template <bool REARRANGED>
  static void preprocessing(
    ConsidersInvalid<REARRANGED>,
    const SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input,
    SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &output)
  {
    //  verticalEdgeFilter(input, output);
    // we use the original images here:
    // images are edge-filtered in advance!
    output = input;
  }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER, int NORMALIZE, bool VARIANT,
          HarrisOutputSelector HOS>
class SinglePreprocessingPolicyHarris
{
public:
  class Parameter
  {
  public:
    Parameter(double k, int binomialFilterCount)
      : k(k), binomialFilterCount(binomialFilterCount)
    {}
    double k;
    int binomialFilterCount;
  };

  // ctor
  SinglePreprocessingPolicyHarris(
    SinglePreprocessingPolicyHarris::Parameter param)
    : k(param.k), binomialFilterCount(param.binomialFilterCount) {};

  void preprocessing(
    ConsidersInvalid<false>,
    const SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input,
    SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &jointMeas)
  {
    SIMDStack<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN, Panorama> jointMeasStk(
      1, input.w, input.h);
    harrisFilterCyclic<NORMALIZE, VARIANT, HOS>(input, jointMeasStk, this->k,
                                                this->binomialFilterCount);
    jointMeas.view(jointMeasStk[0]);
  }

private:
  double k;
  int binomialFilterCount;
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class SingleColumnMeasurePolicyAbs
{
public:
  template <typename T>
  static void columnMeasure(
    ConsidersInvalid<false>,
    const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input,
    SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &colMeas)
  {
    sumOfAbs(input, colMeas);
  }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class SingleColumnMeasurePolicyAbsAccu
{
public:
  template <typename T>
  static void columnMeasure(
    ConsidersInvalid<true>,
    const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input,
    SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &colMeas)
  {
    sumOfAbsAccuInvalid(input, colMeas);
  }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class NormalizePolicyNonRearranged
{
public:
  static const bool isRearranged = false;

  // re-arrange image (makeView instead of pure virtual)
  static void reArrange(
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &imgR,
    const WarpingParameter<SPSType, SIMD_WIDTH> &)
  {
    // TODO: this is fast but may be dangerous, true copy may be advisable
    imgR.view(img);
    // imgR = img;
  }

  // un-re-arrange image (copy instead of pure virtual)
  // typically used for visualization, may be slow
  static void unReArrange(
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &imgR,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img,
    const WarpingParameter<SPSType, SIMD_WIDTH> &)
  {
    // this doesn't work, true copy is needed
    // img.view(imgR);
    img = imgR;
  }

  // re-arrange image (makeView instead of pure virtual)
  static void reArrange(
    const SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &stk,
    SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &stkR,
    const WarpingParameter<SPSType, SIMD_WIDTH> &)
  {
    // TODO: this is fast but may be dangerous, true copy may be advisable
    // TODO: bv: enable following when SIMDStack has view function
    // stkR.view(stk);
    stkR = stk;
  }

  // un-re-arrange image (copy instead of pure virtual)
  // typically used for visualization, may be slow
  static void unReArrange(
    const SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &stkR,
    SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &stk,
    const WarpingParameter<SPSType, SIMD_WIDTH> &)
  {
    // this doesn't work, true copy is needed
    // stk.view(stkR);
    stk = stkR;
  }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class SingleNormalizePolicyAbs
  : public NormalizePolicyNonRearranged<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
{
public:
  static void normalization(
    ConsidersInvalid<false>,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeasure,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &ssMeasure,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &cvMeasure,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &, double postScale,
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &scalePlane,
    const WarpingParameter<SPSType, SIMD_WIDTH> &par)
  {
    computeScalePlaneNormalizedAdd1<false>(jointMeasure, ssMeasure, cvMeasure,
                                           par.wSPS, par.modulo2wpw,
                                           par.shuffle, postScale, scalePlane);
  }

  static void normalization(
    ConsidersInvalid<true>,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeasure,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &combinedMeasure,
    double postScale, SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &scalePlane,
    const WarpingParameter<SPSType, SIMD_WIDTH> &par)
  {
    // here we only use the combined measure
    computeScalePlaneNormalizedAdd1Tilt<false>(
      jointMeasure, combinedMeasure, par.wSPS, par.modulo2wpw, par.shuffle,
      postScale, scalePlane);
  }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class SingleNormalizePolicyAbsRearranged
{
public:
  static const bool isRearranged = true;

  static void normalization(
    ConsidersInvalid<false>,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeasure,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &ssMeasure,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &cvMeasure,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &, double postScale,
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &scalePlane,
    const WarpingParameter<SPSType, SIMD_WIDTH> &par)
  {
    computeScalePlaneNormalizedAdd2<false>(
      jointMeasure, ssMeasure, cvMeasure, par.w, par.wSPS, par.stepAlpha,
      par.nAlpha, par.blkSize, par.modNAlpha, par.shuffleModulo, par.unshuffle,
      postScale, scalePlane);
  }

  static void reArrange(
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &imgR,
    const WarpingParameter<SPSType, SIMD_WIDTH> &par)
  {
    swizzleCyclic(img, par.stepAlpha, numInputSIMDVecs<SPSType, MeasType>(),
                  imgR);
  }

  static void unReArrange(
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &imgR,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img,
    const WarpingParameter<SPSType, SIMD_WIDTH> &par)
  {
    // serial, slow!
    unSwizzleCyclic(imgR, par.stepAlpha, numInputSIMDVecs<SPSType, MeasType>(),
                    img);
  }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class SingleMeasurePolicyNSADA
{
public:
  static void jointMeasure(
    Rearranged<false>, ConsidersInvalid<false>,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input2,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeas,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas1,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &,
    const WarpingParameter<SPSType, SIMD_WIDTH> &)
  {
    // puts("NSADA");
    normSumAbsDiffAbsTerm1(input1, input2, jointMeas);
    colMeas1.view(colMeasRaw1);
    colMeas2.view(colMeasRaw2);
  }

  static void jointMeasure(
    Rearranged<true>, ConsidersInvalid<false>,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input2,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeas,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas1,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &,
    const WarpingParameter<SPSType, SIMD_WIDTH> &par)
  {
    // puts("NSADA");
    normSumAbsDiffAbsTerm2(input1, input2, par.w, par.wSPS, par.stepAlpha,
                           par.nAlpha, par.blkSize, par.modNAlpha, jointMeas);
    colMeas1.view(colMeasRaw1);
    colMeas2.view(colMeasRaw2);
  }

  // NSAD measure is in [0,1]
  static double measureMax() { return 1.0; }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class SingleMeasurePolicyNSADNSADA
{
public:
  static void jointMeasure(
    Rearranged<false>, ConsidersInvalid<false>,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input2,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeas,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas1,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &,
    const WarpingParameter<SPSType, SIMD_WIDTH> &)
  {
    // puts("NSAD-NSADA");
    normSumAvgAbsDiffAbsDiffAbsTerm1(input1, input2, jointMeas);
    colMeas1.view(colMeasRaw1);
    colMeas2.view(colMeasRaw2);
  }

  static void jointMeasure(
    Rearranged<true>, ConsidersInvalid<false>,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input2,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeas,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas1,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &,
    const WarpingParameter<SPSType, SIMD_WIDTH> &par)
  {
    // puts("NSAD-NSADA");
    normSumAvgAbsDiffAbsDiffAbsTerm2(input1, input2, par.w, par.wSPS,
                                     par.stepAlpha, par.nAlpha, par.blkSize,
                                     par.modNAlpha, jointMeas);
    colMeas1.view(colMeasRaw1);
    colMeas2.view(colMeasRaw2);
  }

  // NSAD measure is in [0,1]
  static double measureMax() { return 1.0; }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class SingleMeasurePolicyASCDirect
{
public:
  static void jointMeasure(
    Rearranged<true>, ConsidersInvalid<false>,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input2,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeas,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas1,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &,
    const WarpingParameter<SPSType, SIMD_WIDTH> &par)
  {
    ascDirectTerm(input1, input2, par.w, par.wSPS, par.stepAlpha, par.nAlpha,
                  par.blkSize, par.modNAlpha, jointMeas);
    // colMeas1 = colMeasRaw1;
    // colMeas2 = colMeasRaw2;
    colMeas1.view(colMeasRaw1);
    colMeas2.view(colMeasRaw2);
  }

  // ASCDirect in [0,2]
  static double measureMax() { return 2.0; }
};

// ============================================================================
// Stack
// ============================================================================
// Stack is the base class for multiple preprocessing outputs. It is based on
// former WarpingSPSComputationStack and WarpingSPSComputationSNSAD. It
// implements all virtual functions of WarpingSPSComputation using static
// functions of policy classes given by template parameters:
//
// TODO: maybe split normalization and rearrange
// TODO: Preproc or magnify first?

template <SPS_COMPUTATION_TEMPLATE_PARAMETER, int NPREPROCS,
          class ImageOpsPolicy, class StackPreprocessingPolicy,
          class StackNormalizePolicy, class MeasurePolicy,
          class StackColumnMeasurePolicy>
class Stack : public WarpingSPSComputation<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
{
protected:
  // computeSPS should store input images
  bool storeSPSInput;
  // input images in sps computation
  std::vector<SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama>> ssInputs,
    cvInputs;

  // ===================== non-virtual methods ===========================

  // compute SPS and intermediate input images
  // 27. Feb 18 (rm): vert.res. and horizon taken from images
  void computeSPSInput(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
    const std::vector<double> &scaleFactors, int interpolation,
    const std::vector<double> &pixelScale, double postScale,
    WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps)
  {
    storeSPSInput = true;
    ssInputs.resize(scaleFactors.size());
    cvInputs.resize(scaleFactors.size());
    this->computeSPS(ss, cv, scaleFactors, interpolation, pixelScale, postScale,
                     sps);
    storeSPSInput = false;
  }

public:
  // ======================== ctor, dtor =========================

  Stack(typename StackPreprocessingPolicy::Parameter preproc_param,
        typename StackNormalizePolicy::Parameter normalize_param)
    : WarpingSPSComputation<ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH,
                            SIMD_ALIGN>(),
      storeSPSInput(false),
      stackPreprocessingPolicy(StackPreprocessingPolicy(preproc_param)),
      stackNormalizePolicy(StackNormalizePolicy(normalize_param))
  {}

  virtual ~Stack() {}

  // computation of scale-plane stack
  // (not possible to restrict ImgType to this function since
  // templates may not be virtual)
  void computeSPS(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
    const std::vector<double> &scaleFactors, int interpolation,
    const std::vector<double> &pixelScale, double postScale,
    WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps)
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
      throw SIMDException("WarpingSPSComputationSNSAD::computeSPS", errStr);
    }
    if (ss.addOn != cv.addOn) {
      char errStr[256];
      sprintf(errStr,
              "panorama parameters of ss (%g, %g) and cv (%g, %g) differ",
              ss.addOn.verticalResolution, ss.addOn.horizon,
              cv.addOn.verticalResolution, cv.addOn.horizon);
      throw SIMDException("WarpingSPSComputationSNSAD::computeSPS", errStr);
    }
    // check validity for given image width
    if (!this->validForWidth(w))
      throw SIMDException("WarpingSPSComputationSNSAD::computeSPS",
                          "computation not possible for given image width");
    // get input images
    // convert from input type to preprocessing type
    SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> ssOrigF, cvOrigF;
    ImageOpsPolicy::convertImage(ss, ssOrigF);
    ImageOpsPolicy::convertImage(cv, cvOrigF);
    // preprocess input images
    SIMDStack<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> ssF, cvF;
    stackPreprocessingPolicy.preprocessing(considersInvalid, ssOrigF, ssF);
    stackPreprocessingPolicy.preprocessing(considersInvalid, cvOrigF, cvF);
    // we first process the unmagnified images since the results are needed
    // multiple times
    // 2. Mar 18 (rm): *SR images are the re-arranged versions,
    // earlier in-place re-arrange required heap operations (width differs)
    SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> ssS, cvS, magS, ssSR,
      cvSR, magSR, *cvSPtr, *ssSPtr;
    SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN> ssMeasureS, cvMeasureS,
      cvMeasFinal, ssMeasFinal, combMeas, magMeasureS, jointMeas,
      *cvMeasureSPtr, *ssMeasureSPtr;
    ImageOpsPolicy::copyStack(cvF, pixelScale, cvS);
    StackNormalizePolicy::reArrange(cvS, cvSR, sps.param);
    StackColumnMeasurePolicy::columnMeasure(considersInvalid, cvSR, cvMeasureS);
    ImageOpsPolicy::copyStack(ssF, pixelScale, ssS);
    StackNormalizePolicy::reArrange(ssS, ssSR, sps.param);
    StackColumnMeasurePolicy::columnMeasure(considersInvalid, ssSR, ssMeasureS);

    // create all scale plane images
    for (unsigned int spIndex = 0; spIndex < scaleFactors.size(); spIndex++) {
      // edge images are scaled
      double scale = scaleFactors[spIndex];

      // printf("spIndex %u, scale %g\n", spIndex, scale);
      if (scale > 1.0) {
        // ===== magnify cv by 1/scale and use unmodified ss =====
        SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> magF;
        SIMDStack<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> magFPre;
        magnifyAndScale(cvOrigF, interpolation, 1.0 / scale, 1.0, magF);
        stackPreprocessingPolicy.preprocessing(considersInvalid, magF, magFPre);
        ImageOpsPolicy::copyStack(magFPre, pixelScale, magS);
        StackNormalizePolicy::reArrange(magS, magSR, sps.param);
        cvSPtr = &magSR;
        StackColumnMeasurePolicy::columnMeasure(considersInvalid, magSR,
                                                magMeasureS);
        cvMeasureSPtr = &magMeasureS;
        ssSPtr        = &ssSR;
        ssMeasureSPtr = &ssMeasureS;

      } else if (scale < 1.0) {
        // ===== magnify ss by scale and use unmodified cv =====
        cvSPtr        = &cvSR;
        cvMeasureSPtr = &cvMeasureS;
        SIMDImage<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> magF;
        SIMDStack<ProcType, SIMD_WIDTH, SIMD_ALIGN, Panorama> magFPre;
        magnifyAndScale(ssOrigF, interpolation, scale, 1, magF);
        stackPreprocessingPolicy.preprocessing(considersInvalid, magF, magFPre);
        ImageOpsPolicy::copyStack(magFPre, pixelScale, magS);
        StackNormalizePolicy::reArrange(magS, magSR, sps.param);
        ssSPtr = &magSR;
        StackColumnMeasurePolicy::columnMeasure(considersInvalid, magSR,
                                                magMeasureS);
        ssMeasureSPtr = &magMeasureS;

      } else {
        // ===== scale = 1.0: use unmodified ss and cv =====
        cvSPtr        = &cvSR;
        cvMeasureSPtr = &cvMeasureS;
        ssSPtr        = &ssSR;
        ssMeasureSPtr = &ssMeasureS;
      }
      // store input images if requested
      if (this->storeSPSInput) {
        // 27. Jul 18 (rm)
        StackNormalizePolicy::unReArrange(*ssSPtr, this->ssInputs[spIndex],
                                          sps.param);
        StackNormalizePolicy::unReArrange(*cvSPtr, this->cvInputs[spIndex],
                                          sps.param);
        // this->ssInputs[spIndex] = *ssSPtr;
        // this->cvInputs[spIndex] = *cvSPtr;
      }
      // compute distance term
      MeasurePolicy::jointMeasure(rearranged, considersInvalid, *cvSPtr,
                                  *ssSPtr, jointMeas, sps.param);
      // compute normalized measure
      stackNormalizePolicy.normalization(
        considersInvalid, jointMeas, *ssMeasureSPtr, *cvMeasureSPtr, combMeas,
        postScale, sps.stack[spIndex], sps.param);
    }
    sps.complement();
  }

  // compute a diagram showing all intermediate input images for SPS
  // 27. Feb 18 (rm): vert.res. and hor. taken from images
  void computeSPSInputDiagram(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &ss,
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &cv,
    const std::vector<double> &scaleFactors, int interpolation,
    const std::vector<double> &pixelScale, double postScale,
    WarpingSPS<SPSType, SIMD_WIDTH, SIMD_ALIGN> &sps,
    // outputs:
    MeasType &allMax, MeasType &allMin, MeasType &allRange,
    SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN> &spsInput)
  {
    // TODO: this always shows only the following preprocessing image
    computeSPSInput(ss, cv, scaleFactors, interpolation, pixelScale, postScale,
                    sps);
    unsigned nPlanes = scaleFactors.size();
    // first: determine max and min over all images
    allMax = SIMDTypeInfo<MeasType>::min();
    allMin = SIMDTypeInfo<MeasType>::max();
    for (unsigned planeIdx = 0; planeIdx < nPlanes; planeIdx++) {
      for (int d = 0; d < this->numPreprocs(); d++) {
        MeasType _min, _max;
        _max   = max(ssInputs[planeIdx][d]);
        _min   = min(ssInputs[planeIdx][d]);
        allMax = std::max(allMax, _max);
        allMin = std::min(allMin, _min);
        _max   = max(cvInputs[planeIdx][d]);
        _min   = min(cvInputs[planeIdx][d]);
        allMax = std::max(allMax, _max);
        allMin = std::max(allMin, _min);
      }
    }
    allRange = std::max(allMax, MeasType(-allMin));
    // write into same image
    SIMDImage<SIMDByte, SIMD_WIDTH, SIMD_ALIGN> tmp;
    // asserting all preprocs have same width and height
    int w = ssInputs[0][0].w, h = ssInputs[0][0].h;
    const int margin = SIMDVec<SIMDByte, SIMD_WIDTH>::elements;
    spsInput.resize(this->numPreprocs() * 2 * (w + margin) + margin,
                    nPlanes * (h + margin) + margin);
    spsInput.fill(255);
    double shift = 128.0, scale = 255.0 / allRange;
    int verticalOffset = 0;
    // TODO: replace numPreprocs with NPREPROCS
    for (int d = 0; d < this->numPreprocs(); d++) {
      for (unsigned planeIdx = 0; planeIdx < nPlanes; planeIdx++) {
        copyScaleAndShift(ssInputs[planeIdx][d], scale, shift, tmp);
        insert(tmp, verticalOffset + margin, margin + planeIdx * (h + margin),
               spsInput);
        copyScaleAndShift(cvInputs[planeIdx][d], scale, shift, tmp);
        insert(tmp, verticalOffset + w + 2 * margin,
               margin + planeIdx * (h + margin), spsInput);
      }
      verticalOffset += 2 * w + 2 * margin;
    }
  }

  int numPreprocs() const { return NPREPROCS; }

  std::vector<double> maxDenom(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img)
  {
    return stackPreprocessingPolicy.maxDenom(img);
  };

  // if computed with real numbers, maximal value of the normalized measure
  // (used for WarpingSPSComputation::idealPostScale())
  virtual double measureMax() const { return MeasurePolicy::measureMax(); };

private:
  // check whether input images with given width could be processed
  // this definition should be useful for many cases:
  virtual bool validForWidth(int width) const
  {
    return
      // ProcType->MeasType: uses magnifyImage, copyImage:
      ((width % numSIMDVecsElements<MeasType, ProcType, SIMD_WIDTH>() == 0) &&
       // MeasType->SPSType: uses computeScalePlaneNormalizedAdd:
       (width % numSIMDVecsElements<SPSType, MeasType, SIMD_WIDTH>() == 0));
  }

  static const Rearranged<StackNormalizePolicy::isRearranged> rearranged;
  static const ConsidersInvalid<ImageOpsPolicy::considersInvalidPixels>
    considersInvalid;

  StackPreprocessingPolicy stackPreprocessingPolicy;
  StackNormalizePolicy stackNormalizePolicy;
};

// To prevent clang++ warning for c++98 (no instantiation of variable). Can't be
// in class, because it is static (this is relaxed in c++17). This is meant to
// (most likely) be removed with compiler optimization (empty struct and unused
// variables)
template <SPS_COMPUTATION_TEMPLATE_PARAMETER, int NPREPROCS,
          class ImageOpsPolicy, class StackPreprocessingPolicy,
          class StackNormalizePolicy, class MeasurePolicy,
          class StackColumnMeasurePolicy>
const Rearranged<StackNormalizePolicy::isRearranged>
  Stack<SPS_COMPUTATION_TEMPLATE_ARGUMENTS, NPREPROCS, ImageOpsPolicy,
        StackPreprocessingPolicy, StackNormalizePolicy, MeasurePolicy,
        StackColumnMeasurePolicy>::rearranged = {};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER, int NPREPROCS,
          class ImageOpsPolicy, class StackPreprocessingPolicy,
          class StackNormalizePolicy, class MeasurePolicy,
          class StackColumnMeasurePolicy>
const ConsidersInvalid<ImageOpsPolicy::considersInvalidPixels>
  Stack<SPS_COMPUTATION_TEMPLATE_ARGUMENTS, NPREPROCS, ImageOpsPolicy,
        StackPreprocessingPolicy, StackNormalizePolicy, MeasurePolicy,
        StackColumnMeasurePolicy>::considersInvalid = {};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER, int NPREPROCS,
          HarrisOutputSelector HOS>
class StackPreprocessingPolicyHarris
{
public:
  class Parameter
  {
  public:
    Parameter(double k, int binomialFilterCount)
      : k(k), binomialFilterCount(binomialFilterCount)
    {}
    double k;
    int binomialFilterCount;
  };

  // ctor
  StackPreprocessingPolicyHarris(
    StackPreprocessingPolicyHarris::Parameter param)
    : k(param.k), binomialFilterCount(param.binomialFilterCount) {};

  // preprocessing applied to the input images
  void preprocessing(
    ConsidersInvalid<false>,
    const SIMDImage<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input,
    SIMDStack<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN, Panorama> &output)
  {
    harrisFilterCyclic<2, 0, HOS>(input, output, k, binomialFilterCount);
  }

  std::vector<double> maxDenom(
    const SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &img)
  {
    std::vector<double> md(NPREPROCS);
    // all processing is done in SIMDFloat to avoid overflow
    SIMDImage<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN, Panorama> inputImg;
    SIMDStack<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN, Panorama> preprocStk;
    SIMDImage<SIMDFloat, SIMD_WIDTH, SIMD_ALIGN> sumOfAbsImg;
    inputImg = img;
    harrisFilterCyclic<2, 0, HOS>(inputImg, preprocStk, k, binomialFilterCount);
    for (unsigned int pp = 0; pp < NPREPROCS; pp++) {
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

template <SPS_COMPUTATION_TEMPLATE_PARAMETER, int NPREPROCS>
class StackColumnMeasurePolicyAbs
{
public:
  template <typename T>
  static void columnMeasure(
    ConsidersInvalid<false>,
    const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, Panorama> &in,
    SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN> &jointMeas)
  {
    // jointMeas.h == 1 i.e. one value for each column
    jointMeas.resize(in.numPlanes, in.w, 1);
    // TODO: replace loop with own functions
    for (int i = 0; i < in.numPlanes; i++) { sumOfAbs(in[i], jointMeas[i]); }
  }
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER, int NPREPROCS>
class StackNormalizePolicySigmoid
  : public NormalizePolicyNonRearranged<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
{
public:
  class Parameter
  {
  public:
    Parameter(std::vector<double> w0, std::vector<double> w) : w0(w0), w(w) {}
    std::vector<double> w0;
    std::vector<double> w;
  };

  StackNormalizePolicySigmoid(Parameter param)
    : NormalizePolicyNonRearranged<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>(),
      w0(param.w0), w(param.w)
  {}

  // normalization of jointMeasure and ssMeasure/cvMeasure
  void normalization(
    ConsidersInvalid<false>,
    const SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeasure,
    const SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN> &ssMeasure,
    const SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN> &cvMeasure,
    const SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN> &, double postScale,
    SIMDImage<SPSType, SIMD_WIDTH, SIMD_ALIGN> &scalePlane,
    const WarpingParameter<SPSType, SIMD_WIDTH> &par)
  {
    computeScalePlaneSigmoidNormalizedAdd1<false, NPREPROCS>(
      jointMeasure, ssMeasure, cvMeasure, par.wSPS, par.modulo2wpw, par.shuffle,
      this->w, this->w0, postScale, scalePlane);
  }

private:
  const std::vector<double> w0;
  const std::vector<double> w;
};

template <SPS_COMPUTATION_TEMPLATE_PARAMETER>
class MeasurePolicyNSAD
{
public:
  static void jointMeasure(
    Rearranged<false>, ConsidersInvalid<false>,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input2,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeas,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas1,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &,
    const WarpingParameter<SPSType, SIMD_WIDTH> &)
  {
    // puts("NSAD");
    normSumAbsDiffTerm1(input1, input2, jointMeas);
    colMeas1.view(colMeasRaw1);
    colMeas2.view(colMeasRaw2);
  }

  static void jointMeasure(
    Rearranged<true>, ConsidersInvalid<false>,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input2,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeas,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas1,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &,
    const WarpingParameter<SPSType, SIMD_WIDTH> &par)
  {
    // puts("NSAD");
    normSumAbsDiffTerm2(input1, input2, par.w, par.wSPS, par.stepAlpha,
                        par.nAlpha, par.blkSize, par.modNAlpha, jointMeas);
    colMeas1.view(colMeasRaw1);
    colMeas2.view(colMeasRaw2);
  }

  static void jointMeasure(
    Rearranged<false>, ConsidersInvalid<true>,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &input2,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw1,
    const SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeasRaw2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &output,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas1,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas2,
    SIMDImage<MeasType, SIMD_WIDTH, SIMD_ALIGN> &colMeas12,
    const WarpingParameter<SPSType, SIMD_WIDTH> &)
  {
    // with NSAD match operator
    columnPairMatchingTilt1<MeasType, SIMD_WIDTH, SIMD_ALIGN,
                            MatchAbsDiff<MeasType, SIMD_WIDTH>>(
      input1, input2, output, colMeasRaw1, colMeasRaw2, colMeas1, colMeas2,
      colMeas12);
  }

  // Stack variant
  static void jointMeasure(
    Rearranged<false>, ConsidersInvalid<false>,
    const SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &inputStk1,
    const SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN, Panorama> &inputStk2,
    SIMDStack<MeasType, SIMD_WIDTH, SIMD_ALIGN> &jointMeasStk,
    const WarpingParameter<SPSType, SIMD_WIDTH> &)
  {
    // jointMeasStk.h == input.w
    jointMeasStk.resize(inputStk1.numPlanes, inputStk1.w, inputStk1.w);
    // TODO: remove next line?
    normSumAbsDiffTerm1(inputStk1[0], inputStk2[0], jointMeasStk[0]);
    // TODO: NPREPROCS template instead of .numPlanes
    for (int i = 0; i < inputStk1.numPlanes; i++) {
      // puts("NSAD");
      normSumAbsDiffTerm1(inputStk1[i], inputStk2[i], jointMeasStk[i]);
    }
  }

  // NSAD measure is in [0,1]
  static double measureMax() { return 1.0; }
};

} // namespace modular_sps_computation

} // namespace simd

#endif
