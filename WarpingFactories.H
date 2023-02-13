// ===========================================================================
//
// WarpingFactories.H --
// factory classes for different warping sub-classes
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

// https://
// en.wikibooks.org/wiki/C%2B%2B_Programming/Code/Design_Patterns#Factory

#ifndef _WARPING_FACTORIES_H_
#define _WARPING_FACTORIES_H_

// 18. Apr 18 (rm): removed, not clear why this was needed
// #include "WarpingCompound.H"
// 18. Apr 18 (rm): added, thanks to Benedikt Volkmer
#include "MinWarpingSearchCollection.H"
#include "WarpingSPSComputationCollection.H"
#include "WarpingSPSComputationCollectionModular.H"
#include "WarpingSearchCollection.H"

namespace simd {

// -------------------------------------------------------------------------
// factory for WarpingSPSComputation
// -------------------------------------------------------------------------

// simplifies access to enum constants (no need to pass template
// arguments), class keeps it out of the global scope
class WarpingSPSComputationSelector
{
public:
  enum FirstPhaseType {
    firstPhaseASC1,       // 0  approximated sequential correlation
    firstPhaseNSAD1,      // 1  normalized sum of absolute differences
    firstPhaseASC2,       // 2  swizzled version (faster on AVX)
    firstPhaseNSAD2,      // 3  swizzled version (faster on AVX)
    firstPhaseNSADA1,     // 4  NSADA
    firstPhaseNSADNSADA1, // 5  NSAD-NSADA mix
    firstPhaseNSADA2,     // 6  swizzled version (faster on AVX)
    firstPhaseNSADNSADA2, // 7  swizzled version (faster on AVX)
    firstPhaseASCDirect1, // 8  unused, reserved for ASCDirect1
    firstPhaseASCDirect2, // 9  experimental: ASCDirect2
    firstPhaseNSADTilt,   // 10 NSAD considering invalid pixels
    // Alternative version of NSAD with magnify/preprocess order switched
    // contributed by Benedikt Volkmer
    firstPhaseNSAD2EM, // 11 magnify first
    firstPhaseNSAD2EP, // 12 preprocess first
    // Corner filtered versions of NSADA2
    // contributed by Benedikt Volkmer
    firstPhaseNSAD2CP,   // 13 no normalization, preprocess first
    firstPhaseNSAD2CPN,  // 14 normalized (sqrt-sqrt), preprocess first
    firstPhaseNSAD2CM,   // 15 no normalization, magnify first
    firstPhaseNSAD2CMN,  // 16 normalized (sqrt-sqrt), magnify first
    firstPhaseNSAD2CMNV, // 17 normalized (sqrt-sqrt), magnify first, (1-k)
    firstPhaseNSAD2CMNS, // 18 normalized (sqrt-sqrt), mag. first, scale mag
    firstPhaseNSAD2CPn,  // 19 normalized (sqrt), preprocess first
    firstPhaseNSAD2CMn,  // 20 normalized (sqrt), magnify first
    firstPhaseSNSADEHL1, // 21 sigmoid on NSAD of edge, harris response and
                         //    harris l1
    firstPhaseSNSADE,    // 22 sigmoid on NSAD of edge
    // version of firstPhaseNSADTilt with edge filter (C. Berganski)
    firstPhaseEdgeNSADTilt, // 23 NSAD considering invalid pixels
    // Modular versions (Benedikt Volkmer)
    firstPhaseModularNSAD,        // 24 (-> 1)
    firstPhaseModularNSAD2,       // 25 (-> 3)
    firstPhaseModularNSADA,       // 26 (-> 4)
    firstPhaseModularNSADNSADA,   // 27 (-> 5)
    firstPhaseModularNSADA2,      // 28 (-> 6)
    firstPhaseModularNSADNSADA2,  // 29 (-> 7)
    firstPhaseModularASCDirect2,  // 30 (-> 9)
    firstPhaseModularNSADTilt,    // 31 (-> 10)
    firstPhaseModularNSADCM,      // 32 (-> 16)
    firstPhaseModularSNSAD_EHL1,  // 33 (-> 21)
    firstPhaseModularEdgeNSADTilt // 34 (-> 23)
  };
};

template <typename ImgType, typename ProcType, typename MeasType,
          typename SPSType, int SIMD_WIDTH, int SIMD_ALIGN>
class WarpingSPSComputationFactory : public WarpingSPSComputationSelector
{
public:
  static WarpingSPSComputation<ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH,
                               SIMD_ALIGN>
    *create(int firstPhase, const HarrisParameter &harrisPar,
            const std::vector<double> &sigmoidW,
            const std::vector<double> &sigmoidW0)
  {
    switch (static_cast<FirstPhaseType>(firstPhase)) {
    case firstPhaseASC1:
      return new WarpingSPSComputationEdgeASC1<
        ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN>();
      break;

    case firstPhaseNSAD1:
      return new WarpingSPSComputationEdgeNSAD1<
        ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN>();
      break;

    case firstPhaseASC2:
      return new WarpingSPSComputationEdgeASC2<
        ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN>();
      break;

    case firstPhaseNSAD2:
      return new WarpingSPSComputationEdgeNSAD2<
        ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN>();
      break;

    case firstPhaseNSADA1:
      return new WarpingSPSComputationEdgeNSADA1<
        ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN>();
      break;

    case firstPhaseNSADA2:
      return new WarpingSPSComputationEdgeNSADA2<
        ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN>();
      break;

    case firstPhaseNSADNSADA1:
      return new WarpingSPSComputationEdgeNSADNSADA1<
        ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN>();
      break;

    case firstPhaseNSADNSADA2:
      return new WarpingSPSComputationEdgeNSADNSADA2<
        ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN>();
      break;

    case firstPhaseASCDirect2:
      return new WarpingSPSComputationEdgeASCDirect2<
        ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN>();
      break;

    case firstPhaseNSADTilt:
      return new WarpingSPSComputationNSADTilt<
        ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN>();
      break;

    case firstPhaseNSAD2EM:
      return new WarpingSPSComputationEdgeNSAD2Alt<ImgType, ProcType, MeasType,
                                                   SPSType, SIMD_WIDTH,
                                                   SIMD_ALIGN, true, false>();
      break;

    case firstPhaseNSAD2EP:
      return new WarpingSPSComputationEdgeNSAD2Alt<ImgType, ProcType, MeasType,
                                                   SPSType, SIMD_WIDTH,
                                                   SIMD_ALIGN, false, false>();
      break;

    case firstPhaseNSAD2CP:
      return new WarpingSPSComputationCornerNSAD2<ImgType, MeasType, SPSType,
                                                  SIMD_WIDTH, SIMD_ALIGN, false,
                                                  false, 0, false>(harrisPar);
      break;

    case firstPhaseNSAD2CPn:
      return new WarpingSPSComputationCornerNSAD2<ImgType, MeasType, SPSType,
                                                  SIMD_WIDTH, SIMD_ALIGN, false,
                                                  false, 1, false>(harrisPar);
      break;

    case firstPhaseNSAD2CPN:
      return new WarpingSPSComputationCornerNSAD2<ImgType, MeasType, SPSType,
                                                  SIMD_WIDTH, SIMD_ALIGN, false,
                                                  false, 2, false>(harrisPar);
      break;

    case firstPhaseNSAD2CM:
      return new WarpingSPSComputationCornerNSAD2<ImgType, MeasType, SPSType,
                                                  SIMD_WIDTH, SIMD_ALIGN, true,
                                                  false, 0, false>(harrisPar);
      break;

    case firstPhaseNSAD2CMn:
      return new WarpingSPSComputationCornerNSAD2<ImgType, MeasType, SPSType,
                                                  SIMD_WIDTH, SIMD_ALIGN, true,
                                                  false, 1, false>(harrisPar);
      break;

    case firstPhaseNSAD2CMN:
      return new WarpingSPSComputationCornerNSAD2<ImgType, MeasType, SPSType,
                                                  SIMD_WIDTH, SIMD_ALIGN, true,
                                                  false, 2, false>(harrisPar);
      break;

    case firstPhaseNSAD2CMNV:
      return new WarpingSPSComputationCornerNSAD2<ImgType, MeasType, SPSType,
                                                  SIMD_WIDTH, SIMD_ALIGN, true,
                                                  false, 2, true>(harrisPar);
      break;

    case firstPhaseNSAD2CMNS:
      return new WarpingSPSComputationCornerNSAD2<ImgType, MeasType, SPSType,
                                                  SIMD_WIDTH, SIMD_ALIGN, true,
                                                  true, 2, false>(harrisPar);
      break;

    case firstPhaseSNSADEHL1: {
      const int n = 3;
      if ((sigmoidW.size() != n) || (sigmoidW0.size() != n))
        throw SIMDException("WarpingSPSComputationFactory::create",
                            "sigmoidW/W0 has wrong size");
      // double sigmoidW_a[n], sigmoidW0_a[n];
      // std::copy(sigmoidW.begin(), sigmoidW.end(), sigmoidW_a);
      // std::copy(sigmoidW0.begin(), sigmoidW0.end(), sigmoidW0_a);
      return new WarpingSPSComputationSNSADHarris<
        ImgType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN,
        HOS_VEDGE | HOS_RESPONSE | HOS_L1, n>
        // (harrisPar, sigmoidW_a, sigmoidW0_a);
        (harrisPar, sigmoidW, sigmoidW0);
      break;
    }

    case firstPhaseSNSADE: {
      const int n = 1;
      if ((sigmoidW.size() != n) || (sigmoidW0.size() != n))
        throw SIMDException("WarpingSPSComputationFactory::create",
                            "sigmoidW/W0 has wrong size");
      // double sigmoidW_a[n], sigmoidW0_a[n];
      // std::copy(sigmoidW.begin(), sigmoidW.end(), sigmoidW_a);
      // std::copy(sigmoidW0.begin(), sigmoidW0.end(), sigmoidW0_a);
      return new WarpingSPSComputationSNSADHarris<
        ImgType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN, HOS_VEDGE, n>
        // (harrisPar, sigmoidW_a, sigmoidW0_a);
        (harrisPar, sigmoidW, sigmoidW0);
      break;
    }

    case firstPhaseEdgeNSADTilt:
      return new WarpingSPSComputationEdgeNSADTilt<
        ImgType, ProcType, MeasType, SPSType, SIMD_WIDTH, SIMD_ALIGN>();
      break;

    case firstPhaseModularNSAD: {
      using namespace modular_sps_computation;
      typedef SinglePreprocessingPolicyEdge<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
        PreprocPolicy;
      return new Single<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, false,
        ImageOpsPolicySimple<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>, PreprocPolicy,
        SingleNormalizePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        MeasurePolicyNSAD<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleColumnMeasurePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>>(
        typename PreprocPolicy::Parameter());
      break;
    }

    case firstPhaseModularNSAD2: {
      using namespace modular_sps_computation;
      typedef SinglePreprocessingPolicyEdge<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
        PreprocPolicy;
      return new Single<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, false,
        ImageOpsPolicySimple<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>, PreprocPolicy,
        SingleNormalizePolicyAbsRearranged<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        MeasurePolicyNSAD<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleColumnMeasurePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>>(
        typename PreprocPolicy::Parameter());
      break;
    }

    case firstPhaseModularNSADA: {
      using namespace modular_sps_computation;
      typedef SinglePreprocessingPolicyEdge<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
        PreprocPolicy;
      return new Single<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, false,
        ImageOpsPolicySimple<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>, PreprocPolicy,
        SingleNormalizePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleMeasurePolicyNSADA<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleColumnMeasurePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>>(
        typename PreprocPolicy::Parameter());
      break;
    }

    case firstPhaseModularNSADA2: {
      using namespace modular_sps_computation;
      typedef SinglePreprocessingPolicyEdge<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
        PreprocPolicy;
      return new Single<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, false,
        ImageOpsPolicySimple<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>, PreprocPolicy,
        SingleNormalizePolicyAbsRearranged<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleMeasurePolicyNSADA<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleColumnMeasurePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>>(
        typename PreprocPolicy::Parameter());
      break;
    }

    case firstPhaseModularNSADNSADA: {
      using namespace modular_sps_computation;
      typedef SinglePreprocessingPolicyEdge<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
        PreprocPolicy;
      return new Single<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, false,
        ImageOpsPolicySimple<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>, PreprocPolicy,
        SingleNormalizePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleMeasurePolicyNSADNSADA<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleColumnMeasurePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>>(
        typename PreprocPolicy::Parameter());
      break;
    }

    case firstPhaseModularNSADNSADA2: {
      using namespace modular_sps_computation;
      typedef SinglePreprocessingPolicyEdge<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
        PreprocPolicy;
      return new Single<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, false,
        ImageOpsPolicySimple<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>, PreprocPolicy,
        SingleNormalizePolicyAbsRearranged<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleMeasurePolicyNSADNSADA<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleColumnMeasurePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>>(
        typename PreprocPolicy::Parameter());
      break;
    }

    case firstPhaseModularASCDirect2: {
      using namespace modular_sps_computation;
      typedef SinglePreprocessingPolicyEdge<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
        PreprocPolicy;
      return new Single<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, false,
        ImageOpsPolicySimple<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>, PreprocPolicy,
        SingleNormalizePolicyAbsRearranged<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleMeasurePolicyASCDirect<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleColumnMeasurePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>>(
        typename PreprocPolicy::Parameter());
      break;
    }

    case firstPhaseModularNSADTilt: {
      using namespace modular_sps_computation;
      typedef SinglePreprocessingPolicyNon<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
        PreprocPolicy;
      return new Single<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, false,
        SingleImageOpsPolicyInvalid<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        PreprocPolicy,
        SingleNormalizePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        MeasurePolicyNSAD<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleColumnMeasurePolicyAbsAccu<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>>(
        typename PreprocPolicy::Parameter());
      break;
    }

    case firstPhaseModularEdgeNSADTilt: {
      using namespace modular_sps_computation;
      typedef SinglePreprocessingPolicyEdge<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>
        PreprocPolicy;
      return new Single<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, false,
        SingleImageOpsPolicyInvalid<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        PreprocPolicy,
        SingleNormalizePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        MeasurePolicyNSAD<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleColumnMeasurePolicyAbsAccu<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>>(
        typename PreprocPolicy::Parameter());
      break;
    }

    case firstPhaseModularNSADCM: {
      using namespace modular_sps_computation;
      typedef SinglePreprocessingPolicyHarris<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, 0, false, HOS_RESPONSE>
        PreprocPolicy;
      return new Single<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, false,
        ImageOpsPolicySimple<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>, PreprocPolicy,
        SingleNormalizePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        MeasurePolicyNSAD<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        SingleColumnMeasurePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>>(
        typename PreprocPolicy::Parameter(harrisPar.k,
                                          harrisPar.binomialFilterCount));
      break;
    }

    case firstPhaseModularSNSAD_EHL1: {
      using namespace modular_sps_computation;
      const int n = 3;
      if (sigmoidW.size() != n)
        throw SIMDException("WarpingSPSComputationFactory::create",
                            "sigmoidW has wrong size");
      typedef StackPreprocessingPolicyHarris<SPS_COMPUTATION_TEMPLATE_ARGUMENTS,
                                             n,
                                             HOS_VEDGE | HOS_RESPONSE | HOS_L1>
        PreprocPolicy;
      typedef StackNormalizePolicySigmoid<SPS_COMPUTATION_TEMPLATE_ARGUMENTS, n>
        NormalizePolicy;
      return new Stack<
        SPS_COMPUTATION_TEMPLATE_ARGUMENTS, n,
        ImageOpsPolicySimple<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>, PreprocPolicy,
        NormalizePolicy, MeasurePolicyNSAD<SPS_COMPUTATION_TEMPLATE_ARGUMENTS>,
        StackColumnMeasurePolicyAbs<SPS_COMPUTATION_TEMPLATE_ARGUMENTS, n>>(
        typename PreprocPolicy::Parameter(harrisPar.k,
                                          harrisPar.binomialFilterCount),
        typename NormalizePolicy::Parameter(sigmoidW0, sigmoidW));
      break;
    }

    default:
      throw SIMDException("WarpingSPSComputationFactory::create",
                          "invalid firstPhase argument");
    }
  }
};

// -------------------------------------------------------------------------
// factory for WarpingFull, WarpingPartial
// -------------------------------------------------------------------------

class WarpingSearchSelector
{
public:
  enum SearchMethodType { searchMethodWarpingPARX };
};

template <typename SPSType, typename MatchType, int SIMD_WIDTH, int SIMD_ALIGN>
class WarpingFullFactory : public WarpingSearchSelector
{
public:
  static WarpingFull<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN> *create(
    int searchMethod)
  {
    switch (static_cast<SearchMethodType>(searchMethod)) {
    case searchMethodWarpingPARX:
      return new WarpingFull_PARX<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN>();
      break;
    default:
      throw SIMDException("WarpingFullFactory::create",
                          "invalid searchMethod argument");
    }
  }
};

template <typename SPSType, typename MatchType, int SIMD_WIDTH, int SIMD_ALIGN>
class WarpingPartialFactory : public WarpingSearchSelector
{
public:
  static WarpingPartial<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN> *create(
    int searchMethod)
  {
    switch (static_cast<SearchMethodType>(searchMethod)) {
    case searchMethodWarpingPARX:
      return new WarpingPartial_PARX<SPSType, MatchType, SIMD_WIDTH,
                                     SIMD_ALIGN>();
      break;
    default:
      throw SIMDException("WarpingPartialFactory::create",
                          "invalid searchMethod argument");
    }
  }
};

// -------------------------------------------------------------------------
// factory for MinWarpingFull, MinWarpingPartial
// -------------------------------------------------------------------------

class MinWarpingSearchSelector
{
public:
  enum SearchMethodType { searchMethodMinWarpingXPAY };
};

template <typename SPSType, typename MatchType, int SIMD_WIDTH, int SIMD_ALIGN>
class MinWarpingFullFactory : public MinWarpingSearchSelector
{
public:
  static MinWarpingFull<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN> *create(
    int searchMethod)
  {
    switch (static_cast<SearchMethodType>(searchMethod)) {
    case searchMethodMinWarpingXPAY:
      return new MinWarpingFull_XPAY<SPSType, MatchType, SIMD_WIDTH,
                                     SIMD_ALIGN>();
      break;
    default:
      throw SIMDException("MinWarpingFullFactory::create",
                          "invalid searchMethod argument");
    }
  }
};

template <typename SPSType, typename MatchType, int SIMD_WIDTH, int SIMD_ALIGN>
class MinWarpingPartialFactory : public MinWarpingSearchSelector
{
public:
  static MinWarpingPartial<SPSType, MatchType, SIMD_WIDTH, SIMD_ALIGN> *create(
    int searchMethod)
  {
    switch (static_cast<SearchMethodType>(searchMethod)) {
    case searchMethodMinWarpingXPAY:
      return new MinWarpingPartial_XPAY<SPSType, MatchType, SIMD_WIDTH,
                                        SIMD_ALIGN>();
      break;
    default:
      throw SIMDException("MinWarpingPartialFactory::create",
                          "invalid searchMethod argument");
    }
  }
};

} // namespace simd

#endif
