// ===========================================================================
//
// WarpingTemplate.H --
// template for 2D Warping
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

#ifndef _WARPING_TEMPLATE_H_
#define _WARPING_TEMPLATE_H_

#include "MinWarpingTemplate.H" // uses initializeScaleFactorsAndThresholds
#include "SIMDDefs.H"
#include "SIMDException.H"
#include "SIMDImage.H"

#include <cmath>
#include <limits>
#include <stdint.h>
#include <stdio.h>
#include <vector>

namespace ns_simd {

// =========================================================================
// data types for warping template
// =========================================================================

typedef uint32_t TemplateSPOffset;
typedef int32_t TemplateJY;

typedef SIMD_ATTR_PACKED_STRUCT
{
  TemplateSPOffset sp;
  TemplateJY jy;
}
WarpingTemplateEntry;

// =========================================================================
// Warping template
// =========================================================================

class WarpingTemplate
{
public:
  // width and sps width
  const int w, wSPS;
  // minimal rho, maximal rho
  const double rhoMin, rhoMax;
  // number of rho values between and including rhoMin and rhoMax
  const int nRho;
  // angular with
  const double dTheta;

  // vector storing all thresholds
  std::vector<double> thresholds;

  // rho values
  std::vector<double> rhoTable;
  // template is a vector of vectors (may not be the fastest solution)
  std::vector<std::vector<WarpingTemplateEntry>> tmplt;

  WarpingTemplate(int width, int spsWidth,
                  const std::vector<double> &scaleFactors, double rhoMin,
                  double rhoMax, int nRho)
    : w(width), wSPS(spsWidth), rhoMin(rhoMin), rhoMax(rhoMax), nRho(nRho),
      dTheta((2.0 * M_PI) / w)
  {
    // 10. Nov 11 (rm): we use TemplateSPOffset as offset (sp)
    // in the scale plane stack, this limits the size of the sps (we
    // compute it this way in order not to exceed any limits in this
    // computation, more instructive form of condition is:
    // w*w*nScalePlanes > SCALEPLANEOFFSETTYPE_MAX)
    if (std::numeric_limits<TemplateSPOffset>::max() / scaleFactors.size() <
        (unsigned int) (w * w))
      throw SIMDException("WarpingTemplate::WarpingTemplate",
                          "size of scale-plane stack exceeds offset limit");
    initializeThresholds(scaleFactors);
    createRhoTable();
    createTemplate();
  }

  virtual ~WarpingTemplate() {}

  static double warpFunc(double x, double rho)
  {
    // this special case is actually undefined (could be 0 or PI),
    // we take the intermediate value in order to avoid empty
    // entries in the template (and it is a continuation of the
    // function for rho = 1)
    if ((rho == 1.0) && (x == 0.0)) return M_PI / 2.0;
    // general case: equation (1) in Möller et al. "Three 2D-Warping
    // Schemes for Visual Robot Navigation"
    return atan2(rho * sin(x), 1.0 - rho * cos(x));
  }

  static double sigmaFunc(double x, double rho)
  {
    // equation (4) in Möller et al. "Three 2D-Warping
    // Schemes for Visual Robot Navigation"
    return std::sqrt(1.0 + rho * rho - 2.0 * rho * cos(x));
  }

  // save in human-readable form
  void save(FILE *f) const
  {
    for (int jx = 0; jx < w; jx++) {
      for (int iRho = 0; iRho < nRho; iRho++) {
        fprintf(f, "%4d %2d\t", int(tmplt[iRho][jx].jy),
                int(tmplt[iRho][jx].sp / (wSPS * w)));
      }
      fprintf(f, "\n");
    }
  }

protected:
  // thresholds are chosen in the center between neighboring scale factors
  void initializeThresholds(const std::vector<double> &scaleFactors)
  {
    int nScalePlanes = scaleFactors.size();
    thresholds.resize(nScalePlanes);
    // initialize thresholds: center between scaleFactors
    for (int i = 0; i < nScalePlanes - 1; i++)
      thresholds[i] = 0.5 * (scaleFactors[i] + scaleFactors[i + 1]);
    // upper threshold is fixed
    thresholds[nScalePlanes - 1] = std::numeric_limits<double>::max();
  }

  void createRhoTable()
  {
    if (nRho < 1)
      throw SIMDException("WarpingTemplate::createRhoTable",
                          "nRho >= 1 required");
    rhoTable.resize(nRho);
    if (nRho == 1)
      rhoTable[0] = rhoMin; // we take rhoMin
    else {
      double dRho = (rhoMax - rhoMin) / (nRho - 1), rho = rhoMin;
      for (int i = 0; i < nRho; i++, rho += dRho) rhoTable[i] = rho;
    }
  }

  void createTemplate()
  {
    tmplt.resize(nRho);
    for (int iRho = 0; iRho < nRho; iRho++) {
      tmplt[iRho].resize(w);
      double rho = rhoTable[iRho];
      for (int jx = 0; jx < w; jx++) {
        double x           = jx * dTheta;
        double y           = warpFunc(x, rho);
        int jy             = (w + TemplateJY(::round(y / dTheta))) % w;
        tmplt[iRho][jx].jy = jy;
        double sigma       = sigmaFunc(x, rho);
        if (sigma <= thresholds[thresholds.size() - 1]) {
          bool found = false;
          unsigned t;
          for (t = 0; t < thresholds.size(); t++)
            if (sigma <= thresholds[t]) {
              found = true;
              break;
            }
          if (!found)
            // shouldn't happen (see initializeThresholds())
            throw SIMDException("WarpingTemplate::createTemplate",
                                "scale plane not found");
          tmplt[iRho][jx].sp = (TemplateSPOffset) (t * w * wSPS);
        } else
          // shouldn't happen (see initializeThresholds())
          throw SIMDException("WarpingTemplate::createTemplate",
                              "sigma >= maxThreshold");
      }
    }
  }
};

} // namespace ns_simd

#endif
