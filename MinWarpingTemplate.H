// ===========================================================================
//
// MinWarpingTemplate.H --
// class encoding the minWarping template (x-y arguments)
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

#pragma once
#ifndef MINWARPING_TEMPLATE_H_
#define MINWARPING_TEMPLATE_H_

#include <limits>
#include <stdint.h>
#include <stdio.h>
#include <vector>
// 20. Apr 18 (rm): thanks to Benedikt Volkmer
#include "SIMDDefs.H"
#include "SIMDException.H"
#include "SIMDImage.H"

#include <cmath>

// 20. Sep 22 (Jonas Keller): on Windows M_PI is not defined, so we
// define it here
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace simd {

// =========================================================================
// data types for min-template
// =========================================================================

typedef uint32_t MinTemplateSPOffset;
typedef int32_t MinTemplateIY;

// min_template_col_entry_t contains exactly one entry
// of a min template column.
// sp is an address offset, which has to be added to the start address of the
// scale plane stack and then points to a specific scaleplane.
typedef SIMD_ATTR_PACKED_STRUCT
{
  MinTemplateSPOffset sp;
  // 19. Jun 15 (rm): uint32_t is faster (25600->23500) in this
  // version, but not in CODE3 (there slightly slower)
  // uint16_t iy;
  MinTemplateIY iy;
}
min_templt_col_entry_t;

typedef min_templt_col_entry_t *lp_min_templt_col_entry_t;

// Structure which represents a single column.
// 20. Sep 22 (Jonas Keller):
// used struct name {...}; instead of typedef struct {...} name;
// because compiler warns that type is not C-compatible for linkage
SIMD_ATTR_PACKED_STRUCT min_templt_col_t
{
  // 19. Jun 15 (rm)
  // uint32_t was slower than uint16_t, but with new pointer-based
  // loop, int was better (22500 -> 22000), int32_t even better?
  // uint16_t len;
  // uint32_t len;
  int32_t len;
  // 27. Nov 15 (rm): alternative adapted from suggestion by (as)
  // min_templt_col_entry_t entries[];
  SIMD_INLINE min_templt_col_entry_t *entries()
  {
    return (min_templt_col_entry_t *) ((uint8_t *) (this) + sizeof(len));
  }
};

// Pointer type for the structure above
typedef min_templt_col_t *lp_min_templt_col_t;

// =========================================================================
// MinWarping template
// =========================================================================

class MinWarpingTemplate
{
public:
  // width and sps width (latter needed only to compute plane offset
  // for getImage() and save())
  const int w, wSPS;
  // minimal / maximal rho
  const double rhoMin, rhoMax;
  // angular width
  const double dTheta;

  // vector storing all thresholds
  std::vector<double> thresholds;

  // Pointer on the first column
  lp_min_templt_col_t minTemplate;
  // ptr to the columns of the template mentioned above
  lp_min_templt_col_t *minTemplateColumns;

  MinWarpingTemplate(int width, int spsWidth,
                     const std::vector<double> &scaleFactors,
                     double maxThreshold, double rhoMin, double rhoMax)
    : w(width), wSPS(spsWidth), rhoMin(rhoMin), rhoMax(rhoMax),
      dTheta((2.0 * M_PI) / w), minTemplate(0), minTemplateColumns(0)
  {
    // 10. Nov 11 (rm): we use MinTemplateSPOffset as offset (sp)
    // in the scale plane stack, this limits the size of the sps (we
    // compute it this way in order not to exceed any limits in this
    // computation, more instructive form of condition is:
    // w*w*nScalePlanes > SCALEPLANEOFFSETTYPE_MAX)
    if (std::numeric_limits<MinTemplateSPOffset>::max() / scaleFactors.size() <
        (unsigned int) (w * w))
      throw SIMDException("MinWarpingTemplate::MinWarpingTemplate",
                          "size of scale-plane stack exceeds offset limit");
    initializeThresholds(scaleFactors, maxThreshold);
    createMinTemplate();
  }

  virtual ~MinWarpingTemplate() { destroyMinTemplate(); }

  template <typename ImgType, int SIMD_WIDTH, int SIMD_ALIGN>
  void getImage(SIMDImage<ImgType, SIMD_WIDTH, SIMD_ALIGN> &img, ImgType undef,
                double scale = 1.0, double shift = 0.0) const
  {
    lp_min_templt_col_t tCol;
    int jY, sp;

    img.resize(w, w);
    // filled with negative value for all undefined pixels
    img.fill(undef);
    for (int jX = 0; jX < w; jX++) {
      // get the template column
      tCol = minTemplateColumns[jX];
      if (tCol->len != 0)
        // go through all elements in this template column
        for (int colEntry = 0; colEntry < tCol->len; colEntry++) {
          jY          = tCol->entries()[colEntry].iy;
          sp          = tCol->entries()[colEntry].sp;
          img[jY][jX] = scale * (sp / (wSPS * w)) + shift;
        }
    }
  }

  // write in human-readable form
  void save(FILE *f) const
  {
    for (int jX = 0; jX < w; jX++) {
      // get the template column
      lp_min_templt_col_t tCol = minTemplateColumns[jX];
      fprintf(f, "x = %4d, len = %4d\n", jX, tCol->len);
      if (tCol->len != 0) {
        // go through all elements in this template column
        for (int colEntry = 0; colEntry < tCol->len; colEntry++) {
          int jY = tCol->entries()[colEntry].iy;
          int sp = tCol->entries()[colEntry].sp;
          fprintf(f, "\ty = %4d, sp = %4d\n", jY, sp / (wSPS * w));
        }
      }
    }
  }

  // this is a left-over from earlier experiments; sorts entries
  // according to y value, idea was to improve cache efficiency,
  // didn't help though;
  // if you use it, call this function only once!
  void sortMinTemplateY()
  {
    int ixright = w / 2 + 1;
    // we only need to sort right part of template
    for (int jX = ixright; jX < w; jX++) {
      // get the template column
      lp_min_templt_col_t tCol = minTemplateColumns[jX];
      // we need at least 3 entries to sort: 0 (is left untouched) and 2 more
      if (tCol->len > 2) {
        // go through all elements in this template column
        for (int colEntry1 = 1; colEntry1 < (tCol->len + 1) / 2; colEntry1++) {
          int colEntry2 = tCol->len - colEntry1;
          // swap entries
          min_templt_col_entry_t tmp = tCol->entries()[colEntry1];
          tCol->entries()[colEntry1] = tCol->entries()[colEntry2];
          tCol->entries()[colEntry2] = tmp;
        }
      }
    }
  }

  //--------------------------------------------------------------------------

protected:
  // thresholds are chosen in the center between neighboring scale factors
  void initializeThresholds(const std::vector<double> &scaleFactors,
                            double maxThreshold)
  {
    int nScalePlanes = scaleFactors.size();
    thresholds.resize(nScalePlanes);
    // initialize thresholds: center between scaleFactors
    for (int i = 0; i < nScalePlanes - 1; i++)
      thresholds[i] = 0.5 * (scaleFactors[i] + scaleFactors[i + 1]);
    // upper threshold is given
    thresholds[nScalePlanes - 1] = maxThreshold;
  }

  // creates scale-plane template for min-warping
  // from ScalePlaneStack.C (code by mk)
  // as we know, that we can access linear memory most efficiently, we
  // build our template this way:
  // <2*len0> <sp00> <iy00> ... <spX0> <iyX0>
  // <2*len1> <sp01> <iy01> ... <spX1> <iyX1>
  // ...
  // we must not omit the col_length (although it follows a simple rule)
  // because we need it to mark empty columns as they occur in templates with
  // odd widths!

  void createMinTemplate()
  {
    int ixleft  = (w - 1) / 2;
    int ixright = w / 2 + 1;

    // Reserve enough memory for the min template - there have to be
    // "w" column headers and at least "w * w / 4" column entries
    int templtSpace = sizeof(min_templt_col_t) * w +
                      sizeof(min_templt_col_entry_t) * (w * w / 4);

    // Reserve the complete template memory with malloc, as the
    // memory region is composed of min_templt_col_entry_t and
    // min_templt_col_entry_t - so we can't use "new"
    minTemplate = (lp_min_templt_col_t) (malloc(templtSpace));

    // Reserve the memory for minTemplateColumns - this array holds
    // pointers to the beginning of each column
    minTemplateColumns = new lp_min_templt_col_t[w];

    // Pointer on the current column in the min template
    lp_min_templt_col_t ptr = minTemplate;

    // Skip the first column - set the column pointer to the
    // beginning of the template
    ptr->len              = 0;
    minTemplateColumns[0] = ptr;
    // goto next column (empty "entries" array)
    ptr++;

    // Calculate the left half of the min template
    calculateMinTemplatePart(ptr, 1, ixleft + 1, 1);

    // we insert a gap here if the width is even
    if (!(w & 1)) {
      ptr->len                       = 0;
      minTemplateColumns[ixleft + 1] = ptr;
      // goto next column (empty "entries" array)
      ptr++;
    }

    // Calculate the right half of the min template
    calculateMinTemplatePart(ptr, ixright, w, (-1));

    // Check whether we're still inside the reserved memory region
    // printf("w = %d, %p - %p = %d <= %d ?\n",
    //	  w, ptr, minTemplate, (intptr_t)ptr - (intptr_t)minTemplate,
    //	  templtSpace);
    if ((intptr_t) ptr - (intptr_t) minTemplate > templtSpace)
      throw SIMDException("MinWarpingTemplate::createTemplate",
                          "insufficient reserved memory");
  }

  // calculateMinTemplatePart is used by createMinTemplate to create
  // the left and the right half of the minTemplate.
  //
  // ptr is a reference on the current column pointer
  // xstart is the start value for ix
  // xstop is the end value for ix, where xstop is not reached, ix != xstop
  // iy_dir is the direction into which iy gets incremented
  void calculateMinTemplatePart(lp_min_templt_col_t &ptr, const int xstart,
                                const int xstop, const int iy_dir)
  {
    int ix, iy, j, ej;
    int colHeight;
    unsigned int t;
    double x, y;
    double sigma, sinxpy;
    bool found;
    double rho;

    for (ix = xstart; ix < xstop; ix++) {
      // Set the column pointer
      minTemplateColumns[ix] = ptr;

      x         = ix * dTheta;
      colHeight = (std::abs(2 * ix - w) + 1) / 2;

      iy       = 0;
      ptr->len = 0;
      // entry index (with rho-mechanism: can differ from j!)
      ej = 0;
      for (j = 0; j < colHeight; j++) {
        if (iy < 0) iy += w;

        y      = iy * dTheta;
        sinxpy = sin(x + y);
        sigma  = sin(x) / sinxpy;

        // 16. Jul 12 (rm): we only accept (x,y) pairs within a
        // certain rho-range
        rho = sin(y) / sinxpy;
        // printf("rho = %g\n", rho);
        if ((rho >= rhoMin) && (rho <= rhoMax)) {
          // printf("inside\n");

          if (sigma <= thresholds[thresholds.size() - 1]) {
            found = false;
            for (t = 0; t < thresholds.size(); t++)
              if (sigma <= thresholds[t]) {
                found = true;
                break;
              }
            if (!found)
              throw SIMDException(
                "MinWarpingTemplate::calculateMinTemplatePart",
                "scale plane not found");

            // Increment the column length
            ptr->len++;

            // Calculate the scale-plane address offset (multiples of SPSType)
            // 23. Oct 12 (rm): here we can take any scale plane since
            // the offset is the same for all
            // ptr->entries()[ej].sp =
            //  (uint32_t)(scalePlane[t] - scalePlane[0]);
            // 18. Aug 14 (rm): no access to scalePlane after refactoring
            ptr->entries()[ej].sp = (MinTemplateSPOffset) (t * w * wSPS);

            // Store the row index
            ptr->entries()[ej].iy = iy;

            // 16. Jul 12 (rm) next index
            ej++;

          } else {
            break;
          }
        }

        // Calculate the next row index
        iy += iy_dir;
      }

      // Increment the column pointer position
      // advance by ptr->len elements of type min_templt_col_entry (skip
      // entire entries array) and 1 element of min_templt_col_t
      ptr =
        (lp_min_templt_col_t) ((lp_min_templt_col_entry_t) ptr + ptr->len) + 1;
    }
  }

  void destroyMinTemplate()
  {
    if (minTemplate) {
      free(minTemplate);
      minTemplate = NULL;
    }
    if (minTemplateColumns) {
      delete[] minTemplateColumns;
      minTemplateColumns = 0;
    }
  }
};

} // namespace simd

#endif
