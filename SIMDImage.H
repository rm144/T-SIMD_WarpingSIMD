// ===========================================================================
//
// SIMDImage.H --
// simple image classes used for implementation of SIMD computations
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
#ifndef SIMD_IMAGE_H_
#define SIMD_IMAGE_H_

#include "SIMDException.H"
#include "SharedSIMDPtr.H"

#include <cstddef>
#include <stdio.h>
#include <vector>

#ifdef SIMD_IMAGE_DEBUG_ON
#define SIMD_IMAGE_DEBUG(ACTION) ACTION
#else
#define SIMD_IMAGE_DEBUG(ACTION)
#endif

namespace simd {

// ===========================================================================
// checkWidth
// ===========================================================================

template <typename T, size_t SIMD_WIDTH, size_t SIMD_ALIGN>
void checkWidth(int width)
{
  SIMD_IMAGE_DEBUG(printf("checkWidth(%d)\n", width));
  // width should be a multiple of SIMD_WIDTH
  if ((width * sizeof(T)) % SIMD_WIDTH != 0)
    throw SIMDException("SIMDImage::checkWidth",
                        "width is not a multiple of SIMD_WIDTH");
  // assert that each row start is aligned (alignment of data is
  // guaranteed by SIMDData)
  // (for simplicity, it would also reject a single-row image where
  // this condition is not fulfilled, even though there actually is
  // only a single row)
  // (SIMD_ALIGN is assumed to be a power of two, not checked)
  // (special case: width == 0: accepted)
  if (((width * sizeof(T)) & (SIMD_ALIGN - 1)) != 0)
    throw SIMDException("SIMDImage::checkRowAlignment",
                        "row start is not aligned");
}

// ===========================================================================
// AddOns for SIMDImage
// ===========================================================================

// - AddOns are additional data that can be part of a SIMDImage
// - it is possible to assign any AddOn to DefaultAddOn
// - it is not possible to assign DefaultAddOn to Panorama

// this is an empty AddOn
class DefaultAddOn
{
public:
  DefaultAddOn() {}
  DefaultAddOn(const DefaultAddOn &) {}
  DefaultAddOn &operator=(const DefaultAddOn &) { return *this; }

  // DefaultAddOn accepts all other AddOns (but not the other way round)
  template <class OtherAddOn>
  DefaultAddOn(const OtherAddOn &)
  {}
  template <class OtherAddOn>
  DefaultAddOn &operator=(const OtherAddOn &)
  {
    return *this;
  }

  bool operator==(const DefaultAddOn &) const { return true; }
  bool operator!=(const DefaultAddOn &) const { return false; }

  // properties of this AddOn
  void isInvariantToVerticalReorganization() const {}
  void isInvariantToHorizontalReorganization() const {}
};

// AddOn for panoramic images
class Panorama
{
public:
  double verticalResolution; // [rad/pixel]
  double horizon;            // relates to vertical image coordinate

  // TODO: 3. Jul 18 (rm) this is problematic but not easy to fix:
  // TODO: if e.g. resize is called and no addOn parameter is passed, the
  // TODO: default constr. of Panorama is used, which sets the parameters
  // TODO: to zero; presently this default constr. can't be removed because
  // TODO: it is used in several locations
  Panorama() : verticalResolution(0.0), horizon(0.0) {}

  Panorama(const Panorama &p)
    : verticalResolution(p.verticalResolution), horizon(p.horizon)
  {}
  Panorama(double verticalResolution, double horizon)
    : verticalResolution(verticalResolution), horizon(horizon)
  {}
  Panorama &operator=(const Panorama &p)
  {
    verticalResolution = p.verticalResolution;
    horizon            = p.horizon;
    return *this;
  }
  bool operator==(const Panorama &p) const
  {
    return (verticalResolution == p.verticalResolution) &&
           (horizon == p.horizon);
  }
  bool operator!=(const Panorama &p) const { return !(*this == p); }

  // properties of this AddOn
  void isInvariantToHorizontalReorganization() const {}
};

// ===========================================================================
// SIMDImage
// ===========================================================================

// SIMDImage is derived from SharedSIMDPtr instead of containing a
// SharedSIMDPtr (to avoid another layer of data access:
// myImage.data instead of myImage.sharedPtr.data), but it does not
// behave like a pointer: copy-constructor and operator= allocate new data!

// SIMDImage and SIMDStack are both derived from SharedSIMDPtr,
// therefore they can both be passed to
// SIMDImage::view(SharedSIMDPtr,size_t,size_t)

// AddOns can be added to SIMDImage, e.g. Panorama (with parameters
// of panoramic images)

// simple image class template for single-band images
template <typename T, size_t SIMD_WIDTH, size_t SIMD_ALIGN,
          class AddOn = DefaultAddOn>
class SIMDImage : public SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>
{
public:
  // ***** variables are public, but should not be changed from outside *****
  // width, height
  int w, h;

  AddOn addOn;

  SIMDImage() : SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>(), w(0), h(0), addOn()
  {
    SIMD_IMAGE_DEBUG(puts("SIMDImage()"));
  }

  SIMDImage(int initW, int initH, const AddOn &initAddOn = AddOn())
    : SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>(initW * initH), w(initW),
      h(initH), addOn(initAddOn)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDImage(%d,%d)\n", initW, initH));
    checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
  }

  // 1. Mar 18 (rm): we have to add this, since otherwise it seems that
  // the constructor of SharedSIMDPtr fits better and is used without using
  // a SIMDImage constructor
  SIMDImage(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &si)
    : SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>(si.size), w(si.w), h(si.h),
      addOn(si.addOn)
  {
    SIMD_IMAGE_DEBUG(puts("SIMDImage<...,AddOn>(si&)"));
    checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(si.data);
  }

  // note that copy constructor is called from vector::resize,
  // therefore SharedSIMDPtr must tolerate zero size (copies from
  // empty SIMDImage)
  // TODO: is this copy ctor covered by the one immediately below?
  template <class OtherAddOn>
  SIMDImage(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, OtherAddOn> &si)
    : SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>(si.size), w(si.w), h(si.h),
      addOn(si.addOn)
  {
    SIMD_IMAGE_DEBUG(puts("SIMDImage<...,OtherAddOn>(si&)"));
    checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(si.data);
  }

  // required for things like:
  // SIMDImage<char,16,16> ic;
  // SIMDImage<float,16,16> if = ic;
  // (copy constructor gets called, not operator=)
  template <typename U, size_t OTHER_SIMD_WIDTH, size_t OTHER_SIMD_ALIGN,
            class OtherAddOn>
  SIMDImage(
    const SIMDImage<U, OTHER_SIMD_WIDTH, OTHER_SIMD_ALIGN, OtherAddOn> &si)
    : SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>(si.size), w(si.w), h(si.h),
      addOn(si.addOn)
  {
    SIMD_IMAGE_DEBUG(puts("SIMDImage(si<U,...,OtherAddOn>&)"));
    checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(si.data);
  }

  ~SIMDImage() { SIMD_IMAGE_DEBUG(puts("~SIMDImage()")); }

  void view(const SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN> &sp, size_t extOff,
            int extW, int extH, const AddOn &extAddOn = AddOn())
  {
    SIMD_IMAGE_DEBUG(printf("SIMDImage::view(%lu,%d,%d)\n",
                            (unsigned long) extOff, extW, extH));
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::view(sp, extOff, extW * extH);
    w     = extW;
    h     = extH;
    addOn = extAddOn;
    checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
  }

  void view(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &si)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDImage::view(&si)\n"));
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::view(si, 0, si.size);
    w     = si.w;
    h     = si.h;
    addOn = si.addOn;
    // check width not necessary
  }

  // 13. Sep 19 (rm)
  void protectedView(const SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN> &sp,
                     size_t extOff, int extW, int extH,
                     const AddOn &extAddOn = AddOn())
  {
    SIMD_IMAGE_DEBUG(printf("SIMDImage::protectedView(%lu,%d,%d)\n",
                            (unsigned long) extOff, extW, extH));
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::protectedView(sp, extOff,
                                                            extW * extH);
    w     = extW;
    h     = extH;
    addOn = extAddOn;
    checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
  }

  // 13. Sep 19 (rm)
  void protectedView(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &si)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDImage::protectedView(&si)\n"));
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::protectedView(si, 0, si.size);
    w     = si.w;
    h     = si.h;
    addOn = si.addOn;
    // check width not necessary
  }

  //  1. Mar 18 (rm): this one is needed since otherwise operator= of
  // SharedSIMDPtr is used instead (seems to fit better)
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &operator=(
    const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &si)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDImage::operator=<...,AddOn>\n"));
    resize(si.w, si.h, si.addOn);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(si.data);
    return *this;
  }

  // apparently we need this to disable the default operator=
  // even though it is a special case of the templated version below
  template <class OtherAddOn>
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &operator=(
    const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, OtherAddOn> &si)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDImage::operator=<...,OtherAddOn>\n"));
    resize(si.w, si.h, si.addOn);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(si.data);
    return *this;
  }

  template <typename U, size_t OTHER_SIMD_WIDTH, size_t OTHER_SIMD_ALIGN,
            class OtherAddOn>
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &operator=(
    const SIMDImage<U, OTHER_SIMD_WIDTH, OTHER_SIMD_ALIGN, OtherAddOn> &si)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDImage::operator=<U,...,OtherAddOn>\n"));
    resize(si.w, si.h, si.addOn);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(si.data);
    return *this;
  }

  // resize doesn't clear!
  template <class OtherAddOn>
  void resize(int newW, int newH, const OtherAddOn &newAddOn)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDImage::resize(<...,OtherAddOn>)\n"));
    if (!SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::origOwner) {
      // if we don't own the data, we can only proceed if size is
      // correct (size change is not possible)
      if (!sameSize(newW, newH))
        throw SIMDException("SIMDImage::resize",
                            "data not owned, size change not possible");
    } else {
      // if we own the data, we only have to reallocate if the size
      // changes
      if (!sameSize(newW, newH)) {
        SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::reallocate(newW * newH);
        w     = newW;
        h     = newH;
        addOn = newAddOn;
        checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
      }
    }
  }

  void resize(int newW, int newH)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDImage::resize()\n"));
    resize(newW, newH, AddOn());
  }

  bool sameSize(int otherW, int otherH) const
  {
    return (otherW == w) && (otherH == h);
  }

  template <class OtherAddOn>
  bool sameSize(
    const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, OtherAddOn> &si) const
  {
    return sameSize(si.w, si.h);
  }

  template <class OtherAddOn>
  bool same(const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, OtherAddOn> &si) const
  {
    return sameSize(si) &&
           (SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::data == si.data);
  }

  // access to row j without index check
  // (slow since we have no table with row pointers)
  T *operator[](int j) const
  {
    if (!SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::data)
      throw SIMDException("SIMDImage::operator[]", "empty image");
    return SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::data + j * w;
  }

  void printInfo(FILE *f = stdout)
  {
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::printInfo(f);
    fprintf(f, "w %d, h %d\n", w, h);
  }

  bool validIndex(int x, int y)
  {
    return ((x >= 0) && (x < w) && (y >= 0) && (y < h));
  }
};

// ===========================================================================
// SIMDStack
// ===========================================================================

// TODO: also provide view for SIMDStack as in SIMDImage
// TODO: such that another class can split its data into multiple SIMDStacks?
// TODO: generalize SIMDImage and SIMDStack in a template?
// TODO: test SIMDStack for non-default addon!

// TODO: 10. Jul 19 (rm)
// TODO: - provide view for SIMDStack, may be difficult

// template class for stack of images (of same size)
template <typename T, size_t SIMD_WIDTH, size_t SIMD_ALIGN,
          class AddOn = DefaultAddOn>
class SIMDStack : public SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>
{
public:
  // ***** variables are public, but should not be changed from outside *****
  // number of planes, width and height of each plane
  int numPlanes, w, h;
  AddOn addOn;

  SIMDStack()
    : SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>(), numPlanes(0), w(0), h(0),
      addOn()
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack()\n"));
  }

  /* covered by constructor below 10. Jul 19 (rm)
  SIMDStack(int initNumPlanes, int initW, int initH)
    : SharedSIMDPtr<T,SIMD_WIDTH,SIMD_ALIGN>(initNumPlanes * initW * initH),
      numPlanes(initNumPlanes), w(initW), h(initH)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack(%d,%d,%d)\n",
                            initNumPlanes,initW,initH));
    checkWidth<T,SIMD_WIDTH,SIMD_ALIGN>(w);
    allocatePlanes(AddOn());
  }
  */

  SIMDStack(int initNumPlanes, int initW, int initH,
            const AddOn &initAddOn = AddOn())
    : SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>(initNumPlanes * initW * initH),
      numPlanes(initNumPlanes), w(initW), h(initH), addOn(initAddOn)
  {
    SIMD_IMAGE_DEBUG(
      printf("SIMDStack(%d,%d,%d)\n", initNumPlanes, initW, initH));
    checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
    allocatePlanes(initAddOn);
  }

  SIMDStack(const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &ss)
    : SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>(ss.size),
      numPlanes(ss.numPlanes), w(ss.w), h(ss.h), addOn(ss.addOn)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack(SIMDStack&)\n"));
    checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
    allocatePlanes(ss.addOn);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(ss.data);
  }

  template <class OtherAddOn>
  SIMDStack(const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, OtherAddOn> &ss)
    : SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>(ss.size),
      numPlanes(ss.numPlanes), w(ss.w), h(ss.h), addOn(ss.addOn)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack(SIMDStack&)\n"));
    checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
    allocatePlanes(ss.addOn);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(ss.data);
  }

  template <typename U, size_t OTHER_SIMD_WIDTH, size_t OTHER_SIMD_ALIGN,
            class OtherAddOn>
  SIMDStack(
    const SIMDStack<U, OTHER_SIMD_WIDTH, OTHER_SIMD_ALIGN, OtherAddOn> &ss)
    : SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>(ss.size),
      numPlanes(ss.numPlanes), w(ss.w), h(ss.h), addOn(ss.addOn)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack(SIMDStack<U>&)\n"));
    checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
    allocatePlanes(ss.addOn);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(ss.data);
  }

  ~SIMDStack() { SIMD_IMAGE_DEBUG(puts("~SIMDStack()")); }

  // apparently we need this to disable the default operator=
  // even though it is a special case of the templated version below
  SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &operator=(
    const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &ss)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack::operator=()\n"));
    resize(ss.numPlanes, ss.w, ss.h, ss.addOn);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(ss.data);
    return *this;
  }

  template <class OtherAddOn>
  SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &operator=(
    const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, OtherAddOn> &ss)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack::operator=(T,OtherAddOn)\n"));
    resize(ss.numPlanes, ss.w, ss.h, ss.addOn);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(ss.data);
    return *this;
  }

  template <typename U, size_t OTHER_SIMD_WIDTH, size_t OTHER_SIMD_ALIGN,
            class OtherAddOn>
  SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &operator=(
    const SIMDStack<U, OTHER_SIMD_WIDTH, OTHER_SIMD_ALIGN, OtherAddOn> &ss)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack::operator=(U,OtherAddOn)\n"));
    resize(ss.numPlanes, ss.w, ss.h, ss.addOn);
    SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::copyFrom(ss.data);
    return *this;
  }

  // resize doesn't clear!
  template <class OtherAddOn>
  void resize(int newNumPlanes, int newW, int newH, const OtherAddOn &addOn)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack::resize(OtherAddOn)\n"));
    if (!sameSize(newNumPlanes, newW, newH)) {
      SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::reallocate(newNumPlanes * newW *
                                                           newH);
      numPlanes = newNumPlanes;
      w         = newW;
      h         = newH;
      checkWidth<T, SIMD_WIDTH, SIMD_ALIGN>(w);
      allocatePlanes(addOn);
    }
  }

  void resize(int newNumPlanes, int newW, int newH)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack::resize()\n"));
    resize(newNumPlanes, newW, newH, AddOn());
  }

  // contributed by B. Volkmer
  template <typename U, size_t OTHER_SIMD_WIDTH, size_t OTHER_SIMD_ALIGN,
            class OtherAddOn>
  void resize(SIMDStack<U, OTHER_SIMD_WIDTH, OTHER_SIMD_ALIGN, OtherAddOn> o)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack::resize(U,OtherAddOn)\n"));
    resize(o.numPlanes, o.w, o.h, o.addOn);
  }

  // access to plane i
  inline SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &operator[](int i)
  {
    if (!((i >= 0) && (i < numPlanes)))
      throw SIMDException("SIMDStack::operator[]", "index out of range");
    return planes[i];
  }

  // access to plane i
  inline const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn> &operator[](
    int i) const
  {
    if (!((i >= 0) && (i < numPlanes)))
      throw SIMDException("SIMDStack::operator[]", "index out of range");
    return planes[i];
  }

  // access to all bands as single image
  inline SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &operator()() { return all; }

  // access to all bands as single image
  inline const SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> &operator()() const
  {
    return all;
  }

  bool sameSize(int otherNumPlanes, int otherW, int otherH) const
  {
    return (numPlanes == otherNumPlanes) && (w == otherW) && (h == otherH);
  }

  bool sameSize(const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN> &ss) const
  {
    return sameSize(ss.numPlanes, ss.w, ss.h);
  }

  // contributed by B. Volkmer
  template <class OtherAddOn>
  bool same(const SIMDStack<T, SIMD_WIDTH, SIMD_ALIGN, OtherAddOn> &si) const
  {
    return sameSize(si) &&
           (SharedSIMDPtr<T, SIMD_WIDTH, SIMD_ALIGN>::data == si.data);
  }

  bool validIndex(int plane, int x, int y)
  {
    return ((plane >= 0) && (plane < numPlanes) && (x >= 0) && (x < w) &&
            (y >= 0) && (y < h));
  }

protected:
  // image representing all planes, always DefaultAddOn
  SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN> all;

  // array of images representing the individual planes
  std::vector<SIMDImage<T, SIMD_WIDTH, SIMD_ALIGN, AddOn>> planes;

  void allocatePlanes(const AddOn &addOn)
  {
    SIMD_IMAGE_DEBUG(printf("SIMDStack::allocatePlanes()\n"));
    SIMD_IMAGE_DEBUG(printf("numPlanes %d, %d x %d\n", numPlanes, w, h));
    planes.resize(numPlanes);
    SIMD_IMAGE_DEBUG(
      printf("planes.size = %lu\n", (unsigned long) planes.size()));
    const size_t planeSize = w * h;
    size_t off             = 0;
    for (int i = 0; i < numPlanes; i++, off += planeSize) {
      SIMD_IMAGE_DEBUG(printf("creating plane %i, off = %lu, %d x %d\n", i,
                              (unsigned long) off, w, h));
      // shallow copy in each image!
      // 13. Sep 19 (rm): protected!
      planes[i].protectedView(*this, off, w, h, addOn);
    }
    // 23. Sep 16 (rm): fixed bug: was numPlanes*h x w
    SIMD_IMAGE_DEBUG(
      printf("creating all, off = 0, %d x %d\n", w, numPlanes * h));
    // shallow copy! always DefaultAddOn
    // 13. Sep 19 (rm): protected!
    all.protectedView(*this, 0, w, numPlanes * h);
  }
};

} // namespace simd

#endif
