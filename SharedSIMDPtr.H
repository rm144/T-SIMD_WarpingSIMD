// ===========================================================================
// 
// SharedSIMDPtr.H --
// smart shared pointer class for SIMD data
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

#ifndef _SHARED_SIMD_PTR_H_
#define _SHARED_SIMD_PTR_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "SIMDException.H"

#ifdef SHARED_SIMD_PTR_DEBUG_ON
#define SHARED_SIMD_PTR_DEBUG(ACTION) ACTION
#else
#define SHARED_SIMD_PTR_DEBUG(ACTION)
#endif

// NOTES: 
// - SharedSIMDPtr is NOT thread-safe (no atomic counters in
//   SharedSIMDPtrDesc, no atomic manipulator of pointers etc.)

// TODO: intrusive reference counting (Alexandrescu p.167)
// TODO: (allocate common block for data and desc)
// TODO: (alignment of desc needs to be considered!)
// TODO: thread-safety (atomic counting, atomic pointer manipulation...)

namespace ns_simd {
  
  // -------------------------------------------------------------------------
  // SharedSIMDPtrDesc
  // -------------------------------------------------------------------------

  // use this in .C file to initialize static members
#ifdef SHARED_SIMD_PTR_HEAP_COUNT
#define SHARED_SIMD_PTR_HEAP_COUNT_INIT	     \
  size_t SharedSIMDPtrDesc::totalDescNo = 0; \
  size_t SharedSIMDPtrBase::totalSize = 0;
#endif
  
  class SharedSIMDPtrDesc
  {
    
  protected:

    // number of references to shared data:
    // refCount == 0: no data, not shared ("empty")
    // refCount == 1: data, not shared ("unique")
    // refCount > 1:  data, shared ("shared")
    // (note that in SharedSIMDPtr, the SharedSIMDPtrDesc is deleted if empty,
    //  so the case refCount==0 should never occur)
    int refCount;
    
  public:
    
#ifdef SHARED_SIMD_PTR_HEAP_COUNT
    static size_t totalDescNo;
#endif
    
    SharedSIMDPtrDesc() 
    { 
      SHARED_SIMD_PTR_DEBUG(puts("SharedSIMDPtrDesc(): created desc"));
      init();
#ifdef SHARED_SIMD_PTR_HEAP_COUNT
      totalDescNo++;
#endif
    }
    ~SharedSIMDPtrDesc()
    {
      SHARED_SIMD_PTR_DEBUG(puts("~SharedSIMDPtrDesc(): destroyed desc"));
#ifdef SHARED_SIMD_PTR_HEAP_COUNT
      assert(totalDescNo > 0);
      totalDescNo--;
#endif
    }
    void init() { refCount = 1; }
    void inc()
    {
      refCount++;
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtrDesc::inc() refCount=%d\n",
				   refCount));
    }
    void dec() {
      //  1. Mar 18 (rm): make sure refCount never goes to zero
      assert(refCount > 1);
      refCount--;
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtrDesc::dec() refCount=%d\n",
				   refCount));
    }
    // bool empty()  { return (refCount == 0); }
    bool unique() { return (refCount == 1); }
    bool shared() { return (refCount > 1); }
    int get() { return refCount; }
  };
      
  // -------------------------------------------------------------------------
  // SharedSIMDPtr
  // -------------------------------------------------------------------------
  
  // SharedSIMDPtr behaves (at least partly) like a pointer:
  // - copy-constructor and operator= copy the data pointer, not the data
  //   (implemented through makeView);
  // making a view of a "nil" SharedSIMDPtr just creates a "nil" SharedSIMDPtr;
  // beware if you use a std::vector of SharedSIMDPtr:
  // - std::vector.resize(n) calls SharedSIMDPtr() once (default argument)
  //   and then n times SharedSIMDPtr(const SharedSIMDPtr&),
  //   so we create n independent SharedSIMDPtr's;
  // - std::vector.resize(n, sharedPtr) for a non-nil sharedPtr creates
  //   n views of the same data!

#ifdef SHARED_SIMD_PTR_HEAP_COUNT
  class SharedSIMDPtrBase
  {
  public:
    static size_t totalSize;
  };
#endif

  template <typename T, int SIMD_WIDTH, int SIMD_ALIGN>
  class SharedSIMDPtr
#ifdef SHARED_SIMD_PTR_HEAP_COUNT
    : public SharedSIMDPtrBase
#endif
  {
  protected:
    
    // origData: original data block, may be shared
    // can be NULL pointer
    T *origData;
    // descriptor shared between all shared pointers
    SharedSIMDPtrDesc *desc;
    // origOwner == true indicates that this object is the original owner
    // (creator) of origData
    // this may not be strictly necessary, it just guarantees that only the
    // original owner may reallocate data (not the last one referring to it)
    bool origOwner;
    // 13. Sep 19 (rm)
    // if viewIsProtected is true, then it is not possible to make this object
    // a view of something else (required e.g. for planes of a SIMDStack)
    bool viewIsProtected;
    
  public:
    
    // data: shared data; can be NULL pointer
    T* data;
    // size of the data behind "data"
    size_t size;
    
    // create a shared NULL pointer
    SharedSIMDPtr()
      : viewIsProtected(false)
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr()\n"));
      allocate(0);
    }
    
    // create a shared pointer with data (this is the original owner)
    // must tolerate initSize == 0 (e.g. copy constructor on empty object
    // called from vector::resize)
    SharedSIMDPtr(size_t initSize)
      : viewIsProtected(false)
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr(%ld)\n", (long int)initSize));
      allocate(initSize);
    }

    // TODO: copy constructor creates a view (as operator=), ok?
    SharedSIMDPtr(const SharedSIMDPtr<T,SIMD_WIDTH,SIMD_ALIGN> &ssp)
       : viewIsProtected(false)
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr(SharedSIMDPtr&)\n"));
      makeView(ssp, 0, ssp.size);
    }
    
    // free shared pointer if it is not shared
    virtual ~SharedSIMDPtr()
    {
      SHARED_SIMD_PTR_DEBUG(puts("~SharedSIMDPtr()"));
      release();
    }
    
    // change size
    void resize(size_t newSize)
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr::resize(%d)\n", newSize));
      reallocate(newSize);
    }
    
    // access to the data portion
    // operator T*() const { return data; }
    // not recommended, see Alexandrescu "Modern C++ design" p.172
 
    // -> access
    T* operator->()
    {
      if (!data)
	throw SIMDException("SharedSIMDPtr::operator->",
			    "dereferencing of nil pointer");
      return data;
    }

    // * access
    // TODO: is SharedSIMDPtr::operator*  needed somewhere?
    /*
    T& operator*()
    {
      return *data;
    }
    */
    
    // TODO: methods for comparison of SharedSIMDPtrs?, see Alexandrescu
    
    // different view
    void view(const SharedSIMDPtr<T,SIMD_WIDTH,SIMD_ALIGN> &ssp,
	      size_t portionOffset, size_t portionSize)
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr::view(ssp, %ld, %ld)\n",
				   (long int) portionOffset, 
				   (long int) portionSize));
      makeView(ssp, portionOffset, portionSize);
    }

    // different view, protected
    void protectedView(const SharedSIMDPtr<T,SIMD_WIDTH,SIMD_ALIGN> &ssp,
		       size_t portionOffset, size_t portionSize)
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr::protectedView(ssp, %ld, %ld)\n",
				   (long int) portionOffset, 
				   (long int) portionSize));
      makeProtectedView(ssp, portionOffset, portionSize);
    }
    
    // identical view
    SharedSIMDPtr<T,SIMD_WIDTH,SIMD_ALIGN>&
    operator=(const SharedSIMDPtr<T,SIMD_WIDTH,SIMD_ALIGN> &ssp)
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr::operator=(ssp)\n"));
      makeView(ssp, 0, ssp.size);
      return *this;
    }
    
    void fill(T val)
    {
      if (data)
	for (size_t i = 0; i < size; i++) 
	  data[i] = val;
    }

    void fillByte(int byte)
    {
      if (data) 
	memset(data, byte, size * sizeof(T));
    }

    void clear()
    {
      fillByte(0x00);
    }
    
    bool save(FILE *f) const
    {
      if (!data) return false;
      size_t res = fwrite(data, sizeof(T), size, f);
      return (res == size);
    }
    
    bool load(FILE *f)
    {
      if (!data) return false;
      size_t res = fread(data, sizeof(T), size, f);
      return (res == size);
    }

    // copy data from a source if source and destination exist 
    void copyFrom(const T *const srcData)
    {
      if (data && srcData) 
	memcpy(data, srcData, size * sizeof(T));
    }
    
    // copy data from a source with conversion if source and destination exist
    template <typename U>
    void copyFrom(const U *const srcData)
    {
      if (data && srcData)
	for (size_t i = 0; i < size; i++)
	  data[i] = T(srcData[i]);
    }
    
    // copy data to a destination if source and destination exist
    void copyTo(T *const dstData) const
    {
      if (data && dstData)
	memcpy(dstData, data, size * sizeof(T));
    }
    
    // copy data to a destination with conversion if source and dest. exist
    template <typename U>
    void copyTo(U *const dstData) const
    {
      if (data && dstData)
	for (size_t i = 0; i < size; i++)
	  dstData[i] = U(data[i]);
    }
    
    void printInfo(FILE *f = stdout)
    {
      fprintf(f, "origData %p, data %p, size %ld, origOwner %d, desc %p",
	      (void*) origData, (void*) data, (long int) size, origOwner, 
	      (void*) desc);
      if (desc)
	fprintf(f, " (refCount %d)\n", desc->get());
      else
	fprintf(f, "\n");
    }

  protected:
    
    static void checkSize(size_t sizeToCheck)
    {
      // size must be a multiple of SIMD_WIDTH
      if ((sizeToCheck * sizeof(T)) % SIMD_WIDTH != 0)
	throw SIMDException("SharedSIMDPtr::checksize",
			    "size (in bytes) is not a multiple of SIMD_WIDTH");
    }
    
    static void checkAlignment(const T *const ptr)
    {
      // (SIMD_ALIGN is assumed to be a power of two, not checked)
      if ((((uintptr_t) ptr) & (SIMD_ALIGN - 1)) != 0)
	throw SIMDException("SharedSIMDPtr::checkAlignment",
			    "pointer is not aligned (SIMD_ALIGN)");
    }

    // allocate must tolerate allocSize == 0
    // (called from constructor (2) and reallocate())
    void allocate(size_t allocSize)
    {
      origOwner = true;
      if (allocSize == 0) {
	// empty pointer
	origData = data = 0;
	size = 0;
	desc = 0;
      }
      else {
	// non-empty pointer
	checkSize(allocSize);
	// alignment passed to posix_memalign must be a multiple of
	// sizeof(void*) according to manpage, but we also allow values
	// for SIMD_ALIGN smaller than sizeof(void*) (e.g. 1)
	size_t align = std::max(sizeof(void*), (size_t) SIMD_ALIGN);
	if (align % sizeof(void*) != 0)
	  throw SIMDException("SharedSIMDPtr::allocate",
			      "alignment must be smaller than "
			      "or multiple of sizeof(void*)");
	if (posix_memalign((void**) &origData, align, allocSize * sizeof(T)))
	  throw SIMDException("SharedSIMDPtr::allocate",
			      "posix_memalign failed");
	SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr::allocate(): "
				     "created data %p\n", (void*) origData));
#ifdef SHARED_SIMD_PTR_HEAP_COUNT
	totalSize += allocSize * sizeof(T);
#endif
	desc = new SharedSIMDPtrDesc();
	data = origData;
	size = allocSize;
      }
    }

    // (called from reallocate() (if origOwner) and release() (if last user))
    void deallocate()
    {
      // 1. Mar 18 (rm): just to be on the safe side: make sure we never
      // delete shared data
      assert(!desc || desc->unique());
      if (origData) {
	free(origData);
	SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr::deallocate(): "
				     "destroyed data %p\n", (void*) origData));
#ifdef SHARED_SIMD_PTR_HEAP_COUNT
	assert(totalSize >= size * sizeof(T));
	totalSize -= size * sizeof(T);
#endif
      }
      if (desc)
	delete desc;
    }

    // reallocate() (called from resize())
    //              newSize
    // newSize size ==size origOwner desc action
    // -----------------------------------------------------------------
    // 0       0    true   true      0    nothing (newSize == size)
    // !=0     0    false  true      0    deallocate-, allocate(!=0)
    // 0       !=0  false  true      !=0  shared-chk, deallocate+, allocate(0)
    // 0       !=0  false  false     !=0  exception: not original owner!
    // !=0     !=0  false  true      !=0  shared-chk, deallocate+, allocate(!=0)
    // !=0     !=0  false  false     !=0  exception: not original owner
    // !=0     !=0  true   true      !=0  nothing (newSize == size)
    // !=0     !=0  true   false     !=0  nothing (newSize == size)
    void reallocate(size_t newSize)
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr::reallocate(%ld)\n", 
				   (long int) newSize));
      // no need to resize of the size didn't change
      if (newSize != size) {
	// if we are not the original owner, we are not allowed to reallocate
	if (!origOwner)
	  throw SIMDException("SharedSIMDPtr::reallocate",
			      "only original data owner may reallocate");
	// here we are the original owner
	if (desc)
	  // shared data cannot be reallocated
	  if (desc->shared())
	    throw SIMDException("SharedSIMDPtr::reallocate",
				"shared pointers may not be reallocated");
	// free memory if it was allocated
	deallocate();
	// allocate new memory
	allocate(newSize);
      }
    }
    
    // (called from destructor and makeView())
    void release()
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr::release()\n"));
      // we only free memory if the data is unique (non-NULL and not shared)
      if (desc) {
	if (desc->unique())
	  deallocate();
	// otherwise we just indicate that a reference has been deleted
	else
	  desc->dec();
      }
    }

    void makeView(const SharedSIMDPtr<T,SIMD_WIDTH,SIMD_ALIGN> &ssp,
		  size_t portionOffset, size_t portionSize)
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr::makeView(ssp, %ld, %ld)\n",
				   (long int) portionOffset, 
				   (long int) portionSize));
      // 13. Sep 19 (rm)
      if (viewIsProtected)
	throw SIMDException("SharedSIMDPtr::makeView",
			    "view is protected");
      _makeView(ssp, portionOffset, portionSize);
    }

    // 13. Sep 19 (rm)
    void makeProtectedView(const SharedSIMDPtr<T,SIMD_WIDTH,SIMD_ALIGN> &ssp,
			   size_t portionOffset, size_t portionSize)
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr::makeProtectedView(ssp, %ld, %ld)\n",
				   (long int) portionOffset, 
				   (long int) portionSize));
      if (viewIsProtected)
	throw SIMDException("SharedSIMDPtr::makeProtectedView",
			    "view is protected");
      _makeView(ssp, portionOffset, portionSize);
      viewIsProtected = true;
    }

    // 13. Sep 19 (rm): makeView -> _makeView
    void _makeView(const SharedSIMDPtr<T,SIMD_WIDTH,SIMD_ALIGN> &ssp,
		   size_t portionOffset, size_t portionSize)
    {
      SHARED_SIMD_PTR_DEBUG(printf("SharedSIMDPtr::_makeView(ssp, %ld, %ld)\n",
				   (long int) portionOffset, 
				   (long int) portionSize));
      // 1. Mar 18 (rm): happened in both branches, moved here
      // possibly release memory if we are the last owner
      // 
      release();
      // after release() we have two cases:
      // - the data was deallocated (last user) or
      // - there are other users (shared) or at least one (unique) who will
      //   release the data later
      if (ssp.size == 0) {
	// view of an empty pointer is an empty pointer
	if ((portionOffset != 0) || (portionSize != 0))
	  throw SIMDException("SharedSIMDPtr::makeView",
			      "view of empty object must be empty");
	origData = data = 0;
	size = 0;
	origOwner = true;
	desc = 0;
      }
      else {
	// check size constraints
	if (portionOffset + portionSize > ssp.size)
	  throw SIMDException("SharedSIMDPtr::makeView",
			      "data boundary violation");
	// now make "this" a view of the data in ssp
	origData = ssp.origData;
	data = ssp.data + portionOffset;
	checkAlignment(data);
	size = portionSize;
	origOwner = false;
	desc = ssp.desc;
	desc->inc();
      }
    }
    
  };

}

#endif
