// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyCore.h 
//! \author Cem Yuksel
//! 
//! \brief  Core functions and macros
//! 
//! Core functions and macros for math and other common operations
//! 
//-------------------------------------------------------------------------------
//
// Copyright (c) 2016, Cem Yuksel <cem@cemyuksel.com>
// All rights reserved.
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy 
// of this software and associated documentation files (the "Software"), to deal 
// in the Software without restriction, including without limitation the rights 
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
// copies of the Software, and to permit persons to whom the Software is 
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
// SOFTWARE.
// 
//-------------------------------------------------------------------------------

#ifndef _CY_CORE_H_INCLUDED_
#define _CY_CORE_H_INCLUDED_

//-------------------------------------------------------------------------------

#ifndef _CY_CORE_MEMCPY_LIMIT
#define _CY_CORE_MEMCPY_LIMIT 256
#endif

//-------------------------------------------------------------------------------

#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <type_traits>
#include <limits>

#if !defined(CY_NO_INTRIN_H) && !defined(CY_NO_EMMINTRIN_H)
# include <intrin.h>
#endif

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////
// Compiler compatibility
//////////////////////////////////////////////////////////////////////////

#if defined(__INTEL_COMPILER)
# define _CY_COMPILER_INTEL __INTEL_COMPILER
# define _CY_COMPILER_VER_MEETS(msc,gcc,clang,intel) _CY_COMPILER_INTEL >= intel
# define _CY_COMPILER_VER_BELOW(msc,gcc,clang,intel) _CY_COMPILER_INTEL <  intel
#elif defined(__clang__)
# define _CY_COMPILER_CLANG (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)
# define _CY_COMPILER_VER_MEETS(msc,gcc,clang,intel) _CY_COMPILER_CLANG >= clang
# define _CY_COMPILER_VER_BELOW(msc,gcc,clang,intel) _CY_COMPILER_CLANG <  clang
#elif defined(_MSC_VER)
# define _CY_COMPILER_MSC _MSC_VER
# define _CY_COMPILER_VER_MEETS(msc,gcc,clang,intel) _CY_COMPILER_MSC >= msc
# define _CY_COMPILER_VER_BELOW(msc,gcc,clang,intel) _CY_COMPILER_MSC <  msc
#elif __GNUC__
# define _CY_COMPILER_GCC (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
# define _CY_COMPILER_VER_MEETS(msc,gcc,clang,intel) _CY_COMPILER_GCC >= gcc
# define _CY_COMPILER_VER_BELOW(msc,gcc,clang,intel) _CY_COMPILER_GCC <  gcc
#else
# define _CY_COMPILER_UNKNOWN
# define _CY_COMPILER_VER_MEETS(msc,gcc,clang,intel) false
# define _CY_COMPILER_VER_BELOW(msc,gcc,clang,intel) false
#endif

// constexpr
#ifndef __cpp_constexpr
# if _CY_COMPILER_VER_MEETS(1900,40600,30100,1310)
#  define __cpp_constexpr
# else
#  define constexpr
# endif
#endif

// nullptr
#if _CY_COMPILER_VER_BELOW(1600,40600,20900,1210)
class _cy_nullptr_t {
public:
  template<class T> operator T*() const { return 0; }
  template<class C, class T> operator T C::*() const { return 0; }
private:
  void operator & () const {}
};
static _cy_nullptr_t nullptr;
#endif

// template aliases
#define _CY_TEMPLATE_ALIAS_UNPACK(...) __VA_ARGS__
#if _CY_COMPILER_VER_BELOW(1800,40700,30000,1210)
# define _CY_TEMPLATE_ALIAS(template_name,template_equivalent) class template_name : public _CY_TEMPLATE_ALIAS_UNPACK template_equivalent {}
#else
# define _CY_TEMPLATE_ALIAS(template_name,template_equivalent) using template_name = _CY_TEMPLATE_ALIAS_UNPACK template_equivalent
#endif

// std::is_trivially_copyable
#if _CY_COMPILER_VER_MEETS(1700,50000,30400,1300)
# define _cy_std_is_trivially_copyable 1
#endif

// restrict
#if defined(__INTEL_COMPILER)
//# define restrict restrict
#elif defined(__clang__)
# define restrict __restrict__
#elif defined(_MSC_VER)
# define restrict __restrict
#elif __GNUC__
# define restrict __restrict__
#else
# define restrict
#endif

// alignment
#if _CY_COMPILER_VER_BELOW(1900,40800,30000,1500)
# if defined(_MSC_VER)
#  define alignas(alignment_size) __declspec(align(alignment_size))
# else
#  define alignas(alignment_size) __attribute__((aligned(alignment_size)))
# endif
#endif

// final, override
#if _CY_COMPILER_VER_BELOW(1700,40700,20900,1210)
# define override
# if defined(_MSC_VER)
#  define final sealed
# else
#  define final
# endif
#endif

// static_assert
#if _CY_COMPILER_VER_BELOW(1900,60000,20500,1800)
# define static_assert(condition,message) assert(condition && message)
#endif

// unrestricted unions
#ifndef __cpp_unrestricted_unions
# if _CY_COMPILER_VER_MEETS(1900,40600,30000,1400)
#  define __cpp_unrestricted_unions
# endif
#endif

// nodiscard
#if _CY_COMPILER_VER_MEETS(1901,40800,30000,1500)
# define CY_NODISCARD [[nodiscard]]
#else
# define CY_NODISCARD
#endif

// default and deleted class member functions
#if _CY_COMPILER_VER_MEETS(1800,40400,30000,1200)
# define CY_CLASS_FUNCTION_DEFAULT = default;
# define CY_CLASS_FUNCTION_DELETE  = delete;
#else
# define CY_CLASS_FUNCTION_DEFAULT {}
# define CY_CLASS_FUNCTION_DELETE  { static_assert(false,"Calling deleted method."); }
#endif

//////////////////////////////////////////////////////////////////////////
// Auto Vectorization
//////////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
# if _MSC_VER >= 1700
#  define _CY_IVDEP __pragma(loop(ivdep))
# endif
#elif defined __GNUC__
# if _CY_GCC_VER >= 40900
#  define _CY_IVDEP _Pragma("GCC ivdep");
# endif
#elif defined __clang__
# if _CY_CLANG_VER >= 30500
#  define _CY_IVDEP _Pragma("clang loop vectorize(enable) interleave(enable)");
# endif
#else
//# define _CY_IVDEP _Pragma("ivdep");
# define _CY_IVDEP
#endif

#ifndef _CY_IVDEP
# define _CY_IVDEP
#endif

#define _CY_IVDEP_FOR _CY_IVDEP for

//////////////////////////////////////////////////////////////////////////
// Disabling MSVC's non-standard depreciation warnings
//////////////////////////////////////////////////////////////////////////

#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
# define _CY_CRT_SECURE_NO_WARNINGS     __pragma( warning(push) ) __pragma( warning(disable:4996) )
# define _CY_CRT_SECURE_RESUME_WARNINGS __pragma( warning(pop)  )
#else
# define _CY_CRT_SECURE_NO_WARNINGS
# define _CY_CRT_SECURE_RESUME_WARNINGS
#endif

//////////////////////////////////////////////////////////////////////////
// Math functions
//////////////////////////////////////////////////////////////////////////

//!@name Common math function templates

template <typename T> inline T Max      ( T const &v1, T const &v2 ) { return v1 >= v2 ? v1 : v2; }
template <typename T> inline T Min      ( T const &v1, T const &v2 ) { return v1 <= v2 ? v1 : v2; }
template <typename T> inline T Max      ( T const &v1, T const &v2, T const &v3 ) { return (v1>=v2) ? (v1>=v3?v1:v3) : (v2>=v3?v2:v3); }
template <typename T> inline T Min      ( T const &v1, T const &v2, T const &v3 ) { return (v1<=v2) ? (v1<=v3?v1:v3) : (v2<=v3?v2:v3); }
template <typename T> inline T Clamp    ( T const &v, T minVal=T(0), T maxVal=T(1) ) { return Min(maxVal,Max(minVal,v)); }

template <typename T> inline T ACosSafe ( T const &v ) { return (T) std::acos(Clamp(v,T(-1),T(1))); }
template <typename T> inline T ASinSafe ( T const &v ) { return (T) std::asin(Clamp(v,T(-1),T(1))); }
template <typename T> inline T Sqrt     ( T const &v ) { return (T) std::sqrt(v); }
template <typename T> inline T SqrtSafe ( T const &v ) { return (T) std::sqrt(Max(v,T(0))); }

#ifndef CY_NO_EMMINTRIN_H
template<> inline float  Sqrt    <float> ( float  const &v ) { return _mm_cvtss_f32(_mm_sqrt_ss(_mm_set_ps1(v))); }
template<> inline float  SqrtSafe<float> ( float  const &v ) { return _mm_cvtss_f32(_mm_sqrt_ss(_mm_set_ps1(Max(v,0.0f)))); }
template<> inline double Sqrt    <double>( double const &v ) { __m128d t=_mm_set1_pd(v);          return _mm_cvtsd_f64(_mm_sqrt_sd(t,t)); }
template<> inline double SqrtSafe<double>( double const &v ) { __m128d t=_mm_set1_pd(Max(v,0.0)); return _mm_cvtsd_f64(_mm_sqrt_sd(t,t)); }
#endif

template<typename T> inline T Pi  () { return T(3.141592653589793238462643383279502884197169); }

template <typename T> bool IsFinite( T const &v ) { return std::numeric_limits<T>::is_integer || std::isfinite(v); }

//////////////////////////////////////////////////////////////////////////
// Memory Operations
//////////////////////////////////////////////////////////////////////////

template <typename T>
void MemCopy( T * restrict dest, T const * restrict src, size_t count )
{
#ifdef _cy_std_is_trivially_copyable
	if ( std::is_trivially_copyable<T>() ) {
		memcpy( dest, src, (count)*sizeof(T) );
	} else 
#endif
		for ( size_t i=0; i<count; ++i ) dest[i] = src[i];
}

template <typename T, typename S>
void MemConvert( T * restrict dest, S const * restrict src, size_t count )
{
	for ( size_t i=0; i<count; ++i ) dest[i] = reinterpret_cast<T>(src[i]);
}

template <typename T>
void MemClear( T * dest, size_t count )
{
	memset( dest, 0, count*sizeof(T) );
}

template <typename T> inline void SwapBytes( T &v1, T &v2 ) { char t[sizeof(T)]; memcpy(&t,&v1,sizeof(T)); memcpy(&v1,&v2,sizeof(T)); memcpy(&v2,&t,sizeof(T)); }
template <typename T> inline void Swap     ( T &v1, T &v2 ) { if ( std::is_trivially_copyable<T>::value ) { T t=v1; v1=v2; v2=t; } else SwapBytes(v1,v2); }

/////////////////////////////////////////////////////////////////////////////////
// Sorting functions
/////////////////////////////////////////////////////////////////////////////////

template <bool ascending, typename T>
inline void Sort2( T r[2], T const v[2] )
{
	if ( ascending ) {
		r[0] = Min( v[0], v[1] );
		r[1] = Max( v[1], v[0] );
	} else {
		r[0] = Max( v[0], v[1] );
		r[1] = Min( v[1], v[0] );
	}
}

template <bool ascending, typename T>
void Sort3( T r[3], T const v[3] )
{
	T n01   = Min( v[1], v[0] );
	T x01   = Max( v[0], v[1] );
	T n2x01 = Min( v[2], x01  );
	T r2    = Max( x01,  v[2] );
	T r0    = Min( n2x01, n01 );
	T r1    = Max( n01, n2x01 );
	if ( ascending ) { r[0]=r0; r[1]=r1; r[2]=r2;  }
	else             { r[0]=r2; r[1]=r1; r[2]=r0; }
}

template <bool ascending, typename T>
inline void Sort4( T r[4], T const v[4] )
{
	T n01 = Min( v[1], v[0] );
	T x01 = Max( v[0], v[1] );
	T n23 = Min( v[2], v[3] );
	T x23 = Max( v[3], v[2] );
	T r0  = Min( n01, n23 );
	T x02 = Max( n23, n01 );
	T n13 = Min( x01, x23 );
	T r3  = Max( x23, x01 );
	T r1  = Min( x02, n13 );
	T r2  = Max( n13, x02 );
	if ( ascending ) { r[0]=r0; r[1]=r1; r[2]=r2; r[3]=r3; }
	else             { r[0]=r3; r[1]=r2; r[2]=r1; r[3]=r0; }
}

//////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

#endif

