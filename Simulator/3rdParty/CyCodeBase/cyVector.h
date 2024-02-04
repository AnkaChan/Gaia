// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyVector.h 
//! \author Cem Yuksel
//! 
//! \brief  2D, 3D, 4D, and ND vector classes.
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

#ifndef _CY_VECTOR_H_INCLUDED_
#define _CY_VECTOR_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyCore.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

// Forward declarations
//!	\cond HIDDEN_SYMBOLS
template <typename T> class Vec2;
template <typename T> class Vec3;
template <typename T> class Vec4;
//! \endcond

//-------------------------------------------------------------------------------

//! A general class for N-dimensional vectors.

template <typename T, int N>
class Vec
{
	CY_NODISCARD friend Vec operator - ( T v, Vec const &p ) { return Vec<T,N>(v)-p; }	//!< Subtraction from a constant
	CY_NODISCARD friend Vec operator + ( T v, Vec const &p ) { return p+v; }			//!< Addition with a constant
	CY_NODISCARD friend Vec operator * ( T v, Vec const &p ) { return p*v; }			//!< Multiplication with a constant

public:

	//!@name Components of the vector
	T elem[N];

	//!@name Constructors
	Vec() CY_CLASS_FUNCTION_DEFAULT
	explicit Vec( T const * restrict p ) { MemCopy(elem,p,N); }
	explicit Vec( T v )                  { for ( int i=0; i<N; ++i ) elem[i]=v; }
	template <typename S> explicit Vec( Vec<S,N> const &p ) { MemConvert(elem,p.elem,N); }
	template <int      M> explicit Vec( Vec<T,M> const &p )
	{
		if ( N <= M ) { MemCopy(elem,p.elem,N); }
		else          { MemCopy(elem,p.elem,M); MemClear(elem,N-M); }
	}
	template <typename S, int M> explicit Vec( Vec<S,M> const &p )
	{
		if ( N <= M ) { MemConvert(elem,p.elem,N); }
		else          { MemConvert(elem,p.elem,M); MemClear(elem,N-M); }
	}
	explicit Vec( Vec2<T> const &p );
	explicit Vec( Vec3<T> const &p );
	explicit Vec( Vec4<T> const &p );
	template <typename S> explicit Vec( Vec2<S> const &p );
	template <typename S> explicit Vec( Vec3<S> const &p );
	template <typename S> explicit Vec( Vec4<S> const &p );
	template <typename P> explicit Vec( P       const &p ) { for ( int i=0; i<N; ++i ) elem[i]=(T)p[i]; }

	//!@name Set & Get value methods
	void Zero()                      { MemClear(elem,N); }						//!< Sets the coordinates as zero
	void Get( T * restrict p ) const { MemCopy(p,elem,N); }						//!< Puts the coordinate values into the array
	void Set( T const * restrict p ) { MemCopy(elem,p,N); }						//!< Sets the coordinates using the values in the given array
	void Set( T v )                  { for ( int i=0; i<N; ++i ) elem[i] = v; }	//!< Sets all coordinates using the given value
	template <int M> void CopyData( T * restrict p ) { if ( M <= N ) { MemCopy(p,elem,M); } else { MemCopy(p,elem,N); MemClear(p+N,M-N); }	}
	template <typename S, int M> void ConvertData( S * restrict p ) { if ( M <= N ) { MemConvert(p,elem,M); } else { MemConvert(p,elem,N); MemClear(p+N,M-N); }	}
	void Normalize()                 { *this /= Length(); }						//!< Normalizes the vector, such that its length becomes 1.

	//!@name General methods
	CY_NODISCARD Vec  GetNormalized() const { return *this / Length(); }				//!< Returns a normalized copy of the vector.
	CY_NODISCARD T    LengthSquared() const { Vec p=operator*(*this); return p.Sum(); }	//!< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	CY_NODISCARD T    Length       () const { return cy::Sqrt(LengthSquared()); }		//!< Returns the length of the vector.
	CY_NODISCARD T    Sum          () const { T v=elem[0]; for ( int i=1; i<N; ++i ) v+=elem[i]; return v; }		//!< Returns the sum of its components
	CY_NODISCARD bool IsZero       () const { for ( int i=0; i<N; ++i ) if ( elem[i] != T(0) ) return false; return true; }	//!< Returns true if all components are exactly zero
	CY_NODISCARD T    Min          () const { T m = elem[0]; for ( int i=1; i<N; ++i ) if ( m > elem[i] ) m = elem[i]; return m; }							//!< Returns the minimum component of the vector.
	CY_NODISCARD T    Max          () const { T m = elem[0]; for ( int i=1; i<N; ++i ) if ( m < elem[i] ) m = elem[i]; return m; }							//!< Returns the maximum component of the vector.
	CY_NODISCARD int  MinComp      () const { T m = elem[0]; int ix=0; for ( int i=1; i<N; ++i ) if ( m > elem[i] ) { m = elem[i]; ix = i; } return ix; }	//!< Returns the index of the minimum component of the vector.
	CY_NODISCARD int  MaxComp      () const { T m = elem[0]; int ix=0; for ( int i=1; i<N; ++i ) if ( m < elem[i] ) { m = elem[i]; ix = i; } return ix; }	//!< Returns the index of the maximum component of the vector.
	CY_NODISCARD bool IsFinite     () const { for ( int i=0; i<N; ++i ) if ( ! cy::IsFinite(elem[i]) ) return false; return true; }							//!< Returns true if all components are finite real numbers.
	CY_NODISCARD bool IsUnit       () const { return std::abs(LengthSquared()-T(1)) < T(0.001); }															//!< Returns true if the length of the vector is close to 1.
	CY_NODISCARD Vec  Sqrt         () const { Vec v; for ( int i=0; i<N; ++i ) v.elem[i] = cy::Sqrt(elem[i]); return v; }									//!< Returns the square root of the vector.
	CY_NODISCARD Vec  Abs          () const { Vec v; for ( int i=0; i<N; ++i ) v.elem[i] = std::abs(elem[i]); return v; }									//!< Returns a vector containing the absolute values of all components.

	//!@name Limit methods
	void Clamp   ( T minLimit, T maxLimit ) { ClampMin(minLimit); ClampMax(maxLimit); }		//!< Ensures that all components of the vector are within the given limits.
	void ClampMin( T v ) { for ( int i=0; i<N; ++i ) elem[i] = (elem[i]<v) ? v : elem[i]; }	//!< Ensures that all components of the vector are greater than or equal to the given limit.
	void ClampMax( T v ) { for ( int i=0; i<N; ++i ) elem[i] = (elem[i]>v) ? v : elem[i]; }	//!< Ensures that all components of the vector are smaller than or equal to the given limit.
	void SetAbs  ()      { for ( int i=0; i<N; i++ ) elem[i] = std::abs(elem[i]); }			//!< Converts all negative components to positive values

	//!@name Unary operators
	CY_NODISCARD Vec  operator - () const { Vec r; for ( int i=0; i<N; ++i ) r.elem[i]=-elem[i]; return r; } 

	//!@name Binary operators
	CY_NODISCARD Vec  operator + ( Vec const &p ) const { Vec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] + p.elem[i]; return r; }
	CY_NODISCARD Vec  operator - ( Vec const &p ) const { Vec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] - p.elem[i]; return r; }
	CY_NODISCARD Vec  operator * ( Vec const &p ) const { Vec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] * p.elem[i]; return r; }
	CY_NODISCARD Vec  operator / ( Vec const &p ) const { Vec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] / p.elem[i]; return r; }
	CY_NODISCARD Vec  operator + ( T   const  v ) const { Vec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] + v; return r; }
	CY_NODISCARD Vec  operator - ( T   const  v ) const { Vec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] - v; return r; }
	CY_NODISCARD Vec  operator * ( T   const  v ) const { Vec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] * v; return r; }
	CY_NODISCARD Vec  operator / ( T   const  v ) const { Vec r; for ( int i=0; i<N; ++i ) r.elem[i] = elem[i] / v; return r; }

	//!@name Assignment operators
	Vec const& operator += ( Vec const &p ) { for ( int i=0; i<N; ++i ) elem[i] += p.elem[i]; return *this; }
	Vec const& operator -= ( Vec const &p ) { for ( int i=0; i<N; ++i ) elem[i] -= p.elem[i]; return *this; }
	Vec const& operator *= ( Vec const &p ) { for ( int i=0; i<N; ++i ) elem[i] *= p.elem[i]; return *this; }
	Vec const& operator /= ( Vec const &p ) { for ( int i=0; i<N; ++i ) elem[i] /= p.elem[i]; return *this; }
	Vec const& operator += ( T   const  v ) { for ( int i=0; i<N; ++i ) elem[i] += v; return *this; }
	Vec const& operator -= ( T   const  v ) { for ( int i=0; i<N; ++i ) elem[i] -= v; return *this; }
	Vec const& operator *= ( T   const  v ) { for ( int i=0; i<N; ++i ) elem[i] *= v; return *this; }
	Vec const& operator /= ( T   const  v ) { for ( int i=0; i<N; ++i ) elem[i] /= v; return *this; }

	//!@name Test operators
	CY_NODISCARD bool operator == ( Vec const& p ) const { for ( int i=0; i<N; ++i ) if ( elem[i] != p.elem[i] ) return false; return true; }
	CY_NODISCARD bool operator != ( Vec const& p ) const { for ( int i=0; i<N; ++i ) if ( elem[i] != p.elem[i] ) return true; return false; }

	//!@name Access operators
	CY_NODISCARD T&       operator [] ( int i )       { return Element(i); }
	CY_NODISCARD T        operator [] ( int i ) const { return Element(i); }
	CY_NODISCARD T&       Element     ( int i )       { assert(i>=0 && i<N); return elem[i]; }
	CY_NODISCARD T const& Element     ( int i ) const { assert(i>=0 && i<N); return elem[i]; }
	CY_NODISCARD T*       Elements    ()              { return elem; }
	CY_NODISCARD T const* Elements    ()        const { return elem; }

	//!@name Dot product
	CY_NODISCARD T   Dot        ( Vec const &p ) const { Vec r=operator*(p); return r.Sum(); }	//!< Dot product
	CY_NODISCARD T   operator % ( Vec const &p ) const { return Dot(p); }						//!< Dot product operator
};

//-------------------------------------------------------------------------------

//! 2D vector class

template <typename T>
class Vec2
{
	CY_NODISCARD friend Vec2 operator - ( T v, Vec2 const &p ) { return Vec2<T>(v)-p; }	//!< Subtraction from a constant
	CY_NODISCARD friend Vec2 operator + ( T v, Vec2 const &p ) { return p+v; }			//!< Addition with a constant
	CY_NODISCARD friend Vec2 operator * ( T v, Vec2 const &p ) { return p*v; }			//!< Multiplication with a constant

public:

	//!@name Components of the vector
	union {
		struct { T x, y; };
		T elem[2];	//!< Array-type access to the vector elements x and y
	};

	//!@name Constructors
	Vec2() CY_CLASS_FUNCTION_DEFAULT
	Vec2( T _x, T _y )   : x( _x), y( _y) {}
	explicit Vec2( T v ) : x(v  ), y(v  ) {}
	explicit Vec2( Vec3<T> const &p );
	explicit Vec2( Vec4<T> const &p );
	explicit Vec2( T const * restrict v ) { Set( v ); }
	template <typename S> explicit Vec2( Vec2<S> const &p ) : x(T(p.x)), y(T(p.y)) {}
	template <typename S> explicit Vec2( Vec3<S> const &p );
	template <typename S> explicit Vec2( Vec4<S> const &p );
	template <int N            > explicit Vec2( Vec<T,N> const &p ) { p.template CopyData<2>(elem); }
	template <int N, typename S> explicit Vec2( Vec<S,N> const &p ) { p.template ConvertData<T,2>(elem); }

	//!@name Set & Get value methods
	void Zero()                      { MemClear(elem,2); }				//!< Sets the coordinates as zero.
	void Get( T * restrict p ) const { ((Vec2*)p)->operator=(*this); }	//!< Puts the coordinate values into the array.
	void Set( T const * restrict p ) { operator=(*((Vec2*)p)); }		//!< Sets the coordinates using the values in the given array.
	void Set( T v )                  { x=v; y=v; }						//!< Sets all coordinates using the given value
	void Set( T _x, T _y )           { x=_x; y=_y; }					//!< Sets the coordinates using the given values
	void Normalize()                 { *this /= Length(); }				//!< Normalizes the vector, such that its length becomes 1.

	//!@name General methods
	CY_NODISCARD Vec2 GetNormalized   () const { return *this / Length(); }								//!< Returns a normalized copy of the vector.
	CY_NODISCARD T    LengthSquared   () const { Vec2 p=operator*(*this); return p.Sum(); }				//!< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	CY_NODISCARD T    Length          () const { return cy::Sqrt(LengthSquared()); }					//!< Returns the length of the vector.
	CY_NODISCARD T    Sum             () const { return x+y; }											//!< Returns the sum of its components
	CY_NODISCARD bool IsZero          () const { return x==T(0) && y==T(0); }							//!< Returns true if all components are exactly zero
	CY_NODISCARD T    Min             () const { return elem[MinComp()]; }								//!< Returns the minimum component of the vector.
	CY_NODISCARD T    Max             () const { return elem[MaxComp()]; }								//!< Returns the maximum component of the vector.
	CY_NODISCARD int  MinComp         () const { return x>y; }											//!< Returns the index of the minimum component of the vector.
	CY_NODISCARD int  MaxComp         () const { return x<y; }											//!< Returns the index of the maximum component of the vector.
	CY_NODISCARD bool IsFinite        () const { return cy::IsFinite(x) && cy::IsFinite(y); }			//!< Returns true if all components are finite real numbers.
	CY_NODISCARD bool IsUnit          () const { return std::abs(LengthSquared()-T(1)) < T(0.001); }	//!< Returns true if the length of the vector is close to 1.
	CY_NODISCARD Vec2 Sqrt            () const { return Vec2(cy::Sqrt(x),cy::Sqrt(y)); }				//!< Returns the square root of the vector.
	CY_NODISCARD Vec2 Abs             () const { return Vec2(std::abs(x),std::abs(y)); }				//!< Returns a vector containing the absolute values of all components.
	CY_NODISCARD Vec2 SortAsc         () const { Vec2 v; Sort2<true >( v.elem, elem ); return v; }		//!< Returns a vector with components sorted in ascending order.
	CY_NODISCARD Vec2 SortDesc        () const { Vec2 v; Sort2<false>( v.elem, elem ); return v; }		//!< Returns a vector with components sorted in descending order.
	CY_NODISCARD Vec2 GetPerpendicular() const { return Vec2(-y,x); }									//!< Returns a perpendicular vector (rotated by 90 degrees in counter clockwise direction).

	//!@name Limit methods
	void Clamp   ( T minLimit, T maxLimit ) { ClampMin(minLimit); ClampMax(maxLimit); }	//!< Ensures that all components of the vector are within the given limits.
	void ClampMin( T v ) { x=(x<v)?v:x; y=(y<v)?v:y; }									//!< Ensures that all components of the vector are greater than or equal to the given limit.
	void ClampMax( T v ) { x=(x>v)?v:x; y=(y>v)?v:y; }									//!< Ensures that all components of the vector are smaller than or equal to the given limit.
	void SetAbs  ()      { x=std::abs(x); y=std::abs(y); }								//!< Converts all negative components to positive values

	//!@name Unary operators
	CY_NODISCARD Vec2 operator - () const { Vec2 r; r.x=-x; r.y=-y; return r; } 

	//!@name Binary operators
	CY_NODISCARD Vec2 operator + ( Vec2 const &p ) const { Vec2 r; r.x=x+p.x; r.y=y+p.y; return r; }
	CY_NODISCARD Vec2 operator - ( Vec2 const &p ) const { Vec2 r; r.x=x-p.x; r.y=y-p.y; return r; }
	CY_NODISCARD Vec2 operator * ( Vec2 const &p ) const { Vec2 r; r.x=x*p.x; r.y=y*p.y; return r; }
	CY_NODISCARD Vec2 operator / ( Vec2 const &p ) const { Vec2 r; r.x=x/p.x; r.y=y/p.y; return r; }
	CY_NODISCARD Vec2 operator + ( T    const  v ) const { Vec2 r; r.x=x+v;   r.y=y+v;   return r; }
	CY_NODISCARD Vec2 operator - ( T    const  v ) const { Vec2 r; r.x=x-v;   r.y=y-v;   return r; }
	CY_NODISCARD Vec2 operator * ( T    const  v ) const { Vec2 r; r.x=x*v;   r.y=y*v;   return r; }
	CY_NODISCARD Vec2 operator / ( T    const  v ) const { Vec2 r; r.x=x/v;   r.y=y/v;   return r; }

	//!@name Assignment operators
	Vec2 const& operator += ( Vec2 const &p ) { x+=p.x; y+=p.y; return *this; }
	Vec2 const& operator -= ( Vec2 const &p ) { x-=p.x; y-=p.y; return *this; }
	Vec2 const& operator *= ( Vec2 const &p ) { x*=p.x; y*=p.y; return *this; }
	Vec2 const& operator /= ( Vec2 const &p ) { x/=p.x; y/=p.y; return *this; }
	Vec2 const& operator += ( T    const  v ) { x+=v;   y+=v;   return *this; }
	Vec2 const& operator -= ( T    const  v ) { x-=v;   y-=v;   return *this; }
	Vec2 const& operator *= ( T    const  v ) { x*=v;   y*=v;   return *this; }
	Vec2 const& operator /= ( T    const  v ) { x/=v;   y/=v;   return *this; }

	//!@name Test operators
	CY_NODISCARD bool operator == ( Vec2 const& p ) const { return x==p.x && y==p.y; }
	CY_NODISCARD bool operator != ( Vec2 const& p ) const { return x!=p.x && y!=p.y; }

	//!@name Access operators
	CY_NODISCARD T&       operator [] ( int i )       { return Element(i); }
	CY_NODISCARD T const& operator [] ( int i ) const { return Element(i); }
	CY_NODISCARD T&       Element     ( int i )       { assert(i>=0 && i<2); return elem[i]; }
	CY_NODISCARD T const& Element     ( int i ) const { assert(i>=0 && i<2); return elem[i]; }
	CY_NODISCARD T*       Elements    ()              { return elem; }
	CY_NODISCARD T const* Elements    ()        const { return elem; }

	//!@name Cross product and dot product
	CY_NODISCARD T Cross      ( Vec2 const &p ) const { Vec2 r(-y,x); return r.Dot(p); }	//!< Cross product
	CY_NODISCARD T operator ^ ( Vec2 const &p ) const { return Cross(p); }					//!< Cross product operator
	CY_NODISCARD T Dot        ( Vec2 const &p ) const { return x*p.x + y*p.y; }				//!< Dot product
	CY_NODISCARD T operator % ( Vec2 const &p ) const { return Dot(p); }					//!< Dot product operator

	//!@name Swizzling Methods
	CY_NODISCARD Vec2<T> XX() const { return Vec2<T>(x,x); }
	CY_NODISCARD Vec2<T> XY() const { return *this; }
	CY_NODISCARD Vec2<T> YX() const { return Vec2<T>(y,x); }
	CY_NODISCARD Vec2<T> YY() const { return Vec2<T>(y,y); }
};

//-------------------------------------------------------------------------------

//! 3D vector class

template <typename T>
class Vec3
{
	CY_NODISCARD friend Vec3 operator - ( T v, Vec3 const &p ) { return Vec3<T>(v)-p; }	//!< Subtraction from a constant
	CY_NODISCARD friend Vec3 operator + ( T v, Vec3 const &p ) { return p+v; }			//!< Addition with a constant
	CY_NODISCARD friend Vec3 operator * ( T v, Vec3 const &p ) { return p*v; }			//!< Multiplication with a constant

public:

	//!@name Components of the vector
	union {
		struct { T x, y, z; };
		T elem[3];	//!< Array-type access to the vector elements x, y, and z
	};

	//!@name Constructors
	Vec3() CY_CLASS_FUNCTION_DEFAULT
	Vec3( T _x, T _y, T _z )                  : x( _x), y( _y), z( _z) {}
	explicit Vec3( T v )                      : x(v  ), y(v  ), z(v  ) {}
	explicit Vec3( Vec2<T> const &p, T _z=0 ) : x(p.x), y(p.y), z( _z) {}
	explicit Vec3( Vec4<T> const &p );
	explicit Vec3( T const * restrict v ) { Set( v ); }
	template <typename S> explicit Vec3( Vec3<S> const &p )         : x(T(p.x)), y(T(p.y)), z(T(p.z)) {}
	template <typename S> explicit Vec3( Vec2<S> const &p, T _z=0 ) : x(T(p.x)), y(T(p.y)), z(   _z ) {}
	template <typename S> explicit Vec3( Vec4<S> const &p );
	template <int N            > explicit Vec3( Vec<T,N> const &p ) { p.template CopyData<3>(elem); }
	template <int N, typename S> explicit Vec3( Vec<S,N> const &p ) { p.template ConvertData<T,3>(elem); }

	//!@name Set & Get value methods
	void Zero()                            { MemClear(elem,3); }				//!< Sets the coordinates as zero.
	void Get( T       * restrict p ) const { ((Vec3*)p)->operator=(*this); }	//!< Puts the coordinate values into the array.
	void Set( T const * restrict p )       { operator=(*((Vec3*)p)); }			//!< Sets the coordinates using the values in the given array.
	void Set( T v )                        { x=v; y=v; z=v; }					//!< Sets all coordinates using the given value.
	void Set( T _x, T _y, T _z )           { x= _x; y= _y; z=_z; }				//!< Sets the coordinates using the given values.
	void Set( Vec2<T> const &p, T _z )     { x=p.x; y=p.y; z=_z; }				//!< Sets the coordinates using the given values.
	void Normalize()                       { *this /= Length(); }				//!< Normalizes the vector, such that its length becomes 1.

	//!@name General methods
	CY_NODISCARD Vec3 GetNormalized   () const { return *this / Length(); }								//!< Returns a normalized copy of the vector.
	CY_NODISCARD T    LengthSquared   () const { Vec3 p=operator*(*this); return p.Sum(); }				//!< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	CY_NODISCARD T    Length          () const { return cy::Sqrt(LengthSquared()); }					//!< Returns the length of the vector.
	CY_NODISCARD T    Sum             () const { return x+y+z; }										//!< Returns the sum of its components.
	CY_NODISCARD bool IsZero          () const { return x==T(0) && y==T(0) && z==T(0); }				//!< Returns true if all components are exactly zero.
	CY_NODISCARD T    Min             () const { return elem[MinComp()]; }								//!< Returns the minimum component of the vector.
	CY_NODISCARD T    Max             () const { return elem[MaxComp()]; }								//!< Returns the maximum component of the vector.
	CY_NODISCARD int  MinComp         () const { int yx=y<x; int zx=z<x; int zy=z<y; return (yx|zx)+(zx&zy); }	//!< Returns the index of the minimum component of the vector.
	CY_NODISCARD int  MaxComp         () const { int xy=x<y; int xz=x<z; int yz=y<z; return (xy|xz)+(xz&yz); }	//!< Returns the index of the maximum component of the vector.
	CY_NODISCARD bool IsFinite        () const { return cy::IsFinite(x) && cy::IsFinite(y) && cy::IsFinite(z); }	//!< Returns true if all components are finite real numbers.
	CY_NODISCARD bool IsUnit          () const { return std::abs(LengthSquared()-T(1)) < T(0.001); }	//!< Returns true if the length of the vector is close to 1.
	CY_NODISCARD Vec3 Sqrt            () const { return Vec3(cy::Sqrt(x),cy::Sqrt(y),cy::Sqrt(z)); }	//!< Returns the square root of the vector.
	CY_NODISCARD Vec3 Abs             () const { return Vec3(std::abs(x),std::abs(y),std::abs(z)); }	//!< Returns a vector containing the absolute values of all components.
	CY_NODISCARD Vec3 SortAsc         () const { Vec3 v; Sort3<true >( v.elem, elem ); return v; }		//!< Returns a vector with components sorted in ascending order.
	CY_NODISCARD Vec3 SortDesc        () const { Vec3 v; Sort3<false>( v.elem, elem ); return v; }		//!< Returns a vector with components sorted in descending order.
	CY_NODISCARD Vec3 GetPerpendicular() const { Vec3 v0,v1; GetOrthonormals(v0,v1); return v0; }		//!< Returns a perpendicular vector
	CY_NODISCARD int MaxIndex() const {
		int maxId = x > y ? 0 : 1;
		maxId = elem[maxId] > z ? maxId : 2;
		return maxId;
	}

	void GetOrthonormals ( Vec3 &v0, Vec3 &v1 ) const	//!< Returns two orthogonal vectors to this vector, forming an orthonormal basis
	{
		if ( z >= y ) {
			T const a = T(1)/(1 + z);
			T const b = -x*y*a;
			v0.Set( 1 - x*x*a, b, -x );
			v1.Set( b, 1 - y*y*a, -y );
		} else {
			T const a = T(1)/(1 + y);
			T const b = -x*z*a;
			v0.Set( b, -z, 1 - z*z*a );
			v1.Set( 1 - x*x*a, -x, b );
		}
	}

	//!@name Limit methods
	void Clamp   ( T minLimit, T maxLimit ) { ClampMin(minLimit); ClampMax(maxLimit); }	//!< Ensures that all components of the vector are within the given limits.
	void ClampMin( T v ) { x=(x<v)?v:x; y=(y<v)?v:y; z=(z<v)?v:z; }						//!< Ensures that all components of the vector are greater than or equal to the given limit.
	void ClampMax( T v ) { x=(x>v)?v:x; y=(y>v)?v:y; z=(z>v)?v:z; }						//!< Ensures that all components of the vector are smaller than or equal to the given limit.
	void SetAbs  ()      { x=std::abs(x); y=std::abs(y); z=std::abs(z); }				//!< Converts all negative components to positive values

	//!@name Unary operators
	CY_NODISCARD Vec3 operator - () const { Vec3 r; r.x=-x; r.y=-y; r.z=-z; return r; } 

	//!@name Binary operators
	CY_NODISCARD Vec3 operator + ( Vec3 const &p ) const { Vec3 r; r.x=x+p.x; r.y=y+p.y; r.z=z+p.z; return r; }
	CY_NODISCARD Vec3 operator - ( Vec3 const &p ) const { Vec3 r; r.x=x-p.x; r.y=y-p.y; r.z=z-p.z; return r; }
	CY_NODISCARD Vec3 operator * ( Vec3 const &p ) const { Vec3 r; r.x=x*p.x; r.y=y*p.y; r.z=z*p.z; return r; }
	CY_NODISCARD Vec3 operator / ( Vec3 const &p ) const { Vec3 r; r.x=x/p.x; r.y=y/p.y; r.z=z/p.z; return r; }
	CY_NODISCARD Vec3 operator + ( T    const  v ) const { Vec3 r; r.x=x+v;   r.y=y+v;   r.z=z+v;   return r; }
	CY_NODISCARD Vec3 operator - ( T    const  v ) const { Vec3 r; r.x=x-v;   r.y=y-v;   r.z=z-v;   return r; }
	CY_NODISCARD Vec3 operator * ( T    const  v ) const { Vec3 r; r.x=x*v;   r.y=y*v;   r.z=z*v;   return r; }
	CY_NODISCARD Vec3 operator / ( T    const  v ) const { Vec3 r; r.x=x/v;   r.y=y/v;   r.z=z/v;   return r; }

	//!@name Assignment operators
	Vec3 const& operator += ( Vec3 const &p ) { x+=p.x; y+=p.y; z+=p.z; return *this; }
	Vec3 const& operator -= ( Vec3 const &p ) { x-=p.x; y-=p.y; z-=p.z; return *this; }
	Vec3 const& operator *= ( Vec3 const &p ) { x*=p.x; y*=p.y; z*=p.z; return *this; }
	Vec3 const& operator /= ( Vec3 const &p ) { x/=p.x; y/=p.y; z/=p.z; return *this; }
	Vec3 const& operator += ( T    const  v ) { x+=v;   y+=v;   z+=v;   return *this; }
	Vec3 const& operator -= ( T    const  v ) { x-=v;   y-=v;   z-=v;   return *this; }
	Vec3 const& operator *= ( T    const  v ) { x*=v;   y*=v;   z*=v;   return *this; }
	Vec3 const& operator /= ( T    const  v ) { x/=v;   y/=v;   z/=v;   return *this; }

	//!@name Test operators
	CY_NODISCARD bool operator == ( Vec3 const& p ) const { return x==p.x && y==p.y && z==p.z; }
	CY_NODISCARD bool operator != ( Vec3 const& p ) const { return x!=p.x && y!=p.y && z!=p.z; }

	//!@name Access operators
	CY_NODISCARD T&       operator [] ( int i )       { return Element(i); }
	CY_NODISCARD T const& operator [] ( int i ) const { return Element(i); }
	CY_NODISCARD T&       Element     ( int i )       { assert(i>=0 && i<3); return elem[i]; }
	CY_NODISCARD T const& Element     ( int i ) const { assert(i>=0 && i<3); return elem[i]; }
	CY_NODISCARD T*       Elements    ()              { return elem; }
	CY_NODISCARD T const* Elements    ()        const { return elem; }

	//!@name Cross product and dot product
	CY_NODISCARD Vec3 Cross      ( Vec3 const &p ) const { return Vec3(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x); }	//!< Cross product
	CY_NODISCARD Vec3 operator ^ ( Vec3 const &p ) const { return Cross(p); }										//!< Cross product
	CY_NODISCARD T    Dot        ( Vec3 const &p ) const { return x*p.x + y*p.y + z*p.z; }							//!< Dot product
	CY_NODISCARD T    operator % ( Vec3 const &p ) const { return Dot(p); }											//!< Dot product

	//!@name Swizzling Methods
	CY_NODISCARD Vec2<T> XX() const { return Vec2<T>(x,x); }
	CY_NODISCARD Vec2<T> XY() const { return Vec2<T>(*this); }
	CY_NODISCARD Vec2<T> XZ() const { return Vec2<T>(x,z); }
	CY_NODISCARD Vec2<T> YX() const { return Vec2<T>(y,x); }
	CY_NODISCARD Vec2<T> YY() const { return Vec2<T>(y,y); }
	CY_NODISCARD Vec2<T> YZ() const { return Vec2<T>(y,z); }
	CY_NODISCARD Vec2<T> ZX() const { return Vec2<T>(z,x); }
	CY_NODISCARD Vec2<T> ZY() const { return Vec2<T>(z,y); }
	CY_NODISCARD Vec2<T> ZZ() const { return Vec2<T>(z,z); }

	CY_NODISCARD Vec3<T> XXX() const { return Vec3<T>(x,x,x); }
	CY_NODISCARD Vec3<T> XXY() const { return Vec3<T>(x,x,y); }
	CY_NODISCARD Vec3<T> XXZ() const { return Vec3<T>(x,x,z); }
	CY_NODISCARD Vec3<T> XYX() const { return Vec3<T>(x,y,x); }
	CY_NODISCARD Vec3<T> XYY() const { return Vec3<T>(x,y,y); }
	CY_NODISCARD Vec3<T> XYZ() const { return *this; }
	CY_NODISCARD Vec3<T> XZX() const { return Vec3<T>(x,z,x); }
	CY_NODISCARD Vec3<T> XZY() const { return Vec3<T>(x,z,y); }
	CY_NODISCARD Vec3<T> XZZ() const { return Vec3<T>(x,z,z); }

	CY_NODISCARD Vec3<T> YXX() const { return Vec3<T>(y,x,x); }
	CY_NODISCARD Vec3<T> YXY() const { return Vec3<T>(y,x,y); }
	CY_NODISCARD Vec3<T> YXZ() const { return Vec3<T>(y,x,z); }
	CY_NODISCARD Vec3<T> YYX() const { return Vec3<T>(y,y,x); }
	CY_NODISCARD Vec3<T> YYY() const { return Vec3<T>(y,y,y); }
	CY_NODISCARD Vec3<T> YYZ() const { return Vec3<T>(y,y,z); }
	CY_NODISCARD Vec3<T> YZX() const { return Vec3<T>(y,z,x); }
	CY_NODISCARD Vec3<T> YZY() const { return Vec3<T>(y,z,y); }
	CY_NODISCARD Vec3<T> YZZ() const { return Vec3<T>(y,z,z); }

	CY_NODISCARD Vec3<T> ZXX() const { return Vec3<T>(z,x,x); }
	CY_NODISCARD Vec3<T> ZXY() const { return Vec3<T>(z,x,y); }
	CY_NODISCARD Vec3<T> ZXZ() const { return Vec3<T>(z,x,z); }
	CY_NODISCARD Vec3<T> ZYX() const { return Vec3<T>(z,y,x); }
	CY_NODISCARD Vec3<T> ZYY() const { return Vec3<T>(z,y,y); }
	CY_NODISCARD Vec3<T> ZYZ() const { return Vec3<T>(z,y,z); }
	CY_NODISCARD Vec3<T> ZZX() const { return Vec3<T>(z,z,x); }
	CY_NODISCARD Vec3<T> ZZY() const { return Vec3<T>(z,z,y); }
	CY_NODISCARD Vec3<T> ZZZ() const { return Vec3<T>(z,z,z); }
};

//-------------------------------------------------------------------------------

//! 4D vector class

template <typename T>
class Vec4
{
	CY_NODISCARD friend Vec4 operator - ( T v, Vec4 const &p ) { return Vec4<T>(v)-p; }	//!< Subtraction from a constant
	CY_NODISCARD friend Vec4 operator + ( T v, Vec4 const &p ) { return p+v; }			//!< Addition with a constant
	CY_NODISCARD friend Vec4 operator * ( T v, Vec4 const &p ) { return p*v; }			//!< Multiplication with a constant

public:

	//!@name Components of the vector
	union {
		struct { T x, y, z, w; };
		T elem[4];	//!< Array-type access to the vector elements x, y, z, and w
	};

	//!@name Constructors
	Vec4() CY_CLASS_FUNCTION_DEFAULT
	Vec4( T _x, T _y, T _z, T _w )                    : x( _x), y( _y), z( _z), w( _w) {}
	explicit Vec4( T v )                              : x(v  ), y(v  ), z(v  ), w(v  ) {}
	explicit Vec4( Vec2<T> const &p, T _z=0, T _w=1 ) : x(p.x), y(p.y), z( _z), w( _w) {}
	explicit Vec4( Vec3<T> const &p,         T _w=1 ) : x(p.x), y(p.y), z(p.z), w( _w) {}
	explicit Vec4( T const * restrict v ) { Set( v ); }
	template <typename S> explicit Vec4( Vec2<S> const &p, T _z=0, T _w=1 ) : x(T(p.x)), y(T(p.y)), z(   _z ), w(   _w ) {}
	template <typename S> explicit Vec4( Vec3<S> const &p,         T _w=1 ) : x(T(p.x)), y(T(p.y)), z(T(p.z)), w(   _w ) {}
	template <typename S> explicit Vec4( Vec4<S> const &p )                 : x(T(p.x)), y(T(p.y)), z(T(p.z)), w(T(p.w)) {}
	template <int N            > explicit Vec4( Vec<T,N> const &p ) { p.template CopyData<4>(elem); }
	template <int N, typename S> explicit Vec4( Vec<S,N> const &p ) { p.template ConvertData<T,4>(elem); }

	//!@name Set & Get value methods
	void Zero()                                { MemClear(elem,4); }				//!< Sets the coordinates as zero
	void Get( T       * restrict p ) const     { ((Vec4*)p)->operator=(*this); }	//!< Puts the coordinate values into the array
	void Set( T const * restrict p )           { operator=(*((Vec4*)p)); }			//!< Sets the coordinates using the values in the given array
	void Set( T v )                            { x=v; y=v; z=v; w=v; }				//!< Sets all coordinates using the given value
	void Set( T _x, T _y, T _z, T _w=1 )       { x= _x; y= _y; z= _z; w=_w; }		//!< Sets the coordinates using the given values
	void Set( Vec2<T> const &p, T _z, T _w=1 ) { x=p.x; y=p.y; z= _z; w=_w; }		//!< Sets the coordinates using the given values
	void Set( Vec3<T> const &p,       T _w=1 ) { x=p.x; y=p.y; z=p.z; w=_w; }		//!< Sets the coordinates using the given values
	void Normalize()                           { *this /= Length(); }				//!< Normalizes the vector, such that its length becomes 1.

	//!@name General methods
	CY_NODISCARD Vec4 GetNormalized() const { return *this / Length(); }							//!< Returns a normalized copy of the vector.
	CY_NODISCARD T    LengthSquared() const { Vec4 p=operator*(*this); return p.Sum(); }			//!< Returns the square of the length. Effectively, this is the dot product of the vector with itself.
	CY_NODISCARD T    Length       () const { return cy::Sqrt(LengthSquared()); }					//!< Returns the length of the vector.
	CY_NODISCARD T    Sum          () const { return x+y+z+w; }										//!< Returns the sum of its components
	CY_NODISCARD bool IsZero       () const { return x==T(0) && y==T(0) && z==T(0) && w==T(0); }	//!< Returns true if all components are exactly zero
	CY_NODISCARD T    Min          () const { return elem[MaxComp()]; }								//!< Returns the minimum component of the vector.
	CY_NODISCARD T    Max          () const { return elem[MaxComp()]; }								//!< Returns the maximum component of the vector.
	CY_NODISCARD int  MinComp      () const { int xy=x>y; int zw=(z>w)+2; return elem[xy]<elem[zw]?xy:zw; }	//!< Returns the index of the minimum component of the vector.
	CY_NODISCARD int  MaxComp      () const { int xy=x<y; int zw=(z<w)+2; return elem[xy]>elem[zw]?xy:zw; }	//!< Returns the index of the maximum component of the vector.
	CY_NODISCARD bool IsFinite     () const { return cy::IsFinite(x) && cy::IsFinite(y) && cy::IsFinite(z) && cy::IsFinite(w); }	//!< Returns true if all components are finite real numbers.
	CY_NODISCARD bool IsUnit       () const { return std::abs(LengthSquared()-T(1)) < T(0.001); }				//!< Returns true if the length of the vector is close to 1.
	CY_NODISCARD Vec4 Sqrt         () const { return Vec4(cy::Sqrt(x),cy::Sqrt(y),cy::Sqrt(z),cy::Sqrt(w)); }	//!< Returns the square root of the vector.
	CY_NODISCARD Vec4 Abs          () const { return Vec4(std::abs(x),std::abs(y),std::abs(z),std::abs(w)); }	//!< Returns a vector containing the absolute values of all components.
	CY_NODISCARD Vec4 SortAsc      () const { Vec4 v; Sort4<true >( v.elem, elem ); return v; }		//!< Returns a vector with components sorted in ascending order.
	CY_NODISCARD Vec4 SortDesc     () const { Vec4 v; Sort4<false>( v.elem, elem ); return v; }		//!< Returns a vector with components sorted in descending order.

	//!@name Limit methods
	void Clamp   ( T minLimit, T maxLimit ) { ClampMin(minLimit); ClampMax(maxLimit); }		//!< Ensures that all components of the vector are within the given limits.
	void ClampMin( T v ) { x=(x<v)?v:x; y=(y<v)?v:y; z=(z<v)?v:z; w=(w<v)?v:w; }			//!< Ensures that all components of the vector are greater than or equal to the given limit.
	void ClampMax( T v ) { x=(x>v)?v:x; y=(y>v)?v:y; z=(z>v)?v:z; w=(w>v)?v:w; }			//!< Ensures that all components of the vector are smaller than or equal to the given limit.
	void SetAbs  ()      { x=std::abs(x); y=std::abs(y); z=std::abs(z); w=std::abs(w); }	//!< Converts all negative components to positive values

	//!@name Unary operators
	Vec4 operator - () const { Vec4 r; r.x=-x; r.y=-y; r.z=-z; r.w=-w; return r; } 

	//!@name Binary operators
	CY_NODISCARD Vec4 operator + ( Vec4 const &p ) const { Vec4 r; r.x=x+p.x; r.y=y+p.y; r.z=z+p.z; r.w=w+p.w; return r; }
	CY_NODISCARD Vec4 operator - ( Vec4 const &p ) const { Vec4 r; r.x=x-p.x; r.y=y-p.y; r.z=z-p.z; r.w=w-p.w; return r; }
	CY_NODISCARD Vec4 operator * ( Vec4 const &p ) const { Vec4 r; r.x=x*p.x; r.y=y*p.y; r.z=z*p.z; r.w=w*p.w; return r; }
	CY_NODISCARD Vec4 operator / ( Vec4 const &p ) const { Vec4 r; r.x=x/p.x; r.y=y/p.y; r.z=z/p.z; r.w=w/p.w; return r; }
	CY_NODISCARD Vec4 operator + ( T    const  v ) const { Vec4 r; r.x=x+v;   r.y=y+v;   r.z=z+v;   r.w=w+v;   return r; }
	CY_NODISCARD Vec4 operator - ( T    const  v ) const { Vec4 r; r.x=x-v;   r.y=y-v;   r.z=z-v;   r.w=w-v;   return r; }
	CY_NODISCARD Vec4 operator * ( T    const  v ) const { Vec4 r; r.x=x*v;   r.y=y*v;   r.z=z*v;   r.w=w*v;   return r; }
	CY_NODISCARD Vec4 operator / ( T    const  v ) const { Vec4 r; r.x=x/v;   r.y=y/v;   r.z=z/v;   r.w=w/v;   return r; }

	//!@name Assignment operators
	Vec4 const& operator += ( Vec4 const &p ) { x+=p.x; y+=p.y; z+=p.z; w+=p.w; return *this; }
	Vec4 const& operator -= ( Vec4 const &p ) { x-=p.x; y-=p.y; z-=p.z; w-=p.w; return *this; }
	Vec4 const& operator *= ( Vec4 const &p ) { x*=p.x; y*=p.y; z*=p.z; w*=p.w; return *this; }
	Vec4 const& operator /= ( Vec4 const &p ) { x/=p.x; y/=p.y; z/=p.z; w/=p.w; return *this; }
	Vec4 const& operator += ( T    const  v ) { x+=v;   y+=v;   z+=v;   w+=v;   return *this; }
	Vec4 const& operator -= ( T    const  v ) { x-=v;   y-=v;   z-=v;   w-=v;   return *this; }
	Vec4 const& operator *= ( T    const  v ) { x*=v;   y*=v;   z*=v;   w*=v;   return *this; }
	Vec4 const& operator /= ( T    const  v ) { x/=v;   y/=v;   z/=v;   w/=v;   return *this; }

	//!@name Test operators
	CY_NODISCARD bool operator == ( Vec4 const& p ) const { return x==p.x && y==p.y && z==p.z && w==p.w; }
	CY_NODISCARD bool operator != ( Vec4 const& p ) const { return x!=p.x && y!=p.y && z!=p.z && w!=p.w; }

	//!@name Access operators
	CY_NODISCARD T&       operator [] ( int i )       { return Element(i); }
	CY_NODISCARD T const& operator [] ( int i ) const { return Element(i); }
	CY_NODISCARD T&       Element     ( int i )       { assert(i>=0 && i<4); return elem[i]; }
	CY_NODISCARD T const& Element     ( int i ) const { assert(i>=0 && i<4); return elem[i]; }
	CY_NODISCARD T*       Elements    ()              { return elem; }
	CY_NODISCARD T const* Elements    ()        const { return elem; }

	//!@name Dot product
	CY_NODISCARD T Dot        ( Vec4 const &p ) const { return x*p.x + y*p.y + z*p.z + w*p.w; }	//!< Dot product
	CY_NODISCARD T operator % ( Vec4 const &p ) const { return Dot(p); }						//!< Dot product

	//!@name Swizzling Methods
	CY_NODISCARD Vec2<T> XX() const { return Vec2<T>(x,x); }
	CY_NODISCARD Vec2<T> XY() const { return Vec2<T>(*this); }
	CY_NODISCARD Vec2<T> XZ() const { return Vec2<T>(x,z); }
	CY_NODISCARD Vec2<T> YX() const { return Vec2<T>(y,x); }
	CY_NODISCARD Vec2<T> YY() const { return Vec2<T>(y,y); }
	CY_NODISCARD Vec2<T> YZ() const { return Vec2<T>(y,z); }
	CY_NODISCARD Vec2<T> ZX() const { return Vec2<T>(z,x); }
	CY_NODISCARD Vec2<T> ZY() const { return Vec2<T>(z,y); }
	CY_NODISCARD Vec2<T> ZZ() const { return Vec2<T>(z,z); }

	CY_NODISCARD Vec3<T> XXX() const { return Vec3<T>(x,x,x); }
	CY_NODISCARD Vec3<T> XXY() const { return Vec3<T>(x,x,y); }
	CY_NODISCARD Vec3<T> XXZ() const { return Vec3<T>(x,x,z); }
	CY_NODISCARD Vec3<T> XYX() const { return Vec3<T>(x,y,x); }
	CY_NODISCARD Vec3<T> XYY() const { return Vec3<T>(x,y,y); }
	CY_NODISCARD Vec3<T> XYZ() const { return Vec3<T>(*this); }
	CY_NODISCARD Vec3<T> XZX() const { return Vec3<T>(x,z,x); }
	CY_NODISCARD Vec3<T> XZY() const { return Vec3<T>(x,z,y); }
	CY_NODISCARD Vec3<T> XZZ() const { return Vec3<T>(x,z,z); }

	CY_NODISCARD Vec3<T> YXX() const { return Vec3<T>(y,x,x); }
	CY_NODISCARD Vec3<T> YXY() const { return Vec3<T>(y,x,y); }
	CY_NODISCARD Vec3<T> YXZ() const { return Vec3<T>(y,x,z); }
	CY_NODISCARD Vec3<T> YYX() const { return Vec3<T>(y,y,x); }
	CY_NODISCARD Vec3<T> YYY() const { return Vec3<T>(y,y,y); }
	CY_NODISCARD Vec3<T> YYZ() const { return Vec3<T>(y,y,z); }
	CY_NODISCARD Vec3<T> YZX() const { return Vec3<T>(y,z,x); }
	CY_NODISCARD Vec3<T> YZY() const { return Vec3<T>(y,z,y); }
	CY_NODISCARD Vec3<T> YZZ() const { return Vec3<T>(y,z,z); }

	CY_NODISCARD Vec3<T> ZXX() const { return Vec3<T>(z,x,x); }
	CY_NODISCARD Vec3<T> ZXY() const { return Vec3<T>(z,x,y); }
	CY_NODISCARD Vec3<T> ZXZ() const { return Vec3<T>(z,x,z); }
	CY_NODISCARD Vec3<T> ZYX() const { return Vec3<T>(z,y,x); }
	CY_NODISCARD Vec3<T> ZYY() const { return Vec3<T>(z,y,y); }
	CY_NODISCARD Vec3<T> ZYZ() const { return Vec3<T>(z,y,z); }
	CY_NODISCARD Vec3<T> ZZX() const { return Vec3<T>(z,z,x); }
	CY_NODISCARD Vec3<T> ZZY() const { return Vec3<T>(z,z,y); }
	CY_NODISCARD Vec3<T> ZZZ() const { return Vec3<T>(z,z,z); }

	CY_NODISCARD Vec3<T> GetNonHomogeneous() const { return Vec3<T>(*this)/w; }
};

//-------------------------------------------------------------------------------

// Definitions of the conversion constructors
template <typename T, int N>                       Vec<T,N>::Vec( Vec2<T> const &p ) { if (N<=2) { MemCopy   (elem,&p.x,N); } else { MemCopy   (elem,&p.x,2); MemClear(elem,N-2); } }
template <typename T, int N>                       Vec<T,N>::Vec( Vec3<T> const &p ) { if (N<=3) { MemCopy   (elem,&p.x,N); } else { MemCopy   (elem,&p.x,3); MemClear(elem,N-3); } }
template <typename T, int N>                       Vec<T,N>::Vec( Vec4<T> const &p ) { if (N<=4) { MemCopy   (elem,&p.x,N); } else { MemCopy   (elem,&p.x,4); MemClear(elem,N-4); } }
template <typename T, int N> template <typename S> Vec<T,N>::Vec( Vec2<S> const &p ) { if (N<=2) { MemConvert(elem,&p.x,N); } else { MemConvert(elem,&p.x,2); MemClear(elem,N-2); } }
template <typename T, int N> template <typename S> Vec<T,N>::Vec( Vec3<S> const &p ) { if (N<=3) { MemConvert(elem,&p.x,N); } else { MemConvert(elem,&p.x,3); MemClear(elem,N-3); } }
template <typename T, int N> template <typename S> Vec<T,N>::Vec( Vec4<S> const &p ) { if (N<=4) { MemConvert(elem,&p.x,N); } else { MemConvert(elem,&p.x,4); MemClear(elem,N-4); } }
template <typename T>                              Vec2<T>::Vec2( Vec3<T> const &p ) : x(  p.x ), y(  p.y )            {}
template <typename T>                              Vec2<T>::Vec2( Vec4<T> const &p ) : x(  p.x ), y(  p.y )            {}
template <typename T>                              Vec3<T>::Vec3( Vec4<T> const &p ) : x(  p.x ), y(  p.y ), z(  p.z ) {}
template <typename T>        template <typename S> Vec2<T>::Vec2( Vec3<S> const &p ) : x(T(p.x)), y(T(p.y))            {}
template <typename T>        template <typename S> Vec2<T>::Vec2( Vec4<S> const &p ) : x(T(p.x)), y(T(p.y))            {}
template <typename T>        template <typename S> Vec3<T>::Vec3( Vec4<S> const &p ) : x(T(p.x)), y(T(p.y)), z(T(p.z)) {}

//-------------------------------------------------------------------------------

/// !@name Support functions

template <typename T> inline Vec2<T> Normalize( Vec2<T> const &v ) { return v.GetNormalized(); }
template <typename T> inline Vec3<T> Normalize( Vec3<T> const &v ) { return v.GetNormalized(); }
template <typename T> inline Vec4<T> Normalize( Vec4<T> const &v ) { return v.GetNormalized(); }

//-------------------------------------------------------------------------------

typedef Vec2<float>  Vec2f;	//!< 2D vector class with float type elements
typedef Vec3<float>  Vec3f;	//!< 3D vector class with float type elements
typedef Vec4<float>  Vec4f;	//!< 4D vector class with float type elements

typedef Vec2<double> Vec2d;	//!< 2D vector class with double type elements
typedef Vec3<double> Vec3d;	//!< 3D vector class with double type elements
typedef Vec4<double> Vec4d;	//!< 4D vector class with double type elements

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::Vec2f  cyVec2f;		//!< 2D vector class with float type elements
typedef cy::Vec3f  cyVec3f;		//!< 3D vector class with float type elements
typedef cy::Vec4f  cyVec4f;		//!< 4D vector class with float type elements

typedef cy::Vec2d  cyVec2d;		//!< 2D vector class with double type elements
typedef cy::Vec3d  cyVec3d;		//!< 3D vector class with double type elements
typedef cy::Vec4d  cyVec4d;		//!< 4D vector class with double type elements

//-------------------------------------------------------------------------------

#endif

