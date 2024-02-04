// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyMatrix.h 
//! \author Cem Yuksel
//! 
//! \brief  2x2, 3x3, 3x4, and 4x4 matrix classes
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

#ifndef _CY_MATRIX_H_INCLUDED_
#define _CY_MATRIX_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyVector.h"

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

// Forward declarations
//!	\cond HIDDEN_SYMBOLS
template <typename T> class Matrix3;
template <typename T> class Matrix34;
template <typename T> class Matrix4;
//! \endcond

//-------------------------------------------------------------------------------

#define _CY_VEC_DEFAULT_ERROR_TOLERANCE 0.0001

//-------------------------------------------------------------------------------

// The macros below help the MSVC compiler with auto-vectorization
#define _CY_FORi(start,end,loop) _CY_IVDEP_FOR ( int i=(start); i<(end); ++i ) { loop; }
#define _CY_FORi0(end,loop)      _CY_FORi(0,end,loop)
#define _CY_FOR_4i(loop)         _CY_FORi0(4,loop)
#define _CY_FOR_8i(loop)         _CY_FOR_4i(loop)  _CY_FORi(4,8,loop)
#define _CY_FOR_9i(loop)         _CY_FOR_8i(loop)  { int i=8; loop; }
#define _CY_FOR_12i(loop)        _CY_FOR_8i(loop)  _CY_FORi( 8,12,loop)
#define _CY_FOR_16i(loop)        _CY_FOR_12i(loop) _CY_FORi(12,16,loop)

//-------------------------------------------------------------------------------

//! 2x2 matrix class.
//!
//! Its data stores 4-value array of column-major matrix elements.
//! You can use Matrix2 with Vec2<T> to transform 2D points.

template <typename T>
class Matrix2
{
	CY_NODISCARD friend Matrix2 operator * ( T value, Matrix2 const &right ) { Matrix2 r; _CY_FOR_4i( r.cell[i] = value * right.cell[i] ); return r; }	//!< multiply matrix by a value
	CY_NODISCARD friend Matrix2 operator + ( T value, Matrix2 const &right ) { return Matrix2(value+right.cell[0], right.cell[2], right.cell[1],value+right.cell[3]); }	//!< add a value times identity matrix to a matrix
	CY_NODISCARD friend Matrix2 operator - ( T value, Matrix2 const &right ) { return Matrix2(value-right.cell[0],-right.cell[2],-right.cell[1],value-right.cell[3]); }	//!< subtract matrix from a value times identity matrix
	CY_NODISCARD friend Matrix2 Inverse( Matrix2 const &m ) { return m.GetInverse(); }	//!< return the inverse of the matrix

public:

	//! Elements of the matrix are column-major: \n
	//! | 0  2 | \n
	//! | 1  3 | \n
#ifdef __cpp_unrestricted_unions
	union {
		T       cell[4];
		Vec2<T> column[2];	// column vectors
	};
#else
	T cell[4];
#endif

	//////////////////////////////////////////////////////////////////////////
	//!@name Constructors

	Matrix2() CY_CLASS_FUNCTION_DEFAULT										//!< Default constructor
	template <typename TT> explicit Matrix2<T>( const Matrix2<TT> &matrix ) { MemConvert(cell,matrix.cell,4); }	//!< Copy constructor for different types
	explicit Matrix2( T const * restrict values ) { Set(values); }			//!< Initialize the matrix using an array of 4 values
	explicit Matrix2( T v )                       { SetScale(v); }			//!< Initialize the matrix as identity scaled by v
	explicit Matrix2( Vec2<T> const &x, Vec2<T> const &y ) { Set(x,y); }	//!< Initialize the matrix using two vectors as columns
	explicit Matrix2( Matrix3 <T> const &m );
	explicit Matrix2( Matrix34<T> const &m );
	explicit Matrix2( Matrix4 <T> const &m );

	//! Constructor using row-major order for initialization
	Matrix2( T c00, T c01,
		     T c10, T c11 )
	{
		cell[0] = c00;   cell[2] = c01;
		cell[1] = c10;   cell[3] = c11;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Set & Get Methods

	void Zero    ()       { MemClear(cell,4); }										//!< Set all the values as zero
	bool IsZero  () const { return Column(0).IsZero  () && Column(1).IsZero  (); }	//!< Returns true if the matrix is exactly zero
	bool IsFinite() const { return Column(0).IsFinite() && Column(1).IsFinite(); }	//!< Returns true if all components are finite real numbers.
	void Get( T       * restrict values ) { MemCopy(values,cell,4); }				//!< Copies the matrix cell to the given values array of size 4
	void Set( T const * restrict values ) { MemCopy(cell,values,4); }				//!< Set Matrix using an array of 4 values
	void Set( Vec2<T> const &x, Vec2<T> const &y ) { x.Get(cell); y.Get(cell+2); }	//!< Set Matrix using two vectors as columns
	void SetIdentity()    { SetScale(T(1)); }										//!< Converts the matrix to an identity matrix
	void SetTensorProduct( Vec2<T> const &v0, Vec2<T> const &v1 )					//!< Sets the matrix as the tensor product (outer product) of two vectors
	{
		_CY_IVDEP_FOR ( int i=0; i<2; ++i ) cell[  i] = v0[i] * v1.x;
		_CY_IVDEP_FOR ( int i=0; i<2; ++i ) cell[2+i] = v0[i] * v1.y;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Affine transformations

	//! Sets a uniform scale matrix
	void SetScale( T uniformScale ) { SetScale(uniformScale,uniformScale); }
	//! Sets a scale matrix
	void SetScale( T scaleX, T scaleY ) { cell[0]=scaleX; cell[1]=0; cell[2]=0; cell[3]=scaleY;}
	//! Sets a scale matrix
	void SetScale( Vec2<T> const &scale ) { SetScale(scale.x,scale.y); }
	//! Set a rotation matrix by angle
	void SetRotation( T angle ) { SetRotation( std::sin(angle), std::cos(angle) ); }
	//! Set a rotation matrix by cos and sin of angle
	void SetRotation( T sinAngle, T cosAngle ) { cell[0]=cosAngle; cell[1]=-sinAngle; cell[2]=sinAngle; cell[3]=cosAngle; }
	//! Sets a Cartesian coordinate frame using the given x direction. x must be a unit vector.
	void SetCartesianFrameX( Vec2<T> const &x ) { Set( x, x.GetPerpendicular() ); }
	//! Sets a Cartesian coordinate frame using the given y direction. y must be a unit vector.
	void SetCartesianFrameY( Vec2<T> const &y ) { Set( -y.GetPerpendicular(), y ); }


	//////////////////////////////////////////////////////////////////////////
	//!@name Set Row, Column, or Diagonal

	void SetRow     ( int ri, T x, T y )          { cell[ri]=x; cell[ri+2]=y; }			//!< Sets a row of the matrix
	void SetRow     ( int ri, Vec2<T> const &v )  { SetRow(ri,v.x,v.y); }				//!< Sets a row of the matrix
	void SetColumn  ( int ci, T x, T y )          { Column(ci).Set(x,y); }				//!< Sets a column of the matrix
	void SetColumn  ( int ci, Vec2<T> const &v )  { Column(ci)=v; }						//!< Sets a column of the matrix
	void SetDiagonal( T xx, T yy )                { cell[0]=xx; cell[3]=yy; }			//!< Sets the diagonal values of the matrix
	void SetDiagonal( Vec2<T> const &p )          { SetDiagonal( p.x, p.y ); }			//!< Sets the diagonal values of the matrix
	void SetDiagonal( T const * restrict values ) { SetDiagonal(values[0],values[1]); }	//!< Sets the diagonal values of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Row, Column, or Diagonal

#ifdef __cpp_unrestricted_unions
	CY_NODISCARD Vec2<T>       * Columns()               { return column; }
	CY_NODISCARD Vec2<T> const * Columns()         const { return column; }
	CY_NODISCARD Vec2<T>       & Column ( int ci )       { return column[ci]; }
	CY_NODISCARD Vec2<T> const & Column ( int ci ) const { return column[ci]; }
#else
	CY_NODISCARD Vec2<T>       * Columns()               { return ((Vec2<T>*)cell); }
	CY_NODISCARD Vec2<T> const * Columns()         const { return ((Vec2<T>*)cell); }
	CY_NODISCARD Vec2<T>       & Column ( int ci )       { return Columns()[ci]; }
	CY_NODISCARD Vec2<T> const & Column ( int ci ) const { return Columns()[ci]; }
#endif
	CY_NODISCARD Vec2<T>         GetRow ( int ri ) const { return Vec2<T>( cell[ri], cell[ri+2] ); }		//!< Returns a row of the matrix
	CY_NODISCARD Vec2<T>         GetDiagonal()     const { return Vec2<T>( cell[0], cell[3] ); }			//!< Returns the diagonal of the matrix
	CY_NODISCARD Matrix2         GetRotation()     const { Matrix2 s, r; GetComponents(s,r); return r; }	//!< Returns the rotation portion of the transformation

	//! Returns the scale portion of the transformation.
	//! The returned matrix is symmetric, but not necessarily diagonal, and it can include non-uniform scale.
	CY_NODISCARD Matrix2 GetScale() const
	{
		Matrix2 trns = GetTranspose();
		Matrix2 u2 = *this * trns;
		Vec2<T> v0, v1;
		u2.GetEigenvectors( v0, v1 );
		Matrix2 v(v0,v1);
		Matrix2 vt = v.GetTranspose();
		Matrix2 d2 = vt * (*this) * v;	// the result is a diagonal matrix
		Vec2<T> diag = d2.GetDiagonal();
		Matrix2 d;
		d.SetScale(diag);
		return v * d * vt;
	}

	//! Returns the average scale factor
	CY_NODISCARD T GetAvrgScale() const 
	{
		T det = cell[0]*cell[3]-cell[2]*cell[1];
		T s = Sqrt( std::abs(det) );
		return det >= 0 ? s : -s;
	}

	void GetComponents( Matrix2<T> &scale, Matrix2<T> &rotation ) const { scale = GetScale(); rotation = *this * scale.GetInverse(); }	//!< Returns separate transformation components

	//////////////////////////////////////////////////////////////////////////
	//!@name Comparison Operators

	CY_NODISCARD bool operator == ( Matrix2 const &right ) const { _CY_FOR_4i( if ( cell[i] != right.cell[i] ) return false ); return true;  } //!< compare equal
	CY_NODISCARD bool operator != ( Matrix2 const &right ) const { _CY_FOR_4i( if ( cell[i] != right.cell[i] ) return true  ); return false; } //!< compare not equal


	//////////////////////////////////////////////////////////////////////////
	//!@name Access Operators

	CY_NODISCARD T&        operator () ( int ri, int ci )       { assert( ri>=0 && ri<2 && ci>=0 && ci<2 ); return cell[ ci*2 + ri ]; }	//!< subscript operator
	CY_NODISCARD T const & operator () ( int ri, int ci ) const { assert( ri>=0 && ri<2 && ci>=0 && ci<2 ); return cell[ ci*2 + ri ]; }	//!< constant subscript operator
	CY_NODISCARD T&        operator [] ( int i )                { assert( i>=0 && i<4 ); return cell[i]; }								//!< subscript operator
	CY_NODISCARD T const & operator [] ( int i )          const { assert( i>=0 && i<4 ); return cell[i]; }								//!< constant subscript operator
	
	//////////////////////////////////////////////////////////////////////////
	//!@name Unary and Binary Operators

	// Unary operators
	CY_NODISCARD Matrix2 operator - () const { Matrix2 r; _CY_FOR_4i( r.cell[i] = -cell[i] ); return r; }	//!< negative matrix

	// Binary operators
	CY_NODISCARD Matrix2 operator * ( T       const  value ) const { Matrix2 r; _CY_FOR_4i( r.cell[i] = cell[i] * value         ); return r; }	//!< multiply matrix by a value
	CY_NODISCARD Matrix2 operator / ( T       const  value ) const { Matrix2 r; _CY_FOR_4i( r.cell[i] = cell[i] / value         ); return r; }	//!< divide matrix by a value
	CY_NODISCARD Matrix2 operator + ( Matrix2 const &right ) const { Matrix2 r; _CY_FOR_4i( r.cell[i] = cell[i] + right.cell[i] ); return r; }	//!< add two Matrices
	CY_NODISCARD Matrix2 operator - ( Matrix2 const &right ) const { Matrix2 r; _CY_FOR_4i( r.cell[i] = cell[i] - right.cell[i] ); return r; }	//!< subtract one Matrix2 from another
	CY_NODISCARD Matrix2 operator * ( Matrix2 const &right ) const	//!< multiply a matrix with another
	{
		Matrix2 r;
		r[0] = cell[0] * right.cell[0] + cell[2] * right.cell[1];
		r[1] = cell[1] * right.cell[0] + cell[3] * right.cell[1];
		r[2] = cell[0] * right.cell[2] + cell[2] * right.cell[3];
		r[3] = cell[1] * right.cell[2] + cell[3] * right.cell[3];
		return r;
	}
	CY_NODISCARD Vec2<T> operator * ( Vec2<T> const &p ) const { return Vec2<T>( p.x*cell[0] + p.y*cell[2], p.x*cell[1] + p.y*cell[3] ); }

	CY_NODISCARD Matrix2 operator + ( T value ) const { Matrix2 r=*this; r.cell[0]+=value; r.cell[3]+=value; return r; }	//!< add a value times identity matrix
	CY_NODISCARD Matrix2 operator - ( T value ) const { Matrix2 r=*this; r.cell[0]-=value; r.cell[3]-=value; return r; }	//!< subtract a value times identity matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Assignment Operators

	Matrix2 const & operator += ( Matrix2 const &right ) { _CY_FOR_4i( cell[i] += right.cell[i] ); return *this; }	//!< add two Matrices modify this
	Matrix2 const & operator -= ( Matrix2 const &right ) { _CY_FOR_4i( cell[i] -= right.cell[i] ); return *this; }	//!< subtract one Matrix2 from another matrix and modify this matrix
	Matrix2 const & operator *= ( Matrix2 const &right ) { *this = operator*(right);               return *this; }	//!< multiply a matrix with another matrix and modify this matrix
	Matrix2 const & operator *= ( T       const  value ) { _CY_FOR_4i( cell[i] *= value );         return *this; }	//!< multiply a matrix with a value modify this matrix
	Matrix2 const & operator /= ( T       const  value ) { _CY_FOR_4i( cell[i] /= value );         return *this; }	//!< divide the matrix by a value modify the this matrix
	Matrix2 const & operator += ( T       const  value ) { cell[0]+=value; cell[3]+=value;         return *this; }	//!< add a value times identity matrix
	Matrix2 const & operator -= ( T       const  value ) { cell[0]-=value; cell[3]-=value;         return *this; }	//!< subtract a value times identity matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Other Methods

	void Transpose() { T tmp=cell[0]; cell[0]=cell[3]; cell[3]=tmp; }	//!< Transpose this matrix
	CY_NODISCARD Matrix2 GetTranspose() const							//!< Returns the transpose of this matrix
	{
		Matrix2 m;
		m.cell[0] = cell[0];   m.cell[1] = cell[2];
		m.cell[2] = cell[1];   m.cell[3] = cell[3];
		return m;
	}

	//! Multiply the give vector with the transpose of the matrix
	CY_NODISCARD Vec2<T> TransposeMult( Vec2<T> const &p ) const { return Vec2<T>( p.x*cell[0] + p.y*cell[1], p.x*cell[2] + p.y*cell[3] ); }

	CY_NODISCARD Matrix2 TransposeMult( Matrix2 const & right ) const //!< Multiply a matrix by the transpose of this one (i.e. this^T * right).
	{
		Matrix2 r;
		r[0] = cell[0] * right.cell[0] + cell[1] * right.cell[1];
		r[1] = cell[2] * right.cell[0] + cell[3] * right.cell[1];
		r[2] = cell[0] * right.cell[2] + cell[1] * right.cell[3];
		r[3] = cell[2] * right.cell[2] + cell[3] * right.cell[3];
		return r;
	}
	CY_NODISCARD Matrix2 MultTranspose( Matrix2 const & right ) const //!< Multiply the transpose of a matrix by this one (i.e. this * right^T).
	{
		Matrix2 r;
		r[0] = cell[0] * right.cell[0] + cell[2] * right.cell[2];
		r[1] = cell[1] * right.cell[0] + cell[3] * right.cell[2];
		r[2] = cell[0] * right.cell[1] + cell[2] * right.cell[3];
		r[3] = cell[1] * right.cell[1] + cell[3] * right.cell[3];
		return r;
	}

	CY_NODISCARD Matrix2 TransposeMultSelf() const { return TransposeMult(*this); } //!< Multiply the transpose of this matrix with itself (i.e. this^T * this).
	CY_NODISCARD Matrix2 MultSelfTranspose() const { return MultTranspose(*this); } //!< Multiply the matrix with its transpose (i.e. this * this^T).

	CY_NODISCARD T GetTrace() const { return cell[0] + cell[3]; }	//!< return the Trace of this matrix

	CY_NODISCARD T GetDeterminant() const { return cell[0]*cell[3]-cell[2]*cell[1]; }	//!< Get the determinant of this matrix

	void Invert()	//!< Invert this matrix
	{
		T det = GetDeterminant();
		T cell3 =  cell[0] / det;
		cell[0] =  cell[3] / det;
		cell[1] = -cell[1] / det;
		cell[2] = -cell[2] / det;
		cell[3] =  cell3;
	}
	CY_NODISCARD Matrix2 GetInverse() const	//!< Get the inverse of this matrix
	{
		T det = GetDeterminant();
		Matrix2 inv;
		inv.cell[0] =  cell[3] / det;
		inv.cell[1] = -cell[1] / det;
		inv.cell[2] = -cell[2] / det;
		inv.cell[3] =  cell[0] / det;
		return inv;
	}

	//! Removes the scale component of the matrix by normalizing each column.
	//! The resulting matrix can contain shear, if it originally contained non-uniform scale and rotation.
	void Normalize() { Column(0).Normalize(); Column(1).Normalize(); }

	//! Orthogonalizes the matrix and removes the scale component, preserving the x direction
	void OrthogonalizeX()
	{
		Column(0).Normalize();
		Column(1) -= Column(0) * (Column(1) % Column(0));
		Column(1).Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the y direction
	void OrthogonalizeY()
	{
		Column(1).Normalize();
		Column(0) -= Column(1) * (Column(0) % Column(1));
		Column(0).Normalize();
	}

	//! Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const { return std::abs(cell[0] - T(1)) < tollerance && std::abs(cell[1]) < tollerance && std::abs(cell[2]) < tollerance && std::abs(cell[3] - T(1)) < tollerance; }

	//! Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const { return std::abs(cell[0] - cell[2]) < tollerance; }

	//! Returns if the matrix is diagonal.
	bool IsDiagonal( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const { return std::abs(cell[1]) + std::abs(cell[2]) < tollerance*2; }

	//! Returns the eigenvalues of the matrix.
	//! The eigenvalues are ordered, such that the first one is larger.
	CY_NODISCARD Vec2<T> GetEigenvalues() const
	{
		T t = GetTrace();
		T d = GetDeterminant();
		T a = t*t*T(0.25) - d;
		T s = SqrtSafe<T>(a);
		Vec2<T> lambda;
		lambda.x = t * T(0.5) + s;
		lambda.y = t * T(0.5) - s;
		return lambda;
	}

	//! Returns the eigenvalues and sets the given vectors as the eigenvectors of the matrix.
	//! The eigenvalues are ordered, such that the first one is larger.
	//! The given tollerance is used for checking whether the eigenvalues are the same.
	Vec2<T> GetEigenvectors( Vec2<T> &evec0, Vec2<T> &evec1, T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const
	{
		Vec2<T> lambda = GetEigenvalues();
		if ( std::abs(lambda.x - lambda.y) < tollerance ) {
			evec0 = Column(0);
			evec1 = Column(1);
		} else {
			Matrix2 v0( cell[0]-lambda.y, cell[1], cell[2], cell[3]-lambda.y );
			Matrix2 v1( cell[0]-lambda.x, cell[1], cell[2], cell[3]-lambda.x );
			evec0 = v0.Column(0) + v0.Column(1);
			evec1 = v1.Column(0) + v1.Column(1);
		}
		return lambda;
	}

	//! Singular value decomposition (SVD).
	//! Returns the SVD of the matrix, where U and V are orthogonal matrices and 
	//! S is the diagonal elements of a diagonal matrix (including zeros),
	//! such that this matrix A = U S V^T.
	void SingularValueDecomposition( Matrix2<T> &U, Vec2<T> &S, Matrix2<T> &V )
	{
		Matrix2 AAT = MultSelfTranspose();
		Vec2<T> lambda = AAT.GetEigenvectors( U.Column(0), U.Column(1) );
		S = (lambda.Abs()).Sqrt();
		U.Normalize();
		Matrix2 ATA = TransposeMultSelf();
		AAT.GetEigenvectors( V.Column(0), V.Column(1) );
		V.Normalize();
	}

	//////////////////////////////////////////////////////////////////////////
	//!@name Static Methods

	//! Returns an identity matrix
	CY_NODISCARD static Matrix2 Identity() { T c[] = { 1,0, 0,1 }; return Matrix2(c); }
	//! Returns a rotation matrix by angle in radians
	CY_NODISCARD static Matrix2 Rotation( T angle ) { Matrix2 m; m.SetRotation(angle); return m; }
	//! Returns a rotation matrix by cos and sin of the rotation angle
	CY_NODISCARD static Matrix2 Rotation( T cosAngle, T sinAngle ) { Matrix2 m; m.SetRotation(cosAngle,sinAngle); return m; }
	//! Returns a uniform scale matrix
	CY_NODISCARD static Matrix2 Scale( T uniformScale ) { Matrix2 m; m.SetScale(uniformScale); return m; }
	//! Returns a scale matrix
	CY_NODISCARD static Matrix2 Scale( T scaleX, T scaleY ) { Matrix2 m; m.SetScale(scaleX,scaleY); return m; }
	//! Returns a scale matrix
	CY_NODISCARD static Matrix2 Scale( Vec2<T> const &scale ) { Matrix2 m; m.SetScale(scale); return m; }
	//! Returns the tensor product (outer product) matrix of two vectors
	CY_NODISCARD static Matrix2 TensorProduct( Vec2<T> const &v0, Vec2<T> const &v1 ) { Matrix2 m; m.SetTensorProduct(v0,v1); return m; }

	//////////////////////////////////////////////////////////////////////////
};

//-------------------------------------------------------------------------------

#ifdef CY_NONVECTORIZED_MATRIX3
# define _CY_INIT_MATRIX3_VECTORIZATION const int N = 3; T const *cell_6 = cell + 6;
#else
# define _CY_INIT_MATRIX3_VECTORIZATION const int N = 4; T cell_6[4] = { cell[6], cell[7], cell[8], cell[8] };
#endif

//-------------------------------------------------------------------------------

//! 3x3 matrix class.
//!
//! Its data stores 9-value array of column-major matrix elements.
//! You can use Matrix3 with Vec3<T> to transform 3D points.

template <typename T>
class Matrix3
{
	CY_NODISCARD friend Matrix3 operator * ( T value, Matrix3 const &right ) { Matrix3 r; _CY_FOR_9i( r.cell[i] = value * right.cell[i] ); return r; }	//!< multiply matrix by a value
	CY_NODISCARD friend Matrix3 operator + ( T value, Matrix3 const &right ) { Matrix3 r= right; r.cell[0]+=value; r.cell[4]+=value; r.cell[8]+=value; return r; }	//!< add a value times identity matrix to a matrix
	CY_NODISCARD friend Matrix3 operator - ( T value, Matrix3 const &right ) { Matrix3 r=-right; r.cell[0]+=value; r.cell[4]+=value; r.cell[8]+=value; return r; }	//!< subtract a matrix from a value times identity matrix
	CY_NODISCARD friend Matrix3 Inverse( Matrix3 const &m ) { return m.GetInverse(); }	//!< return the inverse of the matrix

public:

	//! Elements of the matrix are column-major: \n
	//! | 0  3  6 | \n
	//! | 1  4  7 | \n
	//! | 2  5  8 | \n
#ifdef __cpp_unrestricted_unions
	union {
		T       cell[9];
		Vec3<T> column[3];	// column vectors
	};
#else
	T cell[9];
#endif

	//////////////////////////////////////////////////////////////////////////
	//!@name Constructors

	Matrix3() CY_CLASS_FUNCTION_DEFAULT															//!< Default constructor
	template <typename TT> explicit Matrix3<T>( Matrix3<TT> const &matrix ) { MemConvert(cell,matrix.cell,9); }	//!< Copy constructor for different types
	explicit Matrix3( T const * restrict values ) { Set(values); }								//!< Initialize the matrix using an array of 9 values
	explicit Matrix3( T v )                       { SetScale(v); }								//!< Initialize the matrix as identity scaled by v
	explicit Matrix3( Vec3<T> const &x, Vec3<T> const &y, Vec3<T> const &z ) { Set(x,y,z); }	//!< Initialize the matrix using x,y,z vectors as columns
	explicit Matrix3( Matrix2 <T> const &m ) { Column(0).Set(m.Column(0),0); Column(1).Set(m.Column(1),0); Column(2).Set(0,0,1); }
	explicit Matrix3( Matrix34<T> const &m );
	explicit Matrix3( Matrix4 <T> const &m );

	//! Constructor using row-major order for initialization
	Matrix3( T c00, T c01, T c02,
		     T c10, T c11, T c12,
		     T c20, T c21, T c22 )
	{
		cell[0] = c00;   cell[3] = c01;   cell[6] = c02;
		cell[1] = c10;   cell[4] = c11;   cell[7] = c12;
		cell[2] = c20;   cell[5] = c21;   cell[8] = c22;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Set & Get Methods

	void Zero    ()       { MemClear(cell,9); }																	//!< Set all the values as zero
	bool IsZero  () const { return Column(0).IsZero  () && Column(1).IsZero  () && Column(2).IsZero  (); }		//!< Returns true if the matrix is exactly zero
	bool IsFinite() const { return Column(0).IsFinite() && Column(1).IsFinite() && Column(2).IsFinite(); }		//!< Returns true if all components are finite real numbers.
	void Get( T       * restrict values ) { MemCopy(values,cell,9); }											//!< Copies the matrix cell to the given values array of size 9
	void Set( T const * restrict values ) { MemCopy(cell,values,9); }											//!< Set matrix using an array of 9 values
	void Set( Vec3<T> const &x, Vec3<T> const &y, Vec3<T> const &z ) { Column(0)=x; Column(1)=y; Column(2)=z; }	//!< Set matrix using x,y,z vectors as columns
	void SetIdentity()    { SetScale(T(1)); }																	//!< Converts the matrix to an identity matrix
	void SetTensorProduct( Vec3<T> const &v0, Vec3<T> const &v1 )												//!< Sets the matrix as the tensor product (outer product) of two vectors
	{
		_CY_IVDEP_FOR ( int i=0; i<3; ++i ) cell[  i] = v0[i] * v1.x;
		_CY_IVDEP_FOR ( int i=0; i<3; ++i ) cell[3+i] = v0[i] * v1.y;
		_CY_IVDEP_FOR ( int i=0; i<3; ++i ) cell[6+i] = v0[i] * v1.z;
	}
	//! Matrix representation of the cross product ( a x b)
	void SetCrossProd( Vec3<T> const &p ) { cell[0]=T(0); cell[1]=p.z; cell[2]=-p.y; cell[3]=-p.z; cell[4]=T(0); cell[5]=p.x; cell[6]=p.y; cell[7]=-p.x; cell[8]=T(0); }


	//////////////////////////////////////////////////////////////////////////
	//!@name Affine transformations

	//! Sets a uniform scale matrix
	void SetScale( T uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	//! Sets a scale matrix
	void SetScale( T scaleX, T scaleY, T scaleZ )
	{
		cell[0] = scaleX; cell[1] = 0;      cell[2]=0;     
		cell[3] = 0;      cell[4] = scaleY; cell[5]=0;     
		cell[6] = 0;      cell[7] = 0;      cell[8]=scaleZ;
	}
	//! Sets a scale matrix
	void SetScale( Vec3<T> const &scale ) { SetScale(scale.x,scale.y,scale.z); }
	//! Set as rotation matrix around x axis
	void SetRotationX( T angle ) { SetRotationX( std::sin(angle), std::cos(angle) ); }
	//! Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( T sinAngle, T cosAngle )
	{
		cell[0] = T(1);   cell[1] =  T(0);    cell[2] = T(0); 
		cell[3] = T(0);   cell[4] =  cosAngle;   cell[5] = sinAngle;
		cell[6] = T(0);   cell[7] = -sinAngle;   cell[8] = cosAngle;
	}
	//! Set as rotation matrix around y axis
	void SetRotationY( T angle ) { SetRotationY( std::sin(angle), std::cos(angle) ); }
	//! Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( T sinAngle, T cosAngle )
	{
		cell[0] = cosAngle;   cell[1] = T(0);   cell[2] = -sinAngle;
		cell[3] = T(0);       cell[4] = T(1);   cell[5] =  T(0); 
		cell[6] = sinAngle;   cell[7] = T(0);   cell[8] =  cosAngle;
	}
	//! Set as rotation matrix around z axis
	void SetRotationZ( T angle ) { SetRotationZ( std::sin(angle), std::cos(angle) ); }
	//! Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( T sinAngle, T cosAngle )
	{
		cell[0] =  cosAngle;   cell[1] = sinAngle;   cell[2] = T(0);
		cell[3] = -sinAngle;   cell[4] = cosAngle;   cell[5] = T(0);
		cell[6] =  T(0);       cell[7] = T(0);       cell[8] = T(1);
	}
	//! Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
	void SetRotationXYZ( T angleX, T angleY, T angleZ )
	{
		T sx = std::sin(angleX);
		T cx = std::cos(angleX);
		T sy = std::sin(angleY);
		T cy = std::cos(angleY);
		T sz = std::sin(angleZ);
		T cz = std::cos(angleZ);
		cell[0] = cy*cz; 		      cell[1] = cy*sz; 			    cell[2] =-sy;   
		cell[3] = cz*sx*sy - cx*sz;   cell[4] = cx*cz + sx*sy*sz;   cell[5] = cy*sx;
		cell[6] = cx*cz*sy + sx*sz;   cell[7] =-cz*sx + cx*sy*sz;   cell[8] = cx*cy;
	}
	//! Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
	void SetRotationZYX( T angleX, T angleY, T angleZ )
	{
		T sx = std::sin(angleX);
		T cx = std::cos(angleX);
		T sy = std::sin(angleY);
		T cy = std::cos(angleY);
		T sz = std::sin(angleZ);
		T cz = std::cos(angleZ);
		cell[0] =  cy*cz;   cell[1] = cx*sz + sx*sy*cz;   cell[2] = sx*sz - cx*sy*cz;
		cell[3] = -cy*sz;   cell[4] = cx*cz - sx*sy*sz;   cell[5] = sx*cz + cx*sy*sz;
		cell[6] =  sy;      cell[7] = -sx*cy;		      cell[8] = cx*cy;
	}
	//! Set a rotation matrix about the given axis by angle
	void SetRotation( Vec3<T> const &axis, T angle ) { SetRotation(axis,std::sin(angle),std::cos(angle)); }
	//! Set a rotation matrix about the given axis by cos and sin of angle
	void SetRotation( Vec3<T> const &axis, T sinAngle, T cosAngle )
	{
		T t = T(1) - cosAngle;
		Vec3<T> a = t * axis;
		T txy = a.x * axis.y;
		T txz = a.x * axis.z;
		T tyz = a.y * axis.z;
		Vec3<T> s = sinAngle * axis;
		cell[ 0] = a.x * axis.x + cosAngle;   cell[ 1] = txy + s.z;                 cell[ 2] = txz - s.y;
		cell[ 3] = txy - s.z;                 cell[ 4] = a.y * axis.y + cosAngle;   cell[ 5] = tyz + s.x;
		cell[ 6] = txz + s.y;                 cell[ 7] = tyz - s.x;                 cell[ 8] = a.z * axis.z + cosAngle;
	}
	//! Set a rotation matrix that sets [from] unit vector to [to] unit vector
	void SetRotation( Vec3<T> const &from, Vec3<T> const &to )
	{
		assert( from.IsFinite() && to.IsUnit() );
		Vec3<T> axis = from.Cross(to);
		T s = axis.Length();
		if ( s < T(0.000001) ) SetIdentity();
		else {
			T c = from.Dot(to);
			SetRotation(axis/s, s, c);
		}
	}
	//! Set view matrix using position, target and approximate up vector
	void SetView( Vec3<T> const &target, Vec3<T> const &up )
	{
		Vec3<T> f = target;
		f.Normalize();
		Vec3<T> s = f.Cross(up);
		s.Normalize();
		Vec3<T> u = s.Cross(f);
		cell[0] = s.x; cell[1] = u.x; cell[2] = -f.x;
		cell[3] = s.y; cell[4] = u.y; cell[5] = -f.y;
		cell[6] = s.z; cell[7] = u.z; cell[8] = -f.z;
	}
	//! Sets a Cartesian coordinate frame using the given x direction and an approximate y direction. x must be a unit vector.
	void SetCartesianFrameXY( Vec3<T> const &x, Vec3<T> const &y_approx ) { Vec3<T> z = x.Cross(y_approx); z.Normalize(); Vec3<T> y=z.Cross(x); Set(x,y,z); }
	//! Sets a Cartesian coordinate frame using the given x direction and an approximate z direction. x must be a unit vector.
	void SetCartesianFrameXZ( Vec3<T> const &x, Vec3<T> const &z_approx ) { Vec3<T> y = z_approx.Cross(x); y.Normalize(); Vec3<T> z=x.Cross(y); Set(x,y,z); }
	//! Sets a Cartesian coordinate frame using the given y direction and an approximate x direction. y must be a unit vector.
	void SetCartesianFrameYX( Vec3<T> const &y, Vec3<T> const &x_approx ) { Vec3<T> z = x_approx.Cross(y); z.Normalize(); Vec3<T> x=y.Cross(z); Set(x,y,z); }
	//! Sets a Cartesian coordinate frame using the given y direction and an approximate z direction. y must be a unit vector.
	void SetCartesianFrameYZ( Vec3<T> const &y, Vec3<T> const &z_approx ) { Vec3<T> x = y.Cross(z_approx); x.Normalize(); Vec3<T> z=x.Cross(y); Set(x,y,z); }
	//! Sets a Cartesian coordinate frame using the given z direction and an approximate x direction. z must be a unit vector.
	void SetCartesianFrameZX( Vec3<T> const &z, Vec3<T> const &x_approx ) { Vec3<T> y = z.Cross(x_approx); y.Normalize(); Vec3<T> x=y.Cross(z); Set(x,y,z); }
	//! Sets a Cartesian coordinate frame using the given z direction and an approximate y direction. z must be a unit vector.
	void SetCartesianFrameZY( Vec3<T> const &z, Vec3<T> const &y_approx ) { Vec3<T> x = y_approx.Cross(z); x.Normalize(); Vec3<T> y=z.Cross(x); Set(x,y,z); }


	//////////////////////////////////////////////////////////////////////////
	//!@name Set Row, Column, or Diagonal

	void SetRow     ( int ri, T x, T y, T z )     { cell[ri]=x; cell[ri+3]=y; cell[ri+6]=z; }		//!< Sets a row of the matrix
	void SetRow     ( int ri, Vec3<T> const &v )  { SetRow(ri,v.x,v.y,v.z); }						//!< Sets a row of the matrix
	void SetColumn  ( int ci, T x, T y, T z )     { Column(ci).Set(x,y,z); }						//!< Sets a column of the matrix
	void SetColumn  ( int ci, Vec3<T> const &v )  { Column(ci)=v; }									//!< Sets a column of the matrix
	void SetDiagonal( T xx, T yy, T zz )          { cell[0]=xx; cell[4]=yy; cell[8]=zz; }			//!< Sets the diagonal values of the matrix
	void SetDiagonal( Vec3<T> const &p )          { SetDiagonal( p.x, p.y, p.z ); }					//!< Sets the diagonal values of the matrix
	void SetDiagonal( T const * restrict values ) { SetDiagonal(values[0],values[1],values[2]); }	//!< Sets the diagonal values of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Row, Column, or Diagonal
	
#ifdef __cpp_unrestricted_unions
	CY_NODISCARD Vec3<T>       * Columns()               { return column; }
	CY_NODISCARD Vec3<T> const * Columns()         const { return column; }
	CY_NODISCARD Vec3<T>       & Column ( int ci )       { return column[ci]; }
	CY_NODISCARD Vec3<T> const & Column ( int ci ) const { return column[ci]; }
#else
	CY_NODISCARD Vec3<T>       * Columns()               { return ((Vec3<T>*)cell); }
	CY_NODISCARD Vec3<T> const * Columns()         const { return ((Vec3<T>*)cell); }
	CY_NODISCARD Vec3<T>       & Column ( int ci )       { return Columns()[ci]; }
	CY_NODISCARD Vec3<T> const & Column ( int ci ) const { return Columns()[ci]; }
#endif
	CY_NODISCARD Vec3<T>         GetRow ( int ri ) const { return Vec3<T>( cell[ri], cell[ri+3], cell[ri+6] ); }	//!< Returns a row of the matrix
	CY_NODISCARD Vec3<T>         GetDiagonal()     const { return Vec3<T>( cell[0], cell[4], cell[8] ); }			//!< Returns the diagonal of the matrix
	CY_NODISCARD Matrix3         GetRotation()     const { Matrix3 s, r; GetComponents(s,r); return r; }			//!< Returns the rotation portion of the transformation

	//! Returns the scale portion of the transformation.
	//! The returned matrix is symmetric, but not necessarily diagonal, and it can include non-uniform scale.
	CY_NODISCARD Matrix3 GetScale() const
	{
		Matrix3 trns = GetTranspose();
		Matrix3 u2 = *this * trns;
		Vec3<T> v0, v1, v2;
		u2.GetEigenvectors( v0, v1, v2 );
		Matrix3 v(v0,v1,v2);
		Matrix3 vt = v.GetTranspose();
		Matrix3 d2 = vt * (*this) * v;	// the result is a diagonal matrix
		Vec3<T> diag = d2.GetDiagonal();
		Matrix3 d;
		d.SetScale(diag);
		return v * d * vt;
	}

	//! Returns the average scale factor
	CY_NODISCARD T GetAvrgScale() const 
	{
		T det = cell[0] * ( cell[4] * cell[8] - cell[5] * cell[7] ) + 
		        cell[1] * ( cell[5] * cell[6] - cell[3] * cell[8] ) + 
		        cell[2] * ( cell[3] * cell[7] - cell[4] * cell[6] );
		T s = std::pow( std::abs(det), T(1)/T(3) );
		return det >= 0 ? s : -s;
	}

	void GetComponents( Matrix3<T> &scale, Matrix3<T> &rotation ) const { scale = GetScale(); rotation = *this * scale.GetInverse(); }	//!< Returns separate transformation components

	//////////////////////////////////////////////////////////////////////////
	//!@name Get Sub-matrix cell
	
	CY_NODISCARD Matrix2<T> GetSubMatrix2() const { Matrix2<T> m; MemCopy(m.cell,cell,2); MemCopy(m.cell+2,cell+3,2); return m; }	//!< Returns the 2x2 portion of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Comparison Operators

	CY_NODISCARD bool operator == ( Matrix3 const &right ) const { _CY_FOR_9i( if ( cell[i] != right.cell[i] ) return false ); return true;  } //!< compare equal
	CY_NODISCARD bool operator != ( Matrix3 const &right ) const { _CY_FOR_9i( if ( cell[i] != right.cell[i] ) return true  ); return false; } //!< compare not equal


	//////////////////////////////////////////////////////////////////////////
	//!@name Access Operators

	CY_NODISCARD T&        operator () ( int ri, int ci )       { assert( ri>=0 && ri<3 && ci>=0 && ci<3 ); return cell[ ci*3 + ri ]; }	//!< subscript operator
	CY_NODISCARD T const & operator () ( int ri, int ci ) const { assert( ri>=0 && ri<3 && ci>=0 && ci<3 ); return cell[ ci*3 + ri ]; }	//!< constant subscript operator
	CY_NODISCARD T&        operator [] ( int i )                { assert( i>=0 && i<9 ); return cell[i]; }								//!< subscript operator
	CY_NODISCARD T const & operator [] ( int i )          const { assert( i>=0 && i<9 ); return cell[i]; }								//!< constant subscript operator
	

	//////////////////////////////////////////////////////////////////////////
	//!@name Unary and Binary Operators

	// Unary operators
	CY_NODISCARD Matrix3 operator - () const { Matrix3 r; _CY_FOR_9i( r.cell[i] = -cell[i] ); return r; }	//!< negative matrix

	// Binary operators
	CY_NODISCARD Matrix3 operator * ( T       const  value ) const { Matrix3 r; _CY_FOR_9i( r.cell[i] = cell[i] * value         ); return r; }	//!< multiply matrix by a value
	CY_NODISCARD Matrix3 operator / ( T       const  value ) const { Matrix3 r; _CY_FOR_9i( r.cell[i] = cell[i] / value         ); return r; }	//!< divide matrix by a value
	CY_NODISCARD Matrix3 operator + ( Matrix3 const &right ) const { Matrix3 r; _CY_FOR_9i( r.cell[i] = cell[i] + right.cell[i] ); return r; }	//!< add two Matrices
	CY_NODISCARD Matrix3 operator - ( Matrix3 const &right ) const { Matrix3 r; _CY_FOR_9i( r.cell[i] = cell[i] - right.cell[i] ); return r; }	//!< subtract one Matrix3 from another

	CY_NODISCARD Matrix3 operator * ( Matrix3 const &right ) const	//!< multiply a matrix with another
	{
		_CY_INIT_MATRIX3_VECTORIZATION;
		Matrix3 rm;
		for ( int i=0; i<9; i+=3 ) {
			T a[4], b[4], c[4], d[4], r[4];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) a[j] = cell[  j] * right.cell[i  ];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) b[j] = cell[3+j] * right.cell[i+1];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) c[j] = cell_6[j] * right.cell[i+2];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) d[j] = a[j] + b[j];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) r[j] = d[j] + c[j];
			MemCopy( rm.cell+i, r, 3 );
		}
		return rm;
	}
	CY_NODISCARD Vec3<T> operator * ( Vec3<T> const &p ) const
	{
		_CY_INIT_MATRIX3_VECTORIZATION;
		//return Vec3<T>( p.x*cell[0] + p.y*cell[3] + p.z*cell[6], 
		//                p.x*cell[1] + p.y*cell[4] + p.z*cell[7],
		//                p.x*cell[2] + p.y*cell[5] + p.z*cell[8] );
		T a[4], b[4], c[4], d[4], r[4];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) a[i] = p[0] * cell[  i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) b[i] = p[1] * cell[3+i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) c[i] = p[2] * cell_6[i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) d[i] = a[i] + b[i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) r[i] = d[i] + c[i];	
		return Vec3<T>(r);
	}

	CY_NODISCARD Matrix3 operator + ( T value ) const { Matrix3 r=*this; r.cell[0]+=value; r.cell[4]+=value; r.cell[8]+=value; return r; }	//!< add a value times identity matrix
	CY_NODISCARD Matrix3 operator - ( T value ) const { Matrix3 r=*this; r.cell[0]-=value; r.cell[4]-=value; r.cell[8]-=value; return r; }	//!< subtract a value times identity matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Assignment Operators

	Matrix3 const & operator *= ( Matrix3 const &right ) { *this = operator*(right);                       return *this; }	//!< multiply a matrix with another matrix and modify this matrix
	Matrix3 const & operator += ( Matrix3 const &right ) { _CY_FOR_9i( cell[i] += right.cell[i] );         return *this; }	//!< add two Matrices modify this
	Matrix3 const & operator -= ( Matrix3 const &right ) { _CY_FOR_9i( cell[i] -= right.cell[i] );         return *this; }	//!< subtract one Matrix3 from another matrix and modify this matrix
	Matrix3 const & operator *= ( T       const  value ) { _CY_FOR_9i( cell[i] *= value         );         return *this; }	//!< multiply a matrix with a value modify this matrix
	Matrix3 const & operator /= ( T       const  value ) { _CY_FOR_9i( cell[i] /= value         );         return *this; }	//!< divide the matrix by a value modify the this matrix
	Matrix3 const & operator += ( T       const  value ) { cell[0]+=value; cell[4]+=value; cell[8]+=value; return *this; }	//!< add a value times identity matrix
	Matrix3 const & operator -= ( T       const  value ) { cell[0]-=value; cell[4]-=value; cell[8]-=value; return *this; }	//!< subtract a value times identity matrix

	//////////////////////////////////////////////////////////////////////////
	//!@name Other Methods

	//! Adds a diagonal matrix to this matrix and returns the result.
	CY_NODISCARD Matrix3 AddDiagonal( T xx, T yy, T zz ) const
	{
		Matrix3 m;
		m.cell[0] = cell[0] + xx;
		m.cell[1] = cell[1];
		m.cell[2] = cell[2];
		m.cell[3] = cell[3];
		m.cell[4] = cell[4] + yy;
		m.cell[5] = cell[5];
		m.cell[6] = cell[6];
		m.cell[7] = cell[7];
		m.cell[8] = cell[8] + zz;
		return m;
	}
	CY_NODISCARD Matrix3 AddDiagonal( Vec3<T> const &diag ) const { return AddDiagonal(diag.x,diag.y,diag.z); }	//!< Adds a diagonal matrix to this matrix and returns the result.
	CY_NODISCARD Matrix3 AddIdentity( T scale=T(1) )        const { return AddDiagonal( scale, scale, scale ); }	//!< Adds a scaled identity matrix to this matrix and returns the result.

	void Transpose()	//!< Transpose this matrix
	{
		for ( int i = 1; i < 3; ++i ) {
			for ( int j = 0; j < i; j++) {
				T temp = cell[i * 3 + j];
				cell[i * 3 + j] = cell[j * 3 + i];
				cell[j * 3 + i] = temp;
			}
		}
	}
	CY_NODISCARD Matrix3 GetTranspose() const	//!< Return the transpose of this matrix
	{
		Matrix3 m;
		m.cell[0] = cell[0];   m.cell[1] = cell[3];   m.cell[2] = cell[6];
		m.cell[3] = cell[1];   m.cell[4] = cell[4];   m.cell[5] = cell[7];
		m.cell[6] = cell[2];   m.cell[7] = cell[5];   m.cell[8] = cell[8];
		return m;
	}

	//! Multiply the give vector with the transpose of the matrix
	CY_NODISCARD Vec3<T> TransposeMult( Vec3<T> const &p ) const
	{
		return Vec3<T>( p.x*cell[0] + p.y*cell[1] + p.z*cell[2], 
		                p.x*cell[3] + p.y*cell[4] + p.z*cell[5],
		                p.x*cell[6] + p.y*cell[7] + p.z*cell[8] );
	}

	CY_NODISCARD Matrix3 TransposeMult( Matrix3 const & right ) const //!< Multiply a matrix by the transpose of this one (i.e. this^T * right).
	{
		Matrix3 r;
		for ( int i=0, k=0; i<3; ++i ) {
			for ( int j=0; j<3; ++j, ++k ) {
				r.cell[k] = Column(j).Dot( right.Column(i) );
			}
		}
		return r;
	}
	CY_NODISCARD Matrix3 MultTranspose( Matrix3 const & right ) const //!< Multiply the transpose of a matrix by this one (i.e. this * right^T).
	{
		_CY_INIT_MATRIX3_VECTORIZATION;
		Matrix3 rm;
		for ( int i=0; i<3; ++i ) {
			T a[4], b[4], c[4], d[4], r[4];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) a[j] = cell[  j] * right.cell[i  ];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) b[j] = cell[3+j] * right.cell[i+3];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) c[j] = cell_6[j] * right.cell[i+6];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) d[j] = a[j] + b[j];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) r[j] = d[j] + c[j];
			MemCopy( rm.cell+i, r, 3 );
		}
		return rm;
	}

	CY_NODISCARD Matrix3 TransposeMultSelf() const { return TransposeMult(*this); } //!< Multiply the transpose of this matrix with itself (i.e. this^T * this).
	CY_NODISCARD Matrix3 MultSelfTranspose() const { return MultTranspose(*this); } //!< Multiply the matrix with its transpose (i.e. this * this^T).

	CY_NODISCARD T GetTrace() const { return cell[0]+cell[4]+cell[8]; }

	CY_NODISCARD T GetDeterminant() const {	//!< Get the determinant of this matrix
		// 0 (4 8 - 5 7) + 1 (5 6 - 3 8) + 2 (3 7 - 4 6)
		return cell[0] * ( cell[4] * cell[8] - cell[5] * cell[7] ) + 
		       cell[1] * ( cell[5] * cell[6] - cell[3] * cell[8] ) + 
		       cell[2] * ( cell[3] * cell[7] - cell[4] * cell[6] );
	}

	void Invert() { *this = GetInverse(); }	//!< Invert this matrix
	CY_NODISCARD Matrix3 GetInverse() const	//!< Get the inverse of this matrix
	{
		//  ( 4 8 - 5 7    5 6 - 3 8    3 7 - 4 6 ) 
		//  ( 2 7 - 1 8    0 8 - 2 6    1 6 - 0 7 )  / det
		//  ( 1 5 - 2 4    2 3 - 0 5    0 4 - 1 3 ) 

		Matrix3 inverse;

		inverse.cell[0] = (cell[4]*cell[8] - cell[5]*cell[7]);
		inverse.cell[1] = (cell[2]*cell[7] - cell[1]*cell[8]);
		inverse.cell[2] = (cell[1]*cell[5] - cell[2]*cell[4]);

		inverse.cell[3] = (cell[5]*cell[6] - cell[3]*cell[8]);
		inverse.cell[4] = (cell[0]*cell[8] - cell[2]*cell[6]);
		inverse.cell[5] = (cell[2]*cell[3] - cell[0]*cell[5]);

		inverse.cell[6] = (cell[3]*cell[7] - cell[4]*cell[6]);
		inverse.cell[7] = (cell[1]*cell[6] - cell[0]*cell[7]);
		inverse.cell[8] = (cell[0]*cell[4] - cell[1]*cell[3]);

		T det = cell[0] * inverse.cell[0] + cell[1] * inverse.cell[3] + cell[2] * inverse.cell[6];
		return inverse / det;
	}

	//! Removes the scale component of the matrix by normalizing each column.
	//! The resulting matrix can contain shear, if it originally contained non-uniform scale and rotation.
	void Normalize() { Column(0).Normalize(); Column(1).Normalize(); Column(2).Normalize(); }

	//! Orthogonalizes the matrix and removes the scale component, preserving the x direction
	void OrthogonalizeX()
	{
		Column(0).Normalize();
		Column(1) -= Column(0) * (Column(1) % Column(0));
		Column(1).Normalize();
		Column(2) -= Column(0) * (Column(2) % Column(0));
		Column(2) -= Column(1) * (Column(2) % Column(1));
		Column(2).Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the y direction
	void OrthogonalizeY()
	{
		Column(1).Normalize();
		Column(0) -= Column(1) * (Column(0) % Column(1));
		Column(0).Normalize();
		Column(2) -= Column(1) * (Column(2) % Column(1));
		Column(2) -= Column(0) * (Column(2) % Column(0));
		Column(2).Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the z direction
	void OrthogonalizeZ()
	{
		Column(2).Normalize();
		Column(0) -= Column(2) * (Column(0) % Column(2));
		Column(0).Normalize();
		Column(1) -= Column(2) * (Column(1) % Column(2));
		Column(1) -= Column(0) * (Column(1) % Column(0));
		Column(1).Normalize();
	}


	//! Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const
	{
		return std::abs(cell[0]-T(1)) < tollerance && std::abs(cell[1])      < tollerance && std::abs(cell[2])      < tollerance && 
		       std::abs(cell[3])      < tollerance && std::abs(cell[4]-T(1)) < tollerance && std::abs(cell[5])      < tollerance &&
		       std::abs(cell[6])      < tollerance && std::abs(cell[7])      < tollerance && std::abs(cell[8]-T(1)) < tollerance;
	}

	//! Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const { return std::abs(cell[1] - cell[3]) < tollerance && std::abs(cell[2] - cell[6]) < tollerance && std::abs(cell[5] - cell[7]) < tollerance; }

	//! Returns if the matrix is diagonal.
	bool IsDiagonal( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const { return std::abs(cell[1]) + std::abs(cell[2]) + std::abs(cell[3]) + std::abs(cell[5]) + std::abs(cell[6]) + std::abs(cell[7]) < tollerance*6; }

	//! Returns the eigenvalues of the matrix.
	//! The eigenvalues are ordered, such that the first one is the largest.
	//! The given tollerance value is used for checking whether the matrix is diagonal.
	CY_NODISCARD Vec3<T> GetEigenvalues( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const
	{
		Vec3<T> lambda;
		if ( IsDiagonal(tollerance) ) {
			lambda = GetDiagonal();	// diagonal matrix, so the eigenvalues are the same as the diagonal elements.
		} else {
			T q = GetTrace() / 3;
			Matrix3 m = AddIdentity(-q);
			Matrix3 m2 = m * m;
			T p = Sqrt( m2.GetTrace() / T(6) );
			Matrix3 B = m / p;
			T d_2 = B.GetDeterminant() * T(0.5);
			T a = d_2 < T(-1) ? Pi<T>()/3 : ( d_2 > T(1) ? T(0) : ACosSafe<T>(d_2)/3 );	// only guaranteed to work for symmetric matrices
			Vec3<T> b;
			b.x = 2*std::cos( a );
			b.y = 2*std::cos( a + Pi<T>()*(T(2)/T(3)) );
			b.z = 2*std::cos( a + Pi<T>()*(T(4)/T(3)) );
			lambda = b*p + q;
		}
		return lambda;
	}

	//! Returns the eigenvalues and sets the given vector as the eigenvectors of the matrix.
	//! The eigenvalues are ordered, such that the first one is the largest.
	//! The given tollerance is used for checking whether the eigenvalues are the same.
	Vec3<T> GetEigenvectors( Vec3<T> &evec0, Vec3<T> &evec1, Vec3<T> &evec2, T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const
	{
		static auto setVectors = [&]( Vec3<T> &e0, Vec3<T> &e1, Vec3<T> &e2, Matrix3 const &v20, Matrix3 const &v12, Matrix3 const &v01 ) {
			int i = 0;
			Matrix3 v2 = v20 * v12;
			e2 = v2.Column(0) + v2.Column(1) + v2.Column(2);
			if ( (e2 ^ v01.Column(0)).LengthSquared() < tollerance ) {
				e0 = v01.Column(1);
				e1 = v01.Column(2);
			} else {
				e0 = v01.Column(0);
				e1 = v01.Column( ( (e2 ^ v01.Column(1)).LengthSquared() < tollerance ) ? 2 : 1 );
			}
		};

		Vec3<T> lambda = GetEigenvalues(tollerance);
		bool same01 = std::abs(lambda.x - lambda.y) < tollerance;
		bool same12 = std::abs(lambda.y - lambda.z) < tollerance;
		bool same02 = std::abs(lambda.x - lambda.z) < tollerance;
		if ( same01 & same12 ) {
			// All eigenvalues are the same, so, assuming that the matrix is normal, it must be scaled identity.
			evec0 = Column(0);
			evec1 = Column(1);
			evec2 = Column(2);
		} else {
			Matrix3 v12 = AddIdentity( -lambda.x );
			Matrix3 v20 = AddIdentity( -lambda.y );
			Matrix3 v01 = AddIdentity( -lambda.z );
			char same = char(same01) | char(same12)*char(2) | char(same02)*char(3);	// only one of them can be true here
			switch (same) {
				default: case 0: {
					Matrix3 v0 = v01 * v20;
					Matrix3 v1 = v12 * v01;
					Matrix3 v2 = v20 * v12;
					evec0 = v0.Column(0) + v0.Column(1) + v0.Column(2);
					evec1 = v1.Column(0) + v1.Column(1) + v1.Column(2);
					evec2 = v2.Column(0) + v2.Column(1) + v2.Column(2);
				}
				break;
				case 1: setVectors( evec0, evec1, evec2, v20, v12, v01 ); break;
				case 2: setVectors( evec1, evec2, evec0, v01, v20, v12 ); break;
				case 3: setVectors( evec0, evec2, evec1, v12, v01, v20 ); break;
			}
		}
		return lambda;
	}

	//! Singular value decomposition (SVD)
	//! Returns the SVD of the matrix, where U and V are orthogonal matrices and 
	//! S is the diagonal elements of a diagonal matrix (including zeros),
	//! such that this matrix A = U S V^T.
	//! The given tollerance is used for checking whether the eigenvalues are the same.
	void SingularValueDecomposition( Matrix3<T> &U, Vec3<T> &S, Matrix3<T> &V, T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) )
	{
		Matrix3 AAT = MultSelfTranspose();
		Vec3<T> lambda = AAT.GetEigenvectors( U.Column(0), U.Column(1), U.Column(2), tollerance );
		S = (lambda.Abs()).Sqrt();
		U.Normalize();
		Matrix3 ATA = TransposeMultSelf();
		AAT.GetEigenvectors( V.Column(0), V.Column(1), V.Column(2), tollerance );
		V.Normalize();
	}

	//////////////////////////////////////////////////////////////////////////
	//!@name Static Methods

	//! Returns an identity matrix
	CY_NODISCARD static Matrix3 Identity() { T c[] = { 1,0,0, 0,1,0, 0,0,1 }; return Matrix3(c); }
	//! Returns a view matrix using position, target and approximate up vector
	CY_NODISCARD static Matrix3 View( Vec3<T> const &target, Vec3<T> const &up ) { Matrix3 m; m.SetView(target,up); return m; }
	//! Returns a rotation matrix around x axis by angle in radians
	CY_NODISCARD static Matrix3 RotationX( T angle ) { Matrix3 m; m.SetRotationX(angle); return m; }
	//! Returns a rotation matrix around y axis by angle in radians
	CY_NODISCARD static Matrix3 RotationY( T angle ) { Matrix3 m; m.SetRotationY(angle); return m; }
	//! Returns a rotation matrix around z axis by angle in radians
	CY_NODISCARD static Matrix3 RotationZ( T angle ) { Matrix3 m; m.SetRotationZ(angle); return m; }
	//! Returns a rotation matrix around x, y, and then z axes by angle in radians (Rz * Ry * Rx)
	CY_NODISCARD static Matrix3 RotationXYZ( T angleX, T angleY, T angleZ ) { Matrix3 m; m.SetRotationXYZ(angleX,angleY,angleZ); return m; }
	//! Returns a rotation matrix around z, y, and then x axes by angle in radians (Rx * Ry * Rz)
	CY_NODISCARD static Matrix3 RotationZYX( T angleX, T angleY, T angleZ ) { Matrix3 m; m.SetRotationZYX(angleX,angleY,angleZ); return m; }
	//! Returns a rotation matrix about the given axis by angle in radians
	CY_NODISCARD static Matrix3 Rotation( Vec3<T> const &axis, T angle ) { Matrix3 m; m.SetRotation(axis,angle); return m; }
	//! Returns a rotation matrix about the given axis by cos and sin of the rotation angle
	CY_NODISCARD static Matrix3 Rotation( Vec3<T> const &axis, T cosAngle, T sinAngle ) { Matrix3 m; m.SetRotation(axis,cosAngle,sinAngle); return m; }
	//! Returns a rotation matrix that sets [from] unit vector to [to] unit vector
	CY_NODISCARD static Matrix3 Rotation( Vec3<T> const &from, Vec3<T> const &to ) { Matrix3 m; m.SetRotation(from,to); return m; }
	//! Returns a uniform scale matrix
	CY_NODISCARD static Matrix3 Scale( T uniformScale ) { Matrix3 m; m.SetScale(uniformScale); return m; }
	//! Returns a scale matrix
	CY_NODISCARD static Matrix3 Scale( T scaleX, T scaleY, T scaleZ ) { Matrix3 m; m.SetScale(scaleX,scaleY,scaleZ); return m; }
	//! Returns a scale matrix
	CY_NODISCARD static Matrix3 Scale( Vec3<T> const &scale ) { Matrix3 m; m.SetScale(scale); return m; }
	//! Returns the tensor product (outer product) matrix of two vectors
	CY_NODISCARD static Matrix3 TensorProduct( Vec3<T> const &v0, Vec3<T> const &v1 ) { Matrix3 m; m.SetTensorProduct(v0,v1); return m; }
	//! Returns the matrix representation of cross product ( a x b )
	CY_NODISCARD static Matrix3 MatrixCrossProd( Vec3<T> const &a ) { Matrix3 m; m.SetCrossProd(a); return m; }

	//////////////////////////////////////////////////////////////////////////
};

//-------------------------------------------------------------------------------

#ifdef CY_NONVECTORIZED_MATRIX3
# define _CY_INIT_MATRIX34_VECTORIZATION const int N = 3; T const *cell_9 = cell + 9;
#else
# define _CY_INIT_MATRIX34_VECTORIZATION const int N = 4; T cell_9[4] = { cell[9], cell[10], cell[11], cell[11] };
#endif

//-------------------------------------------------------------------------------

//! 3x4 matrix class.
//!
//! Its data stores 12-value array of column-major matrix elements.
//! I chose column-major format to be compatible with OpenGL
//! You can use Matrix34 with Vec3<T> and Vec4<T>
//! to transform 3D and 4D points.

template <typename T>
class Matrix34
{
	CY_NODISCARD friend Matrix34 operator * ( T value, Matrix34 const &right ) { Matrix34 r; _CY_FOR_12i( r.cell[i]=value*right.cell[i] ); return r; }	//!< multiply matrix by a value
	CY_NODISCARD friend Matrix34 operator + ( T value, Matrix34 const &right ) { Matrix34 r= right; r.cell[0]+=value; r.cell[4]+=value; r.cell[8]+=value; return r; }	//!< add a value times identity matrix to a matrix
	CY_NODISCARD friend Matrix34 operator - ( T value, Matrix34 const &right ) { Matrix34 r=-right; r.cell[0]+=value; r.cell[4]+=value; r.cell[8]+=value; return r; }	//!< subtract a matrix from a value times identity matrix
	CY_NODISCARD friend Matrix34 Inverse( Matrix34 const &m ) { return m.GetInverse(); }	//!< return the inverse of the matrix

public:

	//! Elements of the matrix are column-major: \n
	//! | 0   3   6   9 | \n
	//! | 1   4   7  10 | \n
	//! | 2   5   8  11 | \n
#ifdef __cpp_unrestricted_unions
	union {
		T       cell[12];
		Vec3<T> column[4];	// column vectors
	};
#else
	T cell[12];
#endif


	//////////////////////////////////////////////////////////////////////////
	//!@name Constructors

	Matrix34() CY_CLASS_FUNCTION_DEFAULT												//!< Default constructor
	template <typename TT> explicit Matrix34<T>( Matrix34<TT> const &matrix ) { MemConvert(cell,matrix.cell,12); }	//!< Copy constructor for different types
	explicit Matrix34( T const * restrict values ) { Set(values); }						//!< Initialize the matrix using an array of 9 values
	explicit Matrix34( T v )                       { SetScale(v); }						//!< Initialize the matrix as identity scaled by v
	explicit Matrix34( Vec3<T> const &x, Vec3<T> const &y, Vec3<T> const &z, Vec3<T> const &pos ) { Set(x,y,z,pos); }	//!< Initialize the matrix using x,y,z vectors and coordinate center
	explicit Matrix34( Matrix3<T> const &m ) { MemCopy(cell,m.cell,9); Column(3).Zero(); }
	explicit Matrix34( Matrix3<T> const &m, Vec3<T> const &pos ) { MemCopy(cell,m.cell,9); Column(3)=pos; }
	explicit Matrix34( Matrix2<T> const &m ) { Column(0)=Vec3<T>(m.Column(0)); Column(1)=Vec3<T>(m.Column(1)); Column(2).Set(0,0,1); Column(3).Zero(); }
	explicit Matrix34( Matrix4<T> const &m );

	//! Constructor using row-major order for initialization
	Matrix34( T c00, T c01, T c02, T c03,
		      T c10, T c11, T c12, T c13,
		      T c20, T c21, T c22, T c23 )
	{
		cell[ 0] = c00;   cell[ 3] = c01;   cell[ 6] = c02;   cell[ 9] = c03;
		cell[ 1] = c10;   cell[ 4] = c11;   cell[ 7] = c12;   cell[10] = c13;
		cell[ 2] = c20;   cell[ 5] = c21;   cell[ 8] = c22;   cell[11] = c23;
	}

	//////////////////////////////////////////////////////////////////////////
	//!@name Set & Get Methods

	void Zero    ()       { MemClear(cell,12); }																					//!< Set all the values as zero
	bool IsZero  () const { return Column(0).IsZero  () && Column(1).IsZero  () && Column(2).IsZero  () && Column(3).IsZero  (); }	//!< Returns true if the matrix is exactly zero
	bool IsFinite() const { return Column(0).IsFinite() && Column(1).IsFinite() && Column(2).IsFinite() && Column(3).IsFinite(); }	//!< Returns true if all components are finite real numbers.
	void Get( T       * restrict values ) { MemCopy(values,cell,12); }																//!< Copies the matrix cell to the given values array of size 12
	void Set( T const * restrict values ) { MemCopy(cell,values,12); }																//!< Set Matrix using an array of 12 values
	void Set( Vec3<T> const &x, Vec3<T> const &y, Vec3<T> const &z, Vec3<T> const &pos ) { Column(0)=x; Column(1)=y; Column(2)=z; Column(3)=pos; }	//!< Set matrix using x,y,z vectors and coordinate center
	void SetIdentity()    { SetScale(T(1)); }																						//!< Converts the matrix to an identity matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Affine transformations

	//! Sets a uniform scale matrix
	void SetScale( T uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	//! Sets a scale matrix
	void SetScale( T scaleX, T scaleY, T scaleZ )
	{
		Zero();
		cell[0] = scaleX;
		cell[4] = scaleY;
		cell[8] = scaleZ;
	}
	//! Sets a scale matrix
	void SetScale( Vec3<T> const &scale ) { SetScale(scale.x,scale.y,scale.z); }
	//! Set as rotation matrix around x axis
	void SetRotationX( T angle ) { SetRotationX( std::sin(angle), std::cos(angle) ); }
	//! Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( T sinAngle, T cosAngle )
	{
		cell[0] = T(1);   cell[1] = T(0);        cell[2] = T(0); 
		cell[3] = T(0);   cell[4] = cosAngle;    cell[5] = sinAngle;
		cell[6] = T(0);   cell[7] = -sinAngle;   cell[8] = cosAngle;
		MemClear(cell+9,3);
	}
	//! Set as rotation matrix around y axis
	void SetRotationY( T angle ) { SetRotationY( std::sin(angle), std::cos(angle) ); }
	//! Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( T sinAngle, T cosAngle )
	{
		cell[0] = cosAngle;   cell[1] = T(0);   cell[2] = -sinAngle;  
		cell[3] = T(0);       cell[4] = T(1);   cell[5] = T(0);  
		cell[6] = sinAngle;   cell[7] = T(0);   cell[8] = cosAngle; 
		MemClear(cell+9,3);
	}
	//! Set as rotation matrix around z axis
	void SetRotationZ( T angle ) { SetRotationZ( std::sin(angle), std::cos(angle) ); }
	//! Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( T sinAngle, T cosAngle )
	{
		cell[0] =  cosAngle;   cell[1] = sinAngle;   cell[2] = T(0);
		cell[3] = -sinAngle;   cell[4] = cosAngle;   cell[5] = T(0); 
		cell[6] =  T(0);       cell[7] = T(0);       cell[8] = T(1);
		MemClear(cell+9,3);
	}
	//! Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
	void SetRotationXYZ( T angleX, T angleY, T angleZ )
	{
		T sx = std::sin(angleX);
		T cx = std::cos(angleX);
		T sy = std::sin(angleY);
		T cy = std::cos(angleY);
		T sz = std::sin(angleZ);
		T cz = std::cos(angleZ);
		cell[0] = cy*cz; 		      cell[1] = cy*sz; 			    cell[2] =-sy;   
		cell[3] = cz*sx*sy - cx*sz;   cell[4] = cx*cz + sx*sy*sz;   cell[5] = cy*sx; 
		cell[6] = cx*cz*sy + sx*sz;   cell[7] =-cz*sx + cx*sy*sz;   cell[8] = cx*cy;
		MemClear(cell+9,3);
	}
	//! Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
	void SetRotationZYX( T angleX, T angleY, T angleZ )
	{
		T sx = std::sin(angleX);
		T cx = std::cos(angleX);
		T sy = std::sin(angleY);
		T cy = std::cos(angleY);
		T sz = std::sin(angleZ);
		T cz = std::cos(angleZ);
		cell[0] = cy*cz;   cell[1] = cx*sz + sx*sy*cz;   cell[2] = sx*sz - cx*sy*cz;            
		cell[3] =-cy*sz;   cell[4] = cx*cz - sx*sy*sz;   cell[5] = sx*cz + cx*sy*sz;            
		cell[6] = sy;      cell[7] =-sx*cy;			     cell[8] = cx*cy;
		MemClear(cell+9,3);
	}
	//! Set a rotation matrix about the given axis by angle
	void SetRotation( Vec3<T> const &axis, T angle ) { SetRotation(axis,std::sin(angle),std::cos(angle)); }
	//! Set a rotation matrix about the given axis by cos and sin of angle
	void SetRotation( Vec3<T> const &axis, T sinAngle, T cosAngle )
	{
		T t = T(1) - cosAngle;
		Vec3<T> a = t * axis;
		T txy = a.x * axis.y;
		T txz = a.x * axis.z;
		T tyz = a.y * axis.z;
		Vec3<T> s = sinAngle * axis;
		cell[ 0] = a.x * axis.x + cosAngle;   cell[ 1] = txy + s.z;                 cell[ 2] = txz - s.y;
		cell[ 3] = txy - s.z;                 cell[ 4] = a.y * axis.y + cosAngle;   cell[ 5] = tyz + s.x;
		cell[ 6] = txz + s.y;                 cell[ 7] = tyz - s.x;                 cell[ 8] = a.z * axis.z + cosAngle;
		MemClear(cell+9,3);
	}
	//! Set a rotation matrix that sets [from] unit vector to [to] unit vector
	void SetRotation( Vec3<T> const &from, Vec3<T> const &to )
	{
		T c = from.Dot(to);
		if ( c > T(0.9999999) ) SetIdentity();
		else {
			T s = Sqrt(T(1) - c*c);
			Vec3<T> axis = from.Cross(to).GetNormalized();
			SetRotation(axis, s, c);
		}
	}
	//! Sets a translation matrix with no rotation or scale
	void SetTranslation( Vec3<T> const &move ) { Column(0).Set(1,0,0); Column(1).Set(0,1,0); Column(2).Set(0,0,1); Column(3)=move; }
	//! Adds a translation to the matrix
	void AddTranslation( Vec3<T> const &move ) { Column(3) += move; }
	//! Sets the translation component of the matrix
	void SetTranslationComponent( Vec3<T> const &move ) { Column(3) = move; }
	//! Sets the translation component of the matrix to zero
	void SetNoTranslation() { Column(3).Set(0,0,0); }
	//! Set view matrix using position, target and approximate up vector
	void SetView( Vec3<T> const &pos, Vec3<T> const &target, Vec3<T> const &up )
	{
		Vec3<T> f = target - pos;
		f.Normalize();
		Vec3<T> s = f.Cross(up);
		s.Normalize();
		Vec3<T> u = s.Cross(f);
		cell[ 0]=s.x; cell[ 1]=u.x; cell[ 2]=-f.x;
		cell[ 3]=s.y; cell[ 4]=u.y; cell[ 5]=-f.y;
		cell[ 6]=s.z; cell[ 7]=u.z; cell[ 8]=-f.z;
		cell[ 9]= -s % pos;
		cell[10]= -u % pos;
		cell[11]=  f % pos;
	}
	//! Sets a Cartesian coordinate frame using the given x direction, an approximate y direction, and a translation. x must be a unit vector.
	void SetCartesianFrameXY( Vec3<T> const &x, Vec3<T> const &y_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> z = x.Cross(y_approx); z.Normalize(); Vec3<T> y=z.Cross(x); Set(x,y,z,trans); }
	//! Sets a Cartesian coordinate frame using the given x direction, an approximate z direction, and a translation. x must be a unit vector.
	void SetCartesianFrameXZ( Vec3<T> const &x, Vec3<T> const &z_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> y = z_approx.Cross(x); y.Normalize(); Vec3<T> z=x.Cross(y); Set(x,y,z,trans); }
	//! Sets a Cartesian coordinate frame using the given y direction, an approximate x direction, and a translation. y must be a unit vector.
	void SetCartesianFrameYX( Vec3<T> const &y, Vec3<T> const &x_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> z = x_approx.Cross(y); z.Normalize(); Vec3<T> x=y.Cross(z); Set(x,y,z,trans); }
	//! Sets a Cartesian coordinate frame using the given y direction, an approximate z direction, and a translation. y must be a unit vector.
	void SetCartesianFrameYZ( Vec3<T> const &y, Vec3<T> const &z_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> x = y.Cross(z_approx); x.Normalize(); Vec3<T> z=x.Cross(y); Set(x,y,z,trans); }
	//! Sets a Cartesian coordinate frame using the given z direction, an approximate x direction, and a translation. z must be a unit vector.
	void SetCartesianFrameZX( Vec3<T> const &z, Vec3<T> const &x_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> y = z.Cross(x_approx); y.Normalize(); Vec3<T> x=y.Cross(z); Set(x,y,z,trans); }
	//! Sets a Cartesian coordinate frame using the given z direction, an approximate y direction, and a translation. z must be a unit vector.
	void SetCartesianFrameZY( Vec3<T> const &z, Vec3<T> const &y_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> x = y_approx.Cross(z); x.Normalize(); Vec3<T> y=z.Cross(x); Set(x,y,z,trans); }


	//////////////////////////////////////////////////////////////////////////
	//!@name Set Row, Column, or Diagonal

	void SetRow     ( int ri, T x, T y, T z, T w ) { cell[ri]=x; cell[ri+3]=y; cell[ri+6]=z; cell[ri+9]=w; }	//!< Sets a row of the matrix
	void SetRow     ( int ri, Vec4<T> const &v )   { SetRow(ri,v.x,v.y,v.z,v.w); }								//!< Sets a row of the matrix
	void SetColumn  ( int ci, T x, T y, T z )      { Column(ci).Set(x,y,z); }									//!< Sets a column of the matrix
	void SetColumn  ( int ci, Vec3<T> const &v )   { Column(ci)=v; }											//!< Sets a column of the matrix
	void SetDiagonal( T xx, T yy, T zz )           { cell[0]=xx; cell[4]=yy; cell[8]=zz; }						//!< Sets the diagonal values of the matrix
	void SetDiagonal( Vec3<T> const &p )           { SetDiagonal( p.x, p.y, p.z ); }							//!< Sets the diagonal values of the matrix
	void SetDiagonal( T const * restrict values )  { SetDiagonal(values[0],values[1],values[2]); }				//!< Sets the diagonal values of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Row, Column, or Diagonal

#ifdef __cpp_unrestricted_unions
	CY_NODISCARD Vec3<T>       * Columns()               { return column; }
	CY_NODISCARD Vec3<T> const * Columns()         const { return column; }
	CY_NODISCARD Vec3<T>       & Column ( int ci )       { return column[ci]; }
	CY_NODISCARD Vec3<T> const & Column ( int ci ) const { return column[ci]; }
#else
	CY_NODISCARD Vec3<T>       * Columns()               { return ((Vec3<T>*)cell); }
	CY_NODISCARD Vec3<T> const * Columns()         const { return ((Vec3<T>*)cell); }
	CY_NODISCARD Vec3<T>       & Column ( int ci )       { return Columns()[ci]; }
	CY_NODISCARD Vec3<T> const & Column ( int ci ) const { return Columns()[ci]; }
#endif
	CY_NODISCARD Vec4<T>         GetRow ( int ri ) const { return Vec4<T>( cell[ri], cell[ri+3], cell[ri+6], cell[ri+9] ); }	//!< Returns a row of the matrix
	CY_NODISCARD Vec3<T>         GetDiagonal()     const { return Vec3<T>( cell[0], cell[4], cell[8] ); }							//!< Returns the diagonal of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Sub-matrix cell

	CY_NODISCARD Matrix3<T> GetSubMatrix3 () const { Matrix3<T> m; MemCopy(m.cell,cell,9); return m; }	//!< Returns the 3x3 portion of the matrix
	CY_NODISCARD Matrix2<T> GetSubMatrix2 () const { Matrix2<T> m; MemCopy(m.cell,cell,2); MemCopy(m.cell+2,cell+3,2); return m; }	//!< Returns the 2x2 portion of the matrix
	CY_NODISCARD Vec3<T>    GetTranslation() const { return Column(3); }								//!< Returns the translation component of the matrix
	CY_NODISCARD Matrix3<T> GetRotation   () const { Matrix3<T> m(*this); return m.GetRotation(); }		//!< Returns the rotation portion of the transformation
	CY_NODISCARD Matrix3<T> GetScale      () const { Matrix3<T> m(*this); return m.GetScale   (); }		//!< Returns the scale portion of the transformation.

	//! Returns the average scale factor of the 3 by 3 sub-matrix
	T GetAvrgScale() const 
	{
		T det = cell[0] * ( cell[4] * cell[8] - cell[5] * cell[7] ) + 
		        cell[1] * ( cell[5] * cell[6] - cell[3] * cell[8] ) + 
		        cell[2] * ( cell[3] * cell[7] - cell[4] * cell[6] );
		T s = std::pow( std::abs(det), T(1)/T(3) );
		return det >= 0 ? s : -s;
	}

	void GetComponents( Matrix3<T> &scale, Matrix3<T> &rotation, Vec3<T> &translation ) const { Matrix3<T> m(*this); m.GetComponents(scale,rotation); translation=GetTranslation(); }	//!< Returns separate transformation components

	//////////////////////////////////////////////////////////////////////////
	//!@name Comparison Operators

	CY_NODISCARD bool operator == ( Matrix34 const &right ) const { _CY_FOR_12i( if ( cell[i] != right.cell[i] ) return false ); return true;  } //!< compare equal
	CY_NODISCARD bool operator != ( Matrix34 const &right ) const { _CY_FOR_12i( if ( cell[i] != right.cell[i] ) return true  ); return false; } //!< compare not equal


	//////////////////////////////////////////////////////////////////////////
	//!@name Access Operators

	CY_NODISCARD T&        operator () ( int ri, int ci )       { assert( ri>=0 && ri<3 && ci>=0 && ci<4 ); return cell[ ci*3 + ri ]; }	//!< subscript operator
	CY_NODISCARD T const & operator () ( int ri, int ci ) const { assert( ri>=0 && ri<3 && ci>=0 && ci<4 ); return cell[ ci*3 + ri ]; }	//!< constant subscript operator
	CY_NODISCARD T&        operator [] ( int i )                { assert( i>=0 && i<12 ); return cell[i]; }								//!< subscript operator
	CY_NODISCARD T const & operator [] ( int i )          const { assert( i>=0 && i<12 ); return cell[i]; }								//!< constant subscript operator


	//////////////////////////////////////////////////////////////////////////
	//!@name Unary and Binary Operators

	// Unary operators
	CY_NODISCARD Matrix34 operator - () const { Matrix34 r; _CY_FOR_12i( r.cell[i] = -cell[i] ); return r; }	//!< negative matrix

	// Binary operators
	CY_NODISCARD Matrix34 operator * ( T        const  value ) const { Matrix34 r; _CY_FOR_12i( r.cell[i] = cell[i] * value         ); return r; }	//!< multiply matrix by a value
	CY_NODISCARD Matrix34 operator / ( T        const  value ) const { Matrix34 r; _CY_FOR_12i( r.cell[i] = cell[i] / value         ); return r; }	//!< divide matrix by a value
	CY_NODISCARD Matrix34 operator + ( Matrix34 const &right ) const { Matrix34 r; _CY_FOR_12i( r.cell[i] = cell[i] + right.cell[i] ); return r; }	//!< add two Matrices
	CY_NODISCARD Matrix34 operator - ( Matrix34 const &right ) const { Matrix34 r; _CY_FOR_12i( r.cell[i] = cell[i] - right.cell[i] ); return r; }	//!< subtract one Matrix4 from another
	CY_NODISCARD Matrix34 operator * ( Matrix34 const &right ) const	//!< multiply a matrix with another
	{
		_CY_INIT_MATRIX34_VECTORIZATION;
		Matrix34 rm;
		T *rd = rm.cell;
		for ( int i=0; i<12; i+=3, rd+=3 ) {
			T a[4], b[4], c[4], d[4], r[4];
			_CY_IVDEP_FOR ( int k=0; k<N; ++k ) a[k] = cell[  k] * right.cell[i  ];
			_CY_IVDEP_FOR ( int k=0; k<N; ++k ) b[k] = cell[3+k] * right.cell[i+1];
			_CY_IVDEP_FOR ( int k=0; k<N; ++k ) c[k] = cell[6+k] * right.cell[i+2];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) d[j] = a[j] + b[j];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) r[j] = d[j] + c[j];
			MemCopy( rd, r, 3 );
		}
		for ( int j=0; j<3; ++j ) rm.cell[9+j] += cell[9+j];
		return rm;
	}
	CY_NODISCARD Matrix34 operator * ( Matrix3<T> const &right ) const	//!< multiply a matrix with another
	{
		_CY_INIT_MATRIX34_VECTORIZATION;
		Matrix34 rm;
		T *rd = rm.cell;
		for ( int i=0; i<9; i+=3, rd+=3 ) {
			T a[4], b[4], c[4], d[4], r[4];
			_CY_IVDEP_FOR ( int k=0; k<N; ++k ) a[k] = cell[  k] * right.cell[i  ];
			_CY_IVDEP_FOR ( int k=0; k<N; ++k ) b[k] = cell[3+k] * right.cell[i+1];
			_CY_IVDEP_FOR ( int k=0; k<N; ++k ) c[k] = cell[6+k] * right.cell[i+2];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) d[j] = a[j] + b[j];
			_CY_IVDEP_FOR ( int j=0; j<N; ++j ) r[j] = d[j] + c[j];
			MemCopy( rd, r, 3 );
		}
		MemCopy(rm.cell+9,cell+9,3);
		return rm;
	}
	CY_NODISCARD Vec3<T> operator * ( Vec3<T> const &p ) const
	{
		_CY_INIT_MATRIX34_VECTORIZATION;
		//return Vec3<T>( p.x*cell[0] + p.y*cell[3] + p.z*cell[6] + cell[ 9], 
		//                p.x*cell[1] + p.y*cell[4] + p.z*cell[7] + cell[10],
		//                p.x*cell[2] + p.y*cell[5] + p.z*cell[8] + cell[11] );
		T a[4], b[4], c[4], d[4], e[4], r[4];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) a[i] = p.x * cell[  i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) b[i] = p.y * cell[3+i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) c[i] = p.z * cell[6+i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) d[i] = a[i] + b[i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) e[i] = c[i] + cell_9[i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) r[i] = d[i] + e[i];
		return Vec3<T>(r);
	}
	CY_NODISCARD Vec4<T> operator * ( Vec4<T> const &p ) const
	{
		_CY_INIT_MATRIX34_VECTORIZATION;
		//return Vec4<T>( p.x*cell[0] + p.y*cell[3] + p.z*cell[6] + p.w*cell[ 9],
		//                p.x*cell[1] + p.y*cell[4] + p.z*cell[7] + p.w*cell[10],
		//                p.x*cell[2] + p.y*cell[5] + p.z*cell[8] + p.w*cell[11],
		//                0           + 0           + 0           + p.w          );
		Vec4<T> r;
		T a[4], b[4], c[4], d[4], e[4], f[4];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) a[i] = p.x * cell[  i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) b[i] = p.y * cell[3+i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) c[i] = p.z * cell[6+i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) d[i] = p.w * cell_9[i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) e[i] = a[i] + b[i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) f[i] = c[i] + d[i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) r[i] = e[i] + f[3+i];
		r[3] = p.w;
		return r;
	}

	CY_NODISCARD Matrix34 operator + ( T value ) const { Matrix34 r=*this; r.cell[0]+=value; r.cell[4]+=value; r.cell[8]+=value; return r; }	//!< add a value times identity matrix
	CY_NODISCARD Matrix34 operator - ( T value ) const { Matrix34 r=*this; r.cell[0]-=value; r.cell[4]-=value; r.cell[8]-=value; return r; }	//!< subtract a value times identity matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name 3D Vector Transform Methods

	//! Transforms the vector by multiplying it with the matrix, ignoring the translation component.
	CY_NODISCARD Vec3<T> VectorTransform( Vec3<T> const &p ) const
	{
		_CY_INIT_MATRIX34_VECTORIZATION;
		//return Vec3<T>( p.x*cell[0] + p.y*cell[3] + p.z*cell[6], 
		//                p.x*cell[1] + p.y*cell[4] + p.z*cell[7],
		//                p.x*cell[2] + p.y*cell[5] + p.z*cell[8] );
		T a[4], b[4], c[4], d[4], r[4];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) a[i] = p.x * cell[   i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) b[i] = p.y * cell[ 3+i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) c[i] = p.z * cell[ 6+i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) d[i] = a[i] + b[i];
		_CY_IVDEP_FOR ( int i=0; i<N; ++i ) r[i] = d[i] + c[i];
		return Vec3<T>(r);
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Assignment Operators
	
	Matrix34 const & operator += ( Matrix34   const &right ) { _CY_FOR_12i( cell[i] += right.cell[i] );        return *this; }	//!< add two Matrices modify this
	Matrix34 const & operator -= ( Matrix34   const &right ) { _CY_FOR_12i( cell[i] -= right.cell[i] );        return *this; }	//!< subtract one Matrix4 from another matrix and modify this matrix
	Matrix34 const & operator *= ( Matrix34   const &right ) { *this = operator*(right);                       return *this; }	//!< multiply a matrix with another matrix and modify this matrix
	Matrix34 const & operator *= ( Matrix3<T> const &right ) { *this = operator*(right);                       return *this; }	//!< multiply a matrix with another matrix and modify this matrix
	Matrix34 const & operator *= ( T          const  value ) { _CY_FOR_12i( cell[i] *= value );                return *this; }	//!< multiply a matrix with a value modify this matrix
	Matrix34 const & operator /= ( T          const  value ) { _CY_FOR_12i( cell[i] /= value );                return *this; }	//!< divide the matrix by a value modify the this matrix
	Matrix34 const & operator += ( T          const  value ) { cell[0]+=value; cell[4]+=value; cell[8]+=value; return *this; }	//!< add a value times identity matrix
	Matrix34 const & operator -= ( T          const  value ) { cell[0]-=value; cell[4]-=value; cell[8]-=value; return *this; }	//!< subtract a value times identity matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Other Methods

	//! Transpose this matrix
	void Transpose()
	{
		for ( int i=1; i<3; ++i ) {
			for ( int j=0; j<i; j++ ) {
				T temp = cell[i * 3 + j];
				cell[i * 3 + j] = cell[j * 3 + i];
				cell[j * 3 + i] = temp;
			}
		}
		MemClear(cell+9,3);
	}

	CY_NODISCARD Matrix4<T> GetTranspose() const;	//!< Returns the transpose of this matrix

	//! Multiply the give vector with the transpose of the matrix
	CY_NODISCARD Vec4<T> TransposeMult( Vec3<T> const &p ) const
	{
		return Vec4<T>( p.x*cell[ 0] + p.y*cell[ 1] + p.z*cell[ 2], 
		                p.x*cell[ 3] + p.y*cell[ 4] + p.z*cell[ 5],
		                p.x*cell[ 6] + p.y*cell[ 7] + p.z*cell[ 8],
		                p.x*cell[ 9] + p.y*cell[10] + p.z*cell[11] + T(1) );
	}

	//! Multiply the give vector with the transpose of the matrix
	CY_NODISCARD Vec4<T> TransposeMult( Vec4<T> const &p ) const
	{
		return Vec4<T>( p.x*cell[ 0] + p.y*cell[ 1] + p.z*cell[ 2], 
		                p.x*cell[ 3] + p.y*cell[ 4] + p.z*cell[ 5],
		                p.x*cell[ 6] + p.y*cell[ 7] + p.z*cell[ 8],
		                p.x*cell[ 9] + p.y*cell[10] + p.z*cell[11] + p.w );
	}

	CY_NODISCARD T GetDeterminant() const	//!< Get the determinant of this matrix
	{
		// 0 (4 8 - 5 7) + 1 (5 6 - 3 8) + 2 (3 7 - 4 6)
		return cell[0] * ( cell[4] * cell[8] - cell[5] * cell[7] ) + 
		       cell[1] * ( cell[5] * cell[6] - cell[3] * cell[8] ) + 
		       cell[2] * ( cell[3] * cell[7] - cell[4] * cell[6] );
	}
	void Invert() { *this = GetInverse(); }		//!< Invert this matrix
	CY_NODISCARD Matrix34 GetInverse() const	//!< Get the inverse of this matrix
	{
		//  (4 8 - 5 7)    (5 6 - 3 8)    (3 7 - 4 6)    (3 (8 10 - 7 11) + 4 (6 11 - 8  9) + 5 (7  9 - 6 10))
		//  (2 7 - 1 8)    (0 8 - 2 6)    (1 6 - 0 7)    (0 (7 11 - 8 10) + 1 (8  9 - 6 11) + 2 (6 10 -  7 9))    / det
		//  (1 5 - 2 4)    (2 3 - 0 5)    (0 4 - 1 3)    (0 (5 10 - 4 11) + 1 (3 11 - 5  9) + 2 (4  9 - 3 10))

		Matrix34 inverse;

		T data_8_10__7_11 = cell[8] * cell[10] - cell[7] * cell[11];
		T data_6_11__8__9 = cell[6] * cell[11] - cell[8] * cell[ 9];
		T data_7__9__6_10 = cell[7] * cell[ 9] - cell[6] * cell[10];

		inverse.cell[ 0] = (cell[4]*cell[8] - cell[5]*cell[7]);
		inverse.cell[ 1] = (cell[2]*cell[7] - cell[1]*cell[8]);
		inverse.cell[ 2] = (cell[1]*cell[5] - cell[2]*cell[4]);

		inverse.cell[ 3] = (cell[5]*cell[6] - cell[3]*cell[8]);
		inverse.cell[ 4] = (cell[0]*cell[8] - cell[2]*cell[6]);
		inverse.cell[ 5] = (cell[2]*cell[3] - cell[0]*cell[5]);

		inverse.cell[ 6] = (cell[3]*cell[7] - cell[4]*cell[6]);
		inverse.cell[ 7] = (cell[1]*cell[6] - cell[0]*cell[7]);
		inverse.cell[ 8] = (cell[0]*cell[4] - cell[1]*cell[3]);

		inverse.cell[ 9] = cell[3] * data_8_10__7_11 + cell[4] * data_6_11__8__9 + cell[5] * data_7__9__6_10;
		inverse.cell[10] = cell[0] *-data_8_10__7_11 + cell[1] *-data_6_11__8__9 + cell[2] *-data_7__9__6_10;
		inverse.cell[11] = cell[0] * (cell[5] * cell[10] - cell[4] * cell[11]) + 
		                   cell[1] * (cell[3] * cell[11] - cell[5] * cell[ 9]) +
		                   cell[2] * (cell[4] * cell[ 9] - cell[3] * cell[10]);

		T det = cell[0] * inverse.cell[0] + cell[1] * inverse.cell[3] + cell[2] * inverse.cell[6];
		return inverse / det;
	}

	//! Removes the scale component of the matrix by normalizing the first three columns.
	//! The resulting matrix can contain shear, if it originally contained non-uniform scale and rotation.
	void Normalize() { Column(0).Normalize(); Column(1).Normalize(); Column(2).Normalize(); }

	//! Orthogonalizes the matrix and removes the scale component, preserving the x direction
	void OrthogonalizeX()
	{
		Column(0).Normalize();
		Column(1) -= Column(0) * (Column(1) % Column(0));
		Column(1).Normalize();
		Column(2) -= Column(0) * (Column(2) % Column(0));
		Column(2) -= Column(1) * (Column(2) % Column(1));
		Column(2).Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the y direction
	void OrthogonalizeY()
	{
		Column(1).Normalize();
		Column(0) -= Column(1) * (Column(0) % Column(1));
		Column(0).Normalize();
		Column(2) -= Column(1) * (Column(2) % Column(1));
		Column(2) -= Column(0) * (Column(2) % Column(0));
		Column(2).Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the z direction
	void OrthogonalizeZ()
	{
		Column(2).Normalize();
		Column(0) -= Column(2) * (Column(0) % Column(2));
		Column(0).Normalize();
		Column(1) -= Column(2) * (Column(1) % Column(2));
		Column(1) -= Column(0) * (Column(1) % Column(0));
		Column(1).Normalize();
	}

	//! Returns if the matrix is identity within the given error tollerance.
	bool IsIdentity( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const
	{
		return std::abs(cell[0]-T(1)) < tollerance && std::abs(cell[ 1])      < tollerance && std::abs(cell[ 2])      < tollerance && 
			   std::abs(cell[3])      < tollerance && std::abs(cell[ 4]-T(1)) < tollerance && std::abs(cell[ 5])      < tollerance &&
			   std::abs(cell[6])      < tollerance && std::abs(cell[ 7])      < tollerance && std::abs(cell[ 8]-T(1)) < tollerance &&
			   std::abs(cell[9])      < tollerance && std::abs(cell[10])      < tollerance && std::abs(cell[11])      < tollerance;
	}

	//! Returns if the matrix is symmetric within the given error tollerance.
	bool IsSymmetric( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const
	{
		return std::abs(cell[ 1] - cell[3]) < tollerance && 
			   std::abs(cell[ 2] - cell[6]) < tollerance &&
			   std::abs(cell[ 5] - cell[7]) < tollerance &&
			   std::abs(cell[ 9])           < tollerance &&
			   std::abs(cell[10])           < tollerance &&
			   std::abs(cell[11])           < tollerance;
	}

	//! Returns if the matrix is diagonal.
	bool IsDiagonal( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const
	{
		return                      std::abs(cell[ 1]) + std::abs(cell[ 2])
			 + std::abs(cell[ 3])                      + std::abs(cell[ 5])
			 + std::abs(cell[ 6]) + std::abs(cell[ 7])
			 + std::abs(cell[ 9]) + std::abs(cell[10]) + std::abs(cell[11]) < tollerance*9;
	}

	//////////////////////////////////////////////////////////////////////////
	//!@name Static Methods

	//! Returns an identity matrix
	CY_NODISCARD static Matrix34 Identity() { T c[] = { 1,0,0, 0,1,0, 0,0,1, 0,0,0 }; return Matrix34(c); }
	//! Returns a view matrix using position, target and approximate up vector
	CY_NODISCARD static Matrix34 View( Vec3<T> const &pos, Vec3<T> const &target, Vec3<T> const &up ) { Matrix34 m; m.SetView(pos,target,up); return m; }
	//! Returns a rotation matrix around x axis by angle in radians
	CY_NODISCARD static Matrix34 RotationX( T angle ) { Matrix34 m; m.SetRotationX(angle); return m; }
	//! Returns a rotation matrix around y axis by angle in radians
	CY_NODISCARD static Matrix34 RotationY( T angle ) { Matrix34 m; m.SetRotationY(angle); return m; }
	//! Returns a rotation matrix around z axis by angle in radians
	CY_NODISCARD static Matrix34 RotationZ( T angle ) { Matrix34 m; m.SetRotationZ(angle); return m; }
	//! Returns a rotation matrix around x, y, and then z axes by angle in radians (Rz * Ry * Rx)
	CY_NODISCARD static Matrix34 RotationXYZ( T angleX, T angleY, T angleZ ) { Matrix34 m; m.SetRotationXYZ(angleX,angleY,angleZ); return m; }
	//! Returns a rotation matrix around z, y, and then x axes by angle in radians (Rx * Ry * Rz)
	CY_NODISCARD static Matrix34 RotationZYX( T angleX, T angleY, T angleZ ) { Matrix34 m; m.SetRotationZYX(angleX,angleY,angleZ); return m; }
	//! Returns a rotation matrix about the given axis by angle in radians
	CY_NODISCARD static Matrix34 Rotation( Vec3<T> const &axis, T angle ) { Matrix34 m; m.SetRotation(axis,angle); return m; }
	//! Returns a rotation matrix about the given axis by cos and sin of the rotation angle
	CY_NODISCARD static Matrix34 Rotation( Vec3<T> const &axis, T cosAngle, T sinAngle ) { Matrix34 m; m.SetRotation(axis,cosAngle,sinAngle); return m; }
	//! Returns a rotation matrix that sets [from] unit vector to [to] unit vector
	CY_NODISCARD static Matrix34 Rotation( Vec3<T> const &from, Vec3<T> const &to ) { Matrix34 m; m.SetRotation(from,to); return m; }
	//! Returns a uniform scale matrix
	CY_NODISCARD static Matrix34 Scale( T uniformScale ) { Matrix34 m; m.SetScale(uniformScale); return m; }
	//! Returns a scale matrix
	CY_NODISCARD static Matrix34 Scale( T scaleX, T scaleY, T scaleZ ) { Matrix34 m; m.SetScale(scaleX,scaleY,scaleZ); return m; }
	//! Returns a scale matrix
	CY_NODISCARD static Matrix34 Scale( Vec3<T> const &scale ) { Matrix34 m; m.SetScale(scale); return m; }
	//! Returns a translation matrix with no rotation or scale
	CY_NODISCARD static Matrix34 Translation( Vec3<T> const &move ) { Matrix34 m; m.SetTranslation(move); return m; }


	//////////////////////////////////////////////////////////////////////////
};


//-------------------------------------------------------------------------------

//! 4x4 matrix class.
//!
//! Its data stores 16-value array of column-major matrix elements.
//! I chose column-major format to be compatible with OpenGL
//! You can use Matrix4 with Vec3<T> and Vec4<T>
//! to transform 3D and 4D points.

template <typename T>
class Matrix4
{
	CY_NODISCARD friend Matrix4 operator * ( T value, Matrix4 const &right ) { Matrix4 r; _CY_FOR_16i( r.cell[i]=value*right.cell[i] ); return r; }	//!< multiply matrix by a value
	CY_NODISCARD friend Matrix4 Inverse( Matrix4 const &m ) { return m.GetInverse(); }	//!< return the inverse of the matrix

	//! multiply a 4x4 matrix with a 3x4 matrix, treating it as a 4x4 matrix with last row 0,0,0,1
	CY_NODISCARD friend Matrix4 operator * ( Matrix34<T> const &left, Matrix4 const &right )
	{
		Matrix4 rm;
		for ( int i=0; i<16; i+=4 ) {
			T a[4], b[4], c[4], d[4], e[4], f[4], r[4];
			_CY_IVDEP_FOR ( int j=0; j<3; ++j ) a[j] = left.cell[  j] * right.cell[i  ];
			_CY_IVDEP_FOR ( int j=0; j<3; ++j ) b[j] = left.cell[3+j] * right.cell[i+1];
			_CY_IVDEP_FOR ( int j=0; j<3; ++j ) c[j] = left.cell[6+j] * right.cell[i+2];
			_CY_IVDEP_FOR ( int j=0; j<3; ++j ) d[j] = left.cell[9+j] * right.cell[i+3];
			_CY_IVDEP_FOR ( int j=0; j<3; ++j ) e[j] = a[j] + b[j];
			_CY_IVDEP_FOR ( int j=0; j<3; ++j ) f[j] = c[j] + d[j];
			_CY_IVDEP_FOR ( int j=0; j<3; ++j ) r[j] = e[j] + f[j];
			r[3] = right.cell[i+3];
			MemCopy( rm.cell+i, r, 4 );
		}
		return rm;
	}

public:

	//! Elements of the matrix are column-major: \n
	//! | 0   4   8  12 | \n
	//! | 1   5   9  13 | \n
	//! | 2   6  10  14 | \n
	//! | 3   7  11  15 | \n
#ifdef __cpp_unrestricted_unions
	struct ColVec3 {
		Vec3<T> v;
		T s;
		void Set( Vec3<T> const &_v, T _s ) { v=_v; s=_s; }
	};	// column vector plus scalar
	union {
		T       cell [16];
		Vec4<T> column[4];	// column vectors
		ColVec3 col3  [4];	// column vectors plus scalars
	};
#else
	T cell [16];
#endif


	//////////////////////////////////////////////////////////////////////////
	//!@name Constructors

	Matrix4() CY_CLASS_FUNCTION_DEFAULT																					//!< Default constructor
	template <typename TT> explicit Matrix4<T>( Matrix4<TT> const &matrix ) { MemConvert(cell,matrix.cell,16); }		//!< Copy constructor for different types
	explicit Matrix4( T const * restrict values ) { Set(values); }														//!< Initialize the matrix using an array of 9 values
	explicit Matrix4( T v )                      { SetScale(v); }														//!< Initialize the matrix as identity scaled by v
	explicit Matrix4( Vec3<T> const &x, Vec3<T> const &y, Vec3<T> const &z, Vec3<T> const &pos ) { Set(x,y,z,pos); }	//!< Initialize the matrix using x,y,z vectors and coordinate center
	explicit Matrix4( Vec4<T> const &x, Vec4<T> const &y, Vec4<T> const &z, Vec4<T> const &w   ) { Set(x,y,z,w);   }	//!< Initialize the matrix using x,y,z vectors as columns
	explicit Matrix4( Matrix34<T> const &m ) { Column(0).Set(m.Column(0),T(0)); Column(1).Set(m.Column(1),T(0)); Column(2).Set(m.Column(2),T(0)); Column(3).Set(m.Column(3),T(1)); }
	explicit Matrix4( Matrix3 <T> const &m ) { Column(0).Set(m.Column(0),T(0)); Column(1).Set(m.Column(1),T(0)); Column(2).Set(m.Column(2),T(0)); Column(3).Set(0,0,0,1); }
	explicit Matrix4( Matrix2 <T> const &m ) { Column(0).Set(m.Column(0),T(0),T(0)); Column(1).Set(m.Column(1),T(0),T(0)); Column(2).Set(0,0,1,0); Column(3).Set(0,0,0,1); }
	explicit Matrix4( Matrix3 <T> const &m, Vec3<T> const &pos ) { Column(0).Set(m.Column(0),T(0)); Column(1).Set(m.Column(1),T(0)); Column(2).Set(m.Column(2),T(0)); Column(3).Set(pos,T(1)); }

	//! Constructor using row-major order for initialization
	Matrix4( T c00, T c01, T c02, T c03,
		     T c10, T c11, T c12, T c13,
		     T c20, T c21, T c22, T c23,
		     T c30, T c31, T c32, T c33 )
	{
		cell[ 0] = c00;   cell[ 4] = c01;   cell[ 8] = c02;   cell[12] = c03;
		cell[ 1] = c10;   cell[ 5] = c11;   cell[ 9] = c12;   cell[13] = c13;
		cell[ 2] = c20;   cell[ 6] = c21;   cell[10] = c22;   cell[14] = c23;
		cell[ 3] = c30;   cell[ 7] = c31;   cell[11] = c32;   cell[15] = c33;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Set & Get Methods

	void Zero    ()       { MemClear(cell,16); }																					//!< Set all the values as zero
	bool IsZero  () const { return Column(0).IsZero  () && Column(1).IsZero  () && Column(2).IsZero  () && Column(3).IsZero  (); }	//!< Returns true if the matrix is exactly zero
	bool IsFinite() const { return Column(0).IsFinite() && Column(1).IsFinite() && Column(2).IsFinite() && Column(3).IsFinite(); }	//!< Returns true if all components are finite real numbers.
	void Get( T       * restrict values ) { MemCopy(values,cell,16); }																//!< Copies the matrix cell to the given values array of size 16
	void Set( T const * restrict values ) { MemCopy(cell,values,16); }																//!< Set Matrix using an array of 16 values
	void Set( Vec3<T> const &x, Vec3<T> const &y, Vec3<T> const &z, Vec3<T> const &pos ) { Column(0).Set(x,T(0)); Column(1).Set(y,T(0)); Column(2).Set(z,T(0)); Column(3).Set(pos,T(1)); }	//!< Set matrix using x,y,z column vectors and coordinate center
	void Set( Vec4<T> const &x, Vec4<T> const &y, Vec4<T> const &z, Vec4<T> const &w ) { Column(0)=x; Column(1)=y; Column(2)=z; Column(3)=w; }	//!< Set matrix using x,y,z,w column vectors
	void SetIdentity()    { SetScale(T(1)); }																						//!< Converts the matrix to an identity matrix
	void SetTensorProduct( Vec4<T> const &v0, Vec4<T> const &v1 )																	//!< Sets the matrix as the tensor product (outer product) of two vectors
	{
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) cell[   i] = v0[i] * v1.x;	 // cell[0]=v0.x*v1.x;  cell[4]=v0.x*v1.y;  cell[ 8]=v0.x*v1.z;  cell[12]=v0.x*v1.w;
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) cell[ 4+i] = v0[i] * v1.y;	 // cell[1]=v0.y*v1.x;  cell[5]=v0.y*v1.y;  cell[ 9]=v0.y*v1.z;  cell[13]=v0.y*v1.w;
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) cell[ 8+i] = v0[i] * v1.z;	 // cell[2]=v0.z*v1.x;  cell[6]=v0.z*v1.y;  cell[10]=v0.z*v1.z;  cell[14]=v0.z*v1.w;
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) cell[12+i] = v0[i] * v1.w;	 // cell[3]=v0.w*v1.x;  cell[7]=v0.w*v1.y;  cell[11]=v0.w*v1.z;  cell[15]=v0.w*v1.w;
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Affine transformations

	//! Sets a uniform scale matrix
	void SetScale( T uniformScale ) { SetScale(uniformScale,uniformScale,uniformScale); }
	//! Sets a scale matrix
	void SetScale( T scaleX, T scaleY, T scaleZ, T scaleW=T(1) )
	{
		Zero();
		cell[ 0] = scaleX;
		cell[ 5] = scaleY;
		cell[10] = scaleZ;
		cell[15] = scaleW;
	}
	//! Sets a scale matrix
	void SetScale( Vec3<T> const &scale ) { SetScale(scale.x,scale.y,scale.z); }
	//! Set as rotation matrix around x axis
	void SetRotationX( T angle ) { SetRotationX( std::sin(angle), std::cos(angle) ); }
	//! Set as rotation matrix around x axis by cos and sin of angle
	void SetRotationX( T sinAngle, T cosAngle )
	{
		cell[ 0] = T(1);  cell[ 1] =  T(0);       cell[ 2] = T(0);      cell[ 3] = T(0);
		cell[ 4] = T(0);  cell[ 5] =  cosAngle;   cell[ 6] = sinAngle;  cell[ 7] = T(0);
		cell[ 8] = T(0);  cell[ 9] = -sinAngle;   cell[10] = cosAngle;
		MemClear(cell+11,4);
		cell[15] = T(1);
	}
	//! Set as rotation matrix around y axis
	void SetRotationY( T angle ) { SetRotationY( std::sin(angle), std::cos(angle) ); }
	//! Set as rotation matrix around y axis by cos and sin of angle
	void SetRotationY( T sinAngle, T cosAngle )
	{
		cell[ 0] = cosAngle;  cell[ 1] = T(0);  cell[ 2] = -sinAngle;  cell[ 3] = T(0);
		cell[ 4] = T(0);      cell[ 5] = T(1);  cell[ 6] =  T(0);      cell[ 7] = T(0);
		cell[ 8] = sinAngle;  cell[ 9] = T(0);  cell[10] =  cosAngle;
		MemClear(cell+11,4);
		cell[15] = T(1);
	}
	//! Set as rotation matrix around z axis
	void SetRotationZ( T angle ) { SetRotationZ( std::sin(angle), std::cos(angle) ); }
	//! Set as rotation matrix around z axis by cos and sin of angle
	void SetRotationZ( T sinAngle, T cosAngle )
	{
		cell[ 0] =  cosAngle;  cell[ 1] = sinAngle;  cell[ 2] = T(0);  cell[ 3] = T(0);
		cell[ 4] = -sinAngle;  cell[ 5] = cosAngle;  cell[ 6] = T(0);  cell[ 7] = T(0); 
		cell[ 8] =  T(0);      cell[ 9] = T(0);      cell[10] = T(1);
		MemClear(cell+11,4);
		cell[15] = T(1);
	}
	//! Set as rotation matrix around x, y, and then z axes ( Rz * Ry * Rx )
	void SetRotationXYZ( T angleX, T angleY, T angleZ )
	{
		T sx = std::sin(angleX);
		T cx = std::cos(angleX);
		T sy = std::sin(angleY);
		T cy = std::cos(angleY);
		T sz = std::sin(angleZ);
		T cz = std::cos(angleZ);
		cell[ 0] = cy*cz;             cell[ 1] = cy*sz;             cell[ 2] =-sy;     cell[ 3] = T(0);
		cell[ 4] = cz*sx*sy - cx*sz;  cell[ 5] = cx*cz + sx*sy*sz;  cell[ 6] = cy*sx;  cell[ 7] = T(0);
		cell[ 8] = cx*cz*sy + sx*sz;  cell[ 9] =-cz*sx + cx*sy*sz;  cell[10] = cx*cy;
		MemClear(cell+11,4);
		cell[15] = T(1);
	}
	//! Set as rotation matrix around z, y, and then x axes ( Rx * Ry * Rz )
	void SetRotationZYX( T angleX, T angleY, T angleZ )
	{
		T sx = std::sin(angleX);
		T cx = std::cos(angleX);
		T sy = std::sin(angleY);
		T cy = std::cos(angleY);
		T sz = std::sin(angleZ);
		T cz = std::cos(angleZ);
		cell[ 0] = cy*cz;  cell[ 1] = cx*sz + sx*sy*cz;  cell[ 2] = sx*sz - cx*sy*cz;  cell[ 3] = T(0);           
		cell[ 4] =-cy*sz;  cell[ 5] = cx*cz - sx*sy*sz;  cell[ 6] = sx*cz + cx*sy*sz;  cell[ 7] = T(0);           
		cell[ 8] = sy;     cell[ 9] =-sx*cy;             cell[10] = cx*cy;
		MemClear(cell+11,4);
		cell[15] = T(1);
	}
	//! Set a rotation matrix about the given axis by angle
	void SetRotation( Vec3<T> const &axis, T angle ) { SetRotation(axis,std::sin(angle),std::cos(angle)); }
	//! Set a rotation matrix about the given axis by cos and sin of angle
	void SetRotation( Vec3<T> const &axis, T sinAngle, T cosAngle )
	{
		T t = T(1) - cosAngle;
		Vec3<T> a = t * axis;
		T txy = a.x * axis.y;
		T txz = a.x * axis.z;
		T tyz = a.y * axis.z;
		Vec3<T> s = sinAngle * axis;
		cell[ 0] = a.x * axis.x + cosAngle;   cell[ 1] = txy + s.z;                 cell[ 2] = txz - s.y;                 cell[ 3] = T(0);
		cell[ 4] = txy - s.z;                 cell[ 5] = a.y * axis.y + cosAngle;   cell[ 6] = tyz + s.x;                 cell[ 7] = T(0);
		cell[ 8] = txz + s.y;                 cell[ 9] = tyz - s.x;                 cell[10] = a.z * axis.z + cosAngle;
		MemClear(cell+11,4);
		cell[15] = T(1);
	}
	//! Set a rotation matrix that sets [from] unit vector to [to] unit vector
	void SetRotation( Vec3<T> const &from, Vec3<T> const &to )
	{
		T c = from.Dot(to);
		if ( c > T(0.9999999) ) SetIdentity();
		else {
			T s = Sqrt(T(1) - c*c);
			Vec3<T> axis = from.Cross(to).GetNormalized();
			SetRotation(axis, s, c);
		}
	}
	//! Sets a translation matrix with no rotation or scale
	void SetTranslation( Vec3<T> const &move ) { Column(0).Set(1,0,0,0); Column(1).Set(0,1,0,0); Column(2).Set(0,0,1,0); Column(3).Set(move,1); }
	//! Adds a translation to the matrix
	void AddTranslation( Vec3<T> const &move ) { cell[12]+=move.x; cell[13]+=move.y; cell[14]+=move.z; }
	//! Sets the translation component of the matrix
	void SetTranslationComponent( Vec3<T> const &move ) { cell[12]=move.x; cell[13]=move.y; cell[14]=move.z; }
	//! Sets the translation component of the matrix to zero
	void SetNoTranslation() { cell[12]=0; cell[13]=0; cell[14]=0; }
	//! Set view matrix using position, target and approximate up vector
	void SetView( Vec3<T> const &pos, Vec3<T> const &target, Vec3<T> const &up )
	{
		Vec3<T> f = target - pos;
		f.Normalize();
		Vec3<T> s = f.Cross(up);
		s.Normalize();
		Vec3<T> u = s.Cross(f);
		cell[ 0]=s.x; cell[ 1]=u.x; cell[ 2]=-f.x; cell[ 3]=T(0);
		cell[ 4]=s.y; cell[ 5]=u.y; cell[ 6]=-f.y; cell[ 7]=T(0);
		cell[ 8]=s.z; cell[ 9]=u.z; cell[10]=-f.z; cell[11]=T(0);
		cell[12]= -s % pos;
		cell[13]= -u % pos;
		cell[14]=  f % pos;
		cell[15]=T(1);
	}
	//! Sets a Cartesian coordinate frame using the given x direction, an approximate y direction, and a translation. x must be a unit vector.
	void SetCartesianFrameXY( Vec3<T> const &x, Vec3<T> const &y_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> z = x.Cross(y_approx); z.Normalize(); Vec3<T> y=z.Cross(x); Set(x,y,z,trans); }
	//! Sets a Cartesian coordinate frame using the given x direction, an approximate z direction, and a translation. x must be a unit vector.
	void SetCartesianFrameXZ( Vec3<T> const &x, Vec3<T> const &z_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> y = z_approx.Cross(x); y.Normalize(); Vec3<T> z=x.Cross(y); Set(x,y,z,trans); }
	//! Sets a Cartesian coordinate frame using the given y direction, an approximate x direction, and a translation. y must be a unit vector.
	void SetCartesianFrameYX( Vec3<T> const &y, Vec3<T> const &x_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> z = x_approx.Cross(y); z.Normalize(); Vec3<T> x=y.Cross(z); Set(x,y,z,trans); }
	//! Sets a Cartesian coordinate frame using the given y direction, an approximate z direction, and a translation. y must be a unit vector.
	void SetCartesianFrameYZ( Vec3<T> const &y, Vec3<T> const &z_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> x = y.Cross(z_approx); x.Normalize(); Vec3<T> z=x.Cross(y); Set(x,y,z,trans); }
	//! Sets a Cartesian coordinate frame using the given z direction, an approximate x direction, and a translation. z must be a unit vector.
	void SetCartesianFrameZX( Vec3<T> const &z, Vec3<T> const &x_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> y = z.Cross(x_approx); y.Normalize(); Vec3<T> x=y.Cross(z); Set(x,y,z,trans); }
	//! Sets a Cartesian coordinate frame using the given z direction, an approximate y direction, and a translation. z must be a unit vector.
	void SetCartesianFrameZY( Vec3<T> const &z, Vec3<T> const &y_approx, Vec3<T> const &trans=Vec3<T>(T(0),T(0),T(0)) ) { Vec3<T> x = y_approx.Cross(z); x.Normalize(); Vec3<T> y=z.Cross(x); Set(x,y,z,trans); }
	//! Set a project matrix with field of view in radians
	void SetPerspective( T fov, T aspect, T znear, T zfar ) { SetPerspectiveTan(std::tan(fov*T(0.5)),aspect,znear,zfar); }
	//! Set a project matrix with the tangent of the half field of view (tan_fov_2)
	void SetPerspectiveTan( T tan_fov_2, T aspect, T znear, T zfar )
	{
		T yScale = T(1) / tan_fov_2;
		T xScale = yScale / aspect;
		T zdif = znear - zfar;
		Column(0).Set(xScale,0,0,0);
		Column(1).Set(0,yScale,0,0);
		Column(2).Set(0,0,(zfar+znear)/zdif,-1);
		Column(3).Set(0,0,(2*zfar*znear)/zdif,0);
	}


	//////////////////////////////////////////////////////////////////////////
	//!@name Set Row, Column, or Diagonal

	void SetRow     ( int ri, T x, T y, T z, T w ) { cell[ri]=x; cell[ri+4]=y; cell[ri+8]=z; cell[ri+12]=w; }	//!< Sets a row of the matrix
	void SetRow     ( int ri, Vec4<T> const &v )   { SetRow(ri,v.x,v.y,v.z,v.w); }								//!< Sets a row of the matrix
	void SetColumn  ( int ci, T x, T y, T z, T w ) { Column(ci).Set(x,y,z,w); }									//!< Sets a column of the matrix
	void SetColumn  ( int ci, Vec4<T> const &v )   { Column(ci)=v; }											//!< Sets a column of the matrix
	void SetDiagonal( T xx, T yy, T zz, T ww=1 )   { cell[0]=xx; cell[5]=yy; cell[10]=zz; cell[15]=ww; }		//!< Sets the diagonal values of the matrix
	void SetDiagonal( Vec4<T> const &p )           { SetDiagonal( p.x, p.y, p.z, p.w ); }						//!< Sets the diagonal values of the matrix
	void SetDiagonal( Vec3<T> const &p )           { SetDiagonal( p.x, p.y, p.z, T(1) ); }						//!< Sets the diagonal values of the matrix
	void SetDiagonal( T const * restrict values )  { SetDiagonal(values[0],values[1],values[2],values[3]); }	//!< Sets the 4 diagonal values of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Row, Column, or Diagonal

#ifdef __cpp_unrestricted_unions
	CY_NODISCARD Vec4<T>       * Columns()               { return column; }
	CY_NODISCARD Vec4<T> const * Columns()         const { return column; }
	CY_NODISCARD Vec4<T>       & Column ( int ci )       { return column[ci]; }
	CY_NODISCARD Vec4<T> const & Column ( int ci ) const { return column[ci]; }
	CY_NODISCARD Vec3<T>       & Column3( int ci )       { return col3[ci].v; }
	CY_NODISCARD Vec3<T> const & Column3( int ci ) const { return col3[ci].v; }
#else
	CY_NODISCARD Vec4<T>       * Columns()               { return ((Vec4<T>*)cell); }
	CY_NODISCARD Vec4<T> const * Columns()         const { return ((Vec4<T>*)cell); }
	CY_NODISCARD Vec4<T>       & Column ( int ci )       { return Columns()[ci]; }
	CY_NODISCARD Vec4<T> const & Column ( int ci ) const { return Columns()[ci]; }
	CY_NODISCARD Vec3<T>       & Column3( int ci )       { return (Vec3<T>       &)cell[ci*4]; }
	CY_NODISCARD Vec3<T> const & Column3( int ci ) const { return (Vec3<T> const &)cell[ci*4]; }
#endif
	CY_NODISCARD Vec4<T>         GetRow ( int ri ) const { return Vec4<T>( cell[ri], cell[ri+4], cell[ri+8], cell[ri+12] ); }	//!< Returns a row of the matrix
	CY_NODISCARD Vec4<T>         GetDiagonal()     const { return Vec4<T>( cell[0], cell[5], cell[10], cell[15] ); }				//!< Returns the diagonal of the matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Get Sub-matrix cell

	CY_NODISCARD Matrix34<T> GetSubMatrix34() const { Matrix34<T> m; MemCopy(m.cell,cell,3); MemCopy(m.cell+3,cell+4,3); MemCopy(m.cell+6,cell+8,3); MemCopy(m.cell+9,cell+12,3); return m; }	//!< Returns the 3x4 portion of the matrix
	CY_NODISCARD Matrix3<T>  GetSubMatrix3 () const { Matrix3<T>  m; MemCopy(m.cell,cell,3); MemCopy(m.cell+3,cell+4,3); MemCopy(m.cell+6,cell+8,3); return m; }	//!< Returns the 3x3 portion of the matrix
	CY_NODISCARD Matrix2<T>  GetSubMatrix2 () const { Matrix2<T>  m; MemCopy(m.cell,cell,2); MemCopy(m.cell+2,cell+4,2); return m; }	//!< Returns the 2x2 portion of the matrix
	CY_NODISCARD Vec3<T>     GetTranslation() const { return Vec3<T>(cell+12); }						//!< Returns the translation component of the matrix
	CY_NODISCARD Matrix3<T>  GetRotation   () const { Matrix3<T> m(*this); return m.GetRotation(); }	//!< Returns the rotation portion of the transformation
	CY_NODISCARD Matrix3<T>  GetScale      () const { Matrix3<T> m(*this); return m.GetScale   (); }	//!< Returns the scale portion of the transformation.

	//! Returns the average scale factor of the 3 by 3 sub-matrix
	CY_NODISCARD T GetAvrgScale() const 
	{
		T det = cell[0] * ( cell[5] * cell[10] - cell[6] * cell[ 9] ) + 
		        cell[1] * ( cell[6] * cell[ 8] - cell[4] * cell[10] ) + 
		        cell[2] * ( cell[4] * cell[ 9] - cell[5] * cell[ 8] );
		T s = std::pow( std::abs(det), T(1)/T(3) );
		return det >= 0 ? s : -s;
	}

	void GetComponents( Matrix3<T> &scale, Matrix3<T> &rotation, Vec3<T> &translation ) const { Matrix3<T> m(*this); m.GetComponents(scale,rotation); translation=GetTranslation(); }	//!< Returns separate transformation components

	//////////////////////////////////////////////////////////////////////////
	//!@name Comparison Operators
	
	CY_NODISCARD bool operator == ( Matrix4 const &right ) const { _CY_FOR_16i( if ( cell[i] != right.cell[i] ) return false ); return true;  } //!< compare equal
	CY_NODISCARD bool operator != ( Matrix4 const &right ) const { _CY_FOR_16i( if ( cell[i] != right.cell[i] ) return true  ); return false; } //!< compare not equal


	//////////////////////////////////////////////////////////////////////////
	//!@name Access Operators

	CY_NODISCARD T&        operator () ( int ri, int ci )       { assert( ri>=0 && ri<4 && ci>=0 && ci<4 ); return cell[ ci*4 + ri ]; }	//!< subscript operator
	CY_NODISCARD T const & operator () ( int ri, int ci ) const { assert( ri>=0 && ri<4 && ci>=0 && ci<4 ); return cell[ ci*4 + ri ]; }	//!< constant subscript operator
	CY_NODISCARD T&        operator [] ( int i )                { assert( i>=0 && i<16 ); return cell[i]; }								//!< subscript operator
	CY_NODISCARD T const & operator [] ( int i )          const { assert( i>=0 && i<16 ); return cell[i]; }								//!< constant subscript operator


	//////////////////////////////////////////////////////////////////////////
	//!@name Unary and Binary Operators

	// Unary operators
	CY_NODISCARD Matrix4 operator - () const { Matrix4 r; _CY_FOR_16i( r.cell[i] = -cell[i] ); return r; }	//!< negative matrix

	// Binary operators
	CY_NODISCARD Matrix4 operator * ( T       const  value ) const { Matrix4 r; _CY_FOR_16i( r.cell[i] = cell[i] * value         ); return r; }	//!< multiply matrix by a value
	CY_NODISCARD Matrix4 operator / ( T       const  value ) const { Matrix4 r; _CY_FOR_16i( r.cell[i] = cell[i] / value         ); return r; }	//!< divide matrix by a value
	CY_NODISCARD Matrix4 operator + ( Matrix4 const &right ) const { Matrix4 r; _CY_FOR_16i( r.cell[i] = cell[i] + right.cell[i] ); return r; }	//!< add two Matrices
	CY_NODISCARD Matrix4 operator - ( Matrix4 const &right ) const { Matrix4 r; _CY_FOR_16i( r.cell[i] = cell[i] - right.cell[i] ); return r; }	//!< subtract one Matrix4 from another
	CY_NODISCARD Matrix4 operator * ( Matrix4 const &right ) const	//!< multiply a matrix with another
	{
		Matrix4 rm;
		for ( int i=0; i<16; i+=4 ) {
			T a[4], b[4], c[4], d[4], e[4], f[4], r[4];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) a[j] = cell[   j] * right.cell[i  ];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) b[j] = cell[ 4+j] * right.cell[i+1];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) c[j] = cell[ 8+j] * right.cell[i+2];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) d[j] = cell[12+j] * right.cell[i+3];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) e[j] = a[j] + b[j];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) f[j] = c[j] + d[j];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) r[j] = e[j] + f[j];
			MemCopy( rm.cell+i, r, 4 );
		}
		return rm;
	}
	CY_NODISCARD Matrix4 operator * ( Matrix34<T> const &right ) const	//!< multiply a matrix with another
	{
		T a[4], b[4], c[4], d[4], e[4], r[4];
		Matrix4 rm;
		T *rd = rm.cell;
		for ( int i=0; i<9; i+=3, rd+=4 ) {
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) a[k] = cell[  k] * right.cell[i  ];
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) b[k] = cell[4+k] * right.cell[i+1];
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) c[k] = cell[8+k] * right.cell[i+2];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) d[j] = a[j] + b[j];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) r[j] = d[j] + c[j];
			MemCopy( rd, r, 4 );
		}
		_CY_IVDEP_FOR ( int k=0; k<4; ++k ) a[k] = cell[  k] * right.cell[ 9];
		_CY_IVDEP_FOR ( int k=0; k<4; ++k ) b[k] = cell[4+k] * right.cell[10];
		_CY_IVDEP_FOR ( int k=0; k<4; ++k ) c[k] = cell[8+k] * right.cell[11];
		_CY_IVDEP_FOR ( int j=0; j<4; ++j ) d[j] = a[j] + b[j];
		_CY_IVDEP_FOR ( int j=0; j<4; ++j ) e[j] = c[j] + cell[12+j];
		_CY_IVDEP_FOR ( int j=0; j<4; ++j ) r[j] = d[j] + e[j];
		MemCopy( rd, r, 4 );
		return rm;
	}
	CY_NODISCARD Matrix4 operator * ( Matrix3<T> const &right ) const	//!< multiply a matrix with another
	{
		T a[4], b[4], c[4], d[4], r[4];
		Matrix4 rm;
		T *rd = rm.cell;
		for ( int i=0; i<9; i+=3, rd+=4 ) {
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) a[k] = cell[  k] * right.cell[i  ];
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) b[k] = cell[4+k] * right.cell[i+1];
			_CY_IVDEP_FOR ( int k=0; k<4; ++k ) c[k] = cell[8+k] * right.cell[i+2];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) d[j] = a[j] + b[j];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) r[j] = d[j] + c[j];
			MemCopy( rd, r, 4 );
		}
		MemCopy( rm.cell+12, cell+12, 4 );
		return rm;
	}
	CY_NODISCARD Vec4<T> operator * ( Vec3<T> const &p ) const 
	{
		//return Vec4<T>( p.x*cell[0] + p.y*cell[4] + p.z*cell[ 8] + cell[12], 
		//                p.x*cell[1] + p.y*cell[5] + p.z*cell[ 9] + cell[13],
		//                p.x*cell[2] + p.y*cell[6] + p.z*cell[10] + cell[14],
		//                p.x*cell[3] + p.y*cell[7] + p.z*cell[11] + cell[15] );
		Vec4<T> r;
		T a[4], b[4], c[4], d[4], e[4];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) a[i] = p.x * cell[i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) b[i] = p.y * cell[4+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) c[i] = p.z * cell[8+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) d[i] = a[i] + b[i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) e[i] = c[i] + cell[12+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) r[i] = d[i] + e[i];
		return r;
	}
	CY_NODISCARD Vec4<T> operator * ( Vec4<T> const &p ) const 
	{
		//return Vec4<T>( p.x*cell[0] + p.y*cell[4] + p.z*cell[ 8] + p.w*cell[12],
		//                p.x*cell[1] + p.y*cell[5] + p.z*cell[ 9] + p.w*cell[13],
		//                p.x*cell[2] + p.y*cell[6] + p.z*cell[10] + p.w*cell[14],
		//                p.x*cell[3] + p.y*cell[7] + p.z*cell[11] + p.w*cell[15] );
		Vec4<T> r;
		T a[4], b[4], c[4], d[4], e[4], f[4];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) a[i] = p.x * cell[i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) b[i] = p.y * cell[4+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) c[i] = p.z * cell[8+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) d[i] = p.w * cell[12+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) e[i] = a[i] + b[i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) f[i] = c[i] + d[i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) r[i] = e[i] + f[i];
		return r;
	}

	//////////////////////////////////////////////////////////////////////////
	//!@name 3D Vector Transform Methods

	//! Transforms the vector by multiplying it with the matrix, ignoring the translation component.
	CY_NODISCARD Vec4<T> VectorTransform( Vec3<T> const &p ) const
	{
		//return Vec4<T>( p.x*cell[0] + p.y*cell[4] + p.z*cell[ 8], 
		//                p.x*cell[1] + p.y*cell[5] + p.z*cell[ 9],
		//                p.x*cell[2] + p.y*cell[6] + p.z*cell[10],
		//                p.x*cell[3] + p.y*cell[7] + p.z*cell[11] );
		Vec4<T> r;
		T a[4], b[4], c[4], d[4];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) a[i] = p.x * cell[i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) b[i] = p.y * cell[4+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) c[i] = p.z * cell[8+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) d[i] = a[i] + b[i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) r[i] = d[i] + c[i];
		return r;
	}

	//////////////////////////////////////////////////////////////////////////
	//!@name Assignment Operators

	Matrix4 const & operator += ( Matrix4     const &right ) { _CY_FOR_16i( cell[i] += right.cell[i] ); return *this; }	//!< add two Matrices modify this
	Matrix4 const & operator -= ( Matrix4     const &right ) { _CY_FOR_16i( cell[i] -= right.cell[i] ); return *this; }	//!< subtract one Matrix4 from another matrix and modify this matrix
	Matrix4 const & operator *= ( Matrix4     const &right ) { *this = operator*(right);                return *this; }	//!< multiply a matrix with another matrix and modify this matrix
	Matrix4 const & operator *= ( Matrix34<T> const &right ) { *this = operator*(right);                return *this; }	//!< multiply a matrix with another matrix and modify this matrix
	Matrix4 const & operator *= ( Matrix3<T>  const &right ) { *this = operator*(right);                return *this; }	//!< multiply a matrix with another matrix and modify this matrix
	Matrix4 const & operator *= ( T           const  value ) { _CY_FOR_16i( cell[i] *= value );         return *this; }	//!< multiply a matrix with a value modify this matrix
	Matrix4 const & operator /= ( T           const  value ) { _CY_FOR_16i( cell[i] /= value );         return *this; }	//!< divide the matrix by a value modify the this matrix


	//////////////////////////////////////////////////////////////////////////
	//!@name Other Methods

	void Transpose()			//!< Transpose this matrix
	{
		for ( int i = 1; i < 4; ++i ) {
			for ( int j = 0; j < i; j++) {
				T temp = cell[i * 4 + j];
				cell[i * 4 + j] = cell[j * 4 + i];
				cell[j * 4 + i] = temp;
			}
		}
	}
	CY_NODISCARD Matrix4 GetTranspose() const	//!< Return the transpose of this matrix
	{
		Matrix4 m;
		m.cell[ 0] = cell[0];   m.cell[ 1] = cell[4];   m.cell[ 2] = cell[ 8];  m.cell[ 3] = cell[12];
		m.cell[ 4] = cell[1];   m.cell[ 5] = cell[5];   m.cell[ 6] = cell[ 9];  m.cell[ 7] = cell[13];
		m.cell[ 8] = cell[2];   m.cell[ 9] = cell[6];   m.cell[10] = cell[10];  m.cell[11] = cell[14];
		m.cell[12] = cell[3];   m.cell[13] = cell[7];   m.cell[14] = cell[11];  m.cell[15] = cell[15];
		return m;
	}

	//! Multiply the give vector with the transpose of the matrix
	CY_NODISCARD Vec4<T> TransposeMult( Vec3<T> const &p ) const { return TransposeMult( Vec4<T>(p.x, p.y, p.z, T(1)) ); }

	//! Multiply the give vector with the transpose of the matrix
	CY_NODISCARD Vec4<T> TransposeMult( Vec4<T> const &p ) const 
	{
		//return Vec4<T>(	p.x*cell[ 0] + p.y*cell[ 1] + p.z*cell[ 2] + p.w*cell[ 3],
		//						p.x*cell[ 4] + p.y*cell[ 5] + p.z*cell[ 6] + p.w*cell[ 7],
		//						p.x*cell[ 8] + p.y*cell[ 9] + p.z*cell[10] + p.w*cell[11],
		//						p.x*cell[12] + p.y*cell[13] + p.z*cell[14] + p.w*cell[15] );
		T a[4], b[4], c[4], d[4];
		T const *pd = p.Elements();
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) a[i] = pd[i] * cell[   i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) b[i] = pd[i] * cell[ 4+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) c[i] = pd[i] * cell[ 8+i];
		_CY_IVDEP_FOR ( int i=0; i<4; ++i ) d[i] = pd[i] * cell[12+i];
		Vec4<T> rr;
		rr.x = a[0] + a[1] + a[2] + a[3];
		rr.y = b[0] + b[1] + b[2] + b[3];
		rr.z = c[0] + c[1] + c[2] + c[3];
		rr.w = d[0] + d[1] + d[2] + d[3];
		return rr;
	}

	CY_NODISCARD Matrix4 TransposeMult( Matrix4 const & right ) const //!< Multiply a matrix by the transpose of this one (i.e. this^T * right).
	{
		Matrix4 r;
		for ( int i=0, k=0; i<3; ++i ) {
			for ( int j=0; j<3; ++j, ++k ) {
				r.cell[k] = Column(j).Dot( right.Column(i) );
			}
		}
		return r;
	}
	CY_NODISCARD Matrix4 MultTranspose( Matrix4 const & right ) const //!< Multiply the transpose of a matrix by this one (i.e. this * right^T).
	{
		Matrix4 rm;
		T* rd = rm.cell;
		for ( int i=0; i<4; ++i, rd+=4 ) {
			T a[4], b[4], c[4], d[4], e[4], f[4], r[4];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) a[j] = cell[j] * right.cell[i];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) b[j] = cell[4+j] * right.cell[i+ 4];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) c[j] = cell[8+j] * right.cell[i+ 8];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) d[j] = cell[12+j] * right.cell[i+12];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) e[j] = a[j] + b[j];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) f[j] = c[j] + d[j];
			_CY_IVDEP_FOR ( int j=0; j<4; ++j ) r[j] = e[j] + f[j];
			MemCopy( rd, r, 4 );
		}
		return rm;
	}

	CY_NODISCARD Matrix4 TransposeMultSelf() const { return TransposeMult(*this); } //!< Multiply the transpose of this matrix with itself (i.e. this^T * this).
	CY_NODISCARD Matrix4 MultSelfTranspose() const { return MultTranspose(*this); } //!< Multiply the matrix with its transpose (i.e. this * this^T).

	CY_NODISCARD T GetTrace() const { return cell[0]+cell[5]+cell[10]+cell[15]; }

	CY_NODISCARD T GetDeterminant() const	//!< Get the determinant of this matrix
	{
		// 0 (                      5 ( 10 15 - 11 14) + 6 ( 11 13 -  9 15) + 7 (  9 14 - 10 13)) +  
		// 1 ( 4 ( 11 14 - 10 15) +                      6 (  8 15 - 11 12) + 7 ( 10 12 -  8 14)) +  
		// 2 ( 4 (  9 15 - 11 13) + 5 ( 11 12 -  8 15) +                      7 (  8 13 -  9 12)) + 
		// 3 ( 4 ( 10 13 -  9 14) + 5 (  8 14 - 10 12) + 6 (  9 12 -  8 13)) 

		T data_11_14__10_15 = cell[11] * cell[14] - cell[10] * cell[15];
		T data__9_15__11_13 = cell[ 9] * cell[15] - cell[11] * cell[13];
		T data_10_13___9_14 = cell[10] * cell[13] - cell[ 9] * cell[14];
		T data_11_12___8_15 = cell[11] * cell[12] - cell[ 8] * cell[15];
		T data__8_14__10_12 = cell[ 8] * cell[14] - cell[10] * cell[12];
		T data__9_12___8_13 = cell[ 9] * cell[12] - cell[ 8] * cell[13];
		return cell[0] * ( cell[5] * (-data_11_14__10_15 ) + cell[6] * (-data__9_15__11_13 ) + cell[7] * (-data_10_13___9_14 ) ) +  
		       cell[1] * ( cell[4] * ( data_11_14__10_15 ) + cell[6] * (-data_11_12___8_15 ) + cell[7] * (-data__8_14__10_12 ) ) +  
		       cell[2] * ( cell[4] * ( data__9_15__11_13 ) + cell[5] * ( data_11_12___8_15 ) + cell[7] * (-data__9_12___8_13 ) ) + 
		       cell[3] * ( cell[4] * ( data_10_13___9_14 ) + cell[5] * ( data__8_14__10_12 ) + cell[6] * ( data__9_12___8_13 ) );
	}

	void Invert() { *this = GetInverse(); }	//!< Invert this matrix
	CY_NODISCARD Matrix4 GetInverse() const	//!< Get the inverse of this matrix
	{
		//                       5 ( 10 15 - 11 14 ) + 6 ( 11 13 -  9 15 ) + 7 (  9 14 - 10 13 )
		//                       1 ( 11 14 - 10 15 ) + 2 (  9 15 - 11 13 ) + 3 ( 10 13 -  9 14 ) 
		//                       1 (  6 15 -  7 14 ) + 2 (  7 13 -  5 15 ) + 3 (  5 14 -  6 13 )
		//                       1 (  7 10 -  6 11 ) + 2 (  5 11 -  7  9 ) + 3 (  6  9 -  5 10 )
		//					 					   						 					   
		// 4 ( 11 14 - 10 15 ) +                       6 (  8 15 - 11 12 ) + 7 ( 10 12 -  8 14 )
		// 0 ( 10 15 - 11 14 ) +                       2 ( 11 12 -  8 15 ) + 3 (  8 14 - 10 12 )
		// 0 (  7 14 -  6 15 ) +                       2 (  4 15 -  7 12 ) + 3 (  6 12 -  4 14 )      / det
		// 0 (  6 11 -  7 10 ) +                       2 (  8  7 -  4 11 ) + 3 (  4 10 -  6  8 )
		//					 					   						 					   
		// 4 (  9 15 - 11 13 ) + 5 ( 11 12 -  8 15 ) +                       7 (  8 13 -  9 12 ) 
		// 0 ( 11 13 -  9 15 ) + 1 (  8 15 - 11 12 ) +                       3 (  9 12 -  8 13 )
		// 0 (  5 15 -  7 13 ) + 1 (  7 12 -  4 15 ) +                       3 (  4 13 -  5 12 )
		// 0 (  7  9 -  5 11 ) + 1 (  4 11 -  7  8 ) +                       3 (  5  8 -  4  9 )
		//					 					   						 
		// 4 ( 10 13 -  9 14 ) + 5 (  8 14 - 10 12 ) + 6 (  9 12 -  8 13 )
		// 0 (  9 14 - 10 13 ) + 1 ( 10 12 -  8 14 ) + 2 (  8 13 -  9 12 )
		// 0 (  6 13 -  5 14 ) + 1 (  4 14 -  6 12 ) + 2 (  5 12 -  4 13 )
		// 0 (  5 10 -  6  9 ) + 1 (  6  8 -  4 10 ) + 2 (  4  9 -  5  8 )

		Matrix4 inverse;

		T data_11_14__10_15 = cell[11] * cell[14] - cell[10] * cell[15];
		T data_10_15__11_14 = cell[10] * cell[15] - cell[11] * cell[14];
		T data__7_14___6_15 = cell[ 7] * cell[14] - cell[ 6] * cell[15];
		T data__6_11___7_10 = cell[ 6] * cell[11] - cell[ 7] * cell[10];

		T data__9_15__11_13 = cell[ 9] * cell[15] - cell[11] * cell[13];
		T data_11_13___9_15 = cell[11] * cell[13] - cell[ 9] * cell[15];
		T data__5_15___7_13 = cell[ 5] * cell[15] - cell[ 7] * cell[13];
		T data__7__9___5_11 = cell[ 7] * cell[ 9] - cell[ 5] * cell[11];
		
		T data_10_13___9_14 = cell[10] * cell[13] - cell[ 9] * cell[14];
		T data__9_14__10_13 = cell[ 9] * cell[14] - cell[10] * cell[13];
		T data__6_13___5_14 = cell[ 6] * cell[13] - cell[ 5] * cell[14];
		T data__5_10___6__9 = cell[ 5] * cell[10] - cell[ 6] * cell[ 9];
		
		T data_11_12___8_15 = cell[11] * cell[12] - cell[ 8] * cell[15];
		T data__8_15__11_12 = cell[ 8] * cell[15] - cell[11] * cell[12];
		T data__7_12___4_15 = cell[ 7] * cell[12] - cell[ 4] * cell[15];
		T data__4_11___7__8 = cell[ 4] * cell[11] - cell[ 7] * cell[ 8];
		
		T data__8_14__10_12 = cell[ 8] * cell[14] - cell[10] * cell[12];
		T data_10_12___8_14 = cell[10] * cell[12] - cell[ 8] * cell[14];
		T data__4_14___6_12 = cell[ 4] * cell[14] - cell[ 6] * cell[12];
		T data__6__8___4_10 = cell[ 6] * cell[ 8] - cell[ 4] * cell[10];
		
		T data__9_12___8_13 = cell[ 9] * cell[12] - cell[ 8] * cell[13];
		T data__8_13___9_12 = cell[ 8] * cell[13] - cell[ 9] * cell[12];
		T data__5_12___4_13 = cell[ 5] * cell[12] - cell[ 4] * cell[13];
		T data__4__9___5__8 = cell[ 4] * cell[ 9] - cell[ 5] * cell[ 8];

		inverse.cell[ 0] = cell[5] * (-data_11_14__10_15) + cell[6] * (-data__9_15__11_13) + cell[7] * (-data_10_13___9_14);
		inverse.cell[ 1] = cell[1] * (-data_10_15__11_14) + cell[2] * (-data_11_13___9_15) + cell[3] * (-data__9_14__10_13);
		inverse.cell[ 2] = cell[1] * (-data__7_14___6_15) + cell[2] * (-data__5_15___7_13) + cell[3] * (-data__6_13___5_14);
		inverse.cell[ 3] = cell[1] * (-data__6_11___7_10) + cell[2] * (-data__7__9___5_11) + cell[3] * (-data__5_10___6__9);
		
		inverse.cell[ 4] = cell[4] * ( data_11_14__10_15) + cell[6] * (-data_11_12___8_15) + cell[7] * (-data__8_14__10_12);
		inverse.cell[ 5] = cell[0] * ( data_10_15__11_14) + cell[2] * (-data__8_15__11_12) + cell[3] * (-data_10_12___8_14);
		inverse.cell[ 6] = cell[0] * ( data__7_14___6_15) + cell[2] * (-data__7_12___4_15) + cell[3] * (-data__4_14___6_12);
		inverse.cell[ 7] = cell[0] * ( data__6_11___7_10) + cell[2] * (-data__4_11___7__8) + cell[3] * (-data__6__8___4_10);
		
		inverse.cell[ 8] = cell[4] * ( data__9_15__11_13) + cell[5] * ( data_11_12___8_15) + cell[7] * (-data__9_12___8_13);
		inverse.cell[ 9] = cell[0] * ( data_11_13___9_15) + cell[1] * ( data__8_15__11_12) + cell[3] * (-data__8_13___9_12);
		inverse.cell[10] = cell[0] * ( data__5_15___7_13) + cell[1] * ( data__7_12___4_15) + cell[3] * (-data__5_12___4_13);
		inverse.cell[11] = cell[0] * ( data__7__9___5_11) + cell[1] * ( data__4_11___7__8) + cell[3] * (-data__4__9___5__8);

		inverse.cell[12] = cell[4] * ( data_10_13___9_14) + cell[5] * ( data__8_14__10_12) + cell[6] * ( data__9_12___8_13);
		inverse.cell[13] = cell[0] * ( data__9_14__10_13) + cell[1] * ( data_10_12___8_14) + cell[2] * ( data__8_13___9_12);
		inverse.cell[14] = cell[0] * ( data__6_13___5_14) + cell[1] * ( data__4_14___6_12) + cell[2] * ( data__5_12___4_13);
		inverse.cell[15] = cell[0] * ( data__5_10___6__9) + cell[1] * ( data__6__8___4_10) + cell[2] * ( data__4__9___5__8);

		T det = cell[0] * inverse.cell[0] + cell[1] * inverse.cell[4] + cell[2] * inverse.cell[8] + cell[3] * inverse.cell[12];
		return inverse / det;
	}

	//! Removes the scale component of the matrix by normalizing each column of the 3x3 sub-matrix.
	//! The resulting matrix can contain shear, if it originally contained non-uniform scale and rotation.
	void Normalize() { Column3(0).Normalize(); Column3(0).Normalize(); Column3(0).Normalize(); }

	//! Orthogonalizes the matrix and removes the scale component, preserving the x direction
	void OrthogonalizeX()
	{
		Column3(0).Normalize();
		Column3(1) -= Column3(0) * (Column3(1) % Column3(0));
		Column3(1).Normalize();
		Column3(2) -= Column3(0) * (Column3(2) % Column3(0));
		Column3(2) -= Column3(1) * (Column3(2) % Column3(1));
		Column3(2).Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the y direction
	void OrthogonalizeY()
	{
		Column3(1).Normalize();
		Column3(0) -= Column3(1) * (Column3(0) % Column3(1));
		Column3(0).Normalize();
		Column3(2) -= Column3(1) * (Column3(2) % Column3(1));
		Column3(2) -= Column3(0) * (Column3(2) % Column3(0));
		Column3(2).Normalize();
	}
	//! Orthogonalizes the matrix and removes the scale component, preserving the z direction
	void OrthogonalizeZ()
	{
		Column3(2).Normalize();
		Column3(0) -= Column3(2) * (Column3(0) % Column3(2));
		Column3(0).Normalize();
		Column3(1) -= Column3(2) * (Column3(1) % Column3(2));
		Column3(1) -= Column3(0) * (Column3(1) % Column3(0));
		Column3(1).Normalize();
	}

	//! Returns if the matrix is identity within the given error tollerance.
	CY_NODISCARD bool IsIdentity( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const
	{
		return std::abs(cell[ 0]-T(1)) < tollerance && std::abs(cell[ 1])      < tollerance && std::abs(cell[ 2])      < tollerance && std::abs(cell[ 3])      < tollerance && 
			   std::abs(cell[ 4])      < tollerance && std::abs(cell[ 5]-T(1)) < tollerance && std::abs(cell[ 6])      < tollerance && std::abs(cell[ 7])      < tollerance &&
			   std::abs(cell[ 8])      < tollerance && std::abs(cell[ 9])      < tollerance && std::abs(cell[10]-T(1)) < tollerance && std::abs(cell[11])      < tollerance &&
			   std::abs(cell[12])      < tollerance && std::abs(cell[13])      < tollerance && std::abs(cell[14])      < tollerance && std::abs(cell[15]-T(1)) < tollerance;
	}

	//! Returns if the matrix is symmetric within the given error tollerance.
	CY_NODISCARD bool IsSymmetric( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const
	{
		return std::abs(cell[ 1] - cell[ 4]) < tollerance && 
			   std::abs(cell[ 2] - cell[ 8]) < tollerance &&
			   std::abs(cell[ 3] - cell[12]) < tollerance &&
			   std::abs(cell[ 6] - cell[ 9]) < tollerance &&
			   std::abs(cell[ 7] - cell[13]) < tollerance &&
			   std::abs(cell[11] - cell[14]) < tollerance;
	}

	//! Returns if the matrix is diagonal.
	CY_NODISCARD bool IsDiagonal( T tollerance=T(_CY_VEC_DEFAULT_ERROR_TOLERANCE) ) const
	{
		return                      std::abs(cell[ 1]) + std::abs(cell[ 2]) + std::abs(cell[ 3])
			 + std::abs(cell[ 4])                      + std::abs(cell[ 6]) + std::abs(cell[ 7])
			 + std::abs(cell[ 8]) + std::abs(cell[ 9])                      + std::abs(cell[11])
			 + std::abs(cell[12]) + std::abs(cell[13]) + std::abs(cell[14]) < tollerance*12;
	}

	//////////////////////////////////////////////////////////////////////////
	//!@name Static Methods

	//! Returns an identity matrix
	CY_NODISCARD static Matrix4 Identity() { T c[] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 }; return Matrix4(c); }
	//! Returns a view matrix using position, target and approximate up vector
	CY_NODISCARD static Matrix4 View( Vec3<T> const &pos, Vec3<T> const &target, Vec3<T> const &up ) { Matrix4 m; m.SetView(pos,target,up); return m; }
	//! Returns a rotation matrix around x axis by angle in radians
	CY_NODISCARD static Matrix4 RotationX( T angle ) { Matrix4 m; m.SetRotationX(angle); return m; }
	//! Returns a rotation matrix around y axis by angle in radians
	CY_NODISCARD static Matrix4 RotationY( T angle ) { Matrix4 m; m.SetRotationY(angle); return m; }
	//! Returns a rotation matrix around z axis by angle in radians
	CY_NODISCARD static Matrix4 RotationZ( T angle ) { Matrix4 m; m.SetRotationZ(angle); return m; }
	//! Returns a rotation matrix about the given axis by angle in radians
	CY_NODISCARD static Matrix4 Rotation( Vec3<T> const &axis, T angle ) { Matrix4 m; m.SetRotation(axis,angle); return m; }
	//! Returns a rotation matrix about the given axis by cos and sin of the rotation angle
	CY_NODISCARD static Matrix4 Rotation( Vec3<T> const &axis, T cosAngle, T sinAngle ) { Matrix4 m; m.SetRotation(axis,cosAngle,sinAngle); return m; }
	//! Returns a rotation matrix that sets [from] unit vector to [to] unit vector
	CY_NODISCARD static Matrix4 Rotation( Vec3<T> const &from, Vec3<T> const &to ) { Matrix4 m; m.SetRotation(from,to); return m; }
	//! Returns a rotation matrix around x, y, and then z axes by angle in radians (Rz * Ry * Rx)
	CY_NODISCARD static Matrix4 RotationXYZ( T angleX, T angleY, T angleZ ) { Matrix4 m; m.SetRotationXYZ(angleX,angleY,angleZ); return m; }
	//! Returns a rotation matrix around z, y, and then x axes by angle in radians (Rx * Ry * Rz)
	CY_NODISCARD static Matrix4 RotationZYX( T angleX, T angleY, T angleZ ) { Matrix4 m; m.SetRotationZYX(angleX,angleY,angleZ); return m; }
	//! Returns a uniform scale matrix
	CY_NODISCARD static Matrix4 Scale( T uniformScale ) { Matrix4 m; m.SetScale(uniformScale); return m; }
	//! Returns a scale matrix
	CY_NODISCARD static Matrix4 Scale( T scaleX, T scaleY, T scaleZ, T scaleW=T(1) ) { Matrix4 m; m.SetScale(scaleX,scaleY,scaleZ,scaleW); return m; }
	//! Returns a scale matrix
	CY_NODISCARD static Matrix4 Scale( Vec3<T> const &scale ) { Matrix4 m; m.SetScale(scale); return m; }
	//! Returns a translation matrix with no rotation or scale
	CY_NODISCARD static Matrix4 Translation( Vec3<T> const &move ) { Matrix4 m; m.SetTranslation(move); return m; }
	//! Returns a project matrix with field of view in radians
	CY_NODISCARD static Matrix4 Perspective( T fov, T aspect, T znear, T zfar ) { Matrix4 m; m.SetPerspective(fov,aspect,znear,zfar); return m; }
	//! Returns a project matrix with the tangent of the half field of view (tan_fov_2)
	CY_NODISCARD static Matrix4 PerspectiveTan( T tan_fov_2, T aspect, T znear, T zfar ) { Matrix4 m; m.SetPerspectiveTan(tan_fov_2,aspect,znear,zfar); return m; }
	//! Returns the tensor product (outer product) matrix of two vectors
	CY_NODISCARD static Matrix4 TensorProduct( Vec4<T> const &v0, Vec4<T> const &v1 ) { Matrix4 m; m.SetTensorProduct(v0,v1); return m; }

	//////////////////////////////////////////////////////////////////////////
};

//-------------------------------------------------------------------------------

template<typename T> inline Matrix2<T> operator & ( Vec2<T> const &v0, Vec2<T> const &v1 ) { Matrix2<T> r; r.SetTensorProduct(v0,v1); return r; }	//!< tensor product (outer product) of two vectors
template<typename T> inline Matrix3<T> operator & ( Vec3<T> const &v0, Vec3<T> const &v1 ) { Matrix3<T> r; r.SetTensorProduct(v0,v1); return r; }	//!< tensor product (outer product) of two vectors
template<typename T> inline Matrix4<T> operator & ( Vec4<T> const &v0, Vec4<T> const &v1 ) { Matrix4<T> r; r.SetTensorProduct(v0,v1); return r; }	//!< tensor product (outer product) of two vectors

//-------------------------------------------------------------------------------

// Definitions of the conversion constructors
template <typename T>  Matrix2 <T>::Matrix2 ( Matrix3 <T> const &m ) { MemCopy(cell,m.cell,2); MemCopy(cell+2,m.cell+3,2); }
template <typename T>  Matrix2 <T>::Matrix2 ( Matrix34<T> const &m ) { MemCopy(cell,m.cell,2); MemCopy(cell+2,m.cell+3,2); }
template <typename T>  Matrix2 <T>::Matrix2 ( Matrix4 <T> const &m ) { MemCopy(cell,m.cell,2); MemCopy(cell+2,m.cell+4,2); }
template <typename T>  Matrix3 <T>::Matrix3 ( Matrix34<T> const &m ) { MemCopy(cell,m.cell,9); }
template <typename T>  Matrix3 <T>::Matrix3 ( Matrix4 <T> const &m ) { MemCopy(cell,m.cell,3); MemCopy(cell+3,m.cell+4,3); MemCopy(cell+6,m.cell+8,3); }
template <typename T>  Matrix34<T>::Matrix34( Matrix4 <T> const &m ) { MemCopy(cell,m.cell,3); MemCopy(cell+3,m.cell+4,3); MemCopy(cell+6,m.cell+8,3); MemCopy(cell+9,m.cell+12,3); }

template <typename T> inline Matrix4<T> Matrix34<T>::GetTranspose() const
{
	Matrix4<T> m;
	m.cell[ 0] = cell[0];   m.cell[ 1] = cell[3];   m.cell[ 2] = cell[ 6];   m.cell[ 3] = cell[ 9];
	m.cell[ 4] = cell[1];   m.cell[ 5] = cell[4];   m.cell[ 6] = cell[ 7];   m.cell[ 7] = cell[10];
	m.cell[ 8] = cell[2];   m.cell[ 9] = cell[5];   m.cell[10] = cell[ 8];   m.cell[11] = cell[11];
	m.cell[12] = T(0);      m.cell[13] = T(0);      m.cell[14] = T(0);       m.cell[15] = T(1);
	return m;
}

//-------------------------------------------------------------------------------

typedef Matrix2 <float>  Matrix2f;	//!< Single precision (float) 2x2 Matrix class
typedef Matrix3 <float>  Matrix3f;	//!< Single precision (float) 3x3 Matrix class
typedef Matrix34<float>  Matrix34f;	//!< Single precision (float) 3x4 Matrix class
typedef Matrix4 <float>  Matrix4f;	//!< Single precision (float) 4x4 Matrix class

typedef Matrix2 <double> Matrix2d;	//!< Double precision (double) 2x2 Matrix class
typedef Matrix3 <double> Matrix3d;	//!< Double precision (double) 3x3 Matrix class
typedef Matrix34<double> Matrix34d;	//!< Double precision (double) 3x4 Matrix class
typedef Matrix4 <double> Matrix4d;	//!< Double precision (double) 4x4 Matrix class

//-------------------------------------------------------------------------------
} // namespace hf
//-------------------------------------------------------------------------------

typedef cy::Matrix2f  cyMatrix2f;	//!< Single precision (float) 2x2 Matrix class
typedef cy::Matrix3f  cyMatrix3f;	//!< Single precision (float) 3x3 Matrix class
typedef cy::Matrix34f cyMatrix34f;	//!< Single precision (float) 3x4 Matrix class
typedef cy::Matrix4f  cyMatrix4f;	//!< Single precision (float) 4x4 Matrix class

typedef cy::Matrix2d  cyMatrix2d;	//!< Double precision (double) 2x2 Matrix class
typedef cy::Matrix3d  cyMatrix3d;	//!< Double precision (double) 3x3 Matrix class
typedef cy::Matrix34d cyMatrix34d;	//!< Double precision (double) 3x4 Matrix class
typedef cy::Matrix4d  cyMatrix4d;	//!< Double precision (double) 4x4 Matrix class

//-------------------------------------------------------------------------------

#endif
