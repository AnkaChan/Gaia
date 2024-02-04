// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyTriMesh.h 
//! \author Cem Yuksel
//! 
//! \brief  Triangular Mesh class.
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

#ifndef _CY_TRIMESH_H_INCLUDED_
#define _CY_TRIMESH_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyVector.h"
#include <vector>
#include <iostream>

//-------------------------------------------------------------------------------

_CY_CRT_SECURE_NO_WARNINGS

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

//! Triangular Mesh Class

class TriMesh
{
public:
	//! Triangular Mesh Face
	struct TriFace
	{
		unsigned int v[3];	//!< vertex indices
	};

	//! Simple character string
	struct Str
	{
		char *data;	//!< String data
		Str() : data(nullptr) {}							//!< Constructor
		Str( Str const &s ) : data(nullptr) { *this = s; }	//!< Copy constructor
		~Str() { if ( data ) delete [] data; }				//!< Destructor
		operator char const * () { return data; }			//!< Implicit conversion to const char
		void operator = ( Str  const &s ) { *this = s.data; }	//!< Assignment operator
		void operator = ( char const *s ) { if (s) { size_t n=strlen(s); if (data) delete [] data; data=new char[n+1]; strncpy(data,s,n); data[n]='\0'; } else if (data) { delete [] data; data=nullptr; } }	//!< Assignment operator
	};

	//! Material definition
	struct Mtl
	{
		Str   name;		//!< Material name
		float Ka[3];	//!< Ambient color
		float Kd[3];	//!< Diffuse color
		float Ks[3];	//!< Specular color
		float Tf[3];	//!< Transmission color
		float Ns;		//!< Specular exponent
		float Ni;		//!< Index of refraction
		int   illum;	//!< Illumination model
		Str   map_Ka;	//!< Ambient color texture map
		Str   map_Kd;	//!< Diffuse color texture map
		Str   map_Ks;	//!< Specular color texture map
		Str   map_Ns;	//!< Specular exponent texture map
		Str   map_d;	//!< Alpha texture map
		Str   map_bump;	//!< Bump texture map
		Str   map_disp;	//!< Displacement texture map

		//! Constructor sets the default material values
		Mtl()
		{
			Ka[0]=Ka[1]=Ka[2]=0;
			Kd[0]=Kd[1]=Kd[2]=1;
			Ks[0]=Ks[1]=Ks[2]=0;
			Tf[0]=Tf[1]=Tf[2]=0;
			Ns=0;
			Ni=1;
			illum=2;
		}
	};

protected:
	Vec3f   *v;		//!< vertices
	TriFace *f;		//!< faces
	Vec3f   *vn;	//!< vertex normal
	TriFace *fn;	//!< normal faces
	Vec3f   *vt;	//!< texture vertices
	TriFace *ft;	//!< texture faces
	Mtl     *m;		//!< materials
	int     *mcfc;	//!< material cumulative face count

	unsigned int nv;	//!< number of vertices
	unsigned int nf;	//!< number of faces
	unsigned int nvn;	//!< number of vertex normals
	unsigned int nvt;	//!< number of texture vertices
	unsigned int nm;	//!< number of materials

	Vec3f boundMin;	//!< Bounding box minimum bound
	Vec3f boundMax;	//!< Bounding box maximum bound

public:

	//!@name Constructors and Destructor
	TriMesh() : v(nullptr), f(nullptr), vn(nullptr), fn(nullptr), vt(nullptr), ft(nullptr), m(nullptr), mcfc(nullptr)
				, nv(0), nf(0), nvn(0), nvt(0), nm(0),boundMin(1,1,1), boundMax(0,0,0) {}
	TriMesh( TriMesh const &t ) : v(nullptr), f(nullptr), vn(nullptr), fn(nullptr), vt(nullptr), ft(nullptr), m(nullptr), mcfc(nullptr)
				, nv(0), nf(0), nvn(0), nvt(0), nm(0),boundMin(1,1,1), boundMax(0,0,0) { *this = t; }
	virtual ~TriMesh() { Clear(); }

	//!@name Component Access Methods
	Vec3f const &   V (int i) const { return v[i]; }		//!< returns the i^th vertex
	Vec3f&          V (int i)       { return v[i]; }		//!< returns the i^th vertex
	TriFace const & F (int i) const { return f[i]; }		//!< returns the i^th face
	TriFace&        F (int i)       { return f[i]; }		//!< returns the i^th face
	Vec3f const &   VN(int i) const { return vn[i]; }	//!< returns the i^th vertex normal
	Vec3f&          VN(int i)       { return vn[i]; }	//!< returns the i^th vertex normal
	TriFace const & FN(int i) const { return fn[i]; }	//!< returns the i^th normal face
	TriFace&        FN(int i)       { return fn[i]; }	//!< returns the i^th normal face
	Vec3f const &   VT(int i) const { return vt[i]; }	//!< returns the i^th vertex texture
	Vec3f&          VT(int i)       { return vt[i]; }	//!< returns the i^th vertex texture
	TriFace const & FT(int i) const { return ft[i]; }	//!< returns the i^th texture face
	TriFace&        FT(int i)       { return ft[i]; }	//!< returns the i^th texture face
	Mtl const &     M (int i) const { return m[i]; }		//!< returns the i^th material
	Mtl&            M (int i)       { return m[i]; }		//!< returns the i^th material

	unsigned int NV () const { return nv; }		//!< returns the number of vertices
	unsigned int NF () const { return nf; }		//!< returns the number of faces
	unsigned int NVN() const { return nvn; }	//!< returns the number of vertex normals
	unsigned int NVT() const { return nvt; }	//!< returns the number of texture vertices
	unsigned int NM () const { return nm; }		//!< returns the number of materials

	bool HasNormals() const { return NVN() > 0; }			//!< returns true if the mesh has vertex normals
	bool HasTextureVertices() const { return NVT() > 0; }	//!< returns true if the mesh has texture vertices

	//!@name Set Component Count
	void Clear() { SetNumVertex(0); SetNumFaces(0); SetNumNormals(0); SetNumTexVerts(0); SetNumMtls(0); boundMin.Set(1,1,1); boundMax.Zero(); }	//!< Deletes all components of the mesh
	void SetNumVertex  ( unsigned int n ) { Allocate(n,v,nv); }															//!< Sets the number of vertices and allocates memory for vertex positions
	void SetNumFaces   ( unsigned int n ) { Allocate(n,f,nf); if (fn||vn) Allocate(n,fn); if (ft||vt) Allocate(n,ft); }	//!< Sets the number of faces and allocates memory for face data. Normal faces and texture faces are also allocated, if they are used.
	void SetNumNormals ( unsigned int n ) { Allocate(n,vn,nvn); Allocate(n==0?0:nf,fn); }									//!< Sets the number of normals and allocates memory for normals and normal faces.
	void SetNumTexVerts( unsigned int n ) { Allocate(n,vt,nvt); Allocate(n==0?0:nf,ft); }									//!< Sets the number of texture coordinates and allocates memory for texture coordinates and texture faces.
	void SetNumMtls    ( unsigned int n ) { Allocate(n,m,nm); Allocate(n,mcfc); }											//!< Sets the number of materials and allocates memory for material data.
	void operator = ( TriMesh const &t );																					//!< Copies mesh data from the given mesh.

	//!@name Get Property Methods
	bool  IsBoundBoxReady() const { return boundMin.x<=boundMax.x && boundMin.y<=boundMax.y && boundMin.z<=boundMax.z; }	//!< Returns true if the bounding box has been computed.
	Vec3f GetBoundMin() const { return boundMin; }		//!< Returns the minimum values of the bounding box
	Vec3f GetBoundMax() const { return boundMax; }		//!< Returns the maximum values of the bounding box
	Vec3f GetVec     (int faceID, Vec3f const &bc) const { return Interpolate(faceID,v,f,bc); }		//!< Returns the point on the given face with the given barycentric coordinates (bc).
	Vec3f GetNormal  (int faceID, Vec3f const &bc) const { return Interpolate(faceID,vn,fn,bc); }	//!< Returns the the surface normal on the given face at the given barycentric coordinates (bc). The returned vector is not normalized.
	Vec3f GetTexCoord(int faceID, Vec3f const &bc) const { return Interpolate(faceID,vt,ft,bc); }	//!< Returns the texture coordinate on the given face at the given barycentric coordinates (bc).
	int   GetMaterialIndex(int faceID) const;				//!< Returns the material index of the face. This method goes through material counts of all materials to find the material index of the face. Returns a negative number if the face as no material
	int   GetMaterialFaceCount(int mtlID) const { return mtlID>0 ? mcfc[mtlID]-mcfc[mtlID-1] : mcfc[0]; }	//!< Returns the number of faces associated with the given material ID.
	int   GetMaterialFirstFace(int mtlID) const { return mtlID>0 ? mcfc[mtlID-1] : 0; }	//!< Returns the first face index associated with the given material ID. Other faces associated with the same material are placed are placed consecutively.

	//!@name Compute Methods
	void ComputeBoundingBox();						//!< Computes the bounding box
	void ComputeNormals(bool clockwise=false);		//!< Computes and stores vertex normals

	//!@name Load and Save methods
	bool LoadFromFileObj( char const *filename, bool loadMtl=true, std::ostream *outStream=&std::cout );	//!< Loads the mesh from an OBJ file. Automatically converts all faces to triangles.
	bool SaveToFileObj( char const *filename, std::ostream *outStream );									//!< Saves the mesh to an OBJ file with the given name.

private:
	template <class T> void Allocate( unsigned int n, T* &t ) { if (t) delete [] t; if (n>0) t = new T[n]; else t=nullptr; }
	template <class T> bool Allocate( unsigned int n, T* &t, unsigned int &nt ) { if (n==nt) return false; nt=n; Allocate(n,t); return true; }
	template <class T> void Copy( T const *from, unsigned int n, T* &t, unsigned int &nt) { if (!from) n=0; Allocate(n,t,nt); if (t) memcpy(t,from,sizeof(T)*n); }
	template <class T> void Copy( T const *from, unsigned int n, T* &t) { if (!from) n=0; Allocate(n,t); if (t) memcpy(t,from,sizeof(T)*n); }
	static Vec3f Interpolate( int i, Vec3f const *v, TriFace const *f, Vec3f const &bc ) { return v[f[i].v[0]]*bc.x + v[f[i].v[1]]*bc.y + v[f[i].v[2]]*bc.z; }

	// Temporary structures
	struct MtlData
	{
		std::string mtlName;
		unsigned int firstFace;
		unsigned int faceCount;
		MtlData() { faceCount=0; firstFace=0; }
	};
	struct MtlLibName { std::string filename; };
};

//-------------------------------------------------------------------------------

inline void TriMesh::operator = ( TriMesh const &t )
{
	Copy( t.v,  t.nv,  v,  nv  );
	Copy( t.f,  t.nf,  f,  nf  );
	Copy( t.vn, t.nvn, vn, nvn );
	Copy( t.fn, t.nf,  fn );
	Copy( t.vt, t.nvt, vt, nvt );
	Copy( t.ft, t.nf,  ft );
	Allocate(t.nm, m, nm);
	for ( unsigned int i=0; i<nm; i++ ) m[i] = t.m[i];
	Copy( t.mcfc, t.nm,  mcfc );
	boundMin = t.boundMin;
	boundMax = t.boundMax;
}

inline int TriMesh::GetMaterialIndex(int faceID) const
{
	for ( unsigned int i=0; i<nm; i++ ) {
		if ( faceID < mcfc[i] ) return (int) i;
	}
	return -1;
}

inline void TriMesh::ComputeBoundingBox()
{
	if ( nv > 0 ) {
		boundMin=v[0];
		boundMax=v[0];
		for ( unsigned int i=1; i<nv; i++ ) {
			if ( boundMin.x > v[i].x ) boundMin.x = v[i].x;
			if ( boundMin.y > v[i].y ) boundMin.y = v[i].y;
			if ( boundMin.z > v[i].z ) boundMin.z = v[i].z;
			if ( boundMax.x < v[i].x ) boundMax.x = v[i].x;
			if ( boundMax.y < v[i].y ) boundMax.y = v[i].y;
			if ( boundMax.z < v[i].z ) boundMax.z = v[i].z;
		}
	} else {
		boundMin.Set(1,1,1);
		boundMax.Set(0,0,0);
	}
}

inline void TriMesh::ComputeNormals(bool clockwise)
{
	SetNumNormals(nv);
	for ( unsigned int i=0; i<nvn; i++ ) vn[i].Set(0,0,0);	// initialize all normals to zero
	for ( unsigned int i=0; i<nf; i++ ) {
		Vec3f N = (v[f[i].v[1]]-v[f[i].v[0]]) ^ (v[f[i].v[2]]-v[f[i].v[0]]);	// face normal (not normalized)
		if ( clockwise ) N = -N;
		vn[f[i].v[0]] += N;
		vn[f[i].v[1]] += N;
		vn[f[i].v[2]] += N;
		fn[i] = f[i];
	}
	for ( unsigned int i=0; i<nvn; i++ ) vn[i].Normalize();
}

inline bool TriMesh::LoadFromFileObj( char const *filename, bool loadMtl, std::ostream *outStream )
{
	FILE *fp = fopen(filename,"r");
	if ( !fp ) {
		if ( outStream ) *outStream << "ERROR: Cannot open file " << filename << std::endl;
		return false;
	}

	Clear();

	class Buffer
	{
		char data[1024];
		int readLine;
	public:
		int ReadLine(FILE *fp)
		{
			char c = fgetc(fp);
			while ( !feof(fp) ) {
				while ( isspace(c) && ( !feof(fp) || c!='\0' ) ) c = fgetc(fp);	// skip empty space
				if ( c == '#' ) while ( !feof(fp) && c!='\n' && c!='\r' && c!='\0' ) c = fgetc(fp);	// skip comment line
				else break;
			}
			int i=0;
			bool inspace = false;
			while ( i<1024-1 ) {
				if ( feof(fp) || c=='\n' || c=='\r' || c=='\0' ) break;
				if ( isspace(c) ) {	// only use a single space as the space character
					inspace = true;
				} else {
					if ( inspace ) data[i++] = ' ';
					inspace = false;
					data[i++] = c;
				}
				c = fgetc(fp);
			}
			data[i] = '\0';
			readLine = i;
			return i;
		}
		char& operator[](int i) { return data[i]; }
		void ReadVertex( Vec3f &v ) const { v.Zero(); sscanf( data+2, "%f %f %f", &v.x, &v.y, &v.z ); }
		void ReadFloat3( float f[3] ) const { f[2]=f[1]=f[0]=0; int n = sscanf( data+2, "%f %f %f", &f[0], &f[1], &f[2] ); if ( n==1 ) f[2]=f[1]=f[0]; }
		void ReadFloat( float *f ) const { sscanf( data+2, "%f", f ); }
		void ReadInt( int *i, int start ) const { sscanf( data+start, "%d", i ); }
		bool IsCommand( char const *cmd ) const {
			int i=0;
			while ( cmd[i]!='\0' ) {
				if ( cmd[i] != data[i] ) return false;
				i++;
			}
			return (data[i]=='\0' || data[i]==' ');
		}
		char const * Data(int start=0) { return data+start; }
		void Copy( Str &str, int start=0 )
		{
			while ( data[start] != '\0' && data[start] <= ' ' ) start++;
			str = Data(start);
		}
	};
	Buffer buffer;


	struct MtlList {
		std::vector<MtlData> mtlData;
		int GetMtlIndex( char const *mtlName )
		{
			for ( unsigned int i=0; i<mtlData.size(); i++ ) {
				if ( mtlData[i].mtlName == mtlName ) return (int)i;
			}
			return -1;
		}
		int CreateMtl( char const *mtlName, unsigned int firstFace )
		{
			if ( mtlName[0] == '\0' ) return 0;
			int i = GetMtlIndex(mtlName);
			if ( i >= 0 ) return i;
			MtlData m;
			m.mtlName = mtlName;
			m.firstFace = firstFace;
			mtlData.push_back(m);
			return (int)mtlData.size()-1;
		}
	};
	MtlList mtlList;

	std::vector<Vec3f>      _v;		// vertices
	std::vector<TriFace>    _f;		// faces
	std::vector<Vec3f>      _vn;	// vertex normal
	std::vector<TriFace>    _fn;	// normal faces
	std::vector<Vec3f>      _vt;	// texture vertices
	std::vector<TriFace>    _ft;	// texture faces
	std::vector<MtlLibName> mtlFiles;
	std::vector<int>        faceMtlIndex;

	int currentMtlIndex = -1;
	bool hasTextures=false, hasNormals=false;

	while ( int rb = buffer.ReadLine(fp) ) {
		if ( buffer.IsCommand("v") ) {
			Vec3f vertex;
			buffer.ReadVertex(vertex);
			_v.push_back(vertex);
		}
		else if ( buffer.IsCommand("vt") ) {
			Vec3f texVert;
			buffer.ReadVertex(texVert);
			_vt.push_back(texVert);
			hasTextures = true;
		}
		else if ( buffer.IsCommand("vn") ) {
			Vec3f normal;
			buffer.ReadVertex(normal);
			_vn.push_back(normal);
			hasNormals = true;
		}
		else if ( buffer.IsCommand("f") ) {
			int facevert = -1;
			bool inspace = true;
			bool negative = false;
			int type = 0;
			unsigned int index;
			TriFace face, textureFace, normalFace;
			unsigned int nFacesBefore = (unsigned int)_f.size();
			for ( int i=2; i<rb; i++ ) {
				if ( buffer[i] == ' ' ) inspace = true;
				else {
					if ( inspace ) {
						inspace=false;
						negative = false;
						type=0;
						index=0;
						switch ( facevert ) {
							case -1:
								// initialize face
								face.v[0] = face.v[1] = face.v[2] = 0;
								textureFace.v[0] = textureFace.v[1] = textureFace.v[2] = 0;
								normalFace. v[0] = normalFace. v[1] = normalFace. v[2] = 0;
							case 0:
							case 1:
								facevert++;
								break;
							case 2:
								// copy the first two vertices from the previous face
								_f.push_back(face);
								face.v[1] = face.v[2];
								if ( hasTextures ) {
									_ft.push_back(textureFace);
									textureFace.v[1] = textureFace.v[2];
								}
								if ( hasNormals ) {
									_fn.push_back(normalFace);
									normalFace.v[1] = normalFace.v[2];
								}
								faceMtlIndex.push_back(currentMtlIndex);
								break;
						}
					}
					if ( buffer[i] == '/' ) { type++; index=0; }
					if ( buffer[i] == '-' ) negative = true;
					if ( buffer[i] >= '0' && buffer[i] <= '9' ) {
						index = index*10 + (buffer[i]-'0');
						switch ( type ) {
							case 0: face.v       [facevert] = negative ? (unsigned int)_v. size()-index : index-1; break;
							case 1: textureFace.v[facevert] = negative ? (unsigned int)_vt.size()-index : index-1; hasTextures=true; break;
							case 2: normalFace.v [facevert] = negative ? (unsigned int)_vn.size()-index : index-1; hasNormals =true; break;
						}
					}
				}
			}
			_f.push_back(face);
			if ( hasTextures ) _ft.push_back(textureFace);
			if ( hasNormals  ) _fn.push_back(normalFace);
			faceMtlIndex.push_back(currentMtlIndex);
			if ( currentMtlIndex>=0 ) mtlList.mtlData[currentMtlIndex].faceCount += (unsigned int)_f.size() - nFacesBefore;
		}
		else if ( loadMtl ) {
			if ( buffer.IsCommand("usemtl") ) {
				currentMtlIndex = mtlList.CreateMtl(buffer.Data(7), (unsigned int)_f.size());
			}
			if ( buffer.IsCommand("mtllib") ) {
				MtlLibName libName;
				libName.filename = buffer.Data(7);
				mtlFiles.push_back(libName);
			}
		}
		if ( feof(fp) ) break;
	}

	fclose(fp);


	if ( _f.size() == 0 ) return true; // No faces found
	SetNumVertex((unsigned int)_v.size());
	SetNumFaces((unsigned int)_f.size());
	SetNumTexVerts((unsigned int)_vt.size());
	SetNumNormals((unsigned int)_vn.size());
	if ( loadMtl ) SetNumMtls((unsigned int)mtlList.mtlData.size());

	// Copy data
	memcpy(v, _v.data(), sizeof(Vec3f)*_v.size());
	if ( _vt.size() > 0 ) memcpy(vt, _vt.data(), sizeof(Vec3f)*_vt.size());
	if ( _vn.size() > 0 ) memcpy(vn, _vn.data(), sizeof(Vec3f)*_vn.size());

	if ( mtlList.mtlData.size() > 0 ) {
		unsigned int fid = 0;
		for ( int m=0; m<(int)mtlList.mtlData.size(); m++ ) {
			for ( unsigned int i=mtlList.mtlData[m].firstFace, j=0; j<mtlList.mtlData[m].faceCount && i<_f.size(); i++ ) {
				if ( faceMtlIndex[i] == m ) {
					f[fid] = _f[i];
					if ( fn ) fn[fid] = _fn[i];
					if ( ft ) ft[fid] = _ft[i];
					fid++;
					j++;
				}
			}
			mcfc[m] = fid;
		}
		if ( fid <_f.size() ) {
			for ( unsigned int i=0; i<_f.size(); i++ ) {
				if ( faceMtlIndex[i] < 0 ) {
					f[fid] = _f[i];
					if ( fn ) fn[fid] = _fn[i];
					if ( ft ) ft[fid] = _ft[i];
					fid++;
				}
			}
		}
	} else {
		memcpy(f, _f.data(), sizeof(TriFace)*_f.size());
		if ( ft ) memcpy(ft, _ft.data(), sizeof(TriFace)*_ft.size());
		if ( fn ) memcpy(fn, _fn.data(), sizeof(TriFace)*_fn.size());
	}


	// Load the .mtl files
	if ( loadMtl ) {
		// get the path from filename
		char *mtlPathName = nullptr;
		char const *pathEnd = strrchr(filename,'\\');
		if ( !pathEnd ) pathEnd = strrchr(filename,'/');
		if ( pathEnd ) {
			int n = int(pathEnd-filename) + 1;
			mtlPathName = new char[n+1];
			strncpy(mtlPathName,filename,n);
			mtlPathName[n] = '\0';
		}
		for ( unsigned int mi=0; mi<mtlFiles.size(); mi++ ) {
			std::string mtlFilename = ( mtlPathName ) ? std::string(mtlPathName) + mtlFiles[mi].filename : mtlFiles[mi].filename;
			FILE *fp = fopen(mtlFilename.data(),"r");
			if ( !fp ) {
				if ( outStream ) *outStream << "ERROR: Cannot open file " << mtlFilename.c_str() << std::endl;
				continue;
			}
			int mtlID = -1;
			while ( buffer.ReadLine(fp) ) {
				if ( buffer.IsCommand("newmtl") ) {
					mtlID = mtlList.GetMtlIndex(buffer.Data(7));
					if ( mtlID >= 0 ) buffer.Copy( m[mtlID].name, 7 );
				} else if ( mtlID >= 0 ) {
					if ( buffer.IsCommand("Ka") ) buffer.ReadFloat3( m[mtlID].Ka );
					else if ( buffer.IsCommand("Kd") ) buffer.ReadFloat3( m[mtlID].Kd );
					else if ( buffer.IsCommand("Ks") ) buffer.ReadFloat3( m[mtlID].Ks );
					else if ( buffer.IsCommand("Tf") ) buffer.ReadFloat3( m[mtlID].Tf );
					else if ( buffer.IsCommand("Ns") ) buffer.ReadFloat( &m[mtlID].Ns );
					else if ( buffer.IsCommand("Ni") ) buffer.ReadFloat( &m[mtlID].Ni );
					else if ( buffer.IsCommand("illum") ) buffer.ReadInt( &m[mtlID].illum, 5 );
					else if ( buffer.IsCommand("map_Ka"  ) ) buffer.Copy( m[mtlID].map_Ka,   7 );
					else if ( buffer.IsCommand("map_Kd"  ) ) buffer.Copy( m[mtlID].map_Kd,   7 );
					else if ( buffer.IsCommand("map_Ks"  ) ) buffer.Copy( m[mtlID].map_Ks,   7 );
					else if ( buffer.IsCommand("map_Ns"  ) ) buffer.Copy( m[mtlID].map_Ns,   7 );
					else if ( buffer.IsCommand("map_d"   ) ) buffer.Copy( m[mtlID].map_d,    6 );
					else if ( buffer.IsCommand("map_bump") ) buffer.Copy( m[mtlID].map_bump, 9 );
					else if ( buffer.IsCommand("bump"    ) ) buffer.Copy( m[mtlID].map_bump, 5 );
					else if ( buffer.IsCommand("map_disp") ) buffer.Copy( m[mtlID].map_disp, 9 );
					else if ( buffer.IsCommand("disp"    ) ) buffer.Copy( m[mtlID].map_disp, 5 );
				}
			}
			fclose(fp);
		}
		if ( mtlPathName ) delete [] mtlPathName;
	}

	return true;
}

//-------------------------------------------------------------------------------

inline bool TriMesh::SaveToFileObj( char const *filename, std::ostream *outStream )
{
	FILE *fp = fopen(filename,"w");
	if ( !fp ) {
		if ( outStream ) *outStream << "ERROR: Cannot create file " << filename << std::endl;
		return false;
	}

	for ( unsigned int i=0; i<nv; i++ ) {
		fprintf(fp,"v %f %f %f\n",v[i].x, v[i].y, v[i].z);
	}
	for ( unsigned int i=0; i<nvt; i++ ) {
		fprintf(fp,"vt %f %f %f\n",vt[i].x, vt[i].y, vt[i].z);
	}
	for ( unsigned int i=0; i<nvn; i++ ) {
		fprintf(fp,"vn %f %f %f\n",vn[i].x, vn[i].y, vn[i].z);
	}
	int faceFormat = ((nvn>0)<<1) | (nvt>0);
	switch ( faceFormat ) {
	case 0:
		for ( unsigned int i=0; i<nf; i++ ) {
			fprintf(fp,"f %d %d %d\n", f[i].v[0]+1, f[i].v[1]+1, f[i].v[2]+1);
		}
		break;
	case 1:
		for ( unsigned int i=0; i<nf; i++ ) {
			fprintf(fp,"f %d/%d %d/%d %d/%d\n", f[i].v[0]+1, ft[i].v[0]+1, f[i].v[1]+1, ft[i].v[1]+1, f[i].v[2]+1, ft[i].v[2]+1);
		}
		break;
	case 2:
		for ( unsigned int i=0; i<nf; i++ ) {
			fprintf(fp,"f %d//%d %d//%d %d//%d\n", f[i].v[0]+1, fn[i].v[0]+1, f[i].v[1]+1, fn[i].v[1]+1, f[i].v[2]+1, fn[i].v[2]+1);
		}
		break;
	case 3:
		for ( unsigned int i=0; i<nf; i++ ) {
			fprintf(fp,"f %d/%d/%d %d/%d/%d %d/%d/%d\n", f[i].v[0]+1, ft[i].v[0]+1, fn[i].v[0]+1, f[i].v[1]+1, ft[i].v[1]+1, fn[i].v[1]+1, f[i].v[2]+1, ft[i].v[2]+1, fn[i].v[2]+1);
		}
		break;
	}

	fclose(fp);

	return true;
}

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::TriMesh cyTriMesh;	//!< Triangular Mesh Class

//-------------------------------------------------------------------------------

_CY_CRT_SECURE_RESUME_WARNINGS
#endif

