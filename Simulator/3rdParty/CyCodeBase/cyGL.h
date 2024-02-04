// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyGL.h 
//! \author Cem Yuksel
//! 
//! \brief  OpenGL helper classes
//! 
//! The classes in this file are designed to provide convenient interfaces for
//! some OpenGL features. They are not intended to provide the full flexibility
//! of the underlying OpenGL functions, but they greatly simplify the 
//! implementation of some general OpenGL operations.
//!
//-------------------------------------------------------------------------------
//
// Copyright (c) 2017, Cem Yuksel <cem@cemyuksel.com>
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

#ifndef _CY_GL_H_INCLUDED_
#define _CY_GL_H_INCLUDED_

//-------------------------------------------------------------------------------
#if !defined(__gl_h_) && !defined(__GL_H__) && !defined(_GL_H) && !defined(__X_GL_H)
#error gl.h not included before cyGL.h
#endif
#ifndef GL_VERSION_2_0
#error OpenGL 2.0 extensions are required for cyGL.h. You must include an OpenGL extensions header before including cyGL.h.
#endif
#ifndef GL_VERSION_3_0
# define _CY_GL_VERSION_3_0_WARNING "OpenGL version 3 is required for using geometry and tessellation shaders, but the OpenGL extensions header included before cyGL.h does not include OpenGL version 3.0 definitions."
# if _MSC_VER
#  pragma message ("Warning: " _CY_GL_VERSION_3_0_WARNING)
# elif   __GNUC__
#  warning (_CY_GL_VERSION_3_0_WARNING)
# endif
#endif
//-------------------------------------------------------------------------------
#ifndef GL_GEOMETRY_SHADER
#define GL_GEOMETRY_SHADER 0x8DD9
#endif
#ifndef GL_TESS_EVALUATION_SHADER
#define GL_TESS_EVALUATION_SHADER 0x8E87
#endif
#ifndef GL_TESS_CONTROL_SHADER
#define GL_TESS_CONTROL_SHADER 0x8E88
#endif
#ifdef APIENTRY
# define _CY_APIENTRY APIENTRY
#else
# if defined(__MINGW32__) || defined(__CYGWIN__) || (_MSC_VER >= 800) || defined(_STDCALL_SUPPORTED) || defined(__BORLANDC__)
#  define _CY_APIENTRY __stdcall
# else
#  define _CY_APIENTRY
# endif
#endif
//-------------------------------------------------------------------------------

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

//-------------------------------------------------------------------------------

// These includes are needed for checking the OpenGL context.
// If CY_GL_DONT_CHECK_CONTEXT is defined before including this file,
// checking the OpenGL context is disabled and these includes are skipped.
#ifndef CY_GL_DONT_CHECK_CONTEXT
# ifdef _WIN32
#  include <wtypes.h>
#  include <Wingdi.h>
#  define _CY_GL_GET_CONTEXT wglGetCurrentContext()
# elif defined(__APPLE__)
#  include <OpenGL/OpenGL.h>
#  define _CY_GL_GET_CONTEXT CGLGetCurrentContext()
# elif defined(__unix__)
#  include <GL/glx.h>
#  define _CY_GL_GET_CONTEXT glXGetCurrentContext()
# else
#  define _CY_GL_GET_CONTEXT 1
# endif
#else
# define _CY_GL_GET_CONTEXT 1
#endif

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

#define CY_GL_INVALID_ID 0xFFFFFFFF	//!< Invalid ID value

//-------------------------------------------------------------------------------

//! General OpenGL queries
//!
//! This class includes static functions for general OpenGL queries
//! and types used by the other classes in this file.
class GL
{
public:
	//! OpenGL types
	enum Type {
		TYPE_UBYTE = 0,
		TYPE_USHORT,
		TYPE_HALF,
		TYPE_FLOAT,
		TYPE_INT8,
		TYPE_UINT8,
		TYPE_INT16,
		TYPE_UINT16,
		TYPE_INT32,
		TYPE_UINT32,
	};

	static GLenum GetGLType( GLubyte  const * ) { return GL_UNSIGNED_BYTE ; }	//!< Returns the OpenGL type identifier that corresponds to unsigned byte.
	static GLenum GetGLType( GLushort const * ) { return GL_UNSIGNED_SHORT; }	//!< Returns the OpenGL type identifier that corresponds to unsigned short.
	static GLenum GetGLType( GLfloat  const * ) { return GL_FLOAT         ; }	//!< Returns the OpenGL type identifier that corresponds to float.
	static GLenum GetGLType( GLbyte   const * ) { return GL_BYTE          ; }	//!< Returns the OpenGL type identifier that corresponds to byte.
	static GLenum GetGLType( GLshort  const * ) { return GL_SHORT         ; }	//!< Returns the OpenGL type identifier that corresponds to short.
	static GLenum GetGLType( GLint    const * ) { return GL_INT           ; }	//!< Returns the OpenGL type identifier that corresponds to int.
	static GLenum GetGLType( GLuint   const * ) { return GL_UNSIGNED_INT  ; }	//!< Returns the OpenGL type identifier that corresponds to unsigned int.

	static Type GetType( GLubyte  const * ) { return TYPE_UBYTE ; }				//!< Returns the Type that corresponds to unsigned byte.
	static Type GetType( GLushort const * ) { return TYPE_USHORT; }				//!< Returns the Type that corresponds to unsigned short.
	static Type GetType( GLfloat  const * ) { return TYPE_FLOAT ; }				//!< Returns the Type that corresponds to float.
	static Type GetType( GLbyte   const * ) { return TYPE_INT8  ; }				//!< Returns the Type that corresponds to byte.
	static Type GetType( GLshort  const * ) { return TYPE_INT16 ; }				//!< Returns the Type that corresponds to short.
	static Type GetType( GLint    const * ) { return TYPE_INT32 ; }				//!< Returns the Type that corresponds to int.
	static Type GetType( GLuint   const * ) { return TYPE_UINT32; }				//!< Returns the Type that corresponds to unsigned int.

	static GLenum TextureFormat    ( Type type, int numChannels );				//!< Returns the internal OpenGL texture type identifier for the given Type and the number of channels.
	static GLenum TextureDataFormat( Type type, int numChannels );				//!< Returns the OpenGL texture data format identifier for the given Type and the number of channels.

	//! Prints the OpenGL version to the given stream.
	static void PrintVersion(std::ostream *outStream=&std::cout);

	//! Checks all previously triggered OpenGL errors and prints them to the given output stream.
	static void CheckError( char const *sourcefile, int line, char const *call=nullptr, std::ostream *outStream=&std::cout );

	//! Checks if an OpenGL context exists. Returns false if a valid OpenGL context cannot be retrieved.
	//! This is mostly useful for safely deleting previously allocated OpenGL objects.
	static bool CheckContext() { return _CY_GL_GET_CONTEXT ? true : false; }
};

//-------------------------------------------------------------------------------

//! Checks and prints OpenGL error messages to the default output stream.
#define CY_GL_ERROR _CY_GL_ERROR
#define _CY_GL_ERROR cy::GL::CheckError(__FILE__,__LINE__)

//! Checks OpenGL errors before calling the given function,
//! calls the function, and then checks OpenGL errors again.
//! If an error is found, it is printed to the default output stream.
#define CY_GL_ERR(gl_function_call) _CY_GL_ERR(gl_function_call)
#define _CY_GL_ERR(f) cy::GL::CheckError(__FILE__,__LINE__,"a prior call"); f; cy::GL::CheckError(__FILE__,__LINE__,#f)

#ifdef _DEBUG
# define CY_GL_ERROR_D CY_GL_ERROR //!< Checks and prints OpenGL error messages using CY_GL_ERROR in code compiled in debug mode only (with _DEBUG defined).
# define CY_GL_ERR_D   CY_GL_ERR   //!< Checks and prints OpenGL error messages using CY_GL_ERR in code compiled in debug mode only (with _DEBUG defined).
#else
# define CY_GL_ERROR_D			   //!< Checks and prints OpenGL error messages using CY_GL_ERROR in code compiled in debug mode only (with _DEBUG defined).
# define CY_GL_ERR_D			   //!< Checks and prints OpenGL error messages using CY_GL_ERR in code compiled in debug mode only (with _DEBUG defined).
#endif

//-------------------------------------------------------------------------------

#ifdef GL_KHR_debug
#define _CY_GLDebugCallback

//! OpenGL debug callback class.
//!
//! For this class to work, you may need to initialize the OpenGL context in debug mode.
//! This class registers an OpenGL debug callback function, which is called when
//! there is an OpenGL generates a debug message. 
//! The class has no local storage, so it can be safely deleted.
//! Deleting an object of this class, however, does not automatically disable the debug callbacks.
class GLDebugCallback
{
public:
	//! Constructor can register the callback, but only if the OpenGL context is created
	//! before the constructor is called. If the regsiterCallback argument is false,
	//! the other arguments are ignored.
	GLDebugCallback(bool registerCallback=false, bool ignoreNotifications=false, std::ostream *outStream=&std::cout)
		{ if ( registerCallback ) { Register(outStream); IgnoreNotifications(ignoreNotifications); } }

	//! Registers the debug callback function.
	//! The callback function outputs the debug data to the given stream.
	//! Note that if the OpenGL context is not created in debug mode, 
	//! OpenGL may not call the callback function.
	//! If there is a previously registered callback function,
	//! calling this function overwrites the previous callback registration.
	void Register(std::ostream *outStream=&std::cout) { glEnable(GL_DEBUG_OUTPUT); glDebugMessageCallback((GLDEBUGPROC)Callback,outStream); }

	//! Unregisters the OpenGL debug callback function.
	void Unregister() { glDisable(GL_DEBUG_OUTPUT); glDebugMessageCallback(0,0); }

	//! Sets which type of non-critical debug messages should be ignored.
	//! By default, no debug message type is ignored.
	void SetIgnoredTypes( bool deprecated_behavior, bool portability, bool performance, bool other )
	{
		glDebugMessageControl( GL_DONT_CARE, GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR, GL_DONT_CARE, 0, 0, deprecated_behavior );
		glDebugMessageControl( GL_DONT_CARE, GL_DEBUG_TYPE_PORTABILITY,         GL_DONT_CARE, 0, 0, portability );
		glDebugMessageControl( GL_DONT_CARE, GL_DEBUG_TYPE_PERFORMANCE,         GL_DONT_CARE, 0, 0, performance );
		glDebugMessageControl( GL_DONT_CARE, GL_DEBUG_TYPE_OTHER,               GL_DONT_CARE, 0, 0, other );
	}

	//! Sets weather notification messages should be ignored.
	void IgnoreNotifications(bool ignore=true) { glDebugMessageControl( GL_DONT_CARE, GL_DONT_CARE, GL_DEBUG_SEVERITY_NOTIFICATION, 0, 0, !ignore ); }

protected:
	//! This callback function that is called by OpenGL whenever there is a debug message.
	//! See the OpenGL documentation for glDebugMessageCallback for details.
	//! Placing the break point in this function allows easily identifying the
	//! OpenGL call that triggered the debug message (using the call stack).
	static void _CY_APIENTRY Callback( GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, GLchar const *message, void const *userParam );
};

//! Registers the OpenGL callback by ignoring notifications.
//! After this macro is called, the debug messages get printed to the default output stream.
//! The OpenGL context must be created before using this macro.
//! Note that OpenGL may not generate debug messages if it is not created in debug mode.
#define CY_GL_REGISTER_DEBUG_CALLBACK _CY_GL_REGISTER_DEBUG_CALLBACK
#define _CY_GL_REGISTER_DEBUG_CALLBACK { cy::GLDebugCallback callback(true,true); }

#endif // GL_KHR_debug

//-------------------------------------------------------------------------------

//! OpenGL texture base class.
//!
//! This class provides a convenient interface for handling basic texture
//! operations with OpenGL. The template argument TEXTURE_TYPE should be a
//! texture type supported by OpenGL, such as GL_TEXTURE_1D, GL_TEXTURE_2D,
//! GL_TEXTURE_3D, GL_TEXTURE_CUBE_MAP, GL_TEXTURE_RECTANGLE,
//! GL_TEXTURE_1D_ARRAY, GL_TEXTURE_2D_ARRAY, GL_TEXTURE_CUBE_MAP_ARRAY,
//! GL_TEXTURE_BUFFER, GL_TEXTURE_2D_MULTISAMPLE, or 
//! GL_TEXTURE_2D_MULTISAMPLE_ARRAY. This class merely stores the texture id.
//! Note that deleting an object of this class does not automatically delete
//! the texture from the GPU memory. You must explicitly call the Delete() 
//! method to free the texture storage on the GPU.
template <GLenum TEXTURE_TYPE>
class GLTexture
{
private:
	GLuint textureID;	//!< The texture ID

public:
	GLTexture() : textureID(CY_GL_INVALID_ID) {}	//!< Constructor.

	//!@name General Methods

	void   Delete() { if ( textureID != CY_GL_INVALID_ID ) glDeleteTextures(1,&textureID); textureID = CY_GL_INVALID_ID; }	//!< Deletes the texture.
	GLuint GetID () const { return textureID; }													//!< Returns the texture ID.
	bool   IsNull() const { return textureID == CY_GL_INVALID_ID; }								//!< Returns true if the OpenGL texture object is not generated, i.e. the texture id is invalid.
	void   Bind  () const { glBindTexture(TEXTURE_TYPE, textureID); }							//!< Binds the texture to the current texture unit.
	void   Bind  (int textureUnit) const { glActiveTexture(GL_TEXTURE0+textureUnit); Bind(); }	//!< Binds the texture to the given texture unit.
	GLenum Type  () const { return TEXTURE_TYPE; }

	//!@name Texture Creation and Initialization

	//! Generates the texture, only if the texture has not been previously generated.
	//! Initializes the texture sampling parameters
	void Initialize();

#ifdef GL_VERSION_3_0
	//! Builds mipmap levels for the texture. The texture image must be set first.
	void BuildMipmaps() { Bind(); glGenerateMipmap(TEXTURE_TYPE); }
#endif

	//! Sets the texture filtering mode.
	//! The acceptable values are GL_NEAREST and GL_LINEAR.
	//! The minification filter values can also be GL_NEAREST_MIPMAP_NEAREST, GL_LINEAR_MIPMAP_NEAREST, GL_NEAREST_MIPMAP_LINEAR, or GL_LINEAR_MIPMAP_LINEAR.
	//! If the filter argument is zero, the corresponding filter parameter is not changed.
	void SetFilteringMode(GLenum magnificationFilter, GLenum minificationFilter);

#ifdef GL_EXT_texture_filter_anisotropic
	//! Sets the anisotropy level of the texture.
	//! The anisotropy value of 1 disables anisotropic filtering.
	//! Larger values provide provide better anisotropic filtering.
	void SetAnisotropy(float anisotropy) { Bind(); glTexParameterf(TEXTURE_TYPE, GL_TEXTURE_MAX_ANISOTROPY_EXT, anisotropy ); }

	//! Sets anisotropic filtering to the maximum permissible value.
	void SetMaxAnisotropy() { float largest; glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &largest); SetAnisotropy(largest); }

	//! Turns off anisotropic filtering.
	void SetNoAnisotropy() { SetAnisotropy(1.0f); }
#endif
};

//-------------------------------------------------------------------------------

//! OpenGL 1D texture class.
//!
//! This class provides a convenient interface for handling 1D texture
//! operations with OpenGL. The template argument TEXTURE_TYPE should be a
//! 1D texture type supported by OpenGL, such as GL_TEXTURE_1D.
//! This class merely stores the texture id.
//! Note that deleting an object of this class does not automatically delete
//! the texture from the GPU memory. You must explicitly call the Delete() 
//! method to free the texture storage on the GPU.
template <GLenum TEXTURE_TYPE>
class GLTexture1 : public GLTexture<TEXTURE_TYPE>
{
public:
	//! Sets the texture image using the given texture format, data format, and data type.
	void SetImage( GLenum textureFormat, GLenum dataFormat, GLenum dataType, void const *data, GLsizei width, int level=0 ) { GLTexture<TEXTURE_TYPE>::Bind(); glTexImage1D(TEXTURE_TYPE,level,textureFormat,width,0,dataFormat,dataType,data); }

	//! Sets the texture image using the given texture format and data format. The data type is determined by the data pointer type.
	template <typename T> void SetImage( GLenum textureFormat, GLenum dataFormat, T const *data, GLsizei width, int level=0 ) { SetImage(textureFormat,dataFormat,GL::GetGLType(data),data,width,level); }

	//! Sets the texture image using the given texture type. The data format and type are determined by the data pointer type and the number of channels.
	//! The texture format is determined by the texture type and the number of channels.
	template <typename T> void SetImage( GL::Type textureType, T const *data, int numChannels, GLsizei width, int level=0 ) { SetImage(GL::TextureFormat(textureType,numChannels),GL::TextureDataFormat(GL::GetType(data),numChannels),data,width,level); }

	//! Sets the texture image. The texture format uses the matching data pointer type.
	//! If unsigned char is used, the texture uses 8-bit normalized values.
	//! If unsigned short is used, the texture uses 16-bit normalized values.
	//! If float is used, the texture uses non-normalized 32-bit float values.
	//! If char, short, or int is used, the texture uses non-normalized 8-bit, 16-bit, or 32-bit integer values.
	template <typename T> void SetImage( T const *data, int numChannels, GLsizei width, int level=0 ) { SetImage(GL::GetType(data),data,numChannels,width,level); }

	template <typename T> void SetImageRGBA( GL::Type textureType, T const *data, GLsizei width, int level=0 ) { SetImage(textureType,data,4,width,level); }	//!< Sets the texture image with 4 channels.
	template <typename T> void SetImageRGB ( GL::Type textureType, T const *data, GLsizei width, int level=0 ) { SetImage(textureType,data,3,width,level); }	//!< Sets the texture image with 3 channels.
	template <typename T> void SetImageRG  ( GL::Type textureType, T const *data, GLsizei width, int level=0 ) { SetImage(textureType,data,2,width,level); }	//!< Sets the texture image with 2 channels.
	template <typename T> void SetImageR   ( GL::Type textureType, T const *data, GLsizei width, int level=0 ) { SetImage(textureType,data,1,width,level); }	//!< Sets the texture image with 1 channel.

	template <typename T> void SetImageRGBA( T const *data, GLsizei width, int level=0 ) { SetImage(data,4,width,level); }	//!< Sets the texture image with 4 channels.
	template <typename T> void SetImageRGB ( T const *data, GLsizei width, int level=0 ) { SetImage(data,3,width,level); }	//!< Sets the texture image with 3 channels.
	template <typename T> void SetImageRG  ( T const *data, GLsizei width, int level=0 ) { SetImage(data,2,width,level); }	//!< Sets the texture image with 2 channels.
	template <typename T> void SetImageR   ( T const *data, GLsizei width, int level=0 ) { SetImage(data,1,width,level); }	//!< Sets the texture image with 1 channel.

	//! Sets the texture wrapping parameter.
	//! The acceptable values are GL_REPEAT, GL_MIRRORED_REPEAT, GL_CLAMP, and GL_CLAMP_TO_BORDER.
	void SetWrappingMode(GLenum wrapS) { GLTexture<TEXTURE_TYPE>::Bind(); glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_WRAP_S, wrapS); }

};

//-------------------------------------------------------------------------------

//! OpenGL 2D texture class.
//!
//! This class provides a convenient interface for handling 2D texture
//! operations with OpenGL. The template argument TEXTURE_TYPE should be a
//! 2D texture type supported by OpenGL, such as GL_TEXTURE_2D, 
//! GL_TEXTURE_RECTANGLE, or GL_TEXTURE_1D_ARRAY. This class merely stores the texture id.
//! Note that deleting an object of this class does not automatically delete
//! the texture from the GPU memory. You must explicitly call the Delete() 
//! method to free the texture storage on the GPU.
template <GLenum TEXTURE_TYPE>
class GLTexture2 : public GLTexture<TEXTURE_TYPE>
{
public:
	//! Sets the texture image using the given texture format, data format, and data type.
	void SetImage( GLenum textureFormat, GLenum dataFormat, GLenum dataType, void const *data, GLsizei width, GLsizei height, int level=0 ) { GLTexture<TEXTURE_TYPE>::Bind(); glTexImage2D(TEXTURE_TYPE,level,textureFormat,width,height,0,dataFormat,dataType,data); }

	//! Sets the texture image using the given texture format and data format. The data type is determined by the data pointer type.
	template <typename T> void SetImage( GLenum textureFormat, GLenum dataFormat, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(textureFormat,dataFormat,GL::GetGLType(data),data,width,height,level); }

	//! Sets the texture image using the given texture type. The data format and type are determined by the data pointer type and the number of channels.
	//! The texture format is determined by the texture type and the number of channels.
	template <typename T> void SetImage( GL::Type textureType, T const *data, int numChannels, GLsizei width, GLsizei height, int level=0 ) { SetImage(GL::TextureFormat(textureType,numChannels),GL::TextureDataFormat(GL::GetType(data),numChannels),data,width,height,level); }

	//! Sets the texture image. The texture format uses the matching data pointer type.
	//! If GLubyte (unsigned char) is used, the texture uses 8-bit normalized values.
	//! If GLushort (unsigned short) is used, the texture uses 16-bit normalized values.
	//! If GLfloat (float) is used, the texture uses non-normalized 32-bit float values.
	//! If GLbyte, GLshort, GLint, or GLuint is used, the texture uses non-normalized 8-bit, 16-bit, or 32-bit integer values.
	template <typename T> void SetImage( T const *data, int numChannels, GLsizei width, GLsizei height, int level=0 ) { SetImage(GL::GetType(data),data,numChannels,width,height,level); }

	template <typename T> void SetImageRGBA( GL::Type textureType, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(textureType,data,4,width,height,level); }	//!< Sets the texture image with 4 channels.
	template <typename T> void SetImageRGB ( GL::Type textureType, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(textureType,data,3,width,height,level); }	//!< Sets the texture image with 3 channels.
	template <typename T> void SetImageRG  ( GL::Type textureType, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(textureType,data,2,width,height,level); }	//!< Sets the texture image with 2 channels.
	template <typename T> void SetImageR   ( GL::Type textureType, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(textureType,data,1,width,height,level); }	//!< Sets the texture image with 1 channel.

	template <typename T> void SetImageRGBA( T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(data,4,width,height,level); }	//!< Sets the texture image with 4 channels.
	template <typename T> void SetImageRGB ( T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(data,3,width,height,level); }	//!< Sets the texture image with 3 channels.
	template <typename T> void SetImageRG  ( T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(data,2,width,height,level); }	//!< Sets the texture image with 2 channels.
	template <typename T> void SetImageR   ( T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(data,1,width,height,level); }	//!< Sets the texture image with 1 channel.

	//! Sets the texture wrapping parameter.
	//! The acceptable values are GL_REPEAT, GL_MIRRORED_REPEAT, GL_CLAMP, and GL_CLAMP_TO_BORDER.
	//! If the wrap argument is zero, the corresponding wrapping parameter is not changed.
	void SetWrappingMode(GLenum wrapS, GLenum wrapT)
	{
		GLTexture<TEXTURE_TYPE>::Bind();
		if ( wrapS != 0 ) glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_WRAP_S, wrapS);
		if ( wrapT != 0 ) glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_WRAP_T, wrapT);
	}

};

//-------------------------------------------------------------------------------

//! OpenGL 3D texture class.
//!
//! This class provides a convenient interface for handling 3D texture
//! operations with OpenGL. The template argument TEXTURE_TYPE should be a
//! 3D texture type supported by OpenGL, such as GL_TEXTURE_3D or
//! GL_TEXTURE_2D_ARRAY. This class merely stores the texture id.
//! Note that deleting an object of this class does not automatically delete
//! the texture from the GPU memory. You must explicitly call the Delete() 
//! method to free the texture storage on the GPU.
template <GLenum TEXTURE_TYPE>
class GLTexture3 : public GLTexture<TEXTURE_TYPE>
{
public:
	//! Sets the texture image using the given texture format, data format, and data type.
	void SetImage( GLenum textureFormat, GLenum dataFormat, GLenum dataType, void const *data, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { GLTexture<TEXTURE_TYPE>::Bind(); glTexImage3D(TEXTURE_TYPE,level,textureFormat,width,height,depth,0,dataFormat,dataType,data); }

	//! Sets the texture image using the given texture format and data format. The data type is determined by the data pointer type.
	template <typename T> void SetImage( GLenum textureFormat, GLenum dataFormat, T const *data, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { SetImage(textureFormat,dataFormat,GL::GetGLType(data),data,width,height,depth,level); }

	//! Sets the texture image using the given texture type. The data format and type are determined by the data pointer type and the number of channels.
	//! The texture format is determined by the texture type and the number of channels.
	template <typename T> void SetImage( GL::Type textureType, T const *data, int numChannels, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { SetImage(GL::TextureFormat(textureType,numChannels),GL::TextureDataFormat(GL::GetType(data),numChannels),data,width,height,depth,level); }

	//! Sets the texture image. The texture format uses the matching data pointer type.
	//! If unsigned char is used, the texture uses 8-bit normalized values.
	//! If unsigned short is used, the texture uses 16-bit normalized values.
	//! If float is used, the texture uses non-normalized 32-bit float values.
	//! If char, short, or int is used, the texture uses non-normalized 8-bit, 16-bit, or 32-bit integer values.
	template <typename T> void SetImage( T const *data, int numChannels, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { SetImage(GL::GetType(data),data,numChannels,width,height,depth,level); }

	template <typename T> void SetImageRGBA( GL::Type textureType, T const *data, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { SetImage(textureType,data,4,width,height,depth,level); }	//!< Sets the texture image with 4 channels.
	template <typename T> void SetImageRGB ( GL::Type textureType, T const *data, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { SetImage(textureType,data,3,width,height,depth,level); }	//!< Sets the texture image with 3 channels.
	template <typename T> void SetImageRG  ( GL::Type textureType, T const *data, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { SetImage(textureType,data,2,width,height,depth,level); }	//!< Sets the texture image with 2 channels.
	template <typename T> void SetImageR   ( GL::Type textureType, T const *data, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { SetImage(textureType,data,1,width,height,depth,level); }	//!< Sets the texture image with 1 channel.

	template <typename T> void SetImageRGBA( T const *data, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { SetImage(data,4,width,height,depth,level); }	//!< Sets the texture image with 4 channels.
	template <typename T> void SetImageRGB ( T const *data, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { SetImage(data,3,width,height,depth,level); }	//!< Sets the texture image with 3 channels.
	template <typename T> void SetImageRG  ( T const *data, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { SetImage(data,2,width,height,depth,level); }	//!< Sets the texture image with 2 channels.
	template <typename T> void SetImageR   ( T const *data, GLsizei width, GLsizei height, GLsizei depth, int level=0 ) { SetImage(data,1,width,height,depth,level); }	//!< Sets the texture image with 1 channel.

	//! Sets the texture wrapping parameter.
	//! The acceptable values are GL_REPEAT, GL_MIRRORED_REPEAT, GL_CLAMP, and GL_CLAMP_TO_BORDER.
	//! If the wrap argument is zero, the corresponding wrapping parameter is not changed.
	void SetWrappingMode(GLenum wrapS, GLenum wrapT, GLenum wrapR)
	{
		GLTexture<TEXTURE_TYPE>::Bind();
		if ( wrapS != 0 ) glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_WRAP_S, wrapS);
		if ( wrapT != 0 ) glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_WRAP_T, wrapT);
		if ( wrapR != 0 ) glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_WRAP_R, wrapR);
	}

};

//-------------------------------------------------------------------------------

//! Sides of a cube map
enum GLTextureCubeMapSide {
	POSITIVE_X=0,
	NEGATIVE_X,
	POSITIVE_Y,
	NEGATIVE_Y,
	POSITIVE_Z,
	NEGATIVE_Z
};

//-------------------------------------------------------------------------------

//! OpenGL cube map texture class.
//!
//! This class provides a convenient interface for handling cube map texture
//! operations with OpenGL. This class merely stores the texture id.
//! Note that deleting an object of this class does not automatically delete
//! the texture from the GPU memory. You must explicitly call the Delete() 
//! method to free the texture storage on the GPU.
class GLTextureCubeMap : public GLTexture<GL_TEXTURE_CUBE_MAP>
{
public:
	typedef GLTextureCubeMapSide Side;	//!< Sides of the cube map

	//! Sets the texture image using the given texture format, data format, and data type.
	void SetImage( Side side, GLenum textureFormat, GLenum dataFormat, GLenum dataType, void const *data, GLsizei width, GLsizei height, int level=0 ) { GLTexture<GL_TEXTURE_CUBE_MAP>::Bind(); glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X+side,level,textureFormat,width,height,0,dataFormat,dataType,data); }

	//! Sets the texture image using the given texture format and data format. The data type is determined by the data pointer type.
	template <typename T> void SetImage( Side side, GLenum textureFormat, GLenum dataFormat, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(side,textureFormat,dataFormat,GL::GetGLType(data),data,width,height,level); }

	//! Sets the texture image using the given texture type. The data format and type are determined by the data pointer type and the number of channels.
	//! The texture format is determined by the texture type and the number of channels.
	template <typename T> void SetImage( Side side, GL::Type textureType, T const *data, int numChannels, GLsizei width, GLsizei height, int level=0 ) { SetImage(side,GL::TextureFormat(textureType,numChannels),GL::TextureDataFormat(GL::GetType(data),numChannels),data,width,height,level); }

	//! Sets the texture image. The texture format uses the matching data pointer type.
	//! If unsigned char is used, the texture uses 8-bit normalized values.
	//! If unsigned short is used, the texture uses 16-bit normalized values.
	//! If float is used, the texture uses non-normalized 32-bit float values.
	//! If char, short, or int is used, the texture uses non-normalized 8-bit, 16-bit, or 32-bit integer values.
	template <typename T> void SetImage( Side side, T const *data, int numChannels, GLsizei width, GLsizei height, int level=0 ) { SetImage(side,GL::GetType(data),data,numChannels,width,height,level); }

	template <typename T> void SetImageRGBA( Side side, GL::Type textureType, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(side,textureType,data,4,width,height,level); }	//!< Sets the texture image with 4 channels.
	template <typename T> void SetImageRGB ( Side side, GL::Type textureType, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(side,textureType,data,3,width,height,level); }	//!< Sets the texture image with 3 channels.
	template <typename T> void SetImageRG  ( Side side, GL::Type textureType, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(side,textureType,data,2,width,height,level); }	//!< Sets the texture image with 2 channels.
	template <typename T> void SetImageR   ( Side side, GL::Type textureType, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(side,textureType,data,1,width,height,level); }	//!< Sets the texture image with 1 channel.

	template <typename T> void SetImageRGBA( Side side, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(side,data,4,width,height,level); }	//!< Sets the texture image with 4 channels.
	template <typename T> void SetImageRGB ( Side side, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(side,data,3,width,height,level); }	//!< Sets the texture image with 3 channels.
	template <typename T> void SetImageRG  ( Side side, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(side,data,2,width,height,level); }	//!< Sets the texture image with 2 channels.
	template <typename T> void SetImageR   ( Side side, T const *data, GLsizei width, GLsizei height, int level=0 ) { SetImage(side,data,1,width,height,level); }	//!< Sets the texture image with 1 channel.

#ifdef GL_TEXTURE_CUBE_MAP_SEAMLESS
	//! Sets the global seamless cube mapping flag, if supported by the hardware.
	static void SetSeamless(bool enable=true) { if (enable) glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS); else glDisable(GL_TEXTURE_CUBE_MAP_SEAMLESS); }
#endif

};

//-------------------------------------------------------------------------------

#ifdef GL_VERSION_3_0
#define _CY_GLRenderBuffer

//-------------------------------------------------------------------------------

//! OpenGL render buffer
//!
//! This is the base class for helper classes for render to texture in OpenGL.
template <GLenum TEXTURE_TYPE>
class GLRenderBuffer
{
protected:
	GLuint        framebufferID;		//!< The frame-buffer ID
	GLuint        depthbufferID;		//!< The depth-buffer ID
	GLTexture2<TEXTURE_TYPE> texture;	//!< The buffer texture
	GLsizei       bufferWidth;			//!< The width of the frame buffer
	GLsizei       bufferHeight;			//!< The height of the frame buffer
	mutable GLint prevBufferID;			//!< Temporary storage for previous frame-buffer used before binding this buffer
	mutable GLint prevViewport[4];		//!< Temporary storage for the size and position of the previous frame-buffer used before binding this buffer

public:
	GLRenderBuffer() : framebufferID(CY_GL_INVALID_ID), depthbufferID(CY_GL_INVALID_ID), prevBufferID(0) {}	//!< Constructor.
	~GLRenderBuffer() { if ( GL::CheckContext() ) Delete(); }										//!< Destructor.

	//!@name General Methods

	void   Delete    ();														//!< Deletes the render buffer.
	GLuint GetID     () const { return framebufferID; }							//!< Returns the frame buffer ID.
	bool   IsNull    () const { return framebufferID == CY_GL_INVALID_ID; }		//!< Returns true if the render buffer is not initialized, i.e. the render buffer id is invalid.
	void   Bind      () const;													//!< Binds the frame buffer for rendering and adjusts the viewport accordingly.
	void   Unbind    () const;													//!< Binds the frame buffer that was used before this frame buffer was bound and reverts the viewport.
	bool   IsReady   () const { return glIsFramebuffer(framebufferID) > 0; }	//!< Returns true if the frame buffer is ready. This method can be called after initialization.
	bool   IsComplete() const;													//!< Returns true if the render buffer is complete.

	//!@name Texture Methods

	GLuint GetTextureID() const { return texture.GetID(); }						//!< Returns the texture ID.
	void   BindTexture () const { texture.Bind(); }								//!< Binds the texture to the current texture unit.
	void   BindTexture (int textureUnit) const { texture.Bind(textureUnit); }	//!< Binds the texture to the given texture unit.
	void   BuildTextureMipmaps() { texture.BuildMipmaps(); }					//!< Builds mipmap levels for the texture.

	//! Sets the wrapping parameter for the texture.
	//! The acceptable values are GL_REPEAT, GL_MIRRORED_REPEAT, GL_CLAMP, and GL_CLAMP_TO_BORDER.
	//! If the wrap argument is zero, the corresponding wrapping parameter is not changed.
	void SetTextureWrappingMode(GLenum wrapS, GLenum wrapT) { texture.SetWrappingMode(wrapS,wrapT); }

	//! Sets the filtering mode for the texture.
	//! The acceptable values are GL_NEAREST and GL_LINEAR.
	//! The minification filter values can also be GL_NEAREST_MIPMAP_NEAREST, GL_LINEAR_MIPMAP_NEAREST, GL_NEAREST_MIPMAP_LINEAR, or GL_LINEAR_MIPMAP_LINEAR.
	//! If the filter argument is zero, the corresponding filter parameter is not changed.
	void SetTextureFilteringMode(GLenum magnificationFilter=0, GLenum minificationFilter=0) { texture.SetFilteringMode(magnificationFilter,minificationFilter); }

#ifdef GL_EXT_texture_filter_anisotropic
	void SetTextureAnisotropy(float anisotropy) { texture.SetAnisotropy(anisotropy); }	//!< Sets the anisotropy level of the texture.
	void SetTextureMaxAnisotropy() { texture.SetMaxAnisotropy(); }						//!< Sets anisotropic filtering to the maximum permissible value.
	void SetTextureNoAnisotropy () { texture.SetNoAnisotropy(); }						//!< Turns off anisotropic filtering.
#endif

protected:
	void GenerateBuffer();	//!< Generates the frame buffer and initializes the texture
	void SetSize(GLsizei width, GLsizei height) { bufferWidth = width; bufferHeight = height; }	//!< Sets the size of the frame buffer
};

//-------------------------------------------------------------------------------

//! OpenGL render color buffer
//!
//! This class provides a convenient interface for texture rendering in OpenGL with a color texture buffer.
template <GLenum TEXTURE_TYPE>
class GLRenderTexture : public GLRenderBuffer<TEXTURE_TYPE>
{
public:
	//!@name Render Buffer Creation and Initialization

	//! Generates the render buffer.
	//! Returns true if the render buffer is ready.
	bool Initialize( bool useDepthBuffer );

	//! Generates the render buffer and sets its size.
	//! Returns true if the render buffer is ready and complete.
	bool Initialize( bool useDepthBuffer, int numChannels, GLsizei width, GLsizei height, GL::Type type=GL::TYPE_UBYTE ) { return Initialize(useDepthBuffer) ? Resize(numChannels,width,height,type) : false; }

	//! Initializes or changes the size of the render buffer.
	//! Returns true if the buffer is complete.
	bool Resize( int numChannels, GLsizei width, GLsizei height, GL::Type type=GL::TYPE_UBYTE );
};

//-------------------------------------------------------------------------------

//! OpenGL render depth buffer
//!
//! This class provides a convenient interface for texture rendering in OpenGL with a depth texture buffer.
template <GLenum TEXTURE_TYPE>
class GLRenderDepth : public GLRenderBuffer<TEXTURE_TYPE>
{
public:
	//!@name Render Buffer Creation and Initialization

	//! Generates the render buffer.
	//! Returns true if the render buffer is ready.
	//! If depthComparisonTexture is true, initializes the texture for depth comparison.
	bool Initialize( bool depthComparisonTexture=true );

	//! Generates the render buffer and sets its size.
	//! Returns true if the render buffer is ready and complete.
	//! If depthComparisonTexture is true, initializes the texture for depth comparison.
	bool Initialize( bool depthComparisonTexture, GLsizei width, GLsizei height, GLenum depthFormat=GL_DEPTH_COMPONENT ) { return Initialize(depthComparisonTexture) ? Resize(width,height,depthFormat) : false; }

	//! Initializes or changes the size of the render buffer.
	//! Returns true if the buffer is complete.
	bool Resize( GLsizei width, GLsizei height, GLenum depthFormat=GL_DEPTH_COMPONENT );
};

//-------------------------------------------------------------------------------

//! OpenGL render buffer with a cube map texture
//!
//! This class provides a convenient interface for texture rendering in OpenGL with a cube map texture buffer.
template <GLenum ATTACHMENT_TYPE, typename BASE>
class GLRenderTextureCubeBase : public BASE
{
public:
	//! Set the render target to the given side.
	//! Should be called after binding the render texture.
	void SetTarget( GLTextureCubeMapSide side )
	{
		glFramebufferTexture2D( GL_FRAMEBUFFER, ATTACHMENT_TYPE, GL_TEXTURE_CUBE_MAP_POSITIVE_X + side, GLRenderBuffer<GL_TEXTURE_CUBE_MAP>::GetTextureID(), 0 );
	}

#ifdef _CY_MATRIX_H_INCLUDED_
	//! Returns the rotation matrix for the given side.
	static Matrix3f GetRotation( GLTextureCubeMapSide side )
	{
		static Matrix3f m[] = {
			Matrix3f(  0, 0,-1,   0,-1, 0,  -1, 0, 0 ),
			Matrix3f(  0, 0, 1,   0,-1, 0,   1, 0, 0 ),
			Matrix3f(  1, 0, 0,   0, 0, 1,   0,-1, 0 ),
			Matrix3f(  1, 0, 0,   0, 0,-1,   0, 1, 0 ),
			Matrix3f(  1, 0, 0,   0,-1, 0,   0, 0,-1 ),
			Matrix3f( -1, 0, 0,   0,-1, 0,   0, 0, 1 )
		};
		return m[side];
	}

	//! Returns the perspective projection matrix to be used by all sides.
	static Matrix4f GetProjection( float znear, float zfar ) { return Matrix4f::Perspective( Pi<float>()/2, 1, znear, zfar ); }
#endif
};

//-------------------------------------------------------------------------------

#endif // GL_VERSION_3_0

//-------------------------------------------------------------------------------

//! GLSL shader class.
//!
//! This class provides basic functionality for compiling GLSL shaders
//! either from a given source string or a given file.
//! It only stores the shader ID and it can be safely deleted after
//! the shader is used for building (linking) a GLSL program.

class GLSLShader
{
private:
	GLuint shaderID;	//!< The shader ID

public:
	GLSLShader() : shaderID(CY_GL_INVALID_ID) {}					//!< Constructor.
	virtual ~GLSLShader() { if ( GL::CheckContext() ) Delete(); }	//!< Destructor that deletes the shader.

	//!@name General Methods

	void   Delete() { if (shaderID!=CY_GL_INVALID_ID) { glDeleteShader(shaderID); shaderID=CY_GL_INVALID_ID; } }	//!< Deletes the shader.
	GLuint GetID () const { return shaderID; }						//!< Returns the shader ID.
	bool   IsNull() const { return shaderID == CY_GL_INVALID_ID; }	//!< Returns true if the OpenGL shader object is not generated, i.e. the shader id is invalid.

	//!@name Compilation Methods

	//! Compiles the shader using the given file.
	//! If the shader was previously compiled, it is deleted.
	bool CompileFile( char const *filename, GLenum shaderType, std::ostream *outStream=&std::cout ) { return CompileFile(filename,shaderType,0,nullptr,outStream); }

	//! Compiles the shader using the given file.
	//! If the shader was previously compiled, it is deleted.
	//! The prependSource string is added to the beginning of the shader code, so it must begin with the "#version" statement.
	bool CompileFile( char const *filename, GLenum shaderType, char const *prependSource, std::ostream *outStream=&std::cout ) { return CompileFile(filename,shaderType,1,&prependSource,outStream); }

	//! Compiles the shader using the given file.
	//! If the shader was previously compiled, it is deleted.
	//! The prependSources strings are added to the beginning of the shader code, so the first string must begin with "#version" statement.
	bool CompileFile( char const *filename, GLenum shaderType, int prependSourceCount, char const **prependSources, std::ostream *outStream=&std::cout );

	//! Compiles the shader using the given source code.
	//! If the shader was previously compiled, it is deleted.
	bool Compile( char const *shaderSourceCode, GLenum shaderType, std::ostream *outStream=&std::cout ) { return Compile(shaderSourceCode,shaderType,0,nullptr,outStream); }

	//! Compiles the shader using the given source code.
	//! If the shader was previously compiled, it is deleted.
	//! The prependSource string is added to the beginning of the shader code, so it must begin with the "#version" statement.
	bool Compile( char const *shaderSourceCode, GLenum shaderType, char const *prependSource, std::ostream *outStream=&std::cout ) { return Compile(shaderSourceCode,shaderType,1,&prependSource,outStream); }

	//! Compiles the shader using the given source code.
	//! If the shader was previously compiled, it is deleted.
	//! The prependSources strings are added to the beginning of the shader code, so the first string must begin with "#version" statement.
	bool Compile( char const *shaderSourceCode, GLenum shaderType, int prependSourceCount, char const **prependSources, std::ostream *outStream=&std::cout );
};

//-------------------------------------------------------------------------------

//! GLSL program class.
//!
//! This class provides basic functionality for building GLSL programs
//! using vertex and fragment shaders, along with optionally geometry and tessellation shaders.
//! The shader sources can be provides as GLSLShader class objects, source strings, or file names.
//! This class also stores a vector of registered uniform parameter IDs.

class GLSLProgram
{
private:
	GLuint programID;			//!< The program ID
	std::vector<GLint> params;	//!< A list of registered uniform parameter IDs

public:
	GLSLProgram() : programID(CY_GL_INVALID_ID) {}					//!< Constructor
	virtual ~GLSLProgram() { if ( GL::CheckContext() ) Delete(); }	//!< Destructor that deletes the program

	//!@name General Methods

	void   Delete() { if (programID!=CY_GL_INVALID_ID) { glDeleteProgram(programID); programID=CY_GL_INVALID_ID; } }	//!< Deletes the program.
	GLuint GetID () const { return programID; }						//!< Returns the program ID
	bool   IsNull() const { return programID == CY_GL_INVALID_ID; }	//!< Returns true if the OpenGL program object is not generated, i.e. the program id is invalid.
	void   Bind  () const { glUseProgram(programID); }				//!< Binds the program for rendering

	//! Attaches the given shader to the program.
	//! This function must be called before calling Link.
	void CreateProgram() { Delete(); programID = glCreateProgram(); }

	//! Attaches the given shader to the program.
	//! This function must be called before calling Link.
	void AttachShader( GLSLShader const &shader ) { AttachShader(shader.GetID()); }

	//! Attaches the given shader to the program.
	//! This function must be called before calling Link.
	void AttachShader( GLuint shaderID ) { glAttachShader(programID,shaderID); }

	//! Links the program.
	//! The shaders must be attached before calling this function.
	//! Returns true if the link operation is successful.
	//! Writes any error or warning messages to the given output stream.
	bool Link( std::ostream *outStream=&std::cout );

	//!@name Build Methods

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	bool Build( GLSLShader const *vertexShader, 
                GLSLShader const *fragmentShader,
	            GLSLShader const *geometryShader=nullptr,
	            GLSLShader const *tessControlShader=nullptr,
	            GLSLShader const *tessEvaluationShader=nullptr,
	            std::ostream *outStream=&std::cout );

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	bool BuildFiles( char const *vertexShaderFile, 
                     char const *fragmentShaderFile,
	                 char const *geometryShaderFile=nullptr,
	                 char const *tessControlShaderFile=nullptr,
	                 char const *tessEvaluationShaderFile=nullptr,
	                 std::ostream *outStream=&std::cout )
	{ return BuildFiles(vertexShaderFile,fragmentShaderFile,geometryShaderFile,tessControlShaderFile,tessEvaluationShaderFile,0,nullptr,outStream); }

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	//! The prependSource string is added to the beginning of each shader code, so it must begin with the "#version" statement.
	bool BuildFiles( char const *vertexShaderFile, 
                     char const *fragmentShaderFile,
	                 char const *geometryShaderFile,
	                 char const *tessControlShaderFile,
	                 char const *tessEvaluationShaderFile,
	                 char const *prependSource,
	                 std::ostream *outStream=&std::cout )
	{ return BuildFiles(vertexShaderFile,fragmentShaderFile,geometryShaderFile,tessControlShaderFile,tessEvaluationShaderFile,1,&prependSource,outStream); }

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	//! The prependSources strings are added to the beginning of each shader code, so the first string must begin with "#version" statement.
	bool BuildFiles( char const *vertexShaderFile, 
                     char const *fragmentShaderFile,
	                 char const *geometryShaderFile,
	                 char const *tessControlShaderFile,
	                 char const *tessEvaluationShaderFile,
	                 int         prependSourceCount,
	                 char const **prependSource,
	                 std::ostream *outStream=&std::cout );

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	bool BuildSources( char const *vertexShaderSourceCode, 
                       char const *fragmentShaderSourceCode,
	                   char const *geometryShaderSourceCode=nullptr,
	                   char const *tessControlShaderSourceCode=nullptr,
	                   char const *tessEvaluationShaderSourceCode=nullptr,
	                   std::ostream *outStream=&std::cout )
	{ return BuildSources(vertexShaderSourceCode,fragmentShaderSourceCode,geometryShaderSourceCode,tessControlShaderSourceCode,tessEvaluationShaderSourceCode,0,nullptr,outStream); }

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	//! The prependSource string is added to the beginning of each shader code, so it must begin with the "#version" statement.
	bool BuildSources( char const *vertexShaderSourceCode, 
                       char const *fragmentShaderSourceCode,
	                   char const *geometryShaderSourceCode,
	                   char const *tessControlShaderSourceCode,
	                   char const *tessEvaluationShaderSourceCode,
	                   char const *prependSource,
	                   std::ostream *outStream=&std::cout )
	{ return BuildSources(vertexShaderSourceCode,fragmentShaderSourceCode,geometryShaderSourceCode,tessControlShaderSourceCode,tessEvaluationShaderSourceCode,1,&prependSource,outStream); }

	//! Creates a program, compiles the given shaders, and links them.
	//! Returns true if all compilation and link operations are successful.
	//! Writes any error or warning messages to the given output stream.
	//! The prependSources strings are added to the beginning of each shader code, so the first string must begin with "#version" statement.
	bool BuildSources( char const *vertexShaderSourceCode, 
                       char const *fragmentShaderSourceCode,
	                   char const *geometryShaderSourceCode,
	                   char const *tessControlShaderSourceCode,
	                   char const *tessEvaluationShaderSourceCode,
	                   int         prependSourceCount,
	                   char const **prependSource,
	                   std::ostream *outStream=&std::cout );

	//!@name Uniform Parameter Methods

	//! Registers a single uniform parameter.
	//! The index must be unique and the name should match a uniform parameter name in one of the shaders.
	//! The index values for different parameters don't have to be consecutive, but unused index values waste memory.
	void RegisterUniform( unsigned int index, char const *name, std::ostream *outStream=&std::cout );

	//! Registers multiple parameters.
	//! The names should be separated by a space character.
	void RegisterUniforms( char const *names, unsigned int startingIndex=0, std::ostream *outStream=&std::cout );

	//!@{
	//! Sets the value of the uniform parameter with the given index. 
	//! The uniform parameter must be registered before using RegisterUniform() or RegisterUniforms().
	//! The program must be bind by calling Bind() before calling this method.
	void SetUniform (int index, float x)                                { glUniform1f  (params[index],x); }
	void SetUniform (int index, float x, float y)                       { glUniform2f  (params[index],x,y); }
	void SetUniform (int index, float x, float y, float z)              { glUniform3f  (params[index],x,y,z); }
	void SetUniform (int index, float x, float y, float z, float w)     { glUniform4f  (params[index],x,y,z,w); }
	void SetUniform1(int index, int count, float  const *data)          { glUniform1fv (params[index],count,data); }
	void SetUniform2(int index, int count, float  const *data)          { glUniform2fv (params[index],count,data); }
	void SetUniform3(int index, int count, float  const *data)          { glUniform3fv (params[index],count,data); }
	void SetUniform4(int index, int count, float  const *data)          { glUniform4fv (params[index],count,data); }
	void SetUniform (int index, int x)                                  { glUniform1i  (params[index],x); }
	void SetUniform (int index, int x, int y)                           { glUniform2i  (params[index],x,y); }
	void SetUniform (int index, int x, int y, int z)                    { glUniform3i  (params[index],x,y,z); }
	void SetUniform (int index, int x, int y, int z, int w)             { glUniform4i  (params[index],x,y,z,w); }
	void SetUniform1(int index, int count, int    const *data)          { glUniform1iv (params[index],count,data); }
	void SetUniform2(int index, int count, int    const *data)          { glUniform2iv (params[index],count,data); }
	void SetUniform3(int index, int count, int    const *data)          { glUniform3iv (params[index],count,data); }
	void SetUniform4(int index, int count, int    const *data)          { glUniform4iv (params[index],count,data); }
#ifdef GL_VERSION_3_0
	void SetUniform (int index, GLuint x)                               { glUniform1ui (params[index],x); }
	void SetUniform (int index, GLuint x, GLuint y)                     { glUniform2ui (params[index],x,y); }
	void SetUniform (int index, GLuint x, GLuint y, GLuint z)           { glUniform3ui (params[index],x,y,z); }
	void SetUniform (int index, GLuint x, GLuint y, GLuint z, GLuint w) { glUniform4ui (params[index],x,y,z,w); }
	void SetUniform1(int index, int count, GLuint const *data)          { glUniform1uiv(params[index],count,data); }
	void SetUniform2(int index, int count, GLuint const *data)          { glUniform2uiv(params[index],count,data); }
	void SetUniform3(int index, int count, GLuint const *data)          { glUniform3uiv(params[index],count,data); }
	void SetUniform4(int index, int count, GLuint const *data)          { glUniform4uiv(params[index],count,data); }
#endif
#ifdef GL_VERSION_4_0
	void SetUniform (int index, double x)                               { glUniform1d  (params[index],x); }
	void SetUniform (int index, double x, double y)                     { glUniform2d  (params[index],x,y); }
	void SetUniform (int index, double x, double y, double z)           { glUniform3d  (params[index],x,y,z); }
	void SetUniform (int index, double x, double y, double z, double w) { glUniform4d  (params[index],x,y,z,w); }
	void SetUniform1(int index, int count, double const *data)          { glUniform1dv (params[index],count,data); }
	void SetUniform2(int index, int count, double const *data)          { glUniform2dv (params[index],count,data); }
	void SetUniform3(int index, int count, double const *data)          { glUniform3dv (params[index],count,data); }
	void SetUniform4(int index, int count, double const *data)          { glUniform4dv (params[index],count,data); }
#endif

	void SetUniformMatrix2  (int index, float  const *m, int count=1, bool transpose=false) { glUniformMatrix2fv  (params[index],count,transpose,m); }
	void SetUniformMatrix3  (int index, float  const *m, int count=1, bool transpose=false) { glUniformMatrix3fv  (params[index],count,transpose,m); }
	void SetUniformMatrix4  (int index, float  const *m, int count=1, bool transpose=false) { glUniformMatrix4fv  (params[index],count,transpose,m); }
#ifdef GL_VERSION_2_1
	void SetUniformMatrix2x3(int index, float  const *m, int count=1, bool transpose=false) { glUniformMatrix2x3fv(params[index],count,transpose,m); }
	void SetUniformMatrix2x4(int index, float  const *m, int count=1, bool transpose=false) { glUniformMatrix2x4fv(params[index],count,transpose,m); }
	void SetUniformMatrix3x2(int index, float  const *m, int count=1, bool transpose=false) { glUniformMatrix3x2fv(params[index],count,transpose,m); }
	void SetUniformMatrix3x4(int index, float  const *m, int count=1, bool transpose=false) { glUniformMatrix3x4fv(params[index],count,transpose,m); }
	void SetUniformMatrix4x2(int index, float  const *m, int count=1, bool transpose=false) { glUniformMatrix4x2fv(params[index],count,transpose,m); }
	void SetUniformMatrix4x3(int index, float  const *m, int count=1, bool transpose=false) { glUniformMatrix4x3fv(params[index],count,transpose,m); }
#endif
#ifdef GL_VERSION_4_0
	void SetUniformMatrix2  (int index, double const *m, int count=1, bool transpose=false) { glUniformMatrix2dv  (params[index],count,transpose,m); }
	void SetUniformMatrix3  (int index, double const *m, int count=1, bool transpose=false) { glUniformMatrix3dv  (params[index],count,transpose,m); }
	void SetUniformMatrix4  (int index, double const *m, int count=1, bool transpose=false) { glUniformMatrix4dv  (params[index],count,transpose,m); }
	void SetUniformMatrix2x3(int index, double const *m, int count=1, bool transpose=false) { glUniformMatrix2x3dv(params[index],count,transpose,m); }
	void SetUniformMatrix2x4(int index, double const *m, int count=1, bool transpose=false) { glUniformMatrix2x4dv(params[index],count,transpose,m); }
	void SetUniformMatrix3x2(int index, double const *m, int count=1, bool transpose=false) { glUniformMatrix3x2dv(params[index],count,transpose,m); }	
	void SetUniformMatrix3x4(int index, double const *m, int count=1, bool transpose=false) { glUniformMatrix3x4dv(params[index],count,transpose,m); }	
	void SetUniformMatrix4x2(int index, double const *m, int count=1, bool transpose=false) { glUniformMatrix4x2dv(params[index],count,transpose,m); }	
	void SetUniformMatrix4x3(int index, double const *m, int count=1, bool transpose=false) { glUniformMatrix4x3dv(params[index],count,transpose,m); }	
#endif

#ifdef _CY_VECTOR_H_INCLUDED_
	void SetUniform(int index, Vec2<float>  const &p)              { glUniform2fv (params[index],1,    &p.x ); }
	void SetUniform(int index, Vec3<float>  const &p)              { glUniform3fv (params[index],1,    &p.x ); }
	void SetUniform(int index, Vec4<float>  const &p)              { glUniform4fv (params[index],1,    &p.x ); }
	void SetUniform(int index, Vec2<int>    const &p)              { glUniform2iv (params[index],1,    &p.x ); }
	void SetUniform(int index, Vec3<int>    const &p)              { glUniform3iv (params[index],1,    &p.x ); }
	void SetUniform(int index, Vec4<int>    const &p)              { glUniform4iv (params[index],1,    &p.x ); }
	void SetUniform(int index, Vec2<float>  const *p, int count=1) { glUniform2fv (params[index],count,&p->x); }
	void SetUniform(int index, Vec3<float>  const *p, int count=1) { glUniform3fv (params[index],count,&p->x); }
	void SetUniform(int index, Vec4<float>  const *p, int count=1) { glUniform4fv (params[index],count,&p->x); }
	void SetUniform(int index, Vec2<int>    const *p, int count=1) { glUniform2iv (params[index],count,&p->x); }
	void SetUniform(int index, Vec3<int>    const *p, int count=1) { glUniform3iv (params[index],count,&p->x); }
	void SetUniform(int index, Vec4<int>    const *p, int count=1) { glUniform4iv (params[index],count,&p->x); }
# ifdef GL_VERSION_3_0
	void SetUniform(int index, Vec2<GLuint> const &p)              { glUniform2uiv(params[index],1,    &p.x ); }
	void SetUniform(int index, Vec3<GLuint> const &p)              { glUniform3uiv(params[index],1,    &p.x ); }
	void SetUniform(int index, Vec4<GLuint> const &p)              { glUniform4uiv(params[index],1,    &p.x ); }
	void SetUniform(int index, Vec2<GLuint> const *p, int count=1) { glUniform2uiv(params[index],count,&p->x); }
	void SetUniform(int index, Vec3<GLuint> const *p, int count=1) { glUniform3uiv(params[index],count,&p->x); }
	void SetUniform(int index, Vec4<GLuint> const *p, int count=1) { glUniform4uiv(params[index],count,&p->x); }
# endif
# ifdef GL_VERSION_4_0
	void SetUniform(int index, Vec2<double> const &p)              { glUniform2dv (params[index],1,    &p.x ); }
	void SetUniform(int index, Vec3<double> const &p)              { glUniform3dv (params[index],1,    &p.x ); }
	void SetUniform(int index, Vec4<double> const &p)              { glUniform4dv (params[index],1,    &p.x ); }
	void SetUniform(int index, Vec2<double> const *p, int count=1) { glUniform2dv (params[index],count,&p->x); }
	void SetUniform(int index, Vec3<double> const *p, int count=1) { glUniform3dv (params[index],count,&p->x); }
	void SetUniform(int index, Vec4<double> const *p, int count=1) { glUniform4dv (params[index],count,&p->x); }
# endif
#endif

#ifdef _CY_IVECTOR_H_INCLUDED_
	void SetUniform(int index, IVec2<int>    const &p)              { glUniform2iv (params[index],1,    &p.x ); }
	void SetUniform(int index, IVec3<int>    const &p)              { glUniform3iv (params[index],1,    &p.x ); }
	void SetUniform(int index, IVec4<int>    const &p)              { glUniform4iv (params[index],1,    &p.x ); }
	void SetUniform(int index, IVec2<int>    const *p, int count=1) { glUniform2iv (params[index],count,&p->x); }
	void SetUniform(int index, IVec3<int>    const *p, int count=1) { glUniform3iv (params[index],count,&p->x); }
	void SetUniform(int index, IVec4<int>    const *p, int count=1) { glUniform4iv (params[index],count,&p->x); }
# ifdef GL_VERSION_3_0
	void SetUniform(int index, IVec2<GLuint> const &p)              { glUniform2uiv(params[index],1,    &p.x ); }
	void SetUniform(int index, IVec3<GLuint> const &p)              { glUniform3uiv(params[index],1,    &p.x ); }
	void SetUniform(int index, IVec4<GLuint> const &p)              { glUniform4uiv(params[index],1,    &p.x ); }
	void SetUniform(int index, IVec2<GLuint> const *p, int count=1) { glUniform2uiv(params[index],count,&p->x); }
	void SetUniform(int index, IVec3<GLuint> const *p, int count=1) { glUniform3uiv(params[index],count,&p->x); }
	void SetUniform(int index, IVec4<GLuint> const *p, int count=1) { glUniform4uiv(params[index],count,&p->x); }
# endif
#endif

#ifdef _CY_MATRIX_H_INCLUDED_
	void SetUniform(int index, Matrix2 <float>  const &m)              { glUniformMatrix2fv  (params[index],1,    GL_FALSE,m.cell ); }
	void SetUniform(int index, Matrix3 <float>  const &m)              { glUniformMatrix3fv  (params[index],1,    GL_FALSE,m.cell ); }
	void SetUniform(int index, Matrix4 <float>  const &m)              { glUniformMatrix4fv  (params[index],1,    GL_FALSE,m.cell ); }
	void SetUniform(int index, Matrix2 <float>  const *m, int count=1) { glUniformMatrix2fv  (params[index],count,GL_FALSE,m->cell); }
	void SetUniform(int index, Matrix3 <float>  const *m, int count=1) { glUniformMatrix3fv  (params[index],count,GL_FALSE,m->cell); }
	void SetUniform(int index, Matrix4 <float>  const *m, int count=1) { glUniformMatrix4fv  (params[index],count,GL_FALSE,m->cell); }
# ifdef GL_VERSION_2_1
	void SetUniform(int index, Matrix34<float>  const &m)              { glUniformMatrix3x4fv(params[index],1,    GL_FALSE,m.cell ); }
	void SetUniform(int index, Matrix34<float>  const *m, int count=1) { glUniformMatrix3x4fv(params[index],count,GL_FALSE,m->cell); }
# endif
# ifdef GL_VERSION_4_0
	void SetUniform(int index, Matrix2 <double> const &m)              { glUniformMatrix2dv  (params[index],1,    GL_FALSE,m.cell ); }
	void SetUniform(int index, Matrix3 <double> const &m)              { glUniformMatrix3dv  (params[index],1,    GL_FALSE,m.cell ); }
	void SetUniform(int index, Matrix4 <double> const &m)              { glUniformMatrix4dv  (params[index],1,    GL_FALSE,m.cell ); }
	void SetUniform(int index, Matrix34<double> const &m)              { glUniformMatrix3x4dv(params[index],1,    GL_FALSE,m.cell ); }
	void SetUniform(int index, Matrix2 <double> const *m, int count=1) { glUniformMatrix2dv  (params[index],count,GL_FALSE,m->cell); }
	void SetUniform(int index, Matrix3 <double> const *m, int count=1) { glUniformMatrix3dv  (params[index],count,GL_FALSE,m->cell); }
	void SetUniform(int index, Matrix4 <double> const *m, int count=1) { glUniformMatrix4dv  (params[index],count,GL_FALSE,m->cell); }
	void SetUniform(int index, Matrix34<double> const *m, int count=1) { glUniformMatrix3x4dv(params[index],count,GL_FALSE,m->cell); }
# endif
#endif
	//!@}


	//!@{ glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 )
	//! Sets the value of the uniform parameter with the given name, if the uniform parameter is found. 
	//! Since it searches for the uniform parameter first, it is not as efficient as setting the uniform parameter using
	//! a previously registered id. There is no need to bind the program before calling this method.
	void SetUniform (char const *name, float x)                                { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform1f  (id,x); }
	void SetUniform (char const *name, float x, float y)                       { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2f  (id,x,y); }
	void SetUniform (char const *name, float x, float y, float z)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3f  (id,x,y,z); }
	void SetUniform (char const *name, float x, float y, float z, float w)     { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4f  (id,x,y,z,w); }
	void SetUniform1(char const *name, int count, float  const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform1fv (id,count,data); }
	void SetUniform2(char const *name, int count, float  const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2fv (id,count,data); }
	void SetUniform3(char const *name, int count, float  const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3fv (id,count,data); }
	void SetUniform4(char const *name, int count, float  const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4fv (id,count,data); }
	void SetUniform (char const *name, int x)                                  { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform1i  (id,x); }
	void SetUniform (char const *name, int x, int y)                           { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2i  (id,x,y); }
	void SetUniform (char const *name, int x, int y, int z)                    { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3i  (id,x,y,z); }
	void SetUniform (char const *name, int x, int y, int z, int w)             { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4i  (id,x,y,z,w); }
	void SetUniform1(char const *name, int count, int    const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform1iv (id,count,data); }
	void SetUniform2(char const *name, int count, int    const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2iv (id,count,data); }
	void SetUniform3(char const *name, int count, int    const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3iv (id,count,data); }
	void SetUniform4(char const *name, int count, int    const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4iv (id,count,data); }
#ifdef GL_VERSION_3_0
	void SetUniform (char const *name, GLuint x)                               { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform1ui (id,x); }
	void SetUniform (char const *name, GLuint x, GLuint y)                     { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2ui (id,x,y); }
	void SetUniform (char const *name, GLuint x, GLuint y, GLuint z)           { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3ui (id,x,y,z); }
	void SetUniform (char const *name, GLuint x, GLuint y, GLuint z, GLuint w) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4ui (id,x,y,z,w); }
	void SetUniform1(char const *name, int count, GLuint const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform1uiv(id,count,data); }
	void SetUniform2(char const *name, int count, GLuint const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2uiv(id,count,data); }
	void SetUniform3(char const *name, int count, GLuint const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3uiv(id,count,data); }
	void SetUniform4(char const *name, int count, GLuint const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4uiv(id,count,data); }
#endif
#ifdef GL_VERSION_4_0
	void SetUniform (char const *name, double x)                               { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform1d  (id,x); }
	void SetUniform (char const *name, double x, double y)                     { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2d  (id,x,y); }
	void SetUniform (char const *name, double x, double y, double z)           { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3d  (id,x,y,z); }
	void SetUniform (char const *name, double x, double y, double z, double w) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4d  (id,x,y,z,w); }
	void SetUniform1(char const *name, int count, double const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform1dv (id,count,data); }
	void SetUniform2(char const *name, int count, double const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2dv (id,count,data); }
	void SetUniform3(char const *name, int count, double const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3dv (id,count,data); }
	void SetUniform4(char const *name, int count, double const *data)          { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4dv (id,count,data); }
#endif

	void SetUniformMatrix2  (char const *name, float  const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix2fv  (id,count,transpose,m); }
	void SetUniformMatrix3  (char const *name, float  const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3fv  (id,count,transpose,m); }
	void SetUniformMatrix4  (char const *name, float  const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix4fv  (id,count,transpose,m); }
#ifdef GL_VERSION_2_1
	void SetUniformMatrix2x3(char const *name, float  const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix2x3fv(id,count,transpose,m); }
	void SetUniformMatrix2x4(char const *name, float  const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix2x4fv(id,count,transpose,m); }
	void SetUniformMatrix3x2(char const *name, float  const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3x2fv(id,count,transpose,m); }
	void SetUniformMatrix3x4(char const *name, float  const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3x4fv(id,count,transpose,m); }
	void SetUniformMatrix4x2(char const *name, float  const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix4x2fv(id,count,transpose,m); }
	void SetUniformMatrix4x3(char const *name, float  const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix4x3fv(id,count,transpose,m); }
#endif
#ifdef GL_VERSION_4_0
	void SetUniformMatrix2  (char const *name, double const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix2dv  (id,count,transpose,m); }
	void SetUniformMatrix3  (char const *name, double const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3dv  (id,count,transpose,m); }
	void SetUniformMatrix4  (char const *name, double const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix4dv  (id,count,transpose,m); }
	void SetUniformMatrix2x3(char const *name, double const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix2x3dv(id,count,transpose,m); }
	void SetUniformMatrix2x4(char const *name, double const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix2x4dv(id,count,transpose,m); }
	void SetUniformMatrix3x2(char const *name, double const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3x2dv(id,count,transpose,m); }	
	void SetUniformMatrix3x4(char const *name, double const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3x4dv(id,count,transpose,m); }	
	void SetUniformMatrix4x2(char const *name, double const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix4x2dv(id,count,transpose,m); }	
	void SetUniformMatrix4x3(char const *name, double const *m, int count=1, bool transpose=false) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix4x3dv(id,count,transpose,m); }	
#endif

#ifdef _CY_VECTOR_H_INCLUDED_
	void SetUniform(char const *name, Vec2<float>  const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2fv (id,1,    &p.x ); }
	void SetUniform(char const *name, Vec3<float>  const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3fv (id,1,    &p.x ); }
	void SetUniform(char const *name, Vec4<float>  const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4fv (id,1,    &p.x ); }
	void SetUniform(char const *name, Vec2<int>    const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2iv (id,1,    &p.x ); }
	void SetUniform(char const *name, Vec3<int>    const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3iv (id,1,    &p.x ); }
	void SetUniform(char const *name, Vec4<int>    const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4iv (id,1,    &p.x ); }
	void SetUniform(char const *name, Vec2<float>  const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2fv (id,count,&p->x); }
	void SetUniform(char const *name, Vec3<float>  const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3fv (id,count,&p->x); }
	void SetUniform(char const *name, Vec4<float>  const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4fv (id,count,&p->x); }
	void SetUniform(char const *name, Vec2<int>    const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2iv (id,count,&p->x); }
	void SetUniform(char const *name, Vec3<int>    const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3iv (id,count,&p->x); }
	void SetUniform(char const *name, Vec4<int>    const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4iv (id,count,&p->x); }
# ifdef GL_VERSION_3_0
	void SetUniform(char const *name, Vec2<GLuint> const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2uiv(id,1,    &p.x ); }
	void SetUniform(char const *name, Vec3<GLuint> const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3uiv(id,1,    &p.x ); }
	void SetUniform(char const *name, Vec4<GLuint> const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4uiv(id,1,    &p.x ); }
	void SetUniform(char const *name, Vec2<GLuint> const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2uiv(id,count,&p->x); }
	void SetUniform(char const *name, Vec3<GLuint> const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3uiv(id,count,&p->x); }
	void SetUniform(char const *name, Vec4<GLuint> const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4uiv(id,count,&p->x); }
# endif
# ifdef GL_VERSION_4_0
	void SetUniform(char const *name, Vec2<double> const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2dv (id,1,    &p.x ); }
	void SetUniform(char const *name, Vec3<double> const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3dv (id,1,    &p.x ); }
	void SetUniform(char const *name, Vec4<double> const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4dv (id,1,    &p.x ); }
	void SetUniform(char const *name, Vec2<double> const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2dv (id,count,&p->x); }
	void SetUniform(char const *name, Vec3<double> const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3dv (id,count,&p->x); }
	void SetUniform(char const *name, Vec4<double> const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4dv (id,count,&p->x); }
# endif
#endif

#ifdef _CY_IVECTOR_H_INCLUDED_
	void SetUniform(char const *name, IVec2<int>    const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2iv (id,1,    &p.x ); }
	void SetUniform(char const *name, IVec3<int>    const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3iv (id,1,    &p.x ); }
	void SetUniform(char const *name, IVec4<int>    const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4iv (id,1,    &p.x ); }
	void SetUniform(char const *name, IVec2<int>    const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2iv (id,count,&p->x); }
	void SetUniform(char const *name, IVec3<int>    const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3iv (id,count,&p->x); }
	void SetUniform(char const *name, IVec4<int>    const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4iv (id,count,&p->x); }
# ifdef GL_VERSION_3_0
	void SetUniform(char const *name, IVec2<GLuint> const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2uiv(id,1,    &p.x ); }
	void SetUniform(char const *name, IVec3<GLuint> const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3uiv(id,1,    &p.x ); }
	void SetUniform(char const *name, IVec4<GLuint> const &p)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4uiv(id,1,    &p.x ); }
	void SetUniform(char const *name, IVec2<GLuint> const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform2uiv(id,count,&p->x); }
	void SetUniform(char const *name, IVec3<GLuint> const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform3uiv(id,count,&p->x); }
	void SetUniform(char const *name, IVec4<GLuint> const *p, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniform4uiv(id,count,&p->x); }
# endif
#endif

#ifdef _CY_MATRIX_H_INCLUDED_
	void SetUniform(char const *name, Matrix2 <float>  const &m)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix2fv  (id,1,    GL_FALSE,m.cell ); }
	void SetUniform(char const *name, Matrix3 <float>  const &m)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3fv  (id,1,    GL_FALSE,m.cell ); }
	void SetUniform(char const *name, Matrix4 <float>  const &m)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix4fv  (id,1,    GL_FALSE,m.cell ); }
	void SetUniform(char const *name, Matrix2 <float>  const *m, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix2fv  (id,count,GL_FALSE,m->cell); }
	void SetUniform(char const *name, Matrix3 <float>  const *m, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3fv  (id,count,GL_FALSE,m->cell); }
	void SetUniform(char const *name, Matrix4 <float>  const *m, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix4fv  (id,count,GL_FALSE,m->cell); }
# ifdef GL_VERSION_2_1
	void SetUniform(char const *name, Matrix34<float>  const &m)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3x4fv(id,1,    GL_FALSE,m.cell ); }
	void SetUniform(char const *name, Matrix34<float>  const *m, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3x4fv(id,count,GL_FALSE,m->cell); }
# endif
# ifdef GL_VERSION_4_0
	void SetUniform(char const *name, Matrix2 <double> const &m)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix2dv  (id,1,    GL_FALSE,m.cell ); }
	void SetUniform(char const *name, Matrix3 <double> const &m)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3dv  (id,1,    GL_FALSE,m.cell ); }
	void SetUniform(char const *name, Matrix4 <double> const &m)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix4dv  (id,1,    GL_FALSE,m.cell ); }
	void SetUniform(char const *name, Matrix34<double> const &m)              { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3x4dv(id,1,    GL_FALSE,m.cell ); }
	void SetUniform(char const *name, Matrix2 <double> const *m, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix2dv  (id,count,GL_FALSE,m->cell); }
	void SetUniform(char const *name, Matrix3 <double> const *m, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3dv  (id,count,GL_FALSE,m->cell); }
	void SetUniform(char const *name, Matrix4 <double> const *m, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix4dv  (id,count,GL_FALSE,m->cell); }
	void SetUniform(char const *name, Matrix34<double> const *m, int count=1) { glUseProgram(programID); int id = glGetUniformLocation( programID, name ); if ( id >= 0 ) glUniformMatrix3x4dv(id,count,GL_FALSE,m->cell); }
# endif
#endif
	//!@}

	class Param
	{
		friend GLSLProgram;
		GLSLProgram &prog;
		char const *name;
		Param(GLSLProgram &p, char const *n) : prog(p), name(n) {}
	public:
		template <typename T> void operator = ( T const &v ) { prog.SetUniform( name, v ); }
		template <typename T> void Set( T const *v, int count=1 ) { prog.SetUniform( name, v, count ); }
	};

	Param operator [] ( char const *name ) { return Param(*this,name); }


	GLint AttribLocation( char const *name ) const { return glGetAttribLocation( programID, name ); }
	void SetAttribBuffer( char const *name, GLint arrayBufferID, int dimensions, GLenum type=GL_FLOAT, GLboolean normalized=GL_FALSE, GLsizei stride=0, size_t offset=0 )
	{
		glBindBuffer( GL_ARRAY_BUFFER, arrayBufferID );
		GLint a = AttribLocation(name);
		glVertexAttribPointer( a, dimensions, type, normalized, stride, (const void*)offset );
		glEnableVertexAttribArray(a);
	}
	void EnableAttrib ( char const *name ) { glEnableVertexAttribArray ( AttribLocation(name) ); }
	void DisableAttrib( char const *name ) { glDisableVertexAttribArray( AttribLocation(name) ); }
};

//-------------------------------------------------------------------------------
// Implementation of GL
//-------------------------------------------------------------------------------

inline void GL::PrintVersion(std::ostream *outStream)
{
	const GLubyte* version = glGetString(GL_VERSION);
	if ( version ) *outStream << version;
	else {
		int versionMajor=0, versionMinor=0;
#ifdef GL_VERSION_3_0
		glGetIntegerv(GL_MAJOR_VERSION,&versionMajor);
		glGetIntegerv(GL_MINOR_VERSION,&versionMinor);
#endif
		if ( versionMajor > 0 && versionMajor < 100 ) {
			*outStream << versionMajor << "." << versionMinor;
		} else *outStream << "Unknown";
	}
}

inline void GL::CheckError( char const *sourcefile, int line, char const *call, std::ostream *outStream )
{
	GLenum error;
	while ( (error = glGetError()) != GL_NO_ERROR) {
		*outStream << "OpenGL ERROR: " << sourcefile << " (line " << line << "): ";
		if ( call ) *outStream << call << " triggered ";
		*outStream << gluErrorString(error) << std::endl;
	}
}

inline GLenum GL::TextureFormat( GL::Type type, int numChannels )
{
	assert( numChannels > 0 && numChannels <= 4 );
	GLenum const internalFormats[][4] = {
#ifdef CY_GL_TEXTURE_LUMINANCE
		{ GL_LUMINANCE8,   GL_LUMINANCE8_ALPHA8,    GL_RGB8,    GL_RGBA8	},
		{ GL_LUMINANCE16,  GL_LUMINANCE16_ALPHA16,  GL_RGB16,   GL_RGBA16   },
# ifdef GL_VERSION_3_0
		{ GL_LUMINANCE16F_ARB,  GL_LUMINANCE_ALPHA16F_ARB,  GL_RGB16F,  GL_RGBA16F  },
		{ GL_LUMINANCE32F_ARB,  GL_LUMINANCE_ALPHA32F_ARB,  GL_RGB32F,  GL_RGBA32F  },
		{ GL_LUMINANCE8I_EXT,   GL_LUMINANCE_ALPHA8I_EXT,   GL_RGB8I,   GL_RGBA8I   },
		{ GL_LUMINANCE8UI_EXT,  GL_LUMINANCE_ALPHA8UI_EXT,  GL_RGB8UI,  GL_RGBA8UI  },
		{ GL_LUMINANCE16I_EXT,  GL_LUMINANCE_ALPHA16I_EXT,  GL_RGB16I,  GL_RGBA16I  },
		{ GL_LUMINANCE16UI_EXT, GL_LUMINANCE_ALPHA16UI_EXT, GL_RGB16UI, GL_RGBA16UI },
		{ GL_LUMINANCE32I_EXT,  GL_LUMINANCE_ALPHA32I_EXT,  GL_RGB32I,  GL_RGBA32I  },
		{ GL_LUMINANCE32UI_EXT, GL_LUMINANCE_ALPHA32UI_EXT, GL_RGB32UI, GL_RGBA32UI },
# else
		{ GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA },
		{ GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA },
		{ GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA },
		{ GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA },
		{ GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA },
		{ GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA },
		{ GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA },
		{ GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA },
# endif
#else
		{ GL_R8,    GL_RG8,    GL_RGB8,    GL_RGBA8	   },
		{ GL_R16,   GL_RG16,   GL_RGB16,   GL_RGBA16   },
# ifdef GL_VERSION_3_0
		{ GL_R16F,  GL_RG16F,  GL_RGB16F,  GL_RGBA16F  },
		{ GL_R32F,  GL_RG32F,  GL_RGB32F,  GL_RGBA32F  },
		{ GL_R8I,   GL_RG8I,   GL_RGB8I,   GL_RGBA8I   },
		{ GL_R8UI,  GL_RG8UI,  GL_RGB8UI,  GL_RGBA8UI  },
		{ GL_R16I,  GL_RG16I,  GL_RGB16I,  GL_RGBA16I  },
		{ GL_R16UI, GL_RG16UI, GL_RGB16UI, GL_RGBA16UI },
		{ GL_R32I,  GL_RG32I,  GL_RGB32I,  GL_RGBA32I  },
		{ GL_R32UI, GL_RG32UI, GL_RGB32UI, GL_RGBA32UI },
# else
		{ GL_RED, GL_RG, GL_RGB, GL_RGBA },
		{ GL_RED, GL_RG, GL_RGB, GL_RGBA },
		{ GL_RED, GL_RG, GL_RGB, GL_RGBA },
		{ GL_RED, GL_RG, GL_RGB, GL_RGBA },
		{ GL_RED, GL_RG, GL_RGB, GL_RGBA },
		{ GL_RED, GL_RG, GL_RGB, GL_RGBA },
		{ GL_RED, GL_RG, GL_RGB, GL_RGBA },
		{ GL_RED, GL_RG, GL_RGB, GL_RGBA },
# endif
#endif
	};
	
	return internalFormats[type][numChannels-1];
}

inline GLenum GL::TextureDataFormat( GL::Type type, int numChannels )
{
	GLenum const formats[][4] = {
#ifdef CY_GL_TEXTURE_LUMINANCE
		{ GL_LUMINANCE,             GL_LUMINANCE_ALPHA,             GL_RGB,         GL_RGBA         },
# ifdef GL_VERSION_3_0
		{ GL_LUMINANCE_INTEGER_EXT, GL_LUMINANCE_ALPHA_INTEGER_EXT, GL_RGB_INTEGER, GL_RGBA_INTEGER },
# else
		{ GL_LUMINANCE,             GL_LUMINANCE_ALPHA,             GL_RGB,         GL_RGBA         },
# endif
#else
		{ GL_RED,         GL_RG,         GL_RGB,         GL_RGBA         },
# ifdef GL_VERSION_3_0
		{ GL_RED_INTEGER, GL_RG_INTEGER, GL_RGB_INTEGER, GL_RGBA_INTEGER },
# else
		{ GL_RED,         GL_RG,         GL_RGB,         GL_RGBA         },
# endif
#endif
	};
	return formats[type>=GL::TYPE_INT8][numChannels-1];
}

//-------------------------------------------------------------------------------
// Implementation of GLDebugCallback
//-------------------------------------------------------------------------------
#ifdef _CY_GLDebugCallback

inline void _CY_APIENTRY GLDebugCallback::Callback( GLenum source,
                                                    GLenum type,
                                                    GLuint id,
                                                    GLenum severity,
                                                    GLsizei length,
                                                    GLchar const* message,
                                                    void const* userParam )
{
	std::ostream *outStream = (std::ostream*) userParam;

	*outStream << std::endl;
	*outStream << "OpenGL Debug Output:" << std::endl;
	*outStream << "VERSION:  ";
	GL::PrintVersion(outStream);
	*outStream << std::endl;

	*outStream << "SOURCE:   ";
	switch (source) {
		case GL_DEBUG_SOURCE_API:             *outStream << "API";             break;
		case GL_DEBUG_SOURCE_WINDOW_SYSTEM:   *outStream << "Window System";   break;
		case GL_DEBUG_SOURCE_SHADER_COMPILER: *outStream << "Shader Compiler"; break;
		case GL_DEBUG_SOURCE_THIRD_PARTY:     *outStream << "Third Party";     break;
		case GL_DEBUG_SOURCE_APPLICATION:     *outStream << "Application";     break;
		case GL_DEBUG_SOURCE_OTHER:           *outStream << "Other";           break;
		default:                              *outStream << "Unknown";         break;
	}
	*outStream << std::endl;

	*outStream << "TYPE:     ";
	switch (type) {
		case GL_DEBUG_TYPE_ERROR:               *outStream << "Error";               break;
		case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: *outStream << "Deprecated Behavior"; break;
		case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:  *outStream << "Undefined Behavior";  break;
		case GL_DEBUG_TYPE_PORTABILITY:         *outStream << "Portability";         break;
		case GL_DEBUG_TYPE_PERFORMANCE:         *outStream << "Performance";         break;
		case GL_DEBUG_TYPE_OTHER:               *outStream << "Other";               break;
		default:                                *outStream << "Unknown";             break;
	}
	*outStream << std::endl;

	*outStream << "ID:       " << id << std::endl;

	*outStream << "SEVERITY: ";
	switch (severity) {
		case GL_DEBUG_SEVERITY_HIGH:         *outStream << "High";         break;
		case GL_DEBUG_SEVERITY_MEDIUM:       *outStream << "Medium";       break;
		case GL_DEBUG_SEVERITY_LOW:          *outStream << "Low";          break;
		case GL_DEBUG_SEVERITY_NOTIFICATION: *outStream << "Notification"; break;
		default:                             *outStream << "Unknown";      break;
	}
	*outStream << std::endl;

	*outStream << "MESSAGE:  ";
	*outStream << message << std::endl;

	// You can set a breakpoint at the following line. Your debugger will stop the execution,
	// and the call stack will show the OpenGL call causing this callback.
	*outStream << std::endl;
}

#endif
//-------------------------------------------------------------------------------
// GLTexture Implementation
//-------------------------------------------------------------------------------

template <GLenum TEXTURE_TYPE>
inline void GLTexture<TEXTURE_TYPE>::Initialize()
{
	if ( textureID == CY_GL_INVALID_ID ) glGenTextures(1,&textureID);
	glBindTexture(TEXTURE_TYPE, textureID);
	glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

template <GLenum TEXTURE_TYPE>
void GLTexture<TEXTURE_TYPE>::SetFilteringMode(GLenum magnificationFilter, GLenum minificationFilter)
{
	Bind();
	if ( magnificationFilter != 0 ) glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_MAG_FILTER, magnificationFilter);
	if ( minificationFilter  != 0 ) glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_MIN_FILTER, minificationFilter);
}

//-------------------------------------------------------------------------------
// GLRenderBuffer Implementation
//-------------------------------------------------------------------------------
#ifdef _CY_GLRenderBuffer

template <GLenum TEXTURE_TYPE>
inline void GLRenderBuffer<TEXTURE_TYPE>::Delete()
{
	if ( framebufferID != CY_GL_INVALID_ID ) glDeleteFramebuffers (1,&framebufferID); framebufferID = CY_GL_INVALID_ID; 
	if ( depthbufferID != CY_GL_INVALID_ID ) glDeleteRenderbuffers(1,&depthbufferID); depthbufferID = CY_GL_INVALID_ID; 
	texture.Delete();
}

template <GLenum TEXTURE_TYPE>
inline void GLRenderBuffer<TEXTURE_TYPE>::Bind() const
{
	glGetIntegerv(GL_VIEWPORT, prevViewport);
	glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &prevBufferID);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, framebufferID);
	glViewport(0,0,bufferWidth,bufferHeight);
}

template <GLenum TEXTURE_TYPE>
inline void GLRenderBuffer<TEXTURE_TYPE>::Unbind() const
{
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, prevBufferID);
	glViewport(prevViewport[0],prevViewport[1],prevViewport[2],prevViewport[3]);
}

template <GLenum TEXTURE_TYPE>
inline bool GLRenderBuffer<TEXTURE_TYPE>::IsComplete() const
{
	GLint prevbuffer;
	glGetIntegerv(GL_FRAMEBUFFER_BINDING,&prevbuffer);
	glBindFramebuffer(GL_FRAMEBUFFER,framebufferID);
	bool complete = (glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
	glBindFramebuffer(GL_FRAMEBUFFER,prevbuffer);
	return complete;  
}

template <GLenum TEXTURE_TYPE>
inline void GLRenderBuffer<TEXTURE_TYPE>::GenerateBuffer()
{
	GLRenderBuffer<TEXTURE_TYPE>::Delete();
	glGenFramebuffers(1, &framebufferID);
	glBindFramebuffer(GL_FRAMEBUFFER, framebufferID);
	texture.Initialize();
	texture.SetFilteringMode(GL_NEAREST,GL_NEAREST);
	texture.SetWrappingMode(GL_CLAMP_TO_EDGE,GL_CLAMP_TO_EDGE);
}

template <GLenum TEXTURE_TYPE>
inline bool GLRenderTexture<TEXTURE_TYPE>::Initialize( bool useDepthBuffer )
{
	GLint prevBuffer;
	glGetIntegerv( GL_FRAMEBUFFER_BINDING, &prevBuffer );
	GLRenderBuffer<TEXTURE_TYPE>::GenerateBuffer();
	if ( useDepthBuffer ) {
		glGenRenderbuffers(1, &this->depthbufferID);
		glBindRenderbuffer(GL_RENDERBUFFER, this->depthbufferID);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, this->depthbufferID);
	}
	if ( TEXTURE_TYPE != GL_TEXTURE_CUBE_MAP ) {
		glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GLRenderBuffer<TEXTURE_TYPE>::GetTextureID(), 0);
	}
	glDrawBuffer(GL_COLOR_ATTACHMENT0);
	glBindFramebuffer(GL_FRAMEBUFFER, prevBuffer);
	return GLRenderBuffer<TEXTURE_TYPE>::IsReady();
}

template <GLenum TEXTURE_TYPE>
inline bool GLRenderTexture<TEXTURE_TYPE>::Resize( int numChannels, GLsizei width, GLsizei height, GL::Type type )
{
	GLenum textureFormat = GL::TextureFormat(type,numChannels);
	if ( TEXTURE_TYPE == GL_TEXTURE_CUBE_MAP ) {
		for ( int i=0; i<6; ++i ) {
			glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X+i,0,textureFormat,width,height,0,GL_RGBA,GL_UNSIGNED_BYTE,nullptr);
		}
	} else {
		this->texture.SetImage(textureFormat,GL_RGBA,GL_UNSIGNED_BYTE,nullptr,width,height);
	}
	if ( this->depthbufferID != CY_GL_INVALID_ID ) {
		glBindRenderbuffer(GL_RENDERBUFFER, this->depthbufferID);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, width, height);
	}
	GLRenderBuffer<TEXTURE_TYPE>::SetSize(width,height);
	return GLRenderBuffer<TEXTURE_TYPE>::IsComplete();
}

template <GLenum TEXTURE_TYPE>
inline bool GLRenderDepth<TEXTURE_TYPE>::Initialize( bool depthComparisonTexture )
{
	GLint prevBuffer;
	glGetIntegerv( GL_FRAMEBUFFER_BINDING, &prevBuffer );
	GLRenderBuffer<TEXTURE_TYPE>::GenerateBuffer();
	if ( depthComparisonTexture ) {
		glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
		glTexParameteri(TEXTURE_TYPE, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	}
	if ( TEXTURE_TYPE != GL_TEXTURE_CUBE_MAP ) {
		glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GLRenderBuffer<TEXTURE_TYPE>::GetTextureID(), 0);
	}
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, prevBuffer);
	return GLRenderBuffer<TEXTURE_TYPE>::IsReady();
}

template <GLenum TEXTURE_TYPE>
inline bool GLRenderDepth<TEXTURE_TYPE>::Resize( GLsizei width, GLsizei height, GLenum depthFormat )
{
	if ( TEXTURE_TYPE == GL_TEXTURE_CUBE_MAP ) {
		for ( int i=0; i<6; ++i ) {
			glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X+i,0,depthFormat,width,height,0,GL_DEPTH_COMPONENT,GL_FLOAT,nullptr);
		}
	} else {
		this->texture.SetImage(depthFormat,GL_DEPTH_COMPONENT,GL_FLOAT,nullptr,width,height);
	}
	GLRenderBuffer<TEXTURE_TYPE>::SetSize(width,height);
	return GLRenderBuffer<TEXTURE_TYPE>::IsComplete();
}

#endif
//-------------------------------------------------------------------------------
// GLSLShader Implementation
//-------------------------------------------------------------------------------

inline bool GLSLShader::CompileFile( char const *filename, GLenum shaderType, int prependSourceCount, char const **prependSources, std::ostream *outStream )
{
	std::ifstream shaderStream(filename, std::ios::in);
	if(! shaderStream.is_open()) {
		if ( outStream ) *outStream << "ERROR: Cannot open file." << std::endl;
		return false;
	}

	std::string shaderSourceCode((std::istreambuf_iterator<char>(shaderStream)), std::istreambuf_iterator<char>());
	shaderStream.close();

	return Compile( shaderSourceCode.data(), shaderType, prependSourceCount, prependSources, outStream );
}

inline bool GLSLShader::Compile( char const *shaderSourceCode, GLenum shaderType, int prependSourceCount, char const **prependSources, std::ostream *outStream )
{
	Delete();

	shaderID = glCreateShader( shaderType );
	if ( prependSourceCount > 0 ) {
		std::vector<char const*> sources(prependSourceCount+1);
		for ( int i=0; i<prependSourceCount; i++ ) sources[i] = prependSources[i];
		sources[prependSourceCount] = shaderSourceCode;
		glShaderSource(shaderID, prependSourceCount+1, sources.data(), nullptr);
	} else {
		glShaderSource(shaderID, 1, &shaderSourceCode, nullptr);
	}
	glCompileShader(shaderID);

	GLint result = GL_FALSE;
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &result);

	int infoLogLength;
	glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &infoLogLength);
	if ( infoLogLength > 1 ) {
		std::vector<char> compilerMessage(infoLogLength);
		glGetShaderInfoLog( shaderID, infoLogLength, nullptr, compilerMessage.data() );
		if ( outStream ) {
			if ( !result ) *outStream << "ERROR: Cannot compile shader." << std::endl;
			*outStream << "OpenGL Version: ";
			GL::PrintVersion(outStream);
			*outStream << std::endl;
			*outStream << compilerMessage.data() << std::endl;
		}
	}

	if ( result ) {
		GLint stype;
		glGetShaderiv(shaderID, GL_SHADER_TYPE, &stype);
		if ( stype != (GLint)shaderType ) {
			if ( outStream ) *outStream << "ERROR: Incorrect shader type." << std::endl;
			return false;
		}
	}

	return result == GL_TRUE;
}

//-------------------------------------------------------------------------------
// GLSLProgram Implementation
//-------------------------------------------------------------------------------

inline bool GLSLProgram::Link( std::ostream *outStream )
{
	glLinkProgram(programID);

	GLint result = GL_FALSE;
	glGetProgramiv(programID, GL_LINK_STATUS, &result);

	int infoLogLength;
	glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &infoLogLength);
	if ( infoLogLength > 1 ) {
		std::vector<char> compilerMessage(infoLogLength);
		glGetProgramInfoLog( programID, infoLogLength, nullptr, compilerMessage.data() );
		if ( outStream ) *outStream << "ERROR: " << compilerMessage.data() << std::endl;
	}

	return result == GL_TRUE;
}

inline bool GLSLProgram::Build( GLSLShader const *vertexShader, 
                                GLSLShader const *fragmentShader,
	                            GLSLShader const *geometryShader,
	                            GLSLShader const *tessControlShader,
	                            GLSLShader const *tessEvaluationShader,
	                            std::ostream *outStream )
{
	CreateProgram();
	std::stringstream output;
	AttachShader(*vertexShader);
	AttachShader(*fragmentShader);
	if ( geometryShader ) AttachShader(*geometryShader);
	if ( tessControlShader ) AttachShader(*tessControlShader);
	if ( tessEvaluationShader ) AttachShader(*tessEvaluationShader);
	return Link(outStream);
}

inline bool GLSLProgram::BuildFiles( char const *vertexShaderFile, 
                                     char const *fragmentShaderFile,
	                                 char const *geometryShaderFile,
	                                 char const *tessControlShaderFile,
	                                 char const *tessEvaluationShaderFile,
	                                 int         prependSourceCount,
	                                 char const **prependSource,
	                                 std::ostream *outStream )
{
	CreateProgram();
	GLSLShader vs, fs, gs, tcs, tes;
	std::stringstream shaderOutput;
	if ( ! vs.CompileFile(vertexShaderFile, GL_VERTEX_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
		if ( outStream ) *outStream << "ERROR: Failed compiling vertex shader \"" << vertexShaderFile << ".\"" << std::endl << shaderOutput.str();
		return false;
	}
	AttachShader(vs);
	if ( ! fs.CompileFile(fragmentShaderFile, GL_FRAGMENT_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
		if ( outStream ) *outStream << "ERROR: Failed compiling fragment shader \"" << fragmentShaderFile << ".\"" <<  std::endl << shaderOutput.str();
		return false;
	}
	AttachShader(fs);
	if ( geometryShaderFile ) {
		if ( ! gs.CompileFile(geometryShaderFile, GL_GEOMETRY_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling geometry shader \"" << geometryShaderFile << ".\"" <<  std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(gs);
	}
	if ( tessControlShaderFile ) {
		if ( ! tcs.CompileFile(tessControlShaderFile, GL_TESS_CONTROL_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling tessellation control shader \"" << tessControlShaderFile << ".\"" <<  std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(tcs);
	}
	if ( tessEvaluationShaderFile ) {
		if ( ! tes.CompileFile(tessEvaluationShaderFile, GL_TESS_EVALUATION_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling tessellation evaluation shader \"" << tessEvaluationShaderFile << ".\"" <<  std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(tes);
	}
	return Link(outStream);
}


inline bool GLSLProgram::BuildSources( char const *vertexShaderSourceCode, 
                                       char const *fragmentShaderSourceCode,
	                                   char const *geometryShaderSourceCode,
	                                   char const *tessControlShaderSourceCode,
	                                   char const *tessEvaluationShaderSourceCode,
	                                   int         prependSourceCount,
	                                   char const **prependSource,
	                                   std::ostream *outStream )
{
	CreateProgram();
	GLSLShader vs, fs, gs, tcs, tes;
	std::stringstream shaderOutput;
	if ( ! vs.Compile(vertexShaderSourceCode, GL_VERTEX_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
		if ( outStream ) *outStream << "ERROR: Failed compiling vertex shader." << std::endl << shaderOutput.str();
		return false;
	}
	AttachShader(vs);
	if ( ! fs.Compile(fragmentShaderSourceCode, GL_FRAGMENT_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
		if ( outStream ) *outStream << "ERROR: Failed compiling fragment shader." << std::endl << shaderOutput.str();
		return false;
	}
	AttachShader(fs);
	if ( geometryShaderSourceCode ) {
		if ( ! gs.Compile(geometryShaderSourceCode, GL_GEOMETRY_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling geometry shader." << std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(gs);
	}
	if ( tessControlShaderSourceCode ) {
		if ( ! tcs.Compile(tessControlShaderSourceCode, GL_TESS_CONTROL_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling tessellation control shader." << std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(tcs);
	}
	if ( tessEvaluationShaderSourceCode ) {
		if ( ! tes.Compile(tessEvaluationShaderSourceCode, GL_TESS_EVALUATION_SHADER, prependSourceCount, prependSource, &shaderOutput) ) {
			if ( outStream ) *outStream << "ERROR: Failed compiling tessellation evaluation shader." << std::endl << shaderOutput.str();
			return false;
		}
		AttachShader(tes);
	}
	return Link(outStream);
}

inline void GLSLProgram::RegisterUniform( unsigned int index, char const *name, std::ostream *outStream )
{
	if ( params.size() <= index ) params.resize( index+1, -1 );
	params[index] = glGetUniformLocation( programID, name );
	if ( params[index] < 0 ) {
		GLenum error = glGetError();
		GLenum newError;
		while ( (newError = glGetError()) != GL_NO_ERROR ) error = newError; // get the latest error.
		if ( outStream ) {
			*outStream << "ERROR: ";
			switch (error) {
				case GL_INVALID_VALUE:     *outStream << "GL_INVALID_VALUE.";     break;
				case GL_INVALID_OPERATION: *outStream << "GL_INVALID_OPERATION."; break;
			}
			*outStream << " Parameter \"" << name << "\" could not be registered." << std::endl;
		}
	}
}

inline void GLSLProgram::RegisterUniforms( char const *names, unsigned int startingIndex, std::ostream *outStream )
{
	std::stringstream ss(names);
	unsigned int index = startingIndex;
	while ( ss.good() ) {
		std::string name;
		ss >> name;
		RegisterUniform( index++, name.c_str(), outStream );
	}
}

//-------------------------------------------------------------------------------

typedef GLTexture1<GL_TEXTURE_1D       >     GLTexture1D;			//!< OpenGL 1D Texture
typedef GLTexture2<GL_TEXTURE_2D       >     GLTexture2D;			//!< OpenGL 2D Texture
typedef GLTexture3<GL_TEXTURE_3D       >     GLTexture3D;			//!< OpenGL 3D Texture

#ifdef GL_TEXTURE_1D_ARRAY
typedef GLTexture2<GL_TEXTURE_1D_ARRAY >     GLTexture1DArray;		//!< OpenGL 1D Texture Array
typedef GLTexture3<GL_TEXTURE_2D_ARRAY >     GLTexture2DArray;		//!< OpenGL 2D Texture Array
#endif
#ifdef GL_TEXTURE_RECTANGLE
typedef GLTexture2<GL_TEXTURE_RECTANGLE>     GLTextureRect;			//!< OpenGL Rectangle Texture
typedef GLTexture1<GL_TEXTURE_BUFFER   >     GLTextureBuffer;		//!< OpenGL Buffer Texture
#endif

typedef GLRenderTexture<GL_TEXTURE_2D>        GLRenderTexture2D;	//!< OpenGL render color buffer with a 2D texture
typedef GLRenderDepth  <GL_TEXTURE_2D>        GLRenderDepth2D;		//!< OpenGL render depth buffer with a 2D texture
#ifdef GL_TEXTURE_RECTANGLE
typedef GLRenderTexture<GL_TEXTURE_RECTANGLE> GLRenderTextureRect;	//!< OpenGL render color buffer with a rectangle texture
typedef GLRenderDepth  <GL_TEXTURE_RECTANGLE> GLRenderDepthRect;	//!< OpenGL render depth buffer with a rectangle texture
#endif
typedef GLRenderTextureCubeBase< GL_COLOR_ATTACHMENT0, GLRenderTexture<GL_TEXTURE_CUBE_MAP> > GLRenderTextureCube;	//!< OpenGL render color buffer with a cube map texture
typedef GLRenderTextureCubeBase< GL_DEPTH_ATTACHMENT,  GLRenderDepth  <GL_TEXTURE_CUBE_MAP> > GLRenderDepthCube;	//!< OpenGL render depth buffer with a cube map texture

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::GL cyGL;									//!< General OpenGL queries

#ifdef GL_KHR_debug
typedef cy::GLDebugCallback    cyGLDebugCallback;		//!< OpenGL debug callback class
#endif

typedef cy::GLTextureCubeMapSide cyGLTextureCubeMapSide;	//! Sides of a cube map

typedef cy::GLTexture1D        cyGLTexture1D;			//!< OpenGL 1D Texture
typedef cy::GLTexture2D        cyGLTexture2D;			//!< OpenGL 2D Texture
typedef cy::GLTexture3D        cyGLTexture3D;			//!< OpenGL 3D Texture
typedef cy::GLTextureCubeMap   cyGLTextureCubeMap;		//!< OpenGL Cube Map Texture

#ifdef GL_TEXTURE_1D_ARRAY
typedef cy::GLTexture1DArray   cyGLTexture1DArray;		//!< OpenGL 1D Texture Array
typedef cy::GLTexture2DArray   cyGLTexture2DArray;		//!< OpenGL 2D Texture Array
#endif
#ifdef GL_TEXTURE_RECTANGLE
typedef cy::GLTextureRect      cyGLTextureRect;			//!< OpenGL Rectangle Texture
typedef cy::GLTextureBuffer    cyGLTextureBuffer;		//!< OpenGL Buffer Texture
#endif

typedef cy::GLRenderTexture2D   cyGLRenderTexture2D;	//!< OpenGL render color buffer with a 2D texture
typedef cy::GLRenderDepth2D     cyGLRenderDepth2D;		//!< OpenGL render depth buffer with a 2D texture
#ifdef GL_TEXTURE_RECTANGLE
typedef cy::GLRenderTextureRect cyGLRenderTextureRect;	//!< OpenGL render color buffer with a rectangle texture
typedef cy::GLRenderDepthRect   cyGLRenderDepthRect;	//!< OpenGL render depth buffer with a rectangle texture
#endif
typedef cy::GLRenderTextureCube cyGLRenderTextureCube;	//!< OpenGL render color buffer with a cube map texture
typedef cy::GLRenderDepthCube   cyGLRenderDepthCube;	//!< OpenGL render depth buffer with a cube map texture


typedef cy::GLSLShader         cyGLSLShader;			//!< GLSL shader class
typedef cy::GLSLProgram        cyGLSLProgram;			//!< GLSL program class

//-------------------------------------------------------------------------------
#endif
