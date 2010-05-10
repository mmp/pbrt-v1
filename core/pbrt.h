
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef PBRT_PBRT_H
#define PBRT_PBRT_H
// pbrt.h*
// Global Include Files
#if !defined(__APPLE__) && !defined(__OpenBSD__)
#include <malloc.h> // for _alloca, memalign
#endif
#if !defined(WIN32) && !defined(__APPLE__) && !defined(__OpenBSD__)
#include <alloca.h>
#endif
#ifdef __linux__
#include <fpu_control.h>
#endif
#ifdef WIN32
#include <float.h>
#include <float.h>
#define isnan _isnan
#define isinf(f) (!_finite((f)))
#endif
#include <math.h>
#include <stdlib.h>
#define _GNU_SOURCE 1 //NOBOOK
#include <stdio.h>
#include <string.h>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <iostream>
using std::ostream;
#include <algorithm>
using std::min;
using std::max;
using std::swap;
using std::sort;
#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
#include <assert.h>
// Platform-specific definitions
#if defined(WIN32)
#define memalign(a,b) _aligned_malloc(b, a)
#elif defined(__APPLE__)
#define memalign(a,b) valloc(b)
#elif defined(__OpenBSD__)
#define memalign(a,b) malloc(b)
#endif
#ifdef sgi
#define for if (0) ; else for
#endif
#ifdef __APPLE__
#define powf pow
#define sinf sin
#define cosf cos
#define tanf tan
#define asinf asin
#define acosf acos
#define atanf atan
#define atan2f atan2
#define logf log
#define log10f log10
#define expf exp
#define sqrtf sqrt
#if __GNUC__ == 3
extern "C" {
  int isinf(double);
  int isnan(double);
}
#endif // ONLY GCC 3
#endif // __APPLE__

#ifdef WIN32
#pragma warning (disable: 4267 4251 4065 4102)
#endif // WIN32
#ifdef WIN32
#pragma warning( disable: 4190 )
// extern "C" nonsense when returning a template
#endif
#ifdef WIN32
#ifdef CORE_SOURCE
#define COREDLL __declspec(dllexport)
#else
#define COREDLL __declspec(dllimport)
#endif
#define DLLEXPORT __declspec(dllexport)
#else
#define COREDLL
#define DLLEXPORT
#endif
#ifdef WIN32
#define PBRT_PATH_SEP ";"
#else
#define PBRT_PATH_SEP ":"
#endif
#ifdef WIN32
#define isnan _isnan
#define isinf(f) (!_finite((f)))
#endif
// Global Type Declarations
typedef double StatsCounterType;
typedef unsigned char u_char;
typedef unsigned short u_short;
typedef unsigned int u_int;
typedef unsigned long u_long;
// Global Forward Declarations
class Timer;
template<class T, int logBlockSize = 2> class BlockedArray;
struct Matrix4x4;
class ParamSet;
template <class T> struct ParamSetItem;
class TextureParams;
class Shape;
class Scene;
class Vector;
class Point;
class Normal;
class Ray;
class RayDifferential;
class BBox;
class Transform;
struct DifferentialGeometry;
class Primitive;
struct Intersection;
class GeometricPrimitive;
class Spectrum;
class Camera;
class ProjectiveCamera;
class Sampler;
struct Sample;
#define BC_GRID_SIZE 40
typedef vector<int> SampleGrid[BC_GRID_SIZE][BC_GRID_SIZE];
#define GRID(v) (int((v) * BC_GRID_SIZE))
class Filter;
class Film;
class ToneMap;
class BxDF;
class BRDF;
class BTDF;
class Fresnel;
class FresnelConductor;
class FresnelDielectric;
class FresnelNoOp;
class SpecularReflection;
class SpecularTransmission;
class Lambertian;
class OrenNayar;
class Microfacet;
class MicrofacetDistribution;
class BSDF;
class Material;
class TextureMapping2D;
class UVMapping2D;
class SphericalMapping2D;
class CylindricalMapping2D;
class PlanarMapping2D;
class TextureMapping3D;
class IdentityMapping3D;
template <class T> class Texture;
COREDLL float Noise(float x, float y = .5f, float z = .5f);
COREDLL float Noise(const Point &P);
COREDLL float FBm(const Point &P, const Vector &dpdx, const Vector &dpdy,
	float omega, int octaves);
COREDLL float Turbulence(const Point &P, const Vector &dpdx, const Vector &dpdy,
	float omega, int octaves);
class VolumeRegion;
class Light;
struct VisibilityTester;
class AreaLight;
class ShapeSet;
class SurfaceIntegrator;
class Integrator;
class VolumeIntegrator;
// Global Constants
#ifdef WIN32
#define alloca _alloca
#endif
#ifdef M_PI
#undef M_PI
#endif
#define M_PI           3.14159265358979323846f
#define INV_PI  0.31830988618379067154f
#define INV_TWOPI  0.15915494309189533577f
#ifndef INFINITY
#define INFINITY FLT_MAX
#endif
#define PBRT_VERSION 1.05
#define RAY_EPSILON 1e-3f
#define COLOR_SAMPLES 3
// Global Function Declarations
// Setup printf format
#ifdef __GNUG__
#define PRINTF_FUNC __attribute__ \
	((__format__ (__printf__, 1, 2)))
#else
#define PRINTF_FUNC
#endif // __GNUG__
extern COREDLL void Info(const char *, ...) PRINTF_FUNC;
extern COREDLL void Warning(const char *, ...) PRINTF_FUNC;
extern COREDLL void Error(const char *, ...) PRINTF_FUNC;
extern COREDLL void Severe(const char *, ...) PRINTF_FUNC;
extern void StatsPrint(FILE *dest);
extern void StatsCleanup();
COREDLL void *AllocAligned(size_t size);
COREDLL void FreeAligned(void *);
COREDLL bool SolveLinearSystem2x2(const float A[2][2], const float B[2],
	float x[2]);
COREDLL unsigned long genrand_int32(void);
COREDLL extern float genrand_real1(void);
COREDLL extern float genrand_real2(void);
COREDLL Spectrum *ReadImage(const string &name, int *xSize,
	int *ySize);
COREDLL void WriteRGBAImage(const string &name,
	float *pixels, float *alpha, int XRes, int YRes,
	int totalXRes, int totalYRes, int xOffset, int yOffset);
COREDLL void pbrtInit();
COREDLL void pbrtCleanup();
COREDLL Transform Translate(const Vector &delta);
COREDLL Transform Scale(float x, float y, float z);
extern COREDLL Transform RotateX(float angle);
extern COREDLL Transform RotateY(float angle);
extern COREDLL Transform RotateZ(float angle);
extern COREDLL Transform Rotate(float angle, const Vector &axis);
extern COREDLL Transform LookAt(const Point &pos, const Point &look,
                        const Vector &up);
COREDLL Transform Orthographic(float znear, float zfar);
COREDLL
Transform Perspective(float fov, float znear, float zfar);
COREDLL Spectrum FrDiel(float cosi, float cost,
                        const Spectrum &etai,
						const Spectrum &etat);
COREDLL Spectrum FrCond(float cosi,
                        const Spectrum &n,
						const Spectrum &k);
COREDLL
	Spectrum FresnelApproxEta(const Spectrum &intensity);
COREDLL
	Spectrum FresnelApproxK(const Spectrum &intensity);
COREDLL float Lanczos(float, float tau=2);
COREDLL Spectrum EstimateDirect(const Scene *scene, const Light *light, const Point &p,
	const Normal &n, const Vector &wo, BSDF *bsdf,
	const Sample *sample, int lightSamp, int bsdfSamp,
	int bsdfComponent, u_int sampleNum);
extern COREDLL void ComputeStep1dCDF(float *f, int nValues, float *c, float *cdf);
extern COREDLL float SampleStep1d(float *f, float *cdf, float c,
	int nSteps, float u, float *weight);
COREDLL void ConcentricSampleDisk(float u1, float u2, float *dx, float *dy);
COREDLL void UniformSampleTriangle(float ud1, float ud2, float *u, float *v);
COREDLL bool ParseFile(const char *filename);
// Global Classes
struct COREDLL ProgressReporter {
	// ProgressReporter Public Methods
	ProgressReporter(int totalWork, const string &title,
		int barLength=58);
	~ProgressReporter();
	void Update(int num = 1) const;
	void Done() const;
	// ProgressReporter Data
	const int totalPlusses;
	float frequency;
	mutable float count;
	mutable int plussesPrinted;
	mutable Timer *timer;
	FILE *outFile;
	char *buf;
	mutable char *curSpace;
};
class COREDLL StatsCounter {
public:
	// StatsCounter Public Methods
	StatsCounter(const string &category, const string &name);
	void operator++() { ++num; }
	void operator++(int) { ++num; }
	void Max(StatsCounterType val) { num = max(val, num); }
	void Min(StatsCounterType val) { num = min(val, num); }
	operator double() const { return (double)num; }
private:
	// StatsCounter Private Data
	StatsCounterType num;
};
class COREDLL StatsRatio {
public:
	// StatsRatio Public Methods
	StatsRatio(const string &category, const string &name);
	void Add(int a, int b) { na += a; nb += b; }
private:
	// StatsRatio Private Data
	StatsCounterType na, nb;
};
class COREDLL StatsPercentage {
public:
	// StatsPercentage Public Methods
	void Add(int a, int b) { na += a; nb += b; }
	StatsPercentage(const string &category, const string &name);
private:
	// StatsPercentage Private Data
	StatsCounterType na, nb;
};
class COREDLL ReferenceCounted {
public:
	ReferenceCounted() { nReferences = 0; }
	int nReferences;
private:
	ReferenceCounted(const ReferenceCounted &);
	ReferenceCounted &operator=(const ReferenceCounted &);
};
template <class T> class Reference {
public:
	// Reference Public Methods
	Reference(T *p = NULL) {
		ptr = p;
		if (ptr) ++ptr->nReferences;
	}
	Reference(const Reference<T> &r) {
		ptr = r.ptr;
		if (ptr) ++ptr->nReferences;
	}
	Reference &operator=(const Reference<T> &r) {
		if (r.ptr) r.ptr->nReferences++;
		if (ptr && --ptr->nReferences == 0) delete ptr;
		ptr = r.ptr;
		return *this;
	}
	Reference &operator=(T *p) {
		if (p) p->nReferences++;
		if (ptr && --ptr->nReferences == 0) delete ptr;
		ptr = p;
		return *this;
	}
	~Reference() {
		if (ptr && --ptr->nReferences == 0)
			delete ptr;
	}
	T *operator->() { return ptr; }
	const T *operator->() const { return ptr; }
	operator bool() const { return ptr != NULL; }
	bool operator<(const Reference<T> &t2) const {
		return ptr < t2.ptr;
	}
private:
	T *ptr;
};
template <class T> class ObjectArena {
public:
	// ObjectArena Public Methods
	ObjectArena() {
		nAvailable = 0;
	}
	T *Alloc() {
		if (nAvailable == 0) {
			int nAlloc = max((unsigned long)16,
				(unsigned long)(65536/sizeof(T)));
			mem = (T *)AllocAligned(nAlloc * sizeof(T));
			nAvailable = nAlloc;
			toDelete.push_back(mem);
		}
		--nAvailable;
		return mem++;
	}
	operator T *() {
		return Alloc();
	}
	~ObjectArena() { FreeAll(); }
	void FreeAll() {
		for (u_int i = 0; i < toDelete.size(); ++i)
			FreeAligned(toDelete[i]);
		toDelete.erase(toDelete.begin(), toDelete.end());
		nAvailable = 0;
	}
private:
	// ObjectArena Private Data
	T *mem;
	int nAvailable;
	vector<T *> toDelete;
};
class COREDLL MemoryArena {
public:
	// MemoryArena Public Methods
	MemoryArena(u_int bs = 32768) {
		blockSize = bs;
		curBlockPos = 0;
		currentBlock = (char *)AllocAligned(blockSize);
	}
	~MemoryArena() {
		FreeAligned(currentBlock);
		for (u_int i = 0; i < usedBlocks.size(); ++i)
			FreeAligned(usedBlocks[i]);
		for (u_int i = 0; i < availableBlocks.size(); ++i)
			FreeAligned(availableBlocks[i]);
	}
	void *Alloc(u_int sz) {
		// Round up _sz_ to minimum machine alignment
		sz = ((sz + 7) & (~7));
		if (curBlockPos + sz > blockSize) {
			// Get new block of memory for _MemoryArena_
			usedBlocks.push_back(currentBlock);
			if (availableBlocks.size() && sz <= blockSize) {
				currentBlock = availableBlocks.back();
				availableBlocks.pop_back();
			}
			else
				currentBlock = (char *)AllocAligned(max(sz, blockSize));
			curBlockPos = 0;
		}
		void *ret = currentBlock + curBlockPos;
		curBlockPos += sz;
		return ret;
	}
	void FreeAll() {
		curBlockPos = 0;
		while (usedBlocks.size()) {
			availableBlocks.push_back(usedBlocks.back());
			usedBlocks.pop_back();
		}
	}
private:
	// MemoryArena Private Data
	u_int curBlockPos, blockSize;
	char *currentBlock;
	vector<char *> usedBlocks, availableBlocks;
};
template<class T, int logBlockSize> class BlockedArray {
public:
	// BlockedArray Public Methods
	BlockedArray(int nu, int nv, const T *d = NULL) {
		uRes = nu;
		vRes = nv;
		uBlocks = RoundUp(uRes) >> logBlockSize;
		int nAlloc = RoundUp(uRes) * RoundUp(vRes);
		data = (T *)AllocAligned(nAlloc * sizeof(T));
		for (int i = 0; i < nAlloc; ++i)
			new (&data[i]) T();
		if (d)
			for (int v = 0; v < nv; ++v)
				for (int u = 0; u < nu; ++u)
					(*this)(u, v) = d[v * uRes + u];
	}
	int BlockSize() const { return 1 << logBlockSize; }
	int RoundUp(int x) const {
		return (x + BlockSize() - 1) & ~(BlockSize() - 1);
	}
	int uSize() const { return uRes; }
	int vSize() const { return vRes; }
	~BlockedArray() {
		for (int i = 0; i < uRes * vRes; ++i)
			data[i].~T();
		FreeAligned(data);
	}
	int Block(int a) const { return a >> logBlockSize; }
	int Offset(int a) const { return (a & (BlockSize() - 1)); }
	T &operator()(int u, int v) {
		int bu = Block(u), bv = Block(v);
		int ou = Offset(u), ov = Offset(v);
		int offset = BlockSize() * BlockSize() *
		             (uBlocks * bv + bu);
		offset += BlockSize() * ov + ou;
		return data[offset];
	}
	const T &operator()(int u, int v) const {
		int bu = Block(u), bv = Block(v);
		int ou = Offset(u), ov = Offset(v);
		int offset = BlockSize() * BlockSize() * (uBlocks * bv + bu);
		offset += BlockSize() * ov + ou;
		return data[offset];
	}
	void GetLinearArray(T *a) const {
		for (int v = 0; v < vRes; ++v)
			for (int u = 0; u < uRes; ++u)
				*a++ = (*this)(u, v);
	}
private:
	// BlockedArray Private Data
	T *data;
	int uRes, vRes, uBlocks;
};
struct COREDLL Matrix4x4 : public ReferenceCounted {
	// Matrix4x4 Public Methods
	Matrix4x4() {
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				if (i == j) m[i][j] = 1.;
				else m[i][j] = 0.;
	}
	Matrix4x4(float mat[4][4]);
	Matrix4x4(float t00, float t01, float t02, float t03,
	          float t10, float t11, float t12, float t13,
	          float t20, float t21, float t22, float t23,
	          float t30, float t31, float t32, float t33);
	Reference<Matrix4x4> Transpose() const;
	void Print(ostream &os) const {
		os << "[ ";
		for (int i = 0; i < 4; ++i) {
			os << "[ ";
			for (int j = 0; j < 4; ++j)  {
				os << m[i][j];
				if (j != 3) os << ", ";
			}
			os << " ] ";
		}
		os << " ] ";
	}
	static Reference<Matrix4x4>
		Mul(const Reference<Matrix4x4> &m1,
	        const Reference<Matrix4x4> &m2) {
		float r[4][4];
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				r[i][j] = m1->m[i][0] * m2->m[0][j] +
				          m1->m[i][1] * m2->m[1][j] +
				          m1->m[i][2] * m2->m[2][j] +
				          m1->m[i][3] * m2->m[3][j];
		return new Matrix4x4(r);
	}
	Reference<Matrix4x4> Inverse() const;
	float m[4][4];
};
// Global Inline Functions
#ifdef NDEBUG
#define Assert(expr) ((void)0)
#else
#define Assert(expr) \
    ((expr) ? (void)0 : \
		Severe("Assertion \"%s\" failed in %s, line %d", \
               #expr, __FILE__, __LINE__))
#endif // NDEBUG
inline float Lerp(float t, float v1, float v2) {
	return (1.f - t) * v1 + t * v2;
}
inline float Clamp(float val, float low, float high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}
inline int Clamp(int val, int low, int high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}
inline int Mod(int a, int b) {
    int n = int(a/b);
    a -= n*b;
    if (a < 0)
        a += b;
    return a;
}
inline float Radians(float deg) {
	return ((float)M_PI/180.f) * deg;
}
inline float Degrees(float rad) {
	return (180.f/(float)M_PI) * rad;
}
inline float Log2(float x) {
	static float invLog2 = 1.f / logf(2.f);
	return logf(x) * invLog2;
}
inline int Log2Int(float v) {
#if 0
 	return ((*reinterpret_cast<int *>(&v)) >> 23) - 127;
#else
#define _doublemagicroundeps	      (.5-1.4e-11)
	return int(Log2(v) + _doublemagicroundeps);
#endif
}
inline bool IsPowerOf2(int v) {
	return (v & (v - 1)) == 0;
}
inline u_int RoundUpPow2(u_int v) {
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return v+1;
}
#if (defined(__linux__) && defined(__i386__)) || defined(WIN32)
//#define FAST_INT 1
#endif
#define _doublemagicroundeps	      (.5-1.4e-11)
	//almost .5f = .5f - 1e^(number of exp bit)
inline int Round2Int(double val) {
#ifdef FAST_INT
#define _doublemagic			double (6755399441055744.0)
	//2^52 * 1.5,  uses limited precision to floor
	val		= val + _doublemagic;
	return (reinterpret_cast<long*>(&val))[0];
#else
	return int (val+_doublemagicroundeps);
#endif
}
inline int Float2Int(double val) {
#ifdef FAST_INT
	return (val<0) ?  Round2Int(val+_doublemagicroundeps) :
		   Round2Int(val-_doublemagicroundeps);
#else
	return (int)val;
#endif
}
inline int Floor2Int(double val) {
#ifdef FAST_INT
	return Round2Int(val - _doublemagicroundeps);
#else
	return (int)floor(val);
#endif
}
inline int Ceil2Int(double val) {
#ifdef FAST_INT
	return Round2Int(val + _doublemagicroundeps);
#else
	return (int)ceil(val);
#endif
}
inline float RandomFloat();
inline unsigned long RandomUInt();
inline float RandomFloat() {
	return genrand_real2();
}

inline unsigned long RandomUInt() {
	return genrand_int32();
}
inline bool Quadratic(float A, float B, float C, float *t0,
		float *t1) {
	// Find quadratic discriminant
	float discrim = B * B - 4.f * A * C;
	if (discrim < 0.) return false;
	float rootDiscrim = sqrtf(discrim);
	// Compute quadratic _t_ values
	float q;
	if (B < 0) q = -.5f * (B - rootDiscrim);
	else       q = -.5f * (B + rootDiscrim);
	*t0 = q / A;
	*t1 = C / q;
	if (*t0 > *t1) swap(*t0, *t1);
	return true;
}
inline float SmoothStep(float min, float max, float value) {
	float v = Clamp((value - min) / (max - min), 0.f, 1.f);
	return v * v * (-2.f * v  + 3.f);
}
inline float ExponentialAverage(float avg,
                              float val, float alpha) {
	return (1.f - alpha) * val + alpha * avg;
}
#endif // PBRT_PBRT_H
