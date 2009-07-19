
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_TEXTURE_H
#define PBRT_TEXTURE_H
// texture.h*
#include "pbrt.h"
#include "color.h"
#include "geometry.h"
#include "transform.h"
// Texture Declarations
class COREDLL TextureMapping2D {
public:
	// TextureMapping2D Interface
	virtual ~TextureMapping2D() { }
	virtual void Map(const DifferentialGeometry &dg,
		float *s, float *t, float *dsdx, float *dtdx,
		float *dsdy, float *dtdy) const = 0;
};
class COREDLL UVMapping2D : public TextureMapping2D {
public:
	// UVMapping2D Public Methods
	UVMapping2D(float su = 1, float sv = 1,
		float du = 0, float dv = 0);
	void Map(const DifferentialGeometry &dg, float *s, float *t,
		float *dsdx, float *dtdx,
		float *dsdy, float *dtdy) const;
private:
	float su, sv, du, dv;
};
class COREDLL SphericalMapping2D : public TextureMapping2D {
public:
	// SphericalMapping2D Public Methods
	SphericalMapping2D(const Transform &toSph)
		: WorldToTexture(toSph) {
	}
	void Map(const DifferentialGeometry &dg, float *s, float *t,
		float *dsdx, float *dtdx,
		float *dsdy, float *dtdy) const;
private:
	void sphere(const Point &P, float *s, float *t) const;
	Transform WorldToTexture;
};
class
COREDLL CylindricalMapping2D : public TextureMapping2D {
public:
	// CylindricalMapping2D Public Methods
	CylindricalMapping2D(const Transform &toCyl)
		: WorldToTexture(toCyl) {
	}
	void Map(const DifferentialGeometry &dg, float *s, float *t,
		float *dsdx, float *dtdx,
		float *dsdy, float *dtdy) const;
private:
	void cylinder(const Point &P, float *s, float *t) const;
	Transform WorldToTexture;
};
class COREDLL PlanarMapping2D : public TextureMapping2D {
public:
	// PlanarMapping2D Public Methods
	PlanarMapping2D(const Vector &v1, const Vector &v2,
		float du = 0, float dv = 0);
	void Map(const DifferentialGeometry &dg, float *s, float *t,
		float *dsdx, float *dtdx,
		float *dsdy, float *dtdy) const;
private:
	Vector vs, vt;
	float ds, dt;
};
class COREDLL TextureMapping3D {
public:
	// TextureMapping3D Interface
	virtual ~TextureMapping3D() { }
	virtual Point Map(const DifferentialGeometry &dg,
		Vector *dpdx, Vector *dpdy) const = 0;
};
class COREDLL IdentityMapping3D : public TextureMapping3D {
public:
	IdentityMapping3D(const Transform &x)
		: WorldToTexture(x) { }
	Point Map(const DifferentialGeometry &dg, Vector *dpdx,
		Vector *dpdy) const;
private:
	Transform WorldToTexture;
};
template <class T> class Texture : public ReferenceCounted {
public:
	// Texture Interface
	virtual T Evaluate(const DifferentialGeometry &) const = 0;
	virtual ~Texture() { }
};
// ConstantTexture Declarations
template <class T>
class ConstantTexture : public Texture<T> {
public:
	// ConstantTexture Public Methods
	ConstantTexture(const T &v) { value = v; }
	T Evaluate(const DifferentialGeometry &) const {
		return value;
	}
private:
	T value;
};
#endif // PBRT_TEXTURE_H
