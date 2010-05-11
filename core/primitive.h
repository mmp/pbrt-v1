
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

#ifndef PBRT_PRIMITIVE_H
#define PBRT_PRIMITIVE_H
// primitive.h*
#include "pbrt.h"
#include "shape.h"
#include "material.h"
// Primitive Declarations
class COREDLL Primitive : public ReferenceCounted {
public:
	// Primitive Interface
	virtual ~Primitive();
	virtual BBox WorldBound() const = 0;
	virtual bool CanIntersect() const;
	virtual bool Intersect(const Ray &r,
		Intersection *in) const = 0;
	virtual bool IntersectP(const Ray &r) const = 0;
	virtual void
		Refine(vector<Reference<Primitive> > &refined) const;
	void FullyRefine(vector<Reference<Primitive> > &refined)
	const;
	virtual const AreaLight *GetAreaLight() const = 0;
	virtual BSDF *GetBSDF(const DifferentialGeometry &dg,
		const Transform &WorldToObject) const = 0;
};
struct COREDLL Intersection {
	// Intersection Public Methods
	Intersection() { primitive = NULL; }
	BSDF *GetBSDF(const RayDifferential &ray) const;
	Spectrum Le(const Vector &wo) const;
	DifferentialGeometry dg;
	const Primitive *primitive;
	Transform WorldToObject;
};
class COREDLL GeometricPrimitive : public Primitive {
public:
	// GeometricPrimitive Public Methods
	bool CanIntersect() const;
	void Refine(vector<Reference<Primitive> > &refined) const;
	virtual BBox WorldBound() const;
	virtual bool Intersect(const Ray &r,
	                       Intersection *isect) const;
	virtual bool IntersectP(const Ray &r) const;
	GeometricPrimitive(const Reference<Shape> &s,
	                   const Reference<Material> &m,
	                   AreaLight *a);
	const AreaLight *GetAreaLight() const;
	BSDF *GetBSDF(const DifferentialGeometry &dg,
	              const Transform &WorldToObject) const;
private:
	// GeometricPrimitive Private Data
	Reference<Shape> shape;
	Reference<Material> material;
	AreaLight *areaLight;
};
class COREDLL InstancePrimitive : public Primitive {
public:
	// InstancePrimitive Public Methods
	InstancePrimitive(Reference<Primitive> &i,
	                  const Transform &i2w) {
		instance = i;
		InstanceToWorld = i2w;
		WorldToInstance = i2w.GetInverse();
	}
	bool Intersect(const Ray &r, Intersection *in) const;
	bool IntersectP(const Ray &r) const;
	const AreaLight *GetAreaLight() const { return NULL; }
	BSDF *GetBSDF(const DifferentialGeometry &dg,
	              const Transform &WorldToObject) const {
		return NULL;
	}
	BBox WorldBound() const {
		return InstanceToWorld(instance->WorldBound());
	}
private:
	// InstancePrimitive Private Data
	Reference<Primitive> instance;
	Transform InstanceToWorld, WorldToInstance;
};
class COREDLL Aggregate : public Primitive {
public:
	// Aggregate Public Methods
	const AreaLight *GetAreaLight() const;
	BSDF *GetBSDF(const DifferentialGeometry &dg,
	              const Transform &) const;
};
#endif // PBRT_PRIMITIVE_H
