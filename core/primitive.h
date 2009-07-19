
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
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
