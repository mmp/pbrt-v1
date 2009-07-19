
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_SHAPE_H
#define PBRT_SHAPE_H
// shape.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
#include "paramset.h"
// DifferentialGeometry Declarations
struct COREDLL DifferentialGeometry {
	DifferentialGeometry() { u = v = 0.; shape = NULL; }
	// DifferentialGeometry Public Methods
	DifferentialGeometry(const Point &P, const Vector &DPDU,
			const Vector &DPDV, const Vector &DNDU,
			const Vector &DNDV, float uu, float vv,
			const Shape *sh);
	void ComputeDifferentials(const RayDifferential &r) const;
	// DifferentialGeometry Public Data
	Point p;
	Normal nn;
	float u, v;
	const Shape *shape;
	Vector dpdu, dpdv;
	Normal dndu, dndv;
	mutable Vector dpdx, dpdy;
	mutable float dudx, dvdx, dudy, dvdy;
};
// Shape Declarations
class COREDLL Shape : public ReferenceCounted {
public:
	// Shape Interface
	Shape(const Transform &o2w, bool ro);
	virtual ~Shape() { }
	virtual BBox ObjectBound() const = 0;
	virtual BBox WorldBound() const {
		return ObjectToWorld(ObjectBound());
	}
	virtual bool CanIntersect() const { return true; }
	virtual void
		Refine(vector<Reference<Shape> > &refined) const {
		Severe("Unimplemented Shape::Refine() method called");
	}
	virtual bool Intersect(const Ray &ray, float *tHit,
			DifferentialGeometry *dg) const {
		Severe("Unimplemented Shape::Intersect()"
	           "method called");
		return false;
	}
	virtual bool IntersectP(const Ray &ray) const {
		Severe("Unimplemented Shape::IntersectP()"
	           "method called");
		return false;
	}
	virtual void GetShadingGeometry(const Transform &obj2world,
			const DifferentialGeometry &dg,
			DifferentialGeometry *dgShading) const {
		*dgShading = dg;
	}
	virtual float Area() const {
		Severe("Unimplemented Shape::Area() method called");
		return 0.;
	}
	virtual Point Sample(float u1, float u2, Normal *Ns) const {
		Severe("Unimplemented Shape::Sample method called");
		return Point();
	}
	virtual float Pdf(const Point &Pshape) const {
		return 1.f / Area();
	}
	virtual Point Sample(const Point &P,
			float u1, float u2, Normal *Ns) const {
		return Sample(u1, u2, Ns);
	}
	virtual float Pdf(const Point &p, const Vector &wi) const {
		// Intersect sample ray with area light geometry
		DifferentialGeometry dgLight;
		Ray ray(p, wi);
		float thit;
		if (!Intersect(ray, &thit, &dgLight)) return 0.;
		// Convert light sample weight to solid angle measure
		float pdf = DistanceSquared(p, ray(thit)) /
			(AbsDot(dgLight.nn, -wi) * Area());
		if (AbsDot(dgLight.nn, -wi) == 0.f) pdf = INFINITY; // NOBOOK
		return pdf;
	}
	// Shape Public Data
	const Transform ObjectToWorld, WorldToObject;
	const bool reverseOrientation, transformSwapsHandedness;
};
class ShapeSet : public Shape {
public:
	// ShapeSet Public Methods
	Point Sample(float u1, float u2, Normal *Ns) const {
		float ls = RandomFloat();
		u_int sn;
		for (sn = 0; sn < shapes.size()-1; ++sn)
			if (ls < areaCDF[sn]) break;
		return shapes[sn]->Sample(u1, u2, Ns);
	}
	ShapeSet(const vector<Reference<Shape> > &s,
		const Transform &o2w, bool ro)
		: Shape(o2w, ro) {
		shapes = s;
		area = 0;
		vector<float> areas;
		for (u_int i = 0; i < shapes.size(); ++i) {
			float a = shapes[i]->Area();
			area += a;
			areas.push_back(a);
		}
		float prevCDF = 0;
		for (u_int i = 0; i < shapes.size(); ++i) {
			areaCDF.push_back(prevCDF + areas[i] / area);
			prevCDF = areaCDF[i];
		}
	}
	BBox ObjectBound() const {
		BBox ob;
		for (u_int i = 0; i < shapes.size(); ++i)
			ob = Union(ob, shapes[i]->ObjectBound());
		return ob;
	}
	bool CanIntersect() const {
		for (u_int i = 0; i < shapes.size(); ++i)
			if (!shapes[i]->CanIntersect()) return false;
		return true;
	}
	bool Intersect(const Ray &ray, float *t_hitp,
			DifferentialGeometry *dg) const {
		bool anyHit = false;
		for (u_int i = 0; i < shapes.size(); ++i)
			if (shapes[i]->Intersect(ray, t_hitp, dg)) anyHit = true;
		return anyHit;
	}
	void Refine(vector<Reference<Shape> > &refined) const {
		for (u_int i = 0; i < shapes.size(); ++i) {
			if (shapes[i]->CanIntersect())
				refined.push_back(shapes[i]);
			else shapes[i]->Refine(refined);
		}
	
	}
	float Area() const { return area; }
private:
	// ShapeSet Private Data
	float area;
	vector<float> areaCDF;
	vector<Reference<Shape> > shapes;
};
#endif // PBRT_SHAPE_H
