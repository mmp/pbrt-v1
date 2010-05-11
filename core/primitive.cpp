
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

// primitive.cpp*
#include "primitive.h"
#include "light.h"
// Primitive Method Definitions
Primitive::~Primitive() { }

bool Primitive::CanIntersect() const {
	return true;
}

void
Primitive::Refine(vector<Reference<Primitive> > &refined)
const {
	Severe("Unimplemented Primitive::Refine"
	         "method called!");
}
void Primitive::FullyRefine(
		vector<Reference<Primitive> > &refined) const {
	vector<Reference<Primitive> > todo;
	todo.push_back(const_cast<Primitive *>(this));
	while (todo.size()) {
		// Refine last primitive in todo list
		Reference<Primitive> prim = todo.back();
		todo.pop_back();
		if (prim->CanIntersect())
			refined.push_back(prim);
		else
			prim->Refine(todo);
	}
}
const AreaLight *Aggregate::GetAreaLight() const {
	Severe("Aggregate::GetAreaLight() method"
	     "called; should have gone to GeometricPrimitive");
	return NULL;
}
BSDF *Aggregate::GetBSDF(const DifferentialGeometry &,
		const Transform &) const {
	Severe("Aggregate::GetBSDF() method"
	    "called; should have gone to GeometricPrimitive");
	return NULL;
}
// InstancePrimitive Method Definitions
bool InstancePrimitive::Intersect(const Ray &r,
                               Intersection *isect) const {
	Ray ray = WorldToInstance(r);
	if (!instance->Intersect(ray, isect))
		return false;
	r.maxt = ray.maxt;
	isect->WorldToObject = isect->WorldToObject *
		WorldToInstance;
	// Transform instance's differential geometry to world space
	isect->dg.p = InstanceToWorld(isect->dg.p);
	isect->dg.nn = Normalize(InstanceToWorld(isect->dg.nn));
	isect->dg.dpdu = InstanceToWorld(isect->dg.dpdu);
	isect->dg.dpdv = InstanceToWorld(isect->dg.dpdv);
	isect->dg.dndu = InstanceToWorld(isect->dg.dndu);
	isect->dg.dndv = InstanceToWorld(isect->dg.dndv);
	return true;
}
bool InstancePrimitive::IntersectP(const Ray &r) const {
	return instance->IntersectP(WorldToInstance(r));
}
// GeometricPrimitive Method Definitions
BBox GeometricPrimitive::WorldBound() const {
	return shape->WorldBound();
}
bool GeometricPrimitive::IntersectP(const Ray &r) const {
	return shape->IntersectP(r);
}
bool GeometricPrimitive::CanIntersect() const {
	return shape->CanIntersect();
}
void GeometricPrimitive::
        Refine(vector<Reference<Primitive> > &refined)
        const {
	vector<Reference<Shape> > r;
	shape->Refine(r);
	for (u_int i = 0; i < r.size(); ++i) {
		GeometricPrimitive *gp =
		    new GeometricPrimitive(r[i],
			   material, areaLight);
		refined.push_back(gp);
	}
}
GeometricPrimitive::
    GeometricPrimitive(const Reference<Shape> &s,
		const Reference<Material> &m, AreaLight *a)
	: shape(s), material(m), areaLight(a) {
}
bool GeometricPrimitive::Intersect(const Ray &r,
		Intersection *isect) const {
	float thit;
	if (!shape->Intersect(r, &thit, &isect->dg))
		return false;
	isect->primitive = this;
	isect->WorldToObject = shape->WorldToObject;
	r.maxt = thit;
	return true;
}
const AreaLight *GeometricPrimitive::GetAreaLight() const {
	return areaLight;
}
BSDF *
GeometricPrimitive::GetBSDF(const DifferentialGeometry &dg,
		const Transform &WorldToObject) const {
	DifferentialGeometry dgs;
	shape->GetShadingGeometry(WorldToObject.GetInverse(),
		dg, &dgs);
	return material->GetBSDF(dg, dgs);
}
// Intersection Method Definitions
BSDF *Intersection::GetBSDF(const RayDifferential &ray)
		const {
	static StatsCounter pointsShaded("Shading", "Number of points shaded"); // NOBOOK
	++pointsShaded; // NOBOOK
	dg.ComputeDifferentials(ray);
	return primitive->GetBSDF(dg, WorldToObject);
}
Spectrum Intersection::Le(const Vector &w) const {
	const AreaLight *area = primitive->GetAreaLight();
	return area ? area->L(dg.p, dg.nn, w) : Spectrum(0.);
}
