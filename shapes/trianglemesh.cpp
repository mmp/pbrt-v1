
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

// trianglemesh.cpp*
#include "shape.h"
#include "paramset.h"
// TriangleMesh Declarations
class TriangleMesh : public Shape {
public:
	// TriangleMesh Public Methods
	TriangleMesh(const Transform &o2w, bool ro,
	             int ntris, int nverts, const int *vptr,
				 const Point *P, const Normal *N,
				 const Vector *S, const float *uv);
	~TriangleMesh();
	BBox ObjectBound() const;
	BBox WorldBound() const;
	bool CanIntersect() const { return false; }
	void Refine(vector<Reference<Shape> > &refined) const;
	friend class Triangle;
	template <class T> friend class VertexTexture;
protected:
	// TriangleMesh Data
	int ntris, nverts;
	int *vertexIndex;
	Point *p;
	Normal *n;
	Vector *s;
	float *uvs;
};
class Triangle : public Shape {
public:
	// Triangle Public Methods
	Triangle(const Transform &o2w, bool ro,
	         TriangleMesh *m, int n)
			: Shape(o2w, ro) {
		mesh = m;
		v = &mesh->vertexIndex[3*n];
		// Update created triangles stats
		static StatsCounter trisMade("Geometry",
		                             "Triangles created");
		++trisMade;
	}
	BBox ObjectBound() const;
	BBox WorldBound() const;
	bool Intersect(const Ray &ray, float *tHit,
	               DifferentialGeometry *dg) const;
	bool IntersectP(const Ray &ray) const;
	void GetUVs(float uv[3][2]) const;
	float Area() const;
	virtual void GetShadingGeometry(const Transform &obj2world,
			const DifferentialGeometry &dg,
			DifferentialGeometry *dgShading) const {
		if (!mesh->n && !mesh->s) {
			*dgShading = dg;
			return;
		}
		// Initialize _Triangle_ shading geometry with _n_ and _s_
		// Compute barycentric coordinates for point
		float b[3];
		// Initialize _A_ and _C_ matrices for barycentrics
		float uv[3][2];
		GetUVs(uv);
		float A[2][2] =
		    { { uv[1][0] - uv[0][0], uv[2][0] - uv[0][0] },
		      { uv[1][1] - uv[0][1], uv[2][1] - uv[0][1] } };
		float C[2] = { dg.u - uv[0][0], dg.v - uv[0][1] };
		if (!SolveLinearSystem2x2(A, C, &b[1])) {
			// Handle degenerate parametric mapping
			b[0] = b[1] = b[2] = 1.f/3.f;
		}
		else
			b[0] = 1.f - b[1] - b[2];
		// Use _n_ and _s_ to compute shading tangents for triangle, _ss_ and _ts_
		Normal ns;
		Vector ss, ts;
		if (mesh->n) ns = Normalize(obj2world(b[0] * mesh->n[v[0]] +
			b[1] * mesh->n[v[1]] + b[2] * mesh->n[v[2]]));
		else   ns = dg.nn;
		if (mesh->s) ss = Normalize(obj2world(b[0] * mesh->s[v[0]] +
			b[1] * mesh->s[v[1]] + b[2] * mesh->s[v[2]]));
		else   ss = Normalize(dg.dpdu);
		ts = Normalize(Cross(ss, ns));
		ss = Cross(ts, ns);
		Vector dndu, dndv;
		if (mesh->n) {
			// Compute \dndu and \dndv for triangle shading geometry
			float uvs[3][2];
			GetUVs(uvs);
			// Compute deltas for triangle partial derivatives of normal
			float du1 = uvs[0][0] - uvs[2][0];
			float du2 = uvs[1][0] - uvs[2][0];
			float dv1 = uvs[0][1] - uvs[2][1];
			float dv2 = uvs[1][1] - uvs[2][1];
			Vector dn1 = Vector(mesh->n[v[0]] - mesh->n[v[2]]);
			Vector dn2 = Vector(mesh->n[v[1]] - mesh->n[v[2]]);
			float determinant = du1 * dv2 - dv1 * du2;
			if (determinant == 0)
				dndu = dndv = Vector(0,0,0);
			else {
				float invdet = 1.f / determinant;
				dndu = ( dv2 * dn1 - dv1 * dn2) * invdet;
				dndv = (-du2 * dn1 + du1 * dn2) * invdet;
			}
		}
		else
			dndu = dndv = Vector(0,0,0);
		*dgShading = DifferentialGeometry(dg.p, ss, ts,
			ObjectToWorld(dndu), ObjectToWorld(dndv), dg.u, dg.v, dg.shape);
		dgShading->dudx = dg.dudx;  dgShading->dvdx = dg.dvdx; // NOBOOK
		dgShading->dudy = dg.dudy;  dgShading->dvdy = dg.dvdy; // NOBOOK
		dgShading->dpdx = dg.dpdx;  dgShading->dpdy = dg.dpdy; // NOBOOK
	}
	Point Sample(float u1, float u2, Normal *Ns) const;
private:
	// Triangle Data
	Reference<TriangleMesh> mesh;
	int *v;
};
// TriangleMesh Method Definitions
TriangleMesh::TriangleMesh(const Transform &o2w, bool ro,
		int nt, int nv, const int *vi, const Point *P,
		const Normal *N, const Vector *S, const float *uv)
	: Shape(o2w, ro) {
	ntris = nt;
	nverts = nv;
	vertexIndex = new int[3 * ntris];
	memcpy(vertexIndex, vi, 3 * ntris * sizeof(int));
	// Copy _uv_, _N_, and _S_ vertex data, if present
	if (uv) {
		uvs = new float[2*nverts];
		memcpy(uvs, uv, 2*nverts*sizeof(float));
	}
	else uvs = NULL;
	p = new Point[nverts];
	if (N) {
		n = new Normal[nverts];
		memcpy(n, N, nverts*sizeof(Normal));
	}
	else n = NULL;
	if (S) {
		s = new Vector[nverts];
		memcpy(s, S, nverts*sizeof(Vector));
	}
	else s = NULL;
	// Transform mesh vertices to world space
	for (int i  = 0; i < nverts; ++i)
		p[i] = ObjectToWorld(P[i]);
}
TriangleMesh::~TriangleMesh() {
	delete[] vertexIndex;
	delete[] p;
	delete[] s;
	delete[] n;
	delete[] uvs;
}
BBox TriangleMesh::ObjectBound() const {
	BBox bobj;
	for (int i = 0; i < nverts; i++)
		bobj = Union(bobj, WorldToObject(p[i]));
	return bobj;
}
BBox TriangleMesh::WorldBound() const {
	BBox worldBounds;
	for (int i = 0; i < nverts; i++)
		worldBounds = Union(worldBounds, p[i]);
	return worldBounds;
}
void
TriangleMesh::Refine(vector<Reference<Shape> > &refined)
const {
	for (int i = 0; i < ntris; ++i)
		refined.push_back(new Triangle(ObjectToWorld,
		                               reverseOrientation,
                                       (TriangleMesh *)this,
									   i));
}
BBox Triangle::ObjectBound() const {
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];
	return Union(BBox(WorldToObject(p1), WorldToObject(p2)),
		WorldToObject(p3));
}
BBox Triangle::WorldBound() const {
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];
	return Union(BBox(p1, p2), p3);
}
bool Triangle::Intersect(const Ray &ray, float *tHit,
		DifferentialGeometry *dg) const {
	// Initialize triangle intersection statistics
	static
	StatsPercentage triangleHits("Geometry",
	                             "Triangle Ray Intersections");
	// Update triangle tests count
	triangleHits.Add(0, 1);
	// Compute $\VEC{s}_1$
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];
	Vector e1 = p2 - p1;
	Vector e2 = p3 - p1;
	Vector s1 = Cross(ray.d, e2);
	float divisor = Dot(s1, e1);
	if (divisor == 0.)
		return false;
	float invDivisor = 1.f / divisor;
	// Compute first barycentric coordinate
	Vector d = ray.o - p1;
	float b1 = Dot(d, s1) * invDivisor;
	if (b1 < 0. || b1 > 1.)
		return false;
	// Compute second barycentric coordinate
	Vector s2 = Cross(d, e1);
	float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b1 + b2 > 1.)
		return false;
	// Compute _t_ to intersection point
	float t = Dot(e2, s2) * invDivisor;
	if (t < ray.mint || t > ray.maxt)
		return false;
	triangleHits.Add(1, 0); //NOBOOK
	// Fill in _DifferentialGeometry_ from triangle hit
	// Compute triangle partial derivatives
	Vector dpdu, dpdv;
	float uvs[3][2];
	GetUVs(uvs);
	// Compute deltas for triangle partial derivatives
	float du1 = uvs[0][0] - uvs[2][0];
	float du2 = uvs[1][0] - uvs[2][0];
	float dv1 = uvs[0][1] - uvs[2][1];
	float dv2 = uvs[1][1] - uvs[2][1];
	Vector dp1 = p1 - p3, dp2 = p2 - p3;
	float determinant = du1 * dv2 - dv1 * du2;
	if (determinant == 0.f) {
		// Handle zero determinant for triangle partial derivative matrix
		CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
	}
	else {
		float invdet = 1.f / determinant;
		dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
		dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
	}
	// Interpolate $(u,v)$ triangle parametric coordinates
	float b0 = 1 - b1 - b2;
	float tu = b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0];
	float tv = b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1];
	*dg = DifferentialGeometry(ray(t), dpdu, dpdv,
	                           Vector(0,0,0), Vector(0,0,0),
							   tu, tv, this);
	*tHit = t;
	return true;
}
bool Triangle::IntersectP(const Ray &ray) const {
	// Initialize triangle intersection statistics
	static
	StatsPercentage triangleHits("Geometry",
	                             "Triangle Ray Intersections");
	// Update triangle tests count
	triangleHits.Add(0, 1);
	// Compute $\VEC{s}_1$
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];
	Vector e1 = p2 - p1;
	Vector e2 = p3 - p1;
	Vector s1 = Cross(ray.d, e2);
	float divisor = Dot(s1, e1);
	if (divisor == 0.)
		return false;
	float invDivisor = 1.f / divisor;
	// Compute first barycentric coordinate
	Vector d = ray.o - p1;
	float b1 = Dot(d, s1) * invDivisor;
	if (b1 < 0. || b1 > 1.)
		return false;
	// Compute second barycentric coordinate
	Vector s2 = Cross(d, e1);
	float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b1 + b2 > 1.)
		return false;
	// Compute _t_ to intersection point
	float t = Dot(e2, s2) * invDivisor;
	if (t < ray.mint || t > ray.maxt)
		return false;
	triangleHits.Add(1, 0); //NOBOOK
	return true;
}
void Triangle::GetUVs(float uv[3][2]) const {
	if (mesh->uvs) {
		uv[0][0] = mesh->uvs[2*v[0]];
		uv[0][1] = mesh->uvs[2*v[0]+1];
		uv[1][0] = mesh->uvs[2*v[1]];
		uv[1][1] = mesh->uvs[2*v[1]+1];
		uv[2][0] = mesh->uvs[2*v[2]];
		uv[2][1] = mesh->uvs[2*v[2]+1];
	} else {
		uv[0][0] = 0.; uv[0][1] = 0.;
		uv[1][0] = 1.; uv[1][1] = 0.;
		uv[2][0] = 1.; uv[2][1] = 1.;
	}
}
float Triangle::Area() const {
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];
	return 0.5f * Cross(p2-p1, p3-p1).Length();
}
Point Triangle::Sample(float u1, float u2,
		Normal *Ns) const {
	float b1, b2;
	UniformSampleTriangle(u1, u2, &b1, &b2);
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];
	Point p = b1 * p1 + b2 * p2 + (1.f - b1 - b2) * p3;
	Normal n = Normal(Cross(p2-p1, p3-p1));
	*Ns = Normalize(n);
	if (reverseOrientation) *Ns *= -1.f;
	return p;
}
extern "C" DLLEXPORT Shape *CreateShape(const Transform &o2w,
		bool reverseOrientation, const ParamSet &params) {
	int nvi, npi, nuvi, nsi, nni;
	const int *vi = params.FindInt("indices", &nvi);
	const Point *P = params.FindPoint("P", &npi);
	const float *uvs = params.FindFloat("uv", &nuvi);
	if (!uvs) uvs = params.FindFloat("st", &nuvi);
	// XXX should complain if uvs aren't an array of 2...
	if (uvs) {
		if (nuvi < 2 * npi) {
			Error("Not enough of \"uv\"s for triangle mesh.  Expencted %d, "
			      "found %d.  Discarding.\n", 2*npi, nuvi);
			uvs = NULL;
		}
		else if (nuvi > 2 * npi)
			Warning("More \"uv\"s provided than will be used for triangle "
			        "mesh.  (%d expcted, %d found)\n", 2*npi, nuvi);
	}
	if (!vi || !P) return NULL;
	const Vector *S = params.FindVector("S", &nsi);
	if (S && nsi != npi) {
		Error("Number of \"S\"s for triangle mesh must match \"P\"s");
		S = NULL;
	}
	const Normal *N = params.FindNormal("N", &nni);
	if (N && nni != npi) {
		Error("Number of \"N\"s for triangle mesh must match \"P\"s");
		N = NULL;
	}
	if (uvs && N) {
		// if there are normals, check for bad uv's that
		// give degenerate mappings; discard them if so
		const int *vp = vi;
		for (int i = 0; i < nvi; i += 3, vp += 3) {
			float area = .5f * Cross(P[vp[0]]-P[vp[1]], P[vp[2]]-P[vp[1]]).Length();
			if (area < 1e-7) continue; // ignore degenerate tris.
			if ((uvs[2*vp[0]] == uvs[2*vp[1]] &&
				uvs[2*vp[0]+1] == uvs[2*vp[1]+1]) ||
				(uvs[2*vp[1]] == uvs[2*vp[2]] &&
				uvs[2*vp[1]+1] == uvs[2*vp[2]+1]) ||
				(uvs[2*vp[2]] == uvs[2*vp[0]] &&
				uvs[2*vp[2]+1] == uvs[2*vp[0]+1])) {
				Warning("Degenerate uv coordinates in triangle mesh.  Discarding all uvs.");
				uvs = NULL;
				break;
			}
		}
	}
	for (int i = 0; i < nvi; ++i)
		if (vi[i] >= npi) {
			Error("trianglemesh has out of-bounds vertex index %d (%d \"P\" values were given",
				vi[i], npi);
			return NULL;
		}
	return new TriangleMesh(o2w, reverseOrientation, nvi/3, npi, vi, P,
		N, S, uvs);
}
