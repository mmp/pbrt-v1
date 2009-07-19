
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_TRANSFORM_H
#define PBRT_TRANSFORM_H
// transform.h*
#include "pbrt.h"
#include "geometry.h"
// Transform Declarations
class COREDLL Transform {
public:
	// Transform Public Methods
	Transform() {
		m = mInv = new Matrix4x4;
	}
	Transform(float mat[4][4]) {
		m=new Matrix4x4(mat[0][0],mat[0][1],mat[0][2],mat[0][3],
	                	mat[1][0],mat[1][1],mat[1][2],mat[1][3],
	                	mat[2][0],mat[2][1],mat[2][2],mat[2][3],
	                	mat[3][0],mat[3][1],mat[3][2],mat[3][3]);
		mInv = m->Inverse();
	}
	Transform(const Reference<Matrix4x4> &mat) {
		m = mat;
		mInv = m->Inverse();
	}
	Transform(const Reference<Matrix4x4> &mat,
	          const Reference<Matrix4x4> &minv) {
		m = mat;
		mInv = minv;
	}
	friend ostream &operator<<(ostream &, const Transform &);
	Transform GetInverse() const {
		return Transform(mInv, m);
	}
	bool HasScale() const;
	inline Point operator()(const Point &pt) const;
	inline void operator()(const Point &pt,Point *ptrans) const;
	inline Vector operator()(const Vector &v) const;
	inline void operator()(const Vector &v, Vector *vt) const;
	inline Normal operator()(const Normal &) const;
	inline void operator()(const Normal &, Normal *nt) const;
	inline Ray operator()(const Ray &r) const;
	inline void operator()(const Ray &r, Ray *rt) const;
	BBox operator()(const BBox &b) const;
	Transform operator*(const Transform &t2) const;
	bool SwapsHandedness() const;
private:
	// Transform Private Data
	Reference<Matrix4x4> m, mInv;
};
// Transform Inline Functions
inline Point Transform::operator()(const Point &pt) const {
	float x = pt.x, y = pt.y, z = pt.z;
	float xp = m->m[0][0]*x + m->m[0][1]*y + m->m[0][2]*z + m->m[0][3];
	float yp = m->m[1][0]*x + m->m[1][1]*y + m->m[1][2]*z + m->m[1][3];
	float zp = m->m[2][0]*x + m->m[2][1]*y + m->m[2][2]*z + m->m[2][3];
	float wp = m->m[3][0]*x + m->m[3][1]*y + m->m[3][2]*z + m->m[3][3];

	Assert(wp != 0);
	if (wp == 1.) return Point(xp, yp, zp);
	else          return Point(xp, yp, zp)/wp;
}
inline void Transform::operator()(const Point &pt,
		Point *ptrans) const {
	float x = pt.x, y = pt.y, z = pt.z;
	ptrans->x = m->m[0][0]*x + m->m[0][1]*y + m->m[0][2]*z + m->m[0][3];
	ptrans->y = m->m[1][0]*x + m->m[1][1]*y + m->m[1][2]*z + m->m[1][3];
	ptrans->z = m->m[2][0]*x + m->m[2][1]*y + m->m[2][2]*z + m->m[2][3];
	float w   = m->m[3][0]*x + m->m[3][1]*y + m->m[3][2]*z + m->m[3][3];
	if (w != 1.) *ptrans /= w;
}
inline Vector Transform::operator()(const Vector &v) const {
  float x = v.x, y = v.y, z = v.z;
  return Vector(m->m[0][0]*x + m->m[0][1]*y + m->m[0][2]*z,
			    m->m[1][0]*x + m->m[1][1]*y + m->m[1][2]*z,
			    m->m[2][0]*x + m->m[2][1]*y + m->m[2][2]*z);
}
inline void Transform::operator()(const Vector &v,
		Vector *vt) const {
  float x = v.x, y = v.y, z = v.z;
  vt->x = m->m[0][0] * x + m->m[0][1] * y + m->m[0][2] * z;
  vt->y = m->m[1][0] * x + m->m[1][1] * y + m->m[1][2] * z;
  vt->z = m->m[2][0] * x + m->m[2][1] * y + m->m[2][2] * z;
}
inline Normal Transform::operator()(const Normal &n) const {
	float x = n.x, y = n.y, z = n.z;
	return Normal(mInv->m[0][0]*x + mInv->m[1][0]*y + mInv->m[2][0]*z,
                  mInv->m[0][1]*x + mInv->m[1][1]*y + mInv->m[2][1]*z,
                  mInv->m[0][2]*x + mInv->m[1][2]*y + mInv->m[2][2]*z);
}
inline void Transform::operator()(const Normal &n,
		Normal *nt) const {
	float x = n.x, y = n.y, z = n.z;
	nt->x = mInv->m[0][0]*x + mInv->m[1][0]*y + mInv->m[2][0]*z;
	nt->y = mInv->m[0][1]*x + mInv->m[1][1]*y + mInv->m[2][1]*z;
	nt->z = mInv->m[0][2]*x + mInv->m[1][2]*y + mInv->m[2][2]*z;
}
inline Ray Transform::operator()(const Ray &r) const {
	Ray ret;
	(*this)(r.o, &ret.o);
	(*this)(r.d, &ret.d);
	ret.mint = r.mint;
	ret.maxt = r.maxt;
	ret.time = r.time;
	return ret;
}
inline void Transform::operator()(const Ray &r,
                                  Ray *rt) const {
	(*this)(r.o, &rt->o);
	(*this)(r.d, &rt->d);
	rt->mint = r.mint;
	rt->maxt = r.maxt;
	rt->time = r.time;
}
#endif // PBRT_TRANSFORM_H
