
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// volume.cpp*
#include "volume.h"
// Volume Scattering Definitions
COREDLL
float PhaseIsotropic(const Vector &, const Vector &) {
	return 1.f / (4.f * M_PI);
}
COREDLL float PhaseRayleigh(const Vector &w, const Vector &wp) {
	float costheta = Dot(w, wp);
	return  3.f/(16.f*M_PI) * (1 + costheta * costheta);
}
COREDLL float PhaseMieHazy(const Vector &w, const Vector &wp) {
	float costheta = Dot(w, wp);
	return (0.5f + 4.5f * powf(0.5 * (1.f + costheta), 8.f)) / (4.f*M_PI);
}
COREDLL float PhaseMieMurky(const Vector &w, const Vector &wp) {
	float costheta = Dot(w, wp);
	return (0.5f + 16.5f * powf(0.5 * (1.f + costheta), 32.f)) / (4.f*M_PI);
}
COREDLL
float PhaseHG(const Vector &w, const Vector &wp, float g) {
	float costheta = Dot(w, wp);
	return 1.f / (4.f * M_PI) * (1.f - g*g) /
		powf(1.f + g*g - 2.f * g * costheta, 1.5f);
}
COREDLL
float PhaseSchlick(const Vector &w,
                   const Vector &wp, float g) {
	float k = 1.55f * g - .55f * g * g * g;
	float kcostheta = k * Dot(w, wp);
	return 1.f / (4.f * M_PI) * (1.f - k*k) /
		((1.f - kcostheta) * (1.f - kcostheta));
}
Spectrum VolumeRegion::sigma_t(const Point &p,
                               const Vector &w) const {
	return sigma_a(p, w) + sigma_s(p, w);
}
DensityRegion::DensityRegion(const Spectrum &sa,
		const Spectrum &ss, float gg,
	 	const Spectrum &emit,
		const Transform &VolumeToWorld)
	: sig_a(sa), sig_s(ss), le(emit), g(gg)  {
	WorldToVolume = VolumeToWorld.GetInverse();
}
AggregateVolume::
	AggregateVolume(const vector<VolumeRegion *> &r) {
	regions = r;
	for (u_int i = 0; i < regions.size(); ++i)
		bound = Union(bound, regions[i]->WorldBound());
}
Spectrum AggregateVolume::sigma_a(const Point &p,
		const Vector &w) const {
	Spectrum s(0.);
	for (u_int i = 0; i < regions.size(); ++i)
		s += regions[i]->sigma_a(p, w);
	return s;
}
Spectrum AggregateVolume::sigma_s(const Point &p, const Vector &w) const {
	Spectrum s(0.);
	for (u_int i = 0; i < regions.size(); ++i)
		s += regions[i]->sigma_s(p, w);
	return s;
}
Spectrum AggregateVolume::Lve(const Point &p, const Vector &w) const {
	Spectrum L(0.);
	for (u_int i = 0; i < regions.size(); ++i)
		L += regions[i]->Lve(p, w);
	return L;
}
float AggregateVolume::p(const Point &p, const Vector &w, const Vector &wp) const {
	float ph = 0, sumWt = 0;
	for (u_int i = 0; i < regions.size(); ++i) {
		float sigt = regions[i]->sigma_t(p, w).y();
		if (sigt != 0.f) {
			float wt = regions[i]->sigma_s(p, w).y() / sigt;
			sumWt += wt;
			ph += wt * regions[i]->p(p, w, wp);
		}
	}
	return ph / sumWt;
}
Spectrum AggregateVolume::sigma_t(const Point &p, const Vector &w) const {
	Spectrum s(0.);
	for (u_int i = 0; i < regions.size(); ++i)
		s += regions[i]->sigma_t(p, w);
	return s;
}
Spectrum AggregateVolume::Tau(const Ray &ray, float step, float offset) const {
	Spectrum t(0.);
	for (u_int i = 0; i < regions.size(); ++i)
		t += regions[i]->Tau(ray, step, offset);
	return t;
}
bool AggregateVolume::IntersectP(const Ray &ray,
		float *t0, float *t1) const {
	*t0 = INFINITY;
	*t1 = -INFINITY;
	for (u_int i = 0; i < regions.size(); ++i) {
		float tr0, tr1;
		if (regions[i]->IntersectP(ray, &tr0, &tr1)) {
			*t0 = min(*t0, tr0);
			*t1 = max(*t1, tr1);
		}
	}
	return (*t0 < *t1);
}
AggregateVolume::~AggregateVolume() {
	for (u_int i = 0; i < regions.size(); ++i)
		delete regions[i];
}
BBox AggregateVolume::WorldBound() const {
	return bound;
}
Spectrum DensityRegion::Tau(const Ray &r,
		float stepSize, float u) const {
	float t0, t1;
	float length = r.d.Length();
	if (length == 0.f) return 0.f;
	Ray rn(r.o, r.d / length,
	       r.mint * length,
		   r.maxt * length);
	if (!IntersectP(rn, &t0, &t1)) return 0.;
	Spectrum tau(0.);
	t0 += u * stepSize;
	while (t0 < t1) {
		tau += sigma_t(rn(t0), -rn.d);
		t0 += stepSize;
	}
	return tau * stepSize;
}
