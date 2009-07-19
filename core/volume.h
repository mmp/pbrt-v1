
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_VOLUME_H
#define PBRT_VOLUME_H
// volume.h*
#include "pbrt.h"
#include "color.h"
#include "geometry.h"
#include "paramset.h"
#include "transform.h"
#include "transport.h"
// Volume Scattering Declarations
COREDLL float PhaseIsotropic(const Vector &w, const Vector &wp);
COREDLL
float PhaseRayleigh(const Vector &w, const Vector &wp);
COREDLL
float PhaseMieHazy(const Vector &w, const Vector &wp);
COREDLL
float PhaseMieMurky(const Vector &w, const Vector &wp);
COREDLL float PhaseHG(const Vector &w, const Vector &wp, float g);
COREDLL float PhaseSchlick(const Vector &w, const Vector &wp, float g);
class COREDLL VolumeRegion {
public:
	// VolumeRegion Interface
	virtual ~VolumeRegion() { }
	virtual BBox WorldBound() const = 0;
	virtual bool IntersectP(const Ray &ray, float *t0,
		float *t1) const = 0;
	virtual Spectrum sigma_a(const Point &,
		const Vector &) const = 0;
	virtual Spectrum sigma_s(const Point &,
		const Vector &) const = 0;
	virtual
		Spectrum Lve(const Point &, const Vector &) const = 0;
	virtual float p(const Point &, const Vector &,
		const Vector &) const = 0;
	virtual Spectrum sigma_t(const Point &, const Vector &) const;
	virtual Spectrum Tau(const Ray &ray,
		float step = 1.f, float offset = 0.5) const = 0;
};
class COREDLL DensityRegion : public VolumeRegion {
public:
	// DensityRegion Public Methods
	DensityRegion(const Spectrum &sig_a, const Spectrum &sig_s,
		float g, const Spectrum &Le, const Transform &VolumeToWorld);
	virtual float Density(const Point &Pobj) const = 0;
	Spectrum sigma_a(const Point &p, const Vector &) const {
		return Density(WorldToVolume(p)) * sig_a;
	}
	Spectrum sigma_s(const Point &p, const Vector &) const {
		return Density(WorldToVolume(p)) * sig_s;
	}
	Spectrum sigma_t(const Point &p, const Vector &) const {
		return Density(WorldToVolume(p)) * (sig_a + sig_s);
	}
	Spectrum Lve(const Point &p, const Vector &) const {
		return Density(WorldToVolume(p)) * le;
	}
	float p(const Point &p, const Vector &w,
			const Vector &wp) const {
		return PhaseHG(w, wp, g);
	}
	Spectrum Tau(const Ray &r, float stepSize, float offset) const;
protected:
	// DensityRegion Protected Data
	Transform WorldToVolume;
	Spectrum sig_a, sig_s, le;
	float g;
};
class COREDLL AggregateVolume : public VolumeRegion {
public:
	// AggregateVolume Public Methods
	AggregateVolume(const vector<VolumeRegion *> &r);
	~AggregateVolume();
	BBox WorldBound() const;
	bool IntersectP(const Ray &ray, float *t0, float *t1) const;
	Spectrum sigma_a(const Point &, const Vector &) const;
	Spectrum sigma_s(const Point &, const Vector &) const;
	Spectrum Lve(const Point &, const Vector &) const;
	float p(const Point &, const Vector &, const Vector &) const;
	Spectrum sigma_t(const Point &, const Vector &) const;
	Spectrum Tau(const Ray &ray, float, float) const;
private:
	// AggregateVolume Private Data
	vector<VolumeRegion *> regions;
	BBox bound;
};
class VolumeIntegrator : public Integrator {
public:
	virtual Spectrum Transmittance(const Scene *scene,
		const Ray &ray, const Sample *sample,
		float *alpha) const = 0;
};
#endif // PBRT_VOLUME_H
