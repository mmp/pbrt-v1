
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

// exponential.cpp*
#include "volume.h"
// ExponentialDensity Declarations
class ExponentialDensity : public DensityRegion {
public:
	// ExponentialDensity Public Methods
	ExponentialDensity(const Spectrum &sa, const Spectrum &ss,
			float gg, const Spectrum &emit, const BBox &e,
			const Transform &v2w, float aa, float bb,
			const Vector &up)
		: DensityRegion(sa, ss, gg, emit, v2w),
		  extent(e), a(aa), b(bb) {
		upDir = Normalize(up);
	}
	BBox WorldBound() const { return WorldToVolume.GetInverse()(extent); }
	bool IntersectP(const Ray &r, float *t0, float *t1) const {
		Ray ray = WorldToVolume(r);
		return extent.IntersectP(ray, t0, t1);
	}
	float Density(const Point &Pobj) const {
		if (!extent.Inside(Pobj)) return 0;
		float height = Dot(Pobj - extent.pMin, upDir);
		return a * expf(-b * height);
	}
private:
	// ExponentialDensity Private Data
	BBox extent;
	float a, b;
	Vector upDir;
};
// ExponentialDensity Method Definitions
extern "C" DLLEXPORT VolumeRegion *CreateVolumeRegion(const Transform &volume2world,
		const ParamSet &params) {
	// Initialize common volume region parameters
	Spectrum sigma_a = params.FindOneSpectrum("sigma_a", 0.);
	Spectrum sigma_s = params.FindOneSpectrum("sigma_s", 0.);
	float g = params.FindOneFloat("g", 0.);
	Spectrum Le = params.FindOneSpectrum("Le", 0.);
	Point p0 = params.FindOnePoint("p0", Point(0,0,0));
	Point p1 = params.FindOnePoint("p1", Point(1,1,1));
	float a = params.FindOneFloat("a", 1.);
	float b = params.FindOneFloat("b", 1.);
	Vector up = params.FindOneVector("updir", Vector(0,1,0));
	return new ExponentialDensity(sigma_a, sigma_s, g, Le, BBox(p0, p1),
		volume2world, a, b, up);
}
