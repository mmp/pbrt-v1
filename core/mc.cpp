
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

// mc.cpp*
#include "pbrt.h"
#include "geometry.h"
#include "shape.h"
#include "mc.h"
#include "volume.h"
// MC Function Definitions
void ComputeStep1dCDF(float *f, int nSteps, float *c,
		float *cdf) {
	// Compute integral of step function at $x_i$
	int i;
	cdf[0] = 0.;
	for (i = 1; i < nSteps+1; ++i)
		cdf[i] = cdf[i-1] + f[i-1] / nSteps;
	// Transform step function integral into cdf
	*c = cdf[nSteps];
	for (i = 1; i < nSteps+1; ++i)
		cdf[i] /= *c;
}
float SampleStep1d(float *f, float *cdf, float c,
		int nSteps, float u, float *pdf) {
	// Find surrounding cdf segments
	float *ptr = std::upper_bound(cdf, cdf+nSteps+1, u);
	int offset = max(0, (int) (ptr-cdf-1));
	Assert(u >= cdf[offset] && u < cdf[offset+1]);
	// Return offset along current cdf segment
	u = (u - cdf[offset]) / (cdf[offset+1] - cdf[offset]);
	*pdf = f[offset] / c;
	return (offset + u) / nSteps;
}
void RejectionSampleDisk(float *x, float *y) {
	float sx, sy;
	do {
		sx = 1.f - 2.f * RandomFloat();
		sy = 1.f - 2.f * RandomFloat();
	} while (sx*sx + sy*sy > 1.f);
	*x = sx;
	*y = sy;
}
COREDLL Vector UniformSampleHemisphere(float u1, float u2) {
	float z = u1;
	float r = sqrtf(max(0.f, 1.f - z*z));
	float phi = 2 * M_PI * u2;
	float x = r * cosf(phi);
	float y = r * sinf(phi);
	return Vector(x, y, z);
}
COREDLL float UniformHemispherePdf(float theta, float phi) {
	return INV_TWOPI;
}
COREDLL Vector UniformSampleSphere(float u1, float u2) {
	float z = 1.f - 2.f * u1;
	float r = sqrtf(max(0.f, 1.f - z*z));
	float phi = 2.f * M_PI * u2;
	float x = r * cosf(phi);
	float y = r * sinf(phi);
	return Vector(x, y, z);
}
COREDLL float UniformSpherePdf() {
	return 1.f / (4.f * M_PI);
}
COREDLL void UniformSampleDisk(float u1, float u2,
		float *x, float *y) {
	float r = sqrtf(u1);
	float theta = 2.0f * M_PI * u2;
	*x = r * cosf(theta);
	*y = r * sinf(theta);
}
COREDLL void ConcentricSampleDisk(float u1, float u2,
		float *dx, float *dy) {
	float r, theta;
	// Map uniform random numbers to $[-1,1]^2$
	float sx = 2 * u1 - 1;
	float sy = 2 * u2 - 1;
	// Map square to $(r,\theta)$
	// Handle degeneracy at the origin
	if (sx == 0.0 && sy == 0.0) {
		*dx = 0.0;
		*dy = 0.0;
		return;
	}
	if (sx >= -sy) {
		if (sx > sy) {
			// Handle first region of disk
			r = sx;
			if (sy > 0.0)
				theta = sy/r;
			else
				theta = 8.0f + sy/r;
		}
		else {
			// Handle second region of disk
			r = sy;
			theta = 2.0f - sx/r;
		}
	}
	else {
		if (sx <= sy) {
			// Handle third region of disk
			r = -sx;
			theta = 4.0f - sy/r;
		}
		else {
			// Handle fourth region of disk
			r = -sy;
			theta = 6.0f + sx/r;
		}
	}
	theta *= M_PI / 4.f;
	*dx = r*cosf(theta);
	*dy = r*sinf(theta);
}
COREDLL void UniformSampleTriangle(float u1, float u2,
		float *u, float *v) {
	float su1 = sqrtf(u1);
	*u = 1.f - su1;
	*v = u2 * su1;
}
COREDLL float UniformConePdf(float cosThetaMax) {
	return 1.f / (2.f * M_PI * (1.f - cosThetaMax));
}
Vector UniformSampleCone(float u1, float u2,
		float costhetamax) {
	float costheta = Lerp(u1, costhetamax, 1.f);
	float sintheta = sqrtf(1.f - costheta*costheta);
	float phi = u2 * 2.f * M_PI;
	return Vector(cosf(phi) * sintheta,
	              sinf(phi) * sintheta,
		          costheta);
}
COREDLL Vector UniformSampleCone(float u1, float u2, float costhetamax,
		const Vector &x, const Vector &y, const Vector &z) {
	float costheta = Lerp(u1, costhetamax, 1.f);
	float sintheta = sqrtf(1.f - costheta*costheta);
	float phi = u2 * 2.f * M_PI;
	return cosf(phi) * sintheta * x + sinf(phi) * sintheta * y +
		costheta * z;
}
COREDLL Vector SampleHG(const Vector &w, float g,
		float u1, float u2) {
	float costheta;
	if (fabsf(g) < 1e-3)
		costheta = 1.f - 2.f * u1;
	else {
		float sqrTerm = (1.f - g * g) /
				(1.f - g + 2.f * g * u1);
		costheta = (1.f + g * g - sqrTerm * sqrTerm) / (2.f * g);
	}
	float sintheta = sqrtf(max(0.f, 1.f-costheta*costheta));
	float phi = 2.f * M_PI * u2;
	Vector v1, v2;
	CoordinateSystem(w, &v1, &v2);
	return SphericalDirection(sintheta, costheta,
		phi, v1, v2, w);
}
COREDLL float HGPdf(const Vector &w, const Vector &wp,
		float g) {
	return PhaseHG(w, wp, g);
}
