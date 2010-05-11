
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

#ifndef PBRT_MC_H
#define PBRT_MC_H
// mc.h*
// MC Utility Declarations
extern COREDLL void RejectionSampleDisk(float u1, float u2, float *x, float *y);
COREDLL Vector UniformSampleHemisphere(float u1, float u2);
COREDLL float  UniformHemispherePdf(float theta, float phi);
COREDLL Vector UniformSampleSphere(float u1, float u2);
COREDLL float  UniformSpherePdf();
COREDLL Vector UniformSampleCone(float u1, float u2, float costhetamax);
COREDLL Vector UniformSampleCone(float u1, float u2, float costhetamax,
	const Vector &x, const Vector &y, const Vector &z);
COREDLL float  UniformConePdf(float costhetamax);
COREDLL void UniformSampleDisk(float u1, float u2, float *x, float *y);
inline Vector CosineSampleHemisphere(float u1, float u2) {
	Vector ret;
	ConcentricSampleDisk(u1, u2, &ret.x, &ret.y);
	ret.z = sqrtf(max(0.f,
	                  1.f - ret.x*ret.x - ret.y*ret.y));
	return ret;
}
inline float CosineHemispherePdf(float costheta, float phi) {
	return costheta * INV_PI;
}
extern COREDLL Vector SampleHG(const Vector &w, float g, float u1, float u2);
extern COREDLL float HGPdf(const Vector &w, const Vector &wp, float g);
// MC Inline Functions
inline float BalanceHeuristic(int nf, float fPdf, int ng,
		float gPdf) {
	return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}
inline float PowerHeuristic(int nf, float fPdf, int ng,
		float gPdf) {
	float f = nf * fPdf, g = ng * gPdf;
	return (f*f) / (f*f + g*g);
}
#endif // PBRT_MC_H
