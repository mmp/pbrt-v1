
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

// transport.cpp*
#include "transport.h"
#include "scene.h"
// Integrator Method Definitions
Integrator::~Integrator() {
}
// Integrator Utility Functions
COREDLL Spectrum UniformSampleAllLights(const Scene *scene,
		const Point &p, const Normal &n, const Vector &wo,
		BSDF *bsdf, const Sample *sample,
		int *lightSampleOffset, int *bsdfSampleOffset,
		int *bsdfComponentOffset) {
	Spectrum L(0.);
	for (u_int i = 0; i < scene->lights.size(); ++i) {
		Light *light = scene->lights[i];
		int nSamples = (sample && lightSampleOffset) ?
			sample->n2D[lightSampleOffset[i]] : 1;
		// Estimate direct lighting from _light_ samples
		Spectrum Ld(0.);
		for (int j = 0; j < nSamples; ++j)
			Ld += EstimateDirect(scene, light, p, n, wo, bsdf,
				sample, lightSampleOffset[i], bsdfSampleOffset[i],
				bsdfComponentOffset[i], j);
		L += Ld / nSamples;
	}
	return L;
}
COREDLL Spectrum UniformSampleOneLight(const Scene *scene,
		const Point &p, const Normal &n,
		const Vector &wo, BSDF *bsdf, const Sample *sample,
		int lightSampleOffset, int lightNumOffset,
		int bsdfSampleOffset, int bsdfComponentOffset) {
	// Randomly choose a single light to sample, _light_
	int nLights = int(scene->lights.size());
	int lightNum;
	if (lightNumOffset != -1)
		lightNum = Floor2Int(sample->oneD[lightNumOffset][0] *
							 nLights);
	else
		lightNum = Floor2Int(RandomFloat() * nLights);
	lightNum = min(lightNum, nLights-1);
	Light *light = scene->lights[lightNum];
	return (float)nLights *
		EstimateDirect(scene, light, p, n, wo, bsdf, sample,
			lightSampleOffset, bsdfSampleOffset,
			bsdfComponentOffset, 0);
}
COREDLL Spectrum WeightedSampleOneLight(const Scene *scene,
		const Point &p, const Normal &n,
		const Vector &wo, BSDF *bsdf,
		const Sample *sample, int lightSampleOffset,
		int lightNumOffset, int bsdfSampleOffset,
		int bsdfComponentOffset, float *&avgY,
		float *&avgYsample, float *&cdf,
		float &overallAvgY) {
	int nLights = int(scene->lights.size());
	// Initialize _avgY_ array if necessary
	if (!avgY) {
		avgY = new float[nLights];
		avgYsample = new float[nLights];
		cdf = new float[nLights+1];
		for (int i = 0; i < nLights; ++i)
			avgY[i] = avgYsample[i] = 0.;
	}
	Spectrum L(0.);
	if (overallAvgY == 0.) {
		// Sample one light uniformly and initialize luminance arrays
		L = UniformSampleOneLight(scene, p, n,
		    wo, bsdf, sample, lightSampleOffset,
			lightNumOffset, bsdfSampleOffset,
			bsdfComponentOffset);
		float luminance = L.y();
		overallAvgY = luminance;
		for (int i = 0; i < nLights; ++i)
			avgY[i] = luminance;
	}
	else {
		// Choose _light_ according to average reflected luminance
		float c, lightSampleWeight;
		for (int i = 0; i < nLights; ++i)
			avgYsample[i] = max(avgY[i], .1f * overallAvgY);
		ComputeStep1dCDF(avgYsample, nLights, &c, cdf);
		float t = SampleStep1d(avgYsample, cdf, c, nLights,
			sample->oneD[lightNumOffset][0], &lightSampleWeight);
		int lightNum = min(Float2Int(nLights * t), nLights-1);
		Light *light = scene->lights[lightNum];
		L = EstimateDirect(scene, light, p, n, wo, bsdf,
			sample, lightSampleOffset, bsdfSampleOffset,
			bsdfComponentOffset, 0);
		// Update _avgY_ array with reflected radiance due to light
		float luminance = L.y();
		avgY[lightNum] =
			ExponentialAverage(avgY[lightNum], luminance, .99f);
		overallAvgY =
			ExponentialAverage(overallAvgY, luminance, .999f);
		L /= lightSampleWeight;
	}
	return L;
}
Spectrum EstimateDirect(const Scene *scene,
        const Light *light, const Point &p,
		const Normal &n, const Vector &wo,
		BSDF *bsdf, const Sample *sample, int lightSamp,
		int bsdfSamp, int bsdfComponent, u_int sampleNum) {
	Spectrum Ld(0.);
	// Find light and BSDF sample values for direct lighting estimate
	float ls1, ls2, bs1, bs2, bcs;
	if (lightSamp != -1 && bsdfSamp != -1 &&
		sampleNum < sample->n2D[lightSamp] &&
		sampleNum < sample->n2D[bsdfSamp]) {
		ls1 = sample->twoD[lightSamp][2*sampleNum];
		ls2 = sample->twoD[lightSamp][2*sampleNum+1];
		bs1 = sample->twoD[bsdfSamp][2*sampleNum];
		bs2 = sample->twoD[bsdfSamp][2*sampleNum+1];
		bcs = sample->oneD[bsdfComponent][sampleNum];
	}
	else {
		ls1 = RandomFloat();
		ls2 = RandomFloat();
		bs1 = RandomFloat();
		bs2 = RandomFloat();
		bcs = RandomFloat();
	}
	// Sample light source with multiple importance sampling
	Vector wi;
	float lightPdf, bsdfPdf;
	VisibilityTester visibility;
	Spectrum Li = light->Sample_L(p, n,
		ls1, ls2, &wi, &lightPdf, &visibility);
	if (lightPdf > 0. && !Li.Black()) {
		Spectrum f = bsdf->f(wo, wi);
		if (!f.Black() && visibility.Unoccluded(scene)) {
			// Add light's contribution to reflected radiance
			Li *= visibility.Transmittance(scene);
			if (light->IsDeltaLight())
				Ld += f * Li * AbsDot(wi, n) / lightPdf;
			else {
				bsdfPdf = bsdf->Pdf(wo, wi);
				float weight = PowerHeuristic(1, lightPdf, 1, bsdfPdf);
				Ld += f * Li * AbsDot(wi, n) * weight / lightPdf;
			}
		}
	}
	// Sample BSDF with multiple importance sampling
	if (!light->IsDeltaLight()) {
		BxDFType flags = BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
		Spectrum f = bsdf->Sample_f(wo, &wi,
			bs1, bs2, bcs, &bsdfPdf, flags);
		if (!f.Black() && bsdfPdf > 0.) {
			lightPdf = light->Pdf(p, n, wi);
			if (lightPdf > 0.) {
				// Add light contribution from BSDF sampling
				float weight = PowerHeuristic(1, bsdfPdf, 1, lightPdf);
				Intersection lightIsect;
				Spectrum Li(0.f);
				RayDifferential ray(p, wi);
				if (scene->Intersect(ray, &lightIsect)) {
					if (lightIsect.primitive->GetAreaLight() == light)
						Li = lightIsect.Le(-wi);
				}
				else
					Li = light->Le(ray);
				if (!Li.Black()) {
					Li *= scene->Transmittance(ray);
					Ld += f * Li * AbsDot(wi, n) * weight / bsdfPdf;
				}
			}
		}
	}
	return Ld;
}
