
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// scene.cpp*
#include "scene.h"
#include "camera.h"
#include "film.h"
#include "sampling.h"
#include "dynload.h"
#include "volume.h"
// Scene Methods
void Scene::Render() {
	// Allocate and initialize _sample_
	Sample *sample = new Sample(surfaceIntegrator,
	                            volumeIntegrator,
	                            this);
	// Allow integrators to do pre-processing for the scene
	surfaceIntegrator->Preprocess(this);
	volumeIntegrator->Preprocess(this);
	// Trace rays: The main loop
	ProgressReporter progress(sampler->TotalSamples(), "Rendering");
	while (sampler->GetNextSample(sample)) {
		// Find camera ray for _sample_
		RayDifferential ray;
		float rayWeight = camera->GenerateRay(*sample, &ray);
		// Generate ray differentials for camera ray
		++(sample->imageX);
		float wt1 = camera->GenerateRay(*sample, &ray.rx);
		--(sample->imageX);
		++(sample->imageY);
		float wt2 = camera->GenerateRay(*sample, &ray.ry);
		if (wt1 > 0 && wt2 > 0) ray.hasDifferentials = true;
		--(sample->imageY);
		// Evaluate radiance along camera ray
		float alpha;
		Spectrum Ls = 0.f;
		if (rayWeight > 0.f)
			Ls = rayWeight * Li(ray, sample, &alpha);
		// Issue warning if unexpected radiance value returned
		if (Ls.IsNaN()) {
			Error("Not-a-number radiance value returned "
		          "for image sample.  Setting to black.");
			Ls = Spectrum(0.f);
		}
		else if (Ls.y() < -1e-5) {
			Error("Negative luminance value, %g, returned "
		          "for image sample.  Setting to black.", Ls.y());
			Ls = Spectrum(0.f);
		}
		else if (isinf(Ls.y())) {
			Error("Infinite luminance value returned "
		          "for image sample.  Setting to black.");
			Ls = Spectrum(0.f);
		}
		// Add sample contribution to image
		camera->film->AddSample(*sample, ray, Ls, alpha);
		// Free BSDF memory from computing image sample value
		BSDF::FreeAll();
		// Report rendering progress
		static StatsCounter cameraRaysTraced("Camera", "Camera Rays Traced");
		++cameraRaysTraced;
		progress.Update();
	}
	// Clean up after rendering and store final image
	delete sample;
	progress.Done();
	camera->film->WriteImage();
}
Scene::~Scene() {
	delete camera;
	delete sampler;
	delete surfaceIntegrator;
	delete volumeIntegrator;
	delete aggregate;
	delete volumeRegion;
	for (u_int i = 0; i < lights.size(); ++i)
		delete lights[i];
}
Scene::Scene(Camera *cam, SurfaceIntegrator *si,
             VolumeIntegrator *vi, Sampler *s,
             Primitive *accel, const vector<Light *> &lts,
             VolumeRegion *vr) {
	lights = lts;
	aggregate = accel;
	camera = cam;
	sampler = s;
	surfaceIntegrator = si;
	volumeIntegrator = vi;
	volumeRegion = vr;
	if (lts.size() == 0)
		Warning("No light sources defined in scene; "
			"possibly rendering a black image.");
	// Scene Constructor Implementation
	bound = aggregate->WorldBound();
	if (volumeRegion) bound = Union(bound, volumeRegion->WorldBound());
}
const BBox &Scene::WorldBound() const {
	return bound;
}
Spectrum Scene::Li(const RayDifferential &ray,
		const Sample *sample, float *alpha) const {
	Spectrum Lo = surfaceIntegrator->Li(this, ray, sample, alpha);
	Spectrum T = volumeIntegrator->Transmittance(this, ray, sample, alpha);
	Spectrum Lv = volumeIntegrator->Li(this, ray, sample, alpha);
	return T * Lo + Lv;
}
Spectrum Scene::Transmittance(const Ray &ray) const {
	return volumeIntegrator->Transmittance(this, ray, NULL, NULL);
}
