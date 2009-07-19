
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// orthographic.cpp*
#include "camera.h"
#include "film.h"
#include "paramset.h"
// OrthographicCamera Declarations
class OrthoCamera : public ProjectiveCamera {
public:
	// OrthoCamera Public Methods
	OrthoCamera(const Transform &world2cam,
	            const float Screen[4],
		        float hither, float yon,
				float sopen, float sclose,
				float lensr, float focald, Film *film);
	float GenerateRay(const Sample &sample, Ray *) const;
};
// OrthographicCamera Definitions
OrthoCamera::OrthoCamera(const Transform &world2cam,
		const float Screen[4], float hither, float yon,
		float sopen, float sclose, float lensr,
		float focald, Film *f)
	: ProjectiveCamera(world2cam, Orthographic(hither, yon),
		 Screen, hither, yon, sopen, sclose,
		 lensr, focald, f) {
}
float OrthoCamera::GenerateRay(const Sample &sample,
                               Ray *ray) const {
	// Generate raster and camera samples
	Point Pras(sample.imageX, sample.imageY, 0);
	Point Pcamera;
	RasterToCamera(Pras, &Pcamera);
	ray->o = Pcamera;
	ray->d = Vector(0,0,1);
	// Set ray time value
	ray->time = Lerp(sample.time, ShutterOpen, ShutterClose);
	// Modify ray for depth of field
	if (LensRadius > 0.) {
		// Sample point on lens
		float lensU, lensV;
		ConcentricSampleDisk(sample.lensU, sample.lensV,
		                     &lensU, &lensV);
		lensU *= LensRadius;
		lensV *= LensRadius;
		// Compute point on plane of focus
		float ft = (FocalDistance - ClipHither) / ray->d.z;
		Point Pfocus = (*ray)(ft);
		// Update ray for effect of lens
		ray->o.x += lensU;
		ray->o.y += lensV;
		ray->d = Pfocus - ray->o;
	}
	ray->mint = 0.;
	ray->maxt = ClipYon - ClipHither;
	ray->d = Normalize(ray->d);
	CameraToWorld(*ray, ray);
	return 1.f;
}
extern "C" DLLEXPORT Camera *CreateCamera(const ParamSet &params,
		const Transform &world2cam, Film *film) {
	// Extract common camera parameters from _ParamSet_
	float hither = max(1e-4f, params.FindOneFloat("hither", 1e-3f));
	float yon = min(params.FindOneFloat("yon", 1e30f), 1e30f);
	float shutteropen = params.FindOneFloat("shutteropen", 0.f);
	float shutterclose = params.FindOneFloat("shutterclose", 1.f);
	float lensradius = params.FindOneFloat("lensradius", 0.f);
	float focaldistance = params.FindOneFloat("focaldistance", 1e30f);
	float frame = params.FindOneFloat("frameaspectratio",
		float(film->xResolution)/float(film->yResolution));
	float screen[4];
	if (frame > 1.f) {
		screen[0] = -frame;
		screen[1] =  frame;
		screen[2] = -1.f;
		screen[3] =  1.f;
	}
	else {
		screen[0] = -1.f;
		screen[1] =  1.f;
		screen[2] = -1.f / frame;
		screen[3] =  1.f / frame;
	}
	int swi;
	const float *sw = params.FindFloat("screenwindow", &swi);
	if (sw && swi == 4)
		memcpy(screen, sw, 4*sizeof(float));
	return new OrthoCamera(world2cam, screen, hither, yon,
		shutteropen, shutterclose, lensradius, focaldistance,
		film);
}
