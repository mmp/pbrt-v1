
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

// environment.cpp*
#include "camera.h"
#include "film.h"
#include "paramset.h"
// EnvironmentCamera Declarations
class EnvironmentCamera : public Camera {
public:
	// EnvironmentCamera Public Methods
	EnvironmentCamera(const Transform &world2cam, float hither,
		float yon, float sopen, float sclose, Film *film);
	float GenerateRay(const Sample &sample, Ray *) const;
private:
	// EnvironmentCamera Private Data
	Point rayOrigin;
};
// EnvironmentCamera Method Definitions
EnvironmentCamera::
    EnvironmentCamera(const Transform &world2cam,
		float hither, float yon, float sopen, float sclose,
		Film *film)
	: Camera(world2cam, hither, yon, sopen, sclose, film) {
	rayOrigin = CameraToWorld(Point(0,0,0));
}
float EnvironmentCamera::GenerateRay(const Sample &sample,
		Ray *ray) const {
	ray->o = rayOrigin;
	// Generate environment camera ray direction
	float theta = M_PI * sample.imageY / film->yResolution;
	float phi = 2 * M_PI * sample.imageX / film->xResolution;
	Vector dir(sinf(theta) * cosf(phi), cosf(theta),
		sinf(theta) * sinf(phi));
	CameraToWorld(dir, &ray->d);
	// Set ray time value
	ray->time = Lerp(sample.time, ShutterOpen, ShutterClose);
	ray->mint = ClipHither;
	ray->maxt = ClipYon;
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
	(void) lensradius; // don't need this
	(void) focaldistance; // don't need this
	return new EnvironmentCamera(world2cam, hither, yon,
		shutteropen, shutterclose, film);
}
