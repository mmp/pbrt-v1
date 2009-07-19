
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_CAMERA_H
#define PBRT_CAMERA_H
// camera.h*
#include "pbrt.h"
#include "color.h"
#include "sampling.h"
#include "geometry.h"
#include "transform.h"
// Camera Declarations
class COREDLL Camera {
public:
	// Camera Interface
	virtual float GenerateRay(const Sample &sample,
		                      Ray *ray) const = 0;
	virtual ~Camera();
	Camera(const Transform &world2cam, float hither, float yon,
		float sopen, float sclose, Film *film);
	// Camera Public Data
	Film *film;
protected:
	// Camera Protected Data
	Transform WorldToCamera, CameraToWorld;
	float ClipHither, ClipYon;
	float ShutterOpen, ShutterClose;
};
class COREDLL ProjectiveCamera : public Camera {
public:
	// ProjectiveCamera Public Methods
	ProjectiveCamera(const Transform &world2cam,
	    const Transform &proj, const float Screen[4],
		float hither, float yon,
		float sopen, float sclose,
		float lensr, float focald, Film *film);
protected:
	// ProjectiveCamera Protected Data
	Transform CameraToScreen, WorldToScreen, RasterToCamera;
	Transform ScreenToRaster, RasterToScreen;
	float LensRadius, FocalDistance;
};
#endif // PBRT_CAMERA_H
