
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// camera.cpp*
#include "pbrt.h"
#include "camera.h"
#include "film.h"
#include "mc.h"
// Camera Method Definitions
Camera::~Camera() {
	delete film;
}
Camera::Camera(const Transform &world2cam,
               float hither, float yon,
		       float sopen, float sclose, Film *f) {
	WorldToCamera = world2cam;
	CameraToWorld = WorldToCamera.GetInverse();
	ClipHither = hither;
	ClipYon = yon;
	ShutterOpen = sopen;
	ShutterClose = sclose;
	film = f;
	if (WorldToCamera.HasScale())
		Warning("Scaling detected in world-to-camera transformation!\n"
			"The system has numerous assumptions, implicit and explicit,\n"
			"that this transform will have no scale factors in it.\n"
			"Proceed at your own risk; your image may have errors or\n"
			"the system may crash as a result of this.");
}
ProjectiveCamera::ProjectiveCamera(const Transform &w2c,
		const Transform &proj, const float Screen[4],
		float hither, float yon, float sopen,
		float sclose, float lensr, float focald, Film *f)
	: Camera(w2c, hither, yon, sopen, sclose, f) {
	// Initialize depth of field parameters
	LensRadius = lensr;
	FocalDistance = focald;
	// Compute projective camera transformations
	CameraToScreen = proj;
	WorldToScreen = CameraToScreen * WorldToCamera;
	// Compute projective camera screen transformations
	ScreenToRaster = Scale(float(film->xResolution),
	                       float(film->yResolution), 1.f) *
		  Scale(1.f / (Screen[1] - Screen[0]),
				1.f / (Screen[2] - Screen[3]), 1.f) *
		 Translate(Vector(-Screen[0], -Screen[3], 0.f));
	RasterToScreen = ScreenToRaster.GetInverse();
	RasterToCamera =
		CameraToScreen.GetInverse() * RasterToScreen;
}
