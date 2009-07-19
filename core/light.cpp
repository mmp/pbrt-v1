
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// light.cpp*
#include "light.h"
#include "scene.h"
// Light Method Definitions
Light::~Light() {
}
bool VisibilityTester::Unoccluded(const Scene *scene) const {
	// Update shadow ray statistics
	static StatsCounter nShadowRays("Lights",
		"Number of shadow rays traced");
	++nShadowRays;
	return !scene->IntersectP(r);
}
Spectrum VisibilityTester::
	Transmittance(const Scene *scene) const {
	return scene->Transmittance(r);
}
Spectrum Light::Le(const RayDifferential &) const {
	return Spectrum(0.);
}
