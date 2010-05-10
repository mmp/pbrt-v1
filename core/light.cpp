
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
