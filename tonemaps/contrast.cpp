
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

// contrast.cpp*
#include "tonemap.h"
#include "paramset.h"
// ContrastOp Declarations
class ContrastOp : public ToneMap {
public:
	ContrastOp(float day) { displayAdaptationY = day; }
	void Map(const float *y,
	         int xRes, int yRes,
			 float maxDisplayY, float *scale) const;
	float displayAdaptationY;
};
// ContrastOp Method Definitions
void ContrastOp::Map(const float *y, int xRes, int yRes,
		float maxDisplayY, float *scale) const {
	// Compute world adaptation luminance, _Ywa_
	float Ywa = 0.;
	for (int i = 0; i < xRes * yRes; ++i)
		if (y[i] > 0) Ywa += logf(y[i]);
	Ywa = expf(Ywa / (xRes * yRes));
	// Compute contrast-preserving scalefactor, _s_
	float s = powf((1.219f + powf(displayAdaptationY, 0.4f)) /
		(1.219f + powf(Ywa, 0.4f)), 2.5f);
	for (int i = 0; i < xRes*yRes; ++i)
		scale[i] = s;
}
extern "C" DLLEXPORT ToneMap *CreateToneMap(const ParamSet &ps) {
	float day = ps.FindOneFloat("displayadaptationY", 50.f);
	return new ContrastOp(day);
}
