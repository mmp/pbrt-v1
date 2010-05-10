
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

// maxwhite.cpp*
#include "tonemap.h"
// MaxWhiteOp Declarations
class MaxWhiteOp : public ToneMap {
public:
	// MaxWhiteOp Public Methods
	void Map(const float *y, int xRes, int yRes,
	         float maxDisplayY, float *scale) const {
		// Compute maximum luminance of all pixels
		float maxY = 0.;
		for (int i = 0; i < xRes * yRes; ++i)
			maxY = max(maxY, y[i]);
		float s = maxDisplayY / maxY;
		for (int i = 0; i < xRes * yRes; ++i)
			scale[i] = s;
	}
};
// MaxWhiteOp Method Definitions
extern "C" DLLEXPORT ToneMap *CreateToneMap(const ParamSet &ps) {
	return new MaxWhiteOp;
}
