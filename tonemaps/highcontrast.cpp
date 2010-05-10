
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

// highcontrast.cpp*
#include "tonemap.h"
#include "mipmap.h"
// HighContrastOp Declarations
class HighContrastOp : public ToneMap {
public:
	void Map(const float *y,
	         int xRes, int yRes,
			 float maxDisplayY, float *scale) const;
private:
	// HighContrastOp Utility Methods
	static float C(float y) {
		if (y < 0.0034f)
			return y / 0.0014f;
		else if (y < 1)
			return 2.4483f + log10f(y/0.0034f)/0.4027f;
		else if (y < 7.2444f)
			return 16.563f + (y - 1)/0.4027f;
		else
			return 32.0693f + log10f(y / 7.2444f)/0.0556f;
	}
	static float T(float y, float CYmin, float CYmax,
			float maxDisplayY) {
		return maxDisplayY * (C(y) - CYmin) / (CYmax - CYmin);
	}
};
// HighContrastOp Method Definitions
void HighContrastOp::Map(const float *y, int xRes, int yRes,
		float maxDisplayY, float *scale) const {
	// Find minimum and maximum image luminances
	float minY = y[0], maxY = y[0];
	for (int i = 0; i < xRes * yRes; ++i) {
		minY = min(minY, y[i]);
		maxY = max(maxY, y[i]);
	}
	float CYmin = C(minY), CYmax = C(maxY);
	// Build luminance image pyramid
	MIPMap<float> pyramid(xRes, yRes,
	                      y, false,
						  4.f, TEXTURE_CLAMP);
	// Apply high contrast tone mapping operator
	ProgressReporter progress(xRes*yRes, "Tone Mapping"); // NOBOOK
	for (int y = 0; y < yRes; ++y) {
		float yc = (float(y) + .5f) / float(yRes);
		for (int x = 0; x < xRes; ++x) {
			float xc = (float(x) + .5f) / float(xRes);
			// Compute local adaptation luminance at $(x,y)$
			float dwidth = 1.f / float(max(xRes, yRes));
			float maxWidth = 32.f / float(max(xRes, yRes));
			float width = dwidth, prevWidth = 0.f;
			float Yadapt;
			float prevlc = 0.f;
			const float maxLocalContrast = .5f;
			while (1) {
				// Compute local contrast at $(x,y)$
				float b0 = pyramid.Lookup(xc, yc, width,
				                          0.f, 0.f, width);
				float b1 = pyramid.Lookup(xc, yc, 2.f*width,
				                          0.f, 0.f, 2.f*width);
				float lc = fabsf((b0 - b1) / b0);
				// If maximum contrast is exceeded, compute adaptation luminance
				if (lc > maxLocalContrast) {
					float t = (maxLocalContrast - prevlc) / (lc - prevlc);
					float w = Lerp(t, prevWidth, width);
					Yadapt = pyramid.Lookup(xc, yc, w,
					                        0.f, 0.f, w);
					break;
				}
				// Increase search region and prepare to compute contrast again
				prevlc = lc;
				prevWidth = width;
				width += dwidth;
				if (width >= maxWidth) {
					Yadapt = pyramid.Lookup(xc, yc, maxWidth,
					                        0.f, 0.f, maxWidth);
					break;
				}
			}
			// Apply tone mapping based on local adaptation luminance
			scale[x + y*xRes] = T(Yadapt, CYmin, CYmax, maxDisplayY) /
				Yadapt;
		}
		progress.Update(xRes); // NOBOOK
	}
	progress.Done(); // NOBOOK
}
extern "C" DLLEXPORT ToneMap *CreateToneMap(const ParamSet &ps) {
	return new HighContrastOp;
}
