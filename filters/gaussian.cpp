
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

// gaussian.cpp*
#include "sampling.h"
#include "paramset.h"
// Gaussian Filter Declarations
class GaussianFilter : public Filter {
public:
	// GaussianFilter Public Methods
	GaussianFilter(float xw, float yw, float a)
		: Filter(xw, yw) {
		alpha = a;
		expX = expf(-alpha * xWidth * xWidth);
		expY = expf(-alpha * yWidth * yWidth);
	}
	float Evaluate(float x, float y) const;
private:
	// GaussianFilter Private Data
	float alpha;
	float expX, expY;
	// GaussianFilter Utility Functions
	float Gaussian(float d, float expv) const {
		return max(0.f, float(expf(-alpha * d * d) - expv));
	}
};
// Gaussian Filter Method Definitions
float GaussianFilter::Evaluate(float x, float y) const {
	return Gaussian(x, expX) * Gaussian(y, expY);
}
extern "C" DLLEXPORT Filter *CreateFilter(const ParamSet &ps) {
	// Find common filter parameters
	float xw = ps.FindOneFloat("xwidth", 2.);
	float yw = ps.FindOneFloat("ywidth", 2.);
	float alpha = ps.FindOneFloat("alpha", 2.f);
	return new GaussianFilter(xw, yw, alpha);
}
