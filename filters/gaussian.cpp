
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
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
