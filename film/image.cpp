
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

// image.cpp*
#include "pbrt.h"
#include "film.h"
#include "color.h"
#include "paramset.h"
#include "tonemap.h"
#include "sampling.h"
// ImageFilm Declarations
class ImageFilm : public Film {
public:
	// ImageFilm Public Methods
	ImageFilm(int xres, int yres,
	                     Filter *filt, const float crop[4],
		             const string &filename, bool premult,
		             int wf);
	~ImageFilm() {
		delete pixels;
		delete filter;
		delete[] filterTable;
	}
	void AddSample(const Sample &sample, const Ray &ray,
	               const Spectrum &L, float alpha);
	void GetSampleExtent(int *xstart, int *xend,
	                     int *ystart, int *yend) const;
	void WriteImage();
private:
	// ImageFilm Private Data
	Filter *filter;
	int writeFrequency, sampleCount;
	string filename;
	bool premultiplyAlpha;
	float cropWindow[4];
	int xPixelStart, yPixelStart, xPixelCount, yPixelCount;
	struct Pixel {
		Pixel() : L(0.f) {
			alpha = 0.f;
			weightSum = 0.f;
		}
		Spectrum L;
		float alpha, weightSum;
	};
	BlockedArray<Pixel> *pixels;
	float *filterTable;
};
// ImageFilm Method Definitions
ImageFilm::ImageFilm(int xres, int yres,
                     Filter *filt, const float crop[4],
		             const string &fn, bool premult, int wf)
	: Film(xres, yres) {
	filter = filt;
	memcpy(cropWindow, crop, 4 * sizeof(float));
	filename = fn;
	premultiplyAlpha = premult;
	writeFrequency = sampleCount = wf;
	// Compute film image extent
	xPixelStart = Ceil2Int(xResolution * cropWindow[0]);
	xPixelCount =
		max(1, Ceil2Int(xResolution * cropWindow[1]) - xPixelStart);
	yPixelStart =
		Ceil2Int(yResolution * cropWindow[2]);
	yPixelCount =
		max(1, Ceil2Int(yResolution * cropWindow[3]) - yPixelStart);
	// Allocate film image storage
	pixels = new BlockedArray<Pixel>(xPixelCount, yPixelCount);
	// Precompute filter weight table
	#define FILTER_TABLE_SIZE 16
	filterTable =
		new float[FILTER_TABLE_SIZE * FILTER_TABLE_SIZE];
	float *ftp = filterTable;
	for (int y = 0; y < FILTER_TABLE_SIZE; ++y) {
		float fy = ((float)y + .5f) * filter->yWidth /
			FILTER_TABLE_SIZE;
		for (int x = 0; x < FILTER_TABLE_SIZE; ++x) {
			float fx = ((float)x + .5f) * filter->xWidth /
				FILTER_TABLE_SIZE;
			*ftp++ = filter->Evaluate(fx, fy);
		}
	}
}
void ImageFilm::AddSample(const Sample &sample,
		const Ray &ray, const Spectrum &L, float alpha) {
	// Compute sample's raster extent
	float dImageX = sample.imageX - 0.5f;
	float dImageY = sample.imageY - 0.5f;
	int x0 = Ceil2Int (dImageX - filter->xWidth);
	int x1 = Floor2Int(dImageX + filter->xWidth);
	int y0 = Ceil2Int (dImageY - filter->yWidth);
	int y1 = Floor2Int(dImageY + filter->yWidth);
	x0 = max(x0, xPixelStart);
	x1 = min(x1, xPixelStart + xPixelCount - 1);
	y0 = max(y0, yPixelStart);
	y1 = min(y1, yPixelStart + yPixelCount - 1);
	if ((x1-x0) < 0 || (y1-y0) < 0) return;
	// Loop over filter support and add sample to pixel arrays
	// Precompute $x$ and $y$ filter table offsets
	int *ifx = (int *)alloca((x1-x0+1) * sizeof(int));
	for (int x = x0; x <= x1; ++x) {
		float fx = fabsf((x - dImageX) *
		                 filter->invXWidth * FILTER_TABLE_SIZE);
		ifx[x-x0] = min(Floor2Int(fx), FILTER_TABLE_SIZE-1);
	}
	int *ify = (int *)alloca((y1-y0+1) * sizeof(int));
	for (int y = y0; y <= y1; ++y) {
		float fy = fabsf((y - dImageY) *
		                 filter->invYWidth * FILTER_TABLE_SIZE);
		ify[y-y0] = min(Floor2Int(fy), FILTER_TABLE_SIZE-1);
	}
	for (int y = y0; y <= y1; ++y)
		for (int x = x0; x <= x1; ++x) {
			// Evaluate filter value at $(x,y)$ pixel
			int offset = ify[y-y0]*FILTER_TABLE_SIZE + ifx[x-x0];
			float filterWt = filterTable[offset];
			// Update pixel values with filtered sample contribution
			Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
			pixel.L.AddWeighted(filterWt, L);
			pixel.alpha += alpha * filterWt;
			pixel.weightSum += filterWt;
		}
	// Possibly write out in-progress image
	if (--sampleCount == 0) {
		WriteImage();
		sampleCount = writeFrequency;
	}
}
void ImageFilm::GetSampleExtent(int *xstart,
		int *xend, int *ystart, int *yend) const {
	*xstart = Floor2Int(xPixelStart + .5f - filter->xWidth);
	*xend   = Floor2Int(xPixelStart + .5f + xPixelCount  +
		filter->xWidth);
	*ystart = Floor2Int(yPixelStart + .5f - filter->yWidth);
	*yend   = Floor2Int(yPixelStart + .5f + yPixelCount +
		filter->yWidth);
}
void ImageFilm::WriteImage() {
	// Convert image to RGB and compute final pixel values
	int nPix = xPixelCount * yPixelCount;
	float *rgb = new float[3*nPix], *alpha = new float[nPix];
	int offset = 0;
	for (int y = 0; y < yPixelCount; ++y) {
		for (int x = 0; x < xPixelCount; ++x) {
			// Convert pixel spectral radiance to RGB
			float xyz[3];
			(*pixels)(x, y).L.XYZ(xyz);
			const float
				rWeight[3] = { 3.240479f, -1.537150f, -0.498535f };
			const float
				gWeight[3] = {-0.969256f,  1.875991f,  0.041556f };
			const float
				bWeight[3] = { 0.055648f, -0.204043f,  1.057311f };
			rgb[3*offset  ] = rWeight[0]*xyz[0] +
			                  rWeight[1]*xyz[1] +
				              rWeight[2]*xyz[2];
			rgb[3*offset+1] = gWeight[0]*xyz[0] +
			                  gWeight[1]*xyz[1] +
				              gWeight[2]*xyz[2];
			rgb[3*offset+2] = bWeight[0]*xyz[0] +
			                  bWeight[1]*xyz[1] +
				              bWeight[2]*xyz[2];
			alpha[offset] = (*pixels)(x, y).alpha;
			// Normalize pixel with weight sum
			float weightSum = (*pixels)(x, y).weightSum;
			if (weightSum != 0.f) {
				float invWt = 1.f / weightSum;
				rgb[3*offset  ] =
					Clamp(rgb[3*offset  ] * invWt, 0.f, INFINITY);
				rgb[3*offset+1] =
					Clamp(rgb[3*offset+1] * invWt, 0.f, INFINITY);
				rgb[3*offset+2] =
					Clamp(rgb[3*offset+2] * invWt, 0.f, INFINITY);
				alpha[offset] = Clamp(alpha[offset] * invWt, 0.f, 1.f);
			}
			// Compute premultiplied alpha color
			if (premultiplyAlpha) {
				rgb[3*offset  ] *= alpha[offset];
				rgb[3*offset+1] *= alpha[offset];
				rgb[3*offset+2] *= alpha[offset];
			}
			++offset;
		}
	}
	// Write RGBA image
	WriteRGBAImage(filename, rgb, alpha,
		xPixelCount, yPixelCount,
		xResolution, yResolution,
		xPixelStart, yPixelStart);
	// Release temporary image memory
	delete[] alpha;
	delete[] rgb;
}
extern "C" DLLEXPORT Film *CreateFilm(const ParamSet &params, Filter *filter)
{
	string filename = params.FindOneString("filename", "pbrt.exr");
	bool premultiplyAlpha = params.FindOneBool("premultiplyalpha", true);

	int xres = params.FindOneInt("xresolution", 640);
	int yres = params.FindOneInt("yresolution", 480);
	float crop[4] = { 0, 1, 0, 1 };
	int cwi;
	const float *cr = params.FindFloat("cropwindow", &cwi);
	if (cr && cwi == 4) {
		crop[0] = Clamp(min(cr[0], cr[1]), 0., 1.);
		crop[1] = Clamp(max(cr[0], cr[1]), 0., 1.);
		crop[2] = Clamp(min(cr[2], cr[3]), 0., 1.);
		crop[3] = Clamp(max(cr[2], cr[3]), 0., 1.);
	}
	int writeFrequency = params.FindOneInt("writefrequency", -1);

	return new ImageFilm(xres, yres, filter, crop,
		filename, premultiplyAlpha, writeFrequency);
}
