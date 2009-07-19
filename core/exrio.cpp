
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// exrio.cpp*
#ifdef WIN32
#define hypotf hypot // For the OpenEXR headers
#endif

#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>
#include "pbrt.h"
#include "color.h"
using namespace Imf;
using namespace Imath;
// EXR Function Definitions
COREDLL Spectrum *ReadImage(const string &name, int *width, int *height) {
	try {
	InputFile file(name.c_str());
	Box2i dw = file.header().dataWindow();
	*width  = dw.max.x - dw.min.x + 1;
	*height = dw.max.y - dw.min.y + 1;

	half *rgb = new half[3 * *width * *height];

	FrameBuffer frameBuffer;
	frameBuffer.insert("R", Slice(HALF, (char *)rgb,
		3*sizeof(half), *width * 3 * sizeof(half), 1, 1, 0.0));
	frameBuffer.insert("G", Slice(HALF, (char *)rgb+sizeof(half),
		3*sizeof(half), *width * 3 * sizeof(half), 1, 1, 0.0));
	frameBuffer.insert("B", Slice(HALF, (char *)rgb+2*sizeof(half),
		3*sizeof(half), *width * 3 * sizeof(half), 1, 1, 0.0));

	file.setFrameBuffer(frameBuffer);
	file.readPixels(dw.min.y, dw.max.y);

	Spectrum *ret = new Spectrum[*width * *height];
	// XXX should do real RGB -> Spectrum conversion here
	for (int i = 0; i < *width * *height; ++i) {
		float c[3] = { rgb[3*i], rgb[3*i+1], rgb[3*i+2] };
		ret[i] = Spectrum(c);
	}
	delete[] rgb;
	return ret;
	} catch (const std::exception &e) {
		Error("Unable to read image file \"%s\": %s", name.c_str(),
			e.what());
		return NULL;
	}
}

COREDLL void WriteRGBAImage(const string &name, float *pixels,
		float *alpha, int xRes, int yRes,
		int totalXRes, int totalYRes,
		int xOffset, int yOffset) {
	Header header(totalXRes, totalYRes);
	Box2i dataWindow(V2i(xOffset, yOffset), V2i(xOffset + xRes - 1, yOffset + yRes - 1));
	header.dataWindow() = dataWindow;
	header.channels().insert("R", Channel (HALF));
	header.channels().insert("G", Channel (HALF));
	header.channels().insert("B", Channel (HALF));
	header.channels().insert("A", Channel (HALF));

	half *hrgb = new half[3 * xRes * yRes];
	for (int i = 0; i < 3 * xRes * yRes; ++i)
		hrgb[i] = pixels[i];
	half *ha = new half[xRes * yRes];
	for (int i = 0; i < xRes * yRes; ++i)
		ha[i] = alpha[i];

	hrgb -= 3 * (xOffset + yOffset * xRes);
	ha -= (xOffset + yOffset * xRes);

	FrameBuffer fb;
	fb.insert("R", Slice(HALF, (char *)hrgb, 3*sizeof(half),
		3*xRes*sizeof(half)));
	fb.insert("G", Slice(HALF, (char *)hrgb+sizeof(half), 3*sizeof(half),
		3*xRes*sizeof(half)));
	fb.insert("B", Slice(HALF, (char *)hrgb+2*sizeof(half), 3*sizeof(half),
		3*xRes*sizeof(half)));
	fb.insert("A", Slice(HALF, (char *)ha, sizeof(half), xRes*sizeof(half)));

	OutputFile file(name.c_str(), header);
	file.setFrameBuffer(fb);
	try {
		file.writePixels(yRes);
	}
	catch (const std::exception &e) {
		Error("Unable to write image file \"%s\": %s", name.c_str(),
			e.what());
	}

	delete[] (hrgb + 3 * (xOffset + yOffset * xRes));
	delete[] (ha + (xOffset + yOffset * xRes));
}
