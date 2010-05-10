
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

// imagemap.cpp*
#include "pbrt.h"
#include "texture.h"
#include "mipmap.h"
#include "paramset.h"
#include <map>
using std::map;
// ImageMapTexture Declarations
template <class T>
class ImageTexture : public Texture<T> {
public:
	// ImageTexture Public Methods
	ImageTexture(TextureMapping2D *m,
	             const string &filename,
				 bool doTri,
				 float maxAniso,
				 ImageWrap wm);
	T Evaluate(const DifferentialGeometry &) const;
	~ImageTexture();
private:
	// ImageTexture Private Methods
	static MIPMap<T> *GetTexture(const string &filename,
	    bool doTrilinear, float maxAniso, ImageWrap wm);
	static void convert(const Spectrum &from, Spectrum *to) {
		*to = from;
	}
	static void convert(const Spectrum &from, float *to) {
		*to = from.y();
	}
	// ImageTexture Private Data
	MIPMap<T> *mipmap;
	TextureMapping2D *mapping;
};
// ImageMapTexture Method Definitions
template <class T>
ImageTexture<T>::ImageTexture(TextureMapping2D *m,
		const string &filename,
		bool doTrilinear,
		float maxAniso,
		ImageWrap wrapMode) {
	mapping = m;
	mipmap = GetTexture(filename, doTrilinear,
		maxAniso, wrapMode);
}
template <class T> ImageTexture<T>::~ImageTexture() {
	delete mapping;
}
struct TexInfo {
	TexInfo(const string &f, bool dt, float ma, ImageWrap wm)
		: filename(f), doTrilinear(dt), maxAniso(ma), wrapMode(wm) { }
	string filename;
	bool doTrilinear;
	float maxAniso;
	ImageWrap wrapMode;
	bool operator<(const TexInfo &t2) const {
		if (filename != t2.filename) return filename < t2.filename;
		if (doTrilinear != t2.doTrilinear) return doTrilinear < t2.doTrilinear;
		if (maxAniso != t2.maxAniso) return maxAniso < t2.maxAniso;
		return wrapMode < t2.wrapMode;
	}
};
template <class T> MIPMap<T> *ImageTexture<T>::
	GetTexture( const string &filename,
	           bool doTrilinear,
			   float maxAniso,
			   ImageWrap wrap) {
	// Look for texture in texture cache
	static map<TexInfo, MIPMap<T> *> textures;
	TexInfo texInfo(filename, doTrilinear, maxAniso, wrap);
	if (textures.find(texInfo) != textures.end())
		return textures[texInfo];
	static StatsCounter texLoaded("Texture", // NOBOOK
		"Number of image maps loaded"); // NOBOOK
	++texLoaded; // NOBOOK
	int width, height;
	Spectrum *texels = ReadImage(filename, &width, &height);
	MIPMap<T> *ret = NULL;
	if (texels) {
		// Convert texels to type _T_ and create _MIPMap_
		T *convertedTexels = new T[width*height];
		for (int i = 0; i < width*height; ++i)
			convert(texels[i], &convertedTexels[i]);
		ret = new MIPMap<T>(width, height,
							convertedTexels, doTrilinear,
							maxAniso, wrap);
		delete[] texels;
		delete[] convertedTexels;
	}
	else {
		// Create one-valued _MIPMap_
		T *oneVal = new T[1];
		oneVal[0] = 1.;
		ret = new MIPMap<T>(1, 1, oneVal);
		delete[] oneVal;
	}
	textures[texInfo] = ret;
	return ret;
}
template <class T>
T ImageTexture<T>::Evaluate(
		const DifferentialGeometry &dg) const {
	float s, t, dsdx, dtdx, dsdy, dtdy;
	mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
	return mipmap->Lookup(s, t, dsdx, dtdx, dsdy, dtdy);
}
extern "C" DLLEXPORT Texture<float> * CreateFloatTexture(const Transform &tex2world,
		const TextureParams &tp) {
	// Initialize 2D texture mapping _map_ from _tp_
	TextureMapping2D *map = NULL;
	string type = tp.FindString("mapping");
	if (type == "" || type == "uv") {
		float su = tp.FindFloat("uscale", 1.);
		float sv = tp.FindFloat("vscale", 1.);
		float du = tp.FindFloat("udelta", 0.);
		float dv = tp.FindFloat("vdelta", 0.);
		map = new UVMapping2D(su, sv, du, dv);
	}
	else if (type == "spherical") map = new SphericalMapping2D(tex2world.GetInverse());
	else if (type == "cylindrical") map = new CylindricalMapping2D(tex2world.GetInverse());
	else if (type == "planar")
		map = new PlanarMapping2D(tp.FindVector("v1", Vector(1,0,0)),
			tp.FindVector("v2", Vector(0,1,0)),
			tp.FindFloat("udelta", 0.f), tp.FindFloat("vdelta", 0.f));
	else {
		Error("2D texture mapping \"%s\" unknown", type.c_str());
		map = new UVMapping2D;
	}
	// Initialize _ImageTexture_ parameters
	float maxAniso = tp.FindFloat("maxanisotropy", 8.f);
	bool trilerp = tp.FindBool("trilinear", false);
	string wrap = tp.FindString("wrap");
	ImageWrap wrapMode = TEXTURE_REPEAT;
	if (wrap == "" || wrap == "repeat") wrapMode = TEXTURE_REPEAT;
	else if (wrap == "black") wrapMode = TEXTURE_BLACK;
	else if (wrap == "clamp") wrapMode = TEXTURE_CLAMP;
	return new ImageTexture<float>(map, tp.FindString("filename"),
		trilerp, maxAniso, wrapMode);
}

extern "C" DLLEXPORT Texture<Spectrum> * CreateSpectrumTexture(const Transform &tex2world,
		const TextureParams &tp) {
	// Initialize 2D texture mapping _map_ from _tp_
	TextureMapping2D *map = NULL;
	string type = tp.FindString("mapping");
	if (type == "" || type == "uv") {
		float su = tp.FindFloat("uscale", 1.);
		float sv = tp.FindFloat("vscale", 1.);
		float du = tp.FindFloat("udelta", 0.);
		float dv = tp.FindFloat("vdelta", 0.);
		map = new UVMapping2D(su, sv, du, dv);
	}
	else if (type == "spherical") map = new SphericalMapping2D(tex2world.GetInverse());
	else if (type == "cylindrical") map = new CylindricalMapping2D(tex2world.GetInverse());
	else if (type == "planar")
		map = new PlanarMapping2D(tp.FindVector("v1", Vector(1,0,0)),
			tp.FindVector("v2", Vector(0,1,0)),
			tp.FindFloat("udelta", 0.f), tp.FindFloat("vdelta", 0.f));
	else {
		Error("2D texture mapping \"%s\" unknown", type.c_str());
		map = new UVMapping2D;
	}
	// Initialize _ImageTexture_ parameters
	float maxAniso = tp.FindFloat("maxanisotropy", 8.f);
	bool trilerp = tp.FindBool("trilinear", false);
	string wrap = tp.FindString("wrap");
	ImageWrap wrapMode = TEXTURE_REPEAT;
	if (wrap == "" || wrap == "repeat") wrapMode = TEXTURE_REPEAT;
	else if (wrap == "black") wrapMode = TEXTURE_BLACK;
	else if (wrap == "clamp") wrapMode = TEXTURE_CLAMP;
	return new ImageTexture<Spectrum>(map, tp.FindString("filename"),
		trilerp, maxAniso, wrapMode);
}
