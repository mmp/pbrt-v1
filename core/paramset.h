
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

#ifndef PBRT_PARAMSET_H
#define PBRT_PARAMSET_H
// paramset.h*
#include "pbrt.h"
#include "geometry.h"
#include "texture.h"
#include "color.h"
#include <map>
using std::map;
#if (_MSC_VER >= 1400) // NOBOOK
#include <stdio.h>     // NOBOOK
#define snprintf _snprintf // NOBOOK
#endif // NOBOOK
// ParamSet Macros
#define ADD_PARAM_TYPE(T, vec) \
	(vec).push_back(new ParamSetItem<T>(name, (T *)data, nItems))
#define LOOKUP_PTR(vec) \
	for (u_int i = 0; i < (vec).size(); ++i) \
		if ((vec)[i]->name == name) { \
			*nItems = (vec)[i]->nItems; \
			(vec)[i]->lookedUp = true; \
			return (vec)[i]->data; \
		} \
	return NULL
#define LOOKUP_ONE(vec) \
	for (u_int i = 0; i < (vec).size(); ++i) { \
		if ((vec)[i]->name == name && \
			(vec)[i]->nItems == 1) { \
			(vec)[i]->lookedUp = true; \
			return *((vec)[i]->data); \
}		} \
	return d
// ParamSet Declarations
class COREDLL ParamSet {
public:
	// ParamSet Public Methods
	ParamSet() { }
	ParamSet &operator=(const ParamSet &p2);
	ParamSet(const ParamSet &p2);
	void AddFloat(const string &, const float *, int nItems = 1);
	void AddInt(const string &,
	            const int *,
				int nItems = 1);
	void AddBool(const string &,
	             const bool *,
				 int nItems = 1);
	void AddPoint(const string &,
	              const Point *,
				  int nItems = 1);
	void AddVector(const string &,
	               const Vector *,
				   int nItems = 1);
	void AddNormal(const string &,
	               const Normal *,
				   int nItems = 1);
	void AddSpectrum(const string &,
	                const Spectrum *,
					int nItems = 1);
	void AddString(const string &,
	              const string *,
				  int nItems = 1);
	void AddTexture(const string &,
	                const string &);
	bool EraseInt(const string &);
	bool EraseBool(const string &);
	bool EraseFloat(const string &);
	bool ErasePoint(const string &);
	bool EraseVector(const string &);
	bool EraseNormal(const string &);
	bool EraseSpectrum(const string &);
	bool EraseString(const string &);
	bool EraseTexture(const string &);
	float FindOneFloat(const string &, float d) const;
	int FindOneInt(const string &, int d) const;
	bool FindOneBool(const string &, bool d) const;
	Point FindOnePoint(const string &, const Point &d) const;
	Vector FindOneVector(const string &, const Vector &d) const;
	Normal FindOneNormal(const string &, const Normal &d) const;
	Spectrum FindOneSpectrum(const string &,
		const Spectrum &d) const;
	string FindOneString(const string &, const string &d) const;
	string FindTexture(const string &) const;
	const float *FindFloat(const string &, int *nItems) const;
	const int *FindInt(const string &, int *nItems) const;
	const bool *FindBool(const string &, int *nItems) const;
	const Point *FindPoint(const string &, int *nItems) const;
	const Vector *FindVector(const string &, int *nItems) const;
	const Normal *FindNormal(const string &, int *nItems) const;
	const Spectrum *FindSpectrum(const string &,
		int *nItems) const;
	const string *FindString(const string &,
		int *nItems) const;
	void ReportUnused() const;
	~ParamSet() {
		Clear();
	}
	void Clear();
	string ToString() const;
private:
	// ParamSet Data
	vector<ParamSetItem<int> *> ints;
	vector<ParamSetItem<bool> *> bools;
	vector<ParamSetItem<float> *> floats;
	vector<ParamSetItem<Point> *> points;
	vector<ParamSetItem<Vector> *> vectors;
	vector<ParamSetItem<Normal> *> normals;
	vector<ParamSetItem<Spectrum> *> spectra;
	vector<ParamSetItem<string> *> strings;
	vector<ParamSetItem<string> *> textures;
};
template <class T> struct ParamSetItem {
	// ParamSetItem Public Methods
	ParamSetItem<T> *Clone() const {
		return new ParamSetItem<T>(name, data, nItems);
	}
	ParamSetItem(const string &name, const T *val, int nItems = 1);
	~ParamSetItem() {
		delete[] data;
	}
	// ParamSetItem Data
	string name;
	int nItems;
	T *data;
	bool lookedUp;
};
// ParamSetItem Methods
template <class T>
ParamSetItem<T>::ParamSetItem(const string &n,
							  const T *v,
							  int ni) {
	name = n;
	nItems = ni;
	data = new T[nItems];
	for (int i = 0; i < nItems; ++i)
		data[i] = v[i];
	lookedUp = false;
}
// TextureParams Declarations
class COREDLL TextureParams {
public:
	// TextureParams Public Methods
	TextureParams(const ParamSet &geomp, const ParamSet &matp,
			map<string, Reference<Texture<float> > > &ft,
			map<string, Reference<Texture<Spectrum> > > &st)
		: geomParams(geomp),
		  materialParams(matp),
		  floatTextures(ft),
		  spectrumTextures(st) {
	}
	Reference<Texture<Spectrum> > GetSpectrumTexture(const string &name,
			const Spectrum &def) const;
	Reference<Texture<float> > GetFloatTexture(const string &name,
			float def) const;
	float FindFloat(const string &n, float d) const {
		return geomParams.FindOneFloat(n,
			materialParams.FindOneFloat(n, d));
	}
	string FindString(const string &n) const {
	       return geomParams.FindOneString(n, materialParams.FindOneString(n, ""));
	}
	int FindInt(const string &n, int d) const {
	       return geomParams.FindOneInt(n, materialParams.FindOneInt(n, d));
	}
	bool FindBool(const string &n, bool d) const {
	       return geomParams.FindOneBool(n, materialParams.FindOneBool(n, d));
	}
	Point FindPoint(const string &n, const Point &d) const {
	       return geomParams.FindOnePoint(n, materialParams.FindOnePoint(n, d));
	}
	Vector FindVector(const string &n, const Vector &d) const {
	       return geomParams.FindOneVector(n, materialParams.FindOneVector(n, d));
	}
	Normal FindNormal(const string &n, const Normal &d) const {
	       return geomParams.FindOneNormal(n, materialParams.FindOneNormal(n, d));
	}
	Spectrum FindSpectrum(const string &n, const Spectrum &d) const {
	       return geomParams.FindOneSpectrum(n, materialParams.FindOneSpectrum(n, d));
	}
	void ReportUnused() const {
		geomParams.ReportUnused();
		materialParams.ReportUnused();
	}
	const ParamSet &GetGeomParams() const { return geomParams; }
	const ParamSet &GetMaterialParams() const { return materialParams; }
private:
	// TextureParams Private Data
	const ParamSet &geomParams, &materialParams;
	map<string,
	    Reference<Texture<float> > > &floatTextures;
	map<string,
	    Reference<Texture<Spectrum> > > &spectrumTextures;
};
#endif // PBRT_PARAMSET_H
