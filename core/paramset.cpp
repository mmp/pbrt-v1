
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

// paramset.cpp*
#include "paramset.h"
// ParamSet Methods
ParamSet::ParamSet(const ParamSet &p2) {
	*this = p2;
}
ParamSet &ParamSet::operator=(const ParamSet &p2) {
	if (&p2 != this) {
		Clear();
		u_int i;
		for (i = 0; i < p2.ints.size(); ++i)
			ints.push_back(p2.ints[i]->Clone());
		for (i = 0; i < p2.bools.size(); ++i)
			bools.push_back(p2.bools[i]->Clone());
		for (i = 0; i < p2.floats.size(); ++i)
			floats.push_back(p2.floats[i]->Clone());
		for (i = 0; i < p2.points.size(); ++i)
			points.push_back(p2.points[i]->Clone());
		for (i = 0; i < p2.vectors.size(); ++i)
			vectors.push_back(p2.vectors[i]->Clone());
		for (i = 0; i < p2.normals.size(); ++i)
			normals.push_back(p2.normals[i]->Clone());
		for (i = 0; i < p2.spectra.size(); ++i)
			spectra.push_back(p2.spectra[i]->Clone());
		for (i = 0; i < p2.strings.size(); ++i)
			strings.push_back(p2.strings[i]->Clone());
		for (i = 0; i < p2.textures.size(); ++i)
			textures.push_back(p2.textures[i]->Clone());
	}
	return *this;
}
void ParamSet::AddFloat(const string &name,
			            const float *data,
						int nItems) {
	EraseFloat(name);
	floats.push_back(new ParamSetItem<float>(name,
											 data,
											 nItems));
}
void ParamSet::AddInt(const string &name, const int *data, int nItems) {
	EraseInt(name);
	ADD_PARAM_TYPE(int, ints);
}
void ParamSet::AddBool(const string &name, const bool *data, int nItems) {
	EraseInt(name);
	ADD_PARAM_TYPE(bool, bools);
}
void ParamSet::AddPoint(const string &name, const Point *data, int nItems) {
	ErasePoint(name);
	ADD_PARAM_TYPE(Point, points);
}
void ParamSet::AddVector(const string &name, const Vector *data, int nItems) {
	EraseVector(name);
	ADD_PARAM_TYPE(Vector, vectors);
}
void ParamSet::AddNormal(const string &name, const Normal *data, int nItems) {
	EraseNormal(name);
	ADD_PARAM_TYPE(Normal, normals);
}
void ParamSet::AddSpectrum(const string &name, const Spectrum *data, int nItems) {
	EraseSpectrum(name);
	ADD_PARAM_TYPE(Spectrum, spectra);
}
void ParamSet::AddString(const string &name, const string *data, int nItems) {
	EraseString(name);
	ADD_PARAM_TYPE(string, strings);
}
void ParamSet::AddTexture(const string &name, const string &value) {
	EraseTexture(name);
	textures.push_back(new ParamSetItem<string>(name, (string *)&value, 1));
}
bool ParamSet::EraseInt(const string &n) {
	for (u_int i = 0; i < ints.size(); ++i)
		if (ints[i]->name == n) {
			delete ints[i];
			ints.erase(ints.begin() + i);
			return true;
		}
	return false;
}
bool ParamSet::EraseBool(const string &n) {
	for (u_int i = 0; i < bools.size(); ++i)
		if (bools[i]->name == n) {
			delete bools[i];
			bools.erase(bools.begin() + i);
			return true;
		}
	return false;
}
bool ParamSet::EraseFloat(const string &n) {
	for (u_int i = 0; i < floats.size(); ++i)
		if (floats[i]->name == n) {
			delete floats[i];
			floats.erase(floats.begin() + i);
			return true;
		}
	return false;
}
bool ParamSet::ErasePoint(const string &n) {
	for (u_int i = 0; i < points.size(); ++i)
		if (points[i]->name == n) {
			delete points[i];
			points.erase(points.begin() + i);
			return true;
		}
	return false;
}
bool ParamSet::EraseVector(const string &n) {
	for (u_int i = 0; i < vectors.size(); ++i)
		if (vectors[i]->name == n) {
			delete vectors[i];
			vectors.erase(vectors.begin() + i);
			return true;
		}
	return false;
}
bool ParamSet::EraseNormal(const string &n) {
	for (u_int i = 0; i < normals.size(); ++i)
		if (normals[i]->name == n) {
			delete normals[i];
			normals.erase(normals.begin() + i);
			return true;
		}
	return false;
}
bool ParamSet::EraseSpectrum(const string &n) {
	for (u_int i = 0; i < spectra.size(); ++i)
		if (spectra[i]->name == n) {
			delete spectra[i];
			spectra.erase(spectra.begin() + i);
			return true;
		}
	return false;
}
bool ParamSet::EraseString(const string &n) {
	for (u_int i = 0; i < strings.size(); ++i)
		if (strings[i]->name == n) {
			delete strings[i];
			strings.erase(strings.begin() + i);
			return true;
		}
	return false;
}
bool ParamSet::EraseTexture(const string &n) {
	for (u_int i = 0; i < textures.size(); ++i)
		if (textures[i]->name == n) {
			delete textures[i];
			textures.erase(textures.begin() + i);
			return true;
		}
	return false;
}
float ParamSet::FindOneFloat(const string &name,
                             float d) const {
	for (u_int i = 0; i < floats.size(); ++i)
		if (floats[i]->name == name &&
			floats[i]->nItems == 1) {
			floats[i]->lookedUp = true;
			return *(floats[i]->data);
		}
	return d;
}
const float *ParamSet::FindFloat(const string &name,
		int *nItems) const {
	for (u_int i = 0; i < floats.size(); ++i)
		if (floats[i]->name == name) {
			*nItems = floats[i]->nItems;
			floats[i]->lookedUp = true;
			return floats[i]->data;
		}
	return NULL;
}
const int *ParamSet::FindInt(const string &name, int *nItems) const {
	LOOKUP_PTR(ints);
}
const bool *ParamSet::FindBool(const string &name, int *nItems) const {
	LOOKUP_PTR(bools);
}
int ParamSet::FindOneInt(const string &name, int d) const {
	LOOKUP_ONE(ints);
}
bool ParamSet::FindOneBool(const string &name, bool d) const {
	LOOKUP_ONE(bools);
}
const Point *ParamSet::FindPoint(const string &name, int *nItems) const {
	LOOKUP_PTR(points);
}
Point ParamSet::FindOnePoint(const string &name, const Point &d) const {
	LOOKUP_ONE(points);
}
const Vector *ParamSet::FindVector(const string &name, int *nItems) const {
	LOOKUP_PTR(vectors);
}
Vector ParamSet::FindOneVector(const string &name, const Vector &d) const {
	LOOKUP_ONE(vectors);
}
const Normal *ParamSet::FindNormal(const string &name, int *nItems) const {
	LOOKUP_PTR(normals);
}
Normal ParamSet::FindOneNormal(const string &name, const Normal &d) const {
	LOOKUP_ONE(normals);
}
const Spectrum *ParamSet::FindSpectrum(const string &name, int *nItems) const {
	LOOKUP_PTR(spectra);
}
Spectrum ParamSet::FindOneSpectrum(const string &name, const Spectrum &d) const {
	LOOKUP_ONE(spectra);
}
const string *ParamSet::FindString(const string &name, int *nItems) const {
	LOOKUP_PTR(strings);
}
string ParamSet::FindOneString(const string &name, const string &d) const {
	LOOKUP_ONE(strings);
}
string ParamSet::FindTexture(const string &name) const {
	string d = "";
	LOOKUP_ONE(textures);
}
void ParamSet::ReportUnused() const {
#define CHECK_UNUSED(v) \
	for (i = 0; i < (v).size(); ++i) \
		if (!(v)[i]->lookedUp) \
			Warning("Parameter \"%s\" not used", \
				(v)[i]->name.c_str())
	u_int i;
	CHECK_UNUSED(ints);       CHECK_UNUSED(bools);
	CHECK_UNUSED(floats);   CHECK_UNUSED(points);
	CHECK_UNUSED(vectors); CHECK_UNUSED(normals);
	CHECK_UNUSED(spectra); CHECK_UNUSED(strings);
	CHECK_UNUSED(textures);
}
void ParamSet::Clear() {
#define DEL_PARAMS(name) \
	for (u_int i = 0; i < (name).size(); ++i) \
		delete (name)[i]; \
	(name).erase((name).begin(), (name).end())
	DEL_PARAMS(ints);    DEL_PARAMS(bools);
	DEL_PARAMS(floats);  DEL_PARAMS(points);
	DEL_PARAMS(vectors); DEL_PARAMS(normals);
	DEL_PARAMS(spectra); DEL_PARAMS(strings);
	DEL_PARAMS(textures);
#undef DEL_PARAMS
}
string ParamSet::ToString() const {
	string ret;
	u_int i;
	int j;
	string typeString;
	const int bufLen = 48*1024*1024;
	static char *buf = new char[bufLen];
	char *bufEnd = buf + bufLen;
	for (i = 0; i < ints.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<int> *item = ints[i];
		typeString = "integer ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "%d ", item->data[j]);
		ret += buf;
		ret += string("] ");
	}
	for (i = 0; i < bools.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<bool> *item = bools[i];
		typeString = "bool ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "\"%s\" ", item->data[j] ? "true" : "false");
		ret += buf;
		ret += string("] ");
	}
	for (i = 0; i < floats.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<float> *item = floats[i];
		typeString = "float ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "%.8g ", item->data[j]);
		ret += buf;
		ret += string("] ");
	}
	for (i = 0; i < points.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<Point> *item = points[i];
		typeString = "point ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "%.8g %.8g %.8g ", item->data[j].x,
				item->data[j].y, item->data[j].z);
		ret += buf;
		ret += string("] ");
	}
	for (i = 0; i < vectors.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<Vector> *item = vectors[i];
		typeString = "vector ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "%.8g %.8g %.8g ", item->data[j].x,
				item->data[j].y, item->data[j].z);
		ret += buf;
		ret += string("] ");
	}
	for (i = 0; i < normals.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<Normal> *item = normals[i];
		typeString = "normal ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "%.8g %.8g %.8g ", item->data[j].x,
				item->data[j].y, item->data[j].z);
		ret += buf;
		ret += string("] ");
	}
	for (i = 0; i < strings.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<string> *item = strings[i];
		typeString = "string ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "\"%s\" ", item->data[j].c_str());
		ret += buf;
		ret += string("] ");
	}
	for (i = 0; i < textures.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<string> *item = textures[i];
		typeString = "texture ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "\"%s\" ", item->data[j].c_str());
		ret += buf;
		ret += string("] ");
	}
	for (i = 0; i < spectra.size(); ++i) {
		char *bufp = buf;
		*bufp = '\0';
		ParamSetItem<Spectrum> *item = spectra[i];
		typeString = "color ";
		// Print _ParamSetItem_ declaration, determine how many to print
		int nPrint = item->nItems;
		ret += string("\"");
		ret += typeString;
		ret += item->name;
		ret += string("\"");
		ret += string(" [");
		for (j = 0; j < nPrint; ++j)
			bufp += snprintf(bufp, bufEnd - bufp, "%.8g %.8g %.8g ", item->data[j].c[0],
				item->data[j].c[1], item->data[j].c[2]);
		ret += buf;
		ret += string("] ");
	}
	return ret;
}
// TextureParams Method Definitions
Reference<Texture<Spectrum> >
	TextureParams::GetSpectrumTexture(const string &n,
             const Spectrum &def) const {
	string name = geomParams.FindTexture(n);
	if (name == "") name = materialParams.FindTexture(n);
	if (name != "") {
		if (spectrumTextures.find(name) !=
		       spectrumTextures.end())
			return spectrumTextures[name];
		else
			Error("Couldn't find spectrum"
			      "texture named \"%s\"", n.c_str());
	}
	Spectrum val = geomParams.FindOneSpectrum(n,
		materialParams.FindOneSpectrum(n, def));
	return new ConstantTexture<Spectrum>(val);
}
Reference<Texture<float> > TextureParams::GetFloatTexture(const string &n,
		float def) const {
	string name = geomParams.FindTexture(n);
	if (name == "") name = materialParams.FindTexture(n);
	if (name != "") {
		if (floatTextures.find(name) != floatTextures.end())
			return floatTextures[name];
		else
			Error("Couldn't find float texture named \"%s\"", n.c_str());
	}
	float val = geomParams.FindOneFloat(n,
		materialParams.FindOneFloat(n, def));
	return new ConstantTexture<float>(val);
}
