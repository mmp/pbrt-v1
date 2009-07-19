/*
 * The 3D Studio File Format Library
 * Copyright (C) 1996-2001 by J.E. Hoffmann <je-h@gmx.net>
 * All rights reserved.
 *
 * This program is  free  software;  you can redistribute it and/or modify it
 * under the terms of the  GNU Lesser General Public License  as published by 
 * the  Free Software Foundation;  either version 2.1 of the License,  or (at 
 * your option) any later version.
 *
 * This  program  is  distributed in  the  hope that it will  be useful,  but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or  FITNESS FOR A  PARTICULAR PURPOSE.  See the  GNU Lesser General Public  
 * License for more details.
 *
 * You should  have received  a copy of the GNU Lesser General Public License
 * along with  this program;  if not, write to the  Free Software Foundation,
 * Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * $Id: 3ds2lrt.c 1063 2004-05-14 01:23:50Z mmp $
 */
#include <lib3ds/file.h>
#include <lib3ds/vector.h>
#include <lib3ds/matrix.h>
#include <lib3ds/camera.h>
#include <lib3ds/light.h>
#include <lib3ds/material.h>
#include <lib3ds/mesh.h>
#include <lib3ds/node.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
//#include <config.h>
#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif
#include <math.h>

#define max(a,b) ((a)>(b) ? (a) : (b))

/*!
\example 3ds2rib.c

Simple <i>3DS</i> to RIB (RenderMan Interface Bytestream) converter.

\code
Syntax: 3ds2rib [options] filename [options]

Options:\n"
  -h/--help               This help
  -o/--output <filename>  Write output to <filename> instead of stdout
  -c/--camera <name>      Use camera <name> for rendering
  -a/--all                Render all frames
  -f/--frame #            Render frame #
  -d/--downcuts #         Render # frames of animation
\endcode

\author J.E. Hoffmann <je-h@gmx.net>
*/

#define VERSION "MMP-1.0"

static void
help()
{
  fprintf(stderr,
"The 3D Studio File Format Library - 3ds2rib Version " VERSION "\n"
"Copyright (C) 1996-2001 by J.E. Hoffmann <je-h@gmx.net>\n"
"All rights reserved.\n"
"\n"
"Syntax: 3ds2rib [options] filename [options]\n"
"\n"
"Options:\n"
"  -h/--help               This help\n"
"  -o/--output <filename>  Write output to <filename> instead of stdout\n"
"  -c/--camera <name>      Use camera <name> for rendering\n"
"  -a/--all                Render all frames\n"
"  -f/--frame #            Render frame #\n"
"  -d/--downcuts #         Render # frames of animation\n"
"\n"
);
  exit(1);
}


typedef enum _Flags {
  LIB3DS2RIB_ALL  =0x0001
} Flags;

static const char* filename=0;
static const char* output=0;
static Lib3dsDword flags=0;
static float frame=0.0f;
static const char* camera=0;
static int downcuts=0;


static void
parse_args(int argc, char **argv)
{
  int i;
  
  for (i=1; i<argc; ++i) {
    if (argv[i][0]=='-') {
      if ((strcmp(argv[i],"-h")==0) || (strcmp(argv[i],"--help")==0)) {
        help();
      }
      else
      if ((strcmp(argv[i],"-o")==0) || (strcmp(argv[i],"--output")==0)) {
        ++i;
        if (output || (i>=argc)) {
          help();
        }
        output=argv[i];
      }
      else
      if ((strcmp(argv[i],"-a")==0) || (strcmp(argv[i],"--all")==0)) {
        flags|=LIB3DS2RIB_ALL;
      }
      else
      if ((strcmp(argv[i],"-f")==0) || (strcmp(argv[i],"--frame")==0)) {
        ++i;
        if (i>=argc) {
          help();
        }
        frame=(Lib3dsFloat)atof(argv[i]);
      }
      else
      if ((strcmp(argv[i],"-c")==0) || (strcmp(argv[i],"--camera")==0)) {
        ++i;
        if (i>=argc) {
          help();
        }
        camera=argv[i];
      }
      else
      if ((strcmp(argv[i],"-d")==0) || (strcmp(argv[i],"--downcuts")==0)) {
        ++i;
        if (i>=argc) {
          help();
        }
        downcuts=atoi(argv[i]);
      }
      else {
        help();
      }
    }
    else {
      if (filename) {
        help();
      }
      filename=argv[i];
    }
  }
  if (!filename) {
    help();
  }
}


static void
rib_concat_transform(FILE *o, Lib3dsMatrix m)
{
  int i,j;
  
  fprintf(o, "ConcatTransform [");
  for (i=0; i<4; ++i) {
    for (j=0; j<4; ++j) {
      fprintf(o, "%f ", m[i][j]);
    }
  }
  fprintf(o, "]\n");
}


static void
rib_camera(FILE *o, Lib3dsFile *f, Lib3dsMatrix M)
{
  Lib3dsNode *c;
  Lib3dsNode *t;
  const char *name=camera;

  ASSERT(f);
  if (!name) {
    if (f->cameras) {
      name=f->cameras->name;
    }
  }
  if (!name) {
    fprintf(stderr, "***ERROR*** No camera found!\n");
    return;
  }
  c=lib3ds_file_node_by_name(f, name, LIB3DS_CAMERA_NODE);
  t=lib3ds_file_node_by_name(f, name, LIB3DS_TARGET_NODE);
  if (!c || !t) {
    fprintf(stderr, "***ERROR*** Invalid camera/target!\n");
    return;
  }

  lib3ds_matrix_camera(M, c->data.camera.pos, t->data.target.pos, c->data.camera.roll);
  rib_concat_transform(o, M);
  fprintf(o, "Camera \"perspective\" \"float fov\" [%f]\n", c->data.camera.fov);
}


static void
rib_lights(FILE *o, Lib3dsFile *f, Lib3dsMatrix M)
{
  Lib3dsLight *light;
  Lib3dsNode *l;
  Lib3dsNode *s;
  int i=1;

  for (light=f->lights; light; light=light->next) {
    l=lib3ds_file_node_by_name(f, light->name, LIB3DS_LIGHT_NODE);
    s=lib3ds_file_node_by_name(f, light->name, LIB3DS_SPOT_NODE);
    if (!l) {
      fprintf(stderr, "***ERROR*** Invalid light!\n");
      continue;
    }
    if (s) {
      Lib3dsVector pos,spot;
      lib3ds_vector_copy(pos, l->data.light.pos);
      lib3ds_vector_copy(spot, s->data.spot.pos);
      fprintf(o,
        "LightSource "
        "\"spot\" "
        "\"point from\" [%f %f %f] "
        "\"point to\" [%f %f %f]\n",
        pos[0],
        pos[1],
        pos[2],
        spot[0],
        spot[1],
        spot[2]
      );
    }
    else {
      Lib3dsVector pos;
      lib3ds_vector_copy(pos, l->data.light.pos);
      fprintf(o,
        "LightSource "
        "\"point\" "
        "\"from\" [%f %f %f]\n",
        pos[0],
        pos[1],
        pos[2]
      );
    }
    ++i;
  }
}



static void
flushpp(FILE *o, float *P, float *N, float *uv, int nfaces) {
    int i;
    if (!nfaces) return;
    fprintf(o, "Shape \"trianglemesh\" \"integer indices\" [ ");
    for (i = 0; i < nfaces; ++i) fprintf(o, "\t%d %d %d\n", 3*i, 3*i+1, 3*i+2);
    fprintf(o, "] \"point P\" [");
    for (i = 0; i < 9*nfaces; ++i) fprintf(o, "%f%c", P[i], (i%3)==0? '\n':' ');
    fprintf(o, "] \"normal N\" [");
    for (i = 0; i < 9*nfaces; ++i) fprintf(o, "%f%c", N[i], (i%3)==0? '\n':' ');
    fprintf(o, "] \"float st\" [");
    for (i = 0; i < 6*nfaces; ++i) fprintf(o, "%f%c", uv[i], (i%2)==0? '\n':' ');
    fprintf(o, "]\n\n");
}


struct Tex {
    char *map, *mask, *texname;
    int isFloat;
};
/* if more than that, this linear search will be a killer anyway */
#define MAX_TEXTURES 1024
static struct Tex textures[MAX_TEXTURES];
static int nTextures;

char *GetTexture(FILE *o, const char *map, const char *mask, int floattex)
{
    int i, len;
    char *cp, *tp;

    if (!map || !mask) return NULL;
    for (i = 0; i < nTextures; ++i) {
	if (!strcmp(map, textures[i].map) && !strcmp(mask, textures[i].mask) &&
	    textures[i].isFloat == floattex)
	    return textures[i].texname;
    }

    len = strlen(map);
    textures[nTextures].map = malloc(len + 20);
    strcpy(textures[nTextures].map, "textures/");
    tp = rindex(textures[nTextures].map, '/')+1;
    for (i = 0; i < len; ++i)
	tp[i] = tolower(map[i]);
    tp[len] = '\0';
    cp = rindex(tp, '.');
    if (cp) {
	++cp;
	*cp++ = 't';
	*cp++ = 'i';
	*cp++ = 'f';
	*cp = '\0';
    }
    
    textures[nTextures].mask = strdup(mask);
    textures[nTextures].texname = strdup(map);
    if (rindex(textures[nTextures].texname, '.'))
	*rindex(textures[nTextures].texname, '.') = '\0';

    fprintf(o, "Texture \"%s\" \"%s\" \"imagemap\" \"string filename\" [\"%s\"]\n",
	    textures[nTextures].texname, 
	    floattex ? "float" : "color", textures[nTextures].map);
    return textures[nTextures++].texname;
}

static void
create_node(Lib3dsFile *f, Lib3dsNode *node, FILE *o)
{
  Lib3dsMesh *mesh;
  
  if ((node->type==LIB3DS_OBJECT_NODE) && (strcmp(node->name,"$$$DUMMY")!=0)) {
    mesh=lib3ds_file_mesh_by_name(f, node->name);
    ASSERT(mesh);
    if (mesh) {
      Lib3dsObjectData *d=&node->data.object;

      fprintf(o, "\n\n##\n## Object: %s\n##\n", node->name);
      fprintf(o, "AttributeBegin\n");

      {
        Lib3dsMatrix N,M,X;
        lib3ds_matrix_copy(N, node->matrix);
        lib3ds_matrix_translate_xyz(N, -d->pivot[0], -d->pivot[1], -d->pivot[2]);
        lib3ds_matrix_copy(M, mesh->matrix);
        lib3ds_matrix_inv(M);
        lib3ds_matrix_mul(X,N,M);
        rib_concat_transform(o, X);
      }
      {
        unsigned p;
	int i, j;
        Lib3dsVector *normalL=malloc(3*sizeof(Lib3dsVector)*mesh->faces);
        lib3ds_mesh_calculate_normals(mesh, normalL);
	Lib3dsMaterial *lastmat = NULL;

	int nalloc = 256, nfaces = 0;
	float *P = (float *)malloc(nalloc*9*sizeof(float));
	float *N = (float *)malloc(nalloc*9*sizeof(float));
	float *uv = (float *)malloc(nalloc*6*sizeof(float));

        for (p=0; p<mesh->faces; ++p) {
          Lib3dsFace *face=&mesh->faceL[p];
          Lib3dsMaterial *mat=lib3ds_file_material_by_name(f, face->material);
          if (mat && mat != lastmat) {
	      char *Kdmap = NULL, *Omap = NULL;

	      flushpp(o, P, N, uv, nfaces);
	      nfaces = 0;
	      lastmat = mat;

	      if (mat->texture1_map.name[0]) {
		  Kdmap = GetTexture(o, mat->texture1_map.name, mat->texture1_mask.name, 0);
		  char *scale = malloc(strlen(Kdmap)+20);
		  strcpy(scale, Kdmap);
		  strcat(scale, "-scale");
		  fprintf(o, "Texture \"%s\" \"color\" \"scale\" \"texture tex1\" \"%s\" "
			  "\"color tex2\" [.8 .8 .8]\n", scale, Kdmap);
		  Kdmap = scale;
	      }
	      if (mat->opacity_map.name[0])
		  Omap = GetTexture(o, mat->opacity_map.name, mat->opacity_mask.name, 0);

	      char *bumpScale = NULL;
	      if (mat->bump_map.name[0]) {
		  char *Bmap = GetTexture(o, mat->bump_map.name, mat->bump_mask.name, 1);
		  bumpScale = malloc(strlen(Bmap)+20);
		  strcpy(bumpScale, Bmap);
		  strcat(bumpScale, "-scale");
		  fprintf(o, "Texture \"%s\" \"float\" \"scale\" \"texture tex1\" \"%s\" "
			  "\"float tex2\" [.05]\n", bumpScale, Bmap);
	      }

	      fprintf(o, "Material \"uber\" ");
	      if (Kdmap) fprintf(o, "\"texture Kd\" \"%s\" ", Kdmap);
	      else       fprintf(o, "\"color Kd\" [%f %f %f] ",
				 mat->diffuse[0], mat->diffuse[1], mat->diffuse[2]);
	      if (Omap) fprintf(o, "\"texture opacity\" \"%s\" ", Omap);
	      else      fprintf(o, "\"color opacity\" [%f %f %f] ",
				1.f - mat->transparency,
				1.f - mat->transparency,
				1.f - mat->transparency);
	      fprintf(o, "\"color Ks\" [%f %f %f] "
		      "\"float roughness\" [%f] ",
		      mat->specular[0]*mat->shin_strength,
		      mat->specular[1]*mat->shin_strength,
		      mat->specular[2]*mat->shin_strength,
		      mat->shininess
		      );
//CO	      dumptex(o, "tex2", &mat->texture2_map); 
//CO	      dumptex(o, "tex2mask", &mat->texture2_mask);
//CO	      dumptex(o, "specular", &mat->specular_map); 
//CO	      dumptex(o, "specularmask", &mat->specular_mask);
//CO	      dumptex(o, "shininess", &mat->shininess_map); 
//CO	      dumptex(o, "shininessmask", &mat->shininess_mask);
//CO	      dumptex(o, "reflection", &mat->reflection_map); 
//CO	      dumptex(o, "reflectionmask", &mat->reflection_mask);
	      fprintf(o, "\n");

	      if (bumpScale)
		  fprintf(o, "\"float bumpmap\" \"%s\"\n", bumpScale);
          }
	  if (nfaces+1 == nalloc) {
	      nalloc *= 2;
	      P = (float *)realloc(P, nalloc*9*sizeof(float));
	      N = (float *)realloc(N, nalloc*9*sizeof(float));
	      uv = (float *)realloc(uv, nalloc*6*sizeof(float));
	  }

	  for (i = 0; i < 3; ++i) {
	      for (j = 0; j < 3; ++j) {
		  P[9*nfaces+3*i+j] = mesh->pointL[face->points[i]].pos[j];
		  N[9*nfaces+3*i+j] = normalL[3*p+i][j];
		  if (j != 2 && mesh->texelL) 
		      uv[6*nfaces+2*i+j] = mesh->texelL[face->points[i]][j];
	      }
	      if (!mesh->texelL) {
		  uv[6*nfaces+0] = 0;
		  uv[6*nfaces+1] = 0;
		  uv[6*nfaces+2] = 0;
		  uv[6*nfaces+3] = 1;
		  uv[6*nfaces+4] = 1;
		  uv[6*nfaces+5] = 0;
	      }
	  }
	  ++nfaces;
	}

	flushpp(o, P, N, uv, nfaces);
	free(P);
	free(N);
	free(uv);
        free(normalL);
      }
      nTextures = 0;
      fprintf(o, "AttributeEnd\n");
    }
  }
  {
    Lib3dsNode *n;
    for (n=node->childs; n; n=n->next) {
      create_node(f,n,o);
    }
  }
}


static void
create_rib(Lib3dsFile *f, FILE *o, int current)
{
  Lib3dsMatrix M;
  fprintf(o, "Film \"image\" \"integer xresolution\" [400] \"integer yresolution\" [400]\n");
  fprintf(o, "Scale 1 -1 1\n");
  fprintf(o, "Rotate 90 1 0 0\n");

  rib_camera(o,f,M);
  
  fprintf(o, "WorldBegin\n");

  fprintf(o, "AttributeBegin\n");
  fprintf(o, "  CoordSysTransform \"camera\"\n");
  fprintf(o, "#  LightSource \"point\" \"float intensity\" [ 300 ]\n");
  fprintf(o, "AttributeEnd\n\n");

  rib_lights(o,f,M);
  {
    Lib3dsNode *n;
    for (n=f->nodes; n; n=n->next) {
      create_node(f,n,o);
    }
  }
  fprintf(o, "WorldEnd\n");
}


int
main(int argc, char **argv)
{
  Lib3dsFile *f=0;
  FILE *o;

  parse_args(argc, argv);
  f=lib3ds_file_load(filename);
  if (!f) {
    fprintf(stderr, "***ERROR***\nLoading file %s failed\n", filename);
    exit(1);
  }

  if (output) {
    o=fopen(output, "w+");
    if (!o) {
      fprintf(stderr, "***ERROR***\nCan't open %s for writing\n", output);
      exit(1);
    }
  }
  else {
    o=stdout;
  }

  if (flags&LIB3DS2RIB_ALL) {
    int i;
    for (i=f->segment_from; i<=f->segment_to; ++i) {
      lib3ds_file_eval(f,1.0f*i);
      create_rib(f,o,i);
    }
  }
  else
  if (downcuts) {
    int i;
    int delta=f->segment_to-f->segment_from;
    for (i=0; i<downcuts; ++i) {
      float frame=f->segment_from+1.0f*i*delta/(downcuts-1);
      lib3ds_file_eval(f, frame);
      create_rib(f,o,i);
    }
  }
  else {  
    lib3ds_file_eval(f,frame);
    create_rib(f,o, (int)frame);
  }

  if (o!=stdout) {
    fclose(o);
  }
  lib3ds_file_free(f);
  return(0);
}







