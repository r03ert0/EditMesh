/*
 *  editmesh.h
 *  EditMesh
 *
 *  Created by roberto on Sat Sep 13 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */

// Includes
#ifndef __editmesh__
#define __editmesh__

#include "mesh.h"

#define AREAZERO(x)	((T[(x)].a==T[(x)].b) || (T[(x)].b==T[(x)].c) || (T[(x)].c==T[(x)].a))

#define kfire		1
#define kice		2
#define kspectrum	3

#ifndef pi
#define pi 3.14159265358
#endif

// Prototypes
void	em_infoLength(MeshPtr mesh);
void	em_infoNeighbour(MeshPtr mesh);
void	em_infoArea(MeshPtr mesh);
int	em_infoTopology(MeshPtr mesh);

//RGBColor em_colourFromColourmap(int cm_number, float cm_index);

//---------------
// mesh functions
//---------------
void	em_simplify(MeshPtr m);
void	em_simplifyPoint(MeshPtr mesh);
bool	em_deletePoint(MeshPtr mesh,int m);
void	em_take(MeshPtr mesh,int m, int *pstack);

void	em_setFisionThreshold(float fisionthrs);
float	em_getFisionThreshold(void);
bool	em_fision(MeshPtr snake);

void	em_setFusionParam(float fusionthrs, float fusionP);
void	em_getFusionParam(float *fusionthrs, float *fusionP);
void	em_fusion(MeshPtr mesh);
void	em_eraseInsideTriangleFromint2D(int3D t, int sp, char *us, int *np, MeshPtr m);

void	em_splitPoint(MeshPtr m, int np);

void	em_flipEdges(MeshPtr m, int triangle);
void	em_printPointInformation(MeshPtr m, int point);

void	em_setSmooth(float smoothfactor);
float	em_getSmooth(void);
void	em_smooth(MeshPtr mesh);

void	em_setInflate(float scale);
float	em_getInflate(void);
void	em_inflate(MeshPtr mesh);

void	em_setScale(float scale);
float	em_getScale(void);
void	em_scale(MeshPtr m);

void	em_centre(MeshPtr mesh);

void	em_translate(MeshPtr m,float3D t);

void	em_flipTriangles(MeshPtr m);

void	em_pickPoint(MeshPtr mesh);
void	em_drawSelected(MeshPtr mesh);
int 	em_getSelectedPoint(void);

float	em_laplaceFromAngle(float angle);
float	em_angleFromLaplace(float laplace);
void	em_sphereFromTxtr(MeshPtr m, float3D *C);
float3D em_getPointFromSphericalCoordinate(float3D c);

void	em_fusionMeshToMeshEdgeCurve(MeshPtr *mm,MeshPtr m,MeshEdgeCurveRec MEC);
void	em_setFusionMeshToMeshEdgeCurve(MeshPtr *mm,MeshPtr m,MeshEdgeCurveRec MEC,int *npoints,int *ntrian, bool isStoring);

void	em_alignSVD(MeshPtr m);

void	em_changecoordinates(MeshPtr mesh,float *mr);

void    em_adaptMesh(MeshPtr m);

//---------------
//texture functions
//---------------
void	em_pictureToTexture(MeshPtr mesh);
//void	em_openPictureToTextureData(MeshPtr m,FSSpec spec);
//void	em_getTriangleMeanColor(tria2D T,RGBColor *rgb);
void	em_markTriangle(tria2D T);

bool em_textureDepth(MeshPtr m, float *C);
bool em_textureMeanCurvature(MeshPtr m, float *C);
bool em_textureIntegratedMeanCurvature(MeshPtr m, float *C, int iter);
bool em_textureArea(MeshPtr m, float *C);
bool em_textureGaussCurvature(MeshPtr m);

void em_textureSulcalHierarchy(MeshPtr m, float *C);
void em_textureGyralHierarchy(MeshPtr m, float *C);
void em_textureCountSulci(MeshPtr m, float *C);
void em_textureCountGyri(MeshPtr m, float *C);
void fillsulcus(MeshPtr m, float *C, int i, float sul, int *size);
void fillgyrus(MeshPtr m, float *C, int i, float gyr, int *size);

void em_textureSulcalDepth2(MeshPtr m, MeshPtr h, float *C);

void	em_textureLaplaceFD(MeshPtr mesh);
void	em_textureLaplaceFE(MeshPtr m);
void	em_textureCoordinatesFE(MeshPtr m);
void	em_textureHeat(MeshPtr mesh);
void	em_textureConformalMapping(MeshPtr mesh);
void	em_textureTestConformal(MeshPtr mesh);
void	em_textureConformalFE(MeshPtr m);
float3D	em_textureEvalDeriv(float3D aa, float3D bb, float3D cc,int i,double alf);

void em_textureSmoothVerticesData(MeshPtr m, float *C);

void	em_changeCoordinates(MeshPtr mesh);

/*
void	em_drawTalairachBox(x_Rect r, int2D p, float a, int2D *P);
bool	em_clickTalairachBox(Rect *r, int2D *p, float *a, int2D *P);
void	em_talairachCoordinates(MeshPtr mesh);
void	em_setTalairachTransform(float3D *M);
void	em_mapTalairach(MeshPtr m, long *ntxt, char *txtH, int *nlab, float3DPtr *labH);
int	em_mapOpenTalairachCoordinates(long *txtSize, char *txtH);
int	em_mapParseTalairachCoordinates(long txtSize, char *txtH, int *nlab, float3DPtr *labH);
int	em_mapProjectTalairachOnMesh(MeshPtr m,int nlab,float3DPtr labH);
float3D em_getint2DFromTalairachCoordinate(float3D tal);
*/

void em_textureDistortion(MeshPtr m, MeshPtr m0, float *C);
#endif