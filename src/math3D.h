/*
 *  math.h
 *  EditMesh
 *
 *  Created by roberto on Sun Sep 14 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */

// Constants
#ifndef __math3D__
#define __math3D__

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef pi
#define pi 3.14159265358979323
#endif

// structures
typedef struct {
	float	x;
	float	y;
}float2D, **float2DHandle;
typedef struct {
	float	x;
	float	y;
	float	z;
}float3D, *float3DPtr,**float3DHandle;
typedef struct {
	double	x;
	double	y;
}double2D, **double2DHandle;
typedef struct {
	double	x;
	double	y;
	double	z;
}double3D, *double3DPtr,**double3DHandle;

typedef struct {
	float2D	a;
	float2D	b;
	float2D	c;
}tria2D, **tria2DHandle;
typedef struct {
	float3D	a;
	float3D	b;
	float3D	c;
}tria3D, **tria3DHandle;
typedef struct {
	float3D	siz;
	float3D	ide;
}rect3D, *rect3DPtr;

typedef struct {
	int	a;
	int	b;
}int2D, *int2DPtr;
typedef struct{
	int	a;
	int	b;
	int	c;
}int3D,*int3DPtr;
typedef struct{
	int	a;
	int	b;
	int	c;
	int	d;
}int4D,*int4DPtr;

//#define	ABS(x)		(((x)>0)?(x):(-(x)))
#define DIRECTION(x)	sca3D( (x) , 1/norm3D((x)) )

// math Prototypes
double	norm3D(float3D);
double	norm23D(float3D a);
double	dot3D(float3D, float3D);
float3D cross3D(float3D, float3D);
float3D mul3D(float3D, float3D);
float3D sub3D(float3D, float3D);
float3D add3D(float3D, float3D);
float3D sca3D(float3D, float);
int equals3D(float3D a, float3D b);

double2D	add2D(double2D a,double2D b);
double2D	sub2D(double2D a,double2D b);
double2D 	sca2D(double2D a, double);
double 		dot2D(double2D a,double2D b);
double 		cross2D(double2D a,double2D b);
double		norm2D(double2D a);
double		norm22D(double2D a);
double2D 	triCenter2D(double2D a,double2D b,double2D c);

char*		cvector_new(int size);
void		cvector_dispose(char *v);
char**		cmatrix_new(int row,int col);
void		cmatrix_dispose(char **m,int row);
int*		ivector_new(int size);
void		ivector_dispose(int *v);
int**		imatrix_new(int row,int col);
void		imatrix_dispose(int **m,int row);
float*		fvector_new(int size);
void		fvector_dispose(float *v);
void		fvector_save(int n, float v[]);
float**		fmatrix_new(int row,int col);
void		fmatrix_dispose(float **m,int row);
double*		dvector_new(int size);
void		dvector_dispose(double *v);
double**	dmatrix_new(int row,int col);
void		dmatrix_dispose(double **m,int row);
int2D*		i2vector_new(int size);
void		i2vector_dispose(int2D *v);
float2D*	f2vector_new(int size);
void		f2vector_dispose(float2D *v);
double2D*	d2vector_new(int size);
void		d2vector_dispose(double2D *v);
float3D*	f3vector_new(int size);
void		f3vector_dispose(float3D *v);

void	mMat(float *a,float *b,float *c);
void	iMat(float *a,float *b);
float3D	matfloat3D(float *b,float3D a);
float3D	matsolvefloat3D(float *m, float3D v);

int	gausselim(double **A,double *B,int r,int c);
int	backsub(double **A,double *B,int r,int c);

float3D	triPlane(float3D, float3D, float3D);
float3D	triCenter(float3D, float3D, float3D);
double	triArea(float3D, float3D, float3D);
void	triBase(float3D *t, float3D *b);
int	isPtInTriangle(float3D p, float3D *t);

float3D	trianglecoord_to_float3D(float3D *T, float2D c);
int	intersect_VectorTriangle( float3D V, float3D *T, float2D *C, double close_enough);
int	intersect_segments(float3D a[], float3D b[], float3D I[],float c[], double close_enough);
int	intersect_segment_triangle(float3D pe[], float3D ptri[], float3D Itri[], float3D Iseg[],float Ctri[], float Cseg[],int side[2],double close_enough);

double areaOfTriangleInSquare2D(double2D *tr, double2D *sq);
double areaOfTriangleInTriangle2D(double2D *tr0, double2D *tr1);
double triangleArea2D(double2D *tr);
int pointInTriangle2D(double2D pt, double2D *tr,double err);
double areaIntersectingTriangles2D(double2D *tr0, double2D *tr1,double err);
double edgeIntersectEdge2D(int2D a, int2D b, double2D *p,double err);
int intersect_vector_segment(float2D v[2], float2D s[2], float *d);

float3D stereographic(float3D p);
float3D sinusoidal(float3D p);
float3D mercator(float3D p);
void dynstereographic(float *mat, float *src, float *dst);
void dynsinusoidal(float *mat, float *src, float *dst);
void dynmercator(float *mat, float *src, float *dst);

void	fquicksort(int ir, int *arr, float *varr);
void	dquicksort(int ir, int *arr, double *varr);

// SVD
float	pythag(float a, float b);
float SIGN(float a, float b);
void svdcmp(float **a, int m, int n, float w[], float **v);

#endif
