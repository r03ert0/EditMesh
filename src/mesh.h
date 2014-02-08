// Includes
#ifndef __mesh3d__
#define __mesh3d__

#include "math3D.h"

// transform functions
#define SMOOTH(pt,p,f)		(pt) = s3D( (pt) , esc3D( r3D( (p) , (pt) ) , (f) ) )
#define INFLATE(pt,n,f)		(pt) = s3D( (pt) , esc3D( (n) , (f) ) )
#define ANTICONTRACT(pt,p,f)	(pt) = s3D( (pt) , esc3D( (p) , (f) ) )

#define MIN(x,y)    (((x)<(y))?(x):(y))

// Constants
#define SIZESTACK	64

#define kfrontalplan	0
#define klateralplan	1
#define khorizontalplan	2

// Structures
typedef struct {
	int		n;
	int		t[SIZESTACK];
}NTriRec, *NtriPtr;
typedef struct {
	int		n;
	int		e[SIZESTACK];
}NEdgeRec, *NEdgePtr;
typedef struct {
	int			np;
	int			nt;
	float3D		*p;
	int3D		*t;
	
	float3D		center;
	rect3D		bbox;
	char		name[256];

	long		xcontainer;	// extra data container
}MeshRec, *MeshPtr;

// MeshPoints are points defined over the surface of a triangle
typedef struct {
	int		t;		// associated mesh triangle
	float		x;		// x coordinate on the triangle
	float		y;		// y coordinate on the triangle
}MeshPointRec, *MeshPointPtr;
typedef struct {
	int		np;		// number of points
	int		ne;		// number of edges
	MeshPointPtr	p;		// mesh points
	int2DPtr	e;		// edge relations
	float3DPtr	c;		// texture of the meshpoint
	
	char		name[256];	
}MeshCurveRec, *MeshCurvePtr;

// MeshEdgePoints are points defined over the edges of a triangle
typedef struct {
	int	ta;		// 1st associated mesh triangle (ta<0 => tb is a vertex index)
	int	tb;		// 2nd associated mesh triangle
	int	sa;		// side on the 1st triangle
	int	sb;		// side on the 2nd triangle
	float	t;		// edge-point parameter over the 1st triangle (1-t for the 2nd)
}MeshEdgePointRec, *MeshEdgePointPtr;
typedef struct {
	int			np;		// number of points
	int			ne;		// number of edges
	MeshEdgePointPtr	p;		// mesh edge-points
	int2DPtr		e;		// edge relations
	float3DPtr		c;		// texture of the mesh edge-point
	
	char			name[256];	
}MeshEdgeCurveRec, *MeshEdgeCurvePtr;

typedef struct {
	int	nd;		// number of data
	int	sd;		// size in bytes of each datum
	short	id;		// data ID
	char	name[256];	// data name
	
	char	*data;		// data handle
	long	xcontainer;	// extra data container
}ContainerRec,*ContainerPtr;

// IO Prototypes
void msh_parseRawText(MeshRec *mesh, char *path);
int msh_packRawText(MeshRec *m, char **data);

void msh_parseOFFText(MeshRec *mesh, char *path);
int msh_packOFFText(MeshRec *m, char **data);

void msh_parse3DMFText(MeshRec *mesh, char *data);
int msh_pack3DMF(MeshRec *m, char **data);

void msh_parseVerticesDataBin(float **vdat, int *ddim, int np, char *data, int size);
int msh_packVerticesData(float *vdat, int dim, int np, char **data);

int msh_importFSMeshData(MeshRec *mesh, char *path);
int msh_packFSMeshData(MeshRec *mesh, char **data);
void msh_importFSTextureData(float **dat, int np, char *path);
int msh_exportFSTextureData(float *d, int np, char *path);
int msh_importFSMeshAnnotation(float **dat, int np, char *path);

void msh_importBVMeshData(MeshRec *mesh, char *path);

int msh_importVRMLMeshData(MeshRec *mesh, char *path);
int msh_packVRML(MeshRec *m, char **data);

int msh_importPlyMeshData(MeshRec *mesh, char *path);
int msh_packPly(MeshRec *m, float3D *C, char **data);

void exitmesh(void);

// Prototypes
bool msh_new(MeshPtr *m, int np, int nt);
void msh_dispose(MeshPtr mesh);
bool msh_copy(MeshPtr src, MeshPtr *dst);

bool msh_newMeshCurve(MeshCurveRec *MC,int np,int ne);
void msh_disposeMeshCurve(MeshCurveRec MC);
bool msh_newMeshEdgeCurve(MeshEdgeCurveRec *MEC,int np,int ne);
void msh_disposeMeshEdgeCurve(MeshEdgeCurveRec MEC);

void msh_setContainer(ContainerRec *c, int nd, int sd, short id, unsigned char *name);
void msh_addData(MeshPtr m, ContainerRec *c);
void msh_deleteData(MeshPtr m, ContainerRec c);
bool msh_findData(MeshPtr m, ContainerRec *c);

float3D*		msh_getPointsPtr(MeshPtr m);
int3D*			msh_getTrianglesPtr(MeshPtr m);
NEdgeRec*		msh_getNeighborEdgesPtr(MeshPtr m);
NTriRec*		msh_getNeighborTrianglesPtr(MeshPtr m);
int2D*			msh_getEdgesPtr(MeshPtr m);
float3D*		msh_getTexturePtr(MeshPtr m);
char*			msh_getLabelPtr(MeshPtr m);
MeshCurveRec*		msh_getMeshCurvePtr(MeshPtr m);
MeshEdgeCurveRec*	msh_getMeshEdgeCurvePtr(MeshPtr m);
char*			msh_getDataPtr(MeshPtr m,ContainerRec c);

void msh_deleteTexturePtr(MeshPtr m);
void msh_deleteNeighborTrianglesPtr(MeshPtr m);

void update(MeshPtr m);
void msh_setNeighborTriangles(MeshPtr m);
void msh_setNeighborEdges(MeshPtr m);
void msh_setEdges(MeshPtr m);
void msh_storeEdge(int2D *E, int *n, int2D e);

void msh_get_triangle_neighbortriangles(MeshPtr m, int t, int *tn, int *s);
void msh_get_triangle_neighboredges(MeshPtr m, int t, int *e);

void msh_linbcgDN(int n, int en, double D[], double N[], int2D E[], double b[], double x[], int itol, double tol, int itmax, int *iter, double *err);
double msh_snrmDN(int n, double sx[], int itol);
void msh_asolveDN(int n, double D[], double B[], double X[]);
void msh_atimesDN(int n, int en, double D[], double N[], int2D E[], double x[], double r[]);
void msh_SORDN(MeshPtr m, double D[], double N[], double B[], double X[], double tol, int itmax, int *iter, double *err);
void msh_drawDN(int n, int en, double D[], double N[], int2D E[],int dec);

void msh_laplace_fe(MeshPtr m, double *U,int nDir, int *iDir, double *Dir,
                    int nNeu, int *iNeu, double *Neu,
                    int nMFC, int2D *eMFC, double *gMFC,double *hMFC);

int msh_psort(MeshPtr mesh, int m, int *pstack);
int msh_esort(MeshPtr m, int indx, int *pstack, int *estack);

void plane(int meshpoint, MeshPtr mesh, float3D *np, float3D *pp);
float3D tTriPlane(int nt, MeshPtr m);
float3D tTriCenter(int nt, MeshPtr m);
float tTriArea(MeshPtr m,int nt);
float3D tTrinorm3Dal(MeshPtr m, int nt);

void msh_LocalBaseFromPoint(MeshPtr m, int i, float3D *iv, float3D *jv, float3D *kv);
void msh_LocalBaseFromTriangle(MeshPtr m, int t, float3D *iv, float3D *jv, float3D *kv);
float3D msh_VectorFromAngle(float theta, float3D iv, float3D jv);
float3D msh_proj3DtoBase(float3D p, float3D iv, float3D jv, float3D kv);

float determinant(float3D a, float3D b, float3D c);
float volume(float3D *p, int3D *t, int nt);
float surface_area(float3D *p, int3D *t, int nt);

bool trianglesort(int np, MeshPtr mesh, int *pstack);
float3D neighbor3Dmean(MeshPtr mesh, int m);

float3D	msh_meshpoint_to_float3D(MeshPtr m, MeshPointRec mp);
MeshPointRec msh_float3D_to_meshpoint(MeshPtr m, float3D p);
MeshPointRec msh_float3D_to_meshpoint_n(MeshPtr m, float3D x);
MeshPointRec msh_float3D_to_meshpoint_proj(MeshPtr m, float3D x);

float3D	msh_meshedgepoint_to_float3D(MeshPtr m, MeshEdgePointRec mp);
MeshEdgePointRec msh_float3D_to_meshedgepoint(MeshPtr m, float3D p);

#endif
