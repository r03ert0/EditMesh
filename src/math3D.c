/*
 *  math.c
 *  EditMesh
 *
 *  Created by roberto on Sun Sep 14 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */

#include "math3D.h"

#define sqrnorm3D(a)	((a).x*(a).x + (a).y*(a).y + (a).z*(a).z)

#pragma mark -
// vectorial functions
double norm3D(float3D a)
{
    double	xx;
    
    xx= sqrt(pow(a.x,2)+pow(a.y,2)+pow(a.z,2));
    return(xx);
}
double norm23D(float3D a)
{
    double	xx;
    
    xx= pow(a.x,2)+pow(a.y,2)+pow(a.z,2);
    return(xx);
}
double dot3D(float3D a, float3D b)
{
    double xx;
    
    xx=a.x*b.x + a.y*b.y + a.z*b.z;
    return(xx);
}

float3D cross3D(float3D a, float3D b)
{
    float3D	xx;
    
    xx.x = a.y*b.z - a.z*b.y;
    xx.y = -b.z*a.x + b.x*a.z; // SIGNS WERE INVERTED BEFORE!!!
    xx.z = a.x*b.y - a.y*b.x;
    return(xx);
}
float3D mul3D(float3D a, float3D b)
{
    float3D xx;
    
    xx.x=a.x*b.x;
    xx.y=a.y*b.y;
    xx.z=a.z*b.z;
    return(xx);
}
float3D sub3D(float3D a, float3D b)
{
    float3D xx;
    
    xx.x=a.x-b.x;
    xx.y=a.y-b.y;
    xx.z=a.z-b.z;
    return(xx);
}
float3D add3D(float3D a, float3D b)
{
    float3D xx;
    
    xx.x=a.x+b.x;
    xx.y=a.y+b.y;
    xx.z=a.z+b.z;
    return(xx);
}
float3D sca3D(float3D a, float b)
{
    float3D	xx;
    
    xx.x = a.x*b;
    xx.y = a.y*b;
    xx.z = a.z*b;
    return(xx);
}
int equals3D(float3D a, float3D b)
{
    int	x,y,z;
    
    x=(a.x==b.x);
    y=(a.y==b.y);
    z=(a.z==b.z);
    return ((a.x==b.x) && (a.y==b.y) && (a.z==b.z));
}
#pragma mark -
double2D add2D(double2D a,double2D b)
{
    return (double2D){(a).x+(b).x,(a).y+(b).y};
}
double2D sub2D(double2D a,double2D b)
{
    return (double2D){(a).x-(b).x,(a).y-(b).y};
}
double2D sca2D(double2D a, double k)
{
    return (double2D){k*a.x,k*a.y};
}
double	dot2D(double2D a,double2D b)
{
    return ((a).x*(b).x+(a).y*(b).y);
}
double	cross2D(double2D a,double2D b)
{
    return(a.x*b.y-b.x*a.y);
}
double norm2D(double2D a)
{
    return sqrt((a).x*(a).x + (a).y*(a).y);
}
double norm22D(double2D a)
{
    return (a).x*(a).x + (a).y*(a).y;
}
double2D triCenter2D(double2D a,double2D b,double2D c)
{
    return (double2D){((a).x+(b).x+(c).x)/3.0,((a).y+(b).y+(c).y)/3.0};
}
#pragma mark -
char* cvector_new(int size)
{
    char *v;

    v=(char *)calloc(size,sizeof(char));
    if (!v){ printf("allocation failure in cvector_new()\n");}
    return v;
}
void cvector_dispose(char *v)
{
    free(v);
}
char **cmatrix_new(int row,int col)
{
    int		i;
    char	**m;

    m=(char **) calloc(row,sizeof(char*));

    for(i=0;i<row;i++)
    {
        m[i]=(char *) calloc(col,sizeof(char));
    }
    return m;
}
void cmatrix_dispose(char **m,int row)
{
    int i;

    for(i=0;i<row;i++) free(m[i]);
    free(m);
}
int* ivector_new(int size)
{
    int *v;

    v=(int *)calloc(size,sizeof(int));
    if (!v){ printf("allocation failure in ivector_new\n");}
    
    return v;
}
void ivector_dispose(int *v)
{
    free((char *)v);
}
int **imatrix_new(int row,int col)
{
    int		i;
    int		**m;

    m=(int **) calloc(row,sizeof(int*));

    for(i=0;i<row;i++)
    {
            m[i]=(int *) calloc(col,sizeof(int));
    }
    return m;
}
void imatrix_dispose(int **m,int row)
{
    int i;

    for(i=0;i<row;i++) free((char *)m[i]);
    free((char *)m);
}
float* fvector_new(int size)
{
    float *v;

    v=(float *)calloc(size,sizeof(float));
    if (!v){printf("allocation failure in fvector_new\n");}
    
    return v;
}
void fvector_dispose(float *v)
{
    free((char *)v);
}
float **fmatrix_new(int row,int col)
{
    int		i;
    float	**m;

    m=(float **) calloc(row,sizeof(float*));

    for(i=0;i<row;i++)
    {
        m[i]=(float *) calloc(col,sizeof(float));
    }
    return m;
}
void fmatrix_dispose(float **m,int row)
{
    int i;

    for(i=0;i<row;i++) free((char *)m[i]);
    free((char *)m);
}
double* dvector_new(int size)
{
    double *v;

    v=(double *)calloc(size,sizeof(double));
    if (!v){printf("allocation failure in dvector_new\n");}
    
    return v;
}
void dvector_dispose(double *v)
{
    free((char *)v);
}
double **dmatrix_new(int row,int col)
{
    int		i;
    double	**m;

    m=(double **) calloc(row,sizeof(double*));

    for(i=0;i<row;i++)
    {
            m[i]=(double *) calloc(col,sizeof(double));
    }
    return m;
}
void dmatrix_dispose(double **m,int row)
{
    int i;

    for(i=0;i<row;i++) free((char *)m[i]);
    free((char *)m);
}
int2D* i2vector_new(int size)
{
    int2D *v;

    v=(int2D *)calloc(size,sizeof(int2D));
    if (!v){printf("allocation failure in i2vector_new"); }
    
    return v;
}
void i2vector_dispose(int2D *v)
{
    free((char *)v);
}
float2D* f2vector_new(int size)
{
    float2D *v;

    v=(float2D *)calloc(size,sizeof(float2D));
    if (!v){ printf("allocation failure in f2vector_new");}
    
    return v;
}
void f2vector_dispose(float2D *v)
{
    free((char *)v);
}
double2D* d2vector_new(int size)
{
    double2D *v;

    v=(double2D *)calloc(size,sizeof(double2D));
    if (!v){ printf("allocation failure in d2vector_new");}
    
    return v;
}
void d2vector_dispose(double2D *v)
{
    free((char *)v);
}
float3D* f3vector_new(int size)
{
    float3D *v;

    v=(float3D *)calloc(size,sizeof(float3D));
    if (!v){printf("allocation failure in f3vector_new"); }
    
    return v;
}
void f3vector_dispose(float3D *v)
{
    free((char *)v);
}

#pragma mark -
// 3d matrix functions
void mMat(float *a,float *b,float *c)
{
    float	aux[9];
    int		i;

    aux[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
    aux[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
    aux[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
    
    aux[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
    aux[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
    aux[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
    
    aux[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
    aux[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
    aux[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
    
    for(i=0;i<9;i++) c[i]=aux[i];
}
void iMat(float *a,float *b)
// The inverse of the matrix a stored in matrix b
// Input: matrix a[9]
// Output: matrix b[9]
{
    float	det;
    
    det=a[0]*(a[4]*a[8]-a[5]*a[7])
            +a[1]*(a[5]*a[6]-a[3]*a[8])
            +a[2]*(a[3]*a[7]-a[4]*a[6]);
            
    if(det==0)
            b=a;
    else
    {
        b[0]=(a[4]*a[8]-a[5]*a[7])/det;
        b[1]=(a[2]*a[7]-a[1]*a[8])/det;
        b[2]=(a[1]*a[5]-a[2]*a[4])/det;
        
        b[3]=(a[5]*a[6]-a[3]*a[8])/det;
        b[4]=(a[0]*a[8]-a[2]*a[6])/det;
        b[5]=(a[2]*a[3]-a[0]*a[5])/det;
        
        b[6]=(a[3]*a[7]-a[4]*a[6])/det;
        b[7]=(a[1]*a[6]-a[0]*a[7])/det;
        b[8]=(a[0]*a[4]-a[1]*a[3])/det;
    }
}

float3D matfloat3D(float *b,float3D a)
{
    float3D	c;

    c.x = a.x*b[0] + a.y*b[1] + a.z*b[2];
    c.y = a.x*b[3] + a.y*b[4] + a.z*b[5];
    c.z = a.x*b[6] + a.y*b[7] + a.z*b[8];
    
    return c;
}
float3D	matsolvefloat3D(float *m, float3D v)
{
    float3D	p;

    p.x =	(v.x*(m[4]*m[8]-m[5]*m[7]) +
                        v.y*(m[2]*m[7]-m[1]*m[8]) +
                        v.z*(m[1]*m[5]-m[2]*m[4]))
                    /
                    (m[0]*(m[4]*m[8]-m[5]*m[7]) +
                        m[3]*(m[2]*m[7]-m[1]*m[8]) +
                        m[6]*(m[1]*m[5]-m[2]*m[4]));

    p.y =	(v.x*(m[3]*m[8]-m[5]*m[6]) +
                        v.y*(m[2]*m[6]-m[0]*m[8]) +
                        v.z*(m[0]*m[5]-m[2]*m[3]))
                    /
                    (m[1]*(m[3]*m[8]-m[5]*m[6]) +
                        m[4]*(m[2]*m[6]-m[0]*m[8]) +
                        m[7]*(m[0]*m[5]-m[2]*m[3]));

    p.z =	(v.x*(m[3]*m[7]-m[4]*m[6]) +
                        v.y*(m[1]*m[6]-m[0]*m[7]) +
                        v.z*(m[0]*m[4]-m[1]*m[3]))
                    /
                    (m[2]*(m[3]*m[7]-m[4]*m[6]) +
                        m[5]*(m[1]*m[6]-m[0]*m[7]) +
                        m[8]*(m[0]*m[4]-m[1]*m[3]));
    
    return p;
}
int mth_gausselim(double **A,double *B,int r,int c)
{
    int		i,j,k;
    int		max;
    double	temp;
    
    #define SWAP(a,b) {temp=a;a=b;b=temp;}
    for(i=0;i<r;i++)
    {
            // pivot
            max=i;
            for(j=i;j<r;j++)
                    if(A[j][i]>A[max][i])
                            max=j;
            SWAP(B[max],B[i])
            for(j=i;j<c;j++)
                    SWAP(A[max][j],A[i][j])
            
            // norm3Dalize
            B[i]*=1/A[i][i];
            for(k=c-1;k>=i;k--)
                    A[i][k]*=1/A[i][i];
            
            // eliminate
            for(j=i+1;j<r;j++)
            {
                    B[j]-=B[i]*A[j][i];
                    for(k=c-1;k>=0;k--)
                            A[j][k]-=A[i][k]*A[j][i];
            }
    }
    #undef SWAP
    
    return true;
}
int mth_backsub(double **A,double *B,int r,int c)
{
    int	i,j;

    for(i=r-1;i>=0;i--)
            for(j=c-1;j>i;j--)
                    B[i]-=B[j]*A[i][j];
    return true;
}
#pragma mark -
// triangle functions
float3D triCenter(float3D a, float3D b, float3D c)
{
    float3D	p;
    
    p.x=(a.x+b.x+c.x)/3.0;
    p.y=(a.y+b.y+c.y)/3.0;
    p.z=(a.z+b.z+c.z)/3.0;
    
    return(p);
}
double triArea(float3D a, float3D b, float3D c)
{
    double3D	aa=(double3D){a.x-c.x,a.y-c.y,a.z-c.z};
    double3D	bb=(double3D){b.x-c.x,b.y-c.y,b.z-c.z};
    double		p=sqrt(	pow(aa.y*bb.z-aa.z*bb.y,2)+
						pow(aa.z*bb.x-aa.x*bb.z,2)+
						pow(aa.x*bb.y-aa.y*bb.x,2))/2.0;

    return(p);
}
float3D triPlane(float3D a, float3D b, float3D c)
{
    float3D	p,cero={0,0,0};
    
    p= cross3D( sub3D(b,a), sub3D(c,a) );
    if(norm3D(p))
		return(sca3D(p,1.0/norm3D(p)));
    return(cero);		
}
void triBase(float3D *t, float3D *b)
// t[3]: three 3D points of a triangle
// b[3]: an associated 3D base with  origin at t[0]
//	(b[0],b[1],b[2])=(tangent,perpendicular,norm3Dal)
{
    b[0]=sub3D(t[1],t[0]);
    b[0]=sca3D(b[0],1/norm3D(b[0]));
    
    b[2]=cross3D(b[0],sub3D(t[2],t[0]));
    b[2]=sca3D(b[2],1/norm3D(b[2]));
    
    b[1]=cross3D(b[2],b[0]);
}
int mth_isPtInTriangle(float3D p, float3D *t)
// true if the point is inside or in the triangle
{
    float3D	a,b,c;
    float	TEST;
    
    a=cross3D(sub3D(t[1],t[0]),sub3D(p,t[0]));
    b=cross3D(sub3D(t[2],t[1]),sub3D(p,t[1]));
    c=cross3D(sub3D(t[0],t[2]),sub3D(p,t[2]));
    
    TEST=dot3D(a,b);
    TEST=dot3D(b,c);
    TEST=dot3D(c,a);
    
    if(dot3D(a,b)<0 || dot3D(b,c)<0 || dot3D(c,a)<0)
            return false;
    else
            return true;
}
float3D	trianglecoord_to_float3D(float3D *T, float2D c)
// returns a 3D point from triangular coordinates p=T0+s*(T1-T0)+t*(T2-T0)
// input:  triangle T[3], coordinates c
// output: point p
{
    float3D p;
    
    p=add3D(sca3D(sub3D(T[1],T[0]),c.x),sca3D(sub3D(T[2],T[0]),c.y));
    p=add3D(T[0],p);
    
    return p;
}
#define SMALL_NUM  0.0000001 // anything that avoids division overflow
// intersect_VectorTriangle(): intersect a vector with a 3D triangle
//    Input:  a vector V, and a triangle T
//    Output: *C = the coordinates (s,t) of the intersection point (when it exists)
//					I=T0+s*(T1-T0)+t*(T2-T0)
//    Return: -1 = triangle is degenerate (a segment or point)
//             0 = disjoint (no intersect)
//             1 = intersect in unique point I1
//             2 = are in the same plane
// code from:http://geometryalgorithms.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()
int intersect_VectorTriangle( float3D V, float3D *T, float2D *C, double close_enough)
{
    float3D		u, v, n;             // triangle vectors
    float3D		dir, w0, w;          // ray vectors
    float3D		I;					// the intersection point
    double		r, a, b;             // params to calc ray-plane intersect
    double		uu, uv, vv, wu, wv, D;
    double		s, t;
    
   	*C=(float2D){0,0};

    // get triangle edge vectors and plane norm3Dal
    u = sub3D(T[1],T[0]);
    v = sub3D(T[2],T[0]);
    n = cross3D(u,v);      		// cross product
    if (norm3D(n) == 0)    		// triangle is degenerate
        return -1;          	// do not deal with this case

    dir = V;            		// ray direction vector
    w0 =  sca3D(T[0],-1);
    a = dot3D(n,w0);
    b = dot3D(n,dir);
    if (fabs(b) < SMALL_NUM)	// ray is parallel to triangle plane
    {
        if (a == 0)       		// ray lies in triangle plane
        { /*printf("intersect_VectorTriangle: ray lies in triangle plane\n");*/  return 2;}
        else   					// ray disjoint from plane
        { /*printf("intersect_VectorTriangle: ray disjoint from plane\n");*/  return 0;}
    }

    // get intersect point of ray with triangle plane
    r=-a/b;
    if(r<0.0)        			// ray goes away from triangle
   		return 0;				// => no intersect
    // for a segment, also test if (r > 1.0) => no intersect

    I = sca3D(dir,r); 			// intersect point of ray and plane

    // is I inside T?
    uu = dot3D(u,u);
    uv = dot3D(u,v);
    vv = dot3D(v,v);
    w = sub3D(I,T[0]);
    wu = dot3D(w,u);
    wv = dot3D(w,v);
    D = uv * uv - uu * vv;

    // get and test parametric coords
    s = (uv * wv - vv * wu) / D;
	if(fabs(s)<close_enough) s=0;
	if(fabs(1-s)<close_enough) s=1;
 		
    t = (uv * wu - uu * wv) / D;
	if(fabs(t)<close_enough) t=0;
	if(fabs(1-t)<close_enough) t=1;
   
    *C=(float2D){s,t};

    if (s < 0.0 || t < 0.0 || (s + t) > 1.0)  // I is outside T
    	return 0;

    return 1;                      // I is in T
}
int intersect_segments(float3D *L0, float3D *L1, float3D *I,float *cord, double close_enough)
// find the closer_to_intersection point between two segments
// input: a[2],b[2] the 2 points of the first and second segments
// output: I[2], the [0,1] coordinates of the CTI point over the first and second segments
{
    float3D   u = sub3D(L0[1],L0[0]);
    float3D   v = sub3D(L1[1],L1[0]);
    float3D   w = sub3D(L0[0],L1[0]);
    float    a = dot3D(u,u);        // always >= 0
    float    b = dot3D(u,v);
    float    c = dot3D(v,v);        // always >= 0
    float    d = dot3D(u,w);
    float    e = dot3D(v,w);
    float    D = a*c - b*b;         // always >= 0
    float    s,t;
    int		 intersection=0;

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM)
    {
    	// the lines are almost parallel
    	// return the limit in the first segment
    	// contained in between the second segment limits
    	
    	t=dot3D(v,sub3D(L0[0],L1[0]))/c;
    	if(t>=0 && t<=1)
    	{	s=0;	intersection=1;}
    	else
    	{
	    	t=dot3D(v,sub3D(L0[1],L1[0]))/c;
	    	if(t>=0 && t<=1)
	    	{	s=1;	intersection=1;}
    	}
    }
    else
    {
        s = (b*e - c*d)/D;
        if(fabs(s)<close_enough) s=0;
        if(fabs(1-s)<close_enough) s=1;
        
        t = (a*e - b*d)/D;
         if(fabs(t)<close_enough) t=0;
        if(fabs(1-t)<close_enough) t=1;
       
        intersection=(s>=0&&s<=1)&&(t>=0&&t<=1);
    }
    
    I[0]=add3D(sca3D(L0[0],1-s),sca3D(L0[1],s));
    I[1]=add3D(sca3D(L1[0],1-t),sca3D(L1[1],t));
    cord[0]=s;
    cord[1]=t;
    
    return intersection;
}
int intersect_segment_triangle(	float3D pe[], float3D ptri[],
								float3D Itri[], float3D Iseg[],
								float Ctri[], float Cseg[],
								int side[2],double close_enough)
// find the one or two points of closer_to_intersection between a segment and a triangle
// input: pe[2] the 2 points of the segment, ptri[3] the 3 points (a,b,c) of a triangle
// output: Itri[2] the points of intersection over the triangle, Iseg[2] % over the segment
//         Ctri[2] the coordinates of the intersection over the triangle, Cseg[2] % over the segment
//		   points over the triangle, side[2] the #side of the intersections
{
	float3D	sc[2];
	float3D	a[2],b[2],c[2];
	float	ca[2],cb[2],cc[2];
	int		ia,ib,ic;
	int		n=0,sides=0;
	
	sc[0]=ptri[2];sc[1]=ptri[0];

	ia=intersect_segments(ptri,  pe,a,ca, close_enough);
	ib=intersect_segments(ptri+1,pe,b,cb, close_enough);
	ic=intersect_segments(sc,    pe,c,cc, close_enough);
	
	if(ia)
	{	if(ca[0]==1 && (cb[0]==0||cb[0]==1))
		{	Itri[n]=a[0];Iseg[n]=a[1];	Ctri[n]=1;	  Cseg[n]=ca[1];	side[n]=1;	n++;sides-=1;}
		else if(ca[0])
		{	Itri[n]=a[0];Iseg[n]=a[1];	Ctri[n]=ca[0];Cseg[n]=ca[1];	side[n]=1;	n++;}
	}
	if(ib)
	{	if(cb[0]==1 && (cc[0]==0||cc[0]==1))
		{	Itri[n]=b[0];Iseg[n]=b[1];	Ctri[n]=1;	  Cseg[n]=cb[1];	side[n]=2;	n++;sides-=1;}
		else if(cb[0])
		{	Itri[n]=b[0];Iseg[n]=b[1];	Ctri[n]=cb[0];Cseg[n]=cb[1];	side[n]=2;	n++;}
	}
	if(ic)
	{	if(cc[0]==1 && (ca[0]==0||ca[0]==1))
		{	Itri[n]=c[0];Iseg[n]=c[1];	Ctri[n]=1;	  Cseg[n]=cc[1];	side[n]=3;	n++;sides-=1;}
		else if(cc[0])
		{	Itri[n]=c[0];Iseg[n]=c[1];	Ctri[n]=cc[0];Cseg[n]=cc[1];	side[n]=3;	n++;}
	}

	if(n>2)
            printf("ia:%d, ib:%d, ic:%d: intersect_segment_triangle\n",ia,ib,ic);
	
	return ia+ib+ic+sides;
}
#pragma mark -
double areaOfTriangleInSquare2D(double2D *tr, double2D *sq)
{
	double2D tr0[3],tr1[3];
	double  a0,a1;
	
	tr0[0]=sq[0]; tr0[1]=sq[1]; tr0[2]=sq[3];
	tr1[0]=sq[1]; tr1[1]=sq[2]; tr1[2]=sq[3];

	a0=areaOfTriangleInTriangle2D(tr,tr0);
	a1=areaOfTriangleInTriangle2D(tr,tr1);
	
	return a0+a1;
}
double areaOfTriangleInTriangle2D(double2D *tr0, double2D *tr1)
{
	int cs0,cs1,i;
	double  area,err=10e-6;

	cs0=cs1=0;
	for(i=0;i<3;i++)
	{
		cs0+=pointInTriangle2D(tr0[i],tr1,err);
		cs1+=pointInTriangle2D(tr1[i],tr0,err);
	}

	if(cs0==3)			area=triangleArea2D(tr0);
	else if(cs1==3)		area=triangleArea2D(tr1);
	else				area=areaIntersectingTriangles2D(tr0,tr1,err);
	
	//printf("%f\n",area);
	return area;
}
double triangleArea2D(double2D *tr)
{
	double2D c,d;
	
	c=sub2D(tr[1],tr[0]);
	d=sub2D(tr[2],tr[0]);
	
	return fabs(cross2D(c,d)*0.5);
}
int pointInTriangle2D(double2D pt, double2D *tr,double err)
{
	double2D   c,d;
	int		in[3];
	int		i;
	
	for(i=0;i<3;i++)
	{
		c=sub2D(pt,tr[i]);
		d=sub2D(tr[(i+1)%3],tr[i]);
		in[i]=(cross2D(d,c)>-err)?1:0;
	}

	return in[0]*in[1]*in[2];
}
double areaIntersectingTriangles2D(double2D *tr0, double2D *tr1,double err)
{
	int			i,j,iter;
	double		t,area;
	int2D		e[20],e0[20],e1[20];
	int			ne,ne0,ne1;
	double2D	p[20];
	int			np;
	int			sorted[20],n;
	int			A,B,C,D, exit;
	
	//printf("//double\ttr0[3][2]={{%f,%f},{%f,%f},{%f,%f}},",tr0[0].x,tr0[0].y,tr0[1].x,tr0[1].y,tr0[2].x,tr0[2].y);
	//printf("//double\ttr1[3][2]={{%f,%f},{%f,%f},{%f,%f}};\n",tr1[0].x,tr1[0].y,tr1[1].x,tr1[1].y,tr1[2].x,tr1[2].y);
	// 1. fill points and edges vectors
	for(i=0;i<3;i++)
	{   p[i]=tr0[i];
		e0[i]=(int2D){i,(i+1)%3};
	}
	for(i=0;i<3;i++)
	{   p[3+i]=tr1[i];
		e1[i]=(int2D){3+i,3+(i+1)%3};
	}
	np=6;
	ne0=ne1=3;
	
	// 2. intersect edges
	for(i=0;i<ne0;i++) // edges in tr0
	for(j=0;j<ne1;j++) // edges in tr1
	{
		t=edgeIntersectEdge2D(e0[i],e1[j],p,err);
		//printf("%f %i(%i,%i)-%i(%i,%i)\n",t,i,e0[i].a,e0[i].b,3+j,e1[j].a,e1[j].b);
		if(t>0&&t<1)
		{
			p[np].x=(1-t)*p[e0[i].a].x+t*p[e0[i].b].x;
			p[np].y=(1-t)*p[e0[i].a].y+t*p[e0[i].b].y;
	
			e0[ne0]=(int2D){np,e0[i].b};
			e0[i].b=np;
			ne0++;

			e1[ne1]=(int2D){np,e1[j].b};
			e1[j].b=np;
			ne1++;
			
			np++;
		}
	}

	// 3. eliminate edges where 1 or 2 points are not shared
	ne=0;
	for(i=0;i<ne0;i++)
		e[ne++]=e0[i];
	for(i=0;i<ne1;i++)
		e[ne++]=e1[i];
	i=0;
	do
	{
		A=pointInTriangle2D(p[e[i].a],tr0,err);
		B=pointInTriangle2D(p[e[i].b],tr0,err);
		C=pointInTriangle2D(p[e[i].a],tr1,err);
		D=pointInTriangle2D(p[e[i].b],tr1,err);
		//printf("(%i) %i%i\t(%i)%i%i\n",e[i].a,A,C,e[i].b,B,D);
		if(A==0||B==0||C==0||D==0)
		{
			ne--;
			for(j=i;j<ne;j++)
				e[j]=e[j+1];
		}
		else
			i++;
	}
	while(i<ne);
	if(ne==0) return 0;

	// 4. sort the remaining points
	ne0=ne; // variable recycling...
	n=iter=exit=0;
	sorted[n++]=e[0].a;
	sorted[n++]=e[0].b;
	do
	{
		iter++;
		for(i=1;i<ne;i++)
		{
			if(e[i].a==sorted[n-1])
			{
				sorted[n++]=e[i].b;
				if(e[i].b==sorted[0]) exit=1;
				ne--;
				for(j=i;j<ne;j++)
					e[j]=e[j+1];
				break;
			}
			else
			if(e[i].b==sorted[n-1])
			{
				sorted[n++]=e[i].a;
				if(e[i].a==sorted[0]) exit=1;
				ne--;
				for(j=i;j<ne;j++)
					e[j]=e[j+1];
				break;
			}
		}
	}
	while(exit==0 && iter<ne0);
	n--;
	
	if(iter>=ne0) printf("can't close the intersection polygon...\n");
		
	// 5.compute area
	area=0;
	for(i=0;i<n;i++)
		area+=(p[sorted[i]].x+p[sorted[(i+1)%n]].x)*(p[sorted[i]].y-p[sorted[(i+1)%n]].y);
	area=0.5*fabs(area);
	
	return area;
}
double edgeIntersectEdge2D(int2D a, int2D b, double2D *p,double err)
{
	double m0,n0,m1,n1;
	double x,y,s,t;
	double   a10,b10;
	
	a10=(p[a.b].x-p[a.a].x);
	b10=(p[b.b].x-p[b.a].x);

	if(a10==0)
	{
		if(b10==0)
		{
			s=-1000;
			t=-1000;
		}
		else
		{
			m1=(p[b.b].y-p[b.a].y)/b10;
			n1=p[b.a].y-m1*p[b.a].x;
			x=p[a.a].x;
			y=m1*x+n1;
			s=(y-p[a.a].y)/(p[a.b].y-p[a.a].y);
			t=(x-p[b.a].x)/b10;
		}
	}
	else
	{
		m0=(p[a.b].y-p[a.a].y)/a10;
		n0=p[a.a].y-m0*p[a.a].x;
		if(fabs(b10)<err)
		{
			x=p[b.a].x;
			y=m0*x+n0;
			s=(x-p[a.a].x)/a10;
			t=(y-p[b.a].y)/(p[b.b].y-p[b.a].y);
		}
		else
		{
			m1=(p[b.b].y-p[b.a].y)/b10;
			n1=p[b.a].y-m1*p[b.a].x;
			
			if(m0==m1)
			{
				s=-1000;
				t=-1000;
			}
			else
			{
				x=(n1-n0)/(m0-m1);
				y=(m0*n1-m1*n0)/(m0-m1);
				s=(x-p[a.a].x)/a10;
				t=(x-p[b.a].x)/b10;
			}
		}
	}
	if(fabs(s)<err)		s=0;
	if(fabs(s-1)<err)   s=1;
	if(fabs(t)<err)		t=0;
	if(fabs(t-1)<err)   t=1;
	
	if(t<=0 || t>=1) s=t;
	return s;
}
int intersect_vector_segment(float2D v[2], float2D s[2], float *dist)
{
	double  a,b,c,d,e,f;
	double  det,alpha,t;
	
	a=v[1].x;   b=s[0].x-s[1].x;
	c=v[1].y;   d=s[0].y-s[1].y;
	e=s[0].x-v[0].x;
	f=s[0].y-v[0].y;
	
	det=a*d-b*c;
	if(fabs(det)<10e-3)
		return 0;		// 0: lines are parallel
	
	alpha=(e*d-b*f)/det;
	t=(a*f-e*c)/det;
	
	if(t<0 || t>1)
		return 0;		// 0: intersection not in the segment
	*dist=alpha;
	return 1;
}
#pragma mark -
float3D stereographic(float3D p)
// stereographic projection
// input: cartesian vector 0->p
// output: stereographic projection, (-z)-pole open.
{
    float3D	xp,sp;
    float	x,y,z,n;
    float	h,v,delta;
        
    n=norm3D(p);
    xp = sca3D(p,1/n);

    x = acos(xp.x);
    y = acos(xp.y);
    z = acos(xp.z);

    if(z*z<0.000001)	delta=1;
    else		delta = cos(y)/sin(z);
    if(delta<-1) delta=-1;
    if(delta>1) delta=1;

    h = z*sin(acos(delta));
    v = z*delta;
    if(x>pi/2.0)	sp = (float3D){-h,v,n};
    else		sp = (float3D){h,v,n};
    
    return sp;
}
float3D sinusoidal(float3D p)
{
    float	y,a,n;
    float	h,v;
    float3D	xp,sp;
    
    n=norm3D(p);
    xp=sca3D(p,1/n);

    y = acos(xp.y);
    a = atan2(xp.x,xp.z);

    h = sin(y)*a;
    v = pi/2-y;
    
    sp = (float3D){-h,v,n};
    
    return sp;
}
float3D mercator(float3D p)
{
    float	y,a,n;
    float	h,v;
    float3D	xp,sp;
    
    n=norm3D(p);
    xp=sca3D(p,1/n);

    y = acos(xp.z);
    a = atan2(-xp.x,xp.y);

    //h = a;
    h= a-pi/4.0; if(h<-pi) h=h+2*pi;
    v = y-pi/2;
    
    sp = (float3D){h,v,n};
    
    return sp;
}
/*
void dynstereographic(float *mat, float *src, float *dst)
{
    float	xp[3],sp[3];
    float	x,y,z,n;
    float	h,v,delta;
        
    n=sqrt(src[0]*src[0]+src[1]*src[1]+src[2]*src[2]);
    xp[0] =src[0]/n; xp[1] =src[1]/n; xp[2] =src[2]/n;

    x = acos(xp[0]*mat[0]+xp[1]*mat[1]+xp[2]*mat[2]);
    y = acos(xp[0]*mat[3]+xp[1]*mat[4]+xp[2]*mat[5]);
    z = acos(xp[0]*mat[6]+xp[1]*mat[7]+xp[2]*mat[8]);

    if(z*z<0.000001)	delta=1;
    else		delta = cos(y)/sin(z);
    if(delta<-1) delta=-1;
    if(delta>1) delta=1;

    h = z*sin(acos(delta));
    v = z*delta;
    if(x>pi/2.0)
	{	sp[0]=-h; sp[1]=v; sp[2]=n;}
    else
	{   sp[0]=h; sp[1]=v; sp[2]=n;}
    
	dst[0]=mat[0]*sp[0]+mat[3]*sp[1]+mat[6]*sp[2];
	dst[1]=mat[1]*sp[0]+mat[4]*sp[1]+mat[7]*sp[2];
	dst[2]=mat[2]*sp[0]+mat[5]*sp[1]+mat[8]*sp[2];
}
void dynsinusoidal(float *mat, float *src, float *dst)
{
    float	xp[3],sp[3];
    float	x,y,z,a,n;
    float	h,v;
    
    n=sqrt(src[0]*src[0]+src[1]*src[1]+src[2]*src[2]);
    xp[0] =src[0]/n; xp[1] =src[1]/n; xp[2] =src[2]/n;

    x = xp[0]*mat[0]+xp[1]*mat[1]+xp[2]*mat[2];
    y = xp[0]*mat[3]+xp[1]*mat[4]+xp[2]*mat[5];
    z = xp[0]*mat[6]+xp[1]*mat[7]+xp[2]*mat[8];

    y = acos(y);
    a = atan2(x,z);

    h = sin(y)*a;
    v = pi/2-y;
    
    sp[0]=h; sp[1]=v; sp[2]=n;
	
	dst[0]=mat[0]*sp[0]+mat[3]*sp[1]+mat[6]*sp[2];
	dst[1]=mat[1]*sp[0]+mat[4]*sp[1]+mat[7]*sp[2];
	dst[2]=mat[2]*sp[0]+mat[5]*sp[1]+mat[8]*sp[2];
}
void dynmercator(float *mat, float *src, float *dst)
{
    float	xp[3],sp[3];
    float	x,y,z,a,n;
    float	h,v;
    
    n=sqrt(src[0]*src[0]+src[1]*src[1]+src[2]*src[2]);
    xp[0] =src[0]/n; xp[1] =src[1]/n; xp[2] =src[2]/n;

    x = xp[0]*mat[0]+xp[1]*mat[1]+xp[2]*mat[2];
    y = xp[0]*mat[3]+xp[1]*mat[4]+xp[2]*mat[5];
    z = xp[0]*mat[6]+xp[1]*mat[7]+xp[2]*mat[8];

    y = acos(y);
    a = atan2(x,z);

    h = a;
    v = pi/2-y;
    
    sp[0]=h; sp[1]=v; sp[2]=n;
	
	dst[0]=mat[0]*sp[0]+mat[3]*sp[1]+mat[6]*sp[2];
	dst[1]=mat[1]*sp[0]+mat[4]*sp[1]+mat[7]*sp[2];
	dst[2]=mat[2]*sp[0]+mat[5]*sp[1]+mat[8]*sp[2];
}
*/
#pragma mark -
#define SWAP(a,b) temp=(a); (a)=(b); (b)=temp;
#define MM 7
#define NSTACK 50
void fquicksort(int ir, int *arr, float *varr)
// slightly modified Quicksort from Numerical Recipes in C
// input: arr=index vector, varr= float value vector, ir=vector size
// output: arr sorted from min to max
{
	int		i;
	long	j,k,l=1;
	int		jstack=0, istack[NSTACK];
	int		a,temp;

	for(;;)
	{
		if(ir-l<MM)
		{
			for(j=l+1;j<=ir;j++)
			{
				a=arr[j];
				for(i=j-1;i>=1;i--)
				{
					if( varr[arr[i]]<=varr[a] )
						break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if(jstack ==0)
				break;
			ir=istack[jstack--];
			l=istack[jstack--];
		}
		else
		{
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			if( varr[arr[l+1]] > varr[arr[ir]] )
			{
				SWAP(arr[l+1],arr[ir])
			}
			if( varr[arr[l]] > varr[arr[ir]] )
			{
				SWAP(arr[l],arr[ir])
			}
			if( varr[arr[l+1]] > varr[arr[l]] )
			{
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for(;;)
			{
				do i++; while( varr[arr[i]] < varr[a] );
				do j--; while( varr[arr[j]] > varr[a] );
				if(j<i)
					break;
				SWAP(arr[i],arr[j]);
			}
			arr[l]=arr[j];
			arr[j]=a;
			jstack+=2;
			
			if(jstack>NSTACK)
			{
				printf("fquicksort\n");
				//NSTACK too small in sort
				return;
			}
			if(ir-i+1>=j-1)
			{
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			}
			else
			{
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}
void dquicksort(int ir, int *arr, double *varr)
// slightly modified Quicksort from Numerical Recipes in C
// input: arr=index vector, varr= double value vector, ir=vector size
// output: arr sorted from min to max
{
	int		i;
	long	j,k,l=1;
	int		jstack=0, istack[NSTACK];
	int		a,temp;

	for(;;)
	{
		if(ir-l<MM)
		{
			for(j=l+1;j<=ir;j++)
			{
				a=arr[j];
				for(i=j-1;i>=1;i--)
				{
					if( varr[arr[i]]<=varr[a] )
						break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if(jstack ==0)
				break;
			ir=istack[jstack--];
			l=istack[jstack--];
		}
		else
		{
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			if( varr[arr[l+1]] > varr[arr[ir]] )
			{
				SWAP(arr[l+1],arr[ir])
			}
			if( varr[arr[l]] > varr[arr[ir]] )
			{
				SWAP(arr[l],arr[ir])
			}
			if( varr[arr[l+1]] > varr[arr[l]] )
			{
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for(;;)
			{
				do i++; while( varr[arr[i]] < varr[a] );
				do j--; while( varr[arr[j]] > varr[a] );
				if(j<i)
					break;
				SWAP(arr[i],arr[j]);
			}
			arr[l]=arr[j];
			arr[j]=a;
			jstack+=2;
			
			if(jstack>NSTACK)
			{
				printf("dquicksort\n");
				//NSTACK too small in sort
				return;
			}
			if(ir-i+1>=j-1)
			{
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			}
			else
			{
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}
// SVD
#pragma mark -
#pragma mark [   SV Decomposition   ]
#pragma mark -
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define FMAX(a,b)	((a>b)?a:b)
#define IMIN(a,b)	((a<b)?a:b)

float pythag(float a, float b)
//Computes sqrt(a^2 + b^2) without destructive underflow or overflow.
{
    float	absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
float SIGN(float a, float b)
{
    float t;
    if(b>=0.0)	t=fabs(a);
    else	t=-fabs(a);
    return t;
}

void svdcmp(float **a, int m, int n, float w[], float **v)
// Given a matrix a[1..m][1..n], this routine computes its singular value decomposition,
// A = UWV'. The matrix U replaces a on output. The diagonal matrix of singular values W is output
// as a vector w[1..n]. The matrix V (not the transpose V' ) is output as v[1..n][1..n].
{
   int		flag,i,its,j,jj,k,l,nm;
    float	anorm3D,c,f,g,h,s,scale,x,y,z,*rv1;
    
    rv1=fvector_new(n+1);
    g=scale=anorm3D=0.0; //Householder reduction to bidiagonal form.
    for (i=1;i<=n;i++)
    {
        l=i+1;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i <= m)
        {
            for (k=i;k<=m;k++) scale += fabs(a[k][i]);
            if (scale)
            {
                for (k=i;k<=m;k++)
                {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f=a[i][i];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;
                for (j=l;j<=n;j++)
                {
                    for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
                    f=s/h;
                    for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
                }
                for (k=i;k<=m;k++) a[k][i] *= scale;
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if (i <= m && i != n)
        {
            for (k=l;k<=n;k++) scale += fabs(a[i][k]);
            if (scale)
            {
                for (k=l;k<=n;k++)
                {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f=a[i][l];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][l]=f-g;
                for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
                for (j=l;j<=m;j++)
                {
                    for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
                    for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
                }
                for (k=l;k<=n;k++) a[i][k] *= scale;
            }
        }
        anorm3D=FMAX(anorm3D,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n;i>=1;i--)
    { // Accumulation of right-hand transformations.
        if (i < n)
        {
            if (g)
            {
                for (j=l;j<=n;j++) //Double division to avoid possible underflow.
                    v[j][i]=(a[i][j]/a[i][l])/g;
                for (j=l;j<=n;j++)
                {
                    for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
                    for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
                }
            }
            for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
        }
        v[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
    for (i=IMIN(m,n);i>=1;i--)
    { // Accumulation of left-hand transformations.
        l=i+1;
        g=w[i];
        for (j=l;j<=n;j++) a[i][j]=0.0;
        if (g)
        {
            g=1.0/g;
            for (j=l;j<=n;j++)
            {
                for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
                f=(s/a[i][i])*g;
                for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
            }
            for (j=i;j<=m;j++) a[j][i] *= g;
        } else for (j=i;j<=m;j++) a[j][i]=0.0;
        ++a[i][i];
    }
    for (k=n;k>=1;k--)
    { // Diagonalization of the bidiagonal form: Loop over
        // singular values, and over allowed iterations.
        for (its=1;its<=30;its++)
        {
            flag=1;
            for (l=k;l>=1;l--)
            { //Test for splitting.
                nm=l-1; // Note that rv1[1] is always zero.
                if ((float)(fabs(rv1[l])+anorm3D) == anorm3D)
                {
                    flag=0;
                    break;
                }
                if ((float)(fabs(w[nm])+anorm3D) == anorm3D) break;
            }
            if (flag)
            {
                c=0.0; // Cancellation of rv1[l], if l > 1.
                s=1.0;
                for (i=l;i<=k;i++)
                {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if ((float)(fabs(f)+anorm3D) == anorm3D) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=1;j<=m;j++)
                    {
                        y=a[j][nm];
                        z=a[j][i];
                        a[j][nm]=y*c+z*s;
                        a[j][i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l == k)
            { // Convergence.
                if (z < 0.0)
                { // Singular value is made nonnegative.
                    w[k] = -z;
                    for (j=1;j<=n;j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 30) printf("no convergence in 30 svdcmp iterations");
            x=w[l]; // Shift from bottom 2-by-2 minor.
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0; // Next QR transformation:
            for (j=l;j<=nm;j++)
            {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g = g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=1;jj<=n;jj++)
                {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z; // Rotation can be arbitrary if z = 0.
                if (z)
                {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=1;jj<=m;jj++)
                {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
    fvector_dispose(rv1);
}