/* MyMeshView */
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#import <Cocoa/Cocoa.h>
#import "Trackball.h"

#include "editmesh.h"
#include "colourmap.h"

@interface MyMeshView : NSOpenGLView
{
	NSMutableDictionary	*settings;
	
	int		ntris;
	int		nverts;
        
	int			*tris;
	float		*verts;
	GLfloat		*vertsdata;
	GLfloat		*vertscolour;
	GLfloat		*vertsparam;
	GLfloat		*vertsnormals;
	
	GLubyte		*texture;
	NSSize		textureSize;
	GLuint		textureId;
	BOOL		hasTexture;
	
	int			proj;
	GLfloat		*vertsproj;
	
	int			shading;
	
	// cursor
	int			selected;
	float3D		cp[8];
	int3D		ct[12];
	float3D		cc[12];
	float		m[16],invm[16];
	
	Trackball	*m_trackball;
	float		m_rotation[4];	// The main rotation
	float		m_tbRot[4];	// The trackball rotation
	
	float		zoom;		// Zoom = exp(zoom)
	int			cmap_index;
	int			vertsnormals_computed;
	GLuint		shaderTexture;
}
-(void)setSettings:(NSMutableDictionary*)newSettings;

-(float*)verts;
-(float*)vertsproj;
-(int)nverts;
-(void)setVertices:(float *)vert number:(int)n;
-(int*)tris;
-(int)ntris;
-(void)setTriangles:(int *)trian number:(int)n;
-(void)setColourmap:(int)i;
-(void)setMinMaxData;
-(void)setVerticesData:(float *)vertdata;
-(float*)verticesData;
-(void)setVerticesColour:(float *)vertcolour;
-(float*)verticesColour;
-(void)setParam:(float *)param;
-(void)setTexture:(char *)image size:(NSSize)s;
-(void)setTextureActive:(BOOL)isActive;

-(void)setStandardRotation:(int)view;
-(void)addRotation:(float)value toAxis:(int)axis;

-(void)rotateBy:(float *)r;		// trackball method

-(void)getRotationMatrix:(float *)mat;
-(void)setProjection:(int)p;
-(void)projection;
-(void)projection2; // only used when saving a stereographic projection

-(void)setShading:(int)p;

-(void)setZoom:(float)z;

-(void)savePicture:(NSString *)filename;
-(void)getPixels:(char*)baseaddr width:(long)w height:(long)h rowbyte:(long)rb;
-(void)getDepthMap:(float*)baseaddr width:(long)w height:(long)h;

-(void)pickVertex:(NSPoint)mp;
-(void)dragVertex:(NSPoint)mp;
-(void)setSelected:(int)index;
-(int)selected;

-(void)init_vertsnormals:(float*)v;
-(void)recomputeNormals:(float*)v;

void invMat(float *a,float *b);

-(void)spin:(char*)path nframes:(int)nframes;
-(void)morphToMeshAtPath:(char*)mshpath destination:(char*)imgpath nframes:(int)nframes;

-(IBAction)emStandard:(id)sender;
-(IBAction)emProjection:(id)sender;		//
-(IBAction)emColourmap:(id)sender;		//
-(IBAction)emRotate:(id)sender;
-(IBAction)emSaveDepthMap:(id)sender;	//
-(IBAction)emSaveImage:(id)sender;		//
-(IBAction)emShading:(id)sender;		//
-(IBAction)emZoom:(id)sender;
@end
