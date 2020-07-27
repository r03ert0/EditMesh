#import "MyMeshView.h"

@implementation MyMeshView
#pragma mark -
void invMat(float *a,float *b);

void inverse4x4(float b[16],float a[16])
{
	float d=	-a[0]*a[5]*a[10]*a[15]+a[0]*a[5]*a[11]*a[14]
				+a[0]*a[9]*a[6]*a[15]-a[0]*a[9]*a[7]*a[14]
				-a[0]*a[13]*a[6]*a[11]+a[0]*a[13]*a[7]*a[10]
				+a[4]*a[1]*a[10]*a[15]-a[4]*a[1]*a[11]*a[14]
				-a[4]*a[9]*a[2]*a[15]+a[4]*a[9]*a[3]*a[14]
				+a[4]*a[13]*a[2]*a[11]-a[4]*a[13]*a[3]*a[10]
				-a[8]*a[1]*a[6]*a[15]+a[8]*a[1]*a[7]*a[14]
				+a[8]*a[5]*a[2]*a[15]-a[8]*a[5]*a[3]*a[14]
				-a[8]*a[13]*a[2]*a[7]+a[8]*a[13]*a[3]*a[6]
				+a[12]*a[1]*a[6]*a[11]-a[12]*a[1]*a[7]*a[10]
				-a[12]*a[5]*a[2]*a[11]+a[12]*a[5]*a[3]*a[10]
				+a[12]*a[9]*a[2]*a[7]-a[12]*a[9]*a[3]*a[6];
	
	b[0]=-(a[5]*a[10]*a[15]-a[5]*a[11]*a[14]-a[9]*a[6]*a[15]+a[9]*a[7]*a[14]+a[13]*a[6]*a[11]-a[13]*a[7]*a[10])/d;
	b[1]=(a[1]*a[10]*a[15]-a[1]*a[11]*a[14]-a[9]*a[2]*a[15]+a[9]*a[3]*a[14]+a[13]*a[2]*a[11]-a[13]*a[3]*a[10])/d;
	b[2]=-(a[1]*a[6]*a[15]-a[1]*a[7]*a[14]-a[5]*a[2]*a[15]+a[5]*a[3]*a[14]+a[13]*a[2]*a[7]-a[13]*a[3]*a[6])/d;
	b[3]=(a[1]*a[6]*a[11]-a[1]*a[7]*a[10]-a[5]*a[2]*a[11]+a[5]*a[3]*a[10]+a[9]*a[2]*a[7]-a[9]*a[3]*a[6])/d;
	
	b[4]=(a[4]*a[10]*a[15]-a[4]*a[11]*a[14]-a[8]*a[6]*a[15]+a[8]*a[7]*a[14]+a[12]*a[6]*a[11]-a[12]*a[7]*a[10])/d;
	b[5]=(-a[0]*a[10]*a[15]+a[0]*a[11]*a[14]+a[8]*a[2]*a[15]-a[8]*a[3]*a[14]-a[12]*a[2]*a[11]+a[12]*a[3]*a[10])/d;
	b[6]=-(-a[0]*a[6]*a[15]+a[0]*a[7]*a[14]+a[4]*a[2]*a[15]-a[4]*a[3]*a[14]-a[12]*a[2]*a[7]+a[12]*a[3]*a[6])/d;
	b[7]=(-a[0]*a[6]*a[11]+a[0]*a[7]*a[10]+a[4]*a[2]*a[11]-a[4]*a[3]*a[10]-a[8]*a[2]*a[7]+a[8]*a[3]*a[6])/d;
	
	b[8]=-(a[4]*a[9]*a[15]-a[4]*a[11]*a[13]-a[8]*a[5]*a[15]+a[8]*a[7]*a[13]+a[12]*a[5]*a[11]-a[12]*a[7]*a[9])/d;
	b[9]=-(-a[0]*a[9]*a[15]+a[0]*a[11]*a[13]+a[8]*a[1]*a[15]-a[8]*a[3]*a[13]-a[12]*a[1]*a[11]+a[12]*a[3]*a[9])/d;
	b[10]=(-a[0]*a[5]*a[15]+a[0]*a[7]*a[13]+a[4]*a[1]*a[15]-a[4]*a[3]*a[13]-a[12]*a[1]*a[7]+a[12]*a[3]*a[5])/d;
	b[11]=-(-a[0]*a[5]*a[11]+a[0]*a[7]*a[9]+a[4]*a[1]*a[11]-a[4]*a[3]*a[9]-a[8]*a[1]*a[7]+a[8]*a[3]*a[5])/d;

	b[12]=(a[4]*a[9]*a[14]-a[4]*a[10]*a[13]-a[8]*a[5]*a[14]+a[8]*a[6]*a[13]+a[12]*a[5]*a[10]-a[12]*a[6]*a[9])/d;
	b[13]=(-a[0]*a[9]*a[14]+a[0]*a[10]*a[13]+a[8]*a[1]*a[14]-a[8]*a[2]*a[13]-a[12]*a[1]*a[10]+a[12]*a[2]*a[9])/d;
	b[14]=-(-a[0]*a[5]*a[14]+a[0]*a[6]*a[13]+a[4]*a[1]*a[14]-a[4]*a[2]*a[13]-a[12]*a[1]*a[6]+a[12]*a[2]*a[5])/d;
	b[15]=(-a[0]*a[5]*a[10]+a[0]*a[6]*a[9]+a[4]*a[1]*a[10]-a[4]*a[2]*a[9]-a[8]*a[1]*a[6]+a[8]*a[2]*a[5])/d;
}
void v_m(float *r,float *v,float *m)
{
	// v=1x3
	// m=4x4
	// r=1x3
	r[0]=v[0]*m[0*4+0]+v[1]*m[1*4+0]+v[2]*m[2*4+0] + m[3*4+0];
	r[1]=v[0]*m[0*4+1]+v[1]*m[1*4+1]+v[2]*m[2*4+1] + m[3*4+1];
	r[2]=v[0]*m[0*4+2]+v[1]*m[1*4+2]+v[2]*m[2*4+2] + m[3*4+2];
}
// Override NSView's initWithFrame: to specify our pixel format:
- (id) initWithFrame: (NSRect) frame
{
	// 1. Initialize pixel format
    GLuint attribs[] = 
    {
/*
        NSOpenGLPFAOpenGLProfile,
        NSOpenGLProfileVersion3_2Core,
        NSOpenGLPFAColorSize, 24,
        NSOpenGLPFAAlphaSize, 8,
        NSOpenGLPFADoubleBuffer,
        NSOpenGLPFAAccelerated,
*/
            NSOpenGLPFANoRecovery,
            NSOpenGLPFAColorSize, 24,
            NSOpenGLPFAAlphaSize, 8,
            NSOpenGLPFADoubleBuffer,
            NSOpenGLPFAAccelerated,
            NSOpenGLPFADepthSize, 24,
            NSOpenGLPFAStencilSize, 8,
            NSOpenGLPFAAccumSize, 0,
            0
    };

    NSOpenGLPixelFormat* fmt = [[NSOpenGLPixelFormat alloc] initWithAttributes: (NSOpenGLPixelFormatAttribute*) attribs];
    self = [super initWithFrame:frame pixelFormat: [fmt autorelease]];
    if (!fmt)	NSLog(@"No OpenGL pixel format");

    [[self openGLContext] makeCurrentContext];

    
    // 2. Init GL
    //NOTEXTURE glEnable(GL_TEXTURE_2D);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_SMOOTH);
	//NOCULLFACE glEnable(GL_CULL_FACE);
    //NOBACK glCullFace(GL_BACK);
	
    
    // 3. Init default mesh
    tris  = (int *) calloc( 6, sizeof(int) );
    verts = (GLfloat *) calloc( 4*3, sizeof(GLfloat) );
	vertsproj = (GLfloat *) calloc( 4*3, sizeof(GLfloat) );
    vertscolour = (GLfloat *) calloc( 4*3, sizeof(GLfloat) );
    vertsparam = (GLfloat *) calloc( 4*2, sizeof(GLfloat) );
    vertsnormals = (GLfloat *) calloc( 4*3, sizeof(GLfloat) );
    //NOTEXTURE texture = (GLubyte *) calloc( 64*64*3, sizeof(GLubyte) );

    tris[0] = 0; tris[1] = 2; tris[2] = 1;
    tris[3] = 2; tris[4] = 3; tris[5] = 1;
    ntris=2;

    verts[0] = 0.0f;  verts[1] = -2.0f; verts[2] = -2.0f;
    verts[3] = 1.0f;  verts[4] =  0.0f; verts[5] = 0.0f;
    verts[6] = -1.0f; verts[7] =  0.0f; verts[8] = 0.1f;
    verts[9] = 0.0f;  verts[10] = 2.0f; verts[11]= 0.2f;
	nverts=4;

    vertscolour[0] = 1;  vertscolour[1] = 0; vertscolour[2] = 0;
    vertscolour[3] = 1;  vertscolour[4] = 0; vertscolour[5] = 0;
    vertscolour[6] = 1;  vertscolour[7] = 0; vertscolour[8] = 1;
    vertscolour[9] = 1;  vertscolour[10] =0; vertscolour[11]= 0;
    
    vertsparam[0] = 0.0f;  	vertsparam[1] = 0.0f;
    vertsparam[2] = 1.0f;	vertsparam[3] = 0.0f;
    vertsparam[4] = 0.0f; 	vertsparam[5] = 1.0f;
    vertsparam[6] = 1.0f;	vertsparam[7] = 1.0f;

	//---------- Cel-shading texture ----------
	float shaderData[32][3];
	int	i;
	for(i=0;i<32;i++)
	{
		if(i>16)	{shaderData[i][0]=shaderData[i][1]=shaderData[i][2]=1.00;}
		else if(i>8){shaderData[i][0]=shaderData[i][1]=shaderData[i][2]=0.75;}
		else		{shaderData[i][0]=shaderData[i][1]=shaderData[i][2]=0.50;}
		//	shaderData[i][0]=shaderData[i][1]=shaderData[i][2]=i/32.0;
	}
	glGenTextures( 1, &shaderTexture);
	glBindTexture( GL_TEXTURE_1D, shaderTexture);
	glTexParameteri( GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	glTexParameteri( GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
	glTexImage1D( GL_TEXTURE_1D, 0, GL_RGB, 32, 0, GL_RGB , GL_FLOAT, shaderData );
	//-----------------------------------------

    //NOTEXTURE
	/*
	int i, size=64;
	for(i=0;i<size*size*3;i++)
    {	if(i<size*size*3/2)	texture[i]=255*(i%17);
        else			texture[i]=0;
    }
    glGenTextures(1,&textureId);
    glBindTexture(GL_TEXTURE_2D, textureId);
    glTexImage2D(GL_TEXTURE_2D, 0,3, size, size, 0, GL_RGB, GL_UNSIGNED_BYTE,texture);
    free(texture);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	[self setTextureActive:YES];
	*/

    // 4. Initialize the trackball.
    m_trackball = [[Trackball alloc] init];
    m_rotation[0] = m_tbRot[0] = 0.0;
    m_rotation[1] = m_tbRot[1] = 1.0;
    m_rotation[2] = m_tbRot[2] = 0.0;
    m_rotation[3] = m_tbRot[3] = 0.0;
    
	// initialize cursor
	selected=5371;//TEST-1;
	cp[0]=(float3D){1,1,-1};
	cp[1]=(float3D){-1,1,-1};
	cp[2]=(float3D){-1,-1,-1};
	cp[3]=(float3D){1,-1,-1};
	cp[4]=(float3D){1,1,1};
	cp[5]=(float3D){-1,1,1};
	cp[6]=(float3D){-1,-1,1};
	cp[7]=(float3D){1,-1,1};
	ct[0]=(int3D){0,1,5};
	ct[1]=(int3D){0,5,4};
	ct[2]=(int3D){1,2,6};
	ct[3]=(int3D){1,6,5};
	ct[4]=(int3D){2,3,7};
	ct[5]=(int3D){2,7,6};
	ct[6]=(int3D){3,0,4};
	ct[7]=(int3D){3,4,7};
	ct[8]=(int3D){0,2,1};
	ct[9]=(int3D){0,3,2};
	ct[10]=(int3D){4,5,6};
	ct[11]=(int3D){4,6,7};
	cc[0]=(float3D){0,1,0};
	cc[1]=(float3D){0,1,0};
	cc[2]=(float3D){1,0,1};
	cc[3]=(float3D){1,0,1};
	cc[4]=(float3D){1,1,0};
	cc[5]=(float3D){1,1,0};
	cc[6]=(float3D){1,0,0};
	cc[7]=(float3D){1,0,0};
	cc[8]=(float3D){0,1,1};
	cc[9]=(float3D){0,1,1};
	cc[10]=(float3D){0,0,1};
	cc[11]=(float3D){0,0,1};

    // 5. Initialize zoom
	zoom=2;
	
	// 6. Initialize projection to orthographic
	proj=0;
	
	// 7. init colourmap to Jet
	cmap_index=GREY;
	
	// 8. init shading to depth
	shading=0;
	
	vertsnormals_computed=0;
    
    return self;
}
-(BOOL)acceptsFirstResponder
{
	return YES;
}
#pragma mark -
-(float3D)mm :(float3D)V :(float*)M
{
	float3D	r;
	r.x=M[0]*V.x+M[4]*V.y+M[8]*V.z;    // Rotate Around The X Axis
	r.y=M[1]*V.x+M[5]*V.y+M[9]*V.z;    // Rotate Around The Y Axis
	r.z=M[2]*V.x+M[6]*V.y+M[10]*V.z;   // Rotate Around The Z Axis

	return r;
}

- (void) drawRect: (NSRect) rect
{
    float	aspectRatio;
	float3D	tmp;
	int		i;
    
    [self update];

    // init projection
        glViewport(0, 0, (GLsizei) rect.size.width, (GLsizei) rect.size.height);
        glClearColor(1,1,1, 1);
        //glClear(GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT+GL_STENCIL_BUFFER_BIT);
        glClear(GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        aspectRatio = (float)rect.size.width/(float)rect.size.height;
		glOrtho(-aspectRatio*zoom, aspectRatio*zoom, -1.0*zoom, 1.0*zoom, -3000.0, 3000.0);

    // prepare drawing
        glMatrixMode (GL_MODELVIEW);
        glLoadIdentity();
		gluLookAt (0,0,-10, 0,0,0, 0,1,0); // eye,center,updir
        glRotatef(m_tbRot[0],m_tbRot[1], m_tbRot[2], m_tbRot[3]);
        glRotatef(m_rotation[0],m_rotation[1],m_rotation[2],m_rotation[3]);

		glGetFloatv(GL_MODELVIEW_MATRIX, m);
		m[3]=m[7]=m[11]=m[12]=m[13]=m[14]=0;
		inverse4x4(invm,m);

    // draw
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glEnable( GL_POLYGON_SMOOTH );
	
	if(shading==0)	// shading=0 -> Depth shading
	{
		if(proj==0)
		{
			glEnableClientState( GL_VERTEX_ARRAY);
			glVertexPointer( 3, GL_FLOAT, 0, verts );
			if(vertscolour)
			{	glEnableClientState( GL_COLOR_ARRAY );
				glColorPointer( 3, GL_FLOAT, 0, vertscolour );
			}
			if(0)//hasTexture && vertsparam)
			{	glEnableClientState( GL_TEXTURE_COORD_ARRAY );
				glTexCoordPointer( 2, GL_FLOAT, 0, vertsparam );
			}
			glDrawElements( GL_TRIANGLES, ntris*3, GL_UNSIGNED_INT, tris );
			glDisableClientState( GL_VERTEX_ARRAY);
		}
		else
		if(proj>=1&&proj<=3)
		{
			[self projection];
			glBegin(GL_TRIANGLES);
			int i,mx=3;
			float3D	a,b,c;
			
			for(i=0;i<ntris;i++)
			{
				a=*(float3D*)&vertsproj[3*tris[3*i+0]];
				b=*(float3D*)&vertsproj[3*tris[3*i+1]];
				c=*(float3D*)&vertsproj[3*tris[3*i+2]];
				if(norm3D(sub3D(a,b))<mx && norm3D(sub3D(b,c))<mx && norm3D(sub3D(c,a))<mx)
				{
					glColor3fv(&vertscolour[3*tris[3*i+0]]);
					glVertex3fv((float*)&a);
					glColor3fv(&vertscolour[3*tris[3*i+1]]);
					glVertex3fv((float*)&b);
					glColor3fv(&vertscolour[3*tris[3*i+2]]);
					glVertex3fv((float*)&c);
				}
			}
			glEnd();
		}
	}
	else	// Cel or Toon shading
	{
		int		i,j,mx=3;
		float3D	nn,a,b,c;
		float	d,rmat[9];
		
		rmat[0]=m[0]; rmat[1]=m[1]; rmat[2]=m[2];
		rmat[3]=m[4]; rmat[4]=m[5]; rmat[5]=m[6];
		rmat[6]=m[8]; rmat[7]=m[9]; rmat[8]=m[10];
		nn=matfloat3D(rmat,(float3D){0,0,-1});

		[self projection];
		vertsnormals_computed=0;
		[self init_vertsnormals:vertsproj];
		
		if(shading==2) // toon shading
		{
			glEnable( GL_CULL_FACE );
			glPolygonMode( GL_BACK, GL_FILL );
			glCullFace( GL_FRONT );
		}
	
		glEnable( GL_TEXTURE_1D );
		glBindTexture( GL_TEXTURE_1D, shaderTexture);

		glBegin(GL_TRIANGLES);
		for(i=0;i<ntris;i++)
		{
			a=*(float3D*)&vertsproj[3*tris[3*i+0]];
			b=*(float3D*)&vertsproj[3*tris[3*i+1]];
			c=*(float3D*)&vertsproj[3*tris[3*i+2]];
			if(proj>=1&&proj<=3 && (norm3D(sub3D(a,b))>mx || norm3D(sub3D(b,c))>mx || norm3D(sub3D(c,a))>mx))
				continue;
			
			j=3*tris[3*i+0];
			d=dot3D(*(float3D*)&vertsnormals[j],nn);
			if(d<0) d=0;
			glTexCoord1f(d);
			glColor3fv(&vertscolour[j]);
			glVertex3fv((float*)&a);
			
			j=3*tris[3*i+1];
			d=dot3D(*(float3D*)&vertsnormals[j],nn);
			if(d<0) d=0;
			glTexCoord1f(d);
			glColor3fv(&vertscolour[j]);
			glVertex3fv((float*)&b);
			
			j=3*tris[3*i+2];
			d=dot3D(*(float3D*)&vertsnormals[j],nn);
			if(d<0) d=0;
			glTexCoord1f(d);
			glColor3fv(&vertscolour[j]);
			glVertex3fv((float*)&c);
		}
		glEnd();
		glDisable( GL_TEXTURE_1D );

		if(shading==2) // toon shading
		{
			glPolygonMode(GL_FRONT, GL_LINE);
			glLineWidth(3.0);
			glCullFace(GL_BACK);
			glDepthFunc(GL_LESS);
			glColor3f(0,0,0);
			glBegin(GL_TRIANGLES);
			for(i=0;i<ntris;i++)
			{
				a=*(float3D*)&vertsproj[3*tris[3*i+0]];
				b=*(float3D*)&vertsproj[3*tris[3*i+1]];
				c=*(float3D*)&vertsproj[3*tris[3*i+2]];
				if(proj>=1&&proj<=3 && (norm3D(sub3D(a,b))>mx || norm3D(sub3D(b,c))>mx || norm3D(sub3D(c,a))>mx))
					continue;
				
				glVertex3fv((float*)&a);
				glVertex3fv((float*)&b);
				glVertex3fv((float*)&c);
			}
			glEnd();
			glDisable( GL_CULL_FACE );
		}

	}

	// Draw cursor
	glBegin(GL_TRIANGLES);
	if(0)//selected>=0)
	for(i=0;i<12;i++)
	{
		glColor3fv((float*)&cc[i]);
		tmp=add3D(sca3D(cp[ct[i].a],0.02*zoom),((float3D*)verts)[selected]); glVertex3fv((float*)&tmp);
		tmp=add3D(sca3D(cp[ct[i].b],0.02*zoom),((float3D*)verts)[selected]); glVertex3fv((float*)&tmp);
		tmp=add3D(sca3D(cp[ct[i].c],0.02*zoom),((float3D*)verts)[selected]); glVertex3fv((float*)&tmp);
	}
	glEnd();

    [[self openGLContext] flushBuffer];
}
#pragma mark -
-(void)setSettings:(NSMutableDictionary*)newSettings
{
	settings=newSettings;
}
#pragma mark -
-(float*)verts
{
	return verts;
}
-(float*)vertsproj
{
	return vertsproj;
}-(int)nverts
{
	return nverts;
}
- (void) setVertices: (float *) vert number:(int)n
{
	verts=(GLfloat*)vert;
	nverts=n;
	
	free(vertscolour);
	vertscolour=nil;

	free(vertsproj);
	vertsproj = (GLfloat *) calloc( n*3, sizeof(GLfloat) );
	
	free(vertsnormals);
	vertsnormals=(GLfloat *)calloc(n*3,sizeof(GLfloat));
	vertsnormals_computed=0;
}
-(int*)tris
{
	return tris;
}
-(int)ntris
{
	return ntris;
}
- (void) setTriangles: (int *) trian number:(int)n
{
    tris=(int*)trian;
    ntris=n;
}
-(void)init_vertsnormals:(float*)V
{
	if(vertsnormals_computed==1)
		return;
	[self recomputeNormals:V];
}
-(void)recomputeNormals:(float*)V
{
	int	i,j;
	unsigned char *tmp;
	float3D	v,*vv;
	tmp=(unsigned char*)calloc(nverts,sizeof(char));
	
	for(i=0;i<3*nverts;i++)
		vertsnormals[i]=0;
	
	for(i=0;i<ntris;i++)
	{
		v=triPlane(*(float3D*)&V[3*tris[3*i+0]],*(float3D*)&V[3*tris[3*i+1]],*(float3D*)&V[3*tris[3*i+2]]);
		for(j=0;j<3;j++)
		{
			vertsnormals[3*tris[3*i+j]+0]+=v.x;
			vertsnormals[3*tris[3*i+j]+1]+=v.y;
			vertsnormals[3*tris[3*i+j]+2]+=v.z;
			tmp[tris[3*i+j]]++;
			tmp[tris[3*i+j]]++;
			tmp[tris[3*i+j]]++;
		}
	}
	for(i=0;i<nverts;i++)
	{
		vv=(float3D*)&vertsnormals[3*i];
		*vv=sca3D(*vv,1/(float)tmp[i]);
		*vv=sca3D(*vv,1/norm3D(*vv));
	}
	free(tmp);
	
	vertsnormals_computed=1;
}

-(void)smooth:(float*)txr
{
	int		i,j;
	float	*tmp=(float*)calloc(nverts,sizeof(float));
	char	*itmp=(char*)calloc(nverts,sizeof(char));
	
	for(j=0;j<20;j++)
	{
		for(i=0;i<ntris;i++)
		{
			tmp[tris[3*i+0]]+=txr[3*tris[3*i+1]]+txr[3*tris[3*i+2]];
			tmp[tris[3*i+1]]+=txr[3*tris[3*i+2]]+txr[3*tris[3*i+0]];
			tmp[tris[3*i+2]]+=txr[3*tris[3*i+0]]+txr[3*tris[3*i+1]];
			itmp[tris[3*i+0]]+=2;
			itmp[tris[3*i+1]]+=2;
			itmp[tris[3*i+2]]+=2;
		}
		for(i=0;i<nverts;i++)
		{
			txr[3*i+0]=tmp[i]/(float)itmp[i];
			txr[3*i+1]=tmp[i]/(float)itmp[i];
			txr[3*i+2]=tmp[i]/(float)itmp[i];
			tmp[i]=0;
			itmp[i]=0;
		}
	}
	free(tmp);
	free(itmp);
}
-(void)multiplySulcalDepth
{
	int			i;
    float		n,max,*d;
	float3D		p,ce={0,0,0},ide,siz;
	
	// compute sulcal depth
	for(i=0;i<nverts;i++)
	{
		p=*(float3D*)&verts[3*i];
		ce=(float3D){ce.x+p.x,ce.y+p.y,ce.z+p.z};
		
		if(i==0) ide=siz=p;
		
		if(ide.x<p.x) ide.x=p.x;
		if(ide.y<p.y) ide.y=p.y;
		if(ide.z<p.z) ide.z=p.z;
		
		if(siz.x>p.x) siz.x=p.x;
		if(siz.y>p.y) siz.y=p.y;
		if(siz.z>p.z) siz.z=p.z;
	}
	ce=(float3D){ce.x/(float)nverts,ce.y/(float)nverts,ce.z/(float)nverts};

	max=0;
	d=(float*)calloc(nverts,sizeof(float));
    for(i=0;i<nverts;i++)
    {
		p=*(float3D*)&verts[3*i];
        n=	pow(2*(p.x-ce.x)/(ide.x-siz.x),2) +
			pow(2*(p.y-ce.y)/(ide.y-siz.y),2) +
			pow(2*(p.z-ce.z)/(ide.z-siz.z),2);

        d[i] = sqrt(n);
        if(d[i]>max)	max=d[i];
    }
    max*=1.05;	// pure white is not nice...
    for(i=0;i<nverts;i++)
	{
        /*vertscolour[3*i+0]=vertscolour[3*i+0]*(1+d[i]/max)*0.5;
        vertscolour[3*i+1]=vertscolour[3*i+1]*(1+d[i]/max)*0.5;
        vertscolour[3*i+2]=vertscolour[3*i+2]*(1+d[i]/max)*0.5;*/
        vertscolour[3*i+0]=vertscolour[3*i+0]*d[i]/max;
        vertscolour[3*i+1]=vertscolour[3*i+1]*d[i]/max;
        vertscolour[3*i+2]=vertscolour[3*i+2]*d[i]/max;
	}
	free(d);
}
-(void)setColourmap:(int)i
{
	switch(i)
	{
		case 0: cmap_index=JET;		break;
		case 1: cmap_index=HOT;		break;
		case 2: cmap_index=GREY;	break;
		case 3: cmap_index=RED;		break;
		case 4: cmap_index=GREEN;	break;
		case 5: cmap_index=BLUE;	break;
	}
	
	free(vertscolour);
	vertscolour=nil;
	
	if(vertsdata)
		[self setMinMaxData];
	else
		[self setVerticesColour:nil];
}
- (void)setMinMaxData
{
    int				i;
	float			min,max,val;
	unsigned char	c[3],*str;
	
	if(vertsdata==nil)
		return;
	
	if(vertscolour==nil)
	{
		vertscolour=(GLfloat*)calloc(nverts*3,sizeof(GLfloat));
		for(i=0;i<nverts*3;i++)
			vertscolour[i]=1;
	}

	str=(unsigned char*)[[settings objectForKey:@"cmapminmax"] UTF8String];
	printf("%s\n",str);
	sscanf((char*)str," %f , %f ",&min,&max); printf("min:%f, max:%f\n",min,max);
	
	if(1)
	for(i=0;i<nverts;i++)
	{
		val=(vertsdata[i]-min)/(max-min);
		if(val<0) val=0;
		if(val>1) val=1;
		colourmap(val,c,cmap_index);

		vertscolour[3*i+0]=c[0]/255.0;
		vertscolour[3*i+1]=c[1]/255.0;
		vertscolour[3*i+2]=c[2]/255.0;
	}
	
	//[self multiplySulcalDepth];
}
- (void) setVerticesData: (float *) vertdata
{
    int				i,negpos_cm=0;
	float			x,max,min;
	unsigned char	c[3];
	
	if(vertscolour==nil)
	{
		vertscolour=(GLfloat*)calloc(nverts*3,sizeof(GLfloat));
		for(i=0;i<nverts*3;i++)
			vertscolour[i]=1;
	}

	// search max and min values
	for(i=0;i<nverts;i++)
	{
		if(i==0) max=min=vertdata[0];
		
		if(negpos_cm)
		{
			if(fabs(vertdata[i])<min) min=fabs(vertdata[i]);
			if(fabs(vertdata[i])>max) max=fabs(vertdata[i]);
		}
		else
		{
			if(vertdata[i]<min) min=vertdata[i];
			if(vertdata[i]>max) max=vertdata[i];
		}
	}
	[settings setObject:[NSString stringWithFormat:@"%f, %f",min,max] forKey:@"cmapminmax"];
	if(1)
	for(i=0;i<nverts;i++)
	{
		if(negpos_cm)
		{
			c[0]=c[1]=c[2]=0;
			if(vertdata[i]<0)
				c[1]=255*(fabs(vertdata[i])-min)/(max-min);
			else
				c[0]=255*(fabs(vertdata[i])-min)/(max-min);
		}
		else
		{
			x=(vertdata[i]-min)/(max-min);
			colourmap(x,c,cmap_index);
		}

		vertscolour[3*i+0]=c[0]/255.0;
		vertscolour[3*i+1]=c[1]/255.0;
		vertscolour[3*i+2]=c[2]/255.0;
	}
	if(vertsdata!=vertdata)
		free(vertsdata);
	vertsdata=vertdata;
	
	[self multiplySulcalDepth];
}
- (float*)verticesData
{
	return vertsdata;
}
- (void) setVerticesColour: (float *) vertcolour
{
    int		i;
	//float	max,min;
	//unsigned char c[3];
	
	if(vertscolour==nil)
	{
		vertscolour=(GLfloat*)calloc(nverts*3,sizeof(GLfloat));
		for(i=0;i<nverts*3;i++)
			vertscolour[i]=1;
	}

    /*
	// search max and min values
	max=min=0;
	if(vertcolour)
		for(i=0;i<nverts;i++)
		{
			if(i==0) max=min=vertcolour[0];
			if(vertcolour[3*i]<min) min=vertcolour[3*i];
			if(vertcolour[3*i]>max) max=vertcolour[3*i];
		}
	
	for(i=0;i<nverts;i++)
	{
		if(max>min)
			colourmap((vertcolour[3*i]-min)/(max-min),c,cmap_index);
		else
			colourmap(1.0,c,cmap_index);
		vertscolour[3*i+0]=c[0]/255.0;
		vertscolour[3*i+1]=c[1]/255.0;
		vertscolour[3*i+2]=c[2]/255.0;
	}
     */

    for(i=0;i<nverts;i++)
	{
		vertscolour[3*i+0]=vertcolour[3*i+0]/255.0;
		vertscolour[3*i+1]=vertcolour[3*i+1]/255.0;
		vertscolour[3*i+2]=vertcolour[3*i+2]/255.0;
	}

	//[self multiplySulcalDepth];
}
-(float*)verticesColour
{
	return vertscolour;
}
- (void) setParam: (float *)param
{
    vertsparam=(GLfloat*)param;
}
- (void) setTexture: (char *)image size:(NSSize)s
{
    texture=(GLubyte*)image;
    textureSize=s;

    glGenTextures(1,&textureId);
    glBindTexture(GL_TEXTURE_2D, textureId);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0,3, textureSize.width,textureSize.height, 0, GL_RGB, GL_UNSIGNED_BYTE,texture);
}
- (void) setTextureActive:(BOOL)isActive
{
    hasTexture=isActive;
    if(isActive)	glEnable(GL_TEXTURE_2D);
    else			glDisable(GL_TEXTURE_2D);
    [self setNeedsDisplay:YES];
}
#pragma mark -
- (void)setStandardRotation:(int)view
{
    m_rotation[0] = m_tbRot[0] = 0.0;
    m_rotation[1] = m_tbRot[1] = 0.0;
    m_rotation[2] = m_tbRot[2] = 1.0;
    m_rotation[3] = m_tbRot[3] = 0.0;
    
    switch(view)
    {
        case 1:m_rotation[0]=270;	m_rotation[1]=1;m_rotation[2]=0; break; //sup
        case 4:m_rotation[0]= 90;	break; //frn
        case 5:m_rotation[0]=  0;	break; //tmp
        case 6:m_rotation[0]=270;	break; //occ
        case 7:m_rotation[0]=180;	break; //med
        case 9:m_rotation[0]= 90;	m_rotation[1]=1;m_rotation[2]=0; break; //cau
    }
    [self setNeedsDisplay:YES];
}
-(void)addRotation:(float)value toAxis:(int)axis
{
	float	tmp[4]={0,0,0,0};
	tmp[axis]=1;
	tmp[0]=value;
	[m_trackball add:tmp toRotation:m_rotation];
	[self setNeedsDisplay:YES];
}
- (void)rotateBy:(float *)r
{
    m_tbRot[0] = r[0];
    m_tbRot[1] = r[1];
    m_tbRot[2] = r[2];
    m_tbRot[3] = r[3];
}

- (void)mouseDown:(NSEvent *)theEvent
{
	if([theEvent modifierFlags]&NSAlternateKeyMask)
		[self pickVertex:[self convertPoint:[theEvent locationInWindow] fromView:nil]];
	else
		[m_trackball  start:[theEvent locationInWindow] sender:self];
}

- (void)mouseUp:(NSEvent *)theEvent
{
    // Accumulate the trackball rotation
    // into the current rotation.
    [m_trackball add:m_tbRot toRotation:m_rotation];

    m_tbRot[0]=0;
    m_tbRot[1]=1;
    m_tbRot[2]=0;
    m_tbRot[3]=0;
}

- (void)mouseDragged:(NSEvent *)theEvent
{
    [self lockFocus];
	if([theEvent modifierFlags]&NSAlternateKeyMask)
		[self dragVertex:[self convertPoint:[theEvent locationInWindow] fromView:nil]];
	else
		[m_trackball rollTo:[theEvent locationInWindow] sender:self];
    [self unlockFocus];
    [self setNeedsDisplay:YES];
}
-(void)scrollWheel:(NSEvent *)theEvent
{
	float	resolution=0.1;
	zoom = pow(2,log(zoom)/log(2)-[theEvent deltaY]*resolution);
	[self setNeedsDisplay:YES];
}
-(void)keyDown:(NSEvent*)e
{
	printf("key code:%i\n",[e keyCode]);
	if([e keyCode]==76||[e keyCode]==36)	// delete
	{
	}
	else
	if([e keyCode]==76||[e keyCode]==36)	// return and enter
	{
	}
	else
	if([e keyCode]==126)	// up arrow
	{
	}
	else
	if([e keyCode]==125)	// down arrow
	{
	}
	else
	if([e keyCode]==48)	// tab
	{
	}
	else
		[super keyDown:e];
}
#pragma mark -
-(void)setShading:(int)s
{
	shading=s;
}
#pragma mark -
- (void)getRotationMatrix:(float *)mat
{
	glGetFloatv(GL_MODELVIEW_MATRIX,mat);
}
- (void) setProjection:(int)p
{
	proj=p;
	[self setNeedsDisplay:YES];
}
void invMat(float *a,float *b)
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

void dynstereographic(float *mat, float *src, float *dst);
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
	{	sp[0]=-h; sp[1]=v; sp[2]=0;}
    else
	{   sp[0]=h; sp[1]=v; sp[2]=0;}
    
    dst[0]=mat[0]*sp[0]+mat[3]*sp[1]+mat[6]*sp[2];
	dst[1]=mat[1]*sp[0]+mat[4]*sp[1]+mat[7]*sp[2];
	dst[2]=mat[2]*sp[0]+mat[5]*sp[1]+mat[8]*sp[2];
}
void dynstereographic2(float *mat, float *src, float *dst);
void dynstereographic2(float *mat, float *src, float *dst)
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
	{	sp[0]=-h; sp[1]=v; sp[2]=0;}
    else
	{   sp[0]=h; sp[1]=v; sp[2]=0;}
    
    dst[0]=sp[0];
	dst[1]=sp[1];
	dst[2]=sp[2];
}
void dynsinusoidal(float *mat, float *src, float *dst);
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
void dynmercator(float *mat, float *src, float *dst);
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

-(void) projection
{
	float   mat[16],rmat[9],imat[9];
	float   vec[9],ce[3]={0,0,0},tmp[3];
	int		i;
	
	[self getRotationMatrix:mat];
	rmat[0]=mat[0]; rmat[1]=mat[1]; rmat[2]=mat[2];
	rmat[3]=mat[4]; rmat[4]=mat[5]; rmat[5]=mat[6];
	rmat[6]=mat[8]; rmat[7]=mat[9]; rmat[8]=mat[10];
	invMat(rmat,imat);
	for(i=0;i<9;i++) vec[i]=imat[i];
	
	for(i=0;i<nverts;i++)
	{
		ce[0]+=verts[3*i+0];
		ce[1]+=verts[3*i+1];
		ce[2]+=verts[3*i+2];
	}
	ce[0] /=(float)nverts;
	ce[1] /=(float)nverts;
	ce[2] /=(float)nverts;
	
	switch(proj)
	{
		case 0:
			for(i=0;i<3*nverts;i++) vertsproj[i]=verts[i];
			break;
		case 1: // stereographic
			for(i=0;i<nverts;i++)
			{
				tmp[0]=verts[3*i+0]-ce[0];
				tmp[1]=verts[3*i+1]-ce[1];
				tmp[2]=verts[3*i+2]-ce[2];
				dynstereographic(vec,tmp,&vertsproj[3*i]);
			}
			break;
		case 2: // sinusoidal
			for(i=0;i<nverts;i++) dynsinusoidal(vec,&verts[3*i],&vertsproj[3*i]);
			break;
		case 3: // mercator
			for(i=0;i<nverts;i++) dynmercator(vec,&verts[3*i],&vertsproj[3*i]);
			break;
	}
}
-(void) projection2
{
	float   mat[16],rmat[9],imat[9];
	float   vec[9],ce[3]={0,0,0},tmp[3];
	int		i;
	
	[self getRotationMatrix:mat];
	rmat[0]=mat[0]; rmat[1]=mat[1]; rmat[2]=mat[2];
	rmat[3]=mat[4]; rmat[4]=mat[5]; rmat[5]=mat[6];
	rmat[6]=mat[8]; rmat[7]=mat[9]; rmat[8]=mat[10];
	invMat(rmat,imat);
	for(i=0;i<9;i++) vec[i]=imat[i];
	
	for(i=0;i<nverts;i++)
	{
		ce[0]+=verts[3*i+0];
		ce[1]+=verts[3*i+1];
		ce[2]+=verts[3*i+2];
	}
	ce[0] /=(float)nverts;
	ce[1] /=(float)nverts;
	ce[2] /=(float)nverts;
	
    for(i=0;i<nverts;i++)
    {
        tmp[0]=verts[3*i+0]-ce[0];
        tmp[1]=verts[3*i+1]-ce[1];
        tmp[2]=verts[3*i+2]-ce[2];
        dynstereographic2(vec,tmp,&vertsproj[3*i]);
    }
}

-(void)spin:(char*)path nframes:(int)nframes
{
	int		i;
	
	[self display];
	for(i=0;i<nframes;i++)
	{
		NSAutoreleasePool	*pool=[NSAutoreleasePool new];
		[self savePicture:[NSString stringWithFormat:@"%s.%03i.jpg",path,i]];

		m_tbRot[0]+=360.0/(float)nframes;
		m_tbRot[1]=0;
		m_tbRot[2]=1;
		m_tbRot[3]=0;
		
		//zoom-=0.005;
		[self display];
		[pool drain];
	}
}
-(void)morphToMeshAtPath:(char*)mshpath destination:(char*)imgpath nframes:(int)nframes
{
	int		i,j,np;
    float3D *p,*p1,*p0,o={0,0,0};
    FILE    *f;
    char    str[255];
    float   t,g=1/pow(2,5.5);
    
    p=(float3D*)verts;
    
    // store original points
    p0=(float3D*)calloc(nverts,sizeof(float3D));
    for(i=0;i<nverts;i++)
        p0[i]=p[i];

    // load target mesh
    f=fopen(mshpath,"r");
    fgets(str,255,f);
    sscanf(str," %i %*i ",&np);
    p1=(float3D*)calloc(np,sizeof(float3D));
    for(i=0;i<np;i++)
    {
        fscanf(f," %f %f %f ",&(p1[i].x),&(p1[i].y),&(p1[i].z));
        o=add3D(o,p1[i]);
    }
    o=sca3D(o,1/(float)np);
    fclose(f);
    
    // turn &scale mesh
    /*
    for(i=0;i<np;i++)
        p1[i]=(float3D){g*(p1[i].x-o.x),g*(p1[i].y-o.y),-g*(p1[i].z-o.z)};
    */
	
	[self display];
	for(i=0;i<nframes;i++)
	{
		t=i/(float)(nframes-1);
        t=3*t*t-2*t*t*t;
        t=1-t;
        
        NSAutoreleasePool	*pool=[NSAutoreleasePool new];
		[self savePicture:[NSString stringWithFormat:@"%s.%03i.jpg",imgpath,i]];
        
		for(j=0;j<nverts;j++)
            p[j]=add3D(sca3D(p0[j],(1-t)),sca3D(p1[j],t));
        
		[self display];
		[pool drain];
	}
}
#pragma mark -
- (void)setZoom:(float)z
{
    zoom = pow(2,-z);
    [self setNeedsDisplay:YES];
}

#pragma mark -
-(IBAction)emSaveDepthMap:(id)sender
{
    NSRect	frame=[self bounds];
	int		W=frame.size.width;
	int		H=frame.size.height;
    NSBitmapImageRep *bmp=[[[NSBitmapImageRep alloc]
                                initWithBitmapDataPlanes:NULL
                                pixelsWide:W
                                pixelsHigh:H
                                bitsPerSample:32
                                samplesPerPixel:1
                                hasAlpha:NO
                                isPlanar:NO
                                colorSpaceName:NSCalibratedWhiteColorSpace
								bitmapFormat:NSFloatingPointSamplesBitmapFormat
                                bytesPerRow:0
                                bitsPerPixel:32] autorelease];
    float *tmp=(float*)calloc(W*H,sizeof(float));
    float *baseaddr=(float*)[bmp bitmapData];
    NSSavePanel *savePanel;
    int		result,i,j;
	
    [self getDepthMap:tmp width:W height:H];
	for(i=0;i<W;i++)
	for(j=0;j<H;j++)
		baseaddr[W*j+i]=tmp[W*(H-j-1)+i];
	free(tmp);
    
    savePanel = [NSSavePanel savePanel];
    
    [savePanel setAllowedFileTypes:[NSArray arrayWithObject:@"tif"]];
    [savePanel setCanSelectHiddenExtension:YES];
    result=[savePanel runModal];
    if (result == NSOKButton)
    {
        NSString *filename=[[savePanel URL] path];
        [[bmp TIFFRepresentation] writeToFile:filename atomically:YES];
    }
}
-(IBAction)emSaveImage:(id)sender
{
    NSRect	frame=[self bounds];
    NSBitmapImageRep *bmp=[[[NSBitmapImageRep alloc]
                                initWithBitmapDataPlanes:NULL
                                pixelsWide:frame.size.width
                                pixelsHigh:frame.size.height
                                bitsPerSample:8
                                samplesPerPixel:4
                                hasAlpha:YES
                                isPlanar:NO
                                colorSpaceName:NSCalibratedRGBColorSpace
                                bytesPerRow:0
                                bitsPerPixel:0] autorelease];
    NSImage *img;
    unsigned char *baseaddr=[bmp bitmapData];
    NSSavePanel *savePanel;
    int		result;
    
    [self getPixels:(char*)baseaddr width:frame.size.width height:frame.size.height rowbyte:[bmp bytesPerRow]];
    
    img = [[[NSImage alloc] init] autorelease];
    [img addRepresentation:bmp];
    [img setFlipped:YES];
    [img lockFocusOnRepresentation:bmp];
    [img unlockFocus];
    
    savePanel = [NSSavePanel savePanel];
    
    [savePanel setAllowedFileTypes:[NSArray arrayWithObject:@"tif"]];
    [savePanel setCanSelectHiddenExtension:YES];
    result=[savePanel runModal];
    if (result == NSOKButton)
    {
        NSString *filename=[[savePanel URL] path];
        [[img TIFFRepresentation] writeToFile:filename atomically:YES];
    }
}
- (void) savePicture:(NSString *)filename
{
    NSRect				bounds=[self bounds];
	int					i,j,W=bounds.size.width,H=bounds.size.height;
    NSData				*bmp2;
	NSBitmapImageRep	*bmp=[[NSBitmapImageRep alloc]
                                initWithBitmapDataPlanes:NULL
                                pixelsWide:W
                                pixelsHigh:H
                                bitsPerSample:8
                                samplesPerPixel:4
                                hasAlpha:YES
                                isPlanar:NO
                                colorSpaceName:NSDeviceRGBColorSpace
                                bytesPerRow:4*bounds.size.width
                                bitsPerPixel:0];
    unsigned char		*baseaddr=[bmp bitmapData],b[4];

    [self getPixels:(char*)baseaddr width:W height:H rowbyte:4*W];
	
	// flip
	for(i=0;i<bounds.size.width;i++)
		for(j=0;j<bounds.size.height/2;j++)
		{
			b[0]=baseaddr[4*(j*W+i)+0];
			b[1]=baseaddr[4*(j*W+i)+1];
			b[2]=baseaddr[4*(j*W+i)+2];
			b[3]=baseaddr[4*(j*W+i)+3];
			baseaddr[4*(j*W+i)+0]=baseaddr[4*((H-1-j)*W+i)+0];
			baseaddr[4*(j*W+i)+1]=baseaddr[4*((H-1-j)*W+i)+1];
			baseaddr[4*(j*W+i)+2]=baseaddr[4*((H-1-j)*W+i)+2];
			baseaddr[4*(j*W+i)+3]=baseaddr[4*((H-1-j)*W+i)+3];
			baseaddr[4*((H-1-j)*W+i)+0]=b[0];
			baseaddr[4*((H-1-j)*W+i)+1]=b[1];
			baseaddr[4*((H-1-j)*W+i)+2]=b[2];
			baseaddr[4*((H-1-j)*W+i)+3]=b[3];
		}
    
	bmp2=[bmp representationUsingType:NSJPEGFileType
						   properties:[NSDictionary dictionaryWithObject:[NSDecimalNumber numberWithFloat:0.8] forKey:NSImageCompressionFactor]];
	[bmp2 writeToFile:filename atomically:YES];
	[bmp release];
}

-(void) getPixels:(char*)baseaddr width:(long)w height:(long)h rowbyte:(long)rb
{
    glReadPixels(0,0,w,h,GL_RGBA,GL_UNSIGNED_BYTE,baseaddr);
}
-(void) getDepthMap:(float*)baseaddr width:(long)w height:(long)h
{
    glReadPixels(0,0,w,h,GL_DEPTH_COMPONENT,GL_FLOAT,baseaddr);
}
#pragma mark -
-(void)pickVertex:(NSPoint)mp
{
	NSRect	bounds=[self bounds];
	float	dist,min,aspectRatio;
	int		i,indx=-1;
	float3D	r;
	
	aspectRatio = (float)bounds.size.width/(float)bounds.size.height;
	mp=(NSPoint){aspectRatio*zoom*(2*mp.x-bounds.size.width)/bounds.size.width,zoom*(2*mp.y-bounds.size.height)/bounds.size.height};
	min=1000000;
	for(i=0;i<nverts;i++)
	{
		v_m((float*)&r,(float*)&verts[3*i],m);
		dist=pow(r.x-mp.x,2)+pow(r.y-mp.y,2);
		if(dist<min){ min=dist; indx=i;}
	}
	selected=indx;
	if(indx>-1)
	{
		if(vertsdata)
			[settings setValue:[NSString stringWithFormat:@"v%i(%.2f,%.2f,%.2f) = %.2f    v:%i, t:%i",
									indx,
									verts[3*indx+0],verts[3*indx+1],verts[3*indx+2],
									vertsdata[indx],
									nverts,ntris]
								forKey:@"nvertstris"];
		else
			[settings setValue:[NSString stringWithFormat:@"v%i(%.2f,%.2f,%.2f) = 0    v:%i, t:%i",
									indx,
									verts[3*indx+0],verts[3*indx+1],verts[3*indx+2],
									nverts,ntris]
								forKey:@"nvertstris"];
	}
    
    /*
	
	for(i=0;i<nverts;i++)
	{
		if(sqrt(pow(verts[3*selected+0]-verts[3*i+0],2)+pow(verts[3*selected+1]-verts[3*i+1],2)+pow(verts[3*selected+2]-verts[3*i+2],2))<20)
			vertscolour[3*i+0]=vertscolour[3*i+1]=vertscolour[3*i+2]=0.7;
		else
			vertscolour[3*i+0]=vertscolour[3*i+1]=vertscolour[3*i+2]=0;
	}
     
     */
    
    
	[self setNeedsDisplay:YES];
}
-(void)dragVertex:(NSPoint)mp
{
	NSRect	bounds=[self bounds];
	float	aspectRatio;
	float3D	tmp,tp;
	
	v_m((float*)&tp,(float*)&verts[3*selected],m);
	aspectRatio = (float)bounds.size.width/(float)bounds.size.height;
	tmp=(float3D){aspectRatio*zoom*(2*mp.x-bounds.size.width)/bounds.size.width,zoom*(2*mp.y-bounds.size.height)/bounds.size.height,tp.z};
	v_m((float*)&verts[3*selected],(float*)&tmp,invm);

	[self setNeedsDisplay:YES];
}
-(void)setSelected:(int)index
{
	selected=index;
}
-(int)selected
{
	return selected;
}
#pragma mark -
-(IBAction)emStandard:(id)sender
{
	int	tag=[[sender selectedCell] tag];
	[self setStandardRotation:tag];
}
-(IBAction)emProjection:(id)sender
{
	[self setProjection:[sender indexOfSelectedItem]];
	[self setNeedsDisplay:YES];
}
-(IBAction)emColourmap:(id)sender
{
	[self setColourmap:[sender indexOfSelectedItem]];
	[self setNeedsDisplay:YES];
}
-(IBAction)emRotate:(id)sender
{
	int	tag=[[sender selectedCell] tag];
	int	value=[sender intValue];
	
	[self addRotation:(value-0.5)*60 toAxis:tag];	
}
-(IBAction)emShading:(id)sender
{
	[self setShading:[sender indexOfSelectedItem]];
	[self setNeedsDisplay:YES];
}
-(IBAction)emZoom:(id)sender
{
	float	z=[sender floatValue];
	[self setZoom:z];
}

@end
