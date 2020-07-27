//
//  MyDocument.m
//  EditMesh
//
//  Created by rOBERTO tORO on 13/11/2005.
//  Copyright __MyCompanyName__ 2005 . All rights reserved.
//

#import "MyDocument.h"

@implementation MyDocument

- (id)init
{
    self = [super init];
    if (self) {
        M.p=0;
		M.t=0;
		M.np=0;
		M.nt=0;
    }
    return self;
}

- (NSString *)windowNibName
{
    // Override returning the nib file name of the document
    // If you need to use a subclass of NSWindowController or if your document supports multiple NSWindowControllers, you should remove this method and override -makeWindowControllers instead.
    return @"MyDocument";
}

-(void)configureMesh
{
	float  	S,V,logabsg;
	
	S=surface_area(M.p, M.t, M.nt);
	V=volume(M.p, M.t, M.nt);
	logabsg=log(S)-2*log(V)/3.0-log(36*pi)/3.0;
	
	[[settings content] setValue:[NSString stringWithFormat:@"v:%i, t:%i, Euler: %i, S=%.2f, V=%.2f, log(G)=%.2f",M.np,M.nt,M.np-M.nt/2,S,V,logabsg] forKey:@"nvertstris"];
	
	[view setVertices:(float*)M.p number:M.np];
	[view setTriangles:(int*)M.t number:M.nt];
	update(&M);

    float	*C=(float*)calloc(M.np,sizeof(float));
	em_textureDepth(&M, C);
	[view setVerticesData:C];
	
    [[settings content] setValue:[NSString stringWithFormat:@"%f,%f,%f\n",M.center.x,M.center.y,M.center.z] forKey:@"centre"];
	[view setNeedsDisplay:YES];
}
- (void)windowControllerDidLoadNib:(NSWindowController *) aController
{
    [super windowControllerDidLoadNib:aController];
	[text setApp:self];
	[view setSettings:[settings content]];
    if(M.p&&M.t)
	{
		[self configureMesh];
	}
}
- (NSData *)dataRepresentationOfType:(NSString *)aType
{
	NSData	*data=nil;
	char	*bytes;
	int		n;

	if([aType isEqualTo:@"OFFType"])
	{
		n=msh_packOFFText(&M, &bytes);
		data=[NSData dataWithBytes:bytes length:n];
	}
	else
	if([aType isEqualTo:@"3DMFType"])
	{
		n=msh_pack3DMF(&M, &bytes);
		data=[NSData dataWithBytes:bytes length:n];
	}
	else
	if([aType isEqualTo:@"TextMeshType"])
	{
		n=msh_packRawText(&M, &bytes);
		data=[NSData dataWithBytes:bytes length:n];
	}
	else
	if([aType isEqualTo:@"FSSurfType"])
	{
		n=msh_packFreesurferMeshData(&M, &bytes);
		data=[NSData dataWithBytes:bytes length:n];
	}
	else
	if([aType isEqualTo:@"VRMLType"])
	{
		n=msh_packVRML(&M, &bytes);
		data=[NSData dataWithBytes:bytes length:n];
	}
	else
	if([aType isEqualTo:@"PlyMeshType"])
	{
		n=msh_packPly(&M,(float3D*)[view verticesColour],&bytes);
		data=[NSData dataWithBytes:bytes length:n];
	}
    else
    if([aType isEqualTo:@"VTKMeshType"])
    {
        n=msh_packVTK(&M,(float3D*)[view verticesColour],&bytes);
        data=[NSData dataWithBytes:bytes length:n];
    }
	
    return data;
}

- (BOOL)readFromURL:(NSURL*)url ofType:(NSString *)aType error:(NSError**)outError
{
	if([aType isEqualTo:@"OFFType"])
		msh_readOFFText(&M, (char*)[[url path] UTF8String]);
	else
	if([aType isEqualTo:@"3DMFType"])
		msh_parse3DMFText(&M, (char*)[[NSData dataWithContentsOfURL:url] bytes]);
	else
	if([aType isEqualTo:@"TextMeshType"])
		msh_readRawText(&M, (char*)[[url path] UTF8String]);
	else
	if([aType isEqualTo:@"FSSurfType"])
		msh_readFreesurferMesh(&M, (char*)[[url path] UTF8String]);
	else
	if([aType isEqualTo:@"VRMLType"])
		msh_importVRMLMeshData(&M, (char*)[[url path] UTF8String]);
	else
	if([aType isEqualTo:@"BVMeshType"])
		msh_importBVMeshData(&M, (char*)[[url path] UTF8String]);
    else
    if([aType isEqualTo:@"PlyMeshType"])
        msh_importPlyMeshData(&M, (char*)[[url path] UTF8String]);
    else
    if([aType isEqualTo:@"VTKMeshType"])
        msh_importVTKMeshData(&M, (char*)[[url path] UTF8String]);
	
	return YES;
}
-(IBAction)exportVertexData:(id)sender
{
    NSSavePanel *savePanel;
    int		result;
	
    savePanel = [NSSavePanel savePanel];    
    result=[savePanel runModal];
    if (result == NSOKButton)
    {
        NSString *filename=[[savePanel URL] path];
        msh_exportFSTextureData([view verticesData], M.np, (char*)[filename UTF8String]);
    }
}
-(IBAction)importVertexData:(id)sender
{
	NSOpenPanel *open=[NSOpenPanel openPanel];
	int			result,i,imax,ddim;
	float		*tmp;
	float3D		*c;
	NSString	*path,*ext;
	NSData		*data;
	
    [open setAllowedFileTypes:[NSArray arrayWithObjects:@"annot",@"area",@"avg_curv",@"curv",@"sulc",@"target",@"thickness",@"float",@"inflated",@"txt",@"txt1",nil]];
	result=[open runModal];
    //result=[open runModalForDirectory:nil file:nil types:];
	if (result!=NSOKButton)
		return;
	
	path=[[[open URLs] objectAtIndex:0] path];
	ext=[path pathExtension];
	
	if([ext isEqualTo:@"annot"])
    {
        float3D *C=msh_getTexturePtr(&M);
        msh_importFSMeshAnnotation(C,M.np, (char*)[path UTF8String]);
        [view setVerticesColour:(float*)C];
    }
    else
    if([ext isEqualTo:@"float"])
	{
		data=[NSData dataWithContentsOfFile:path];
		msh_parseVerticesDataBin(&tmp, &ddim, M.np, (char *)[data bytes], [data length]);
		//for(i=0;i<M.np;i++) swapfloat(&(tmp[i]));
        
        if(ddim==1)
        {
            [view setVerticesData:tmp];
            c=msh_getTexturePtr(&M);
            for(i=0;i<M.np;i++)
                c[i].x=c[i].y=c[i].z=tmp[i];
            imax=0;
            for(i=0;i<M.np;i++)
                if(tmp[imax]<tmp[i])
                    imax=i;
            printf("Maximum value:%f at vertex index:%i\n",tmp[imax],imax);
            [view setSelected:imax];
        }
        else
        if(ddim==3)
        {
            float3D *C=msh_getTexturePtr(&M);
            for(i=0;i<M.np;i++)
                C[i]=((float3D*)tmp)[i];
            [view setVerticesColour:tmp];
            [view setNeedsDisplay:YES];
        }
	}
	else
	if([ext isEqualTo:@"txt"]||[ext isEqualTo:@"txt1"])
	{
		char	str[256];
		FILE	*f=fopen([path UTF8String],"r");
		int		np;
		
		fgets(str,256,f);
		sscanf(str," %i ", &np);
		tmp=(float*)calloc(np,sizeof(float));
		for(i=0;i<np;i++)
			fscanf(f," %f ",&tmp[i]);
		[view setVerticesData:tmp];
		c=msh_getTexturePtr(&M);
		imax=0;
		for(i=0;i<M.np;i++)
		{
			c[i].x=c[i].y=c[i].z=tmp[i];
			if(tmp[imax]<tmp[i])
				imax=i;
		}
		printf("Maximum value:%f at vertex index:%i\n",tmp[imax],imax);
		[view setSelected:imax];
	}
	else
	if([ext isEqualTo:@"inflated"])
	{
		c=msh_getTexturePtr(&M);
		msh_readFreesurferMesh(&M, (char*)[path UTF8String]);
		[view setVertices:(float*)M.p number:M.np];
		[view setTriangles:(int*)M.t number:M.nt];
		[view setVerticesColour:(float*)c];
	}
	else
	{
		msh_importFSTextureData(&tmp,M.np, (char*)[path UTF8String]);
		c=msh_getTexturePtr(&M);
		[view setVerticesData:tmp];

		c=msh_getTexturePtr(&M);
		imax=0;
		for(i=0;i<M.np;i++)
		{
			c[i].x=c[i].y=c[i].z=tmp[i];
			if(tmp[imax]<tmp[i])
				imax=i;
		}
		printf("Maximum value:%f at vertex index:%i\n",tmp[imax],imax);
		[view setSelected:imax];
	}

	[view setNeedsDisplay:YES];
}
//TEST
/*
-(IBAction)importVertexData:(id)sender
{
	NSOpenPanel *open=[NSOpenPanel openPanel];
	int			result;
	float3D		*c,*sulcaldepth=msh_getTexturePtr(&M);
	NSString	*path;
	float		val,*vol; // vol is a colin27-like volume, ie. float 181*217*181, colin's origin
	FILE		*f;
	int			i;
	
	result=[open runModalForDirectory:nil file:nil types:[NSArray arrayWithObjects:@"img",nil]];
	if (result!=NSOKButton)
		return;
	
	path=[[open filenames] objectAtIndex:0];
	f=fopen([path cString],"r");
	vol=(float*)calloc(181*217*181,sizeof(float));
	fread(vol,181*217*181,sizeof(float),f);
	fclose(f);
	c=(float3D*)calloc(M.np,sizeof(float3D));
	for(i=0;i<M.np;i++)
	{
		val=vol[((int)(M.p[i].z+72))*217*181+((int)(M.p[i].y+126))*181+(int)(M.p[i].x+90)];
		if(val>2.0)
			c[i].x=c[i].y=val;
			
		else
			c[i]=sulcaldepth[i];
	}

	[view setVerticesColour:(float*)c];
	[view setNeedsDisplay:YES];
}
*/

-(IBAction)importVertexParameterization:(id)sender
{
	NSOpenPanel *open=[NSOpenPanel openPanel];
	int			result;
	NSString	*path;
	FILE		*f;
	char		str[256];
	float		*param;
	
	[open setAllowedFileTypes:[NSArray arrayWithObjects:@"float3D",nil]];
    result=[open runModal];
	if (result!=NSOKButton)
		return;
	
	path=[[[open URLs] objectAtIndex:0] path];
	f=fopen([path UTF8String],"r");
	fgets(str,255,f);
	param=(float*)calloc(atoi(str)*3,sizeof(float));
	fread(param,atoi(str)*3,sizeof(float),f);
	fclose(f);

	[view setParam:(float*)param];
}
-(IBAction)importTextureImage:(id)sender
{
	NSOpenPanel			*open=[NSOpenPanel openPanel];
	int					result;
	NSString			*path;
	NSImage				*img;
	NSBitmapImageRep	*bmp;
	char				*data;
	int					i,j,k;
	
	[open setAllowedFileTypes:[NSArray arrayWithObjects:@"tif",@"png",@"psd",nil]];
    result=[open runModal];
	if (result!=NSOKButton)
		return;
	
	path=[[[open URLs] objectAtIndex:0] path];
	
	img=[[NSImage alloc] initWithContentsOfFile:path];
	bmp=[[img representations] objectAtIndex:0];
	data=calloc([bmp pixelsHigh]*[bmp pixelsWide]*[bmp samplesPerPixel],1);
	for(j=0;j<[bmp pixelsHigh];j++)
	for(i=0;i<[bmp pixelsWide];i++)
	for(k=0;k<[bmp samplesPerPixel];k++)
		data[(j*[bmp pixelsWide]+i)*[bmp samplesPerPixel]+k]=[bmp bitmapData][j*[bmp bytesPerRow]+[bmp samplesPerPixel]*i+k];

	[view setTexture:data size:(NSSize){[bmp pixelsWide],[bmp pixelsHigh]}];
	[view setNeedsDisplay:YES];
}
//
#pragma mark -
-(IBAction)emParam2d:(id)sender
{
	float3D	*C;

	// compute laplace coordinates
	em_textureCoordinatesFE(&M);
	
	// convert mesh to sphere
	C=msh_getTexturePtr(&M);
	if(1)
		em_sphereFromTxtr(&M,C);
	
	if(0)
	{
		int	i;
		float	*vc=[view verticesColour];
		for(i=0;i<M.np;i++)
		{
			vc[3*i+0]=C[i].x;
			vc[3*i+1]=C[i].y;
			vc[3*i+2]=C[i].z;
		}
	}
	if(0)
	{
		// save equivalent sphere points
		FILE *f=fopen("/Users/roberto/Desktop/sphere.txt","w");
		int	i;
		for(i=0;i<M.np;i++)
			fprintf(f,"%f %f %f\n",M.p[i].x,M.p[i].y,M.p[i].z);
		fclose(f);
	}

	[view setNeedsDisplay:YES];
}
-(IBAction)emSmooth:(id)sender
{
	em_smooth(&M);
	[view recomputeNormals:(float*)M.p];
	update(&M);
	[view setNeedsDisplay:YES];
}
-(IBAction)emCurvature:(id)sender
{
	float	*C=(float*)calloc(M.np,sizeof(float));
	int		i;

	em_textureMeanCurvature(&M,C);
	for(i=0;i<M.np;i++)
		C[i]=(C[i]-0.5)/0.5;

	[view setVerticesData:C];
}
-(IBAction)emDepth:(id)sender
{
	float	*C=(float*)calloc(M.np,sizeof(float));
	em_textureDepth(&M, C);
	[view setVerticesData:C];
	[view setNeedsDisplay:YES];
}
-(void)emIntegratedCurvature:(int)iter
{
	float	*C=(float*)calloc(M.np,sizeof(float));
	int		i;
	
	em_textureIntegratedMeanCurvature(&M,C,iter);
	for(i=0;i<M.np;i++)
		C[i]=(C[i]-0.5)/0.5;
	[view setVerticesData:C];
	[view setNeedsDisplay:YES];
}
-(IBAction)emRostrocaudal:(id)sender
{
	float	*d=[view verticesData];
	int		i,j,*n,N=100,x0,x1;
	float	*s;
	
	n=(int*)calloc(N,sizeof(int));
	s=(float*)calloc(N,sizeof(float));
	
	x0=x1=0;
	for(i=0;i<M.np;i++)
	{
		if(M.p[i].y<M.p[x0].y)	x0=i;
		if(M.p[i].y>M.p[x1].y)	x1=i;
	}
	
	for(i=0;i<M.np;i++)
	{
		j=(int)((N-1)*(M.p[i].y-M.p[x0].y)/(M.p[x1].y-M.p[x0].y));
		n[j]++;
		//s[j]+=d[i];
		s[j]=MAX(s[j],d[i]);
	}
	
	for(i=0;i<N;i++)
		//printf("%f ",s[i]/(float)n[i]);
		printf("%f ",s[i]);
	printf("\n");
	
	[view setSelected:x0];
	[view setNeedsDisplay:YES];
}
-(IBAction)emChangeCentre:(id)sender
{
	if(M.p==nil)
		return;
		
	float			x,y,z;
	unsigned char	*str;
	
	str=(unsigned char*)[[[settings content] objectForKey:@"centre"] UTF8String];
	sscanf((char*)str," %f , %f , %f ",&x,&y,&z); printf("x:%f, y:%f, z:%f\n",x,y,z);
	em_translate(&M,(float3D){x,y,z});
	[[settings content] setValue:[NSString stringWithFormat:@"%.1f,%.1f,%.1f\n",M.center.x,M.center.y,M.center.z] forKey:@"centre"];

	[view setNeedsDisplay:YES];
}
-(IBAction)emFlipTriangles:(id)sender
{
	em_flipTriangles(&M);
	[view setNeedsDisplay:YES];
}
-(IBAction)emChangeCMapMinMax:(id)sender
{
	[view setMinMaxData];
	[view setNeedsDisplay:YES];
}
-(void)maxVertex
{
	float	*d=[view verticesData];
	int		i,imax=0;
	
	for(i=0;i<M.np;i++)
		if(d[imax]<d[i])
			imax=i;
	[text insertText:[NSString stringWithFormat:@"\nvertex index: %i\nvertex value:%f\n",imax,d[imax]]];
}
-(void)valueAtVertex:(int)i
{
	float	*d=[view verticesData];
	[text insertText:[NSString stringWithFormat:@"\nvalue:%f, coords=(%f,%f,%f)\n",d[i],M.p[i].x,M.p[i].y,M.p[i].z]];
}
double	sum;
int		*tmark,icmax,ncverts;
-(void)averageClusterAtVertex:(int)ip threshold:(float)thr data:(float*)D nt:(NTriRec*)NT
{
	int		i,j,it;
	int		inside;
	int		*tt;
	
	tmark[ip]=1;
	
	ncverts++;
	
	if(D[ip]>D[icmax])
		icmax=ip;
	
	sum+=D[ip];

	for(i=0;i<=NT[ip].n;i++)
	{
		it=NT[ip].t[i];		
		tt=(int*)&(M.t[it]);
		inside=0;
		for(j=0;j<3;j++)
			if(D[tt[j]]>=thr && tmark[tt[j]]==0)
				[self averageClusterAtVertex:tt[j] threshold:thr data:D nt:NT];
	}
}
-(void)averageClusterValue:(float)thr
{
	float	*d=[view verticesData];
	int		n;
	int		i;
	NTriRec	*NT;
	
	NT=(NTriRec*)calloc(M.np,sizeof(NTriRec));		
	for(i=0;i<M.nt;i++)
	{
		NT[M.t[i].a].t[NT[M.t[i].a].n++] = i;
		NT[M.t[i].b].t[NT[M.t[i].b].n++] = i;
		NT[M.t[i].c].t[NT[M.t[i].c].n++] = i;
	}
	
	n=1;
	tmark=(int*)calloc(M.np,sizeof(int));
	for(i=0;i<M.np;i++)
	if(d[i]>=thr && tmark[i]==0)
	{
		icmax=i;
		ncverts=0;
		sum=0;
		[self averageClusterAtVertex:i threshold:thr data:d nt:NT];
		printf("cluster %i. max(vertex,value,coords)=(%i,%f,(%.2f,%.2f,%.2f)), avrg(vertices,value)=(%i,%f)\n",n,icmax,d[icmax],M.p[icmax].x,M.p[icmax].y,M.p[icmax].z,ncverts,sum/(float)ncverts);
		n++;
	}
		
	free(tmark);
	free(NT);
}
-(void)paintVertex:(float)x :(float)y :(float)z
{
	float	*vc=[view verticesColour];
	int		i=[view selected];
	
	vc[3*i+0]=x;
	vc[3*i+1]=y;
	vc[3*i+2]=z;
	
	[view setNeedsDisplay:YES];
}
-(void)foldLength
{
    
    int     nt=M.nt;
    float3D *p=M.p;
    int3D   *t=M.t;
    float   *data=(float*)calloc(M.np,sizeof(float));
    int     i,j;
    float   length=0,a,x;
    float3D p0[3];
    
    em_textureMeanCurvature(&M,data);
	for(i=0;i<M.np;i++)
		data[i]=(data[i]-0.5)/0.5;
    
    for(i=0;i<nt;i++)
    {
        j=0;
        if(data[t[i].a]*data[t[i].b]<0)
        {
            a=fabs(data[t[i].a]);
            x=a/(a+fabs(data[t[i].b]));
            p0[j++]=add3D(sca3D(p[t[i].a],1-x),sca3D(p[t[i].b],x));
        }
        if(data[t[i].b]*data[t[i].c]<0)
        {
            a=fabs(data[t[i].b]);
            x=a/(a+fabs(data[t[i].c]));
            p0[j++]=add3D(sca3D(p[t[i].b],1-x),sca3D(p[t[i].c],x));
        }
        if(data[t[i].c]*data[t[i].a]<0)
        {
            a=fabs(data[t[i].c]);
            x=a/(a+fabs(data[t[i].a]));
            p0[j++]=add3D(sca3D(p[t[i].c],1-x),sca3D(p[t[i].a],x));
        }
        if(j==2)
            length+=norm3D(sub3D(p0[0],p0[1]));
    }
    free(data);
    printf("foldLength: %f\n",length/2.0);
}
-(void)protrusions
{
	// get spherical parametrisation
	// vertex colour as a function of original/spherical edge length ratio
	float3D	*p0,*p1,*C;
	int3D	*t;
	MeshRec	M1;
	int	i;
	float	*num,*den,n[3],d[3];
	float	*D;
	
	M1=M;
	M1.p=(float3D*)calloc(M.np, sizeof(float3D));
	for(i=0;i<M1.np;i++)
		M1.p[i]=M.p[i];
	
	p0=msh_getPointsPtr(&M);
	p1=msh_getPointsPtr(&M1);
	t=msh_getTrianglesPtr(&M);

	// compute laplace coordinates
	em_textureCoordinatesFE(&M1);
	
	// convert mesh to sphere
	C=msh_getTexturePtr(&M);
	em_sphereFromTxtr(&M1,C);
	free(C);
	
	// compute deformation ratio
	num=(float*)calloc(M.np,sizeof(float));
	den=(float*)calloc(M.np,sizeof(float));
	for(i=0;i<M.nt;i++)
	{
		n[0]=norm3D(sub3D(p1[t[i].a],p1[t[i].b]));
		n[1]=norm3D(sub3D(p1[t[i].b],p1[t[i].c]));
		n[2]=norm3D(sub3D(p1[t[i].c],p1[t[i].a]));
		d[0]=norm3D(sub3D(p0[t[i].a],p0[t[i].b]));
		d[1]=norm3D(sub3D(p0[t[i].b],p0[t[i].c]));
		d[2]=norm3D(sub3D(p0[t[i].c],p0[t[i].a]));
		num[t[i].a]+=n[0]+n[2];
		num[t[i].b]+=n[1]+n[0];
		num[t[i].c]+=n[2]+n[1];
		den[t[i].a]+=d[0]+d[2];
		den[t[i].b]+=d[1]+d[0];
		den[t[i].c]+=d[2]+d[1];
	}
	free(M1.p);
	
	D=(float*)calloc(M.np,sizeof(float));
	for(i=0;i<M.np;i++)
		D[i]=den[i]/MAX(num[i],0.00001);
	free(num);
	free(den);
	[view setVerticesData:D];

	[view setNeedsDisplay:YES];
}
-(void)scale:(float)x
{
	float3D	*p0;
	int	i;
	
	p0=msh_getPointsPtr(&M);
	for(i=0;i<M.np;i++)
		p0[i]=sca3D(p0[i],x);
	
	[view setNeedsDisplay:YES];
}
-(void)deleteVertices:(float)thr
{
	int		i,sum,np,nt;
	float	*D=[view verticesData];
	float3D	*p=msh_getPointsPtr(&M);
	int3D	*t=msh_getTrianglesPtr(&M);
	int		*lu;
	
	lu=(int*)calloc(M.np,sizeof(int));
	np=0;
	for(i=0;i<M.np;i++)
	{
		if(D[i]>=thr)
			lu[i]=-1;
		else
			lu[i]=np++;
	}
	
	np=0;
	for(i=0;i<M.np;i++)
	if(lu[i]>=0)
		p[np++]=p[i];
	
	nt=0;
	for(i=0;i<M.nt;i++)
	{
		sum=0;
		sum+=lu[t[i].a]>=0;
		sum+=lu[t[i].b]>=0;
		sum+=lu[t[i].c]>=0;
		
		if(sum==3)
		{
			t[nt].a=lu[t[i].a];
			t[nt].b=lu[t[i].b];
			t[nt].c=lu[t[i].c];
			nt++;
		}
	}
	free(lu);

	M.np=np;
	M.nt=nt;
	[[settings content] setValue:[NSString stringWithFormat:@"v:%i, t:%i, Euler: %i",M.np,M.nt,M.np-M.nt/2] forKey:@"nvertstris"];
	[view setVertices:(float*)p number:np];
	[view setTriangles:(int*)t number:nt];
	update(&M);
	em_textureDepth(&M,D);
	[view setVerticesData:D];
	msh_setNeighborTriangles(&M);

	if(0)
	{
	for(i=0;i<np;i++)
		printf("%f %f %f\n",p[i].x,p[i].y,p[i].z);
	for(i=0;i<nt;i++)
		printf("%i %i %i\n",t[i].a,t[i].b,t[i].c);
	}

	[view setNeedsDisplay:YES];
}
/*
3 5 4 5		c1: e1.a<e2.a					r -1
4 5 3 5		c2: e1.a>e2.a					r  1
3 5 3 7		c3: e1.a==e2.a, e1.b<e2.b		r -1
3 7 3 5		c4: e1.a==e2.a, e1.b>e2.b		r  1
			c5: e1.a==e2.a, e1.b==e2.b		r  0
*/

int compareEdges (const void *a, const void *b)
{
	int2D	e1=*(int2D*)a;
	int2D	e2=*(int2D*)b;

	if(e1.a==e2.a)
	{
		if(e1.b==e2.b)
			return 0;
		else
		if(e1.b<e2.b)
			return -1;
		else
			return 1;
	}
	else
	{
		if(e1.a<e2.a)
			return -1;
		else
			return	1;
	}
}
void addTriangle(MeshRec *M,int3D *T, int *p1, int t1, int i1, int j1)
{
	if((T[t1].a==i1&&T[t1].b==j1)&&(T[t1].b==i1&&T[t1].c==j1)&&(T[t1].c==i1&&T[t1].a==j1))
		T[M->nt++]=(int3D){p1[j1],p1[i1],M->np};	// add triangle j1,i1,l0
	else
		T[M->nt++]=(int3D){p1[i1],p1[j1],M->np};	// add triangle i1,j1,l0
}
-(void)nonManifold
{
	int		i,j,n;
	float3D	*p=msh_getPointsPtr(&M);
	int		*t=(int*)msh_getTrianglesPtr(&M);
	int3D	*T,*e;
	int3D	*nme;
	int2D	*m;
	int		*p0,*p1,nn,s0,nl,sz;
	int		i1,j1,k,found,t0,t1;
	int		newp,newt;
	float3D	l0;
	
	T=(int3D*)t;

	// make a list of all edges
	e=(int3D*)calloc(M.nt*3,sizeof(int3D));
	for(i=0;i<M.nt;i++)
	for(j=0;j<3;j++)
	{
		if(t[3*i+j]<t[3*i+(j+1)%3])
		{
			e[3*i+j].a=t[3*i+j];
			e[3*i+j].b=t[3*i+(j+1)%3];
			e[3*i+j].c=i;
		}
		else
		{
			e[3*i+j].a=t[3*i+(j+1)%3];
			e[3*i+j].b=t[3*i+j];
			e[3*i+j].c=i;
		}
	}
	
	// sort edges
	qsort(e,M.nt*3,sizeof(int3D),compareEdges);
	
	// count nonmanifold edges
	i=0;
	j=0;
	n=0;
	do
	{
		if(e[j].a==e[j+1].a && e[j].b==e[j+1].b)
			j+=2;
		else
		{
			j++;
			n++;
		}
	}
	while(j<M.nt*3);
	printf("nonmanifold edges:%i\n",n);
	
	// store nonmanifold edges
	nme=(int3D*)calloc(n,sizeof(int3D));
	i=0;
	j=0;
	n=0;
	do
	{
		if(e[j].a==e[j+1].a && e[j].b==e[j+1].b)
			j+=2;
		else
		{
			nme[n++]=e[j];
			j++;
		}
	}
	while(j<M.nt*3);
	
	// recode vertices
	p0=(int*)calloc(M.np,sizeof(int));	// p0[#orig]=#new
	for(i=0;i<M.np;i++)
		p0[i]=-1;
	p1=(int*)calloc(n,sizeof(int));		// p1[#new]=#orig
	nn=0;
	for(i=0;i<n;i++)
	{
		if(p0[nme[i].a]<0)
		{
			p0[nme[i].a]=nn;
			p1[nn]=nme[i].a;
			nn++;
		}
		if(p0[nme[i].b]<0)
		{
			p0[nme[i].b]=nn;
			p1[nn]=nme[i].b;
			nn++;
		}
	}
	
	// make a non-manifold edge matrix
	m=(int2D*)calloc(n*n,sizeof(int2D));
	for(i=0;i<n;i++)
	{
		m[p0[nme[i].b]*n+p0[nme[i].a]]=(int2D){-1,i};
		m[p0[nme[i].a]*n+p0[nme[i].b]]=(int2D){-1,i};
	}
	
	// find loop intersections (if any)
	for(i=0;i<n;i++)
	{
		s0=0;
		for(j=0;j<n;j++)
			if(m[i*n+j].a<0)
				s0++;
		if(s0>2)
		{
			printf("loop intersection: %i\n",s0);
			for(j=0;j<n;j++)
				if(m[i*n+j].a<0)
					m[i*n+j].a=m[j*n+i].a=-2;
		}
	}
	
	// find and fill loops
	// WARNING: MESH MEMORY HAS TO BE EXPANDED BEFORE ADDING NEW TRIANGLES AND VERTICES!!
	newp=M.np;
	newt=M.nt;
	nl=1;
	for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	if(m[i*n+j].a==-1)
	{
		i1=i;
		j1=j;
		t0=nme[m[i*n+j].b].c;
		addTriangle(&M,T, p1, t0, i1, j1);
		l0=p[p1[i1]];
		sz=1;
		printf("loop %i: ",nl);
		do
		{
			m[i1*n+j1].a=nl;
			m[j1*n+i1].a=nl;
			found=0;
			for(k=0;k<n;k++)
			{
				if(m[i1*n+k].a==-1)
				{
					j1=k;
					t1=nme[m[i1*n+k].b].c;
					addTriangle(&M,T, p1, t1, i1, j1);
					l0=add3D(l0,p[p1[k]]);
					sz++;
					found=1;
					break;
				}
				if(m[k*n+j1].a==-1)
				{
					i1=k;
					t1=nme[m[k*n+j1].b].c;
					addTriangle(&M,T, p1, t1, i1, j1);
					l0=add3D(l0,p[p1[k]]);
					sz++;
					found=1;
					break;
				}
			}
			//if(found) printf(", %i(%i)",p1[k],t1);
		}
		while(found);
		l0=sca3D(l0,1/(float)sz);
		p[M.np++]=l0;
		printf(" %i vertices\n",sz);
		nl++;
	}
	free(e);
	
	printf("newp = %i\nnewt = %i\n",M.np-newp,M.nt-newt);
	printf("v:%i, t:%i, Euler: %i",M.np,M.nt,M.np-M.nt/2);

	[[settings content] setValue:[NSString stringWithFormat:@"v:%i, t:%i, Euler: %i",M.np,M.nt,M.np-M.nt/2] forKey:@"nvertstris"];
	[view setVertices:(float*)p number:M.np];
	[view setTriangles:(int*)t number:M.nt];
	update(&M);
	float	*D=(float*)calloc(M.np,sizeof(float));
	em_textureDepth(&M,D);
	[view setVerticesData:D];
	msh_setNeighborTriangles(&M);
}
-(void)emSmoothData
{
	float	*C;
	int		i;
	
	C=[view verticesData];
	for(i=0;i<100;i++)
		em_textureSmoothVerticesData(&M,C);
	[view setVerticesData:C];
	[view setNeedsDisplay:YES];
}
-(void)laplaceSmoothWithLambda:(float)lambda iterations:(int)N
{
	int	i;
	
	em_setSmooth(lambda);
	for(i=0;i<N;i++)
		em_smooth(&M);
	update(&M);
}
-(void)taubinSmoothWithLambda:(float)lambda mu:(float)mu iterations:(int)N
{
    printf(">>\n");
	int	i;
	
	for(i=0;i<N;i++)
		if(i%2==0)
		{
			em_setSmooth(lambda);
			em_smooth(&M);
		}
		else
		{
			em_setSmooth(mu);
			em_smooth(&M);
		}
	update(&M);
}
-(void)swapData
{
	float	*C;
	char 	*by;
	char	sw[4];
	int		i;

	C=[view verticesData];

	for(i=0;i<M.np;i++)
	{
		by=(char*)&(C[i]);
		sw[0]=by[3];
		sw[1]=by[2];
		sw[2]=by[1];
		sw[3]=by[0];
		C[i]=*(float*)sw;
	}
	[view setVerticesData:C];
}
-(void)adaptMesh
{
	em_adaptMesh(&M);

	float  	S,V,logabsg;
    float3D *C=msh_getTexturePtr(&M);
	
	S=surface_area(M.p, M.t, M.nt);
	V=volume(M.p, M.t, M.nt);
	logabsg=log(S)-2*log(V)/3.0-log(36*pi)/3.0;
	
	[[settings content] setValue:[NSString stringWithFormat:@"v:%i, t:%i, Euler: %i, S=%.2f, V=%.2f, log(G)=%.2f",M.np,M.nt,M.np-M.nt/2,S,V,logabsg] forKey:@"nvertstris"];
	
	[view setVertices:(float*)M.p number:M.np];
	[view setTriangles:(int*)M.t number:M.nt];
    [view setVerticesColour:(float*)C];
	update(&M);
	
    [[settings content] setValue:[NSString stringWithFormat:@"%f,%f,%f\n",M.center.x,M.center.y,M.center.z] forKey:@"centre"];
	[view setNeedsDisplay:YES];
}
#pragma mark -
#pragma mark Surface ratio test
double	sum;
char	*ttmark;

-(void)sumneighbours:(int)ip mesh:(MeshRec*)m nt:(NTriRec*)NT
{
	int		i,j,it;
	int		inside;
	int		*tt;
	float	d;
	double	R=20;
	
	for(i=0;i<=NT[ip].n;i++)
	{
		it=NT[ip].t[i];
		if(ttmark[it]==1)
			continue;
		
		tt=(int*)&((*m).t[it]);
		inside=0;
		for(j=0;j<3;j++)
		{
			d=norm3D(sub3D((*m).p[[view selected]],(*m).p[tt[j]]));
			if(d<R)
				inside++;
		}
		if(inside>1)
		{
			ttmark[it]=1;
			sum+=triArea((*m).p[tt[0]],(*m).p[tt[1]],(*m).p[tt[2]]);
			//printf("sum=%f\n",sum);
			for(j=0;j<3;j++)
				if(tt[j]!=ip)
					[self sumneighbours:tt[j] mesh:m nt:NT];
		}
	}
}
-(IBAction)emSurfaceRatio:(id)sender
{
	int		i;
	float	lf,R=20;
	MeshRec	m;
	float3D	*C;
	NTriRec	*NT;
	int		np,nt,n=0,*ip;
	float	*verts=[view verts];
	int		*tris=[view tris];
	float3D	*p;
	int3D	*t;
	FILE	*f;
	
	m.p=(float3D*)verts;
	m.t=(int3D*)tris;
	m.np=[view nverts];
	m.nt=[view ntris];

	NT=(NTriRec*)calloc(m.np,sizeof(NTriRec));		
	for(i=0;i<m.nt;i++)
	{
		NT[m.t[i].a].t[NT[m.t[i].a].n++] = i;
		NT[m.t[i].b].t[NT[m.t[i].b].n++] = i;
		NT[m.t[i].c].t[NT[m.t[i].c].n++] = i;
	}
	
	ttmark=(char*)calloc(m.nt,1);
	
	sum=0;
	[self sumneighbours:[view selected] mesh:&m nt:NT];
	lf=sum/(pi*R*R);

	C=(float3D*)calloc(m.np,sizeof(float3D));
	nt=0;
	for(i=0;i<m.nt;i++)
		if(ttmark[i])
		{
			C[m.t[i].a]=C[m.t[i].b]=C[m.t[i].c]=(float3D){1,0,0};
			nt++;
		}
		else
			C[m.t[i].a]=C[m.t[i].b]=C[m.t[i].c]=(float3D){0.5,0.5,0.5};
	[view setVerticesColour:(float*)C];
	free(C);
	free(NT);
	printf("vertex: %i, surface:%f, surface ratio:%f\n",[view selected],sum,lf);
	
	// save suface patch
	t=(int3D*)calloc(nt,sizeof(int3D));
	for(i=0;i<m.nt;i++)
	if(ttmark[i])
		t[n++]=*(int3D*)&tris[3*i];
	free(ttmark);
	
	ip=(int*)calloc(m.np,sizeof(int));
	for(i=0;i<m.np;i++)
		ip[i]=-1;
	np=0;
	for(i=0;i<nt;i++)
	{
		if(ip[t[i].a]==-1)
			ip[t[i].a]=np++;
		if(ip[t[i].b]==-1)
			ip[t[i].b]=np++;
		if(ip[t[i].c]==-1)
			ip[t[i].c]=np++;
		t[i].a=ip[t[i].a];
		t[i].b=ip[t[i].b];
		t[i].c=ip[t[i].c];
	}
	p=(float3D*)calloc(np,sizeof(float3D));
	for(i=0;i<m.np;i++)
	if(ip[i]>-1)
		p[ip[i]]=*(float3D*)&verts[3*i];
	
	f=fopen("/Users/roberto/Desktop/patch.txt","w");
	fprintf(f,"%i %i\n",np,nt);
	for(i=0;i<np;i++)
		fprintf(f,"%f,%f,%f\n",p[i].x,p[i].y,p[i].z);
	for(i=0;i<nt;i++)
		fprintf(f,"%i %i %i\n",t[i].a,t[i].b,t[i].c);
	fclose(f);

}
#pragma mark -
-(void)applyScript:(NSString*)s
{
	char		*cs,cmd[64];
	int			n;
	char		str[256];
	float		x,y,z;
	int			a;
	
	cs=(char*)[s UTF8String];
	n=0;
	while(cs[n]!=(char)0 && cs[n]!=' ' && cs[n]!='\r' && n<64)
    {
		cmd[n]=cs[n];
        n++;
    }
	cmd[n]=(char)0;
	printf("[%s]\n",cmd);

	if(strcmp(cmd,"savePicture")==0)
    {
        n=sscanf(cs," savePicture %s ",str);
        if(n==1)
            [view savePicture:[NSString stringWithUTF8String:str]];
    }
    else
    if(strcmp(cmd,"saveData")==0)
	{
		n=sscanf(cs," saveData %s ",str);
		if(n==1)
		{
			printf("write data to '%s'\n",str);
			float	*C=[view verticesData];
			FILE	*f=fopen(str,"w");
			fprintf(f,"%i\n",M.np);
			int	i;
			for(i=0;i<M.np;i++)
				fprintf(f,"%f\n",C[i]);
			fclose(f);
		}
	}
	else
	if(strcmp(cmd,"smoothMesh")==0)
		[self emSmooth:self];
	else
    if(strcmp(cmd,"adaptMesh")==0)
        [self adaptMesh];
	else
	if(strcmp(cmd,"smoothData")==0)
		[self emSmoothData];
	else
	if(strcmp(cmd,"nonManifold")==0)
		[self nonManifold];
	else
	if(strcmp(cmd,"depth")==0)
		[self emDepth:self];
	else
	if(strcmp(cmd,"curvature")==0)
		[self emCurvature:self];
	else
	if(strcmp(cmd,"icurvature")==0)
	{
		n=sscanf(cs," icurvature %i ",&a);
		if(n==1)
			[self emIntegratedCurvature:a];
	}
	else
	if(strcmp(cmd,"maxVertex")==0)
		[self maxVertex];
	else
	if(strcmp(cmd,"valueAtVertex")==0)
	{
		n=sscanf(cs," valueAtVertex %i ",&a);
		if(n==1)
			[self valueAtVertex:a];
	}
	else
    if(strcmp(cmd,"rostrocaudal")==0)
        [self emRostrocaudal:self];
    else
    if(strcmp(cmd,"foldLength")==0)
        [self foldLength];
    else
	if(strcmp(cmd,"protrusions")==0)
		[self protrusions];
	else
	if(strcmp(cmd,"swapData")==0)
		[self swapData];
	else
	if(strcmp(cmd,"averageClusterValue")==0)
	{
		n=sscanf(cs," averageClusterValue %f ",&x);
		if(n==1)
			[self averageClusterValue:x];
	}
	else
	if(strcmp(cmd,"laplaceSmooth")==0)
	{
		n=sscanf(cs," laplaceSmooth %f %i",&x,&a);	// lambda, iterations
		if(n==2)
			[self laplaceSmoothWithLambda:x iterations:a];
	}
	else
	if(strcmp(cmd,"taubinSmooth")==0)
	{
		n=sscanf(cs," taubinSmooth %f %f %i",&x,&y,&a);	// lambda, mu, iterations
		if(n==3)
			[self taubinSmoothWithLambda:x mu:y iterations:a];
	}
	else
	if(strcmp(cmd,"deleteVertices")==0)
	{
		n=sscanf(cs," deleteVertices %f ",&x);
		if(n==1)
			[self deleteVertices:x];
	}
	else
	if(strcmp(cmd,"paintVertex")==0)
	{
		n=sscanf(cs," paintVertex %f , %f , %f ",&x,&y,&z);
		if(n==3)
			[self paintVertex:x:y:z];
	}
	else
	if(strcmp(cmd,"scale")==0)
	{
		n=sscanf(cs," scale %f ",&x);
		if(n==1)
			[self scale:x];
	}
	else
	if(strcmp(cmd,"spin")==0)
	{
		n=sscanf(cs," spin %s %f",str,&x);
		if(n==2)
			[view spin:str nframes:x];
	}
	else
    if(strcmp(cmd,"morph")==0)
    {
        char    str1[512];
        n=sscanf(cs," morph %s %s %f ",str,str1,&x);
        if(n==3)
            [view morphToMeshAtPath:str destination:str1 nframes:x];
    }
    else
    if(strcmp(cmd,"addMesh")==0)
    {
        n=sscanf(cs," addMesh %s ",str);
        if(n==1)
            [self addMesh:str];
    }
    else
    if(strcmp(cmd,"applyRotation")==0)
    {
        [self applyRotation];
    }
    else
    if(strcmp(cmd,"saveStereographic")==0)
    {
        n=sscanf(cs," saveStereographic %s ",str);
        if(n==1)
            [self saveStereographic:str];
    }
	/*
	if(strcmp(cmd,"setCentreMesh")==0)
	{
		n=sscanf(cs," setCentreMesh %s ",str);
		if(n==1)
		{
			if(strcmp(str,"on")==0)
				[view setCentreMesh:1];
			else
			if(strcmp(str,"off")==0)
				[view setCentreMesh:0];
		}
	}
	if(strcmp(cmd,"setHemisphere")==0)
	{
		n=sscanf(cs," setHemisphere %s ",str);
		if(n==1)
			[self sendMessage:str];
	}
	*/
	[view setNeedsDisplay:YES];
}
-(void)emAddMesh:(id)sender
{
	NSOpenPanel *open=[NSOpenPanel openPanel];
	int			result;
	NSString	*path;
	NSString    *ext;
    MeshRec     Mtmp;
	
	result=[open runModal];
	if (result!=NSOKButton)
		return;
	path=[[[open URLs] objectAtIndex:0] path];
    ext=[path pathExtension];

    if([ext isEqualToString:@"txt"])
    {
        msh_readRawText(&Mtmp, (char*)[path UTF8String]);
    }
    else
    if([ext isEqualToString:@"inflated"] ||
       [ext isEqualToString:@"orig"] ||
       [ext isEqualToString:@"pial"] ||
       [ext isEqualToString:@"reg"] ||
       [ext isEqualToString:@"smoothwm"] ||
       [ext isEqualToString:@"sphere"] ||
       [ext isEqualToString:@"white"])
    {
        msh_readFreesurferMesh(&Mtmp, (char*)[path UTF8String]);
    }
    else
    {
        printf("ERROR: can't read this file type\n");
        return;
    }

	int        i;
    float3D	*newverts;
	int3D	*newtris;
    
	newverts=(float3D*)calloc(M.np+Mtmp.np,sizeof(float3D));
	newtris=(int3D*)calloc(M.nt+Mtmp.nt,sizeof(int3D));
	for(i=0;i<M.np;i++)
		newverts[i]=M.p[i];
	for(i=0;i<M.nt;i++)
		newtris[i]=M.t[i];
	for(i=0;i<Mtmp.np;i++)
		newverts[i+M.np]=Mtmp.p[i];
	for(i=0;i<Mtmp.nt;i++)
		newtris[i+M.nt]=(int3D){Mtmp.t[i].a+M.np,Mtmp.t[i].b+M.np,Mtmp.t[i].c+M.np};
	free(M.p);
	free(M.t);
	
	M.p=newverts;
	M.t=newtris;
	M.np+=Mtmp.np;
	M.nt+=Mtmp.nt;
    
    msh_dispose(&Mtmp);
	
	[self configureMesh];
}
-(void)applyRotation
{
    float   mat[16];
    float   rmat[9],imat[9];
    int		i;
    float3D p,q;
    
    [view getRotationMatrix:mat];
    rmat[0]=mat[0]; rmat[1]=mat[1]; rmat[2]=mat[2];
    rmat[3]=mat[4]; rmat[4]=mat[5]; rmat[5]=mat[6];
    rmat[6]=mat[8]; rmat[7]=mat[9]; rmat[8]=mat[10];
    invMat(rmat,imat);
    printf("%g %g %g\n%g %g %g\n%g %g %g\n",
           imat[0],imat[1],imat[2],
           imat[3],imat[4],imat[5],
           imat[6],imat[7],imat[8]);
    for(i=0;i<M.np;i++)
    {
        p=M.p[i];
        q.x=imat[0]*p.x+imat[1]*p.y+imat[2]*p.z;
        q.y=imat[3]*p.x+imat[4]*p.y+imat[5]*p.z;
        q.z=imat[6]*p.x+imat[7]*p.y+imat[8]*p.z;
        M.p[i]=q;
    }
}
-(void)addMesh:(char*)path
{
    FILE	*f;
	int		i,np,nt;
	float3D	*p1;
	int3D	*t1;
	char	s[256];
	float3D	*newverts;
	int3D	*newtris;
	
	// read new mesh
	f=fopen(path,"r");
	if(f==nil)
		return;
	fgets(s,256,f);
	sscanf(s," %i %i ",&np,&nt);
	p1=(float3D*)calloc(np,sizeof(float3D));
	t1=(int3D*)calloc(nt,sizeof(int3D));
	for(i=0;i<np;i++)
	{
		fgets(s,256,f);
		sscanf(s," %f %f %f ",&(p1[i]).x,&(p1[i]).y,&(p1[i]).z);
	}
	for(i=0;i<nt;i++)
	{
		fgets(s,256,f);
		sscanf(s," %i %i %i ",&t1[i].a,&t1[i].b,&t1[i].c);
	}
	fclose(f);
	
	newverts=(float3D*)calloc(M.np+np,sizeof(float3D));
	newtris=(int3D*)calloc(M.nt+nt,sizeof(int3D));
	for(i=0;i<M.np;i++)
		newverts[i]=M.p[i];
	for(i=0;i<M.nt;i++)
		newtris[i]=M.t[i];
	for(i=0;i<np;i++)
		newverts[i+M.np]=p1[i];
	for(i=0;i<nt;i++)
		newtris[i+M.nt]=(int3D){t1[i].a+M.np,t1[i].b+M.np,t1[i].c+M.np};
	free(M.p);
	free(M.t);
	
	M.p=newverts;
	M.t=newtris;
	M.np+=np;
	M.nt+=nt;
	
	[self configureMesh];
}
-(void)saveStereographic:(char*)path
{
    printf("[saveStereographic]\n");
    FILE      *f;
    float3D   *p,o={0,0,0},v;
    int3D     *t;
    int       *T;
    int       i,j,np,nt,flag,n;
    float     R=0;
	
    //[view projection2];
    [view projection];

    p=(float3D*)[view vertsproj];
    t=(int3D*)[view tris];
    np=[view nverts];
    nt=[view ntris];
    
    // find centre
    for(i=0;i<np;i++)
        o=add3D(o,p[i]);
    o=sca3D(o,1/(float)np);
    for(i=0;i<np;i++)
        if(R<norm3D(sub3D(p[i],o)))
            R=norm3D(sub3D(p[i],o));
    
    // count number of triangles to conserve
    n=0;
    for(i=0;i<nt;i++)
    {
        T=(int*)&(t[i]);
        flag=0;
        for(j=0;j<3;j++)
        {
            v=sub3D(p[T[j]],o);
            if(norm3D(v)>=R*0.8)    // save only triangles within 80% of the radius
                flag=1;
        }
        if(flag==0)
            n++;
    }
    
    // save mesh
    f=fopen(path,"w");
    fprintf(f,"%i %i\n",np,n);
    for(i=0;i<np;i++)
    {
        v=sub3D(p[i],o);
        fprintf(f,"%f %f %f\n",v.x,v.y,v.z);
    }
    for(i=0;i<nt;i++)
    {
        T=(int*)&(t[i]);
        flag=0;
        for(j=0;j<3;j++)
        {
            v=sub3D(p[T[j]],o);
            if(norm3D(v)>=R*0.8)    // save only triangles within 80% of the radius
                flag=1;
        }
        if(flag==0)
            fprintf(f,"%i %i %i\n",t[i].a,t[i].b,t[i].c);
    }
    fclose(f);
}
@end
