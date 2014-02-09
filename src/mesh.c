
/*
 *  mesh.c
 *  EditMesh
 *
 *  Created by roberto on Sun Sep 14 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include "util.h"
#include <string.h>

// IO globals
// Handle	g_msh_txtH;
// long		g_msh_txtSize;

void swapint(int *n)
{
	char 	*by=(char*)n;
	char	sw[4]={by[3],by[2],by[1],by[0]};
	
	*n=*(int*)sw;
}
void swapfloat(float *n)
{
	char 	*by=(char*)n;
	char	sw[4]={by[3],by[2],by[1],by[0]};
	
	*n=*(float*)sw;
}
void msh_readRawText(MeshRec *mesh, char *path)
{
	FILE		*f;
	int			i;
	char		str[256];

	f=fopen(path,"r");
	fgets(str,255,f);
	sscanf(str,"%i %i",&(*mesh).np,&(*mesh).nt);
	
	(*mesh).p=(float3D*)calloc((*mesh).np,sizeof(float3D));
	(*mesh).t=(int3D*)calloc((*mesh).nt,sizeof(int3D));
	for(i=0;i<(*mesh).np;i++)
	{	
		fgets(str,255,f);
		sscanf(str,"%f %f %f",&(*mesh).p[i].x,&(*mesh).p[i].y,&(*mesh).p[i].z);
	}
	for(i=0;i<(*mesh).nt;i++)
	{	
		fgets(str,255,f);
		sscanf(str,"%i %i %i",&(*mesh).t[i].a,&(*mesh).t[i].b,&(*mesh).t[i].c);
	}
	
	mesh->xcontainer=0;
    
    msh_getNeighborTrianglesPtr(mesh);

    update(mesh);
}
int msh_packRawText(MeshRec *m, char **data)
{
    char	str[255];
    int		i,sum;
    float3D	*p;
    int3D	*T;
    
    p=msh_getPointsPtr(m);
    T=msh_getTrianglesPtr(m);
    
    sum=0;
    sum+=sprintf(str,"%i %i\n",(*m).np,(*m).nt);
    for(i=0;i<(*m).np;i++)
        sum+=sprintf(str,"%f %f %f\n",p[i].x,p[i].y,p[i].z);
    for(i=0;i<(*m).nt;i++)
        sum+=sprintf(str,"%i %i %i\n",T[i].a,T[i].b,T[i].c);
    
    *data=calloc(sum,1);
    sum=0;
    sum+=sprintf(&(*data)[sum],"%i %i\n",(*m).np,(*m).nt);
    for(i=0;i<(*m).np;i++)
        sum+=sprintf(&(*data)[sum],"%f %f %f\n",p[i].x,p[i].y,p[i].z);
    for(i=0;i<(*m).nt;i++)
        sum+=sprintf(&(*data)[sum],"%i %i %i\n",T[i].a,T[i].b,T[i].c);

    return sum;
}
void msh_readOFFText(MeshRec *mesh, char *path)
{
	int			i;
	FILE		*f;
	char		str[512];
	
	f=fopen(path,"r");

	fgets(str,512,f);	// skip OFF
	fscanf(f," %i %i %*i ",&(*mesh).np,&(*mesh).nt);
	(*mesh).p=(float3D*)calloc((*mesh).np,sizeof(float3D));
	(*mesh).t=(int3D*)calloc((*mesh).nt,sizeof(int3D));
	for(i=0;i<(*mesh).np;i++)
		fscanf(f," %f %f %f ",&(*mesh).p[i].x,&(*mesh).p[i].y,&(*mesh).p[i].z);
	for(i=0;i<(*mesh).nt;i++)
		fscanf(f," %*i %i %i %i ",&(*mesh).t[i].a,&(*mesh).t[i].b,&(*mesh).t[i].c);
	
	(*mesh).xcontainer=0;
    
    msh_getNeighborTrianglesPtr(mesh);

    update(mesh);
}
int msh_packOFFText(MeshRec *m, char **data)
{
    char	str[255];
    int		i,sum;
    float3D	*p;
    int3D	*T;
    
    p=msh_getPointsPtr(m);
    T=msh_getTrianglesPtr(m);
    
    sum=0;
    sum+=sprintf(str,"%i %i 0\n",(*m).np,(*m).nt);
    for(i=0;i<(*m).np;i++)
        sum+=sprintf(str,"%f %f %f\n",p[i].x,p[i].y,p[i].z);
    for(i=0;i<(*m).nt;i++)
        sum+=sprintf(str,"3 %i %i %i\n",T[i].a,T[i].b,T[i].c);
    
    *data=calloc(sum,1);
    sum=0;
    sum+=sprintf(&(*data)[sum],"%i %i 0\n",(*m).np,(*m).nt);
    for(i=0;i<(*m).np;i++)
        sum+=sprintf(&(*data)[sum],"%f %f %f\n",p[i].x,p[i].y,p[i].z);
    for(i=0;i<(*m).nt;i++)
        sum+=sprintf(&(*data)[sum],"3 %i %i %i\n",T[i].a,T[i].b,T[i].c);

    return sum;
}
void msh_parse3DMFText(MeshRec *mesh, char *data)
{
    int	i,j;
    int3D	*T;
    float3D	*p;
    char	stop;

    i=0;

    // npoints
    getnumfromtxt(&i,data,&stop);	// jump 3
    getnumfromtxt(&i,data,&stop);	// jump 1
    getnumfromtxt(&i,data,&stop);	// jump 5
    getnumfromtxt(&i,data,&stop);	// jump 0
    
	(*mesh).np = getnumfromtxt(&i,data,&stop);
    (*mesh).p = (float3D*)calloc((*mesh).np,sizeof(float3D));
    if((*mesh).p==NULL) printf("Not enough memory: mesh points\n");
    p = (*mesh).p;	
    
    // points
    for(j=0;j<(*mesh).np;j++)
    {
            p[j].x= getnumfromtxt(&i,data,&stop);
            p[j].y= getnumfromtxt(&i,data,&stop);
            p[j].z= getnumfromtxt(&i,data,&stop);
    }
    
    // ntrian
    (*mesh).nt = getnumfromtxt(&i,data,&stop);
    (*mesh).t = (int3DPtr)calloc((*mesh).nt,sizeof(int3D));
    if((*mesh).t==NULL) printf("Not enough memory: mesh triangles\n");
    T = (*mesh).t;
    
    // ncontour
    getnumfromtxt(&i,data,&stop);
    
    // triangles
    for(j=0;j<(*mesh).nt;j++)
    {
            getnumfromtxt(&i,data,&stop);
            T[j].a= getnumfromtxt(&i,data,&stop);
            T[j].b= getnumfromtxt(&i,data,&stop);
            T[j].c= getnumfromtxt(&i,data,&stop);
    }
	
	(*mesh).xcontainer=0;
    
    msh_getNeighborTrianglesPtr(mesh);

    update(mesh);
}
int msh_pack3DMF(MeshRec *m, char **data)
{
    char	str[255];
    int		i,sum;
    float3D	*p;
    int3D	*T;
    
    p=msh_getPointsPtr(m);
    T=msh_getTrianglesPtr(m);
    
    sum=0;
    sum+=sprintf(str,"3DMetafile ( 1 5 norm3Dal tableofcontents0> )\nMesh\n(\n");
    sum+=sprintf(str,"%i\n",(*m).np);
    for(i=0;i<(*m).np;i++)
        sum+=sprintf(str,"%f %f %f\n",p[i].x,p[i].y,p[i].z);
    sum+=sprintf(str,"%i\n",(*m).nt);
    sum+=sprintf(str,"%i\n",0);
    for(i=0;i<(*m).nt;i++)
        sum+=sprintf(str,"3 %i %i %i\n",T[i].a,T[i].b,T[i].c);
    sum+=sprintf(str,")\n");
    
    *data=calloc(sum,1);
    sum=0;
    sum+=sprintf(&(*data)[sum],"3DMetafile ( 1 5 norm3Dal tableofcontents0> )\nMesh\n(\n");
    sum+=sprintf(&(*data)[sum],"%i\n",(*m).np);
    for(i=0;i<(*m).np;i++)
        sum+=sprintf(&(*data)[sum],"%f %f %f\n",p[i].x,p[i].y,p[i].z);
    sum+=sprintf(&(*data)[sum],"%i\n",(*m).nt);
    sum+=sprintf(&(*data)[sum],"%i\n",0);
    for(i=0;i<(*m).nt;i++)
        sum+=sprintf(&(*data)[sum],"3 %i %i %i\n",T[i].a,T[i].b,T[i].c);
    sum+=sprintf(&(*data)[sum],")\n");

    return sum;
}

void msh_parseVerticesDataBin(float **vdat, int *ddim, int np, char *data, int size)
{
    // OUT	vdat:	pointer to received data
    // IN	nd:	data dimension
    // 		np:	data points
    // 		data:	raw data
    // 		size:	size of raw data

    int		i,j=0,n;
    char	b4[4],swb4[4],str[1024];
    int		np_file,nd_file,version_file,version=1;
	float	max,min;

    while(data[j]!='\n'&&data[j]!='\r')
        str[j++]=data[j];
    str[j++]=(char)0;
    n=sscanf(str," %i %i %i ",&np_file, &nd_file, &version_file);
    *ddim=1;
    if(n==1)
        version=1;
    else
    if(n==2)
        version=2;
    else
    if(n==3)
    {
        if(version_file!=3)
        {
            printf("ERROR: look more closely at this data file: it looks like version 3 but it's not...");
            exit(1);
        }
        version=3;
        *ddim=nd_file;
    }
    if(np_file!=np)
        printf("error: texture dimensions\n");
    
    *vdat=(float*)calloc(*ddim*np,sizeof(float));
    for(i=0;i<*ddim*np;i++)
    {
        if(version<2)
        {
            *(float*)b4 = ((float*)(data+j))[i];
            swb4[0]=b4[3];
            swb4[1]=b4[2];
            swb4[2]=b4[1];
            swb4[3]=b4[0];
            (*vdat)[i] = *((float*)(swb4));
        }
        else
            (*vdat)[i]=((float*)(data+j))[i];
        
        if(i==0)                max=min=(*vdat)[i];
        else if((*vdat)[i]<min) min=(*vdat)[i];
        else if((*vdat)[i]>max) max=(*vdat)[i];
        
        if(!((*vdat)[i]==(*vdat)[i]))
            printf("Nan\n");
    }
    printf("[min,max]=[%f,%f]\n",min,max);        
}

int msh_packVerticesData(float *vdat, int dim, int np, char **data)
{
    char	str[255];
    int		sum;
    
    sum=0;
    sum+=sprintf(str,"%i\n",np);
    sum+=dim*sizeof(float)*np;
    
    *data=calloc(sum,1);
    sum=0;
    sum+=sprintf(&(*data)[sum],"%i\n",np);
    memcpy(&(*data)[sum],vdat,dim*sizeof(float)*np);
    sum+=dim*sizeof(float)*np;

    return sum;
}

int msh_readFreesurferMesh(MeshRec *mesh, char *path)
{
    FILE	*f;
	int		i,j;
    int3D	*T;
    float3D	*p;
    int		id,a,b,c,d;
    char	date[256],info[256];
	char	byte12[12];

	char testInt[]={1,0,0,0};
	int endianness=*(int*)testInt;

    f=fopen(path,"r");

	// read triangle/quad identifier: 3 bytes
    a=(unsigned char)fgetc(f);
	a=a<<16;
    b=(unsigned char)fgetc(f);
	b=b<<8;
    c=(unsigned char)fgetc(f);
    id=a+b+c;
    
    if(id==16777214)	// triangle mesh
    {
		printf("FS id (16777214 triangle) %i\n",id);
        // get creation date text line
        j=0;
		do
		{
			date[j]=fgetc(f);
		}
        while(date[j++]!=(char)10);
        date[j-1]=(char)0;
		printf("FS date %s\n",date);
        // get info text line
        j=0;
		do
		{
			info[j]=fgetc(f);
		}
        while(info[j++]!=(char)10);
		info[j-1]=(char)0;
		printf("FS info %s\n",info);
        
        // get number of vertices
        a=((int)(u_int8_t)fgetc(f))<<24;
        b=((int)(u_int8_t)fgetc(f))<<16;
        c=((int)(u_int8_t)fgetc(f))<<8;
        d=(u_int8_t)fgetc(f);
        (*mesh).np=a+b+c+d;
		printf("FS #points %i\n",(*mesh).np);
        (*mesh).p = (float3D*)calloc((*mesh).np,sizeof(float3D));
        if((*mesh).p==0) printf("Not enough memory: mesh points\n");
        p = (*mesh).p;
        	
        // get number of triangles
        a=((int)(u_int8_t)fgetc(f))<<24;
        b=((int)(u_int8_t)fgetc(f))<<16;
        c=((int)(u_int8_t)fgetc(f))<<8;
        d=(u_int8_t)fgetc(f);
        (*mesh).nt=a+b+c+d;
		printf("FS #triangles %i\n",(*mesh).nt);
        (*mesh).t = (int3D*)calloc((*mesh).nt,sizeof(int3D));
        if((*mesh).t==0) printf("Not enough memory: mesh triangles\n");
        T = (*mesh).t;
       
		if(endianness==1) // Intel endianness
		{
			for(j=0;j<(*mesh).np;j++)
			{
				for(i=0;i<12;i++) byte12[i/4*4+3-i%4]=fgetc(f);
				p[j] = *(float3D*)byte12;
			}
			for(j=0;j<(*mesh).nt;j++)
			{
				for(i=0;i<12;i++) byte12[i/4*4+3-i%4]=fgetc(f);
				T[j] = *(int3D*)byte12;
				T[j]=(int3D){T[j].a,T[j].c,T[j].b}; // flip triangle
			}
		}
		else				// Motorola endianness
		{
			// read vertices
			fread(p,mesh->np,sizeof(float3D),f);

			// read triangles
			for(j=0;j<(*mesh).nt;j++)
			{
				fread(&T[j],1,sizeof(int3D),f);
				T[j]=(int3D){T[j].a,T[j].c,T[j].b}; // flip triangle
			}
		}
        
        mesh->xcontainer=0;
        
        msh_getNeighborTrianglesPtr(mesh);

        update(mesh);
    }
	printf("FSSurf finished\n");
	return 1;
}
int msh_packFreesurferMeshData(MeshRec *m, char **data)
{
    int		sum,i,j,k;
    int		id=16777214,a,b,c;
    char	date[256]="EMPTY",info[256]="EMPTY";
    float	*p;
    int		*T;
	int		np,nt;

	char testInt[]={1,0,0,0};
	int endianness=*(int*)testInt;
    
    p=(float*)msh_getPointsPtr(m);
    T=(int*)msh_getTrianglesPtr(m);
	np=(*m).np;
	nt=(*m).nt;

	printf("[FreeSurfer_save_mesh]\n");

	sum=3+6+6;
	sum+=sizeof(int)+sizeof(int);
	sum+=3*np*sizeof(float)+3*nt*sizeof(int);
	
    *data=calloc(sum,1);

	i=0;
    a=id>>16;
    b=(id&0xff00)>>8;
    c=(id&0xff);
	(*data)[i++]=(char)a;
	(*data)[i++]=(char)b;
	(*data)[i++]=(char)c;

	date[strlen(date)]=(char)10;
	info[strlen(info)]=(char)10;
	
	memcpy(&((*data)[i]),date,6);
	i+=strlen(date);
	memcpy(&((*data)[i]),info,6);
	i+=strlen(info);
	
	if(endianness==1)	// Intel
	{
		for(j=0;j<4;j++)
			(*data)[i+j]=((char*)&np)[3-j];
		i+=sizeof(int);
		for(j=0;j<4;j++)
			(*data)[i+j]=((char*)&nt)[3-j];
		i+=sizeof(int);
		for(k=0;k<3*np;k++)
			for(j=0;j<4;j++)
				(*data)[i+k*4+j]=((char*)&(p[k]))[3-j];
		i+=3*np*sizeof(float);
		for(k=0;k<3*nt;k++)
			for(j=0;j<4;j++)
				(*data)[i+k*4+j]=((char*)&(T[k]))[3-j];
		i+=3*nt*sizeof(int);
	}
	else				// Motorola
	{
		memcpy(&((*data)[i]),&np,sizeof(int));
		i+=sizeof(int);
		memcpy(&((*data)[i]),&nt,sizeof(int));
		i+=sizeof(int);
		memcpy(&((*data)[i]),p,3*np*sizeof(float));
		i+=3*np*sizeof(float);
		memcpy(&((*data)[i]),T,3*nt*sizeof(int));
		i+=3*nt*sizeof(int);
	}

	return sum;
}
int msh_exportFSTextureData(float *d, int np, char *path)
{
	FILE	*f;
    int		id=16777215,a,b,c,FaceCount=0,ValsPerVertex=1,i,n;
    float	x;

    f=fopen(path,"w");
	
	if(f==NULL)
		return 1;
	
	// write data identifier: 3 bytes
    a=id>>16;
    b=(id&0xff00)>>8;
    c=(id&0xff);
    fputc((char)a,f);
    fputc((char)b,f);
    fputc((char)c,f);
    
	n=np;
	if(1)//(endianness==kINTEL)
	{
		swapint(&n);
		swapint(&FaceCount);
		swapint(&ValsPerVertex);
	}

    // write number of vertices
	fwrite(&n,1,sizeof(int),f);
	
	// write empty FaceCount and ValsPerVertex
	fwrite(&FaceCount,1,sizeof(int),f);
	fwrite(&ValsPerVertex,1,sizeof(int),f);

	if(1)//(endianness==kINTEL)
	{
		for(i=0;i<np;i++)
		{
			x=d[i];
			swapfloat(&x);
			fwrite(&x,1,sizeof(float),f);
		}
	}
	else
		fwrite(d,np,sizeof(float),f);
	fclose(f);

	return 0;
}
void msh_importFSTextureData(float **dat, int np, char *path)
{
    FILE	*f;
	int		i,j;
    float	v;
    int		id,a,b,c,n,nan=0;
	char	byte4[4];

	char testInt[]={1,0,0,0};
	int endianness=*(int*)testInt;

    f=fopen(path,"r");

	// read triangle/quad identifier: 3 bytes
    a=((int)(u_int8_t)fgetc(f))<<16;
    b=((int)(u_int8_t)fgetc(f))<<8;
    c=(u_int8_t)fgetc(f);
    id=a+b+c;
    
    if(id==16777215)	// triangle mesh
    {
		printf("FS id (16777215 vertex data) %i\n",id);
        // get number of vertices
		if(endianness==1)	// Intel
			for(i=0;i<4;i++) byte4[3-i]=fgetc(f);
		else
			fread(byte4,4,sizeof(char),f);
		n=*(int*)byte4;
		printf("FS #vertex_data %i\n",n);
		if(n!=np)
			printf("ERROR: data vertices (%i) different from mesh vertices (%i)\n",n,np);
        
		*dat=(float*)calloc(np,sizeof(float));
		
		// disregard FaceCount and ValsPerVertex
		fgetc(f);fgetc(f);fgetc(f);fgetc(f);
		fgetc(f);fgetc(f);fgetc(f);fgetc(f);
		
        // read vertex data
        for(j=0;j<n;j++)
		{
			if(endianness==1)	// Intel
				for(i=0;i<4;i++) byte4[3-i]=fgetc(f);
			else
				fread(byte4,4,sizeof(char),f);
            (*dat)[j]=*(float*)byte4;

			if(!(v==v))
				nan++;
		}
    }
	if(nan)
		printf("ERROR: Nan count:%i\n",nan);
	printf("FSVertex finished\n");
	
	fclose(f);
}
int msh_importFSMeshAnnotation(float **dat, int np, char *path)
{
    FILE	*f;
	int		i,l;
    int		n;
	char	*tmp;
    
	char testInt[]={1,0,0,0};
	int endianness=*(int*)testInt;
    
    f=fopen(path,"r");
	if(f==NULL)
		return 0;
	
	fread(&n,1,sizeof(int),f);	if(endianness==1) swapint(&n);
    if(n!=np) printf("Annotation file corrupted. points:%i annotations:%i [msh_importFSMeshAnnotation]\n",np,n);
    else
    {
        tmp=calloc(np,2*sizeof(int));
        *dat=(float*)calloc(3*np,sizeof(float));
        if(tmp==NULL) printf("Cannot allocate memory for tmp [msh_importFSMeshAnnotation]\n");
        else
        {
            fread(tmp,n,2*sizeof(int),f);
            for(i=0;i<MIN(np,n);i++)
            {
                l=((int*)tmp)[2*i+1]; if(endianness==1) swapint(&l);
                (*dat)[3*i+0]=(l&0xff);
                (*dat)[3*i+1]=((l>>8)&0xff);
                (*dat)[3*i+2]=((l>>16)&0xff);
            }
        }
        free(tmp);
    }
	fclose(f);

	return 1;
}
void msh_importBVMeshData(MeshRec *mesh, char *path)
{
    FILE	*f;
	char	tmp[6];
    int		i;
    int		endian,ignore;
    int3D	*T;
    float3D	*p;
    
	if(mesh==NULL)
        mesh=(MeshRec*)calloc(1,sizeof(MeshRec));

    f=fopen(path,"r");

	// get format (ascii, binar)
    for(i=0;i<5;i++) tmp[i]=fgetc(f);
    tmp[5]=(char)0;
    
    if(strcmp(tmp,"ascii")==0)
    {
    	fscanf(f," %*s ");					// ignore VOID
    	fscanf(f," %*i %*i %*i ");			// ignore 3 integers
    	
    	// READ 3-D COORDINATES
    	fscanf(f," %i ",&((*mesh).np));
		(*mesh).p = (float3D*)calloc((*mesh).np,sizeof(float3D));
    	for(i=0;i<(*mesh).np;i++)
    		fscanf(f," ( %f , %f , %f ) ", &(*mesh).p[i].x,&(*mesh).p[i].y,&(*mesh).p[i].z);
    	
    	fscanf(f," %i ",&ignore);			// ignore number of normal vectors
		if(ignore==(*mesh).np)
			for(i=0;i<(*mesh).np;i++)		// ignore normal vectors
				fscanf(f," ( %*f , %*f , %*f ) ");
    	fscanf(f," %*i ");					// ignore an integer
    	
    	// READ TRIANGLES
    	fscanf(f," %i ",&(*mesh).nt);
		(*mesh).t = (int3D*)calloc((*mesh).nt,sizeof(int3D));
    	for(i=0;i<(*mesh).nt;i++)
    		fscanf(f," ( %i , %i , %i ) ", &(*mesh).t[i].a,&(*mesh).t[i].b,&(*mesh).t[i].c);
    }
    else
    if(strcmp(tmp,"binar")==0)
    {
        for(i=0;i<4;i++) tmp[i]=fgetc(f);
        tmp[4]=(char)0;
        
        endian=-1;
        if(strcmp(tmp,"ABCD")==0)	endian=0;	
        if(strcmp(tmp,"DCBA")==0)	endian=1;
        if(endian==-1){ printf("not ABCD nor DCBA order...exit.\n"); return;}

        for(i=0;i<4;i++) tmp[i]=fgetc(f);	// ignore "VOID" string length
        for(i=0;i<4;i++) tmp[i]=fgetc(f);	// ignore "VOID" string

		ignore=fgetint(f,endian);			// verify number of vertices per polygon
        if(ignore!=3){ printf("Only able to read triangle meshes. This mesh has %i vertices per polygon.\n",ignore); return;}
        fgetint(f,endian);			// ignore time steps
        fgetint(f,endian);			// ignore time step index
    
        (*mesh).np=fgetint(f,endian);		// read number of polygons
		(*mesh).p = (float3DPtr)calloc((*mesh).np,sizeof(float3D));
		if((*mesh).p==NULL) printf("Not enough memory: mesh points\n");
		p = (*mesh).p;	
        printf("%i vertices\n", (*mesh).np);
        for(i=0;i<(*mesh).np;i++)
        {
            p[i].x=fgetfloat(f,endian);
            p[i].y=fgetfloat(f,endian);
            p[i].z=fgetfloat(f,endian);
        }
        
        ignore=fgetint(f,endian);		// ignore number of norm3Dal vectors
		if(ignore==(*mesh).np)
			for(i=0;i<(*mesh).np;i++)	// ignore norm3Dal vectors
			{
				fgetfloat(f,endian);
				fgetfloat(f,endian);
				fgetfloat(f,endian);
			}
        
        fgetint(f,endian);			// ignore number of texture coordinates
        
        (*mesh).nt=fgetint(f,endian);		// read number of triangles
		(*mesh).t = (int3DPtr)calloc((*mesh).nt,sizeof(int3D));
		if((*mesh).t==NULL) printf("Not enough memory: mesh triangles\n");
		T = (*mesh).t;
        printf("%i triangles\n", (*mesh).nt);
        for(i=0;i<(*mesh).nt;i++)		// read triangles
        {
            T[i].a=fgetint(f,endian);
            T[i].b=fgetint(f,endian);
            T[i].c=fgetint(f,endian);
        }
    }
	fclose(f);
    
    return;
}
/*
{
    FILE	*f;
	char	type[6],word[256];
    int		i,ii=0;
    int		endian,ignore;
	MeshPtr m=*mesh;
    int3D	*T;
    float3D	*p;
    char	stop;
    
	if(m==NULL)
        *mesh = m = (MeshPtr)calloc(1,sizeof(MeshRec));

    f=fopen(path,"r");

	// get format (ascii, binar)
    for(i=0;i<5;i++) type[i]=data[ii++];
    type[5]=(char)0;
    
    if(strcmp(type,"ascii")==0)
    {
        // read points
		getwordfromtxt(&ii, data, word);   // ignore "VOID"
		getnumfromtxt(&ii, data, &stop);
		getnumfromtxt(&ii, data, &stop);
		getnumfromtxt(&ii, data, &stop);
		(*m).np=getnumfromtxt(&ii, data, &stop);
		(*m).p = (float3DPtr)calloc((*m).np,sizeof(float3D));
		if((*m).p==NULL) printf("Not enough memory: mesh points\n");
		p = (*m).p;	
        printf("%i vertices\n", (*m).np);
        for(i=0;i<(*m).np;i++)
        {
			p[i].x=getnumfromtxt(&ii, data, &stop);
			p[i].y=getnumfromtxt(&ii, data, &stop);
			p[i].z=getnumfromtxt(&ii, data, &stop);
        }
        
        // ignore norm3Dals
        getnumfromtxt(&ii, data, &stop);
        for(i=0;i<(*m).np;i++)
		{
			getnumfromtxt(&ii, data, &stop);
			getnumfromtxt(&ii, data, &stop);
			getnumfromtxt(&ii, data, &stop);
		}
        
        // ignore texture
		getnumfromtxt(&ii, data, &stop);
        
        // read triangles
        (*m).nt=getnumfromtxt(&ii, data, &stop);
		(*m).t = (int3DPtr)calloc((*m).nt,sizeof(int3D));
		if((*m).t==NULL) printf("Not enough memory: mesh triangles\n");
		T = (*m).t;
        printf("%i triangles\n",(*m).nt);
        for(i=0;i<(*m).nt;i++)
        {
			T[i].a=getnumfromtxt(&ii, data, &stop);
			T[i].b=getnumfromtxt(&ii, data, &stop);
			T[i].c=getnumfromtxt(&ii, data, &stop);
        }
    }
    else
    if(strcmp(type,"binar")==0)
    {
        for(i=0;i<4;i++) type[i]=data[ii++];
        type[4]=(char)0;
        
        endian=-1;
        if(strcmp(type,"ABCD")==0)	endian=0;	
        if(strcmp(type,"DCBA")==0)	endian=1;
        if(endian==-1){ printf("not ABCD nor DCBA order...exit.\n"); return;}

        for(i=0;i<4;i++) data[ii++];		// ignore "VOID" string length
        ii+=4;								// ignore "VOID"

        ignore=getint(&ii,data,endian);	// verify vertices in polygon
        if(ignore!=3){ printf("only able to read triangle meshes...exit.\n"); return;}
        getint(&ii,data,endian);			// ignore time steps
        getint(&ii,data,endian);			// ignore time step index
    
        (*m).np=getint(&ii,data,endian);		// read number of polygons
		(*m).p = (float3DPtr)calloc((*m).np,sizeof(float3D));
		if((*m).p==NULL) printf("Not enough memory: mesh points\n");
		p = (*m).p;	
        printf("%i vertices\n", (*m).np);
        for(i=0;i<(*m).np;i++)
        {
            p[i].x=getfloat(&ii,data,endian);
            p[i].y=getfloat(&ii,data,endian);
            p[i].z=getfloat(&ii,data,endian);
        }
        
        getint(&ii,data,endian);		// ignore number of norm3Dal vectors
        for(i=0;i<(*m).np;i++)			// ignore norm3Dal vectors
        {
            getfloat(&ii,data,endian);
            getfloat(&ii,data,endian);
            getfloat(&ii,data,endian);
        }
        
        getint(&ii,data,endian);		// ignore number of texture coordinates
        
        (*m).nt=getint(&ii,data,endian);		// read number of triangles
		(*m).t = (int3DPtr)calloc((*m).nt,sizeof(int3D));
		if((*m).t==NULL) printf("Not enough memory: mesh triangles\n");
		T = (*m).t;
        printf("%i triangles\n", (*m).nt);
        for(i=0;i<(*m).nt;i++)		// read triangles
        {
            T[i].a=getint(&ii,data,endian);
            T[i].b=getint(&ii,data,endian);
            T[i].c=getint(&ii,data,endian);
        }
    }
    
    return;
}
*/
int msh_importVRMLMeshData(MeshRec *mesh, char *path)
{
    FILE	*f;
	int		i,n,loop;
	char	str[256],*tmp;

    f=fopen(path,"r");

	(*mesh).np=0;
	(*mesh).nt=0;
	
	loop=1;
	while(loop)
	{
		fgets(str,255,f);
		if(strstr(str,"point"))
			loop=0;
	}
	loop=1;
	while(loop)
	{
		fgets(str,255,f);
		if(strstr(str,"]"))
			loop=0;
		else
		if(strstr(str,"[")==NULL)
		{
			tmp=str;
			while(tmp=strstr(tmp+1,","),tmp)
			{
				(*mesh).np++;
				tmp++;
			}
		}
	}
	(*mesh).np++;
	loop=1;
	while(loop)
	{
		fgets(str,255,f);
		if(strstr(str,"coordIndex"))
			loop=0;
	}
	loop=1;
	while(loop)
	{
		fgets(str,255,f);
		if(strstr(str,"]"))
			loop=0;
		else
		{
			tmp=str;
			while(tmp=strstr(tmp+1,"-1"),tmp)
			{
				(*mesh).nt++;
				tmp+=2;
			}
		}
	}

	(*mesh).p = (float3D*)calloc((*mesh).np,sizeof(float3D));
	(*mesh).t = (int3D*)calloc((*mesh).nt,sizeof(int3D));
	fseek(f,0,SEEK_SET);
   
	loop=1;
	while(loop)
	{
		fgets(str,255,f);
		if(strstr(str,"point"))
			loop=0;
	}
	loop=1;
	i=0;
	while(loop)
	{
		fgets(str,255,f);
		if(strstr(str,"]"))
			loop=0;
		else
		if(strstr(str,"[")==NULL)
		{
			tmp=str;
			do
			{
				n=sscanf(tmp," %f %f %f ",&(*mesh).p[i].x,&(*mesh).p[i].y,&(*mesh).p[i].z);
				tmp=strstr(tmp,",");
				if(tmp)
					tmp++;
				if(n>0)
					i++;
			}
			while(tmp);
		}
	}
	loop=1;
	while(loop)
	{
		fgets(str,255,f);
		if(strstr(str,"coordIndex"))
			loop=0;
	}
	loop=1;
	i=0;
	while(loop)
	{
		fgets(str,255,f);
		if(strstr(str,"]"))
			loop=0;
		else
		if(strstr(str,"[")==NULL)
		{
			tmp=str;
			do
			{
				n=sscanf(tmp," %i%*[ ,\t]%i%*[ ,\t]%i ",&(*mesh).t[i].a,&(*mesh).t[i].c,&(*mesh).t[i].b);
				tmp=strstr(tmp,"-1,");
				if(tmp)
					tmp+=3;
				if(n>0)
					i++;
			}
			while(tmp);
		}
	}
	fclose(f);

	msh_getNeighborTrianglesPtr(mesh);

	update(mesh);
	printf("VRML finished\n");
	return 1;
}
int msh_packVRML(MeshRec *m, char **data)
{
    char	str[255];
    int		i,sum;
    float3D	*p;
    int3D	*T;
    
    p=msh_getPointsPtr(m);
    T=msh_getTrianglesPtr(m);
    
    sum=0;
    sum+=sprintf(str,"#VRML V1.0 ascii\n");
    sum+=sprintf(str,"Separator {\n");
    sum+=sprintf(str,"ShapeHints {\n");
    sum+=sprintf(str,"vertexOrdering COUNTERCLOCKWISE\n");
    sum+=sprintf(str,"faceType CONVEX\n");
    sum+=sprintf(str,"}\n");
    sum+=sprintf(str,"Coordinate3 {\n");
    sum+=sprintf(str,"point [\n");
    for(i=0;i<(*m).np;i++)
        sum+=sprintf(str,"%f %f %f,\n",p[i].x,p[i].y,p[i].z);
	sum+=sprintf(str,"]\n");
	sum+=sprintf(str,"}\n");
	sum+=sprintf(str,"IndexedFaceSet {\n");
	sum+=sprintf(str,"coordIndex [\n");
    for(i=0;i<(*m).nt;i++)
        sum+=sprintf(str,"%i,%i,%i,-1\n",T[i].a,T[i].b,T[i].c);
	sum+=sprintf(str,"]\n");
	sum+=sprintf(str,"}\n");
	sum+=sprintf(str,"}\n");
    
    *data=calloc(sum,1);
    sum=0;
    sum+=sprintf(&(*data)[sum],"#VRML V1.0 ascii\n");
    sum+=sprintf(&(*data)[sum],"Separator {\n");
    sum+=sprintf(&(*data)[sum],"ShapeHints {\n");
    sum+=sprintf(&(*data)[sum],"vertexOrdering COUNTERCLOCKWISE\n");
    sum+=sprintf(&(*data)[sum],"faceType CONVEX\n");
    sum+=sprintf(&(*data)[sum],"}\n");
    sum+=sprintf(&(*data)[sum],"Coordinate3 {\n");
    sum+=sprintf(&(*data)[sum],"point [\n");
    for(i=0;i<(*m).np;i++)
        sum+=sprintf(&(*data)[sum],"%f %f %f,\n",p[i].x,p[i].y,p[i].z);
	sum+=sprintf(&(*data)[sum],"]\n");
	sum+=sprintf(&(*data)[sum],"}\n");
	sum+=sprintf(&(*data)[sum],"IndexedFaceSet {\n");
	sum+=sprintf(&(*data)[sum],"coordIndex [\n");
    for(i=0;i<(*m).nt;i++)
        sum+=sprintf(&(*data)[sum],"%i,%i,%i,-1\n",T[i].a,T[i].b,T[i].c);
	sum+=sprintf(&(*data)[sum],"]\n");
	sum+=sprintf(&(*data)[sum],"}\n");
	sum+=sprintf(&(*data)[sum],"}\n");

    return sum;
}
int msh_importPlyMeshData(MeshRec *mesh, char *path)
{
    FILE	*f;
    int		i,x;
    char	str[512],str1[256],str2[256];
    
    f=fopen(path,"r");
	if(f==NULL){printf("ERROR: Cannot open file\n");return 1;}
	
	// READ HEADER
	mesh->np=mesh->nt=0;
	do
	{
		fgets(str,511,f);
		sscanf(str," %s %s %i ",str1,str2,&x);
		if(strcmp(str1,"element")==0&&strcmp(str2,"vertex")==0)
			mesh->np=x;
		else
		if(strcmp(str1,"element")==0&&strcmp(str2,"face")==0)
			mesh->nt=x;
	}
	while(strcmp(str1,"end_header")!=0 && !feof(f));
	if(mesh->np*mesh->nt==0)
	{
		printf("ERROR: Bad Ply file header format\n");
		return 1;
	}
	(*mesh).p = (float3D*)calloc((*mesh).np,sizeof(float3D));
	(*mesh).t = (int3D*)calloc((*mesh).nt,sizeof(int3D));
	// READ VERTICES
	if(mesh->p==NULL){printf("ERROR: Not enough memory for mesh vertices\n");return 1;}
	for(i=0;i<mesh->np;i++)
	{
		fgets(str,512,f);
		sscanf(str," %f %f %f ",&((*mesh).p[i].x),&((*mesh).p[i].y),&((*mesh).p[i].z));
	}
	printf("Read %i vertices\n",mesh->np);
	
	// READ TRIANGLES
	if(mesh->t==NULL){printf("ERROR: Not enough memory for mesh triangles\n"); return 1;}
	for(i=0;i<mesh->nt;i++)
		fscanf(f," 3 %i %i %i ",&((*mesh).t[i].a),&((*mesh).t[i].b),&((*mesh).t[i].c));
	printf("Read %i triangles\n",mesh->nt);
	
	fclose(f);
    
    return 0;
}
int msh_packPly(MeshRec *m,float3D *C,char **data)
{
    char	str[255];
    int		i,sum;
    float3D	*p;
    int3D	*t;
    
    p=msh_getPointsPtr(m);
    t=msh_getTrianglesPtr(m);
    
    sum=0;
	// WRITE HEADER
	sum+=sprintf(str,"ply\n");
	sum+=sprintf(str,"format ascii 1.0\n");
	sum+=sprintf(str,"comment meshconvert, R. Toro 2010\n");
	sum+=sprintf(str,"element vertex %i\n",m->np);
	sum+=sprintf(str,"property float x\n");
	sum+=sprintf(str,"property float y\n");
	sum+=sprintf(str,"property float z\n");
	sum+=sprintf(str,"element face %i\n",m->nt);
	sum+=sprintf(str,"property list uchar int vertex_indices\n");
	sum+=sprintf(str,"end_header\n");
	// WRITE VERTICES
	for(i=0;i<m->np;i++)
		//sum+=sprintf(str,"%f %f %f %i %i %i\n",p[i].x,p[i].y,p[i].z,(int)(255*C[i].x),(int)(255*C[i].y),(int)(255*C[i].z));	
		sum+=sprintf(str,"%f %f %f \n",p[i].x,p[i].y,p[i].z);	
	// WRITE TRIANGLES
	for(i=0;i<m->nt;i++)
		sum+=sprintf(str,"3 %i %i %i\n",t[i].a,t[i].b,t[i].c);
    
    *data=calloc(sum,1);
    sum=0;
	sum+=sprintf(&(*data)[sum],"ply\n");
	sum+=sprintf(&(*data)[sum],"format ascii 1.0\n");
	sum+=sprintf(&(*data)[sum],"comment meshconvert, R. Toro 2010\n");
	sum+=sprintf(&(*data)[sum],"element vertex %i\n",m->np);
	sum+=sprintf(&(*data)[sum],"property float x\n");
	sum+=sprintf(&(*data)[sum],"property float y\n");
	sum+=sprintf(&(*data)[sum],"property float z\n");
	sum+=sprintf(&(*data)[sum],"element face %i\n",m->nt);
	sum+=sprintf(&(*data)[sum],"property list uchar int vertex_indices\n");
	sum+=sprintf(&(*data)[sum],"end_header\n");
	// WRITE VERTICES
	for(i=0;i<m->np;i++)
		//sum+=sprintf(&(*data)[sum],"%f %f %f %i %i %i\n",p[i].x,p[i].y,p[i].z,(int)(255*C[i].x),(int)(255*C[i].y),(int)(255*C[i].z));
		sum+=sprintf(&(*data)[sum],"%f %f %f\n",p[i].x,p[i].y,p[i].z);	
	// WRITE TRIANGLES
	for(i=0;i<m->nt;i++)
		sum+=sprintf(&(*data)[sum],"3 %i %i %i\n",t[i].a,t[i].b,t[i].c);
	
    return sum;
}

#pragma mark -
//
// mesh functions
//
bool msh_new(MeshRec **m, int np, int nt)
{
    bool	isGood = true;
    int		i;
    char	name[]=" untitled mesh";
    
    *m = (MeshPtr)calloc(1,sizeof(MeshRec));
    if(*m==NULL)	isGood=false;

    (**m).np = np;
    (**m).nt = nt;
    (**m).p = (float3DPtr)calloc(np,sizeof(float3D));
    if((**m).p==NULL)	isGood=false;
    (**m).t = (int3DPtr)calloc(nt,sizeof(int3D));
    if((**m).t==NULL)	isGood=false;
    
    (**m).center = (float3D){0,0,0};
    (**m).bbox = (rect3D){{0,0,0},{0,0,0}};
    name[0]=(char)14;
    for(i=0;i<=name[0];i++)
        (**m).name[i]=name[i];
    
    (**m).xcontainer=0;
    
    return(isGood);
}

void msh_dispose(MeshRec *m)
{
    ContainerPtr	xcH;
    long		temp;

    free((*m).p);
    free((*m).t);
    
    if((*m).xcontainer)
    {
        xcH=(ContainerPtr)(*m).xcontainer;
        while((long)xcH)
        {
            temp=(*xcH).xcontainer;

            free((*xcH).data);
            free(xcH);
            
            if(temp)	xcH=(ContainerPtr)temp;
            else		xcH=NULL;
        }
    }
}
bool msh_copy(MeshRec *src, MeshRec **dst)
{
	int	i;
	bool	isGood = true;
	
	*dst = (MeshPtr)calloc(1,sizeof(MeshRec));
	if(*dst==NULL)	isGood=false;

	(**dst).np = (*src).np;
	(**dst).nt = (*src).nt;
	(**dst).p = (float3DPtr)calloc((*src).np,sizeof(float3D));
	if((**dst).p==NULL)	isGood=false;
	(**dst).t = (int3DPtr)calloc((*src).nt,sizeof(int3D));
	if((**dst).t==NULL)	isGood=false;
	
	for(i=0;i<(*src).np;i++)
		(**dst).p[i]=(*src).p[i];
	for(i=0;i<(*src).nt;i++)
		(**dst).t[i]=(*src).t[i];
		
	(**dst).center	=(*src).center;
	(**dst).bbox	= (*src).bbox;
        for(i=0;i<=(*src).name[0];i++)
            (**dst).name[i]=(*src).name[i];
	
	(**dst).xcontainer=0;
	
	return(isGood);
}
bool msh_newMeshCurve(MeshCurveRec *MC,int np,int ne)
{
	bool	isGood = true;
        char	name[]=" untitled meshCurve";
        int	i;
	
	(*MC).p = (MeshPointPtr)calloc(np,sizeof(MeshPointRec));
	if((*MC).p==NULL)	isGood=false;

	(*MC).e = (int2DPtr)calloc(ne,sizeof(int2D));
	if((*MC).e==NULL)	isGood=false;

	(*MC).c = (float3DPtr)calloc(np,sizeof(float3D));
	if((*MC).c==NULL)	isGood=false;
	
	(*MC).np = np;
	(*MC).ne = ne;
	
	name[0]=18;
        for(i=0;i<=name[0];i++)
            (*MC).name[i]=name[i];
	
	return(isGood);
}
void msh_disposeMeshCurve(MeshCurveRec MC)
{
	free(MC.p);
        free(MC.e);
	free(MC.c);
}
bool msh_newMeshEdgeCurve(MeshEdgeCurveRec *MEC,int np,int ne)
{
	bool	isGood = true;
        char	name[]=" untitled MeshEdgeCurve";
        int	i;
	
	(*MEC).p = (MeshEdgePointPtr)calloc(np,sizeof(MeshEdgePointRec));
	if((*MEC).p==NULL)	isGood=false;

	(*MEC).e = (int2DPtr)calloc(ne,sizeof(int2D));
	if((*MEC).e==NULL)	isGood=false;

	(*MEC).c = (float3DPtr)calloc(np,sizeof(float3D));
	if((*MEC).c==NULL)	isGood=false;
	
	(*MEC).np = np;
	(*MEC).ne = ne;
	
	name[0]=22;
        for(i=0;i<=name[0];i++)
            (*MEC).name[i]=name[i];
	
	return(isGood);
}
void msh_disposeMeshEdgeCurve(MeshEdgeCurveRec MEC)
{
	free(MEC.p);
	free(MEC.e);
	free(MEC.c);
}

#pragma mark _
void msh_setContainer(ContainerRec *c, int nd, int sd, short id, unsigned char *name)
{
        int	i;
        
	(*c).nd=nd;				// number of data objects
	(*c).sd=sd;				// size in bytes of each data object
	(*c).id=id;				// identification number
	for(i=0;i<=name[0];i++)
            (*c).name[i]=name[i];		// name
	(*c).xcontainer=0;
}
void msh_addData(MeshPtr m, ContainerRec *c)
{
	ContainerPtr	newxcH,xcH;
	long		temp;
	
	newxcH=(ContainerPtr)calloc(1,sizeof(ContainerRec));
	
	(*c).data=calloc((*c).nd,(*c).sd);
	if((*c).data==NULL)	printf("not enough memory");
	
	if((*m).xcontainer)
	{
		xcH=(ContainerPtr)(*m).xcontainer;
		while((long)xcH)
		{
			if((*xcH).id>=(*c).id)	(*c).id++;

			temp=(*xcH).xcontainer;
			if(temp)	xcH=(ContainerPtr)temp;
			else		break;
		}
		(*xcH).xcontainer=(long)newxcH;
	}
	else
		(*m).xcontainer=(long)newxcH;
	
	*newxcH=*c;
}
void msh_deleteData(MeshPtr m, ContainerRec c)
{
    ContainerPtr	xcH,parent;
    long			temp;
    
    if((*m).xcontainer)
    {
        xcH=(ContainerPtr)(*m).xcontainer;
        if((*xcH).id==c.id)
        {
            (*m).xcontainer=(*xcH).xcontainer;
            
            free((*xcH).data);
            free(xcH);
        }
        else	
        while((long)xcH)
        {
            if((*xcH).id==c.id)
            {
                (*parent).xcontainer=(*xcH).xcontainer;
                
                free((*xcH).data);
                free(xcH);
                break;
            }

            parent=xcH;
            temp=(*xcH).xcontainer;
            if(temp)	xcH=(ContainerPtr)temp;
            else		break;
        }
    }
}
bool msh_findData(MeshPtr m, ContainerRec *c)
{
	ContainerPtr	xcH;
	long		temp;
	bool		found=false;
        bool		cmp;
        int		i;
	
	if((*m).xcontainer)
	{
		xcH=(ContainerPtr)(*m).xcontainer;
		while((long)xcH)
		{
			if((*xcH).name[0]==(*c).name[0])
                        {
                            cmp=true;
                            for(i=0;i<=(*c).name[0];i++)
                            if((*xcH).name[i]!=(*c).name[i])
                            {
                                cmp=false;
                                break;
                            }
                            if(	((*c).id>0 && (*xcH).id==(*c).id) ||
                                    ((*c).name[0]>0 && cmp))
                            {
                                    *c=*xcH;
                                    found=true;
                                    break;
                            }
                        }

			temp=(*xcH).xcontainer;
			if(temp)	xcH=(ContainerPtr)temp;
			else		break;
		}
	}
	
	return found;
}
#pragma mark _
float3D* msh_getPointsPtr(MeshPtr m)
{
	return (*m).p;
}
int3D* msh_getTrianglesPtr(MeshPtr m)
{
	return (*m).t;
}
int2D* msh_getEdgesPtr(MeshPtr m)
{
	ContainerRec	c;

	msh_setContainer(&c,(*m).nt*3/2,sizeof(int2D),0,(unsigned char *)"edge");
	if(msh_findData(m,&c)==false){	msh_addData(m,&c); msh_setEdges(m);}
	
	return (int2D*)c.data;
}
NTriRec* msh_getNeighborTrianglesPtr(MeshPtr m)
{
	ContainerRec	c;

	msh_setContainer(&c,(*m).np,sizeof(NTriRec),0,"\pntri");
	if(msh_findData(m,&c)==false)
        {   msh_addData(m,&c);
            msh_setNeighborTriangles(m);}
	
	return (NTriRec*)c.data;
}
NEdgeRec* msh_getNeighborEdgesPtr(MeshPtr m)
{
	ContainerRec	c;

	msh_setContainer(&c,(*m).np,sizeof(NEdgeRec),0,"\pnedg");
	if(msh_findData(m,&c)==false){	msh_addData(m,&c); msh_setNeighborEdges(m);}
	
	return (NEdgeRec*)c.data;
}
float3D* msh_getTexturePtr(MeshPtr m)
{
    ContainerRec	c;

    msh_setContainer(&c,(*m).np,sizeof(float3D),0,"\ptxtr");
    if(msh_findData(m,&c)==false)	msh_addData(m,&c);
    
    return (float3D*)c.data;
}
char* msh_getLabelPtr(MeshPtr m)
{
	ContainerRec	c;

	msh_setContainer(&c,(*m).np,sizeof(char),0,"\plabl");
	if(msh_findData(m,&c)==false)	msh_addData(m,&c);
	
	return (char*)c.data;
}
MeshCurveRec* msh_getMeshCurvePtr(MeshPtr m)
{
	ContainerRec	c;

	msh_setContainer(&c,1,sizeof(MeshCurveRec),0,"\pmeshcurve");
	if(msh_findData(m,&c)==false)	msh_addData(m,&c);
	
	return (MeshCurveRec*)c.data;
}
MeshEdgeCurveRec* msh_getMeshEdgeCurvePtr(MeshPtr m)
{
	ContainerRec	c;

	msh_setContainer(&c,1,sizeof(MeshEdgeCurveRec),0,"\pmeshedgecurve");
	if(msh_findData(m,&c)==false)	msh_addData(m,&c);
	
	return (MeshEdgeCurveRec*)c.data;
}
char* msh_getDataPtr(MeshPtr m,ContainerRec c)
{
    if(msh_findData(m,&c)==false)	msh_addData(m,&c);
    return (char*)c.data;
}
#pragma mark _
void msh_deleteTexturePtr(MeshPtr m)
{
    ContainerRec c;
    
    msh_setContainer(&c,(*m).np,sizeof(float3D),0,"\ptxtr");
    if(msh_findData(m,&c))
        msh_deleteData(m,c);
}
void msh_deleteNeighborTrianglesPtr(MeshPtr m)
{
    ContainerRec c;
    
    msh_setContainer(&c,(*m).np,sizeof(NTriRec),0,"\pntri");
    if(msh_findData(m,&c))
        msh_deleteData(m,c);
}
#pragma mark -
void update(MeshPtr mesh)
{
    int		i;
    float3D		*p;
    
    p = (*mesh).p;
    
    (*mesh).bbox.siz=p[0];
    (*mesh).bbox.ide=p[0];
    for(i=0;i<(*mesh).np;i++)
    {
            if(p[i].x<(*mesh).bbox.siz.x)
                    (*mesh).bbox.siz.x=p[i].x;
            if(p[i].y<(*mesh).bbox.siz.y)
                    (*mesh).bbox.siz.y=p[i].y;
            if(p[i].z<(*mesh).bbox.siz.z)
                    (*mesh).bbox.siz.z=p[i].z;
            
            if(p[i].x>(*mesh).bbox.ide.x)
                    (*mesh).bbox.ide.x=p[i].x;
            if(p[i].y>(*mesh).bbox.ide.y)
                    (*mesh).bbox.ide.y=p[i].y;
            if(p[i].z>(*mesh).bbox.ide.z)
                    (*mesh).bbox.ide.z=p[i].z;
    }
    (*mesh).center = sca3D(add3D((*mesh).bbox.ide,(*mesh).bbox.siz),0.5);
}
void msh_setNeighborTriangles(MeshPtr m)
{
	int		i;
	int3D		*T;
	NTriRec		*NT;
	
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	
	for(i=0;i<(*m).np;i++)
		NT[i].n = 0;
		
	for(i=0;i<(*m).nt;i++)
	{
		NT[T[i].a].t[NT[T[i].a].n++] = i;
		NT[T[i].b].t[NT[T[i].b].n++] = i;
		NT[T[i].c].t[NT[T[i].c].n++] = i;
	}
	
	//TEST
	for(i=0;i<(*m).np;i++)
		if(NT[i].n>=SIZESTACK)
			printf("more neighbors than what we can handle");
}
void msh_setNeighborEdges(MeshPtr m)
{
	int			i;
	int2D		*E;
	NEdgeRec	*NE;
	
	E = msh_getEdgesPtr(m);
	NE = msh_getNeighborEdgesPtr(m);
	
	for(i=0;i<(*m).np;i++)
		NE[i].n = 0;
		
	for(i=0;i<(*m).nt*3/2;i++)
	{
		NE[E[i].a].e[NE[E[i].a].n++] = i;
		NE[E[i].b].e[NE[E[i].b].n++] = i;
	}
}
void msh_setEdges(MeshPtr m)
{
	int		i,n;
	int		a,b,c;
	int3D	*T;
	int2D	*E;

	T = msh_getTrianglesPtr(m);
	E = msh_getEdgesPtr(m);
	
	n=0;
	for(i=0;i<(*m).nt;i++)
	{
		a=T[i].a;	b=T[i].b;	c=T[i].c;
		
		if(a<b)	msh_storeEdge(E,&n,(int2D){a,b});
		else	msh_storeEdge(E,&n,(int2D){b,a});
		
		if(b<c)	msh_storeEdge(E,&n,(int2D){b,c});
		else	msh_storeEdge(E,&n,(int2D){c,b});
		
		if(c<a)	msh_storeEdge(E,&n,(int2D){c,a});
		else	msh_storeEdge(E,&n,(int2D){a,c});
		
		if(n>(*m).nt*3/2)
			printf("more edges than what we can handle");
	}
}
void msh_storeEdge(int2D *E, int *n, int2D e)
{
	int	i,max,min;
	int	j;
	bool	ismin=false,ismax=false,isequal=false;
	
	max=(*n);min=0;i=(*n)/2;
	if((*n)>0)
		do
		{
			if(E[i].a<e.a){	min=i;	ismin=true;}
			else
			if(E[i].a>e.a){	max=i;	ismax=true;}
			else
			if(E[i].b<e.b){	min=i;	ismin=true;}
			else
			if(E[i].b>e.b){	max=i;	ismax=true;}
			else
			{				isequal=true;break;}
			
			i = (max+min)/2;
		}
		while(((max-min)>1 || !ismin) && max>0);

	if(isequal==false)
	{
		i=max;
		for(j=*n;j>i;j--)	E[j]=E[j-1];
		E[i]=e;
		(*n)++;
	}
}
void msh_get_triangle_neighbortriangles(MeshPtr m, int t, int *tn, int *ts)
// find the 3 neighbor triangles of the triangle t in the mesh m
// input: mesh m, triangle index t
// output: tn[3] the triangle index of the neighbors for the 1,2,3 #side,
//		   ts[3] the #side of the associated neighbor triangle
{
	int3D			*T;
	NTriRec			*NT;
	int			i;
	int3D			t0,t1;
	
	T =msh_getTrianglesPtr(m);
	NT=msh_getNeighborTrianglesPtr(m);
	
	t0=T[t];
	
	for(i=0;i<NT[t0.a].n;i++)
	{
		t1=T[NT[t0.a].t[i]];
		
		if(t0.a==t1.a && t0.b==t1.c){	tn[0]=NT[t0.a].t[i];	ts[0]=3;	break;}
		if(t0.a==t1.b && t0.b==t1.a){	tn[0]=NT[t0.a].t[i];	ts[0]=1;	break;}
		if(t0.a==t1.c && t0.b==t1.b){	tn[0]=NT[t0.a].t[i];	ts[0]=2;	break;}
	}

	for(i=0;i<NT[t0.c].n;i++)
	{
		t1=T[NT[t0.c].t[i]];
		
		if(t0.b==t1.a && t0.c==t1.c){	tn[1]=NT[t0.c].t[i];	ts[1]=3;	}
		if(t0.b==t1.b && t0.c==t1.a){	tn[1]=NT[t0.c].t[i];	ts[1]=1;	}
		if(t0.b==t1.c && t0.c==t1.b){	tn[1]=NT[t0.c].t[i];	ts[1]=2;	}

		if(t0.c==t1.a && t0.a==t1.c){	tn[2]=NT[t0.c].t[i];	ts[2]=3;	}
		if(t0.c==t1.b && t0.a==t1.a){	tn[2]=NT[t0.c].t[i];	ts[2]=1;	}
		if(t0.c==t1.c && t0.a==t1.b){	tn[2]=NT[t0.c].t[i];	ts[2]=2;	}
	}
	
}
void msh_get_triangle_neighboredges(MeshPtr m, int t, int *e)
// get the 3 edges that compose the triangle t in the mesh m
// input: mesh m, triangle index t
// output: e[3] edge index for the 1,2,3 #side,
{
	int3D			*T;
	int2D			*E;
	NEdgeRec		*NE;
	int				i;
	int3D			t0;
	int2D			e0;
	
	T =msh_getTrianglesPtr(m);
	E =msh_getEdgesPtr(m);
	NE=msh_getNeighborEdgesPtr(m);
	
	t0=T[t];
	
	for(i=0;i<NE[t0.a].n;i++)
	{
		e0=E[NE[t0.a].e[i]];
		
		if((t0.a==e0.a && t0.b==e0.b) || (t0.b==e0.a && t0.a==e0.b))
		{
			e[0]=NE[t0.a].e[i];
			break;
		}
	}
	for(i=0;i<NE[t0.c].n;i++)
	{
		e0=E[NE[t0.c].e[i]];

		if((t0.b==e0.a && t0.c==e0.b) || (t0.c==e0.a && t0.b==e0.b))
			e[1]=NE[t0.c].e[i];
		if((t0.c==e0.a && t0.a==e0.b) || (t0.a==e0.a && t0.c==e0.b))
			e[2]=NE[t0.c].e[i];
	}
}

#pragma mark -
#pragma mark __________Conjugated gradients
#define EPS 1.0e-18
void msh_linbcgDN(int n,int en, double D[], double N[], int2D *E, double b[], double x[], int itol, double tol, int itmax, int *iter, double *err)
// Solves A * x = b for x[1..n], given b[1..n], by the iterative biconjugate gradient method.
// On input x[1..n] should be set to an initial guess of the solution (or all zeros); itol is 1,2,3,
// or 4, specifying which convergence test is applied (see text); itmax is the maximum number
// of allowed iterations; and tol is the desired convergence tolerance. On output, x[1..n] is
// reset to the improved solution, iter is the number of iterations actually taken, and err is the
// estimated error. The matrix A is referenced only through the user-supplied routines atimes,
// which computes the product of either A or its transpose on a vector; and asolve, which solves
// Ì*x = b or Ì^T*x = b for some preconditioner matrix Ì (possibly the trivial diagonal part of A).
//
// RT:
// In our PDEs, for each entry of the A=D matrix just the neighbors are non-zero
// So, data in the D matrix is organized in the NT matrix.
{
	unsigned	long j;
	double		ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	double		*p,*pp,*r,*rr,*z,*zz;//Double precision is a good idea in this routine.

	p=dvector_new(n);
	pp=dvector_new(n);
	r=dvector_new(n);
	rr=dvector_new(n);
	z=dvector_new(n);
	zz=dvector_new(n);
	
	//Calculate initial residual.
	*iter=0;
	msh_atimesDN(n,en,D,N,E,x,r);	// not transposed
									// Input to atimes is x[1..n], output is r[1..n];
									// the final 0 indicates that the matrix (not its
									// transpose) is to be used.
	for (j=0;j<n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
// atimes(n,r,rr,0);				// Uncomment this line to get the "minimum resid-
									// ual" variant of the algorithm.
	if (itol == 1) {
		bnrm=msh_snrmDN(n,b,itol);
		msh_asolveDN(n,D,r,z);		// not transposed
									// Input to asolve is r[1..n], output is z[1..n];
									// the final 0 indicates that the matrix Ì (not
									// its transpose) is to be used.
	}
	else if (itol == 2) {
		msh_asolveDN(n,D,b,z);		// not transposed
		bnrm=msh_snrmDN(n,z,itol);
		msh_asolveDN(n,D,r,z);		// not transposed
	}
	else if (itol == 3 || itol == 4) {
		msh_asolveDN(n,D,b,z);		// not transposed
		bnrm=msh_snrmDN(n,z,itol);
		msh_asolveDN(n,D,r,z);		// not transposed
		znrm=msh_snrmDN(n,z,itol);
	} else printf("illegal itol in linbcg");
	while (*iter <= itmax) {		// Main loop.
		++(*iter);
		msh_asolveDN(n,D,rr,zz);		// transposed
									// Final 1 indicates use of transpose matrix Ì^T.
		for (bknum=0.0,j=0;j<n;j++) bknum += z[j]*rr[j];
		// Calculate coefficient bk and direction vectors p and pp.
		if (*iter == 1) {
			for (j=0;j<n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else {
			bk=bknum/bkden;
			for (j=0;j<n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;				// Calculate coefficient ak, new iterate x, and new
									// residuals r and rr.
		msh_atimesDN(n,en,D,N,E,p,z);		// not transposed
		for (akden=0.0,j=0;j<n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		msh_atimesDN(n,en,D,N,E,pp,zz);	// transposed
		for (j=0;j<n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		msh_asolveDN(n,D,r,z);		// not transposed
									// Solve Ì*z = r and check stopping criterion.
		if (itol == 1)
			*err=msh_snrmDN(n,r,itol)/bnrm;
		else if (itol == 2)
			*err=msh_snrmDN(n,z,itol)/bnrm;
		else if (itol == 3 || itol == 4) {
			zm1nrm=znrm;
			znrm=msh_snrmDN(n,z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*msh_snrmDN(n,p,itol);
				*err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				*err=znrm/bnrm;		// Error may not be accurate, so loop again.
				continue;
			}
			xnrm=msh_snrmDN(n,x,itol);
			if (*err <= 0.5*xnrm) *err /= xnrm;
			else {
				*err=znrm/bnrm;		// Error may not be accurate, so loop again.
				continue;
			}
		}
		if (*err <= tol) break;
	}
	
	/*psadd(s,"\pCG: iter=");
	NumToString(*iter,s1);psadd(s,s1);
	psadd(s,"\p err=");
	fnum_to_string(*err,20,s1);
	psadd(s,s1);
	printstring(s);*/

	dvector_dispose(p);
	dvector_dispose(pp);
	dvector_dispose(r);
	dvector_dispose(rr);
	dvector_dispose(z);
	dvector_dispose(zz);
}

double msh_snrmDN(int n, double sx[], int itol)
// Compute one of two norm3Ds for a vector sx[1..n], as signaled by itol. Used by linbcg.
{
	unsigned long i,isamax;
	double ans;
	if (itol <= 3) {
		ans = 0.0;
		for (i=0;i<n;i++) ans += sx[i]*sx[i];		// Vector magnitude norm3D.
			return sqrt(ans);
	} else {
		isamax=1;
		for (i=0;i<n;i++) {							// Largest component norm3D.
			if (sx[i] > fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}
void msh_asolveDN(int n, double D[], double b[], double x[])
{
	unsigned long i;
	for(i=0;i<n;i++) x[i]=(D[i] != 0.0 ? b[i]/D[i] : b[i]);
	// The matrix Ì is the diagonal part of A, stored in the first n elements of sa. Since the
	// transpose matrix has the same diagonal, the ag itrnsp is not used.
}
void msh_atimesDN(int n, int en, double D[], double N[], int2D E[], double x[], double r[])
{
	int	i;
	
	for(i=0;i<n;i++)
		r[i]=D[i]*x[i];
	
	for(i=0;i<en;i++)
	{
		r[E[i].a]+= N[i]*x[E[i].b];
		r[E[i].b]+= N[i]*x[E[i].a];
	}
}
#pragma mark __________SOR
void msh_SORDN(MeshPtr m, double D[], double N[], double B[], double X[], double tol, int itmax, int *iter, double *err)
{
	int				i,j,indx;
	double			resid,omega=1;
	int2D			*E;
	NEdgeRec		*NE;
	
	E = msh_getEdgesPtr(m);
	NE = msh_getNeighborEdgesPtr(m);

	for((*iter)=0;(*iter)<itmax;(*iter)++)
	{
		(*err)=0.0;
		for(i=0;i<(*m).np;i++)
		{
			resid=B[i]-D[i]*X[i];
			for(j=0;j<NE[i].n;j++)
			{
				if(E[NE[i].e[j]].a!=i)	indx=E[NE[i].e[j]].a;
				else					indx=E[NE[i].e[j]].b;
				
				resid-=N[NE[i].e[j]]*X[indx];
			}
			X[i]=X[i]+omega*resid/D[i];

			(*err)+=fabs(resid);
		}
		if(*err<=tol)	break;
	}
}
#pragma mark ____________ draw matrix
/*void msh_drawDN(int n, int en, double D[], double N[], int2D E[],int dec)
{
	int			i;
	Rect		w,r;
	int			rw,cw,sz=25,fz=11,e=3;
	char		s[256];
	RGBColor	gray={0xdfff,0xdfff,0xdfff};
	
	GetWindowPortBounds(FrontWindow(),&w);
        EraseRect(&w);

	TextSize(fz);

	rw=(w.bottom/(fz+1)>n)?n:(w.bottom/(fz+1));
	cw=(w.right/sz>n)?n:(w.right/sz);
	
	RGBForeColor(&gray);
	// background
	for(i=0;i<((rw<cw)?rw:cw);i+=10)
	{
		SetRect(&r,0,fz+(i-1)*(fz+e),w.right,fz+(i+4)*(fz+e));
		PaintRect(&r);
	}
	ForeColor(blackColor);
	//diagonal
	for(i=0;i<((rw<cw)?rw:cw);i++)
	{
		MoveTo(i*sz,fz+i*(fz+e));
		
		if(D[i]!=0)	fnum_to_string(D[i],dec,s);
		else		pscopy(s,"\p-");
		
		DrawString(s);
	}
	
	//neighbors
	for(i=0;i<en;i++)
	{
		if(E[i].a<rw && E[i].b<cw)
		{
			MoveTo(E[i].b*sz,fz+E[i].a*(fz+e));
			
			if(N[i]!=0)	fnum_to_string(N[i],dec,s);
			else		pscopy(s,"\p-");
			
			DrawString(s);
		}
		if(E[i].b<rw && E[i].a<cw)
		{
			MoveTo(E[i].a*sz,fz+E[i].b*(fz+e));
			
			if(N[i]!=0)	fnum_to_string(N[i],dec,s);
			else		pscopy(s,"\p-");
			
			DrawString(s);
		}
	}
}*/
#pragma mark -
void msh_laplace_fe(MeshPtr m, double *U,int nDir, int *iDir, double *Dir,
											int nNeu, int *iNeu, double *Neu,
											int nMFC, int2D *eMFC, double *gMFC,double *hMFC)
//	Solves the Laplace equation over a mesh.
//	Accepts Dirichlet, Neumann and inhomoegneous Multifreedom boundary conditions
//	input:	m: the mesh, U:vector for the resulting value in each point,
//			nDir: number of Dirichlet conditions, iDir[]: vector with the mesh point index for the Dirichlet conditions,
//			Dir[]: vector with the value of the Dirichlet condition,
//			nNeu: number of Neumann conditions, iNeu[]: vector with the mesh point index for the Neumann conditions,
//			Neu[]: vector with the value of the Neumann condition,
//			nMFC: number of Multifreedom conditions, eMFC[]: vector with the mesh edge index for the Multifreedom conditions,
//			mMFC[]: vector with the multiplicative value of the Multifreedom condition,
//			nMFC[]: vector with the additive value of the Multifreedom condition,
{
	float3D		*p;
	int2D		*E,*tempE;
	NEdgeRec	*NE;
	int			i,j,k,indx;
	float3D		a,b,c;
	float3D		ir,jr,is,js;
	int			pstack[SIZESTACK],estack[SIZESTACK];
	double		coef;
	double		*D,*N,*B;
	char		*M;
	int			iter=0;
	double		err=0;
	bool		isDirichlet;
	double		g,h;
	int			ms,sl;

	if(m==NULL)	return;
	
	p = msh_getPointsPtr(m);		// points of the mesh
	E = msh_getEdgesPtr(m);			// Edges of the mesh
	NE = msh_getNeighborEdgesPtr(m);// Neighbor Edges for each point

	D=dvector_new((*m).np);		// Diagonal
	N=dvector_new((*m).nt*3/2);	// Neighbor (non diagonal)
	B=dvector_new((*m).np);		// Boundary
	
	if(nMFC>0)
		M=cvector_new((*m).np);	// Mark to avoid double masters or slaves in MFC

	// 1. prepare the matrix
	// ---------------------
	//	1.1. border conditions
	//		1.1.a Dirichlet
			for(i=0;i<nDir;i++)
			{
				B[iDir[i]]=Dir[i];
				if(nMFC>0)
					M[iDir[i]]=1;
			}
	//		1.1.b Neumann
			for(i=0;i<nNeu;i++)
				B[iNeu[i]]+=Neu[i];

	//	1.2. coefficient matrix
			#define cot(x,y) dot3D(x,y)/norm3D(cross3D(x,y))
			for(i=0;i<(*m).np;i++)
			{
				isDirichlet=false;
				for(k=0;k<nDir;k++)
                                    if(i==iDir[k]){	isDirichlet=true;	break;	}
				
				if(isDirichlet==false)
				{
					msh_esort(m,i,pstack,estack);
					D[i]=0;
					for(j=0;j<NE[i].n;j++)
					{
						a=p[pstack[(j+2)%NE[i].n]];
						b=p[pstack[(j+1)%NE[i].n]];
						c=p[pstack[j]];
						
						if(E[estack[(j+1)%NE[i].n]].a!=i)	indx=E[estack[(j+1)%NE[i].n]].a;
						else					indx=E[estack[(j+1)%NE[i].n]].b;
						
						ir=sub3D(p[i],a);	ir=sca3D(ir,1/norm3D(ir));
						jr=sub3D(b,a);	jr=sca3D(jr,1/norm3D(jr));
						
						is=sub3D(b,c);	is=sca3D(is,1/norm3D(is));
						js=sub3D(p[i],c);	js=sca3D(js,1/norm3D(js));

						coef = 0.5*(cot(ir,jr)+cot(is,js));
						
						D[i]+=coef;

						isDirichlet=false;
						for(k=0;k<nDir;k++)
							if(indx==iDir[k]){	isDirichlet=true;	break;	}
						if(isDirichlet==false)
							N[estack[(j+1)%NE[i].n]]=-coef;
						else
							B[i]+=B[indx]*coef;
					}
				}
				else
					D[i]=1;
			}
			#undef cot

	//	1.3. multifreedom constraints of the form Usl=g*Ums+h
			if(nMFC>0)
			{
				tempE = i2vector_new((*m).nt*3/2);
				for(i=0;i<(*m).nt*3/2;i++)
					tempE[i]=E[i];

				for(i=0;i<nMFC;i++)
				if(M[eMFC[i].a]+M[eMFC[i].b]==0)	// (avoid modifying twice)
				{
					g=gMFC[i];
					h=hMFC[i];
					sl=eMFC[i].a;
					ms=eMFC[i].b;
					// 1st. Master Diagonal and Boundary absorves slave's.
					//      Slave Diagonal is set to 1
					D[ms]+=g*g*D[sl];
					B[ms]+=g*h*B[sl];
					D[sl]=1;
					// 2nd. Slave Neighbors adjoined to master's. Inhomogeneous
					//      conditions are substracted
					for(j=0;j<NE[sl].n;j++)
					{
						if(E[NE[sl].e[j]].a==sl)
						{	E[NE[sl].e[j]].a=ms;
							B[E[NE[sl].e[j]].b]-=h*N[NE[sl].e[j]];
						}
						else
						{	E[NE[sl].e[j]].b=ms;
							B[E[NE[sl].e[j]].a]-=h*N[NE[sl].e[j]];
						}
						
						N[NE[sl].e[j]]*=g;
					}
					// 3rd. Master and slave are marked out
					M[ms]=1;
					M[sl]=1;
				}
			}

	// 2. solve
	// ---------
	//	2.1 conjugate gradients in the mesh sparse matrix
		msh_linbcgDN((*m).np,(*m).nt*3/2, D,N,E,B, U, 2, 1e-10, 5*(*m).np, &iter, &err);
	//	2.2 if there are multifreedom constraints, solve for the slaves
		if(nMFC>0)
		{
			// solve for the Slaves
			for(i=0;i<nMFC;i++)
			{
				g=gMFC[i];
				h=hMFC[i];
				sl=eMFC[i].a;
				ms=eMFC[i].b;
			
				U[sl]=g*U[ms]+h;
			}
			
			// restore the Edges vector and dispose the temporal
			for(i=0;i<(*m).nt*3/2;i++)
				E[i]=tempE[i];
			i2vector_dispose(tempE);
		}	
			
	// 3. Dispose
	// ----------
		dvector_dispose(D);
		dvector_dispose(N);
		dvector_dispose(B);
		
		if(nMFC>0)
			cvector_dispose(M);
}

#pragma mark -
//
// point and triangle functions
//
int msh_psort(MeshPtr m, int indx, int *pstack)
{
	int			i,j,n;
	int			nstack;
	int3D		*T;
	float3D		*p;
	NTriRec		*NT;
	int			stack[SIZESTACK];
	int			val;
	bool		isNest;
	int			old_n;
	
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	
	//1. stack the first two vertices of the first triangle
	//2. stack further points to the right
	
	nstack = NT[indx].n;
	for(i=0;i<nstack;i++)
		stack[i]=NT[indx].t[i];

	//1.
	n=0;
	
	if(T[stack[0]].a==indx)
	{
		pstack[n++]=T[stack[0]].b;
		pstack[n++]=T[stack[0]].c;
	}
	else if(T[stack[0]].b==indx)
	{
		pstack[n++]=T[stack[0]].c;
		pstack[n++]=T[stack[0]].a;
	}
	else if(T[stack[0]].c==indx)
	{
		pstack[n++]=T[stack[0]].a;
		pstack[n++]=T[stack[0]].b;
	}

	for(j=0;j<nstack-1;j++)
		stack[j]=stack[j+1];
	
	//2.
	do
	{
		isNest = false;
		old_n = n;
		for(i=0;i<nstack;i++)
		{
			if(T[stack[i]].a==pstack[n-1])
			{
				val=T[stack[i]].b;

				for(j=0;j<n;j++)
					if(val==pstack[j])	isNest=true;
				
				if(!isNest)	pstack[n++]=val;
				break;
			}
			else if(T[stack[i]].b==pstack[n-1])
			{
				val=T[stack[i]].c;

				for(j=0;j<n;j++)
					if(val==pstack[j])	isNest=true;
				
				if(!isNest)	pstack[n++]=val;
				break;
			}
			else if(T[stack[i]].c==pstack[n-1])
			{
				val=T[stack[i]].a;

				for(j=0;j<n;j++)
					if(val==pstack[j])	isNest=true;
				
				if(!isNest)	pstack[n++]=val;
				break;
			}
		}

		for(j=i;j<nstack-1;j++)
			stack[j]=stack[j+1];
		nstack--;
			
	}while(n>old_n && nstack>2);		// n==old_n => nest, nstack==1 => first_p==last_p
	
	return(n);
}
int msh_esort(MeshPtr m, int indx, int *pstack, int *estack)
{
	int			i,j,n;
	int			nt,ne;
	float3D		*p;
	int3D		*T;
	NTriRec		*NT;
	int2D		*E;
	NEdgeRec	*NE;
	int			tstack[SIZESTACK];
	int			val;
	bool		isNest;
	int			old_n;
	
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	E = msh_getEdgesPtr(m);
	NE = msh_getNeighborEdgesPtr(m);
	
	nt = NT[indx].n;
	for(i=0;i<nt;i++)
		tstack[i]=NT[indx].t[i];

	//1. stack the first two vertices of the first triangle
	//	 find and store the first egde
	n=0;
	ne=0;
	
	if(T[tstack[0]].a==indx){	pstack[n++]=T[tstack[0]].b;
								pstack[n++]=T[tstack[0]].c;	}
	else
	if(T[tstack[0]].b==indx){	pstack[n++]=T[tstack[0]].c;
								pstack[n++]=T[tstack[0]].a;	}
	else
	if(T[tstack[0]].c==indx){	pstack[n++]=T[tstack[0]].a;
								pstack[n++]=T[tstack[0]].b;	}

	for(j=0;j<nt-1;j++)
		tstack[j]=tstack[j+1];
	
	for(j=0;j<NE[indx].n;j++)
	if(E[NE[indx].e[j]].a==pstack[0] || E[NE[indx].e[j]].b==pstack[0])
	{	estack[ne++]=NE[indx].e[j];	break;}
	for(j=0;j<NE[indx].n;j++)
	if(E[NE[indx].e[j]].a==pstack[1] || E[NE[indx].e[j]].b==pstack[1])
	{	estack[ne++]=NE[indx].e[j];	break;}
	
	//2. stack further points to the right
	do
	{
		isNest = false;
		old_n = n;
		for(i=0;i<nt;i++)
		{
			if(T[tstack[i]].a==pstack[n-1])
			{	val=T[tstack[i]].b;
				for(j=0;j<n;j++)	if(val==pstack[j])	isNest=true;
				if(!isNest)	pstack[n++]=val;
				break;
			}
			else
			if(T[tstack[i]].b==pstack[n-1])
			{	val=T[tstack[i]].c;
				for(j=0;j<n;j++)	if(val==pstack[j])	isNest=true;
				if(!isNest)	pstack[n++]=val;
				break;
			}
			else
			if(T[tstack[i]].c==pstack[n-1])
			{	val=T[tstack[i]].a;
				for(j=0;j<n;j++)	if(val==pstack[j])	isNest=true;
				if(!isNest)	pstack[n++]=val;
				break;
			}
		}
		for(j=i;j<nt-1;j++)
			tstack[j]=tstack[j+1];
		nt--;
		
		for(j=0;j<NE[indx].n;j++)
		if(E[NE[indx].e[j]].a==pstack[n-1] || E[NE[indx].e[j]].b==pstack[n-1])
		{	estack[ne++]=NE[indx].e[j];	break;}

	}while(n>old_n && nt>2);		// n==old_n => nest, nstack==1 => first_p==last_p
	
	return(n);
}
#pragma mark _
float3D tTriPlane(int nt,MeshPtr m)
{
	float3D	*p,xp;
	int3D	*T;
	float3D	cero={0,0,0};
	
	p = (*m).p;
	T = (*m).t;
	
	xp= cross3D( sub3D(p[T[nt].a],p[T[nt].b]), sub3D(p[T[nt].a],p[T[nt].c]) );
	
	if(norm3D(xp))
		return(sca3D(xp,1.0/norm3D(xp)));
	return(cero);		
}

float3D tTriCenter(int nt, MeshPtr m)
{
	float3D	*p,xp;
	int3D	*t;
	
	p = (*m).p;
	t = (*m).t;
	
	xp.x=(p[t[nt].a].x+p[t[nt].b].x+p[t[nt].c].x)/3.0;
	xp.y=(p[t[nt].a].y+p[t[nt].b].y+p[t[nt].c].y)/3.0;
	xp.z=(p[t[nt].a].z+p[t[nt].b].z+p[t[nt].c].z)/3.0;
	
	return(xp);
}

float tTriArea(MeshPtr m, int nt)
{
	float3D	*p;
	int3D	*t;
	float	aa;
	
	p = (*m).p;
	t = (*m).t;
	
	aa=norm3D( cross3D(sub3D(p[t[nt].a],p[t[nt].c]),sub3D(p[t[nt].b],p[t[nt].a])) )/2.0;

	return(aa);
}
float3D tTrinorm3Dal(MeshPtr m, int nt)
{
	float3D n;
	float3D	*p;
	int3D	*T;
	
	p = (*m).p;
	T = (*m).t;
	
	n=cross3D(sub3D(p[T[nt].b],p[T[nt].a]),sub3D(p[T[nt].c],p[T[nt].a]));
	n=sca3D(n,1/norm3D(n));
	
	return n;
}
void plane(int indx, MeshPtr m, float3D *np, float3D *pp)
{
	int			i;
	float		area=0;
	float3D		zero={0,0,0};
	NTriRec		*NT;
	int3D		*T;
	
	NT = msh_getNeighborTrianglesPtr(m);
	T = msh_getTrianglesPtr(m);
	
	*np=zero;
	*pp=zero;
	
	for(i=0;i<NT[indx].n;i++)
	{
		*np=add3D(*np,tTriPlane(NT[indx].t[i],m));
		*pp=add3D(*pp,
					sca3D(	tTriCenter(NT[indx].t[i],m),
							tTriArea(m,NT[indx].t[i])	  ));
		area+=tTriArea(m,NT[indx].t[i]);
	}

	if(norm3D(*np))
		*np=sca3D(*np,1.0/norm3D(*np));
	if(area)
		*pp=sca3D(*pp,1.0/area);
}
float3D neighbor3Dmean(MeshPtr m, int indx)
{
	int			i;
	float3D		*p;
	NTriRec		*NT;
	int3D		*T;
	float3D		zero = {0,0,0};
	float3D		sum;
	
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	
	sum = zero;
	for(i=0;i<NT[indx].n;i++)
	{
		if(T[NT[indx].t[i]].a!=indx)	sum = add3D(sum,p[T[NT[indx].t[i]].a]);
		if(T[NT[indx].t[i]].b!=indx)	sum = add3D(sum,p[T[NT[indx].t[i]].b]);
		if(T[NT[indx].t[i]].c!=indx)	sum = add3D(sum,p[T[NT[indx].t[i]].c]);
	}
	
	return sca3D(sum , 1/(float)(2*NT[indx].n));
}
void msh_LocalBaseFromPoint(MeshPtr m, int i, float3D *iv, float3D *jv, float3D *kv)
{
	float3D	*p;
	int3D	*T;
	NTriRec	*NT;
	
	int		j;
	float3D	zero={0,0,0};

	//----------------
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	//----------------
	
	*iv=*jv=*kv=zero;

	//kv = norm3Dal
		for(j=0;j<NT[i].n;j++)
			*kv = add3D(*kv,tTriPlane(NT[i].t[j],m));
		*kv = sca3D(*kv,1/norm3D(*kv));
	//iv = tangent
		if(T[NT[i].t[0]].a!=i)
			*iv = sub3D(p[T[NT[i].t[0]].a],p[i]);
		else
			*iv = sub3D(p[T[NT[i].t[0]].b],p[i]);
		*iv = sca3D(*iv,1/norm3D(*iv));
	//jv = iv x kv
		*jv = cross3D(*kv,*iv);
}
void msh_LocalBaseFromTriangle(MeshPtr m, int t, float3D *iv, float3D *jv, float3D *kv)
{
	float3D	*p;
	int3D	*T;
	NTriRec	*NT;
	float3D	zero={0,0,0};

	//----------------
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	//----------------
	
	*iv=*jv=*kv=zero;

	//kv = norm3Dal
		*kv = add3D(*kv,tTriPlane(t,m));
		*kv = sca3D(*kv,1/norm3D(*kv));
	//iv = tangent
		*iv = sub3D( p[T[t].a] , tTriCenter(t,m) );
		*iv = sca3D(*iv,1/norm3D(*iv));
	//jv = iv x kv
		*jv = cross3D(*kv,*iv);
}
float3D msh_VectorFromAngle(float theta, float3D iv, float3D jv)
{
	return add3D(sca3D(iv,cos(theta)),sca3D(jv,sin(theta)));
}
float3D msh_proj3DtoBase(float3D p, float3D iv, float3D jv, float3D kv)
{
	float3D	pp;
	
	pp.x = dot3D(p,iv);
	pp.y = dot3D(p,jv);
	pp.z = dot3D(p,kv);
	
	return pp;
}
float determinant(float3D a, float3D b, float3D c)
{
	float	D=	a.x*(b.y*c.z-c.y*b.z)+
				a.y*(b.z*c.x-c.z*b.x)+
				a.z*(b.x*c.y-c.x*b.y);
	return D;
}
float volume(float3D *p, int3D *t, int nt) // Code By S Melax from http://www.melax.com/volint/
{
	float	vol=0;
	int		i;
	
	for(i=0;i<nt;i++)
		vol += determinant(p[t[i].a],p[t[i].b],p[t[i].c]); //divide by 6 later for efficiency
	return fabs(vol)/6.0f;  // since the determinant give 6 times tetra volume
}
float surface_area(float3D *p, int3D *t, int nt)
{
	int		i;
	float	area;
	
	area=0;
	for(i=0;i<nt;i++)
		area+=triArea(p[t[i].a],p[t[i].b],p[t[i].c]);
	return area;
}

#pragma mark -
float3D	msh_meshpoint_to_float3D(MeshPtr m, MeshPointRec mp)
// if t<0, the the point is the mesh vertex -t-1,
// on the contrary, t is a triangle index and x,y the triangle coordinates for the point
{
	float3D		*p;
	int3D		*T;
	float3D		x;
	
	p=msh_getPointsPtr(m);
	T=msh_getTrianglesPtr(m);
	
	if(mp.t<0)
		x=p[-mp.t-1];
	else
		x = add3D(add3D(sca3D(sub3D(p[T[mp.t].b],p[T[mp.t].a]),mp.x),
					sca3D(sub3D(p[T[mp.t].c],p[T[mp.t].a]),mp.y)),
				p[T[mp.t].a]);

	return x;
}
MeshPointRec msh_float3D_to_meshpoint(MeshPtr m, float3D x)
// Find the closest representation of the point _x_ (in cartesian coords)
// as a mesh point (a triangle plus natural coordinates over the triangle)
// Input: a cartesian point _x_ with center at {0,0,0}; a spherical mesh _m_
// Output: the meshPoint coordinates, i.e., triangle plus natural coords
{
	float3D			*p,tr[3];
	int3D			*T;
	NTriRec			*NT;
	float3D			m0;
	float2D			cp;
	MeshPointRec	mp={-1,0,0};
	int				i,j,pindex,t;
	
	p=msh_getPointsPtr(m);
	T=msh_getTrianglesPtr(m);
	NT=msh_getNeighborTrianglesPtr(m);
	
	// find the closest point in the spherical representation
	pindex=0;
	for(i=0;i<(*m).np;i++)
	{
		if( norm3D(sub3D(x,sca3D(sub3D(p[  i   ],(*m).center),1/norm3D(sub3D(p[i],(*m).center))))) <
			norm3D(sub3D(x,sca3D(sub3D(p[pindex],(*m).center),1/norm3D(sub3D(p[i],(*m).center))))) )
		{
			pindex = i;
		}
	}
	
	t=-1;
	for(j=0;j<NT[pindex].n;j++)
	{
		m0=sub3D(p[T[NT[pindex].t[j]].a],(*m).center);	tr[0]=sca3D(m0,1/norm3D(m0));
		m0=sub3D(p[T[NT[pindex].t[j]].b],(*m).center);	tr[1]=sca3D(m0,1/norm3D(m0));
		m0=sub3D(p[T[NT[pindex].t[j]].c],(*m).center);	tr[2]=sca3D(m0,1/norm3D(m0));

		if(intersect_VectorTriangle(x,tr,&cp,0.00001))
		{
			t = NT[pindex].t[j];
			break;
		}
	}
	
	if(t>=0)
		mp=(MeshPointRec){t,cp.x,cp.y};
	else
		printf("triangle not found for meshPoint\n");

	return mp;
}
MeshPointRec msh_float3D_to_meshpoint_n(MeshPtr m, float3D x)
// Find the closest representation of the point _x_ (in cartesian coords)
// as a mesh point (a triangle plus natural coordinates over the triangle).
// The triangle selected is the one with norm3Dal vector more parallel to the
// _x_ vector.
// Input: a cartesian point _x_ with center at {0,0,0}; a spherical mesh _m_
// Output: the meshPoint coordinates, i.e., triangle plus natural coords
{
	float3D			*p,tr[3];
	int3D			*T;
	NTriRec			*NT;
	float3D			m0;
	float2D			cp;
	MeshPointRec	mp={-1,0,0};
	int				i,pindex,t;
	
	p=msh_getPointsPtr(m);
	T=msh_getTrianglesPtr(m);
	NT=msh_getNeighborTrianglesPtr(m);
	
	// find the closest point in the spherical representation
	pindex=0;
	for(i=0;i<(*m).np;i++)
	{
		if( norm3D(sub3D(x,sca3D(sub3D(p[  i   ],(*m).center),1/norm3D(sub3D(p[  i   ],(*m).center))))) <
			norm3D(sub3D(x,sca3D(sub3D(p[pindex],(*m).center),1/norm3D(sub3D(p[pindex],(*m).center))))) )
		{
			pindex = i;
		}
	}
	
	t=NT[pindex].t[0];
	/*for(j=0;j<NT[pindex].n;j++)
	{
		if(	fabs(dot3D(x,tTrinorm3Dal(m,NT[pindex].t[j]))) <
			fabs(dot3D(x,tTrinorm3Dal(m,t              ))))
		{
			t = NT[pindex].t[j];
		}
	}*/
	
	m0=sub3D(p[T[t].a],(*m).center);	tr[0]=sca3D(m0,1/norm3D(m0));
	m0=sub3D(p[T[t].b],(*m).center);	tr[1]=sca3D(m0,1/norm3D(m0));
	m0=sub3D(p[T[t].c],(*m).center);	tr[2]=sca3D(m0,1/norm3D(m0));
	intersect_VectorTriangle(x,tr,&cp,0.00001);
	mp=(MeshPointRec){t,cp.x,cp.y};

	return mp;
}
MeshPointRec msh_float3D_to_meshpoint_proj(MeshPtr m, float3D x)
// Find the closest representation of the point _x_ (in cartesian coords)
// as a mesh point (a triangle plus natural coordinates over the triangle).
// First the point is projected over the plane of the triangle, then
// the inclusion is verifyied, and finally we find its coordinates.
// Input: a cartesian point _x_ with center at {0,0,0}; a spherical mesh _m_
// Output: the meshPoint coordinates, i.e., triangle plus natural coords
{
	float3D			*p,tr[3],xx;
	float3D			iv,jv,kv;
	float3D			pa,pb,pc;
	float2D			c;
	float			A,B,C,d,e;
	int3D			*T;
	NTriRec			*NT;
	float3D			m0;
	float2D			cp;
	MeshPointRec	mp={-1,0,0};
	int				i,j,pindex,t,tt;
	
	p=msh_getPointsPtr(m);
	T=msh_getTrianglesPtr(m);
	NT=msh_getNeighborTrianglesPtr(m);
	
	// find the closest point in the spherical representation
	pindex=0;
	for(i=0;i<(*m).np;i++)
	{
		if( norm3D(sub3D(x,sca3D(sub3D(p[  i   ],(*m).center),1/norm3D(sub3D(p[  i   ],(*m).center))))) <
			norm3D(sub3D(x,sca3D(sub3D(p[pindex],(*m).center),1/norm3D(sub3D(p[pindex],(*m).center))))) )
		{
			pindex = i;
		}
	}
	
	tt=-1;
	for(j=0;j<NT[pindex].n;j++)
	{
		t=NT[pindex].t[j];
		
		//ForeColor(greenColor);
		//MoveTo3D(p[T[t].a]);LineTo3D(p[T[t].b]);LineTo3D(p[T[t].c]);LineTo3D(p[T[t].a]);
		//ForeColor(redColor);
		//FrameRect3D(p[pindex],1);
		
		pa=sub3D(p[T[t].a],(*m).center);
		pb=sub3D(p[T[t].b],(*m).center);
		pc=sub3D(p[T[t].c],(*m).center);
		
		iv=sub3D(pb,pa);iv=sca3D(iv,1/norm3D(iv));
		kv=cross3D(iv,sub3D(pc,pa));kv=sca3D(kv,1/norm3D(kv));
		jv=cross3D(kv,iv);
		
		d=dot3D(pa,kv);
		e=dot3D(x,kv);
		if(e>0)
		{
			xx=sca3D(x,d/e);
			c=(float2D){dot3D(sub3D(xx,pa),iv),dot3D(sub3D(xx,pa),jv)};
			
			A=dot3D(cross3D(sub3D(pb,pa),sub3D(xx,pa)),kv);
			B=dot3D(cross3D(sub3D(pc,pb),sub3D(xx,pb)),kv);
			C=dot3D(cross3D(sub3D(pa,pc),sub3D(xx,pc)),kv);
			
			if(!(A<0||B<0||C<0))
			{
				tt=t;
				break;
			}
		}
	}
	
	if(tt<0)
		for(t=0;t<(*m).nt;t++)
		{
			//ForeColor(greenColor);
			//MoveTo3D(p[T[t].a]);LineTo3D(p[T[t].b]);LineTo3D(p[T[t].c]);LineTo3D(p[T[t].a]);
			//ForeColor(redColor);
			//FrameRect3D(p[pindex],1);
			
			pa=sub3D(p[T[t].a],(*m).center);
			pb=sub3D(p[T[t].b],(*m).center);
			pc=sub3D(p[T[t].c],(*m).center);
			
			iv=sub3D(pb,pa);iv=sca3D(iv,1/norm3D(iv));
			kv=cross3D(iv,sub3D(pc,pa));kv=sca3D(kv,1/norm3D(kv));
			jv=cross3D(kv,iv);
			
			d=dot3D(pa,kv);
			e=dot3D(x,kv);
			if(e>0)
			{
				xx=sca3D(x,d/e);
				c=(float2D){dot3D(sub3D(xx,pa),iv),dot3D(sub3D(xx,pa),jv)};
				
				A=dot3D(cross3D(sub3D(pb,pa),sub3D(xx,pa)),kv);
				B=dot3D(cross3D(sub3D(pc,pb),sub3D(xx,pb)),kv);
				C=dot3D(cross3D(sub3D(pa,pc),sub3D(xx,pc)),kv);
				
				if(!(A<0||B<0||C<0))
				{
					tt=t;
					break;
				}
			}
		}
	if(tt<0)
		printf("FUNAWER!!!!!\n");
	
	t=tt;
	m0=sub3D(p[T[t].a],(*m).center);	tr[0]=sca3D(m0,1/norm3D(m0));
	m0=sub3D(p[T[t].b],(*m).center);	tr[1]=sca3D(m0,1/norm3D(m0));
	m0=sub3D(p[T[t].c],(*m).center);	tr[2]=sca3D(m0,1/norm3D(m0));
	intersect_VectorTriangle(x,tr,&cp,0.00001);
	mp=(MeshPointRec){t,cp.x,cp.y};

	return mp;
}
float3D	msh_meshedgepoint_to_float3D(MeshPtr m, MeshEdgePointRec mep)
{
	float3D		*p;
	int3D		*T;
	float3D		x;
	
	p=msh_getPointsPtr(m);
	T=msh_getTrianglesPtr(m);
	
	if(mep.ta<0)	// if ta<0 the point is the vertex point tb
		x=p[mep.tb];
	else
	{
		switch(mep.sa)
		{
			case 1: x = add3D( sca3D(p[T[mep.ta].a],1-mep.t), sca3D(p[T[mep.ta].b],mep.t) );	break;
			case 2: x = add3D( sca3D(p[T[mep.ta].b],1-mep.t), sca3D(p[T[mep.ta].c],mep.t) );	break;
			case 3: x = add3D( sca3D(p[T[mep.ta].c],1-mep.t), sca3D(p[T[mep.ta].a],mep.t) );	break;
		}

		/*switch(mep.sb)
		{	case 1: x = add3D( sca3D(p[T[mep.tb].a],mep.t), sca3D(p[T[mep.tb].b],1-mep.t) );	break;
			case 2: x = add3D( sca3D(p[T[mep.tb].b],mep.t), sca3D(p[T[mep.tb].c],1-mep.t) );	break;
			case 3: x = add3D( sca3D(p[T[mep.tb].c],mep.t), sca3D(p[T[mep.tb].a],1-mep.t) );	break;	}*/
	}
	return x;
}
MeshEdgePointRec msh_float3D_to_meshedgepoint(MeshPtr m, float3D x)
{
	float3D				*p,tr[3];
	int3D				*T;
	NTriRec				*NT;
	float3D				m0;
	float2D				cp;
	MeshEdgePointRec	mep;
	float				min,xmin;
	int					i,j,pindex,t;
	
	p=msh_getPointsPtr(m);
	T=msh_getTrianglesPtr(m);
	NT=msh_getNeighborTrianglesPtr(m);
	
	// find the closest point in the spherical representation
	min=100000;
	for(i=0;i<(*m).np;i++)
	{
		m0=sub3D(p[i],(*m).center);
		m0=sca3D(m0,1/norm3D(m0));
		xmin=norm3D(sub3D(x,m0));
		if( xmin<min )
		{
			pindex = i;
			min = xmin;
		}
	}
	
	t=-1;
	for(j=0;j<NT[pindex].n;j++)
	{
		m0=sub3D(p[T[NT[pindex].t[j]].a],(*m).center);	tr[0]=sca3D(m0,1/norm3D(m0));
		m0=sub3D(p[T[NT[pindex].t[j]].b],(*m).center);	tr[1]=sca3D(m0,1/norm3D(m0));
		m0=sub3D(p[T[NT[pindex].t[j]].c],(*m).center);	tr[2]=sca3D(m0,1/norm3D(m0));

		if(intersect_VectorTriangle(x,tr,&cp,0.1))
		{
			t = NT[pindex].t[j];
			break;
		}
	}
	
	//--------------------------------------
	// WARNING
	//_____________________________________
	
	// this function is incomplete.
	// here I have the triangle for the point and its coordinates {t,cp.x,cp.y}.
	// Now I need to get i. the closer point over one of the triangle sides
	// and the neighbor triangle that shares that side
	// or ii. a vertex if it is closer.
	
	mep=(MeshEdgePointRec){t,t,1,1,cp.x};	// not valid

	return mep;
}
