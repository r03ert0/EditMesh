/*
 *  util.c
 *  EditMesh
 *
 *  Created by roberto on Thu Sep 18 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */

#include "util.h"

float getnumfromtxt(int *x, char *t, char *stop)
{
    char		ch;
    char		s[256],ms[256],es[256];
    int			i=*x,sig=1;
    double		num,m=0;
    long		n,e=0;
    
    ch = *(t+i++);
    while((ch<'0' || ch>'9') && ch!='+' && ch!='-')
    {
        if(ch=='#')
        {
                do ch = *(t+i++);
                while(ch!='\n' && ch!='\r');
        }
        ch = *(t+i++);
    }
    s[0]=(char)0;
    
    if(ch=='-')
            sig=-1;
    do
    {
        s[++s[0]]=ch;
        ch = *(t+i++);
        
        if(ch=='.'||ch==',')
        {
            ms[0]=(char)0;
            ch = *(t+i++);
            do
            {
                    ms[ms[0]]=ch;
                    ch = *(t+i++);
            }while(ch>='0' && ch<='9');
            
            ms[(int)ms[0]+1]=' ';
            n=atoi(ms+1);
            m = n/pow(10,ms[0]);
        }
        if(ch=='e')
        {
            es[0]=(char)0;
            ch = *(t+i++);
            do
            {
                    es[++es[0]]=ch;
                    ch = *(t+i++);
            }while(ch>='0' && ch<='9');
            
            es[(int)es[0]+1]=' ';
            e=atoi(es+1);
        }
    }while(ch>='0' && ch<='9');

    s[(int)s[0]+1]=' ';
    n=atoi(s+1);
    num = sig * (n*sig + m) * pow(10,e);
    *x = i;
    *stop=ch;
    
    return(num);
}
int getcharfromtxt(int *x, char *t, char *stop)
{
    char		ch;
    int			i=*x,type=0;
    
    do
	{
		ch = *(t+i);
		if((ch>='0' && ch<='9') || ch=='+' || ch=='-')
			type=1;
		else
		if((ch>='a' && ch<='z') || (ch>='A' && ch<='Z') || ch=='_')
			type=2;
		else
		if(ch=='=')
			type=3;
		else
		if(ch=='#')
		{
			do ch = *(t+i++);
			while(ch!='\n' && ch!='\r');
		}
		else
			i++;
	}
	while(type==0);
	
	*stop=ch;
	*x=i;
	return type;
}
void getwordfromtxt(int *x, char *t, char *word)
{
    char		ch;
    int			i=*x;
    
    word[0]=(char)0;
	do  ch = *(t+i++);
	while(ch==' ');
	
	do
	{
		word[++word[0]]=ch;
		ch = *(t+i++);
	}
	while((ch>='0' && ch<='9') || (ch>='a' && ch<='z') || (ch>='A' && ch<='Z') || ch=='_');
	
	word[++word[0]]=(char)0;
	
	*x=i;
}
void flookforword(int *x, char *t,char *word, int length)
{
	int i,equal=1;
	char	a;
	do
	{
		i=0;
		do
		{
			a=t[(*x)++];
			equal=(a==word[i++]);
			if(word[i]==(char)0)
				break;
		}while(equal);
		if(equal)
			break;
	}while((*x)<length);
} 
#pragma mark -
bool	configured=false;
float	cm_gray[3][256];		// grayscale colour map
float	cm_gray_discret[3][256];	// discret grayscale colour map
float	cm_redgreen[3][256];		// <0.5 red, >0.5 green colour map
float	cm_redgreen_discret[3][256];	// discret <0.5 red, >0.5 green colour map
float	cm_rgb[3][256];			// rgb colour map
float	cm_rgb_discret[3][256];		// discret rgb colour map
float	cm_rygcb[3][256];		// monochromatic colour map
float	cm_rygcb_discret[3][256];	// discret monochromatic colour map
float	cm_colour[3][256];		// monochromatic colour map
float	cm_colour_discret[3][256];	// discret monochromatic colour map
float	cm_colour_flat[3][256];		// flat monochromatic colour map
void initColourmaps(void)
{
    int	i;

    // continous colourmaps
    for(i=0;i<128;i++)
    {
        cm_rgb[2][i]=cm_rgb[1][i+128]=1-i/127.0;
        cm_rgb[1][i]=cm_rgb[0][i+128]=i/127.0;
    }
    for(i=0;i<256;i++)
    {
        cm_gray[0][i]=cm_gray[1][i]=cm_gray[2][i]=i/255.0;
        cm_colour[0][i]=i/255.0;	// init to RED
        cm_colour_flat[0][i]=1;		// init to RED
    }
    for(i=0;i<128;i++)
    {
        cm_redgreen[0][i+128]=i/127.0;
        cm_redgreen[1][i]=1-i/127.0;
    }
    for(i=0;i<64;i++)
    {
        cm_rygcb[0][i+192]=cm_rygcb[1][i+64]=cm_rygcb[1][i+128]=cm_rygcb[2][i]=1;
        cm_rygcb[1][i+192]=cm_rygcb[2][i+64]=1-i/63.0;
        cm_rygcb[0][i+128]=cm_rygcb[1][i]=i/63.0;
    }
    
    // discret colourmaps
    for(i=0;i<86;i++)
    {
        cm_gray_discret[0][i+85]=cm_gray_discret[1][i+85]=cm_gray_discret[2][i+85]=0.5;
        cm_gray_discret[0][i+170]=cm_gray_discret[1][i+170]=cm_gray_discret[2][i+170]=1;
        
        cm_redgreen_discret[1][i]=
        cm_redgreen_discret[0][i+170]=1; 
        
        cm_rgb_discret[2][i]=
        cm_rgb_discret[1][i+85]=
        cm_rgb_discret[0][i+170]=1;
        
        cm_colour_discret[0][i+85]=0.5;
        cm_colour_discret[0][i+170]=1;
    }
    for(i=0;i<51;i++)
        cm_rygcb_discret[0][i+153]=cm_rygcb_discret[0][i+204]=
        cm_rygcb_discret[1][i+51]=cm_rygcb_discret[1][i+102]=cm_rygcb_discret[1][i+153]=
        cm_rygcb_discret[2][i]=cm_rygcb_discret[2][i+51]=1;
    cm_rygcb_discret[0][255]=1;
    
    configured=true;
}
void colourFromColourmap(float index, float *c, int cm)
{
    int	i;
    
    if(configured==false) initColourmaps();
    switch(cm)
    {
        case 0: for(i=0;i<3;i++) c[i]=cm_gray[i][(int)(255*index)];		break;
        case 1: for(i=0;i<3;i++) c[i]=cm_gray_discret[i][(int)(255*index)];	break;
        case 2: for(i=0;i<3;i++) c[i]=cm_redgreen[i][(int)(255*index)];		break;
        case 3: for(i=0;i<3;i++) c[i]=cm_redgreen_discret[i][(int)(255*index)];	break;
        case 4: for(i=0;i<3;i++) c[i]=cm_rgb[i][(int)(255*index)];		break;
        case 5: for(i=0;i<3;i++) c[i]=cm_rgb_discret[i][(int)(255*index)];	break;
        case 6: for(i=0;i<3;i++) c[i]=cm_rygcb[i][(int)(255*index)];		break;
        case 7: for(i=0;i<3;i++) c[i]=cm_rygcb_discret[i][(int)(255*index)];	break;
        case 8: for(i=0;i<3;i++) c[i]=cm_colour[i][(int)(255*index)];		break;
        case 9: for(i=0;i<3;i++) c[i]=cm_colour_discret[i][(int)(255*index)];	break;
        case 10:for(i=0;i<3;i++) c[i]=cm_colour_flat[i][(int)(255*index)];	break;
    }
}
void configureMonochromaticColourmap(float *c)
{
    int	i;
    
    for(i=0;i<256;i++)
    {
        cm_colour[0][i]=c[0]*i/255.0;
        cm_colour[1][i]=c[1]*i/255.0;
        cm_colour[2][i]=c[2]*i/255.0;
        
        cm_colour_flat[0][i]=c[0];
        cm_colour_flat[1][i]=c[1];
        cm_colour_flat[2][i]=c[2];
    }
    for(i=0;i<86;i++)
    {
        cm_colour_discret[0][i]=0;
        cm_colour_discret[0][i+85]=0.5*c[0];
        cm_colour_discret[0][i+170]=c[0];
        
        cm_colour_discret[1][i]=0;
        cm_colour_discret[1][i+85]=0.5*c[1];
        cm_colour_discret[1][i+170]=c[1];
        
        cm_colour_discret[2][i]=0;
        cm_colour_discret[2][i+85]=0.5*c[2];
        cm_colour_discret[2][i+170]=c[2];
    }
}
int getint(int *x, char *b, int endian)
{
    char	n[4];
	char	tmp[4]={0,0,0,1};
	int		hostEndian=(*(int*)tmp==1);
    
    if(endian==hostEndian)
    {
        n[3]=b[(*x)++];
        n[2]=b[(*x)++];
        n[1]=b[(*x)++];
        n[0]=b[(*x)++];
    }
    else
    {
        n[0]=b[(*x)++];
        n[1]=b[(*x)++];
        n[2]=b[(*x)++];
        n[3]=b[(*x)++];
    }
    return *(int*)n;
}
float getfloat(int *x, char *b, int endian)
{
    char n[4];
	char	tmp[4]={0,0,0,1};
	int		hostEndian=(*(int*)tmp==1);
    
    if(endian==hostEndian)
    {
        n[3]=b[(*x)++];
        n[2]=b[(*x)++];
        n[1]=b[(*x)++];
        n[0]=b[(*x)++];
    }
    else
    {
        n[0]=b[(*x)++];
        n[1]=b[(*x)++];
        n[2]=b[(*x)++];
        n[3]=b[(*x)++];
    }
    return *(float*)n;
}
int fgetint(FILE *f, int endian)
{
    char	n[4],b[4];
	int		i,result;
	char	tmp[4]={0,0,0,1};
	int		hostEndian=(*(int*)tmp==1);
	
	for(i=0;i<4;i++) b[i]=fgetc(f);
    
    if(endian==hostEndian)
    {
        for(i=0;i<4;i++)
            n[3-i]=b[i];
        result=*(int*)n;
    }
    else
		result=*(int*)b;
    return result;
}
float fgetfloat(FILE *f, int endian)
{
    char	n[4],b[4];
	int		i;
	float	result;
	char	tmp[4]={0,0,0,1};
	int		hostEndian=(*(int*)tmp==1);
	
	for(i=0;i<4;i++) b[i]=fgetc(f);
    
    if(endian==hostEndian)
    {
        for(i=0;i<4;i++)
            n[3-i]=b[i];
        result=*(float*)n;
    }
    else
		result=*(float*)b;
    return result;
}