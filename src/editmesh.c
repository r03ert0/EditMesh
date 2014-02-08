/*
 *  editmesh.c
 *  EditMesh
 *
 *  Created by roberto on Sat Sep 13 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */

#include "editmesh.h"

// transform globals
float		g_em_scale=2.0;
float		g_em_simplify=0.2;
float		g_em_smooth=0.5;
float		g_em_inflate=0.01;

// fision / fusion globals
float		g_em_fisionthrs=4.0;
int		g_em_maxpoint=75000;
int		g_em_maxtrian=150000;

float		g_em_lengthfusionthrs=0.5;
float		g_em_areafusionthrs=0.5;
float		g_em_fusionP=0.0;

// standard rotations
float		g_em_viewZ[9] =		{ 1,0, 0,	0,1, 0,	 0, 0, 1};
float		g_em_view_Z[9] =	{-1,0, 0,	0,1, 0,	 0, 0,-1};
float		g_em_viewY[9] = 	{ 1,0, 0,	0,0, 1,	 0,-1, 0};
float		g_em_view_Y[9] = 	{ 1,0, 0,	0,0,-1,	 0, 1, 0};
float		g_em_viewX[9] = 	{ 0,0, 1,	0,1, 0,	-1, 0, 0};
float		g_em_view_X[9] = 	{ 0,0,-1,	0,1, 0,	 1, 0, 0};

// selected point for editing
bool		g_em_isSelected=false;
int		g_em_selectedpoint;

// colour maps
float		g_em_fire_r[] = {0,0,1,25,49,73,98,122,146,162,173,184,195,207,217,229,240,252,255,255,255,255,255,255,255,255,255,255,255,255,255,255};
float		g_em_fire_g[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,14,35,57,79,101,117,133,147,161,175,190,205,219,234,248,255,255,255,255};
float		g_em_fire_b[] = {31,61,96,130,165,192,220,227,210,181,151,122,93,64,35,5,0,0,0,0,0,0,0,0,0,0,0,35,98,160,223,255};

float		g_em_ice_r[] = {0,0,0,0,0,0,19,29,50,48,79,112,134,158,186,201,217,229,242,250,250,250,250,251,250,250,250,250,251,251,243,230};
float		g_em_ice_g[] = {156,165,176,184,190,196,193,184,171,162,146,125,107,93,81,87,92,97,95,93,93,90,85,69,64,54,47,35,19,0,4,0};
float		g_em_ice_b[] = {140,147,158,166,170,176,209,220,234,225,236,246,250,251,250,250,245,230,230,222,202,180,163,142,123,114,106,94,84,64,26,27};

//texture: hierarchical curvature global
//int			g_em_smoothiter = 10;//30;

// Talairach transformation Matrix
float3D			g_tal_M[4];
float3DHandle	g_tal_labH;
int				g_tal_nlabels=150;

float g_em_harmonics[512]={	1.00000000000000, 0.96017605125000, 0.92035210250000, 0.88052815375000, 0.84070420500000,
							0.80088025625000, 0.76105630750000, 0.72123235875000, 0.68140841000000, 0.63897514000000,
							0.62566871000000, 0.61236228000000, 0.59905585000000, 0.58574942000000, 0.57244299000000,
							0.55913656000000, 0.54583013000000, 0.53564651666667, 0.52546290333333, 0.51527929000000,
							0.50628711750000, 0.49729494500000, 0.48830277250000, 0.47931060000000, 0.47138139500000,
							0.46345219000000, 0.45792816500000, 0.45240414000000, 0.44688011500000, 0.44135609000000,
							0.43426554000000, 0.42717499000000, 0.42022759000000, 0.41542266000000, 0.41061773000000,
							0.40581280000000, 0.40003754666667, 0.39426229333333, 0.38848704000000, 0.38426593000000,
							0.38222253000000, 0.37795896500000, 0.37369540000000, 0.36945697333333, 0.36521854666667,
							0.36098012000000, 0.35658319000000, 0.35218626000000, 0.34723115000000, 0.34416249000000,
							0.33990001500000, 0.33563754000000, 0.33207006333333, 0.32850258666667, 0.32493511000000,
							0.32192370000000, 0.32146204000000, 0.31658140000000, 0.31170076000000, 0.30911565666667,
							0.30653055333333, 0.30394545000000, 0.30085295500000, 0.29776046000000, 0.29439834000000,
							0.29217690500000, 0.28995547000000, 0.28685711500000, 0.28375876000000, 0.28073308000000,
							0.27770740000000, 0.27468172000000, 0.27407864000000, 0.26926082000000, 0.26799724000000,
							0.26432067000000, 0.26166349500000, 0.25900632000000, 0.25827509000000, 0.25490630000000,
							0.25299427000000, 0.25133610000000, 0.24879822000000, 0.24626034000000, 0.24356134000000,
							0.24126829000000, 0.23876484500000, 0.23626140000000, 0.23521683000000, 0.23280239000000,
							0.22945222000000, 0.22870363000000, 0.22493730000000, 0.22374172500000, 0.22254615000000,
							0.21962409000000, 0.21842812000000, 0.21609493000000, 0.21331905000000, 0.21119159000000,
							0.20906413000000, 0.20712914000000, 0.20519415000000, 0.20333605000000, 0.20297825000000,
							0.20077564000000, 0.19827932000000, 0.19584453000000, 0.19501199000000, 0.19313936000000,
							0.19126673000000, 0.18919508000000, 0.18781886000000, 0.18564902000000, 0.18381105000000,
							0.18305048000000, 0.18121995000000, 0.17867383000000, 0.17791472000000, 0.17537352000000,
							0.17465691000000, 0.17195807000000, 0.17067376000000, 0.16900049000000, 0.16778342000000,
							0.16517833000000, 0.16402040000000, 0.16255434000000, 0.16149811000000, 0.15997708000000,
							0.15716885000000, 0.15617490500000, 0.15518096000000, 0.15354230000000, 0.15145667000000,
							0.14999968000000, 0.14925636000000, 0.14708392000000, 0.14561006000000, 0.14431795000000,
							0.14213172000000, 0.14051592000000, 0.13889498000000, 0.13764985000000, 0.13640159000000,
							0.13446414000000, 0.13425836000000, 0.13162099000000, 0.13037430000000, 0.12912761000000,
							0.12793604000000, 0.12593231000000, 0.12499724000000, 0.12374388000000, 0.12301020000000,
							0.12118306000000, 0.11972722000000, 0.11814676000000, 0.11703368000000, 0.11491345000000,
							0.11397808000000, 0.11221533000000, 0.11134135000000, 0.11043078000000, 0.10835493000000,
							0.10684173000000, 0.10554413000000, 0.10516541000000, 0.10329884000000, 0.10216175000000,
							0.10074545000000, 0.09907018000000, 0.09862638000000, 0.09737607000000, 0.09577937000000,
							0.09430736000000, 0.09296854000000, 0.09255683000000, 0.09039403000000, 0.08934227000000,
							0.08876959000000, 0.08712226000000, 0.08530658000000, 0.08420894000000, 0.08377215000000,
							0.08165161000000, 0.08028018000000, 0.07969181000000, 0.07872528000000, 0.07686730000000,
							0.07597684000000, 0.07494645000000, 0.07289722000000, 0.07188097000000, 0.07035693000000,
							0.06941156000000, 0.06860622000000, 0.06716582000000, 0.06617359000000, 0.06479642000000,
							0.06378989000000, 0.06308818000000, 0.06191645000000, 0.06032592000000, 0.05920655000000,
							0.05760932000000, 0.05678081500000, 0.05595231000000, 0.05417852000000, 0.05361692000000,
							0.05191611000000, 0.05103591000000, 0.04998863000000, 0.04810179000000, 0.04666962000000,
							0.04629311000000, 0.04525017000000, 0.04362446000000, 0.04296560000000, 0.04160564000000,
							0.04027436000000, 0.03863633000000, 0.03823184000000, 0.03715204000000, 0.03549206000000,
							0.03502371000000, 0.03373144000000, 0.03238219000000, 0.03074061000000, 0.02994172000000,
							0.02917591000000, 0.02711291000000, 0.02648273500000, 0.02585256000000, 0.02458071000000,
							0.02330169000000, 0.02200439000000, 0.02121164000000, 0.01963169000000, 0.01805174000000,
							0.01803836000000, 0.01686032000000, 0.01560804000000, 0.01412341000000, 0.01310045000000,
							0.01207749000000, 0.01105453000000, 0.01003157000000, 0.00900861000000, 0.00780362000000,
							0.00650301666667, 0.00520241333333, 0.00390181000000, 0.00260120666667, 0.00130060333333,
							               0,-0.000750301111111, -0.00130060333333, -0.00260120666667, -0.00390181000000,
  							-0.00520241333333, -0.00650301666667, -0.00780362000000, -0.00900861000000, -0.01003157000000,
  							-0.01105453000000, -0.01207749000000, -0.01310045000000, -0.01412341000000, -0.01560804000000,
  							-0.01686032000000, -0.01803836000000, -0.01805174000000, -0.01963169000000, -0.02121164000000,
  							-0.02200439000000, -0.02330169000000, -0.02458071000000, -0.02585256000000, -0.02648273500000,
  							-0.02711291000000, -0.02917591000000, -0.02994172000000, -0.03074061000000, -0.03238219000000,
  							-0.03373144000000, -0.03502371000000, -0.03549206000000, -0.03715204000000, -0.03823184000000,
  							-0.03863633000000, -0.04027436000000, -0.04160564000000, -0.04296560000000, -0.04362446000000,
  							-0.04525017000000, -0.04629311000000, -0.04666962000000, -0.04810179000000, -0.04998863000000,
  							-0.05103591000000, -0.05191611000000, -0.05361692000000, -0.05417852000000, -0.05595231000000,
  							-0.05678081500000, -0.05760932000000, -0.05920655000000, -0.06032592000000, -0.06191645000000,
  							-0.06308818000000, -0.06378989000000, -0.06479642000000, -0.06617359000000, -0.06716582000000,
  							-0.06860622000000, -0.06941156000000, -0.07035693000000, -0.07188097000000, -0.07289722000000,
  							-0.07494645000000, -0.07597684000000, -0.07686730000000, -0.07872528000000, -0.07969181000000,
  							-0.08028018000000, -0.08165161000000, -0.08377215000000, -0.08420894000000, -0.08530658000000,
  							-0.08712226000000, -0.08876959000000, -0.08934227000000, -0.09039403000000, -0.09255683000000,
  							-0.09296854000000, -0.09430736000000, -0.09577937000000, -0.09737607000000, -0.09862638000000,
  							-0.09907018000000, -0.10074545000000, -0.10216175000000, -0.10329884000000, -0.10516541000000,
  							-0.10554413000000, -0.10684173000000, -0.10835493000000, -0.11043078000000, -0.11134135000000,
  							-0.11221533000000, -0.11397808000000, -0.11491345000000, -0.11703368000000, -0.11814676000000,
  							-0.11972722000000, -0.12118306000000, -0.12301020000000, -0.12374388000000, -0.12499724000000,
  							-0.12593231000000, -0.12793604000000, -0.12912761000000, -0.13037430000000, -0.13162099000000,
  							-0.13425836000000, -0.13446414000000, -0.13640159000000, -0.13764985000000, -0.13889498000000,
  							-0.14051592000000, -0.14213172000000, -0.14431795000000, -0.14561006000000, -0.14708392000000,
  							-0.14925636000000, -0.14999968000000, -0.15145667000000, -0.15354230000000, -0.15518096000000,
  							-0.15617490500000, -0.15716885000000, -0.15997708000000, -0.16149811000000, -0.16255434000000,
  							-0.16402040000000, -0.16517833000000, -0.16778342000000, -0.16900049000000, -0.17067376000000,
  							-0.17195807000000, -0.17465691000000, -0.17537352000000, -0.17791472000000, -0.17867383000000,
  							-0.18121995000000, -0.18305048000000, -0.18381105000000, -0.18564902000000, -0.18781886000000,
  							-0.18919508000000, -0.19126673000000, -0.19313936000000, -0.19501199000000, -0.19584453000000,
  							-0.19827932000000, -0.20077564000000, -0.20297825000000, -0.20333605000000, -0.20519415000000,
  							-0.20712914000000, -0.20906413000000, -0.21119159000000, -0.21331905000000, -0.21609493000000,
  							-0.21842812000000, -0.21962409000000, -0.22254615000000, -0.22374172500000, -0.22493730000000,
  							-0.22870363000000, -0.22945222000000, -0.23280239000000, -0.23521683000000, -0.23626140000000,
  							-0.23876484500000, -0.24126829000000, -0.24356134000000, -0.24626034000000, -0.24879822000000,
  							-0.25133610000000, -0.25299427000000, -0.25490630000000, -0.25827509000000, -0.25900632000000,
  							-0.26166349500000, -0.26432067000000, -0.26799724000000, -0.26926082000000, -0.27407864000000,
  							-0.27468172000000, -0.27770740000000, -0.28073308000000, -0.28375876000000, -0.28685711500000,
  							-0.28995547000000, -0.29217690500000, -0.29439834000000, -0.29776046000000, -0.30085295500000,
  							-0.30394545000000, -0.30653055333333, -0.30911565666667, -0.31170076000000, -0.31658140000000,
  							-0.32146204000000, -0.32192370000000, -0.32493511000000, -0.32850258666667, -0.33207006333333,
  							-0.33563754000000, -0.33990001500000, -0.34416249000000, -0.34723115000000, -0.35218626000000,
  							-0.35658319000000, -0.36098012000000, -0.36521854666667, -0.36945697333333, -0.37369540000000,
  							-0.37795896500000, -0.38222253000000, -0.38426593000000, -0.38848704000000, -0.39426229333333,
  							-0.40003754666667, -0.40581280000000, -0.41061773000000, -0.41542266000000, -0.42022759000000,
  							-0.42717499000000, -0.43426554000000, -0.44135609000000, -0.44688011500000, -0.45240414000000,
  							-0.45792816500000, -0.46345219000000, -0.47138139500000, -0.47931060000000, -0.48830277250000,
  							-0.49729494500000, -0.50628711750000, -0.51527929000000, -0.52546290333333, -0.53564651666667,
  							-0.54583013000000, -0.55913656000000, -0.57244299000000, -0.58574942000000, -0.59905585000000,
  							-0.61236228000000, -0.62566871000000, -0.63897514000000, -0.68140841000000, -0.72123235875000,
  							-0.76105630750000, -0.80088025625000, -0.84070420500000, -0.88052815375000, -0.92035210250000,
  							-0.96017605125000, -1.00000000000000};

#pragma mark ____________________get info
/*void em_infoLength(MeshPtr m)
{
	float3D		*p;
	int3D		*T;
	NTriRec		*NT;
	int		i,j;
	int		h[100],maxh;
	float		x,min,max;
	Rect		r,er;
	Str255		s;
	
	int			barsize=5;
	int			ndecimal=2;

	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);	
	
	min=1000;
	max=0;
	for(i=0;i<(*m).np;i++)
		for(j=0;j<NT[i].n;j++)
		{
			if(T[NT[i].t[j]].a!=i)
			{
				x = norm3D(sub3D(p[T[NT[i].t[j]].a],p[i]));
				min = (x<min)?x:min;
				max = (x>max)?x:max;
			}

			if(T[NT[i].t[j]].b!=i)
			{
				x = norm3D(sub3D(p[T[NT[i].t[j]].b],p[i]));
				min = (x<min)?x:min;
				max = (x>max)?x:max;
			}

			if(T[NT[i].t[j]].c!=i)
			{
				x = norm3D(sub3D(p[T[NT[i].t[j]].c],p[i]));
				min = (x<min)?x:min;
				max = (x>max)?x:max;
			}
		}

	for(i=0;i<100;i++)	h[i]=0;

	for(i=0;i<(*m).np;i++)
		for(j=0;j<NT[i].n;j++)
		{
			if(T[NT[i].t[j]].a!=i)
			{
				x = norm3D(sub3D(p[T[NT[i].t[j]].a],p[i]));
				h[(int)(99*(x-min)/(max-min))]++;
			}

			if(T[NT[i].t[j]].b!=i)
			{
				x = norm3D(sub3D(p[T[NT[i].t[j]].b],p[i]));
				h[(int)(99*(x-min)/(max-min))]++;
			}

			if(T[NT[i].t[j]].c!=i)
			{
				x = norm3D(sub3D(p[T[NT[i].t[j]].c],p[i]));
				h[(int)(100*(x-min)/(max-min))]++;
			}
		}

	maxh=0;
	for(i=0;i<100;i++)
		maxh = (h[i]>maxh)?h[i]:maxh;
	
	EraseRect(GetWindowPortBounds(FrontWindow(),&er));
	TextSize(10);
	
	ForeColor(greenColor);
	for(i=0;i<100;i++)
	{
		SetRect(&r, 10+i*barsize,200-150*(h[i]/(float)maxh),10+barsize*(i+1)-1,200);
		PaintRect(&r);
	}
	
	ForeColor(blackColor);
	MoveTo(10,200);LineTo(20+100*barsize,200);
	for(i=0;i<=100;i++)
		if(i%10==0)
		{
			MoveTo(10+barsize*(i+0.5),205);LineTo(10+barsize*(i+0.5),200);
			
			fnum_to_string(i*(max-min)/100.0+min,ndecimal,s);
			MoveTo(10+i*barsize,220);DrawString(s);
		}
		else
		{
			MoveTo(10+barsize*(i+0.5),202);LineTo(10+barsize*(i+0.5),200);
		}
	MoveTo(10,200);LineTo(10,40);
	MoveTo(5,50);LineTo(10,50);
	fnum_to_string(maxh/(float)((*m).nt*2)*100,ndecimal,s);
	MoveTo(15,50);DrawString(s);DrawString("\p%");
}
void em_infoNeighbour(MeshPtr m)
{
	NTriRec		*NT;
	int			i;
	int			h[30],maxh;
	Rect		r,er;
	Str255		s;
	
	int			v=30;
	int			barsize=20;
	int			ndecimal=2;

	NT = msh_getNeighborTrianglesPtr(m);
	
	for(i=0;i<v;i++)	h[i]=0;

	for(i=0;i<(*m).np;i++)
		h[NT[i].n]++;

	maxh=0;
	for(i=0;i<v;i++)
		maxh = (h[i]>maxh)?h[i]:maxh;
	
	EraseRect(GetWindowPortBounds(FrontWindow(),&er));
	TextSize(10);
	
	ForeColor(greenColor);
	for(i=0;i<v;i++)
	{
		SetRect(&r, 10+i*barsize,200-150*(h[i]/(float)maxh),10+barsize*(i+1)-1,200);
		PaintRect(&r);
	}
	
	ForeColor(blackColor);
	MoveTo(10,200);LineTo(20+100*barsize,200);
	for(i=0;i<=v;i++)
	{
		MoveTo(10+i*barsize,205);LineTo(10+i*barsize,200);
		
		NumToString(i,s);
		MoveTo(10+i*barsize,220);DrawString(s);
	}
	MoveTo(10,200);LineTo(10,40);
	MoveTo(5,50);LineTo(10,50);
	fnum_to_string(maxh/(float)((*m).nt*2)*100,ndecimal,s);
	MoveTo(15,50);DrawString(s);DrawString("\p%");
	
	//
	for(i=0;i<v;i++)
	{
		MoveTo(10+barsize*v+70*(i>7),40+(i-8*(i>7))*12);
		NumToString(i,s);
		DrawString(s);	DrawString("\p : ");
		NumToString(h[i],s);
		DrawString(s);
	}
		
}
void em_infoArea(MeshPtr m)
// constructs an histogram of the triangle's area
// Gives also the total area of the surface
{
	float3D		*p;
	int3D		*T;
	NTriRec		*NT;
	int			i;
	float		sum;
	int			h[100],maxh;
	float		x,min,max;
	Rect		r,er;
	Str255		s;
	
	int			barsize=5;
	int			ndecimal=2;

	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);	
	
	// 1. triangles area histogram
	min=1000;
	max=0;
	for(i=0;i<(*m).nt;i++)
	{
		x = tTriArea(m,i);
		min = (x<min)?x:min;
		max = (x>max)?x:max;
	}

	for(i=0;i<100;i++)	h[i]=0;

	for(i=0;i<(*m).nt;i++)
	{
		x = tTriArea(m,i);
		h[(int)(99*(x-min)/(max-min))]++;
	}

	maxh=0;
	for(i=0;i<100;i++)
		maxh = (h[i]>maxh)?h[i]:maxh;
	
	EraseRect(GetWindowPortBounds(FrontWindow(),&er));
	TextSize(10);
	
	ForeColor(greenColor);
	for(i=0;i<100;i++)
	{
		SetRect(&r, 10+i*barsize,200-150*(h[i]/(float)maxh),10+barsize*(i+1)-1,200);
		PaintRect(&r);
	}
	
	ForeColor(blackColor);
	MoveTo(10,200);LineTo(20+100*barsize,200);
	for(i=0;i<=100;i++)
		if(i%10==0)
		{
			MoveTo(10+barsize*(i+0.5),205);LineTo(10+barsize*(i+0.5),200);
			
			fnum_to_string(i*(max-min)/100.0+min,ndecimal,s);
			MoveTo(10+i*barsize,220);DrawString(s);
		}
		else
		{
			MoveTo(10+barsize*(i+0.5),202);LineTo(10+barsize*(i+0.5),200);
		}
	MoveTo(10,200);LineTo(10,40);
	MoveTo(5,50);LineTo(10,50);
	fnum_to_string(maxh/(float)((*m).nt*2)*100,ndecimal,s);
	MoveTo(15,50);DrawString(s);DrawString("\p%");

	// 2. number of triangles with area=0
	sum=0;
	for(i=0;i<(*m).nt;i++)
		sum+=(tTriArea(m,i)==0);
		
	MoveTo(20,160);
	DrawString("\pWith zero area : ");
	NumToString(sum,s);
	DrawString(s);		

	// 3. total surface area
	sum=0;
	for(i=0;i<(*m).nt;i++)
		sum+=tTriArea(m,i);
		
	MoveTo(20,180);
	DrawString("\pTotal surface area : ");
	fnum_to_string(sum,4,s);
	DrawString(s);
	
	// 4. Bounding box size
	updatebbox(m);
	
	MoveTo(20,110);DrawString("\pBounding box:");
	MoveTo(20,120);DrawString("\p X=");
	fnum_to_string((*m).bbox.siz.x-(*m).bbox.ide.x,4,s);	DrawString(s);
	MoveTo(20,130);DrawString("\p Y=");
	fnum_to_string((*m).bbox.siz.y-(*m).bbox.ide.y,4,s);	DrawString(s);
	MoveTo(20,140);DrawString("\p Z=");
	fnum_to_string((*m).bbox.siz.z-(*m).bbox.ide.z,4,s);	DrawString(s);
}

int em_infoTopology(MeshPtr m)
// return the number of invalid points with non-flat neighbor triangles
// (that we cannot sort in one closing list)
// input: the mesh _m_
// output: the number of errors
{
	int		i,n,t;
	NTriRec	*NT;
	int		pstack[SIZESTACK];
	Str255	s;
        Rect	er;
	
	NT = msh_getNeighborTrianglesPtr(m);
	
	EraseRect(GetWindowPortBounds(FrontWindow(),&er));
	
	t=0;
	for(i=0;i<(*m).np;i++)
	{
		n = msh_psort(m,i,pstack);
		
		if(n!=NT[i].n)
		{
			t++;
			NumToString(i,s);
			printstring(s);
		}
	}
	
	MoveTo(10,10);
	DrawString("\perror points : ");
	NumToString(t,s);
	DrawString(s);
	
	return t;	
}

RGBColor em_colourFromColorMap(int cm_number, float cm_index)
{
	RGBColor	rgb={0,0,0};
	int			i0,i1;
	float		alfa;
	
	if(cm_index>2 || cm_index<0)
		return rgb;
	
	i0 = cm_index*31/2.0;
	i1 = i0+1;
	
	alfa = (cm_index*31/2.0-i0)/(i1-i0);
	
	switch(cm_number)
	{
		case kfire:
			rgb.red   = 0xff*g_em_fire_r[i0]*(1-alfa) + 0xff*g_em_fire_r[i1]*alfa;
			rgb.green = 0xff*g_em_fire_g[i0]*(1-alfa) + 0xff*g_em_fire_g[i1]*alfa;
			rgb.blue  = 0xff*g_em_fire_b[i0]*(1-alfa) + 0xff*g_em_fire_b[i1]*alfa;
			break;
		case kice:
			rgb.red   = 0xff*g_em_ice_r[i0]*(1-alfa) + 0xff*g_em_ice_r[i1]*alfa;
			rgb.green = 0xff*g_em_ice_g[i0]*(1-alfa) + 0xff*g_em_ice_g[i1]*alfa;
			rgb.blue  = 0xff*g_em_ice_b[i0]*(1-alfa) + 0xff*g_em_ice_b[i1]*alfa;
			break;
		case kspectrum:
			break;
	}
	
	return rgb;
}
*/
#pragma mark -
#pragma mark ____________________mesh operations
void em_simplify(MeshPtr m)
{
    int		i,iter;
    NTriRec	*NT;
    
    NT = msh_getNeighborTrianglesPtr(m);

    printf(">>Simplify\n");
	for(iter=0;iter<5;iter++)
    {
        em_fusion(m);

        for(i=0;i<(*m).np;i++)
        {
			if(NT[i].n>=12)
			{	em_splitPoint(m,i);
					msh_setNeighborTriangles(m);
			}
        }
		printf("\t%i.%i\n",iter,(*m).np);
    }
}
void em_simplifyPoint(MeshPtr m)
{
	int			i,j;
	
	// randomly choose points to simplify
	j=0;
	do
	{
		i= (*m).np*(rand()/65535.0+0.5);
		
		if(em_deletePoint(m,i))
			j++;
					
	}while(j<(*m).np*g_em_simplify);
}
bool em_deletePoint(MeshPtr m,int indx)
{
	int			i,k,n;
	float3D		*p;
	int3D		*T;
	NTriRec		*NT;
	int			pstack[SIZESTACK];
	int			nstack;
	bool		used=false;

	// 2.sort neighbour points
	// 3.if there are no nested triangles, delete the implied triangles
	// 4.remake the triangles
	// 5.delete the point

	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	
	//2.
	nstack = msh_psort(m, indx,pstack);
	
	//3.
	if(nstack==NT[indx].n)
	{
		used=true;
		
		for(n=0;n<NT[indx].n;n++)
		{
			// delete the triangle
			for(i=NT[indx].t[n];i<(*m).nt;i++)
				T[i] = T[i+1];
			(*m).nt--;

			// update the stack
			for(i=n+1;i<nstack;i++)
				if(NT[indx].t[i]>NT[indx].t[n])
					NT[indx].t[i]--;
		}					
		//4.
		em_take(m,indx,pstack);
		msh_setNeighborTriangles(m);

		//5.
		for(i=indx;i<(*m).np-1;i++)
		{
			p[i] = p[i+1];
			
			for(k=0;k<NT[i+1].n;k++)
			{
				if(T[NT[i+1].t[k]].a==i+1)
					T[NT[i+1].t[k]].a -= 1;
				if(T[NT[i+1].t[k]].b==i+1)
					T[NT[i+1].t[k]].b -= 1;
				if(T[NT[i+1].t[k]].c==i+1)
					T[NT[i+1].t[k]].c -= 1;
			}
		}
		(*m).np--;
		
		msh_setNeighborTriangles(m);
	}
	return(used);
}

void em_take(MeshPtr m,int indx, int *pstack)
{
	int		i,nstack;
	int		xstack[SIZESTACK];
	int		n;
	int3D	*T;
	NTriRec	*NT;
	
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	nstack = NT[indx].n;
	
	n=0;
	while(nstack>2)
	{
		xstack[n++]=pstack[0];
		
		for(i=0;i<nstack-1;i+=2)
		{
			T[(*m).nt].a=pstack[i];
			T[(*m).nt].b=pstack[i+1];
			T[(*m).nt].c=pstack[(i+2)%nstack];
			
			(*m).nt++;
			
			if(i+2<nstack)
				xstack[n++]=pstack[i+2];
		}
		
		for(i=0;i<n;i++)
			pstack[i]=xstack[i];
		nstack=n;
		n=0;
	}
}

#pragma mark _
void em_setFisionThreshold(float fisionthrs)
{
	g_em_fisionthrs = fisionthrs;
}
float em_getFisionThreshold(void)
{
	return g_em_fisionthrs;
}
bool em_fision(MeshPtr m)
{
	int		i,j,k,l;
	float3D	*p;
	NTriRec	*NT;
	int3D	*T;
	int		ps[SIZESTACK];
	bool	ok= true;
	
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);

	for(i=0;i<(*m).np;i++)
	{
		
		if(NT[i].n>SIZESTACK)
		{
                    printf("tip overflow");
                    break;
		}
		
		if(msh_psort(m,i,ps))
		for(j=0;j<NT[i].n;j++)
		{
                    if(	norm3D(sub3D(p[i],p[ps[j]]))>=g_em_fisionthrs )
                    {
                        p[(*m).np]=sca3D(add3D(p[i],p[ps[j]]),0.5);
                        NT[(*m).np].n=0;
                        
                        for(k=0;k<NT[i].n;k++)
                        {
                            if(T[NT[i].t[k]].a==i && T[NT[i].t[k]].b==ps[j])
                            {
                                //1
                                // new triangles
                                T[NT[i].t[k]].b=(*m).np;
                                
                                T[(*m).nt].a=(*m).np;
                                T[(*m).nt].b=ps[j];
                                T[(*m).nt].c=T[NT[i].t[k]].c;
                                
                                // reset the triangle list of the implied points
                                NT[(*m).np].t[NT[(*m).np].n++]=NT[i].t[k];
                                NT[(*m).np].t[NT[(*m).np].n++]=(*m).nt;
                                
                                NT[T[NT[i].t[k]].c].t[NT[T[NT[i].t[k]].c].n++]=(*m).nt;
                                
                                for(l=0;l<NT[ps[j]].n;l++)
                                        if(NT[ps[j]].t[l]==NT[i].t[k])
                                                NT[ps[j]].t[l]=(*m).nt;
                                
                                (*m).nt++;
                            }
                            else if(T[NT[i].t[k]].b==i && T[NT[i].t[k]].c==ps[j])
                            {
                                //2
                                T[NT[i].t[k]].c=(*m).np ;
                                
                                T[(*m).nt].a=T[NT[i].t[k]].a;
                                T[(*m).nt].b=(*m).np;
                                T[(*m).nt].c=ps[j];
                                
                                NT[(*m).np].t[NT[(*m).np].n++]=NT[i].t[k];
                                NT[(*m).np].t[NT[(*m).np].n++]=(*m).nt;
                                
                                NT[T[NT[i].t[k]].a].t[NT[T[NT[i].t[k]].a].n++]=(*m).nt;
                                
                                for(l=0;l<NT[ps[j]].n;l++)
                                        if(NT[ps[j]].t[l]==NT[i].t[k])
                                                NT[ps[j]].t[l]=(*m).nt;
                                
                                (*m).nt++;
                            }
                            else if(T[NT[i].t[k]].c==i && T[NT[i].t[k]].a==ps[j])
                            {
                                //3
                                T[NT[i].t[k]].a=(*m).np;
                                
                                T[(*m).nt].a=ps[j];
                                T[(*m).nt].b=T[NT[i].t[k]].b;
                                T[(*m).nt].c=(*m).np;
                                
                                NT[(*m).np].t[NT[(*m).np].n++]=NT[i].t[k];
                                NT[(*m).np].t[NT[(*m).np].n++]=(*m).nt;
                                
                                NT[T[NT[i].t[k]].b].t[NT[T[NT[i].t[k]].b].n++]=(*m).nt;
                                
                                for(l=0;l<NT[ps[j]].n;l++)
                                        if(NT[ps[j]].t[l]==NT[i].t[k])
                                                NT[ps[j]].t[l]=(*m).nt;
                                
                                (*m).nt++;
                            }
                            else if(T[NT[i].t[k]].a==ps[j] && T[NT[i].t[k]].b==i)
                            {
                                //4
                                T[NT[i].t[k]].a=(*m).np;
                                
                                T[(*m).nt].a=ps[j];
                                T[(*m).nt].b=(*m).np;
                                T[(*m).nt].c=T[NT[i].t[k]].c;
                                
                                NT[(*m).np].t[NT[(*m).np].n++]=NT[i].t[k];
                                NT[(*m).np].t[NT[(*m).np].n++]=(*m).nt;
                                
                                NT[T[NT[i].t[k]].c].t[NT[T[NT[i].t[k]].c].n++]=(*m).nt;
                                
                                for(l=0;l<NT[ps[j]].n;l++)
                                        if(NT[ps[j]].t[l]==NT[i].t[k])
                                                NT[ps[j]].t[l]=(*m).nt;
                                
                                (*m).nt++;
                            }
                            else if(T[NT[i].t[k]].b==ps[j] && T[NT[i].t[k]].c==i)
                            {
                                //5
                                T[NT[i].t[k]].b=(*m).np;
                                
                                T[(*m).nt].a=T[NT[i].t[k]].a;
                                T[(*m).nt].b=ps[j];
                                T[(*m).nt].c=(*m).np;
                                
                                NT[(*m).np].t[NT[(*m).np].n++]=NT[i].t[k];
                                NT[(*m).np].t[NT[(*m).np].n++]=(*m).nt;
                                
                                NT[T[NT[i].t[k]].a].t[NT[T[NT[i].t[k]].a].n++]=(*m).nt;
                                
                                for(l=0;l<NT[ps[j]].n;l++)
                                        if(NT[ps[j]].t[l]==NT[i].t[k])
                                                NT[ps[j]].t[l]=(*m).nt;
                                
                                (*m).nt++;
                            }
                            else if(T[NT[i].t[k]].c==ps[j] && T[NT[i].t[k]].a==i)
                            {
                                //6
                                T[NT[i].t[k]].c=(*m).np;
                                
                                T[(*m).nt].a=(*m).np;
                                T[(*m).nt].b=T[NT[i].t[k]].b;
                                T[(*m).nt].c=ps[j];
                                
                                NT[(*m).np].t[NT[(*m).np].n++]=NT[i].t[k];
                                NT[(*m).np].t[NT[(*m).np].n++]=(*m).nt;
                                
                                NT[T[NT[i].t[k]].b].t[NT[T[NT[i].t[k]].b].n++]=(*m).nt;
                                
                                for(l=0;l<NT[ps[j]].n;l++)
                                        if(NT[ps[j]].t[l]==NT[i].t[k])
                                                NT[ps[j]].t[l]=(*m).nt;
                                
                                (*m).nt++;
                            }

                            if(NT[(*m).np].n>=SIZESTACK)
                            {	ok=false;	break;	}

                            if(NT[T[NT[i].t[k]].a].n>=SIZESTACK)
                            {	ok=false;	break;	}
                            if(NT[T[NT[i].t[k]].b].n>=SIZESTACK)
                            {	ok=false;	break;	}
                            if(NT[T[NT[i].t[k]].c].n>=SIZESTACK)
                            {	ok=false;	break;	}

                            if((*m).nt>g_em_maxtrian-50)
                            {	ok=false;	break;	}

                        }

                        (*m).np++;
                        if((*m).np>g_em_maxpoint-50)
                        {	ok=false;	break;	}
                        }
		}
	}
	if(!ok)
	{
		printf("neighbour or trian or point overflow");
	}
	return(ok);
}

#pragma mark _
void em_setFusionParam(float fusionthrs, float fusionP)
{
	g_em_lengthfusionthrs = fusionthrs;
	g_em_fusionP = fusionP;
}
void em_getFusionParam(float *fusionthrs, float *fusionP)
{
	*fusionthrs = g_em_lengthfusionthrs;
	*fusionP = g_em_fusionP;
}
void em_fusion(MeshPtr m)
{
	int			indx,w,n,i,j;
	float3D		*p;
	int3D		*T;
	NTriRec		*NT;
	
	char		*usedP;
	char		*us, xpointMk=1,pointMk=2;
	char		*newpP;
	int		*np;
	
	bool		found,loop;
	int			i0,i1,a0,a1;
	
	int			i0_pstack[SIZESTACK],i1_pstack[SIZESTACK];
	int			xi,xj;
	//float		xnorm3D,xmax,ymax;
	int			x,y;//xp;
	
	//1. select a triangle side to colapse
	//2. select the triangle that shares the same side
	//3. colapse the side points
	//4. update the point list
	//5. update the triangle list
	
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	
	usedP = calloc((*m).nt,sizeof(char));
	if(usedP==NULL) printf("Out of memory:simplify_line");
	us = usedP;

	newpP = calloc((*m).np,sizeof(int));
	if(newpP==NULL) printf("Out of memory:simplify_line");
	np = (int *)newpP;
	for(i=0;i<(*m).np;i++)
            np[i]=i;

		
	for(indx=0;indx<(*m).nt;indx++)
	{
		//1.
		found = true;
		if(norm3D(sub3D(p[T[indx].a],p[T[indx].b]))<g_em_lengthfusionthrs)
				n = 0;
		else if(norm3D(sub3D(p[T[indx].b],p[T[indx].c]))<g_em_lengthfusionthrs)
				n = 1;
		else if(norm3D(sub3D(p[T[indx].c],p[T[indx].a]))<g_em_lengthfusionthrs)
				n = 2;
		else
		if(tTriArea(m,indx)<g_em_areafusionthrs)
		{
			if(norm3D(sub3D(p[T[indx].a],p[T[indx].b]))<=norm3D(sub3D(p[T[indx].b],p[T[indx].c]))
				&&norm3D(sub3D(p[T[indx].a],p[T[indx].b]))<=norm3D(sub3D(p[T[indx].c],p[T[indx].a])))
					n=0;
			if(norm3D(sub3D(p[T[indx].b],p[T[indx].c]))<=norm3D(sub3D(p[T[indx].a],p[T[indx].b]))
				&&norm3D(sub3D(p[T[indx].b],p[T[indx].c]))<=norm3D(sub3D(p[T[indx].c],p[T[indx].a])))
					n=1;
			if(norm3D(sub3D(p[T[indx].c],p[T[indx].a]))<=norm3D(sub3D(p[T[indx].a],p[T[indx].b]))
				&&norm3D(sub3D(p[T[indx].c],p[T[indx].a]))<=norm3D(sub3D(p[T[indx].b],p[T[indx].c])))
					n=2;
		}
		else
			found=false;
		
		if( found && (rand()/65535.0+0.5)>g_em_fusionP && (!AREAZERO(indx)))
		{
			//2.
			found = false;
			
			switch(n)
			{
				case 0:if(us[T[indx].a] + us[T[indx].b] == 0)
						{	i0 = T[indx].a;
							i1 = T[indx].b;
							a0 = T[indx].c;
							found = true;	}	break;
				case 1:if(us[T[indx].b] + us[T[indx].c] == 0)
						{	i0 = T[indx].b;
							i1 = T[indx].c;
							a0 = T[indx].a;
							found = true;	}	break;
				case 2:if(us[T[indx].c] + us[T[indx].a] == 0)
						{	i0 = T[indx].c;
							i1 = T[indx].a;
							a0 = T[indx].b;
							found = true;	}	break;
			}

			if(found && us[a0]==0)
			{
				for(i=0;i<NT[i0].n;i++)
				{
					if(		NT[i0].t[i] != indx && T[NT[i0].t[i]].a == i1)
					{	w=	NT[i0].t[i];	a1 = T[w].c;	break;}
					if(		NT[i0].t[i] != indx && T[NT[i0].t[i]].b == i1)
					{	w=	NT[i0].t[i];	a1 = T[w].a;	break;}
					if(		NT[i0].t[i] != indx && T[NT[i0].t[i]].c == i1)
					{	w=	NT[i0].t[i];	a1 = T[w].b;	break;}
				}
				
				if(us[a1]==0)
				{
					msh_psort(m,i0, i0_pstack);
					msh_psort(m,i1, i1_pstack);
					
					xi=0;	while(i0_pstack[xi]!=i1)	xi++;
					xj=0;	while(i1_pstack[xj]!=i0)	xj++;
					
					//maxtrian x, y
					x=a0;	y=a1;
					//xmax=0;	ymax=0;		
					for(i=2;i<NT[i0].n-1;i++)
						for(j=2;j<NT[i1].n-1;j++)
							if(np[i0_pstack[(xi+i)%NT[i0].n]] == np[i1_pstack[(xj+j)%NT[i1].n]])
							{
								found = false;
								/*	xp = i0_pstack[(xi+i)%NT[i0].n];
									xnorm3D = norm3D(sub3D(p[xp],sca3D(add3D(p[i0],p[i1]),0.5)));
									if(dot3D(	cross3D(sub3D(p[xp],p[i1]),sub3D(p[i0],p[i1])) ,
												cross3D(sub3D(p[a0],p[i1]),sub3D(p[i0],p[i1])) ) > 0)
									{
										if(xnorm3D>xmax && xp!=a1)
										{	xmax = xnorm3D;
											x = xp;			}
									}
									else
									{
										if(xnorm3D>ymax && xp!=a0)
										{	ymax = xnorm3D;
											y = xp;			}
								}*/
							}
					if(found/*us[x] + us[y] == 0*/)
					{
					//3.
						p[i0] = sca3D(add3D(p[i0],p[i1]),0.5);
							
						p[i1] = p[i0];
						np[i1] = i0;
						
						us[i0] = pointMk;
						us[i1] = xpointMk;
						us[x] = pointMk;
						us[y] = pointMk;
					//i0
						NT[i0].n -= 2;
						n = 0;
						for(i=0;i<NT[i0].n;i++)
						{
							while(NT[i0].t[i+n]==indx || NT[i0].t[i+n]==w)
								n++;
							NT[i0].t[i] = NT[i0].t[i+n];
						}
						for(i=0;i<NT[i1].n;i++)
							if(NT[i1].t[i]!=indx && NT[i1].t[i]!=w)	NT[i0].t[NT[i0].n++] = NT[i1].t[i];

						/*if(us[a0]==0)
							eraseinside_tri_from_point((int3D){i0,i1,x}, a0, us,np,m);
						if(us[a1]==0)
							eraseinside_tri_from_point((int3D){i0,y,i1}, a1, us,np,m);*/
					}
				}
			}
		}
		
	}
	
	//4.
	i=0;
	n=0;
	do
	{
		while(us[i+n]==xpointMk)
			n++;
		
		np[i+n]=i;
		i++;
	}while(i+n<(*m).np );
	
	for(i=0;i<(*m).np;i++)	if(us[i]==xpointMk)
									np[i] = np[np[i]];
	i=0;
	n=0;
	do
	{
		while(us[i+n]==xpointMk)
			n++;
		
		p[i]=p[i+n];
		i++;
	}while(i+n<(*m).np);
	(*m).np -= n;
	
	//5.
	i=0;
	n=0;
	do
	{
		do
		{
			loop=false;

			T[i+n].a = np[T[i+n].a];
			T[i+n].b = np[T[i+n].b];
			T[i+n].c = np[T[i+n].c];
		
			if(AREAZERO(i+n))
			{
				n++;
				
				if(i+n<(*m).nt)
					loop=true;
				else
					break;
			}
		}while(loop);
		
		if(i+n<(*m).nt)
		{
			T[i]=T[i+n];
			i++;
		}
	}while(i+n<(*m).nt);
	(*m).nt -= n;
	
	free(usedP);
	free(newpP);
	
	msh_setNeighborTriangles(m);
}

void em_eraseInsideTriangleFromPoint(int3D t, int sp, char *us, int *np, MeshPtr m)
{
	NTriRec		*NT;
	int		i;
	char		xpointMk=1;
	int		pstack[SIZESTACK];
	
	NT = msh_getNeighborTrianglesPtr(m);

	us[sp] = xpointMk;
	np[sp] = t.a;

	msh_psort(m, sp, pstack);
	for(i=0;i<NT[sp].n;i++)
		if(us[pstack[i]]==0)
			em_eraseInsideTriangleFromPoint(t, pstack[i], us, np, m);
}

#pragma mark _
void em_splitPoint(MeshPtr m, int np)
{
	int			i,n;
	float3D		*p;
	int3D		*T;
	NTriRec		*NT;
	int			nstack;
	int			pstack[SIZESTACK];
	bool		isupdatable;
	
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	
	if(NT[np].n<6)
		return;
	
	nstack = msh_psort(m, np, pstack);
	
	if(nstack!=NT[np].n)
		return;
	
	// one new point, plus a changed one
	p[np] = add3D(sca3D(p[pstack[(int)(nstack*0.75)]],1/3.0),sca3D(p[pstack[(int)(nstack*0.25)]],2/3.0));
	p[(*m).np] = add3D(sca3D(p[pstack[(int)(nstack*0.75)]],2/3.0),sca3D(p[pstack[(int)(nstack*0.25)]],1/3.0));
	
	// two new triangles
	T[(*m).nt].a = pstack[0];
	T[(*m).nt].b = np;
	T[(*m).nt].c = (*m).np;
	(*m).nt += 1;
	
	T[(*m).nt].a = pstack[(int)(nstack*0.5)];
	T[(*m).nt].b = (*m).np;
	T[(*m).nt].c = np;
	(*m).nt += 1;
	
	// update the triangles of the side of the new point
	for(i=0;i<NT[np].n;i++)
	{
		if(	T[NT[np].t[i]].a==np )
		{
			isupdatable = false;
			for(n=(int)(nstack*0.5);n<nstack;n++)	if(T[NT[np].t[i]].b==pstack[n])
													{
														isupdatable=true;
														break;
													}
			if(isupdatable)
				T[NT[np].t[i]].a = (*m).np;
		}
		else
		if(	T[NT[np].t[i]].b==np )
		{
			isupdatable = false;
			for(n=(int)(nstack*0.5);n<nstack;n++)	if(T[NT[np].t[i]].c==pstack[n])
													{
														isupdatable=true;
														break;
													}
			if(isupdatable)
				T[NT[np].t[i]].b = (*m).np;
		}
		else
		if(	T[NT[np].t[i]].c==np )
		{
			isupdatable = false;
			for(n=(int)(nstack*0.5);n<nstack;n++)	if(T[NT[np].t[i]].a==pstack[n])
													{
														isupdatable=true;
														break;
													}
			if(isupdatable)
				T[NT[np].t[i]].c = (*m).np;
		}
	}
	
	(*m).np += 1;
}
#pragma mark _

void em_flipEdges(MeshPtr m, int triangle)
{
	int		i=triangle;
	int		tn[3],ts[3];
	float3D	*p;
	int3D	*T;
	NTriRec	*NT;
	float	R,x;
	bool	changed;
	
	R=fabs((*m).bbox.siz.x-(*m).center.x);
	x=2*R*sqrt(pi/sqrt(3)/(*m).nt);
	
	p=msh_getPointsPtr(m);
	T=msh_getTrianglesPtr(m);
	NT=msh_getNeighborTrianglesPtr(m);

	msh_get_triangle_neighbortriangles(m,i, tn, ts);
	
	changed=false;
	if(NT[T[i].a].n>3 && NT[T[i].b].n>3)
	switch(ts[0])
	{
		case 1:	if(norm3D(sub3D(p[T[i].c],p[T[tn[0]].c]))<norm3D(sub3D(p[T[i].a],p[T[i].b])))
				{	T[i].b=T[tn[0]].c;
					T[tn[0]].b=T[i].c;
					changed=true;
				}
				break;
		case 3:	if(norm3D(sub3D(p[T[i].c],p[T[tn[0]].b]))<norm3D(sub3D(p[T[i].a],p[T[i].b])))
				{	T[i].b=T[tn[0]].b;
					T[tn[0]].a=T[i].c;
					changed=true;
				}
				break;
		case 2:	if(norm3D(sub3D(p[T[i].c],p[T[tn[0]].a]))<norm3D(sub3D(p[T[i].a],p[T[i].b])))
				{	T[i].b=T[tn[0]].a;
					T[tn[0]].c=T[i].c;
					changed=true;
				}
				break;
	}
	if(!changed)
	if(NT[T[i].b].n>3 && NT[T[i].c].n>3)
	switch(ts[1])
	{
		case 1:	if(norm3D(sub3D(p[T[i].a],p[T[tn[1]].c]))<norm3D(sub3D(p[T[i].b],p[T[i].c])))
				{	T[i].c=T[tn[1]].c;
					T[tn[1]].b=T[i].a;
					changed=true;
				}
				break;
		case 2:	if(norm3D(sub3D(p[T[i].a],p[T[tn[1]].a]))<norm3D(sub3D(p[T[i].b],p[T[i].c])))
				{	T[i].c=T[tn[1]].a;
					T[tn[1]].c=T[i].a;
					changed=true;
				}
				break;
		case 3:	if(norm3D(sub3D(p[T[i].a],p[T[tn[1]].b]))<norm3D(sub3D(p[T[i].b],p[T[i].c])))
				{	T[i].c=T[tn[1]].b;
					T[tn[1]].a=T[i].a;
					changed=true;
				}
				break;
	}
	if(!changed)
	if(NT[T[i].c].n>3 && NT[T[i].a].n>3)
	switch(ts[2])
	{
		case 1:	if(norm3D(sub3D(p[T[i].b],p[T[tn[2]].c]))<norm3D(sub3D(p[T[i].c],p[T[i].a])))
				{	T[i].a=T[tn[2]].c;
					T[tn[2]].b=T[i].b;
					changed=true;
				}
				break;
		case 2:	if(norm3D(sub3D(p[T[i].b],p[T[tn[2]].a]))<norm3D(sub3D(p[T[i].c],p[T[i].a])))
				{	T[i].a=T[tn[2]].a;
					T[tn[2]].c=T[i].b;
					changed=true;
				}
				break;
		case 3:	if(norm3D(sub3D(p[T[i].b],p[T[tn[2]].b]))<norm3D(sub3D(p[T[i].c],p[T[i].a])))
				{	T[i].a=T[tn[2]].b;
					T[tn[2]].a=T[i].b;
					changed=true;
				}
				break;
	}

	if(changed)
		msh_setNeighborTriangles(m);
}
/*void em_printPointInformation(MeshPtr m, int point)
{
	int		i;
	Str255	s,s1;
	float3D	*p;
	int3D	*T;
	NTriRec	*NT;
	
	p=msh_getPointsPtr(m);
	T=msh_getTrianglesPtr(m);
	NT=msh_getNeighborTrianglesPtr(m);
	
	printstring("\p>> ");
	NumToString(point,s);DrawString(s);

	for(i=0;i<NT[point].n;i++)
	{
		NumToString(NT[point].t[i],s);
		psadd(s,"\p: ");
		NumToString(T[NT[point].t[i]].a,s1);
		psadd(s,s1);psadd(s,"\p,");
		NumToString(T[NT[point].t[i]].b,s1);
		psadd(s,s1);psadd(s,"\p,");
		NumToString(T[NT[point].t[i]].c,s1);
		psadd(s,s1);
		
		printstring(s);
	}
}*/
#pragma mark _

void em_setSmooth(float smoothfactor)
{
    g_em_smooth = smoothfactor;
}
float em_getSmooth(void)
{
    return g_em_smooth;
}

void em_smooth(MeshPtr m) // Laplace smooth
{
	float3D	*tmp,x,dx;
	int		*n;
	int		i;
	float3D	*p;
    int3D	*t;
	
    p = msh_getPointsPtr(m);
    t = msh_getTrianglesPtr(m);
    tmp=(float3D*)calloc(m->np,sizeof(float3D));
    n=(int*)calloc(m->np,sizeof(int));
    for(i=0;i<m->nt;i++)
    {
    	tmp[t[i].a]=add3D(tmp[t[i].a],add3D(p[t[i].b],p[t[i].c]));
    	tmp[t[i].b]=add3D(tmp[t[i].b],add3D(p[t[i].c],p[t[i].a]));
    	tmp[t[i].c]=add3D(tmp[t[i].c],add3D(p[t[i].a],p[t[i].b]));
    	n[t[i].a]+=2;
    	n[t[i].b]+=2;
    	n[t[i].c]+=2;
    }
    for(i=0;i<m->np;i++)
    {
    	x=sca3D(tmp[i],1/(float)n[i]);
    	dx=sub3D(x,p[i]);
    	p[i]=add3D(p[i],sca3D(dx,g_em_smooth));	// p=p+l(x-p)
    }
    free(tmp);
    free(n);
}
#pragma mark _
void em_setInflate(float t)
{
    g_em_inflate=t;
}
float em_getInflate(void)
{
    return g_em_inflate;
}
void em_inflate(MeshPtr m)
{
	int		i,j;
	float3D	*p,xp,nn,zero={0,0,0};
	NTriRec	*NT;
	
	p = msh_getPointsPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	
	for(i=0;i<(*m).np;i++)
	{
		nn = zero;
		for(j=0;j<NT[i].n;j++)
			nn = add3D(nn,tTriPlane(NT[i].t[j],m));
		
			xp = sca3D(nn,1/norm3D(nn));
		
		p[i] = add3D(p[i],sca3D(xp,g_em_inflate));
	}
}

#pragma mark _
void em_setScale(float scale)
{
    g_em_scale = scale;
}
float em_getScale(void)
{
    return g_em_scale;
}
void em_scale(MeshPtr m)
{
    int			i;
    float3D		*p;
    int3D		*T;
    
    p = (*m).p;
    T = (*m).t;
    
    for(i=0;i<(*m).np;i++)
            p[i]= add3D( (*m).center, sca3D( sub3D(p[i],(*m).center) , g_em_scale ) );
    
    (*m).bbox.siz=p[0];
    (*m).bbox.ide=p[0];
    for(i=0;i<(*m).np;i++)
    {
            if(p[i].x<(*m).bbox.siz.x)
                    (*m).bbox.siz.x=p[i].x;
            if(p[i].y<(*m).bbox.siz.y)
                    (*m).bbox.siz.y=p[i].y;
            if(p[i].z<(*m).bbox.siz.z)
                    (*m).bbox.siz.z=p[i].z;
            
            if(p[i].x>(*m).bbox.ide.x)
                    (*m).bbox.ide.x=p[i].x;
            if(p[i].y>(*m).bbox.ide.y)
                    (*m).bbox.ide.y=p[i].y;
            if(p[i].z>(*m).bbox.ide.z)
                    (*m).bbox.ide.z=p[i].z;
    }
}
#pragma mark _
void em_centre(MeshPtr m)
{
    int			i;
    float3D		*p;
    
    p = (*m).p;
    
    update(m);
	for(i=0;i<(*m).np;i++)
            p[i]= sub3D(p[i],(*m).center);
    update(m);
}
void em_translate(MeshPtr m,float3D t)
{
    int			i;
    float3D		*p;
    
    p = (*m).p;
    
    update(m);
	for(i=0;i<(*m).np;i++)
		p[i]= sub3D(p[i],t);
    update(m);
}
void em_flipTriangles(MeshPtr m)
{
    int			i;
    int3D		*t;
    
    t = (*m).t;
    
	for(i=0;i<(*m).nt;i++)
            t[i]= (int3D){t[i].a,t[i].c,t[i].b};
}
void em_adaptMesh(MeshPtr m)
{
    int			i,j,a,b,c;
    int         np1,nt1;
    float3D     *p,*p1,*C,*C1,pab,pbc,pca,pabc;
    int3D		*t,*t1;
    int         *T,flag;
    MeshPtr     m1;
    
	p = msh_getPointsPtr(m);
	t = msh_getTrianglesPtr(m);
	C = msh_getTexturePtr(m);
    
    // count number of new points and triangles
	np1=m->np;
    nt1=m->nt;
    for(i=0;i<m->nt;i++)
    {
        T=(int*)&(t[i]);
        flag=0;
        for(j=0;j<3;j++)
            if(!equals3D(C[T[j]],C[T[(j+1)%3]]))
                flag++;
        
        if(flag==2) // triangle contains a frontier between two classes
        {
            np1+=4;
            nt1+=2;
        }
        else
        if(flag==3) // triangle contains a frontier between three classes
        {
            np1+=9;
            nt1+=6;
        }
    }
    msh_new(&m1, np1, nt1);
	p1 = msh_getPointsPtr(m1);
	t1 = msh_getTrianglesPtr(m1);
	C1 = msh_getTexturePtr(m1);
    for(i=0;i<m->np;i++)
    {
        p1[i]=p[i];
        C1[i]=C[i];
    }
    for(i=0;i<m->nt;i++)
        t1[i]=t[i];

    // create adapted mesh
    np1=m->np;
    nt1=m->nt;
	for(i=0;i<m->nt;i++)
    {
        T=(int*)&(t[i]);
        
        // check 3-frontier
        if(!equals3D(C[T[0]],C[T[1]]) && !equals3D(C[T[1]],C[T[2]]) && !equals3D(C[T[2]],C[T[0]]))
        {
            // new vertex positions
            pab=sca3D(add3D(p[T[0]],p[T[1]]),0.5);
            pbc=sca3D(add3D(p[T[1]],p[T[2]]),0.5);
            pca=sca3D(add3D(p[T[2]],p[T[0]]),0.5);
            pabc=sca3D(add3D(add3D(p[T[0]],p[T[1]]),p[T[2]]),1/3.0);
            
            // new vertices
            p1[np1+0]=pab;
            p1[np1+1]=pab;
            p1[np1+2]=pbc;
            p1[np1+3]=pbc;
            p1[np1+4]=pca;
            p1[np1+5]=pca;
            p1[np1+6]=pabc;
            p1[np1+7]=pabc;
            p1[np1+8]=pabc;
            
            // new vertex values
            C1[np1+0]=C[T[0]];
            C1[np1+1]=C[T[1]];
            C1[np1+2]=C[T[1]];
            C1[np1+3]=C[T[2]];
            C1[np1+4]=C[T[2]];
            C1[np1+5]=C[T[0]];
            C1[np1+6]=C[T[0]];
            C1[np1+7]=C[T[1]];
            C1[np1+8]=C[T[2]];
            
            // modified and new triangles
            t1[i]=(int3D){T[0],np1+0,np1+6};
            t1[nt1+0]=(int3D){T[0],np1+6,np1+4};
            t1[nt1+1]=(int3D){T[1],np1+2,np1+7};
            t1[nt1+2]=(int3D){T[1],np1+7,np1+1};
            t1[nt1+3]=(int3D){T[2],np1+5,np1+8};
            t1[nt1+4]=(int3D){T[2],np1+8,np1+3};
            
            np1+=9;
            nt1+=6;

            continue;
        }
        
        // check 2-frontier
        flag=0;
        for(j=0;j<3;j++)
            if(!equals3D(C[T[j]],C[T[(j+1)%3]]))
            {
                flag=1;
                break;
            }
        if(flag==0)
            continue;
        if(equals3D(C[T[(j+1)%3]],C[T[(j+2)%3]]))
        {
            // single point is j: assign it to a
            a=T[j];
            b=T[(j+1)%3];
            c=T[(j+2)%3];
            
        }
        else
        {
            // single point is j+1: assign it to a
            a=T[(j+1)%3];
            b=T[(j+2)%3];
            c=T[j];
        }
        
        // new points: p_ab and p_ac (these points are duplicated [which will change the surface's topology])
        pab=sca3D(add3D(p[a],p[b]),0.5);
        pca=sca3D(add3D(p[a],p[c]),0.5);
        p1[np1+0]=pab;
        p1[np1+1]=pab;
        p1[np1+2]=pca;
        p1[np1+3]=pca;
        C1[np1+0]=C[a];
        C1[np1+1]=C[b];
        C1[np1+2]=C[a];
        C1[np1+3]=C[b];
        
        // modified triangle {a,pab,pac} with C=C[a], and new triangles {pab,b,c} and {c,pac,pab} with C=C[b] (=C[c])
        t1[i]=(int3D){a,np1,np1+2};
        t1[nt1]=(int3D){np1+1,b,c};
        t1[nt1+1]=(int3D){c,np1+3,np1+1};
        
        np1+=4;
        nt1+=2;
    }

    free(m->p);
    free(m->t);
    msh_deleteTexturePtr(m);
    
    m->np=np1;
    m->nt=nt1;
    m->p=m1->p;
    m->t=m1->t;
    m->xcontainer=m1->xcontainer;

}
#pragma mark -
/*
void em_pickPoint(MeshPtr m)
{
	int3D		*T;
	float3D		*p;
	NTriRec		*NT;
	int			*arr;
	int			i;
	Point		mouse,x;
	Rect		r;
	
	if(g_em_arrH==NULL)
		return;

	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	arr = (int *)*g_em_arrH;
	
	GetMouse(&mouse);
	
	for(i=(*m).nt;i>=0;i--)
	{
		x=float3DTo2D(p[T[arr[i]].a]);
		SetRect(&r,x.h-2,x.v-2,x.h+2,x.v+2);
		if(PtInRect(mouse,&r))
		{
			g_em_isSelected = true;
			g_em_selectedpoint = T[arr[i]].a;
			break;
		}
		
		x=float3DTo2D(p[T[arr[i]].b]);
		SetRect(&r,x.h-2,x.v-2,x.h+2,x.v+2);
		if(PtInRect(mouse,&r))
		{
			g_em_isSelected = true;
			g_em_selectedpoint = T[arr[i]].b;
			break;
		}
		
		x=float3DTo2D(p[T[arr[i]].c]);
		SetRect(&r,x.h-2,x.v-2,x.h+2,x.v+2);
		if(PtInRect(mouse,&r))
		{
			g_em_isSelected = true;
			g_em_selectedpoint = T[arr[i]].c;
			break;
		}
	}
	
	drawselected(m);
}
void em_drawSelected(MeshPtr m)
{
	int3D		*T;
	float3D		*p;
	NTriRec		*NT;
	int			j;
	Point		x;
	Rect		r;
	Str255		s;

	if(m==NULL)
		return;
	
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	
	ForeColor(greenColor);
	for(j=0;j<NT[g_em_selectedpoint].n;j++)
	{
		if(T[NT[g_em_selectedpoint].t[j]].a==g_em_selectedpoint)
		{
			MoveTo3D(p[g_em_selectedpoint]);	LineTo3D(p[T[NT[g_em_selectedpoint].t[j]].b]);
			MoveTo3D(p[g_em_selectedpoint]);	LineTo3D(p[T[NT[g_em_selectedpoint].t[j]].c]);
		}
		if(T[NT[g_em_selectedpoint].t[j]].b==g_em_selectedpoint)
		{
			MoveTo3D(p[g_em_selectedpoint]);	LineTo3D(p[T[NT[g_em_selectedpoint].t[j]].a]);
			MoveTo3D(p[g_em_selectedpoint]);	LineTo3D(p[T[NT[g_em_selectedpoint].t[j]].c]);
		}
		if(T[NT[g_em_selectedpoint].t[j]].c==g_em_selectedpoint)
		{
			MoveTo3D(p[g_em_selectedpoint]);	LineTo3D(p[T[NT[g_em_selectedpoint].t[j]].a]);
			MoveTo3D(p[g_em_selectedpoint]);	LineTo3D(p[T[NT[g_em_selectedpoint].t[j]].b]);
		}
	}
	ForeColor(redColor);
	x=float3DTo2D(p[g_em_selectedpoint]);
	SetRect(&r,x.h-2,x.v-2,x.h+2,x.v+2);
	PaintRect(&r);
	ForeColor(whiteColor);
	NumToString(g_em_selectedpoint,s);
	MoveTo(x.h,x.v);DrawString(s);
}
*/
int	em_getSelectedPoint(void)
{
	return g_em_selectedpoint;
}
#pragma mark _
float em_laplaceFromAngle(float angle)
{
	/*float	l0,l1,laplace;
	int		a0,a1;
	
	angle*=511/pi;
	a0=(int)angle;
	a1=a0+1;
	
	l0=g_em_harmonics[a0];
	l1=g_em_harmonics[a1];
	
	laplace=l0+(l1-l0)*(angle-a0)/(a1-a0);*/
	float	laplace;
	float	rpos,rneg;
	float	ex=0.39;
	
	rpos=sqrt(pow(sin(angle),2)+pow(1-cos(angle),2));
	rneg=sqrt(pow(sin(angle),2)+pow(1+cos(angle),2));
	
	laplace= -(pow(rpos,ex)-pow(rneg,ex))/(pow(rpos,ex)+pow(rneg,ex));

	return laplace;
}
float em_angleFromLaplace(float laplace)
{
	int		i;
	float	l0,l1;
	int		a0,a1;
	float	angle;
	
	float	rpos,rneg;
	float	ex=0.39;
	
	a0=510;a1=511;
	for(i=0;i<511;i++)
		if(g_em_harmonics[i]<=laplace)
		{
			a0=i-1;
			a1=i;
			break;
		}
	
	l0=g_em_harmonics[a0];
	l1=g_em_harmonics[a1];
	
	//TEST
	rpos=sqrt(pow(sin(a0*pi/511.0),2)+pow(1-cos(a0*pi/511.0),2));
	rneg=sqrt(pow(sin(a0*pi/511.0),2)+pow(1+cos(a0*pi/511.0),2));
	l0= -(pow(rpos,ex)-pow(rneg,ex))/(pow(rpos,ex)+pow(rneg,ex));
	rpos=sqrt(pow(sin(a1*pi/511.0),2)+pow(1-cos(a1*pi/511.0),2));
	rneg=sqrt(pow(sin(a1*pi/511.0),2)+pow(1+cos(a1*pi/511.0),2));
	l1= -(pow(rpos,ex)-pow(rneg,ex))/(pow(rpos,ex)+pow(rneg,ex));

	angle=a0*pi/511.0 + (a1-a0)*(pi/511.0)*(laplace-l0)/(l1-l0);
	
	return angle;
}
void em_sphereFromTxtr(MeshPtr m, float3D *C)
{
	int			i;
	float3D		*p;
	int3D		*T;
	float		area,R;
	int			maxy,maxx,minz;
	float3D		a0,a1,a2,b0,b1,b2;

	p=msh_getPointsPtr(m);
	T=msh_getTrianglesPtr(m);
	
	// The following code is being programed to roughly align the chimpanzee
	// surfaces in /Users/roberto/Documents/2007_10Primates/2008_10Chimpanzees-Mangin/leschimp/01 LmeshChimp/
	// from the configurations before and after spherical deformation.
	// The alignment is performed based on the occipital pole (+Y in the data), the medial-most tip (+X),
	// and the top of the dorsal margin (-Z)
	maxy=maxx=minz=0;
	for(i=0;i<(*m).np;i++)
	{
		if(p[maxy].y<p[i].y)
			maxy=i;
		if(p[maxx].x<p[i].x)
			maxx=i;
		if(p[minz].z>p[i].z)
			minz=i;
	}
	a0=p[maxy];
	a1=p[maxx];
	a2=p[minz];
	
	area=0;
	for(i=0;i<(*m).nt;i++)
		area+=tTriArea(m,i);
	area=area/3.0;
	R=sqrt(area/pi);
	for(i=0;i<(*m).np;i++)
	{
            p[i]=em_getPointFromSphericalCoordinate(C[i]);
            p[i]=sca3D((float3D){-p[i].x,p[i].y,-p[i].z},R);
            p[i]=add3D(p[i],(*m).center);
	}
	b0=p[maxy];
	b1=p[maxx];
	b2=p[minz];
	
	
	// Compute the rotation necessary to pass from b0,b1,b2 to a0,a1,a2.
	// First, find orthonormal basis for as and bs
	// Then, the rotation matrix is A=B*r, r=B'*A
	float3D	A0,A1,A2,B0,B1,B2;
	float	r[9];
	
	A0=sub3D(a1,a0);
	A0=sca3D(A0,1/norm3D(A0));
	A2=cross3D(A0,sub3D(a2,a0));
	A2=sca3D(A2,1/norm3D(A2));
	A1=cross3D(A2,A0);
	B0=sub3D(b1,b0);
	B0=sca3D(B0,1/norm3D(B0));
	B2=cross3D(B0,sub3D(b2,b0));
	B2=sca3D(B2,1/norm3D(B2));
	B1=cross3D(B2,B0);
	r[0]=B0.x*A0.x+B1.x*A1.x+B2.x*A2.x;
	r[1]=B0.x*A0.y+B1.x*A1.y+B2.x*A2.y;
	r[2]=B0.x*A0.z+B1.x*A1.z+B2.x*A2.z;
	r[3]=B0.y*A0.x+B1.y*A1.x+B2.y*A2.x;
	r[4]=B0.y*A0.y+B1.y*A1.y+B2.y*A2.y;
	r[5]=B0.y*A0.z+B1.y*A1.z+B2.y*A2.z;
	r[6]=B0.z*A0.x+B1.z*A1.x+B2.z*A2.x;
	r[7]=B0.z*A0.y+B1.z*A1.y+B2.z*A2.y;
	r[8]=B0.z*A0.z+B1.z*A1.z+B2.z*A2.z;
	for(i=0;i<(*m).np;i++)
		p[i]=matfloat3D(r,p[i]);
}
float3D em_getPointFromSphericalCoordinate(float3D c)
{
    float	n;
    float3D	p;

    p.x=cos(em_angleFromLaplace(c.x/2.0));
    p.y=cos(em_angleFromLaplace(c.y/2.0));
    p.z=cos(em_angleFromLaplace(c.z/2.0));
    
    n=sqrt(p.x*p.x+p.y*p.y+p.z*p.z); //n=1; -> not norm3Dalized
    p=sca3D(p,1/n);
    
    return p;
}
#pragma mark _
void em_fusionMeshToMeshEdgeCurve(MeshPtr *mm,MeshPtr m,MeshEdgeCurveRec MEC)
{
	int		npoints,ntrian;
	
	// get the number of necessary extra points and triangles
	em_setFusionMeshToMeshEdgeCurve(mm,m,MEC,&npoints,&ntrian,false);
	em_setFusionMeshToMeshEdgeCurve(mm,m,MEC,&npoints,&ntrian,true);
}
void em_setFusionMeshToMeshEdgeCurve(MeshPtr *mm,MeshPtr m,MeshEdgeCurveRec MEC,int *npoints,int *ntrian,bool isStoring)
{
	int			i;
	int			t,side,a,b,c,x0,x1,s[4];
	MeshEdgePointRec	*mep,mep0,mep1;
	float3D			*p;
	int3D			*T;
	int			*index;
	
	mep=MEC.p;
	index=ivector_new(MEC.np);

	if(isStoring)
	{
		msh_new(mm,(*m).np+*npoints,(*m).nt+*ntrian);
		p=msh_getPointsPtr(*mm);
		T=msh_getTrianglesPtr(*mm);

		for(i=0;i<(*m).np;i++)	p[i]=((*m).p)[i];
		for(i=0;i<(*m).nt;i++)	T[i]=((*m).t)[i];
	}
	
	// 1. create the new points, save an index.
	for(i=0;i<MEC.np;i++)
	{
		if(isStoring)
			p[(*m).np+i]=msh_meshedgepoint_to_float3D(m,mep[i]);
		index[i]=(*m).np+i;
	}
	*npoints=MEC.np;
	
	// TEST
	//if(isStoring)
	//{
	//	(***mm).np=(*m).np;
	//	topologyinfo(*mm);
	//}
	
	// 2. create the new triangles, modify the old
	*ntrian=0;
	for(i=0;i<MEC.ne;i++)
	{
		mep0=mep[MEC.e[i].a];
		mep1=mep[MEC.e[i].b];
		
		if(mep0.ta<0 && mep1.ta<0)
		{
			// the two points are vertices: keep just the coordinates
		}
		else
		if(mep0.ta<0 && mep1.ta>=0)
		{
			// 1st a vertex, then a side: create 1 triangle, modify 1 triangle
			// (dividing the side of mep1, including the opposite point)
			if(isStoring)
			{
				// determine the triangle to modify
				t=-1;
				if(T[mep1.ta].a==mep0.tb || T[mep1.ta].b==mep0.tb || T[mep1.ta].c==mep0.tb){	t=mep1.ta; side=mep1.sa;}
				if(T[mep1.tb].a==mep0.tb || T[mep1.tb].b==mep0.tb || T[mep1.tb].c==mep0.tb){	t=mep1.tb; side=mep1.sb;}
				
				if(t>=0)
				{
					// create & modify
                                    T[(*m).nt+*ntrian]=T[t];
                                    switch(side)
                                    {
                                        case 1:	T[t].b=T[(*m).nt+*ntrian].a = index[(MEC.e)[i].b];	break;
                                        case 2:	T[t].c=T[(*m).nt+*ntrian].b = index[(MEC.e)[i].b];	break;
                                        case 3:	T[t].a=T[(*m).nt+*ntrian].c = index[(MEC.e)[i].b];	break;
                                    }
				}
				else
                                    printf("A fusion_mesh_meshedgecurve:EDGE OUTSIDE TRIANGLE\n");
			}
			(*ntrian)++;
		}
		else
		if(mep0.ta>=0 && mep1.ta<0)
		{
			// 1st a side, then a vertex: create 1 triangle, modify 1 triangle
			if(isStoring)
			{
                            // determine the triangle to modify
                            t=-1;
                            if(T[mep0.ta].a==mep1.tb || T[mep0.ta].b==mep1.tb || T[mep0.ta].c==mep1.tb){	t=mep0.ta;	side=mep0.sa;}
                            if(T[mep0.tb].a==mep1.tb || T[mep0.tb].b==mep1.tb || T[mep0.tb].c==mep1.tb){	t=mep0.tb;	side=mep0.sb;}
                            
                            if(t>=0)
                            {
                                // create & modify
                                T[(*m).nt+*ntrian]=T[t];
                                switch(side)
                                {
                                    case 1:	T[t].b=T[(*m).nt+*ntrian].a = index[(MEC.e)[i].a];	break;
                                    case 2:	T[t].c=T[(*m).nt+*ntrian].b = index[(MEC.e)[i].a];	break;
                                    case 3:	T[t].a=T[(*m).nt+*ntrian].c = index[(MEC.e)[i].a];	break;
                                }
                            }
                            else
                                printf("B fusion_mesh_meshedgecurve:EDGE OUTSIDE TRIANGLE\n");
			}
			(*ntrian)++;
		}
		else
		{
			// two sides: create 2 triangles, modify 1 triangle
			if(isStoring)
			{
				// determine the triangle to modify
				t=-1;
				if(mep0.ta==mep1.ta || mep0.ta==mep1.tb){	t=mep0.ta;}
				if(mep0.tb==mep1.ta || mep0.tb==mep1.tb){	t=mep0.tb;}
				
				if(t>=0)
				{
					s[1]=s[2]=s[3]=-1;
					if(mep0.ta==t) s[mep0.sa]=index[(MEC.e)[i].a];	else	s[mep0.sb]=index[(MEC.e)[i].a];
					if(mep1.ta==t) s[mep1.sa]=index[(MEC.e)[i].b];	else	s[mep1.sb]=index[(MEC.e)[i].b];
					
					// create & modify
					/*if(*ntrian>=780)
					{
						PenSize(3,3);
						ForeColor(greenColor);MoveTo3D(p[T[t].a]);LineTo3D(p[T[t].b]);LineTo3D(p[T[t].c]);LineTo3D(p[T[t].a]);
						ForeColor(blackColor);FrameRect3D(msh_meshedgepoint_to_float3D(m,mep[(*MEC.e)[i].a]),0.5);
											  FrameRect3D(msh_meshedgepoint_to_float3D(m,mep[(*MEC.e)[i].b]),0.5);
					}*/
					T[(*m).nt+*ntrian+1]=T[(*m).nt+*ntrian]=T[t];
					side=1*(s[1]<0)+2*(s[2]<0)+3*(s[3]<0);
					switch(side)
					{
						case 1:	a=T[t].a;b=T[t].b;c=T[t].c;		x0=s[3]; x1=s[2];	break;
						case 2:	a=T[t].b;b=T[t].c;c=T[t].a;		x0=s[1]; x1=s[3];	break;
						case 3:	a=T[t].c;b=T[t].a;c=T[t].b;		x0=s[2]; x1=s[1];	break;
					}
					if(((s[1]<0)+(s[2]<0)+(s[3]<0))>=2)
						printf("C fusion_mesh_meshedgecurve:REPEATED EDGE-POINT\n");
	
					T[(*m).nt+*ntrian+1] =(int3D){a,  b,  x1};
					T[(*m).nt+*ntrian]	  =(int3D){a,  x1, x0};
					T[t]				  =(int3D){x0, x1,  c};

					/*if(*ntrian>=780)
					{
						PenSize(1,1);
						ForeColor(redColor);MoveTo3D(p[T[t].a]);LineTo3D(p[T[t].b]);LineTo3D(p[T[t].c]);LineTo3D(p[T[t].a]);
						ForeColor(yellowColor);MoveTo3D(p[T[(*m).nt+*ntrian].a]);LineTo3D(p[T[(*m).nt+*ntrian].b]);LineTo3D(p[T[(*m).nt+*ntrian].c]);LineTo3D(p[T[(*m).nt+*ntrian].a]);
						ForeColor(magentaColor);MoveTo3D(p[T[(*m).nt+*ntrian+1].a]);LineTo3D(p[T[(*m).nt+*ntrian+1].b]);LineTo3D(p[T[(*m).nt+*ntrian+1].c]);LineTo3D(p[T[(*m).nt+*ntrian+1].a]);
					}*/
				}
				else
					printf("D fusion_mesh_meshedgecurve:EDGE OUTSIDE TRIANGLE\n");
				
			}
			(*ntrian)+=2;
		}
		//ForeColor(redColor);
		//MoveTo3D(msh_meshedgepoint_to_float3D(m,mep0));
		//LineTo3D(msh_meshedgepoint_to_float3D(m,mep1));
		
	}
	//ForeColor(blackColor);
	
	// TEST
	if(isStoring)
	{
		//topologyinfo(*mm);
	}
}
void em_alignSVD(MeshPtr m)
{
    float3D *p,cp={0,0,0},min,max,eve[3];
	int		maxp=100;
	int		i,j,step=(*m).np/maxp;
    float	**a,*w,**v;
    float	eva[3];
    int		ie[3];
    
    p=msh_getPointsPtr(m);

	// centroid
    for(i=0;i<(*m).np;i++)
		cp=add3D(cp,p[i]);
	cp=sca3D(cp,1/(float)(*m).np);

    // orientation
    a=fmatrix_new(3+1,maxp+1);				// allocate memory
    w=fvector_new(maxp+1);
    v=fmatrix_new(maxp+1,maxp+1);
    
    j=1;
	for(i=0;i<(*m).np;i+=step)		// mean-centred data points
    {
        a[1][j]=p[i].x-cp.x;
        a[2][j]=p[i].y-cp.y;
		a[3][j]=p[i].z-cp.z;
		j++;
    }
    
    svdcmp(a,3,maxp,w,v);			// SVD
    
    ie[0]=ie[1]=ie[2]=-1;
    eva[0]=eva[1]=eva[2]=0;			// find the three largest eigenvectors
    for(i=1;i<=maxp;i++)
    {
        if(fabs(w[i])>fabs(eva[2]))
		{   if(fabs(w[i])>fabs(eva[1]))
            {   if(fabs(w[i])>fabs(eva[0]))
				{   eva[2]=eva[1];  ie[2]=ie[1];
					eva[1]=eva[0];  ie[1]=ie[0];
					eva[0]=w[i];	ie[0]=i;
				}
				else
				{   eva[2]=eva[1];  ie[2]=ie[1];
					eva[1]=w[i];	ie[1]=i;
				}
            }
            else
            {	eva[2]=w[i];
                ie[2]=i;
            }
		}
    }

	for(i=0;i<3;i++)
	{
		eve[i]=(float3D){a[1][ie[i]],a[2][ie[i]],a[3][ie[i]]};
		eve[i]=sca3D(eve[i],1/norm3D(eve[i]));
	}
	min=max=p[0];
	for(i=0;i<(*m).np;i++)
	{
		p[i]=(float3D){dot3D(p[i],eve[0]),dot3D(p[i],eve[1]),dot3D(p[i],eve[2])};
		if(min.x>p[i].x) min.x=p[i].x;
		if(min.y>p[i].y) min.y=p[i].y;
		if(min.z>p[i].z) min.z=p[i].z;
		if(max.x<p[i].x) max.x=p[i].x;
		if(max.y<p[i].y) max.y=p[i].y;
		if(max.z<p[i].z) max.z=p[i].z;
	}
	printf("bbox:\n\tx=[%f,%f]\n\ty=[%f,%f]\n\tz=[%f,%f]\n",min.x,max.x,min.y,max.y,min.z,max.z);
	printf("centre: {%f,%f,%f}\nmean: {%f,%f,%f}\n",
				cp.x,cp.y,cp.z,
				(min.x+max.x)/2.0,
				(min.y+max.y)/2.0,
				(min.z+max.z)/2.0);
	printf("c-m={%f,%f,%f}\n",
				cp.x-(min.x+max.x)/2.0,
				cp.y-(min.y+max.y)/2.0,
				cp.z-(min.z+max.z)/2.0);

    fmatrix_dispose(a,3+1);				// free memory
    fvector_dispose(w);
    fmatrix_dispose(v,maxp+1);
}
void em_changecoordinates(MeshPtr m,float *mr)
{
	int i;
	float3D *p,x=*(float3D*)mr,y=((float3D*)mr)[1],z=((float3D*)mr)[2];
	
	p=msh_getPointsPtr(m);
	for(i=0;i<(*m).np;i++)
		p[i]=(float3D){dot3D(p[i],x),dot3D(p[i],y),dot3D(p[i],z)};
}
#pragma mark -
#pragma mark ____________________texture operations
/*
void em_pictureToTexture(MeshPtr m)
{
	FSSpec		spec;
	OSErr		err;
	
	err=nav_OpenDialog(&spec);
	
	if(err==noErr)
            open_pict2texture_data(m, spec);
}
void em_openPictureToTextureData(MeshPtr m,FSSpec spec)
{
	PicHandle	pict;
	long		pictSize,hdr=512;
	short		picRefNum;
	char		pictHdr[512];
	
	GWorldPtr	gw;
	GWorldPtr	origPort;
	GDHandle	origDev;
	
	int			i;
	tria2D		t;
	RGBColor	rgb,black={0,0,0},white={0xffff,0xffff,0xffff};
	int		*arr;
	float3D		*p;
	int3D		*T;
	float3D		*C;
        Rect		fwr;
        PicHandle	picrender;
        
        GetWindowPortBounds(FrontWindow(),&fwr);
	
	// open texture picture data
	if(FSpOpenDF(&spec, fsRdPerm, &picRefNum)!=noErr)
		return;

	GetEOF(picRefNum, &pictSize);
	FSRead(picRefNum,&hdr, pictHdr);
	pictSize -= hdr;
	pict = (PicHandle)NewHandle(pictSize);
	HLock((Handle)pict);
	FSRead(picRefNum, &pictSize, *pict);
	HUnlock((Handle)pict);
	FSClose(picRefNum);
	
	// setup a graphic world to use as mask
	GetGWorld(&origPort, &origDev);
	NewGWorld(&gw,32,&fwr,NULL,NULL,NULL);
	SetGWorld(gw, NULL);
	RGBForeColor(&black);
	PaintRect(&fwr);

	// i) render to sort the triangles and ii) draw the picture to project
	SetGWorld(origPort, origDev);
	render(m,1,FrontWindow(),&picrender);
	DrawPicture(pict,&(**pict).picFrame);
	KillPicture(pict);
	
	// a render from near to far that takes the picture as texture
	arr = (int *)*g_em_arrH;
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	C = msh_getTexturePtr(m);
	
	BeginPaint3D();
	for(i=(*m).nt;i>=1;i--)
	{
		t = Transform3D((tria3D){p[T[arr[i]].a],p[T[arr[i]].b],p[T[arr[i]].c]});
		
		// verify that the region is not already used
		SetGWorld(gw, NULL);
		gettrianglemeancolour(t,&rgb);
		
		// if not used, take it as texture
		if(rgb.red+rgb.green+rgb.blue==0)
		{
			SetGWorld(origPort, origDev);
			gettrianglemeancolour(t,&rgb);

			C[T[arr[i]].a].x = rgb.red/(float)0x7fff;
			C[T[arr[i]].a].y = rgb.green/(float)0x7fff;
			C[T[arr[i]].a].z = rgb.blue/(float)0x7fff;
			
			C[T[arr[i]].b].x = rgb.red/(float)0x7fff;
			C[T[arr[i]].b].y = rgb.green/(float)0x7fff;
			C[T[arr[i]].b].z = rgb.blue/(float)0x7fff;
			
			C[T[arr[i]].c].x = rgb.red/(float)0x7fff;
			C[T[arr[i]].c].y = rgb.green/(float)0x7fff;
			C[T[arr[i]].c].z = rgb.blue/(float)0x7fff;

			// mark the region as used
			SetGWorld(gw, NULL);
			RGBForeColor(&white);
			marktriangle(t);
		}
	
		if(Button())
			break;
	}
	EndPaint3D();

	// dispose the gworld
	DisposeGWorld(gw);
}
void em_getTriangleMeanColor(tria2D T, RGBColor *rgb)
{
	Point		min,med,max;
	float		x,y;
	int			i,l0,l1,temp;
	RGBColor	c;
	int			n;
	
	*rgb=(RGBColor){0,0,0};
	n=0;
	
	if(T.a.x<=T.b.x && T.b.x<=T.c.x)	{	min.h=T.a.x;min.v=T.a.y;	med.h=T.b.x;med.v=T.b.y;	max.h=T.c.x;max.v=T.c.y;	}
	if(T.a.x<=T.c.x && T.c.x<=T.b.x)	{	min.h=T.a.x;min.v=T.a.y;	med.h=T.c.x;med.v=T.c.y;	max.h=T.b.x;max.v=T.b.y;	}
	if(T.b.x<=T.a.x && T.a.x<=T.c.x)	{	min.h=T.b.x;min.v=T.b.y;	med.h=T.a.x;med.v=T.a.y;	max.h=T.c.x;max.v=T.c.y;	}
	if(T.b.x<=T.c.x && T.c.x<=T.a.x)	{	min.h=T.b.x;min.v=T.b.y;	med.h=T.c.x;med.v=T.c.y;	max.h=T.a.x;max.v=T.a.y;	}
	if(T.c.x<=T.a.x && T.a.x<=T.b.x)	{	min.h=T.c.x;min.v=T.c.y;	med.h=T.a.x;med.v=T.a.y;	max.h=T.b.x;max.v=T.b.y;	}
	if(T.c.x<=T.b.x && T.b.x<=T.a.x)	{	min.h=T.c.x;min.v=T.c.y;	med.h=T.b.x;med.v=T.b.y;	max.h=T.a.x;max.v=T.a.y;	}
	

	if(	min.h>g_m3d_portRect.left && max.h<g_m3d_portRect.right  &&
		min.v>g_m3d_portRect.top  && min.v<g_m3d_portRect.bottom &&
		med.v>g_m3d_portRect.top  && med.v<g_m3d_portRect.bottom &&
		max.v>g_m3d_portRect.top  && max.v<g_m3d_portRect.bottom)
	{
			for(x=min.h+1; x<=max.h;x++)
			{
				l0 = min.v + (x-min.h)*(max.v-min.v)/(float)(max.h-min.h);
				if(x<med.h)
					l1 = ((med.h-min.h)?(min.v + (x-min.h)*(med.v-min.v)/(float)(med.h-min.h)):(med.v));
				else
					l1 = ((max.h-med.h)?(med.v + (x-med.h)*(max.v-med.v)/(float)(max.h-med.h)):(med.v));
				
				if(l0>l1)
				{	temp=l0;l0=l1;l1=temp;}
				y = l0;					
				i = (l1-l0);
				do
				{
					GetCPixel(x,y,&c);
					y++;
					
					*rgb = (RGBColor){	((*rgb).red*n	+ c.red)/(float)(1+n),
										((*rgb).green*n	+ c.green)/(float)(1+n),
										((*rgb).blue*n	+ c.blue)/(float)(1+n)	};
					n++;
				}
				while(i-->0+1);
			}
	}
}
void em_markTriangle(tria2D T)
{
	Point		min,med,max;
	float		x,y,temp;
	int			i,l0,l1;	
	
	if(T.a.x<=T.b.x && T.b.x<=T.c.x)	{	min.h=T.a.x;min.v=T.a.y;	med.h=T.b.x;med.v=T.b.y;	max.h=T.c.x;max.v=T.c.y;	}
	if(T.a.x<=T.c.x && T.c.x<=T.b.x)	{	min.h=T.a.x;min.v=T.a.y;	med.h=T.c.x;med.v=T.c.y;	max.h=T.b.x;max.v=T.b.y;	}
	if(T.b.x<=T.a.x && T.a.x<=T.c.x)	{	min.h=T.b.x;min.v=T.b.y;	med.h=T.a.x;med.v=T.a.y;	max.h=T.c.x;max.v=T.c.y;	}
	if(T.b.x<=T.c.x && T.c.x<=T.a.x)	{	min.h=T.b.x;min.v=T.b.y;	med.h=T.c.x;med.v=T.c.y;	max.h=T.a.x;max.v=T.a.y;	}
	if(T.c.x<=T.a.x && T.a.x<=T.b.x)	{	min.h=T.c.x;min.v=T.c.y;	med.h=T.a.x;med.v=T.a.y;	max.h=T.b.x;max.v=T.b.y;	}
	if(T.c.x<=T.b.x && T.b.x<=T.a.x)	{	min.h=T.c.x;min.v=T.c.y;	med.h=T.b.x;med.v=T.b.y;	max.h=T.a.x;max.v=T.a.y;	}
	

	if(	min.h>g_m3d_portRect.left && max.h<g_m3d_portRect.right  &&
		min.v>g_m3d_portRect.top  && min.v<g_m3d_portRect.bottom &&
		med.v>g_m3d_portRect.top  && med.v<g_m3d_portRect.bottom &&
		max.v>g_m3d_portRect.top  && max.v<g_m3d_portRect.bottom)
	{
			for(x=min.h+1; x<=max.h;x++)
			{
				l0 = min.v + (x-min.h)*(max.v-min.v)/(float)(max.h-min.h);
				if(x<med.h)
					l1 = ((med.h-min.h)?(min.v + (x-min.h)*(med.v-min.v)/(float)(med.h-min.h)):(med.v));
				else
					l1 = ((max.h-med.h)?(med.v + (x-med.h)*(max.v-med.v)/(float)(max.h-med.h)):(med.v));
				
				if(l0>l1)
				{	temp=l0;l0=l1;l1=temp;}
				y = l0;					
				i = (l1-l0);
				do
				{
					MoveTo(x,y);LineTo(x,y);
					y++;
				}while(i-->0+1);
			}
	}
}
*/
#pragma mark -
bool em_textureDepth(MeshPtr m, float *C)
{
    int			i;
    float		n,max;
    float3D		*p=msh_getPointsPtr(m);
    float       dx,dy,dz;

    dx=((*m).bbox.ide.x-(*m).bbox.siz.x);
    dy=((*m).bbox.ide.y-(*m).bbox.siz.y);
    dz=((*m).bbox.ide.z-(*m).bbox.siz.z);
    max=0;
    for(i=0;i<(*m).np;i++)
    {
        n=	pow(2*(p[i].x-(*m).center.x)/(dx?dx:1),2) +
            pow(2*(p[i].y-(*m).center.y)/(dy?dy:1),2) +
            pow(2*(p[i].z-(*m).center.z)/(dz?dz:1),2);
        C[i] = sqrt(n);
        if(C[i]>max)	max=C[i];
    }
    max*=1.05;	// pure white is not nice...
    for(i=0;i<(*m).np;i++)
        C[i]=C[i]/max;
    return true;
}
bool em_textureMeanCurvature(MeshPtr m, float *C)
{
    float3D		nn,sd,zero={0,0,0};
    float		absmax;
    int			i,j;
    float3D		*p = msh_getPointsPtr(m);
    int3D		*T = msh_getTrianglesPtr(m);
    NTriRec		*NT = msh_getNeighborTrianglesPtr(m);

    for(i=0;i<(*m).np;i++)
    {
        //norm3Dal
        nn = zero;
        for(j=0;j<NT[i].n;j++)
            nn = add3D(nn,tTriPlane(NT[i].t[j],m));
        nn = sca3D(nn,1/norm3D(nn));

        //smoothing direction
        sd = zero;
        for(j=0;j<NT[i].n;j++)
        {
            if(T[NT[i].t[j]].a!=i)	sd = add3D(sd,p[T[NT[i].t[j]].a]);
            if(T[NT[i].t[j]].b!=i)	sd = add3D(sd,p[T[NT[i].t[j]].b]);
            if(T[NT[i].t[j]].c!=i)	sd = add3D(sd,p[T[NT[i].t[j]].c]);
        }
        sd = sca3D(sd,1/(float)(2*NT[i].n));
        sd = sub3D(sd,p[i]);

        C[i] = -dot3D(nn,sd);
    }
    absmax=-1;
    for(i=0;i<(*m).np;i++)
        absmax = (fabs(C[i])>absmax)?fabs(C[i]):absmax;
	//printf("absmax:%f\n",absmax);
	//absmax=0.05;
    
    // there are too big curvatures I
    absmax*=0.25;
    for(i=0;i<(*m).np;i++)
    {
        C[i] = 0.5+ 0.5*C[i]/absmax;
        
        // there are too big curvatures II
        if(C[i]>1)	C[i]=1;
        if(C[i]<0)	C[i]=0;
    }
    
    return true;
}

bool em_textureIntegratedMeanCurvature(MeshPtr m, float *C, int iter)
{
    float3D		nn,sd,zero={0,0,0};
    float		absmax;
    int			i,j,k;
    float3D		*p;
    int3D		*T;
    NTriRec		*NT;

    //----------------
    p = msh_getPointsPtr(m);
    T = msh_getTrianglesPtr(m);
    NT = msh_getNeighborTrianglesPtr(m);
    //----------------

    for(i=0;i<(*m).np;i++)
        C[i]=0;

    printf("integrated mean curvature\n");
    for(k=0;k<iter;k++)
    {
        printf("%i/%i\n",k+1,iter/*g_em_smoothiter*/);
        for(i=0;i<(*m).np;i++)
        {
            //norm3Dal
            nn = zero;
            for(j=0;j<NT[i].n;j++)
                    nn = add3D(nn,tTriPlane(NT[i].t[j],m));
            nn = sca3D(nn,1/norm3D(nn));

            //smoothing direction
            sd = zero;
            for(j=0;j<NT[i].n;j++)
            {
                if(T[NT[i].t[j]].a!=i)	sd = add3D(sd,p[T[NT[i].t[j]].a]);
                if(T[NT[i].t[j]].b!=i)	sd = add3D(sd,p[T[NT[i].t[j]].b]);
                if(T[NT[i].t[j]].c!=i)	sd = add3D(sd,p[T[NT[i].t[j]].c]);
            }
            sd = sca3D(sd,1/(float)(2*NT[i].n));
            sd = sub3D(sd,p[i]);
    
            C[i] += -dot3D(nn,sd)*(1+k/50.0);
        }
        
        em_smooth(m);
        em_smooth(m);
        em_smooth(m);
        em_smooth(m);
        em_smooth(m);
        em_smooth(m);
        em_smooth(m);
        em_smooth(m);
        em_smooth(m);
        em_smooth(m);
    }
    absmax=-1;
    for(i=0;i<(*m).np;i++)
        absmax = (fabs(C[i])>absmax)?fabs(C[i]):absmax;
    for(i=0;i<(*m).np;i++)
        C[i] = 0.5+0.5*C[i]/absmax;
    
    return true;
}
bool em_textureArea(MeshPtr m, float *C)
{
    double		a,max;
    int			i;
    int3D		*T = msh_getTrianglesPtr(m);
    NTriRec		*NT = msh_getNeighborTrianglesPtr(m);

    for(i=0;i<(*m).np;i++) C[i]=0;
	for(i=0;i<(*m).nt;i++)
    {
		a=tTriArea(m,i);
		C[T[i].a]+=a;
		C[T[i].b]+=a;
		C[T[i].c]+=a;
	}
	max=-1;
	for(i=0;i<(*m).np;i++)
	{
		C[i]=C[i]/(double)NT[i].n;
		if(max<C[i])	max=C[i];
	}
	for(i=0;i<(*m).np;i++)
		C[i]=C[i]/max;
    
    return true;
}
bool em_textureGaussCurvature(MeshPtr m)
{
	int			i,j;
	double		sum,angle;
	double		mx,mn;
	int			pstack[SIZESTACK];
	float3D		*p;
	int2D		*E;
	NEdgeRec	*NE;
	float3D		*C;
	char		*M;


	//----------------
	p = msh_getPointsPtr(m);
	E = msh_getEdgesPtr(m);
	NE = msh_getNeighborEdgesPtr(m);
	C = msh_getTexturePtr(m);
	M = msh_getLabelPtr(m);
	//----------------
	
	mx=mn=0;
	for(i=0;i<(*m).np;i++)
	{
            msh_psort(m,i,pstack);
            sum=0;
            for(j=0;j<NE[i].n;j++)
            {
                angle = dot3D(sub3D(p[pstack[j]],p[i]),sub3D(p[pstack[(j+1)%NE[i].n]],p[i]))
                                /norm3D(sub3D(p[pstack[j]],p[i]))
                                /norm3D(sub3D(p[pstack[(j+1)%NE[i].n]],p[i]));
                angle = acos(angle);
                
                sum += angle;
            }
            sum = sum-2*pi;
    
            if(sum>mx)	mx=sum;
            if(sum<mn)	mn=sum;
    
            C[i].z = sum;
            M[i] = 0;
	}

	for(i=0;i<(*m).np;i++)
	{
		C[i].x = (C[i].z>0)?( 20*C[i].z):0;
		C[i].y = (C[i].z<0)?(-20*C[i].z):0;
		C[i].z = 0;
		
		if(C[i].x>2) C[i].x=2;
		if(C[i].y>2) C[i].y=2;
	}
	
	return true;
}
#pragma mark _
void em_textureSulcalHierarchy(MeshPtr m, float *C)
{
	float		sul=0.5;
	int			i,n=0,size,maxsize;
	float		*CC;

	float		a,sulcal,total;
	int3D		*T;
	
	em_textureMeanCurvature(m, C);
	em_textureSmoothVerticesData(m, C);
	
	CC=fvector_new((*m).np);
	for(i=0;i<(*m).np;i++) CC[i]=C[i];

	// search a sulcus
	maxsize=0;
	for(i=0;i<(*m).np;i++)
	{
		if(C[i]>sul)
		{
			size=0;
			fillsulcus(m, C, i, sul,&size);

			if(size>maxsize) maxsize=size;

			if(size>25)
			{
				size=0;
				fillsulcus(m,CC,i,sul,&size);
				n++;
			}
		}
	}
	//printf("%i\n",n);

	// paint
	for(i=0;i<(*m).np;i++)
		if(CC[i]<-0.1)	C[i]=1;
		else			C[i]=0.5;

	// compute sulcal/total area
	T=msh_getTrianglesPtr(m);
	sulcal=total=0;
	for(i=0;i<(*m).nt;i++)
	{
		a=tTriArea(m,i);
		total+=a;
		if(C[T[i].a]==1||C[T[i].b]==1||C[T[i].c]==1)
			sulcal+=a;
	}
	printf("%f %f %f;\n",sulcal,total,sulcal/total);

	fvector_dispose(CC);
}
void em_textureGyralHierarchy(MeshPtr m, float *C)
{
	float		gyr=0.5;
	int			i,n=0,size,maxsize;
	float		*CC;
	
	em_textureMeanCurvature(m, C);
	em_textureSmoothVerticesData(m, C);
	
	CC=fvector_new((*m).np);
	for(i=0;i<(*m).np;i++) CC[i]=C[i];

	// search a sulcus
	maxsize=0;
	for(i=0;i<(*m).np;i++)
	{
		if(C[i]<gyr && C[i]>=0)
		{
			size=0;
			fillgyrus(m, C, i, gyr,&size);
			if(size>maxsize) maxsize=size;
			
			if(size>25)
			{
				size=0;
				fillgyrus(m,CC,i,gyr,&size);
				n++;
			}
		}
	}
	printf("%i\n",n);

	// paint
	for(i=0;i<(*m).np;i++)
		if(CC[i]<0) C[i]=0;
		else		C[i]=0.5;
}
void em_textureCountSulci(MeshPtr m, float *C)
{
	float		sul=0.5;
	int			i,n=0,size,maxsize;
	
	// search a sulcus
	maxsize=0;
	for(i=0;i<(*m).np;i++)
	{
		if(C[i]>sul)
		{
			size=0;
			fillsulcus(m, C, i, sul,&size);
			if(size>maxsize) maxsize=size;
			n++;
		}
	}
	printf("%i (%i)\n",n,maxsize);

	// paint
	for(i=0;i<(*m).np;i++)
		if(C[i]<0)  C[i]=1;
		else		C[i]=0.5;
}
void em_textureCountGyri(MeshPtr m, float *C)
{
	float		gyr=0.01;
	int			i,n=0,size,maxsize;
	
	// search a sulcus
	maxsize=0;
	for(i=0;i<(*m).np;i++)
	{
		if(C[i]<gyr && C[i]>=0)
		{
			size=0;
			fillgyrus(m, C, i, gyr,&size);
			if(size>maxsize) maxsize=size;
			n++;
		}
	}
	printf("%i\n",n);

	// paint
	for(i=0;i<(*m).np;i++)
		if(C[i]<0)  C[i]=0;
		else		C[i]=0.5;
}
void fillsulcus(MeshPtr m, float *C, int i, float sul, int *size)
{
	NEdgeRec	*NE;
	int2D		*E;
	int			ii,j;
	
	C[i]*=-1;
	(*size)++;

	NE=msh_getNeighborEdgesPtr(m);
	E = msh_getEdgesPtr(m);
	
	for(j=0;j<NE[i].n;j++)
	{
		if(E[NE[i].e[j]].a==i)  ii=E[NE[i].e[j]].b;
		else					ii=E[NE[i].e[j]].a;
		
		if(C[ii]>sul)   fillsulcus(m,C,ii,sul,size);
	}
}
void fillgyrus(MeshPtr m, float *C, int i, float gyr, int *size)
{
	NEdgeRec	*NE;
	int2D		*E;
	int			ii,j;
	
	C[i]-=10;
	(*size)++;

	NE=msh_getNeighborEdgesPtr(m);
	E = msh_getEdgesPtr(m);
	
	for(j=0;j<NE[i].n;j++)
	{
		if(E[NE[i].e[j]].a==i)  ii=E[NE[i].e[j]].b;
		else					ii=E[NE[i].e[j]].a;
		
		if(C[ii]<gyr&&C[ii]>=0)   fillgyrus(m,C,ii,gyr,size);
	}
}
#pragma mark _
void em_textureSulcalDepth2(MeshPtr m, MeshPtr h, float *C)
{
	int i,j,result;
	float   d,max=0;
	float2D Coord;
	float3D V,t[3],p;
	int3D   *T;
	
	T=msh_getTrianglesPtr(h);
	for(i=0;i<(*m).np;i++)
	{
		if(i%1000==0) printf("%i\n",i);
		d=-1;
		for(j=0;j<(*h).nt;j++)
		{
			V=sca3D(sub3D((*m).p[i],(*m).center),10);
			t[0]=sub3D((*h).p[T[j].a],(*m).center);
			t[1]=sub3D((*h).p[T[j].b],(*m).center);
			t[2]=sub3D((*h).p[T[j].c],(*m).center);
			result=intersect_VectorTriangle(V, t, &Coord, 10e-4);
			if(result==1)
			{
				p=add3D(t[0],add3D(sca3D(sub3D(t[1],t[0]),Coord.x),sca3D(sub3D(t[2],t[0]),Coord.y)));
				d=norm3D(sub3D(p,(*m).p[i]));
				C[i]=d;
				if(d>max) max=d;
				//printf("%f\n",d);
				break;
			}
		}
		if(d<0) printf(">>>>>>>>>>>>>>>>> >>>>>>>>>>>>> oops.\n");
	}
	printf("%f\n",max);
	for(i=0;i<(*m).np;i++)
		C[i]=1-C[i]/max;
}
#pragma mark _
void em_textureLaplaceFD(MeshPtr m)
{
	float3D			*p;
	int3D			*T;
	int2D			*E;
	NTriRec			*NT;
	NEdgeRec		*NE;
	float3D			*C;
	char			*M;
	int				i,j,indx;
	int				ixmax,ixmin;
	double			area;
	double			*D,*N,*R,*xR;
	ContainerRec	Dcont,Ncont;
	int				iter=0;
	double			err=0;
	int				pstack[SIZESTACK],estack[SIZESTACK];
	int				a,b;
	float			px,nx;

	if(m==NULL)		return;
	
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT= msh_getNeighborTrianglesPtr(m);
	E = msh_getEdgesPtr(m);
	NE= msh_getNeighborEdgesPtr(m);
	C = msh_getTexturePtr(m);
	M = msh_getLabelPtr(m);
	
	//0. find a pole (the frontal pole)
	ixmax=0;
	ixmin=0;
	for(i=0;i<(*m).np;i++)
	{
		if(p[ixmax].x<p[i].x)	ixmax=i;
		if(p[ixmin].x>p[i].x)	ixmin=i;
	}
	M[ixmax]=1;
	M[ixmin]=1;
	
	//-----------------------
	// Laplace equation uxx=0
	// (finite differencing)
	//-----------------------

	// 1.1. border conditions
	R=dvector_new((*m).np);
	xR=dvector_new((*m).np);

	R[ixmax]= 1;
	R[ixmin]=-1;

	// 1.2. coefficient matrix
	msh_setContainer(&Dcont,(*m).np,sizeof(double),0,"\pdiag");
	D=(double*)msh_getDataPtr(m,Dcont);
	msh_setContainer(&Ncont,(*m).nt*3/2+1,sizeof(double),0,"\pneig");
	N=(double*)msh_getDataPtr(m,Ncont);	

	for(i=0;i<(*m).np;i++)
	{
		msh_esort(m,i,pstack,estack);
		if(i!=ixmax && i!=ixmin)
		{
			D[i]=0;
			for(j=0;j<NE[i].n;j++)
			{
				a=(j+1)%NE[i].n;
				b=(j+2)%NE[i].n;
				
				if(E[estack[a]].a!=i)	indx=E[estack[a]].a;
				else					indx=E[estack[a]].b;
				
				area=	triArea(p[pstack[j]],p[pstack[a]],p[i]) +
                                        triArea(p[pstack[a]],p[pstack[b]],p[i]);
				
				area=area/500.0;
				
				D[i]+=area;
				// dirichlet boundary
				if(indx!=ixmax && indx!=ixmin)		
					N[estack[a]]=-area;
				else
					R[i]+=R[indx]*area;
			}
		}
		else
			// dirichlet boundary
			D[i]=1;
	}

	// 1.3. solve: conjugate gradients in the mesh sparse matrix
	msh_linbcgDN((*m).np,(*m).nt*3/2, D,N,E,R, xR, 2, 1e-10, 10*(*m).np, &iter, &err);
	
	// 2.1. calculate positive v/s negative area
	px=nx=0;
	for(i=0;i<(*m).np;i++)
	{
		area=0;
		for(j=0;j<NT[i].n;j++)
			area+=tTriArea(m,NT[i].t[j]);
		
		if(xR[i]>0)	px+=area;
		if(xR[i]<0)	nx+=area;
	}
        printf("x: +a/-a = %f / %f", px, nx);

	/*
        for(i=0;i<(*m).np;i++)
	{
		ForeColor(fabs(xR[i])>1?redColor:blackColor);
		MoveTo(10+acos(dot3D(sca3D(sub3D(p[i],(*m).center),1/norm3D(sub3D(p[i],(*m).center))),(float3D){1,0,0}))*100,300);
		LineTo(10+acos(dot3D(sca3D(sub3D(p[i],(*m).center),1/norm3D(sub3D(p[i],(*m).center))),(float3D){1,0,0}))*100,300-xR[i]*100);
		
		ForeColor(greenColor);
		MoveTo(10+i*100*pi/(float)(*m).np,300-100*cos(i*pi/(float)(*m).np));
		LineTo(10+i*100*pi/(float)(*m).np,300-100*cos(i*pi/(float)(*m).np));
	}
        */

	// 2.2. set texture as coordinates
	for(i=0;i<(*m).np;i++)
	{
		C[i].x=2*xR[i];
		C[i].y=0;
		C[i].z=0;
	}
	
	msh_deleteData(m,Dcont);
	msh_deleteData(m,Ncont);
	dvector_dispose(R);
	dvector_dispose(xR);
}
void em_textureLaplaceFE(MeshPtr m)
{
	float3D		*p;
	NTriRec		*NT;
	int2D		*E;
	float3D		*C;
	int		i,j;
	int		ixmax,ixmin;
	int		nDir,iDir[2];
	double		Dir[2];
	double		*X;
	double		sum,px,nx;
	
	if(m==NULL)	return;
	
	p = msh_getPointsPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	E = msh_getEdgesPtr(m);
	C = msh_getTexturePtr(m);
	
	X=dvector_new((*m).np);

	// Solve Laplace equation with the X poles as Dirichlet conditions
	ixmax=ixmin=0;
	for(i=0;i<(*m).np;i++)	{	if(p[ixmax].x<p[i].x)	ixmax=i;
								if(p[ixmin].x>p[i].x)	ixmin=i;
							}
	nDir=2;
	iDir[0]=ixmax;	iDir[1]=ixmin;
	Dir[0]= 1;		Dir[1]=-1;
	msh_laplace_fe(m,X, nDir,iDir,Dir, 0,0,0, 0,0,0,0);

	// Calculate positive v/s negative area
	px=nx=0;
	for(i=0;i<(*m).np;i++)
	{
		sum=0;
		for(j=0;j<NT[i].n;j++)	sum+=tTriArea(m,NT[i].t[j]);
		
		if(X[i]>0)	px+=sum/3.0;
		if(X[i]<0)	nx+=sum/3.0;
	}
        printf("+a/-a = %f / %f", px,nx);
        printf("+X/-X = %f / %f", X[ixmax], X[ixmin]);

        /*
	SH=fvector_new(520);
	max=1;
	for(i=0;i<(*m).np;i++)
	{
		theta = i*pi/(float)(*m).np;
		ForeColor(fabs(X[i])>1?redColor:blueColor);
		MoveTo(10+acos(dot3D(sca3D(sub3D(p[i],(*m).center),1/norm3D(sub3D(p[i],(*m).center))),(float3D){1,0,0}))*100,300);
		LineTo(10+acos(dot3D(sca3D(sub3D(p[i],(*m).center),1/norm3D(sub3D(p[i],(*m).center))),(float3D){1,0,0}))*100,300-X[i]*100);
		
		ForeColor(greenColor);
		MoveTo(10+100*theta,300-100*(2.5*pow(cos(theta),3)-1.5*cos(theta)+cos(theta)));
		LineTo(10+100*theta,300-100*(2.5*pow(cos(theta),3)-1.5*cos(theta)+cos(theta)));

		ForeColor(magentaColor);
		ee=0.875;
		MoveTo(10+100*theta,300+100*((sqrt(pow(sin(theta),ee)+pow(1-cos(theta),ee))-sqrt(pow(sin(theta),ee)+pow(1+cos(theta),ee)))/(sqrt(pow(sin(theta),ee)+pow(1-cos(theta),ee))+sqrt(pow(sin(theta),ee)+pow(1+cos(theta),ee)))));
		LineTo(10+100*theta,300+100*((sqrt(pow(sin(theta),ee)+pow(1-cos(theta),ee))-sqrt(pow(sin(theta),ee)+pow(1+cos(theta),ee)))/(sqrt(pow(sin(theta),ee)+pow(1-cos(theta),ee))+sqrt(pow(sin(theta),ee)+pow(1+cos(theta),ee)))));
		if(fabs(X[i])>max)	max=fabs(X[i]);
		
		SH[(int)(512*acos(dot3D(sca3D(sub3D(p[i],(*m).center),1/norm3D(sub3D(p[i],(*m).center))),(float3D){1,0,0}))/pi)]=X[i];
	}
	//fvector_save(513,SH);
	fvector_dispose(SH);
	pscopy(s,"\pmax |X| = ");fnum_to_string(max,4,s1);psadd(s,s1);
	printstring(s);
	*/

	// Set Laplace equation solution as texture
	for(i=0;i<(*m).np;i++)
		C[i]=(float3D){8*em_angleFromLaplace(X[i])/(2*pi),0,0};

	dvector_dispose(X);
}

void em_textureCoordinatesFE(MeshPtr m)
{
	float3D		*p;
	NTriRec		*NT;
	int2D		*E;
	float3D		*C;
	char		*M;
	int		i,j;
	int		ixmax,ixmin,iymax,iymin,izmax,izmin;
	
	int		nDir,iDir[2];
	double		Dir[2];
	
	int		nMFC,ms,sl;
	int2D		*eMFC;
	double		*gMFC,*hMFC,t;
	
	
	double		*X,*Y,*Z;
	int		*U;

	double		area,sum,X0,Y0,Z0;
	double		px,nx,py,ny,pz,nz;
        
	double		closeto0;

	if(m==NULL)	return;
	
	p = msh_getPointsPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	E = msh_getEdgesPtr(m);
	C = msh_getTexturePtr(m);
	M = msh_getLabelPtr(m);
	
	X=dvector_new((*m).np);
	Y=dvector_new((*m).np);
	Z=dvector_new((*m).np);

	U=ivector_new((*m).np);

	// 0. calculate the total area
	area=0;
	for(i=0;i<(*m).np;i++)
	{
		for(j=0;j<NT[i].n;j++)
			area+=tTriArea(m,NT[i].t[j])/3.0;
	}
	
	//__________________________________________________________________________
	//1. Great circle X
	//__________________________________________________________________________
	//	1.1 find preliminary X poles
		printf("computing X great circle\n");
		ixmax=ixmin=0;
		for(i=0;i<(*m).np;i++)	{	if(p[ixmax].x<p[i].x)	ixmax=i;
                                                if(p[ixmin].x>p[i].x)	ixmin=i;
                                        }
	//	1.2 solve Laplace equation
		nDir=2;
		iDir[0]=ixmax;	iDir[1]=ixmin;
		Dir[0]= 1;		Dir[1]=-1;
		msh_laplace_fe(m,X, nDir,iDir,Dir, 0,0,0, 0,0,0,0);
	//	1.3. Find the great circle equidistant to the poles xmax,xmin
		for(i=0;i<(*m).np;i++)	U[i]=i;
		dquicksort((*m).np,U-1,X);
		sum=0;
		for(i=0;i<(*m).np;i++)
		{	for(j=0;j<NT[U[i]].n;j++) sum+=tTriArea(m,NT[U[i]].t[j])/3.0;
			if(sum>area/2.0){	X0=X[U[i]];	break;	}
		}

	//__________________________________________________________________________
	//2. Great circle Z
	//__________________________________________________________________________
		printf("computing Z great circle\n");
	//	2.1 find preliminary Z poles
                closeto0=((*m).bbox.ide.y-(*m).bbox.siz.y)*0.1; // 10% of the total height
		izmax=izmin=-1;
		for(i=0;i<(*m).nt*3/2;i++)
                    if((X[E[i].a]-X0)*(X[E[i].b]-X0)<=0 && (fabs(p[E[i].a].y-(*m).center.y)<closeto0))
                    {
                        if(izmax==-1)	izmax=izmin=E[i].a;
                        
                        if(p[E[i].a].z>p[izmax].z)	izmax=E[i].a;
                        if(p[E[i].b].z>p[izmax].z)	izmax=E[i].b;
                        
                        if(p[E[i].a].z<p[izmin].z)	izmin=E[i].a;
                        if(p[E[i].b].z<p[izmin].z)	izmin=E[i].b;
                    }
	//	2.2 solve Laplace equation
		nDir=2;
		iDir[0]=izmax;	iDir[1]=izmin;
		Dir[0]= 1;		Dir[1]=-1;
		msh_laplace_fe(m,Z, nDir,iDir,Dir, 0,0,0, 0,0,0,0);
	//	2.3. Find the great circle equidistant to the poles zmax,zmin
		for(i=0;i<(*m).np;i++)	U[i]=i;
		dquicksort((*m).np,U-1,Z);
		sum=0;
		for(i=0;i<(*m).np;i++)
		{	for(j=0;j<NT[U[i]].n;j++)	sum+=tTriArea(m,NT[U[i]].t[j])/3.0;
			if(sum>area/2.0){	Z0=Z[U[i]];	break;	}
		}

	//__________________________________________________________________________
	//3. Great circle Y
	//__________________________________________________________________________
		printf("computing Y great circle\n");
	//	3.1 Find antipodal pair Y
		iymax=iymin=-1;
		for(i=0;i<(*m).nt*3/2;i++)
		{
			if( (X[E[i].a]-X0)*(X[E[i].b]-X0)<=0 &&
				(Z[E[i].a]-Z0)*(Z[E[i].b]-Z0)<=0 )
			{
				if(iymax==-1)	iymax=iymin=E[i].a;
				
				if(p[E[i].a].y>p[iymax].y)	iymax=E[i].a;
				if(p[E[i].b].y>p[iymax].y)	iymax=E[i].b;
				
				if(p[E[i].a].y<p[iymin].y)	iymin=E[i].a;
				if(p[E[i].b].y<p[iymin].y)	iymin=E[i].b;
			}
		}
	//	3.2 solve Laplace equation
		nDir=2;
		iDir[0]=iymax;	iDir[1]=iymin;
		Dir[0]= 1;		Dir[1]=-1;
		msh_laplace_fe(m,Y, nDir,iDir,Dir, 0,0,0, 0,0,0,0);
	//	3.3. Find the great circle equidistant to the poles zmax,zmin
		for(i=0;i<(*m).np;i++)	U[i]=i;
		dquicksort((*m).np,U-1,Y);
		sum=0;
		for(i=0;i<(*m).np;i++)
		{	for(j=0;j<NT[U[i]].n;j++)	sum+=tTriArea(m,NT[U[i]].t[j])/3.0;
			if(sum>area/2.0){	Y0=Y[U[i]];	break;	}
		}

	//__________________________________________________________________________
	// 4. Find X,Z antipodal pairs
	//__________________________________________________________________________
		ixmax=ixmin=izmax=izmin=-1;
		for(i=0;i<(*m).nt*3/2;i++)
		{
			if( (Y[E[i].a]-Y0)*(Y[E[i].b]-Y0)<=0 &&
				(Z[E[i].a]-Z0)*(Z[E[i].b]-Z0)<=0 )
			{
				if(ixmax==-1)	ixmax=ixmin=E[i].a;
				
				if(p[E[i].a].x>p[ixmax].x)	ixmax=E[i].a;
				if(p[E[i].b].x>p[ixmax].x)	ixmax=E[i].b;
				
				if(p[E[i].a].x<p[ixmin].x)	ixmin=E[i].a;
				if(p[E[i].b].x<p[ixmin].x)	ixmin=E[i].b;
			}
			if( (X[E[i].a]-X0)*(X[E[i].b]-X0)<=0 &&
				(Y[E[i].a]-Y0)*(Y[E[i].b]-Y0)<=0 )
			{
				if(izmax==-1)	izmax=izmin=E[i].a;
				
				if(p[E[i].a].z>p[izmax].z)	izmax=E[i].a;
				if(p[E[i].b].z>p[izmax].z)	izmax=E[i].b;
				
				if(p[E[i].a].z<p[izmin].z)	izmin=E[i].a;
				if(p[E[i].b].z<p[izmin].z)	izmin=E[i].b;
			}
		}
	//__________________________________________________________________________
	// 5. Set the Laplace coordinate system with
	//		Dirichlet boundaries at the X,Y,Z antipodal pairs
	//		and Multifreedom constrains at the X,Y,Z great circles
	//__________________________________________________________________________
	//	5.1 X coordinates
		printf("computing X coordinates\n");
		// multifreedom constraints
		nMFC=0;
		for(i=0;i<(*m).nt*3/2;i++)
			if((X[E[i].a]-X0)*(X[E[i].b]-X0)<=0)	nMFC++;
		eMFC=i2vector_new(nMFC);
		gMFC=dvector_new(nMFC);
		hMFC=dvector_new(nMFC);
		nMFC=0;
		for(i=0;i<(*m).nt*3/2;i++)
			if((X[E[i].a]-X0)*(X[E[i].b]-X0)<=0)
			{
				// slave is the point with coordinate closer to the intersection
				// to avoid an eventual zero division
				// t is the parametric intersection point from slave to master
				// (if t=0, this turns to be a Dirichlet condition)
				// the edge is stored slave (t=0) first
				if(fabs(X[E[i].a]-X0)>fabs(X[E[i].b]-X0)){	ms=E[i].a;	sl=E[i].b;}
				else									 {	ms=E[i].b;	sl=E[i].a;}
				t=(X[sl]-X0)/(X[sl]-X[ms]);
				
				eMFC[nMFC]=(int2D){sl,ms};
				gMFC[nMFC]=t/(t-1);
				nMFC++;
			}
		// dirichlet boundary
		nDir=0;
		iDir[nDir]=ixmax;	Dir[nDir]= 1;	nDir++;
		iDir[nDir]=ixmin;	Dir[nDir]=-1;	nDir++;
		// solve Laplace
		msh_laplace_fe(m,X,nDir,iDir,Dir,0,0,0,nMFC,eMFC,gMFC,hMFC);
		// dispose
		i2vector_dispose(eMFC);
		dvector_dispose(gMFC);
		dvector_dispose(hMFC);

	//	5.2 Y coordinates
		printf("computing Y coordinates\n");
		// multifreedom constraints
		nMFC=0;
		for(i=0;i<(*m).nt*3/2;i++)
			if((Y[E[i].a]-Y0)*(Y[E[i].b]-Y0)<=0)	nMFC++;
		eMFC=i2vector_new(nMFC);
		gMFC=dvector_new(nMFC);
		hMFC=dvector_new(nMFC);
		nMFC=0;
		for(i=0;i<(*m).nt*3/2;i++)
			if((Y[E[i].a]-Y0)*(Y[E[i].b]-Y0)<=0)
			{
				if(fabs(Y[E[i].a]-Y0)>fabs(Y[E[i].b]-Y0)){	ms=E[i].a;	sl=E[i].b;}
				else									 {	ms=E[i].b;	sl=E[i].a;}
				t=(Y[sl]-Y0)/(Y[sl]-Y[ms]);
				
				eMFC[nMFC]=(int2D){sl,ms};
				gMFC[nMFC]=t/(t-1);
				nMFC++;
			}
		// dirichlet boundary
		nDir=0;
		iDir[nDir]=iymax;	Dir[nDir]= 1;	nDir++;
		iDir[nDir]=iymin;	Dir[nDir]=-1;	nDir++;
		// solve Laplace
		msh_laplace_fe(m,Y,nDir,iDir,Dir,0,0,0,nMFC,eMFC,gMFC,hMFC);
		// dispose
		i2vector_dispose(eMFC);
		dvector_dispose(gMFC);
		dvector_dispose(hMFC);
		
	//	5.3 Z coordinates
		printf("computing Z coordinates\n");
		// multifreedom constraints
		nMFC=0;
		for(i=0;i<(*m).nt*3/2;i++)
			if((Z[E[i].a]-Z0)*(Z[E[i].b]-Z0)<=0)	nMFC++;
		eMFC=i2vector_new(nMFC);
		gMFC=dvector_new(nMFC);
		hMFC=dvector_new(nMFC);
		nMFC=0;
		for(i=0;i<(*m).nt*3/2;i++)
			if((Z[E[i].a]-Z0)*(Z[E[i].b]-Z0)<=0)
			{
				if(fabs(Z[E[i].a]-Z0)>fabs(Z[E[i].b]-Z0)){	ms=E[i].a;	sl=E[i].b;}
				else									 {	ms=E[i].b;	sl=E[i].a;}
				t=(Z[sl]-Z0)/(Z[sl]-Z[ms]);
				
				eMFC[nMFC]=(int2D){sl,ms};
				gMFC[nMFC]=t/(t-1);
				nMFC++;
			}
		// dirichlet boundary
		nDir=0;
		iDir[nDir]=izmax;	Dir[nDir]= 1;	nDir++;
		iDir[nDir]=izmin;	Dir[nDir]=-1;	nDir++;	
		// solve Laplace
		msh_laplace_fe(m,Z,nDir,iDir,Dir,0,0,0,nMFC,eMFC,gMFC,hMFC);
		// dispose
		i2vector_dispose(eMFC);
		dvector_dispose(gMFC);
		dvector_dispose(hMFC);
		
		X0=Y0=Z0=0;

	//__________________________________________________________________________
	// 6. Display information
	//__________________________________________________________________________
	//	6.1 draw the X,Z,Y antipodal pairs
		M[ixmax]=M[ixmin]=M[iymax]=M[iymin]=M[izmax]=M[izmin]=6;

	//	6.2 calculate positive v/s negative area
		px=nx=py=ny=pz=nz=0;
		for(i=0;i<(*m).np;i++)
		{
			sum=0;
			for(j=0;j<NT[i].n;j++)	sum+=tTriArea(m,NT[i].t[j]);
			
			if(X[i]-X0>0)	px+=sum/3.0;
			if(X[i]-X0<0)	nx+=sum/3.0;

			if(Y[i]-Y0>0)	py+=sum/3.0;
			if(Y[i]-Y0<0)	ny+=sum/3.0;

			if(Z[i]-Z0>0)	pz+=sum/3.0;
			if(Z[i]-Z0<0)	nz+=sum/3.0;
		}
                printf("x: +a/-a = %f / %f\n", px, nx);
                printf("y: +a/-a = %f / %f\n", py, ny);
                printf("z: +a/-a = %f / %f\n", pz, nz);

	//	6.3. Set texture as coordinates
		for(i=0;i<(*m).np;i++)
		{
			C[i].x=2*X[i];
			C[i].y=2*Y[i];
			C[i].z=2*Z[i];
/*			C[i].x=48*angle_from_laplace(X[i])/(2*pi);//2*X[i];
			C[i].y=48*angle_from_laplace(Y[i])/(2*pi);//2*Y[i];
			C[i].z=48*angle_from_laplace(Z[i])/(2*pi);//2*Z[i];*/
		}
	
	//__________________________________________________________________________
	//7. Dispose
	//__________________________________________________________________________
	dvector_dispose(X);
	dvector_dispose(Y);
	dvector_dispose(Z);

	ivector_dispose(U);
}
/*
void em_textureHeat(MeshPtr m)
{
	ContainerRec	c;
	float3D		*p,*C;
	int3D		*T;
	NTriRec		*NT;
	double3D		*DC;
	char		*M;
	int			ii,i,j,k,l,indx;
	int			ixmax,ixmin, iymax,iymin;
	double		area,c0,c1;
	//Rect		r={0,0,30,300};
	double		tau=0.01;
	double		gxsum,gXsum;
	double		gysum,gYsum;

	if(m==NULL)		return;
	
	msh_setContainer(&c,(*m).np,sizeof(double3D),0,"\ptemp");
	DC = (double3D*)msh_getDataPtr(m,c);

	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	C = msh_getTexturePtr(m);
	M = msh_getLabelPtr(m);
	
	// 1. compute point cell area and find max,min mesh points in X,Y
	ixmax=ixmin=iymax=iymin=0;
	for(i=0;i<(*m).np;i++)
	{
		area=0;
		for(j=0;j<NT[i].n;j++)
			area+=triArea(p[T[NT[i].t[j]].a],p[T[NT[i].t[j]].b],p[T[NT[i].t[j]].c]);

		DC[i].x=0;
		DC[i].y=0;
		DC[i].z=area;
		
		#define DIST(X,V) norm3D(sub3D(DIRECTION(sub3D(X,g_m3d_center)),V))
		if(DIST(p[ixmax],((float3D){ 1,0,0})) <DIST(p[i],((float3D){ 1,0,0})))	{ixmax=i;}
		if(DIST(p[ixmin],((float3D){-1,0,0})) <DIST(p[i],((float3D){-1,0,0})))	{ixmin=i;}
		
		if(DIST(p[iymax],((float3D){0, 1,0})) <DIST(p[i],((float3D){0, 1,0})))	{iymax=i;}
		if(DIST(p[iymin],((float3D){0,-1,0})) <DIST(p[i],((float3D){0,-1,0})))	{iymin=i;}
		#undef DIST
	}
		
	//----------------------
	// heat equation ut=uxx
	// monte carlo finite difference iterations
	//----------------------
	
	l=0;k=0;
	do
	{
		MoveTo(l*3,k);LineTo(l*3+1,k);k++;
		if(k>300)// FrontWindow()->portRect.bottom)
		{
			k=0;l++;
			gxsum=gXsum=gysum=gYsum=0;
			for(i=0;i<(*m).np;i++)
			{
				if(DC[i].x>0)	gXsum+= DC[i].z*DC[i].x;
				if(DC[i].x<0)	gxsum+=-DC[i].z*DC[i].x;
				
				if(DC[i].y>0)	gYsum+= DC[i].z*DC[i].y;
				if(DC[i].y<0)	gysum+=-DC[i].z*DC[i].y;

				C[i].x=DC[i].x;
				C[i].y=DC[i].y;
				C[i].z=0;
				
			}
			MoveTo(2,12);
			EraseRect(&r);
			NumToString(gxsum,s);DrawString("\p·x=");DrawString(s);
			NumToString(gXsum,s);DrawString("\p ·X=");DrawString(s);
			NumToString(gysum,s);DrawString("\p·y=");DrawString(s);
			NumToString(gYsum,s);DrawString("\p ·Y=");DrawString(s);
			render(m,1,FrontWindow(),&pict);
		}

	
		// x,y diffusion
		for(ii=0;ii<(*m).np;ii++)
		{
			i=(*m).np*(0.5+Random()/65536.0);
			j=NT[i].n*(0.5+Random()/65536.0);
			if(T[NT[i].t[j]].a!=i)	indx=T[NT[i].t[j]].a;
			if(T[NT[i].t[j]].b!=i)	indx=T[NT[i].t[j]].b;
			if(T[NT[i].t[j]].c!=i)	indx=T[NT[i].t[j]].c;
			
			DC[ixmax].x+=-0.001/DC[ixmax].z;//	DC[ixmax].y=0;
			DC[ixmin].x+= 0.001/DC[ixmin].z;//	DC[ixmin].y=0;

			DC[iymax].y+=-0.001/DC[iymax].z;//	DC[iymax].x=0;
			DC[iymin].y+= 0.001/DC[iymin].z;//	DC[iymin].x=0;

			c0=DC[i].x;	c1=DC[indx].x;
			DC[i].x 	+= tau*(c1-c0)*DC[indx].z/(DC[indx].z+DC[i].z);
			DC[indx].x 	+= tau*(c0-c1)*DC[i].z/(DC[indx].z+DC[i].z);

			c0=DC[i].y;	c1=DC[indx].y;
			DC[i].y 	+= tau*(c1-c0)*DC[indx].z/(DC[indx].z+DC[i].z);
			DC[indx].y 	+= tau*(c0-c1)*DC[i].z/(DC[indx].z+DC[i].z);

		}
	}
	while(!Button());
	for(i=0;i<(*m).np;i++)
	{	C[i].x = DC[i].x;C[i].y=DC[i].y;C[i].z=0;}
	
	msh_deleteData(m,c);
}
*/
/*
void em_textureConformalMapping(MeshPtr m)
{
	ContainerRec	c;
	float3D		*p,zero={0,0,0};
	int3D		*T;
	NTriRec		*NT;
	float3D		*C,*C1;
	char		*M;
	int			i,j,k,l;
	int			ixmax,ixmin, iymax,iymin;
	float		area;
	float3D		du0,dv0,du1,dv1;
	float3D		u0,v0,u1,v1;
	float3D		gt[3],gb[3];
	float		mx[9],b0,b1;
	
	if(m==NULL)		return;

	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	C = msh_getTexturePtr(m);
	M = msh_getLabelPtr(m);
	
	msh_setContainer(&c,(*m).np,sizeof(float3D),0,"\ptemp");
	C1 = (float3D*)msh_getDataPtr(m,c);
	
	// 1. compute point cell area and find max,min mesh points in X,Y
	ixmax=ixmin=iymax=iymin=0;
	for(i=0;i<(*m).np;i++)
	{
		area=0;
		for(j=0;j<NT[i].n;j++)
			area+=tTriArea(m,NT[i].t[j]);

		C1[i]=zero;
		C[i].z=area;
		
		#define DIST(X,V) norm3D(sub3D(DIRECTION(sub3D(X,g_m3d_center)),V))
		if(DIST(p[ixmax],((float3D){ 1,0,0})) <DIST(p[i],((float3D){ 1,0,0})))	{ixmax=i;}
		if(DIST(p[ixmin],((float3D){-1,0,0})) <DIST(p[i],((float3D){-1,0,0})))	{ixmin=i;}
		
		if(DIST(p[iymax],((float3D){0, 1,0})) <DIST(p[i],((float3D){0, 1,0})))	{iymax=i;}
		if(DIST(p[iymin],((float3D){0,-1,0})) <DIST(p[i],((float3D){0,-1,0})))	{iymin=i;}
		#undef DIST
	}
		
	// 2. diffusion	
	l=0;k=0;
	do
	{
		MoveTo(l*3,k);LineTo(l*3+1,k);k++;
		if(k>20) //FrontWindow()->portRect.bottom)
		{
			k=0;l++;
			render(m,1,FrontWindow(),&pict);
		}

		// x,y diffusion
		//for(i=0;i<(**mesh).np;i++)
		//{
		//	C[ixmax].x=-1;C[ixmax].y= 0;	C[ixmin].x=1;C[ixmin].y= 0;
		//	C[iymax].y= 0;C[iymax].y=-1;	C[iymin].x=0;C[iymin].y=+1;
                //
		//	xsum=ysum=0;	area=0;
		//	for(j=0;j<NT[i].n;j++)
		//	{
		//		if(T[NT[i].t[j]].a!=i)	indx=T[NT[i].t[j]].a;
		//		if(T[NT[i].t[j]].b!=i)	indx=T[NT[i].t[j]].b;
		//		if(T[NT[i].t[j]].c!=i)	indx=T[NT[i].t[j]].c;
				
		//		xsum += C[indx].x*C[indx].z;
		//		ysum += C[indx].y*C[indx].z;
		//		area += C[indx].z;
		//	}
		//	C1[i].x = xsum/area;
		//	C1[i].y = ysum/area;
		//}
		//for(i=0;i<(**mesh).np;i++)
		//{
		//	C[i].x = C1[i].x;
		//	C[i].y = C1[i].y;
		//}
		
		// conformality
		for(i=0;i<(*m).np;i++)
			C1[i]=zero;
		for(i=0;i<(*m).nt;i++)
		{
			// local base from the triangle
			gt[0]=p[T[i].a];	gt[1]=p[T[i].b];	gt[2]=p[T[i].c];
			triBase(gt,gb);

			mx[0]=0;								mx[1]=0;								mx[2]=1;
			mx[3]=dot3D(sub3D(gt[1],gt[0]),gb[0]);	mx[4]=dot3D(sub3D(gt[1],gt[0]),gb[1]);	mx[5]=1;
			mx[6]=dot3D(sub3D(gt[2],gt[0]),gb[0]);	mx[7]=dot3D(sub3D(gt[2],gt[0]),gb[1]);	mx[8]=1;
			
			// u coordinates local base
			u0=(float3D){C[T[i].a].x,C[T[i].b].x,C[T[i].c].x};
			du0=mth_matsolvefloat3D(mx,u0);
			// v coordinates local base
			v0=(float3D){C[T[i].a].y,C[T[i].b].y,C[T[i].c].y};
			dv0= mth_matsolvefloat3D(mx,v0);
			
			//affine coordinate vectors
			du1=sca3D((float3D){du0.x+dv0.y,du0.y-dv0.x,0},0.5);
			dv1=sca3D((float3D){dv0.x-du0.y,du0.x+dv0.y,0},0.5);
			
			b0=(u0.x+u0.y+u0.z)/3.0;
			u1 = mth_matfloat3D(mx,(float3D){du1.x,du1.y,0});
			b1=(u1.x+u1.y+u1.z)/3.0;
			u1 = (float3D){b0+(u1.x-b1),b0+(u1.y-b1),b0+(u1.z-b1)};
			C1[T[i].a].x+=u1.x-u0.x;
			C1[T[i].b].x+=u1.y-u0.y;
			C1[T[i].c].x+=u1.z-u0.z;

			b0=(v0.x+v0.y+v0.z)/3.0;
			v1 = mth_matfloat3D(mx,(float3D){dv1.x,dv1.y,0});
			b1=(v1.x+v1.y+v1.z)/3.0;
			v1 = (float3D){b0+(v1.x-b1),b0+(v1.y-b1),b0+(v1.z-b1)};
			C1[T[i].a].y+=v1.x-v0.x;
			C1[T[i].b].y+=v1.y-v0.y;
			C1[T[i].c].y+=v1.z-v0.z;			
		}
		for(i=0;i<(*m).np;i++)
		{
			C1[i]=sca3D(C1[i],1/(float)NT[i].n);
			C[i].x+=C1[i].x*0.01;
			C[i].y+=C1[i].y*0.01;
		}
	}
	while(!Button());
	for(i=0;i<(*m).np;i++)
		C[i].z = 0;
	
	msh_deleteData(m,c);
}
void em_textureTestConformal(MeshPtr m)
{
	// for each vertex, compute the gradient of the texture coordinates u,v
	// then show how far from orthogonal they are
	ContainerRec	c;
	float3D		*p;
	int3D		*T;
	NTriRec		*NT;
	float3D		*C,*C1;
	char		*M;
	int		i,j;
	float3D		gt[3],gb[3];			// geometric triangle and orthonorm3Dal base
	float3D		du0,dv0;
	float3D		u0,v0;
	float		mx[9];
	float		dot,kcoord=20;
	int		*arr;

	arr = (int *)*g_em_arrH;
	
	if(m==NULL)
		return;

	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	C = msh_getTexturePtr(m);
	M = msh_getLabelPtr(m);
	
	msh_setContainer(&c,(*m).np,sizeof(float3D),0,"\ptemp");
	C1 = (float3D*)msh_getDataPtr(m,c);
	
	for(i=0;i<(*m).np;i++)
		C1[i]=(float3D){0,0,0};
	for(i=0;i<(*m).nt;i++)
	{
		// local base from the triangle
		gt[0]=p[T[i].a];	gt[1]=p[T[i].b];	gt[2]=p[T[i].c];
		triBase(gt,gb);

		mx[0]=0;								mx[1]=0;								mx[2]=1;
		mx[3]=dot3D(sub3D(gt[1],gt[0]),gb[0]);	mx[4]=dot3D(sub3D(gt[1],gt[0]),gb[1]);	mx[5]=1;
		mx[6]=dot3D(sub3D(gt[2],gt[0]),gb[0]);	mx[7]=dot3D(sub3D(gt[2],gt[0]),gb[1]);	mx[8]=1;
		
		// u coordinates local base
		u0 = sca3D((float3D){C[T[i].a].x,C[T[i].b].x,C[T[i].c].x},kcoord);
		du0 = mth_matsolvefloat3D(mx,u0);
		
		// v coordinates local base
		v0 = sca3D((float3D){C[T[i].a].y,C[T[i].b].y,C[T[i].c].y},kcoord);
		dv0 = mth_matsolvefloat3D(mx,v0);
		
		//conformality
		dot=dot3D(du0,dv0);

		C1[T[i].a].z+=dot/(float)NT[T[i].a].n/20.0;
		C1[T[i].b].z+=dot/(float)NT[T[i].b].n/20.0;
		C1[T[i].c].z+=dot/(float)NT[T[i].c].n/20.0;
	}

	// transform, render, detransform
	for(i=0;i<(*m).np;i++)
	{	C1[i].x=C[i].x;						C1[i].y=C[i].y;
		C[i].x=(C1[i].z>0)?(C1[i].z):1;		C[i].y=(C1[i].z<0)?(-C1[i].z):1;
		C1[i].z=C[i].z;						C[i].z=0;
	}
	render(m,1,FrontWindow(),&pict);
	for(j=0;j<(*m).nt*g_em_cutplane/100.0;j++)
	{
		i=arr[j];
		// local base from the triangle
		gt[0]=p[T[i].a];	gt[1]=p[T[i].b];	gt[2]=p[T[i].c];
		triBase(gt,gb);

		mx[0]=0;								mx[1]=0;								mx[2]=1;
		mx[3]=dot3D(sub3D(gt[1],gt[0]),gb[0]);	mx[4]=dot3D(sub3D(gt[1],gt[0]),gb[1]);	mx[5]=1;
		mx[6]=dot3D(sub3D(gt[2],gt[0]),gb[0]);	mx[7]=dot3D(sub3D(gt[2],gt[0]),gb[1]);	mx[8]=1;
		
		// u coordinates local base
		u0 = sca3D((float3D){C[T[i].a].x,C[T[i].b].x,C[T[i].c].x},kcoord);
		du0 = mth_matsolvefloat3D(mx,u0);
		
		// v coordinates local base
		v0 = sca3D((float3D){C[T[i].a].y,C[T[i].b].y,C[T[i].c].y},kcoord);
		dv0 = mth_matsolvefloat3D(mx,v0);

		// draw vectors		
		ForeColor(redColor);
		MoveTo3D(gt[0]);LineTo3D(add3D(gt[0],add3D(sca3D(gb[0],du0.x),sca3D(gb[1],du0.y))));
		ForeColor(greenColor);
		MoveTo3D(gt[0]);LineTo3D(add3D(gt[0],add3D(sca3D(gb[0],dv0.x),sca3D(gb[1],dv0.y))));
	}
	for(i=0;i<(*m).np;i++){	C[i]=C1[i];C[i].z=0;}
	
	msh_deleteData(m,c);
}
float3D em_textureEvalDeriv(float3D aa, float3D bb, float3D cc,int i,double alf)
{
	float3D	zero={0,0,0},z,y;
	float3D	a,b,c;
	float2D	B,C;
	double	A2;
	
	switch(i)
	{
		case 1:	b=sub3D(bb,aa);	c=sub3D(cc,aa);	a=zero;	break;
		case 2:	b=sub3D(cc,bb);	c=sub3D(aa,bb);	a=zero;	break;
		case 3:	b=sub3D(aa,cc);	c=sub3D(bb,cc);	a=zero;	break;
	}
	bb=sub3D(bb,aa);
	cc=sub3D(cc,aa);
	aa=zero;
	switch(0)
	{
		case 0:
			z=(float3D){0,sin(alf),cos(alf)};
			y=(float3D){0,cos(alf),-sin(alf)};
			break;
		case 1:	
			z=sca3D(bb,1/norm3D(bb));
			y=sub3D(cc,sca3D(bb,dot3D(cc,bb)/norm3D(bb)));
			y=sca3D(y,1/norm3D(y));
			break;
	}
	
	B=(float2D){dot3D(b,z),dot3D(b,y)};
	C=(float2D){dot3D(c,z),dot3D(c,y)};
	
	A2=B.x*C.y-B.y*C.x;
	
	return (float3D){(B.y-C.y)/A2,(C.x-B.x)/A2,1};
	
}
void em_textureConformalFE(MeshPtr m)
{
	float3D		*p;
	int3D		*T;
	NTriRec		*NT;
	int2D		*E;
	NEdgeRec	*NE;
	float3D		*C;
	char		*M;
	int			i,j,indx;
	int			ixmax,ax,bx,cx;
	int			ixmin,an,bn,cn;
	float3D		a,b,c;
	float3D		ir,jr,is,js;
	int			pstack[SIZESTACK],estack[SIZESTACK];
	double		coef;
	double		*D,*N,*R,*I,*xR,*xI;
	ContainerRec	Dcont,Ncont;
	int			iter=0;
	double		err=0;

	if(m==NULL)	return;
	
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	E = msh_getEdgesPtr(m);
	NE = msh_getNeighborEdgesPtr(m);
	C = msh_getTexturePtr(m);
	M = msh_getLabelPtr(m);
	
	//0. find a pole (the frontal pole)
	ixmax=ixmin=0;
	for(i=0;i<(*m).np;i++)
	{
		if(p[ixmax].x>p[i].x)	ixmax=i;
		if(p[ixmin].x<p[i].x)	ixmin=i;
	}

	//---------------------
	//1. Border conditions
	//---------------------
	R=dvector_new((*m).np);	xR=dvector_new((*m).np);	
	I=dvector_new((*m).np);	xI=dvector_new((*m).np);
	
	// 1.1 neuman condition for the MAX pole
	ax=T[NT[ixmax].t[0]].a;
	bx=T[NT[ixmax].t[0]].b;
	cx=T[NT[ixmax].t[0]].c;
	
	R[ax]=txtr_evalderiv(p[ax],p[bx],p[cx],1,0).x;
	R[bx]=txtr_evalderiv(p[ax],p[bx],p[cx],2,0).x;
	R[cx]=txtr_evalderiv(p[ax],p[bx],p[cx],3,0).x;
	
	I[ax]=txtr_evalderiv(p[ax],p[bx],p[cx],1,0).y;
	I[bx]=txtr_evalderiv(p[ax],p[bx],p[cx],2,0).y;
	I[cx]=txtr_evalderiv(p[ax],p[bx],p[cx],3,0).y;
	
	// 1.2 neuman condition for the MIN pole
	an=T[NT[ixmin].t[0]].a;
	bn=T[NT[ixmin].t[0]].b;
	cn=T[NT[ixmin].t[0]].c;
	
	R[an]=txtr_evalderiv(p[an],p[bn],p[cn],1,0).x;
	R[bn]=txtr_evalderiv(p[an],p[bn],p[cn],2,0).x;
	R[cn]=txtr_evalderiv(p[an],p[bn],p[cn],3,0).x;
	
	I[an]=txtr_evalderiv(p[an],p[bn],p[cn],1,0).y;
	I[bn]=txtr_evalderiv(p[an],p[bn],p[cn],2,0).y;
	I[cn]=txtr_evalderiv(p[an],p[bn],p[cn],3,0).y;

	//-----------------------
	// 2. Coefficient matrix
	//-----------------------
	msh_setContainer(&Dcont,(*m).np,sizeof(double),0,"\pdiag");
	D=(double*)msh_getDataPtr(m,Dcont);
	msh_setContainer(&Ncont,(*m).nt*3/2,sizeof(double),0,"\pneig");
	N=(double*)msh_getDataPtr(m,Ncont);
	#define cot(x,y) dot3D(x,y)/norm3D(cross3D(x,y))
	for(i=0;i<(*m).np;i++)
	{
		msh_esort(m,i,pstack,estack);
		D[i]=0;
		for(j=0;j<NT[i].n;j++)
		{
			a=p[pstack[(j+2)%NT[i].n]];
			b=p[pstack[(j+1)%NT[i].n]];
			c=p[pstack[j]];
			
			if(E[estack[(j+1)%NT[i].n]].a!=i)	indx=E[estack[(j+1)%NT[i].n]].a;
			else								indx=E[estack[(j+1)%NT[i].n]].b;
			
			ir=sub3D(p[i],a);	ir=sca3D(ir,1/norm3D(ir));
			jr=sub3D(b,a);	jr=sca3D(jr,1/norm3D(jr));
			
			is=sub3D(b,c);	is=sca3D(is,1/norm3D(is));
			js=sub3D(p[i],c);	js=sca3D(js,1/norm3D(js));

			coef = 0.5*(cot(ir,jr)+cot(is,js));

			D[i]+=coef;

			N[estack[(j+1)%NT[i].n]]=-coef;
		}
	}
	#undef cot

	//----------
	// 3. Solve
	//----------
	
	// conjugate gradients in the mesh sparse matrix
	msh_linbcgDN((*m).np,(*m).nt*3/2, D,N,E,R, xR, 1, 0.000001, (*m).np, &iter, &err);
	msh_linbcgDN((*m).np,(*m).nt*3/2, D,N,E,I, xI, 1, 0.000001, (*m).np, &iter, &err);
	msh_drawDN((*m).np,0,xR,N,E,6);
	
	//4. Set texture as coordinates
	for(i=0;i<(*m).np;i++)
	{
		C[i].x=2*xR[i];
		C[i].y=2*xI[i];
		C[i].z=0;
	}
	
	dvector_dispose(R);
	dvector_dispose(I);
	
	msh_deleteData(m,Dcont);
	msh_deleteData(m,Ncont);
	
	printstring("\pmax pole value:");
	fnum_to_string(R[ixmax],10,s);DrawString(s);DrawString("\p,");fnum_to_string(I[ixmax],10,s);DrawString(s);
	printstring("\pmax pole value:");
	fnum_to_string(R[ixmin],10,s);DrawString(s);DrawString("\p,");fnum_to_string(I[ixmin],10,s);DrawString(s);
	
}
*/
#pragma mark _
/*
void em_changeCoordinates(MeshPtr m)
{
	
	int			i;
	float3D		*p;
	NTriRec		*NT;
	int3D		*T;
	float		c[9]={0,0,0, 0,0,0, 0,0,0};
	
	p = msh_getPointsPtr(m);
	T = msh_getTrianglesPtr(m);
	NT = msh_getNeighborTrianglesPtr(m);
	
	for(i=0;i<(*m).np;i++)
	{
		c[0]=p[i].x-(*m).center.x;
		c[1]=p[i].y-(*m).center.y;
		c[2]=p[i].z-(*m).center.z;
		
		mMat(c,g_m3d_xyz,c);

		p[i] = (float3D){c[0]+g_m3d_center.x,c[1]+g_m3d_center.y,c[2]+g_m3d_center.z};
	}
	
	updatebbox(m);
	(*m).center = sca3D(add3D((*m).bbox.ide,(*m).bbox.siz),0.5);
}
*/
void em_textureSmoothVerticesData(MeshPtr m, float *C)
{
	float3D		*p;
	int2D		*E;
	NEdgeRec	*NE;
	int			i,j,indx;
	float3D		a,b,c;
	float3D		ir,jr,is,js;
	int			pstack[SIZESTACK],estack[SIZESTACK];
	double		coef;
	double		*N;//*D;
	double		tau=0.5;
	int			iter,MaxIter=1;
	float		*CC;
	double		sum,max;

	if(m==NULL)	return;
	
	p = msh_getPointsPtr(m);		// points of the mesh
	E = msh_getEdgesPtr(m);			// Edges of the mesh
	NE = msh_getNeighborEdgesPtr(m);// Neighbor Edges for each point

	//D=dvector_new((*m).np);			// Diagonal
	N=dvector_new((*m).nt*3/2);	// Neighbor (non diagonal)

	// 1. Matrix for laplace-beltrami in finite elements
	#define cot(x,y) dot3D(x,y)/norm3D(cross3D(x,y))
	for(i=0;i<(*m).np;i++)
	{
		msh_esort(m,i,pstack,estack);
		//D[i]=0;
		for(j=0;j<NE[i].n;j++)
		{
			a=p[pstack[(j+2)%NE[i].n]];
			b=p[pstack[(j+1)%NE[i].n]];
			c=p[pstack[j]];
			
			if(E[estack[(j+1)%NE[i].n]].a!=i)
				indx=E[estack[(j+1)%NE[i].n]].a;
			else
				indx=E[estack[(j+1)%NE[i].n]].b;
			
			ir=sub3D(p[i],a);	ir=sca3D(ir,1/norm3D(ir));
			jr=sub3D(b,a);	jr=sca3D(jr,1/norm3D(jr));
			
			is=sub3D(b,c);	is=sca3D(is,1/norm3D(is));
			js=sub3D(p[i],c);	js=sca3D(js,1/norm3D(js));

			coef = 0.5*(cot(ir,jr)+cot(is,js));
			
			//D[i]+=coef;
			N[estack[(j+1)%NE[i].n]]=-coef;
		}
	}
	#undef cot
	
	CC=fvector_new((*m).np);
	// 2. Finite difference iterations (diffusion equation)
	for(iter=0;iter<MaxIter;iter++)
	{
		max=C[0];
		for(i=0;i<(*m).np;i++)		// ying: compute on C, store on CC
		{
			sum=0;
			coef=0;
			for(j=0;j<NE[i].n;j++)
			{
				if(E[NE[i].e[j]].a==i)
					sum+=N[NE[i].e[j]]*C[E[NE[i].e[j]].b];
				else
					sum+=N[NE[i].e[j]]*C[E[NE[i].e[j]].a];
				coef+=N[NE[i].e[j]];
			}
			CC[i]=C[i]-tau*(C[i]-sum/coef);
			
			if(fabs(CC[i])>max) max=fabs(CC[i]);
		}
		for(i=0;i<(*m).np;i++)		// yang: compute on CC, store on C
		{
			sum=0;
			coef=0;
			for(j=0;j<NE[i].n;j++)
			{
				if(E[NE[i].e[j]].a==i)
					sum+=N[NE[i].e[j]]*CC[E[NE[i].e[j]].b];
				else
					sum+=N[NE[i].e[j]]*CC[E[NE[i].e[j]].a];
				coef+=N[NE[i].e[j]];
			}
			C[i]=CC[i]-tau*(CC[i]-sum/coef);
			
			C[i]*=1/max;
		}
	}
	fvector_dispose(CC);
	dvector_dispose(N);
}
#pragma mark -
void em_textureDistortion(MeshPtr m, MeshPtr m0, float *C)
{
    int3D	*T=msh_getTrianglesPtr(m);
    NTriRec	*NT=msh_getNeighborTrianglesPtr(m);
    int	i;
    float a,a0;
    double d,mind,maxd;
    
    for(i=0;i<(*m).np;i++) C[i]=0;
    for(i=0;i<(*m).nt;i++)
    {
        a=tTriArea(m,i);
        a0=tTriArea(m0,i);
        d=a/a0;
        
        if(i==0) mind=maxd=d;
        else
        if(d>maxd)	maxd=d;
        else
        if(d<mind)	mind=d;
        
        C[T[i].a]+=d;
        C[T[i].b]+=d;
        C[T[i].c]+=d;
    }
    
    maxd=5.0;
    for(i=0;i<(*m).np;i++)
        C[i]=(C[i]/(double)NT[i].n-mind)/(maxd-mind);
}