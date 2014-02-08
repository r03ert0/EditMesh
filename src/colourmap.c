#include "colourmap.h"
void colourmap(float val, unsigned char *c, int index)
{
	int		n=192,i=val*(n-1);
	float	*cm;
	
	switch(index)
	{
		case AUTUMN:	cm=(float*)autumn;		break;
		case BONE:		cm=(float*)bone;		break;
		case WINTER:	cm=(float*)winter;		break;
		case HOT:		cm=(float*)hot;			break;
		case WATER:		cm=(float*)water;		break;
	
		case JET:		cm=(float*)jet;			break;
		case NEGPOS:	cm=(float*)negpos;		break;
		case GREY:		cm=(float*)gray;		break;

		case RED:		cm=(float*)red;			break;
		case GREEN:		cm=(float*)green;		break;
		case BLUE:		cm=(float*)blue;		break;
	}
	c[0]=255*(cm[    i]+(val-i/(float)(n-1))*(cm[    i+1]-cm[    i]));
	c[1]=255*(cm[  n+i]+(val-i/(float)(n-1))*(cm[  n+i+1]-cm[  n+i]));
	c[2]=255*(cm[2*n+i]+(val-i/(float)(n-1))*(cm[2*n+i+1]-cm[2*n+i]));
}