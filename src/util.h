/*
 *  util.h
 *  EditMesh
 *
 *  Created by roberto on Thu Sep 18 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __util__
#define __util__

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

float getnumfromtxt(int *x, char *t, char *stop);
int getcharfromtxt(int *x, char *t, char *stop);
void getwordfromtxt(int *x, char *t, char *word);
void flookforword(int *x, char *t,char *word, int length);

void initColourmaps(void);
void colourFromColourmap(float index, float *c, int cm);
void configureMonochromaticColourmap(float *c);

int getint(int *x, char *b, int endian);
float getfloat(int *x, char *b, int endian);
int fgetint(FILE *f, int endian);
float fgetfloat(FILE *f, int endian);

#endif