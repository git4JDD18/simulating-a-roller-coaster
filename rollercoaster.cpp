#include "stdafx.h"
#include <pic.h>
#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <GL/glu.h>
#include <GL/glut.h>
#include <cmath>

using namespace std;


int splineIndex = 0;	// Index of the spline being viewed by gluLookAt()
int pointIndex = 0;		// Index of the point on the spline being viewed by gluLookAt()
double u=0;

int g_iMenuId;
int g_vMousePos[2] = {0, 0};
int g_iLeftMouseButton = 0;    /* 1 if pressed, 0 if not */
int g_iMiddleMouseButton = 0;
int g_iRightMouseButton = 0;

typedef enum { ROTATE, TRANSLATE, SCALE } CONTROLSTATE;
CONTROLSTATE g_ControlState = ROTATE;



/* state of the world */
float g_vLandRotate[3] = {0.0, 0.0, 0.0};
float g_vLandTranslate[3] = {0.0, 0.0, 0.0};
float g_vLandScale[3] = {1.0, 1.0, 1.0};


Pic *cubeTexUp;	// Sky texture image
Pic *cubeTexSide1;	// left texture image
Pic *cubeTexSide2;	// right texture image
Pic *cubeTexSide3;	// Sky texture image
Pic *cubeTexSide4;	// Sky texture image
Pic *groundTex;	// Ground texture image
Pic *woodTex1;	// Metal texture image
Pic *woodTex2;
GLuint cubeTexIdUp; //The id of the texture
GLuint cubeTexIdSide1; //The id of the texture
GLuint cubeTexIdSide2; //The id of the texture
GLuint cubeTexIdSide3; //The id of the texture
GLuint cubeTexIdSide4; //The id of the texture

GLuint groundTexId; //The id of the texture
GLuint woodTex1Id;	// ID of metal texture
//GLuint woodTex2Id;	// ID of 2nd metal texture

// Variable for writing filename(incremented everytime a file is saved)
int filename_count = 0;
char filename_in[100];	// input file name


// Write a screenshot to the specified filename 
void saveScreenshot (char *filename)
{
  int i;
  Pic *in = NULL;

  if (filename == NULL)
    return;

  //Allocate a picture buffer 
  in = pic_alloc(640, 480, 3, NULL);

  printf("File to save to: %s\n", filename);

  for (i=479; i>=0; i--) {
    glReadPixels(0, 479-i, 640, 1, GL_RGB, GL_UNSIGNED_BYTE,
                 &in->pix[i*in->nx*in->bpp]);
  }

  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);
}



/* represents one control point along the spline */
struct point {
	double x;
	double y;
	double z;
};

struct point currentLookAt;// Stores the current Point for gluLookAt()
struct point controlPointOfSegment[4]; //array of 4 control point of a given segment

// Structure for a vector (Used for tangents, normals , binormals)
struct vector {
	double x;
	double y;
	double z;
};

struct vector tangent, normal, binormal, prevBinormal;	// Store the tangent, normal, binormal values (for gluLookAt)

/* spline struct which contains how many control points, and an array of control points */
struct spline {
	int numControlPoints;
	struct point *points;
};



/* the spline array */
struct spline *g_Splines;

/* total number of splines */
int g_iNumOfSplines;


int loadSplines(char *argv) {
	char *cName = (char *)malloc(128 * sizeof(char));
	FILE *fileList;
	FILE *fileSpline;
	int iType, i = 0, j, iLength;

	/* load the track file */
	 fopen_s(&fileList,argv, "r");
	if (fileList == NULL) {
		printf ("can't open file\n");
		exit(1);
	}
  
	/* stores the number of splines in a global variable */
	fscanf_s(fileList, "%d", &g_iNumOfSplines);

	g_Splines = (struct spline *)malloc(g_iNumOfSplines * sizeof(struct spline));

	/* reads through the spline files */
	for (j = 0; j < g_iNumOfSplines; j++) {
		i = 0;
		fscanf_s(fileList, "%s",cName,128);
		fopen_s(&fileSpline,cName, "r");

		if (fileSpline == NULL) {
			printf ("can't open file\n");
			exit(1);
		}

		/* gets length for spline file */
		fscanf_s(fileSpline, "%d %d", &iLength, &iType);

		/* allocate memory for all the points */
		g_Splines[j].points = (struct point *)malloc(iLength * sizeof(struct point));
		g_Splines[j].numControlPoints = iLength;

		/* saves the data to the struct */
		while (fscanf_s(fileSpline, "%lf %lf %lf", 
			&g_Splines[j].points[i].x, 
			&g_Splines[j].points[i].y, 
			&g_Splines[j].points[i].z) != EOF) {
			i++;
		}
	}

	free(cName);
	//free(g_Splines);
	//for (j = 0; j < g_iNumOfSplines; j++) 
	//{
	//	free(g_Splines[j].points);
	//}
	return 0;
}


// Loads the texture pointed by the image into memory
GLuint loadTexture(Pic* image) {
	GLuint textureId;
	glGenTextures(1, &textureId); 
	glBindTexture(GL_TEXTURE_2D, textureId); 
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,image->nx, image->ny,0,GL_RGB,GL_UNSIGNED_BYTE,image->pix);               
	return textureId; 
}


void cubeTexture()
{
	cubeTexUp=jpeg_read("sky.jpg",NULL);
	cubeTexIdUp=loadTexture(cubeTexUp);
	delete cubeTexUp;

	cubeTexSide1=jpeg_read("snow.jpg",NULL);
	cubeTexIdSide1=loadTexture(cubeTexSide1);
	delete cubeTexSide1;

	cubeTexSide2=jpeg_read("snow.jpg",NULL);
	cubeTexIdSide2=loadTexture(cubeTexSide1);
	delete cubeTexSide1;

	cubeTexSide3=jpeg_read("snow.jpg",NULL);
	cubeTexIdSide2=loadTexture(cubeTexSide2);
	delete cubeTexSide2;

	cubeTexSide4=jpeg_read("snow.jpg",NULL);
	cubeTexIdSide3=loadTexture(cubeTexSide2);
	delete cubeTexSide2;


	groundTex=jpeg_read("grass.jpg",NULL);
	groundTexId = loadTexture(groundTex);
	delete groundTex;

	woodTex1=jpeg_read("rail.jpg",NULL);
	woodTex1Id = loadTexture(woodTex1);
	delete woodTex1;

}


void myinit()	// initialize function
{
	/* setup gl view here */
	glClearColor(0.0, 0.0, 0.0, 0.0);
	// Enable smooth shading
	glShadeModel(GL_SMOOTH);
	
	// Enable GL_DEPTH_TEST
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING); //Enable lighting
	glEnable(GL_LIGHT0); 
	glEnable(GL_LIGHT1); 
	cubeTexture();
	
}


void loadSplinePoints(struct point *pCPList,int splineIndex, int pointIndex)
{
	// load the four points of the spline curve
	pCPList[0].x = g_Splines[splineIndex].points[pointIndex].x;
	pCPList[0].y = g_Splines[splineIndex].points[pointIndex].y;
	pCPList[0].z = g_Splines[splineIndex].points[pointIndex].z;

	pCPList[1].x = g_Splines[splineIndex].points[pointIndex + 1].x;
	pCPList[1].y = g_Splines[splineIndex].points[pointIndex + 1].y;
	pCPList[1].z = g_Splines[splineIndex].points[pointIndex + 1].z;

	pCPList[2].x = g_Splines[splineIndex].points[pointIndex + 2].x;
	pCPList[2].y = g_Splines[splineIndex].points[pointIndex + 2].y;
	pCPList[2].z = g_Splines[splineIndex].points[pointIndex + 2].z;

	pCPList[3].x = g_Splines[splineIndex].points[pointIndex + 3].x;
	pCPList[3].y = g_Splines[splineIndex].points[pointIndex + 3].y;
	pCPList[3].z = g_Splines[splineIndex].points[pointIndex + 3].z;

}


void computeCurrentLookAt(struct point *pCPList, double u)
{
	

	//double u2=u*u;
	//double u3=u*u*u;

	// Compute the xyz values of the point for the current value of u
	float x = (-0.5*u*u*u+u*u-0.5*u)*pCPList[0].x+(1.5*u*u*u-2.5*u*u+1)*pCPList[1].x+(-1.5*u*u*u+2*u*u+0.5*u)*pCPList[2].x+(0.5*u*u*u-0.5*u*u)*pCPList[3].x;
	float y = (-0.5*u*u*u+u*u-0.5*u)*pCPList[0].y+(1.5*u*u*u-2.5*u*u+1)*pCPList[1].y+(-1.5*u*u*u+2*u*u+0.5*u)*pCPList[2].y+(0.5*u*u*u-0.5*u*u)*pCPList[3].y;
	float z = (-0.5*u*u*u+u*u-0.5*u)*pCPList[0].z+(1.5*u*u*u-2.5*u*u+1)*pCPList[1].z+(-1.5*u*u*u+2*u*u+0.5*u)*pCPList[2].z+(0.5*u*u*u-0.5*u*u)*pCPList[3].z;
	// Store these values in currLookAt
	currentLookAt.x=x;
	currentLookAt.y=y;
	currentLookAt.z=z;

}


void TNB_Calculation(struct point *pCPList,double u)
{
	// If first point, then compute the normal and binormal using an arbitrary vector (in this case (0,1,0))
	if (u == 0 && pointIndex == 0 )
	{
		//Tangent Calculation
		float tx0=-0.5*pCPList[0].x+0.5*pCPList[2].x;
		float ty0=-0.5*pCPList[0].y+0.5*pCPList[2].y;
		float tz0=-0.5*pCPList[0].z+0.5*pCPList[2].z;

		float tangentMagnitude = sqrt(tx0*tx0+ty0*ty0+tz0*tz0);
		tx0=tx0/tangentMagnitude;
		ty0=ty0/tangentMagnitude;
		tz0=tz0/tangentMagnitude;

		tangent.x=tx0;
		tangent.y=ty0;
		tangent.z=tz0;


		//Normal Calcuation
		float arbitrary_x0 = 0; 
		float arbitrary_y0 = 1; 
		float arbitrary_z0 = 0;
		

		float nx0=(ty0*arbitrary_z0 - tz0*arbitrary_y0);
		float ny0=(tz0*arbitrary_x0 - tx0*arbitrary_z0);
		float nz0=(tx0*arbitrary_y0 - ty0*arbitrary_x0);
		float normalMagnitude = sqrt(nx0*nx0+ny0*ny0+nz0*nz0);

		nx0=nx0/normalMagnitude;
		ny0=ny0/normalMagnitude;
		nz0=nz0/normalMagnitude;

		normal.x=nx0;
		normal.y=ny0;
		normal.z=nz0;

		//Binormal Calculation

		float bx0=(ty0*nz0-tz0*ny0);
		float by0=(tz0*nx0-tx0*nz0);
		float bz0=(tx0*ny0-ty0*nx0);

		float binormalMagnitude = sqrt(bx0*bx0+by0*by0+bz0*bz0);
		bx0=bx0/binormalMagnitude;
		by0=by0/binormalMagnitude;
		bz0=bz0/binormalMagnitude;

		binormal.x=bx0;
		binormal.y=by0;
		binormal.z=bz0;

	}

	// For other points, compute the normal and binormal using the value of the binormal at the previous point
	else
	{

			float bx0=prevBinormal.x;
			float by0=prevBinormal.y;
			float bz0=prevBinormal.z;

			float tx = (-0.5*3*u*u+2*u-0.5*1)*pCPList[0].x+(1.5*3*u*u-2.5*2*u)*pCPList[1].x+(-1.5*3*u*u+2*2*u+0.5*1)*pCPList[2].x+(0.5*3*u*u-0.5*2*u)*pCPList[3].x;
			float ty = (-0.5*3*u*u+2*u-0.5*1)*pCPList[0].y+(1.5*3*u*u-2.5*2*u)*pCPList[1].y+(-1.5*3*u*u+2*2*u+0.5*1)*pCPList[2].y+(0.5*3*u*u-0.5*2*u)*pCPList[3].y;
			float tz = (-0.5*3*u*u+2*u-0.5*1)*pCPList[0].z+(1.5*3*u*u-2.5*2*u)*pCPList[1].z+(-1.5*3*u*u+2*2*u+0.5*1)*pCPList[2].z+(0.5*3*u*u-0.5*2*u)*pCPList[3].z;

			float tangentMagnitude = sqrt(tx*tx+ty*ty+tz*tz);

			// Normalize the values
			tx=tx/tangentMagnitude;
			ty=ty/tangentMagnitude;
			tz=tz/tangentMagnitude;

			float nx=(by0*tz-bz0*ty);
			float ny=(bz0*tx-bx0*tz);
			float nz=(bx0*ty-by0*tx);
			float normalMagnitude = sqrt(nx*nx+ny*ny+nz*nz);

			// Normalize the values
			nx=nx/normalMagnitude;
			ny=ny/normalMagnitude;
			nz=nz/normalMagnitude;

			float bx=(ty*nz-tz*ny);
			float by=(tz*nx-tx*nz);
			float bz=(tx*ny-ty*nx);
			float binormalMagnitude = sqrt(bx*bx+by*by+bz*bz);
			// Normalize the values
			bx=bx/binormalMagnitude;
			by=by/binormalMagnitude;
			bz=bz/binormalMagnitude;

			// Store the values of tangent, normal, binormal in a global structure so that those can be used by the gluLookAt function
			tangent.x=tx;
			tangent.y=ty;
			tangent.z=tz;
		
			normal.x=nx;
			normal.y=ny;
			normal.z=nz;

			binormal.x=bx;
			binormal.y=by;
			binormal.z=bz;

	}
	// Store the value of value of binormal as the value of previous binormal for the next point
	prevBinormal.x=binormal.x;
	prevBinormal.y=binormal.y;
	prevBinormal.z=binormal.z;

}


void drawSpline(int splineNum)
{

	float avX=0,avY=1,avZ=0;	// Arbitrary vector for normal computation 
	float local_x,local_y,local_z;	// xyz coordinates of the point
	float local_tx,local_ty,local_tz,local_magT,local_nx,local_ny,local_nz,local_magN,local_bx,local_by,local_bz,local_magB;	// Tangent, normal and binormal values
	float step = 0.02;	// Stepsize for incrementing u
	float uNext,xNext,yNext,zNext,txNext,tyNext,tzNext,magTNext,nxNext,nyNext,nzNext,magN_Next,bxNext,byNext,bzNext,magBNext;	// Values for the next point on the spline

	float centralH1=0.015,centralH2=0.025,centralW=0.01;	// Height and width of the central support bar
	float connectorH=0.01,connectorW=0.008;

	float railHeight=0.0025, railW1=0.032,railW2=0.028;
	float barH1=0.0025,barH2=0.005,barW=0.04,barL=0.0075;

	struct point l1,l2,l3,l4,r1,r2,r3,r4,f1,f2,f3,f4,back1,back2,back3,back4,lp1,lp2,lp3,lp4,rp1,rp2,rp3,rp4,ctl,ctr;
	struct point cpList[4];

	// For each control point
	for(int i=0; i< g_Splines[splineNum].numControlPoints;i++)
	{
		// Load the four control points
		loadSplinePoints(cpList,splineNum,i);
	
		glEnable(GL_TEXTURE_2D);
		// Enable textures
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,GL_REPLACE);
	
	
		
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, woodTex1Id);

		glBegin(GL_QUADS);
		// Change u from 0 to 1-step
		for(float u=0;u<=1-step;u+=step)
		{
			// Compute the x,y,z values for the current point
			local_x = (-0.5*u*u*u+u*u-0.5*u)*cpList[0].x+(1.5*u*u*u-2.5*u*u+1)*cpList[1].x+(-1.5*u*u*u+2*u*u+0.5*u)*cpList[2].x+(0.5*u*u*u-0.5*u*u)*cpList[3].x;
			local_y = (-0.5*u*u*u+u*u-0.5*u)*cpList[0].y+(1.5*u*u*u-2.5*u*u+1)*cpList[1].y+(-1.5*u*u*u+2*u*u+0.5*u)*cpList[2].y+(0.5*u*u*u-0.5*u*u)*cpList[3].y;
			local_z = (-0.5*u*u*u+u*u-0.5*u)*cpList[0].z+(1.5*u*u*u-2.5*u*u+1)*cpList[1].z+(-1.5*u*u*u+2*u*u+0.5*u)*cpList[2].z+(0.5*u*u*u-0.5*u*u)*cpList[3].z;
			// Compute the tangent values for the current point
			local_tx = (-0.5*3*u*u+2*u-0.5*1)*cpList[0].x+(1.5*3*u*u-2.5*2*u)*cpList[1].x+(-1.5*3*u*u+2*2*u+0.5*1)*cpList[2].x+(0.5*3*u*u-0.5*2*u)*cpList[3].x;
			local_ty = (-0.5*3*u*u+2*u-0.5*1)*cpList[0].y+(1.5*3*u*u-2.5*2*u)*cpList[1].y+(-1.5*3*u*u+2*2*u+0.5*1)*cpList[2].y+(0.5*3*u*u-0.5*2*u)*cpList[3].y;
			local_tz = (-0.5*3*u*u+2*u-0.5*1)*cpList[0].z+(1.5*3*u*u-2.5*2*u)*cpList[1].z+(-1.5*3*u*u+2*2*u+0.5*1)*cpList[2].z+(0.5*3*u*u-0.5*2*u)*cpList[3].z;
			// Normalize the tangent values
			local_magT=sqrt(local_tx*local_tx+local_ty*local_ty+local_tz*local_tz);
			local_tx = local_tx/local_magT;
			local_ty = local_ty/local_magT;
			local_tz = local_tz/local_magT;

			// For the first point point, compute the normal using an arbitray vector
			if(u==0 && i==0 )
			{
				local_nx=(local_ty*0-local_tz*1);
				local_ny=(local_tz*0-local_tx*0);
				local_nz=(local_tx*1-local_ty*0);
				local_magN = sqrt(local_nx*local_nx+local_ny*local_ny+local_nz*local_nz);

				local_nx = local_nx/local_magN;
				local_ny = local_ny/local_magN;
				local_nz = local_nz/local_magN;			
			}
			// For other points, compute the normal using the binormal at the previosu point
			else
			{
				local_nx=(local_by*local_tz-local_bz*local_ty);
				local_ny=(local_bz*local_tx-local_bx*local_tz);
				local_nz=(local_bx*local_ty-local_by*local_tx);
				local_magN = sqrt(local_nx*local_nx+local_ny*local_ny+local_nz*local_nz);
	
				local_nx=local_nx/local_magN;
				local_ny=local_ny/local_magN;
				local_nz=local_nz/local_magN;
			}
			
			// Compute the binormal, normalize its components
			local_bx=(local_ty*local_nz-local_tz*local_ny);
			local_by=(local_tz*local_nx-local_tx*local_nz);
			local_bz=(local_tx*local_ny-local_ty*local_nx);
			local_magB = sqrt(local_bx*local_bx+local_by*local_by+local_bz*local_bz);
			local_bx=local_bx/local_magB;
			local_by=local_by/local_magB;
			local_bz=local_bz/local_magB;

			//Top left point of left rail
			l1.x=local_x+local_nx*railHeight-local_bx*railW1;
			l1.y=local_y+local_ny*railHeight-local_by*railW1;
			l1.z=local_z+local_nz*railHeight-local_bz*railW1;

			//Bottom left point of left rail
			l2.x=local_x-local_nx*railHeight-local_bx*railW1;
			l2.y=local_y-local_ny*railHeight-local_by*railW1;
			l2.z=local_z-local_nz*railHeight-local_bz*railW1;

			//Top Right point of left rail
			l3.x=local_x+local_nx*railHeight-local_bx*railW2;
			l3.y=local_y+local_ny*railHeight-local_by*railW2;
			l3.z=local_z+local_nz*railHeight-local_bz*railW2;

			//Bottom right point of left rail
			l4.x=local_x-local_nx*railHeight-local_bx*railW2;
			l4.y=local_y-local_ny*railHeight-local_by*railW2;
			l4.z=local_z-local_nz*railHeight-local_bz*railW2;

			//Top right point of right rail
			r1.x=local_x+local_nx*railHeight+local_bx*railW1;
			r1.y=local_y+local_ny*railHeight+local_by*railW1;
			r1.z=local_z+local_nz*railHeight+local_bz*railW1;

			//Bottom right point of right railo
			r2.x=local_x-local_nx*railHeight+local_bx*railW1;
			r2.y=local_y-local_ny*railHeight+local_by*railW1;
			r2.z=local_z-local_nz*railHeight+local_bz*railW1;

			//Top left point of right rail
			r3.x=local_x+local_nx*railHeight+local_bx*railW2;
			r3.y=local_y+local_ny*railHeight+local_by*railW2;
			r3.z=local_z+local_nz*railHeight+local_bz*railW2;

			//Bottom left point of right rail
			r4.x=local_x-local_nx*railHeight+local_bx*railW2;
			r4.y=local_y-local_ny*railHeight+local_by*railW2;
			r4.z=local_z-local_nz*railHeight+local_bz*railW2;

			//Top left point of the front face of central support tube
			f1.x=local_x-local_nx*centralH1-local_bx*centralW;
			f1.y=local_y-local_ny*centralH1-local_by*centralW;
			f1.z=local_z-local_nz*centralH1-local_bz*centralW;

			//Bottom left point of the front face of central support tube
			f2.x=local_x-local_nx*centralH2-local_bx*centralW;
			f2.y=local_y-local_ny*centralH2-local_by*centralW;
			f2.z=local_z-local_nz*centralH2-local_bz*centralW;

			//Top Right point of the front face of central support tube
			f3.x=local_x-local_nx*centralH1+local_bx*centralW;
			f3.y=local_y-local_ny*centralH1+local_by*centralW;
			f3.z=local_z-local_nz*centralH1+local_bz*centralW;

			//Bottom right point of the front face of central support tube
			f4.x=local_x-local_nx*centralH2+local_bx*centralW;
			f4.y=local_y-local_ny*centralH2+local_by*centralW;
			f4.z=local_z-local_nz*centralH2+local_bz*centralW;

			// Left point of the central connector
			ctl.x=local_x-local_nx*connectorH-local_bx*connectorW;
			ctl.y=local_y-local_ny*connectorH-local_by*connectorW;
			ctl.z=local_z-local_nz*connectorH-local_bz*connectorW;

			// Right point of the central connector
			ctr.x=local_x-local_nx*connectorH+local_bx*connectorW;
			ctr.y=local_y-local_ny*connectorH+local_by*connectorW;
			ctr.z=local_z-local_nz*connectorH+local_bz*connectorW;
			


			glColor3f(0.5,0.5,0.5);
			
			//glVertex3f(x+bx*0.05,y+by*0.05,z+bz*0.05);
			//glVertex3f(x-bx*0.05,y-by*0.05,z-bz*0.05);



			
			// Repeat the procedure for the next point on the spline
			uNext = u+step;
			xNext = (-0.5*uNext*uNext*uNext+uNext*uNext-0.5*uNext)*cpList[0].x+(1.5*uNext*uNext*uNext-2.5*uNext*uNext+1)*cpList[1].x+(-1.5*uNext*uNext*uNext+2*uNext*uNext+0.5*uNext)*cpList[2].x+(0.5*uNext*uNext*uNext-0.5*uNext*uNext)*cpList[3].x;
			yNext = (-0.5*uNext*uNext*uNext+uNext*uNext-0.5*uNext)*cpList[0].y+(1.5*uNext*uNext*uNext-2.5*uNext*uNext+1)*cpList[1].y+(-1.5*uNext*uNext*uNext+2*uNext*uNext+0.5*uNext)*cpList[2].y+(0.5*uNext*uNext*uNext-0.5*uNext*uNext)*cpList[3].y;
			zNext = (-0.5*uNext*uNext*uNext+uNext*uNext-0.5*uNext)*cpList[0].z+(1.5*uNext*uNext*uNext-2.5*uNext*uNext+1)*cpList[1].z+(-1.5*uNext*uNext*uNext+2*uNext*uNext+0.5*uNext)*cpList[2].z+(0.5*uNext*uNext*uNext-0.5*uNext*uNext)*cpList[3].z;
				// Compute the tangent values for the current point
			txNext = (-0.5*3*uNext*uNext+2*uNext-0.5*1)*cpList[0].x+(1.5*3*uNext*u-2.5*2*uNext)*cpList[1].x+(-1.5*3*uNext*uNext+2*2*uNext+0.5*1)*cpList[2].x+(0.5*3*uNext*uNext-0.5*2*uNext)*cpList[3].x;
			tyNext = (-0.5*3*uNext*uNext+2*uNext-0.5*1)*cpList[0].y+(1.5*3*uNext*u-2.5*2*uNext)*cpList[1].y+(-1.5*3*uNext*uNext+2*2*uNext+0.5*1)*cpList[2].y+(0.5*3*uNext*uNext-0.5*2*uNext)*cpList[3].y;
			tzNext = (-0.5*3*uNext*uNext+2*uNext-0.5*1)*cpList[0].z+(1.5*3*uNext*u-2.5*2*uNext)*cpList[1].z+(-1.5*3*uNext*uNext+2*2*uNext+0.5*1)*cpList[2].z+(0.5*3*uNext*uNext-0.5*2*uNext)*cpList[3].z;
			// Normalize the tangent values
			magTNext=sqrt(txNext*txNext+tyNext*tyNext+tzNext*tzNext);
			txNext=txNext/magTNext;
			tyNext=tyNext/magTNext;
			tzNext=tzNext/magTNext;

			nxNext=(local_by*tzNext-local_bz*tyNext);
			nyNext=(local_bz*txNext-local_bx*tzNext);
			nzNext=(local_bx*tyNext-local_by*txNext);
			magN_Next = sqrt(nxNext*nxNext+nyNext*nyNext+nzNext*nzNext);
	
			nxNext=nxNext/magN_Next;
			nyNext=nyNext/magN_Next;
			nzNext=nzNext/magN_Next;
			
			
			// Compute the binormal, normalize its components
			bxNext=(tyNext*nzNext-tzNext*nyNext);
			byNext=(tzNext*nxNext-txNext*nzNext);
			bzNext=(txNext*nyNext-tyNext*nxNext);
			magBNext = sqrt(bxNext*bxNext+byNext*byNext+bzNext*bzNext);
			bxNext=bxNext/magBNext;
			byNext=byNext/magBNext;
			bzNext=bzNext/magBNext;


			back1.x=xNext-local_nx*centralH1-local_bx*centralW;
			back1.y=yNext-local_ny*centralH1-local_by*centralW;
			back1.z=zNext-local_nz*centralH1-local_bz*centralW;

			back2.x=xNext-local_nx*centralH2-local_bx*centralW;
			back2.y=yNext-local_ny*centralH2-local_by*centralW;
			back2.z=zNext-local_nz*centralH2-local_bz*centralW;

			
			back3.x=xNext-local_nx*centralH1+local_bx*centralW;
			back3.y=yNext-local_ny*centralH1+local_by*centralW;
			back3.z=zNext-local_nz*centralH1+local_bz*centralW;

			back4.x=xNext-local_nx*centralH2+local_bx*centralW;
			back4.y=yNext-local_ny*centralH2+local_by*centralW;
			back4.z=zNext-local_nz*centralH2+local_bz*centralW;


			lp1.x=xNext+local_nx*railHeight-local_bx*railW1+local_tx*0.005;
			lp1.y=yNext+local_ny*railHeight-local_by*railW1+local_ty*0.005;
			lp1.z=zNext+local_nz*railHeight-local_bz*railW1+local_tz*0.005;

			lp2.x=xNext-local_nx*railHeight-local_bx*railW1+local_tx*0.005;
			lp2.y=yNext-local_ny*railHeight-local_by*railW1+local_ty*0.005;
			lp2.z=zNext-local_nz*railHeight-local_bz*railW1+local_tz*0.005;

			lp3.x=xNext+local_nx*railHeight-local_bx*railW2+local_tx*0.005;
			lp3.y=yNext+local_ny*railHeight-local_by*railW2+local_ty*0.005;
			lp3.z=zNext+local_nz*railHeight-local_bz*railW2+local_tz*0.005;

			lp4.x=xNext-local_nx*railHeight-local_bx*railW2+local_tx*0.005;
			lp4.y=yNext-local_ny*railHeight-local_by*railW2+local_ty*0.005;
			lp4.z=zNext-local_nz*railHeight-local_bz*railW2+local_tz*0.005;


			rp1.x=xNext+local_nx*railHeight+local_bx*railW1+local_tx*0.005;
			rp1.y=yNext+local_ny*railHeight+local_by*railW1+local_ty*0.005;
			rp1.z=zNext+local_nz*railHeight+local_bz*railW1+local_tz*0.005;

			rp2.x=xNext-local_nx*railHeight+local_bx*railW1+local_tx*0.005;
			rp2.y=yNext-local_ny*railHeight+local_by*railW1+local_ty*0.005;
			rp2.z=zNext-local_nz*railHeight+local_bz*railW1+local_tz*0.005;

			rp3.x=xNext+local_nx*railHeight+local_bx*railW2+local_tx*0.005;
			rp3.y=yNext+local_ny*railHeight+local_by*railW2+local_ty*0.005;
			rp3.z=zNext+local_nz*railHeight+local_bz*railW2+local_tz*0.005;

			rp4.x=xNext-local_nx*railHeight+local_bx*railW2+local_tx*0.005;
			rp4.y=yNext-local_ny*railHeight+local_by*railW2+local_ty*0.005;
			rp4.z=zNext-local_nz*railHeight+local_bz*railW2+local_tz*0.005;

			
	
			// Left rail
			
			glTexCoord2f(0.0f, 0.5f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(l1.x,l1.y,l1.z);
			glTexCoord2f(0.0f, 1.0f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(lp1.x,lp1.y,lp1.z);
			glTexCoord2f(0.5f, 1.0f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(lp2.x,lp2.y,lp2.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(l2.x,l2.y,l2.z);
			
			
			
			glTexCoord2f(0.0f, 0.5f);
			glNormal3f(local_bx,local_by,local_bz);
			glVertex3f(l3.x,l3.y,l3.z);
			glTexCoord2f(0.0f, 1.0f);
			glNormal3f(local_bx,local_by,local_bz);
			glVertex3f(lp3.x,lp3.y,lp3.z);
			glTexCoord2f(0.5f, 1.0f);
			glNormal3f(local_bx,local_by,local_bz);
			glVertex3f(lp4.x,lp4.y,lp4.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(local_bx,local_by,local_bz);
			glVertex3f(l4.x,l4.y,l4.z);
			
			
			glTexCoord2f(0.0f, 0.5f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(l1.x,l1.y,l1.z);
			glTexCoord2f(0.0f, 1.0f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(lp1.x,lp1.y,lp1.z);
			glTexCoord2f(0.5f, 1.0f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(lp3.x,lp3.y,lp3.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(l3.x,l3.y,l3.z);
			
			glTexCoord2f(0.0f, 0.5f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(l2.x,l2.y,l2.z);
			glTexCoord2f(0.0f, 1.0f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(lp2.x,lp2.y,lp2.z);
			glTexCoord2f(0.5f, 1.0f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(lp4.x,lp4.y,lp4.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(l4.x,l4.y,l4.z);
			
			
			// Right rail					
			glTexCoord2f(0.0f, 0.5f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(r1.x,r1.y,r1.z);
			glTexCoord2f(0.0f, 1.0f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(rp1.x,rp1.y,rp1.z);
			glTexCoord2f(0.5f, 1.0f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(rp2.x,rp2.y,rp2.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(r2.x,r2.y,r2.z);
			
			glTexCoord2f(0.0f, 0.5f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(r3.x,r3.y,r3.z);
			glTexCoord2f(0.0f, 1.0f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(rp3.x,rp3.y,rp3.z);			
			glTexCoord2f(0.5f, 1.0f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(rp4.x,rp4.y,rp4.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(r4.x,r4.y,r4.z);
			
			
			glTexCoord2f(0.0f, 0.5f);
			glNormal3f(local_bx,local_by,local_bz);
			glVertex3f(r1.x,r1.y,r1.z);
			glTexCoord2f(0.0f, 1.0f);
			glNormal3f(local_bx,local_by,local_bz);
			glVertex3f(rp1.x,rp1.y,rp1.z);
			glTexCoord2f(0.5f, 1.0f);
			glNormal3f(local_bx,local_by,local_bz);
			glVertex3f(rp3.x,rp3.y,rp3.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(local_bx,local_by,local_bz);
			glVertex3f(r3.x,r3.y,r3.z);
			
			glTexCoord2f(0.0f, 0.5f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(r2.x,r2.y,r2.z);
			glTexCoord2f(0.0f, 1.0f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(rp2.x,rp2.y,rp2.z);
			glTexCoord2f(0.5f, 1.0f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(rp4.x,rp4.y,rp4.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(r4.x,r4.y,r4.z);
			
	
			
		// Central Support tube
			
	
			glColor3f(1.0,0.4,0.4);

			glTexCoord2f(0.5f, 0.5f);
			glVertex3f(f1.x,f1.y,f1.z);
			glTexCoord2f(0.0f, 0.5f);
			glVertex3f(f3.x,f3.y,f3.z);
			glTexCoord2f(0.5f, 0.0f);
			glVertex3f(f4.x,f4.y,f4.z);
			glTexCoord2f(0.0f, 0.0f);
			glVertex3f(f2.x,f2.y,f2.z);
			
			glTexCoord2f(0.0f, 0.5f);
			glVertex3f(back1.x,back1.y,back1.z);
			glTexCoord2f(0.5f, 0.5f);
			glVertex3f(back3.x,back3.y,back3.z);
			glTexCoord2f(0.5f, 0.0f);
			glVertex3f(back4.x,back4.y,back4.z);
			glTexCoord2f(0.0f, 0.0f);
			glVertex3f(back2.x,back2.y,back2.z);
			
			glTexCoord2f(0.0f, 0.5f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(f1.x,f1.y,f1.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(back1.x,back1.y,back1.z);
			glTexCoord2f(0.5f, 0.0f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(back2.x,back2.y,back2.z);
			glTexCoord2f(0.0f, 0.0f);
			glNormal3f(-local_bx,-local_by,-local_bz);
			glVertex3f(f2.x,f2.y,f2.z);
			
			
			
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(local_bx,local_by,local_bz);
			glVertex3f(f3.x,f3.y,f3.z);
			glNormal3f(local_bx,local_by,local_bz);
			glTexCoord2f(0.0f, 0.5f);
			glVertex3f(back3.x,back3.y,back3.z);
			glTexCoord2f(0.5f, 0.0f);
			glNormal3f(local_bx,local_by,local_bz);
			glVertex3f(back4.x,back4.y,back4.z);
			glTexCoord2f(0.0f, 0.0f);
			glNormal3f(local_bx,local_by,local_bz);
			glVertex3f(f4.x,f4.y,f4.z);
			
			
			glTexCoord2f(0.0f, 0.5f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(f1.x,f1.y,f1.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(back1.x,back1.y,back1.z);			
			glTexCoord2f(0.5f, 0.0f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(back3.x,back3.y,back3.z);
			glTexCoord2f(0.0f, 0.0f);
			glNormal3f(local_nx,local_ny,local_nz);
			glVertex3f(f3.x,f3.y,f3.z);
			
			
			glTexCoord2f(0.0f, 0.5f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(f2.x,f2.y,f2.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(f4.x,f4.y,f4.z);
			glTexCoord2f(0.5f, 0.0f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(back2.x,back2.y,back2.z);
			glTexCoord2f(0.0f, 0.0f);
			glNormal3f(-local_nx,-local_ny,-local_nz);
			glVertex3f(back4.x,back4.y,back4.z);
			
			// Connectors between rail and central support
			glTexCoord2f(0.5f, 0.0f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(l4.x,l4.y,l4.z);
			glTexCoord2f(1.0f, 0.0f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(l2.x,l2.y,l2.z);
			glTexCoord2f(1.0f, 0.5f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(f2.x,f2.y,f2.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(ctl.x,ctl.y,ctl.z);

			

			glTexCoord2f(0.5f, 0.0f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(r4.x,r4.y,r4.z);
			glTexCoord2f(1.0f, 0.0f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(r2.x,r2.y,r2.z);
			glTexCoord2f(1.0f, 0.5f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(f4.x,f4.y,f4.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(ctr.x,ctr.y,ctr.z);

			glTexCoord2f(0.5f, 0.0f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(f1.x,f1.y,f1.z);
			glTexCoord2f(1.0f, 0.0f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(ctl.x,ctl.y,ctl.z);
			glTexCoord2f(1.0f, 0.5f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(ctr.x,ctr.y,ctr.z);
			glTexCoord2f(0.5f, 0.5f);
			glNormal3f(-local_tx,-local_ty,-local_tz);
			glVertex3f(f3.x,f3.y,f3.z);
			
			
			
			

			
		}
		glEnd();
		glDisable(GL_TEXTURE_2D);

	}
}


void doIdle(int value)
{	
	// Compute the xyz, tangent, normal, binormal values at the current point
	loadSplinePoints(controlPointOfSegment,splineIndex,pointIndex);
	computeCurrentLookAt(controlPointOfSegment,u);
	TNB_Calculation(controlPointOfSegment,u);
	// Increment the value of u
	u = u +0.05;
	// if u reaches 1, make u 0 and increment numpoint(use next four points on the spline for computation)
	if(u >= 1)
	{
		u = 0;
		pointIndex++;
	}

	// If all control points of the current spline have been used, render next spline.
	if(pointIndex == g_Splines[splineIndex].numControlPoints-4)
	{
		splineIndex++;
		pointIndex = 0;
	}

	if(splineIndex == g_iNumOfSplines)
		splineIndex = 0;

	glutPostRedisplay();
	
	// Call a delay of 65ms
	glutTimerFunc(65,doIdle,0);
}


void menufunc(int value)
{
  switch (value)
  {
    case 0:
      exit(0);
      break;
  }
}


void keyboard(unsigned char key, int x, int y)
{
	// if s is pressed save the file
	if(key == 's'|| key == 'S')
	{
		sprintf_s(filename_in,"%d.jpg",filename_count);
		filename_count++; // increment filename count everytime
		saveScreenshot(filename_in);

	}

}



void reshape(int width, int height)
 
{
	// Enable perspective projection
   glViewport(0, 0, width, height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluPerspective(60.0, (float)width / height, 0.01, 1000);
   glMatrixMode(GL_MODELVIEW);
}



void display()	// display function
{
  
	// Clear the depth buffer and the color buffer
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	
	//Add ambient light
	GLfloat ambientColor[] = {0.2f, 0.2f, 0.2f, 1.0f}; //Color (0.2, 0.2, 0.2)
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);
	
	//Add positioned light
	GLfloat lightColor1[] = {0.2f, 0.2f, 0.7f, 1.0f}; //Color (0.2, 0.2, 0.7)
	GLfloat lightPos1[] = {20.0f, 20.0f, 30.0f, 1.0f}; //Positioned at (20,20,30)
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor1);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos1);
	
	//Add directed light
	GLfloat lightColor2[] = {0.2f, 0.2f, 0.2f, 1.0f}; //Color (0.2, 0.2, 0.2)
	//Coming from the direction (0.5, 0.5, 0.5)
	GLfloat lightPos2[] = {0.5f, 0.5f, 0.5f, 0.0f};
	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor2);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPos2);
	

	
	

	for(int i=0;i<g_iNumOfSplines;i++)
	{
		drawSpline(i);
	}

	
	glEnable(GL_TEXTURE_2D);
	
	//top face
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,GL_REPLACE);
	glBindTexture(GL_TEXTURE_2D, cubeTexIdUp);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, 0.0f);
	glVertex3f(-50.0f, -50.0f, 95.0f);
	glTexCoord2f(1.0f, 0.0f);
	glVertex3f(50.0f, -50.0f, 95.0f);
	glTexCoord2f(1.0f, 1.0f);
	glVertex3f(50.0f, 50.0f, 95.0f);
	glTexCoord2f(0.0f, 1.0f);
	glVertex3f(-50.0f, 50.0f, 95.0f);

	glEnd();

	//side1 face
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,GL_REPLACE);
	glBindTexture(GL_TEXTURE_2D, cubeTexIdSide1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, 0.0f);
	glVertex3f(-50.0f, -50.0f, 95.0f);
	glTexCoord2f(1.0f, 0.0f);
	glVertex3f(50.0f, -50.0f, 95.0f);
	glTexCoord2f(1.0f, 1.0f);
	glVertex3f(50.0f, -50.0f, -5.0f);
	glTexCoord2f(0.0f, 1.0f);
	glVertex3f(-50.0f, -50.0f, -5.0f);
	glEnd();

	//side2 face
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,GL_REPLACE);
	glBindTexture(GL_TEXTURE_2D, cubeTexIdSide2);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, 0.0f);
	glVertex3f(-50.0f, 50.0f, 95.0f);
	glTexCoord2f(1.0f, 0.0f);
	glVertex3f(50.0f, 50.0f, 95.0f);
	glTexCoord2f(1.0f, 1.0f);
	glVertex3f(50.0f, 50.0f, -5.0f);
	glTexCoord2f(0.0f, 1.0f);
	glVertex3f(-50.0f, 50.0f, -5.0f);
	glEnd();
	
	//side3 face
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,GL_REPLACE);
	glBindTexture(GL_TEXTURE_2D, cubeTexIdSide3);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, 0.0f);
	glVertex3f(-50.0f, -50.0f, 95.0f);
	glTexCoord2f(1.0f, 0.0f);
	glVertex3f(-50.0f, 50.0f, 95.0f);
	glTexCoord2f(1.0f, 1.0f);
	glVertex3f(-50.0f, 50.0f, -5.0f);
	glTexCoord2f(0.0f, 1.0f);
	glVertex3f(-50.0f, -50.0f, -5.0f);
	//glEnd();


	//side4 face
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,GL_REPLACE);
	//glBindTexture(GL_TEXTURE_2D, cubeTexIdSide4);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	
	glBegin(GL_QUADS);
	glTexCoord2f(1.0f, 1.0f);
	glVertex3f(50.0f, 50.0f, -5.0f);
	glTexCoord2f(0.0f, 1.0f);
	glVertex3f(50.0f, -50.0f, -5.0f);
	glTexCoord2f(0.0f, 0.0f);
	glVertex3f(50.0f, -50.0f, 95.0f);
	glTexCoord2f(1.0f, 0.0f);
	glVertex3f(50.0f, 50.0f, 95.0f);
	glEnd();
	
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,GL_REPLACE);
	glBindTexture(GL_TEXTURE_2D, groundTexId);
	
	//Ground
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);


	glBegin(GL_QUADS);
	glColor3f(1.0,1.0,1.0);
	glNormal3f(0.0, 1.0f, -1.0f);
	glTexCoord2f(1.0f, 0.0f);
	glVertex3f(-50.0f, -50.0f, -5.0f);
	glTexCoord2f(1.0f, 1.0f);
	glVertex3f(-50.0f, 50.0f, -5.0f);
	glTexCoord2f(0.0f, 1.0f);
	glVertex3f(50.0f, 50.0f, -5.0f);
	glTexCoord2f(0.0f, 0.0f);
	glVertex3f(50.0f, -50.0f, -5.0f);
	glEnd();

	

	glDisable(GL_TEXTURE_2D);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(currentLookAt.x+0.02*normal.x,currentLookAt.y+0.02*normal.y,currentLookAt.z+0.02*normal.z,currentLookAt.x+tangent.x,currentLookAt.y+tangent.y,currentLookAt.z+tangent.z,normal.x,normal.y,normal.z);

		  
	// Use double buffering
	glutSwapBuffers();
}


int _tmain(int argc, _TCHAR* argv[])
{
	// I've set the argv[1] to track.txt.
	// To change it, on the "Solution Explorer",
	// right click "assign1", choose "Properties",
	// go to "Configuration Properties", click "Debugging",
	// then type your track file name for the "Command Arguments"
	if (argc<2)
	{  
		printf ("usage: %s <trackfile>\n", argv[0]);
		exit(0);
	}

	loadSplines(argv[1]);
	
	loadSplinePoints(controlPointOfSegment,0,0);

	glutInit(&argc,(char**)argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
	
	// Create a window of size 640x480
	glutInitWindowSize(640,480);
	glutInitWindowPosition(200, 100);
	glutCreateWindow("CR-spline");
	
	/* tells glut to use a particular display function to redraw */
	glutDisplayFunc(display);

	// Reshape function call
	glutReshapeFunc (reshape);

  
	/* allow the user to quit using the right mouse button menu */
	g_iMenuId = glutCreateMenu(menufunc);
	glutSetMenu(g_iMenuId);
	glutAddMenuEntry("Quit",0);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

  
	/* replace with any animate code */
	//glutIdleFunc(doIdle);
	glutTimerFunc(66, doIdle, 0);

	/* callback for keyboard */
	glutKeyboardFunc(keyboard);

	
	
	/* do initialization */
	myinit();
	
	glutMainLoop();


	
	return 0;
}
