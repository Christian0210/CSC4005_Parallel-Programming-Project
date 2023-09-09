#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include "time.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <math.h>
#include <pthread.h>

#define X_RESN 200
#define Y_RESN 200
#define iteration 5000
#define num_threads 4
#define legal(x,n) ((x)>=0 && (x) < (n))
bool PIC = false;
int x, y;



Window         win;                            /* initialization for a window */
unsigned int   width, height,                  /* window size */
			border_width,                   /* border width in pixels */
			idth, display_height,  		  /* size of screen */
			screen;                         /* which screen */

GC             gc;
unsigned long  valuemask = 0;
XGCValues      values;
Display        *display;
XSizeHints     size_hints;
Pixmap         bitmap;
FILE           *fp, *fopen ();	
Colormap		default_cmap;
XColor		color[256];


typedef struct TemperatureField
{
	int x, y;
	double **t;
	double *storage;
}TemperatureField;


TemperatureField *field;
TemperatureField *tempField, *swapField;
int dx[4] = {0,-1,0,1};
int dy[4] = {1,0,-1,0};


int count = 0;
int local_count = 0;


void *slave(void *id)
{
	int tid = *((int *)id);
	int subtask = field->y/num_threads;
	int i,j,d;
	for(i = tid * subtask; i < (tid + 1)*subtask; i++)
	{
		for(j = 0; j < field->y; j++)
		{
			int cnt = 0;
			tempField->t[i][j] = 0;
			for (d = 0; d < 4; ++d)
			{
				if (legal(i+dx[d], field->x) && legal(j+dy[d], field->y))
				{
					tempField->t[i][j] += field->t[i+dx[d]][j+dy[d]];
					++cnt;
				}
			}
			tempField->t[i][j] /= cnt;
		}
	}
}

void XWindow_Init(TemperatureField *field)
{    
     XSetWindowAttributes attr[1];       
       
     /* connect to Xserver */
	display = XOpenDisplay (NULL);
    	if (display == NULL)
	{
        	fprintf (stderr, "Cannot connect to X server \n");
        	exit (-1);
    	}
        
	/* get screen size */

	screen = DefaultScreen (display);

	/* set window size *///XFlush (display);

	width = field->y;
	height = field->x;

	/* create opaque window */

	border_width = 4;
	win = XCreateSimpleWindow (display, RootWindow (display, screen),
						width, height, width, height, border_width, 
						BlackPixel (display, screen), WhitePixel (display, screen));

	size_hints.flags = USPosition|USSize;
	size_hints.x = 0;
	size_hints.y = 0;
	size_hints.width = width;
	size_hints.height = height;
	size_hints.min_width = 300;
	size_hints.min_height = 300;
	
	XSetNormalHints (display, win, &size_hints);
	

	/* create graphics context */

	gc = XCreateGC (display, win, valuemask, &values);

	default_cmap = DefaultColormap(display, screen);
	XSetBackground (display, gc, WhitePixel (display, screen));
	XSetForeground (display, gc, BlackPixel (display, screen));
	XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);

	attr[0].backing_store = Always;
	attr[0].backing_planes = 1;
	attr[0].backing_pixel = BlackPixel(display, screen);

	XChangeWindowAttributes(display, win, CWBackingStore | CWBackingPlanes | CWBackingPixel, attr);

	XMapWindow (display, win);
	XSync(display, 0); 

	for (int i=0; i<20; ++i)
	{
	    double w = double(i)/20;
	    color[i].green = 66*64+150*64*(1-w);
	    color[i].red = 65535*w;
	    color[i].blue = 65535*(1-w);
	    color[i].flags = DoRed | DoGreen | DoBlue;
	    XAllocColor(display, default_cmap, &color[i]);
	}
}


int temperatue_to_color_pixel(double t)
{
	return color[(int)(t/5.0f)].pixel;
}


void XRedraw(TemperatureField *field)
{
    int i, j;
    for (i=0; i<field->x; ++i)
        for (j=0; j<field->y; ++j)
	{
		XSetForeground(display, gc, temperatue_to_color_pixel(field->t[i][j]));
	        XDrawPoint (display, win, gc, j, i);
	}
    XFlush (display);
}


void temperatue_iterate(TemperatureField* field)
{
	int i;
	int resu = field->x;
	pthread_t thread[num_threads];
	int indexes[num_threads];
	for (i = 0; i <num_threads;i++)
	{
		indexes[i] = i;
		pthread_create(&thread[i],NULL,slave,(void*) & (indexes[i]));
	}

	for (i = 0; i <num_threads;i++)
		pthread_join(thread[i], NULL);	

	for (i = Y_RESN/2-Y_RESN/12; i < Y_RESN/2+Y_RESN/12; i++)
          tempField->t[0][i] = 100.0f;
}


void deleteField(TemperatureField *field)
{
	free(field->t);
	free(field->storage);
}


void newField(TemperatureField *field, int x, int y, int sourceX, int sourceY)
{
	TemperatureField temp = *field;
	field->storage = (double *)malloc( sizeof(double) * x * y);
	field->t = (double **)malloc( sizeof(double*) * x);
	field->x = x;
	field->y = y;
	int i, j;
	for (i=0; i<x; ++i)
		field->t[i] = &field->storage[i*y];
	if (sourceX)
	{
		double scaleFactorX = (double)sourceX/x;
		double scaleFactorY = (double)sourceY/y;
		for (i=0; i<x; ++i)
			for (j=0; j<y; ++j)
				field->t[i][j] = temp.t[(int)(i*scaleFactorX)][(int)(j*scaleFactorY)];
		deleteField(&temp);
	}
	else memset(field->storage, 0, sizeof(double)*x*y);
}


void initField(TemperatureField *field)
{
	int i, j;
	for (i=0; i<field->x; ++i)
		for (j=0; j<field->y; ++j)
			field->t[i][j] = 20.0;
}


void HeatSimulationPthread(int argc, char **argv)
{
	x = X_RESN;
	y = Y_RESN;
	
	struct timeval timeStart, timeEnd;
	double runTime=0;

	gettimeofday(&timeStart, NULL );
	field = (TemperatureField *)malloc(sizeof(TemperatureField));
	tempField = (TemperatureField *)malloc(sizeof(TemperatureField));
	newField(field, x, y, 0, 0);
	newField(tempField,x,y,0,0);
	initField(field);
	int i ,j;
	for (i = 0; i < field->x; i++)
		for(j = 0; j < field->y; j++)
			field->t[i][j] = 0;
		
	if (PIC)
		XWindow_Init(field);
	
	int iter;
	for(iter = 0; iter < iteration; iter++)
	{
		temperatue_iterate(field);
		swapField = field;
		field = tempField;
		tempField = swapField;
		
		if(PIC)
			XRedraw(field);
	}
	gettimeofday( &timeEnd, NULL );
	runTime = (timeEnd.tv_sec - timeStart.tv_sec ) + (double)(timeEnd.tv_usec -timeStart.tv_usec)/1000000;
	printf("Name: Cui Yuncong\n");
	printf("Student ID: 118010045\n");
	printf("Assignment 4, Heat Simulation, Pthread Implementation.\n");
	printf("runTime is %lfs\n", runTime);	
}


int main(int argc, char *argv[])
{
    HeatSimulationPthread(argc, argv);
    return 0;
}
