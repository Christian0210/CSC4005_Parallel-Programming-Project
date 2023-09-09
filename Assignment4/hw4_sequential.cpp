#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include "time.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <math.h>


#define RES_X 200
#define RES_Y 200
#define ITER 5000
#define legal(x,n) ((x)>=0 && (x) < (n))
bool PIC = false;

int total = RES_X * RES_Y;
int iteration = ITER;
int x = RES_X;
int y = RES_Y;
int dx[4] = {0, -1, 0, 1};
int dy[4] = {1, 0, -1, 0};

double *map;
double *tempmap;


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
FILE           *fp, *fopen();	
Colormap		default_cmap;
XColor		color[256];


int temperatue_to_color_pixel(double t)
{
	return color[(int)(t/5.0f)].pixel;
}


void Xwindow_Init()
{
	XSetWindowAttributes attr[1];

	/* connect to Xserver */
	display = XOpenDisplay (NULL);
	if (display == NULL) {
		fprintf (stderr, "Cannot connect to X server \n");
		exit (-1);
	}

	/* get screen size */
	screen = DefaultScreen (display);

	/* set window size */
	width = RES_X;
	height = RES_Y;

	/* set window position */
	x = 10;
	y = 10;

	/* create opaque window */
	border_width = 4;
	win = XCreateSimpleWindow (display, RootWindow (display, screen),
						x, y, width, height, border_width,
						BlackPixel (display, screen), WhitePixel (display, screen));
	size_hints.flags = USPosition|USSize;
	size_hints.x = 0;
	size_hints.y = 0;
	size_hints.width = width;
	size_hints.height = height;
	size_hints.min_width = 300;
	size_hints.min_height = 300;
		
	XSetNormalHints (display, win, &size_hints);

	//XStoreName(display, win, window_name);

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

	/* create color */
	int i;
	for (i=0; i<20; i++)
	{
	    double w = double(i)/20;
	    color[i].green = 66*64+150*64*(1-w);
	    color[i].red = 65535*w;
	    color[i].blue = 65535*(1-w);
	    color[i].flags = DoRed | DoGreen | DoBlue;
	    XAllocColor(display, default_cmap, &color[i]);
	}
}


void XRedraw(double *map)
{
	int i, j;
	for (i=0; i<RES_X; ++i)
		for (j=0; j<RES_Y; ++j)
		{
			XSetForeground(display, gc, temperatue_to_color_pixel(map[i * RES_X + j]));
			XDrawPoint (display, win, gc, j, i);
		}
	XFlush (display);
}


void temperatue_iterate(double *map)
{
	int i, j, d;
	for (i = 0; i < RES_X; i++)
	{
		for (j = 0; j < RES_X; j++)
		{
			int cnt = 0;
			tempmap[i * RES_X + j] = 0;
			for (d = 0; d < 4; d++)
			{
				if (legal(i+dx[d], RES_X) && legal(j+dy[d], RES_X))
				{
					tempmap[i*RES_X + j] += map[(i+dx[d])*RES_X + (j+dy[d])];
					cnt ++;
				}
			}
			tempmap[i*RES_X + j] /= cnt;
		}
	}
	
	for (j = RES_Y/2-RES_Y/12; j < RES_Y/2+RES_Y/12; j++)
		tempmap[j] = 100.0;
}


void HeatSimulationSequential(int argc, char **argv)
{
	clock_t start, endtime;
	float time_sequential = 0;
	int i,j;
	int BUFFER_SIZE, sub_total;
	start = clock();
	map = (double *)malloc(total*sizeof(double));
	tempmap = (double *)malloc(total*sizeof(double));

	if (PIC)
		Xwindow_Init();
	for (i = 0; i < RES_X; i++)
	{
		int cur = i*RES_X;
		tempmap[cur] = 20.0f;
		tempmap[cur + RES_X - 1] = 20.0;
		tempmap[i] = 20.0f;
		tempmap[(RES_X-1)*RES_Y+i] = 20.0; 
	}
	int iter;
	for (iter = 0; iter < iteration; iter++)
	{
		temperatue_iterate(map);
		map = tempmap;
		if (PIC)
			XRedraw(map);
	}

	endtime = clock();
	time_sequential += (float)(endtime - start) / CLOCKS_PER_SEC;

	printf("Name: Cui Yuncong\n");
	printf("Student ID: 118010045\n");
	printf("Assignment 4, Heat Simulation, Sequential Implementation.\n");
	printf("runTime is %fs\n", time_sequential);
}


int main(int argc, char *argv[])
{
    HeatSimulationSequential(argc, argv);
    return 0;

}
