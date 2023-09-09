#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdbool.h>
#include <mpi.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <omp.h>

#define MASTER 0
#define X_RESN 200
#define Y_RESN 200
#define ITER 5000
#define legal(x,n) ((x)>=0 && (x) < (n))
#define EPSILON 1

bool  PIC = false;

int dx[4] = {0,-1,0,1};
int dy[4] = {1,0,-1,0};

int total = X_RESN * Y_RESN;
int iteration = ITER;
int x = X_RESN;
int y = Y_RESN;


Window          win;                            /* initialization for a window */
unsigned int    width, height,                  /* window size */
			    border_width,                   /* border width in pixels */
			    idth, display_height,  		    /* size of screen */
			    screen;                         /* which screen */
GC              gc;
unsigned long   valuemask = 0;
XGCValues       values;
Display         *display;
XSizeHints      size_hints;
Pixmap          bitmap;
FILE            *fp, *fopen();	
Colormap		default_cmap;
XColor		    color[256];


int temperatue_to_color_pixel(double t)
{
    return color[(int)(t/5.0)].pixel;
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
    width = X_RESN;
    height = Y_RESN;

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


void XRedraw(double *map)
{
    int i, j;
    #pragma omp for schedule(static)
    for (i=0; i<X_RESN; i++)
        for (j=0; j<Y_RESN; j++)
        {
            XSetForeground(display, gc, temperatue_to_color_pixel(map[i * X_RESN + j]));
            XDrawPoint (display, win, gc, j, i);
        }
    XFlush (display);
}


void temperature_iterate(double *sbuf, double *map, int my_rank, int BUFFER_SIZE)
{
    int i,j,d;
    int cur = my_rank * BUFFER_SIZE;
    int cur_end = (my_rank+1)*BUFFER_SIZE;
    #pragma omp for schedule(static)
    for (i = cur; i < cur_end; i++)
    {
        for (j = 0; j <X_RESN; j++)
        {
            int cnt = 0;
            sbuf[(i-cur) * X_RESN + j] = 0;
            for (d = 0; d < 4; d++)
            {
                if (legal(i+dx[d], X_RESN) && legal(j+dy[d], X_RESN))
                {
                    sbuf[(i-cur)*X_RESN + j] += map[(i+dx[d])*X_RESN + (j+dy[d])];
                    cnt ++;
                }
            }
            sbuf[(i-cur)*X_RESN + j] /= cnt;
        }
    }
}


void HeatSimulationHybrid(int argc, char **argv)
{
    int comm_sz,my_rank,len;
    int iter;
    double *sbuf;
    double *map;
    double *tempmap;
    double start_time;
    double cur_time;
    double end_time;
    int i,j;
    int BUFFER_SIZE,sub_total;
    char hostname[MPI_MAX_PROCESSOR_NAME];

    int NUM_THREADS = 4;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Get_processor_name(hostname, &len);
    MPI_Status status;

    BUFFER_SIZE = X_RESN / comm_sz;

    sub_total = BUFFER_SIZE * X_RESN; //one for map, one for temp

    sbuf = (double *)malloc(sub_total*sizeof(double));
    map = (double *)malloc(total*sizeof(double));
    tempmap = (double *)malloc(total*sizeof(double));

    if(my_rank == 0)
    {
        if (PIC)
            Xwindow_Init();
        start_time = MPI_Wtime();
    }

    MPI_Bcast(&map[0], total, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (i = 0; i < X_RESN; i++)
    {
        int cur = i * X_RESN;
        tempmap[cur] = 20.0f;
        tempmap[cur+X_RESN-1] = 20.0;
        tempmap[i] = 20.f;
        tempmap[(X_RESN-1)*Y_RESN+i] = 20.0;
    }

    for (iter = 0; iter < iteration; iter++)
    {
        omp_set_num_threads(NUM_THREADS);

        int cur = my_rank * BUFFER_SIZE;
        temperature_iterate(sbuf, map, my_rank, BUFFER_SIZE);
        MPI_Gather(sbuf, sub_total, MPI_DOUBLE, tempmap, sub_total, MPI_DOUBLE,0, MPI_COMM_WORLD);
        if (my_rank == 0)
        {            
            for (j = Y_RESN/2-Y_RESN/12; j < Y_RESN/2+Y_RESN/12; j++)
                tempmap[j] = 100.0;
            map = tempmap;
            if (PIC)
                XRedraw(map);
            for (i = 1; i < comm_sz; i++)
            {
                if (i != comm_sz-1)
                    MPI_Send(&map[(i*BUFFER_SIZE-1)*X_RESN],(BUFFER_SIZE+2)*X_RESN,MPI_DOUBLE,i, 1, MPI_COMM_WORLD);
                else
                    MPI_Send(&map[(i*BUFFER_SIZE-1)*X_RESN],(BUFFER_SIZE+1)*X_RESN,MPI_DOUBLE,i, 1, MPI_COMM_WORLD);
            }
        }

        if (my_rank != 0 && my_rank != comm_sz -1)
            MPI_Recv(&map[(my_rank*BUFFER_SIZE-1)*X_RESN],(BUFFER_SIZE+2)*X_RESN,MPI_DOUBLE,0, 1, MPI_COMM_WORLD, &status);

        if (my_rank == comm_sz - 1)
            MPI_Recv(&map[(my_rank*BUFFER_SIZE-1)*X_RESN],(BUFFER_SIZE+1)*X_RESN,MPI_DOUBLE,0, 1, MPI_COMM_WORLD, &status);   
    } 

    end_time = MPI_Wtime();
    if (my_rank == 0)
    {
        printf("Name: Cui Yuncong\n");
        printf("Student ID: 118010045\n");
        printf("Assignment 4, Heat Simulation, MPI_OpenMP Implementation.\n");
        printf("runTime is = %fs\n",end_time-start_time);
    }
    MPI_Finalize();
}


int main(int argc, char *argv[])
{
    HeatSimulationHybrid(argc, argv);
    return 0;
}









