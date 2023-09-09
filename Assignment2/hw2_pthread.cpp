#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdbool.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <pthread.h>

#define NUM_THREADS 3
double real_min;
double real_max;
double imag_min;
double imag_max;
int width;
int height;
bool PIC = true;
pthread_mutex_t mutex;
pthread_cond_t mutex_threshold;

typedef struct complextype
{
    float real, imag;
}Compl;

struct Local
{
    double xStart, xEnd;
    double yStart, yEnd;
    int width, height;
    int* output;
    int threadId;
    int numThreads;
};

GC gc;
Display *display;
Window window;
int screen;

void initGraph(int width,int height)
{
    display = XOpenDisplay(NULL);
    if(display == NULL) 
    {
        fprintf(stderr, "cannot open display\n");
        exit(1);
    }

    screen = DefaultScreen(display);

    int x = 0, y = 0, border_width = 0;

    window = XCreateSimpleWindow(display, RootWindow(display, screen), x, y, width, height, border_width, BlackPixel(display, screen), WhitePixel(display, screen));

    XGCValues values;
    long valuemask = 0;

    gc = XCreateGC(display, window, valuemask, &values);
    XSetForeground (display, gc, BlackPixel (display, screen));
    XSetBackground(display, gc, 0X0000FF00);
    XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);

    XMapWindow(display, window);
    XSync(display, 0);
}

int cal_pixel(Compl c)
{
    Compl z;
    int count, max;
    double temp, lengthsq;
    count = 0;
    max = 100;

    z.real = 0.0;
    z.imag = 0.0;

    do{
        temp = z.real * z.real - z.imag * z.imag + c.real;
        z.imag = 2 * z.real * z.imag + c.imag;
        z.real = temp;
        lengthsq = z.real * z.real + z.imag * z.imag;
        count ++;
    }while((lengthsq < 4.0) && (count < max));

    return count;
}

void *localExecution(void *idp)
{
    Local *local = (Local*)idp;
    Compl z, c;
    int repeats;
    double scale_real = (local->xEnd - local->xStart ) / local->width;
    double scale_imag = (local->yEnd - local->yStart) / local->height;
    double x, y;

    int local_width = local->width / local->numThreads;
    if (local->width % local->numThreads != 0)
        local_width = local_width + 1;

    int start = local_width * local->threadId;
    int end;
    if ((local_width * local->threadId + local_width) > local->width)
        end = local->width-1;
    else
        end = local_width * local->threadId + local_width-1;

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < local->height; j++)
        {
            x= local->xStart + i * scale_real;
            y= local->yStart + j * scale_imag;
            z.real = 0.0;
            z.imag = 0.0;
            c.real= x;
            c.imag = y;
            repeats = cal_pixel(c);
            int index = (i * local->height + j);
            local->output[index] = repeats;
       }
    }
}

void MandelbrotPthread(double real_min, double real_max, double imag_min,
                       double imag_max, int width, int height, int *output)
{
    int rc;
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    pthread_t threads[NUM_THREADS];
    struct timeval timeStart, timeEnd, timeSystemStart;
    double runTime=0, systemRunTime;
    gettimeofday(&timeStart, NULL );

    Local *local = (Local *)malloc(sizeof(Local) * NUM_THREADS);

    for (int i = 0; i < NUM_THREADS; i++)
    {
        local[i].xStart = real_min;
        local[i].xEnd = real_max;
        local[i].yStart = imag_min;
        local[i].yEnd = imag_max;
        local[i].width = width;
        local[i].height = height;
        local[i].output = output;
        local[i].threadId = i;
        local[i].numThreads = NUM_THREADS;
    }

    for (int i = 0; i < NUM_THREADS; i++)
        rc = pthread_create(&threads[i], NULL, localExecution, &local[i]);

    for (int i = 0;i < NUM_THREADS; i++)
        pthread_join(threads[i], NULL);
   
    gettimeofday( &timeEnd, NULL );
    runTime = (timeEnd.tv_sec - timeStart.tv_sec ) + (double)(timeEnd.tv_usec
               - timeStart.tv_usec) / 1000000;
    printf("runingTime: %lf\n", runTime);

    pthread_attr_destroy(&attr);
}

int main(int argc, char *argv[])
{
    double real_min = -2;
    double imag_min = -2;
    double real_max = 2;
    double imag_max = 2;
    int width = 200, height = 200;

    printf("Name: Cui Yuncong\nStudent ID: 118010045\nAssignment 2, Mandelbrot Set, Pthread implementation\n");
    int *results = (int *)malloc(sizeof(int) * width * height);

    if(PIC)
        initGraph(width,height);

    MandelbrotPthread(real_min, real_max, imag_min, imag_max, width, height,results);

    if(PIC)
    {
        for(int i = 0; i < width; i++)
            for(int j = 0; j <height; j++)
            {
                XSetForeground (display, gc,  1024 * 1024 * (results[i*height + j] % 256));
                XDrawPoint (display, window, gc, i, j);
            }
        XFlush(display);
        sleep(5);
    }

    return 0;
}
 