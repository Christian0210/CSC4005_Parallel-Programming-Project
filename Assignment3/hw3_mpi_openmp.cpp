#include <mpi.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <time.h>
#include <omp.h>
using namespace std;

const double e = 0.5;
const int MASTER = 0;
const double loss = -0.8;
const double G = 6.67;
const int Xborder = 200;
const int Yborder = 200;
const double delta_t = 0.01;
const int Iteration = 1500;
const int body_num = 200;
const int minX = 0.125 * double(Xborder);     
const int maxX = 0.875 * double(Xborder);    
const int minY = 0.125 * double(Yborder);   
const int maxY = 0.875 * double(Yborder);

bool PIC = false;

Window          win;
Display         *display;
GC              gc;
unsigned long   valuemask = 0;
XGCValues       values;
XSizeHints      size_hints;
XSetWindowAttributes   attr[1];
int             width, height,                  /* window size */
                x, y,                           /* window position */
                border_width,                   /*border width in pixels */
                display_width, display_height,  /* size of screen */
                screen;                         /* which screen */

typedef struct Body{
    double x_pos, y_pos;
    double x_v, y_v;
    double x_f, y_f; 
    double mass;
} body_t;


body_t* initialization(body_t *body, int body_num)
{
    srand(time(0) + rand());
    for (int i = 0; i < body_num; i++)
    { 
        body[i].mass = rand() % 1000 + 3000;
        body[i].x_v = 0;
        body[i].y_v = 0;
        body[i].x_f = 0;
        body[i].y_f = 0;
        body[i].x_pos = rand() % (maxX - minX) + minX; 
        body[i].y_pos = rand() % (maxY - minY) + minY;
    }
    return body;
}


void calculation(body_t *body, int body_num, int localStart, int localEnd) 
{
    int j;
    #pragma omp parallel private(j)
    {
        #pragma omp for schedule(static)
        for (int i = localStart; i < localEnd; i++)
        {   
            double x1, x2, y1, y2, m1, m2, xv1, xv2, yv1, yv2, r, f_x, f_y;
            body[i].x_f = 0;
            body[i].y_f = 0;

            for (int j = 0; j < body_num; j++)
            {
                x1 = body[i].x_pos;
                x2 = body[j].x_pos;
                y1 = body[i].y_pos;
                y2 = body[j].y_pos;
                m1 = body[i].mass;
                m2 = body[j].mass;
                xv1 = body[i].x_v;
                xv2 = body[j].x_v;
                yv1 = body[i].y_v;
                yv2 = body[j].y_v;

                if (x1 != x2 || y1 != y2)
                {
                    r = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2)); 
                    f_x = G * m1 * m2 * (x2 - x1) / pow(r, 3);     
                    f_y = G * m1 * m2 * (y2 - y1) / pow(r, 3);   
                    
                    if (r > 0.5)
                    {
                        body[i].x_f += f_x;
                        body[i].y_f += f_y;
                    }
                    else
                    {
                        body[i].x_v = ((m1 - m2) / (m1 + m2)) * xv1 + (2 * m2 / (m1 + m2)) * xv2;
                        body[j].x_v = (2 * m1 / (m1 + m2)) * xv1 + ((m2 - m1) / (m1 + m2)) * xv2;
                        body[i].y_v = ((m1 - m2) / (m1 + m2)) * yv1 + (2 * m2 / (m1 + m2)) * yv2;
                        body[j].y_v = (2 * m1 / (m1 + m2)) * yv1 + ((m2 - m1) / (m1 + m2)) * yv2;
                    }
                } 
            }
        }
    }
}


void move(body_t *body, int localStart, int localEnd)
{
    #pragma omp for schedule(static)
    for (int i = localStart; i < localEnd; i++)
    {
        double m, x_a, y_a, x_v0, y_v0, x_d, y_d;
        m = body[i].mass;
        x_a = body[i].x_f / m;
        y_a = body[i].y_f / m;
        x_v0 = body[i].x_v;
        y_v0 = body[i].y_v;
        x_d = x_v0 * delta_t + 0.5*x_a*pow(delta_t, 2);
        y_d = y_v0 * delta_t + 0.5*y_a*pow(delta_t, 2);

        body[i].x_pos += x_d;
        body[i].y_pos += y_d;

        if (body[i].x_pos < 0 || body[i].x_pos > Xborder)
            body[i].x_v = -body[i].x_v;

        if (body[i].y_pos < 0 || body[i].y_pos > Yborder)
            body[i].y_v = -body[i].y_v;
    }
}


void collision(body_t *body, int body_num, int localStart, int localEnd)
{
    int j;
    #pragma omp parallel private(j)
    {
        #pragma omp for schedule(static)
        for (int i = localStart; i < localEnd; i++)
        {
            double dx, dy, r, m1, m2, xv1, xv2, yv1, yv2;
                
            for(int j = 0; j < body_num; j++)
            {
                if (j != i)
                {
                    dx = body[j].x_pos - body[i].x_pos;
                    dy = body[j].y_pos - body[i].y_pos;
                    r = sqrt(pow(dx, 2) + pow(dy, 2));
                    if (r < 0.5)
                    {
                        m1 = body[i].mass;
                        m2 = body[j].mass;
                        xv1 = body[i].x_v;
                        xv2 = body[j].x_v;
                        yv1 = body[i].y_v;
                        yv2 = body[j].y_v;                        
                        body[i].x_v = ((m1 - m2) / (m1 + m2)) * xv1 + (2 * m2 / (m1 + m2)) * xv2;
                        body[j].x_v = (2 * m1 / (m1 + m2)) * xv1 + ((m2 - m1) / (m1 + m2)) * xv2;
                        body[i].y_v = ((m1 - m2) / (m1 + m2)) * yv1 + (2 * m2 / (m1 + m2)) * yv2;
                        body[j].y_v = (2 * m1 / (m1 + m2)) * yv1 + ((m2 - m1) / (m1 + m2)) * yv2;
                    }
                }
            }
            if (body[i].x_pos < minX)
            { 
                if (body[i].x_v < 0)
                    body[i].x_v = loss * body[i].x_v;
                body[i].x_pos = minX;
            }
            if (body[i].x_pos > maxX)
            {  
                if (body[i].x_v > 0)
                    body[i].x_v = loss * body[i].x_v;
                body[i].x_pos = maxX;
            }
            if (body[i].y_pos < minY)
            {  
                if (body[i].y_v < 0)
                    body[i].y_v = loss * body[i].y_v;
                body[i].y_pos = minY;
            }
            if (body[i].y_pos > maxY)
            {  
                if (body[i].y_v > 0)
                    body[i].y_v = loss * body[i].y_v;
                body[i].y_pos = maxY;
            }
        }
    }
}


void draw_body_t(Display *display, Window win, GC gc, body_t *body, int body_num)
{
     XClearWindow(display, win);
     for (int i = 0; i < body_num; i++)
         XDrawPoint(display, win, gc, (int)(body[i].x_pos), (int)(body[i].y_pos));
     XFlush(display);
}


void NbodyMPIOpenMP(int argc, char **argv)
{
    int i, k, comm_sz, my_rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int NUM_THREADS = atoi(argv[1]);

    struct timeval timeStart, timeEnd;
    double runTime=0, systemRunTime;
    double start_time,end_time;
    MPI_Status status;
    MPI_Datatype MPI_body_t;
    MPI_Type_contiguous(7, MPI_DOUBLE, &MPI_body_t);
    MPI_Type_commit(&MPI_body_t);
    
    int local_body_num = body_num / comm_sz;
    int last_local_body_num = local_body_num + body_num % comm_sz;
    int localStart = my_rank * local_body_num;
    int localEnd = localStart + local_body_num;
    body_t *local_body = (body_t *)malloc(sizeof(body_t)*local_body_num);
    body_t *body = (body_t *)malloc(sizeof(body_t)*body_num);
    int *displacement = (int *)malloc(comm_sz * sizeof(int));
    int *send_count = (int *)malloc(comm_sz * sizeof(int));

    for (i = 0; i < comm_sz; i++)
        displacement[i] = i * local_body_num;
    for (i = 0; i < comm_sz - 1; i++)
        send_count[i] = local_body_num;
    send_count[comm_sz-1] = last_local_body_num; 
   
    if (my_rank == MASTER)
    {
        printf("Name: Cui Yuncong\nStudent ID: 118010045\nAssignment 3, N-body simulation, MPI-OpenMP Implementation\n");

        if (PIC)
        {
            display = XOpenDisplay (NULL);
            if (display == NULL)
            {
                fprintf (stderr, "Cannot connect to X server \n");
                exit (-1);
            }
            
            screen = DefaultScreen(display);
            display_width = DisplayWidth(display, screen);
            display_height = DisplayHeight(display, screen);

            width = Xborder;
            height = Yborder;
            x = 0;
            y = 0;

            border_width = 4;
            win = XCreateSimpleWindow(display, RootWindow(display, screen),
                                    x, y, width, height, border_width,
                                    WhitePixel(display, screen), BlackPixel(display, screen));

            size_hints.flags = USPosition|USSize;
            size_hints.x = x;
            size_hints.y = y;
            size_hints.width = width;
            size_hints.height = height;
            size_hints.min_width = 300;
            size_hints.min_height = 300;

            XSetNormalHints (display, win, &size_hints);

            gc = XCreateGC(display, win, valuemask, &values);
            XSetBackground(display, gc, BlackPixel(display, screen));
            XSetForeground(display, gc, WhitePixel(display, screen));
            XSetLineAttributes(display, gc, 1, LineSolid, CapRound, JoinRound);

            attr[0].backing_store = Always;
            attr[0].backing_planes = 1;
            attr[0].backing_pixel = WhitePixel(display, screen);

            XChangeWindowAttributes(display, win, CWBackingStore | CWBackingPlanes | CWBackingPixel, attr);

            XMapWindow (display, win);
            XSync(display, 0);
        
            XFlush(display);
        }

        body = initialization(body, body_num);

        for (i = 1; i < comm_sz; i++)
            MPI_Send(body, body_num, MPI_body_t, i, i, MPI_COMM_WORLD);
    }

    omp_set_num_threads(NUM_THREADS);

    if (my_rank == comm_sz - 1)
    {
        local_body_num = last_local_body_num;
        localEnd = body_num;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    start_time = MPI_Wtime();
    for (k = 0; k < Iteration; k++)
    {
        if (my_rank != MASTER)
        {
            MPI_Status status;
            MPI_Recv(body, body_num, MPI_body_t, MASTER, my_rank, MPI_COMM_WORLD, &status);
        }
        calculation(body, body_num, localStart, localEnd);
        move(body, localStart, localEnd);
        collision(body, body_num, localStart, localEnd);

        for (i = 0; i < local_body_num; i++)
            local_body[i] = body[my_rank * (body_num/comm_sz) + i];

        MPI_Gatherv(local_body, send_count[my_rank], MPI_body_t, body, send_count, displacement, MPI_body_t, MASTER, MPI_COMM_WORLD);

        if (my_rank == MASTER)
        {
            if (PIC)
                draw_body_t(display, win, gc, body, body_num);
            for (i = 1; i < comm_sz; i++)
                MPI_Send(body, body_num, MPI_body_t, i, i, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    runTime = end_time - start_time;

    if (my_rank == MASTER)
        printf("runTime is  %lf\n", runTime);

    MPI_Finalize();
}


int main(int argc, char *argv[])
{
    NbodyMPIOpenMP(argc, argv);
    return 0;
}
