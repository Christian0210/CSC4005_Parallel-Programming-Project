#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <cmath>
#include <sys/time.h>
#include <omp.h>

pthread_mutex_t mutex;

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

typedef struct body{
  double x_pos, y_pos;
  double x_v, y_v;
  double x_f, y_f;
  double mass;
}body_t;

struct Local{
    int body_num;
    int local_body_num;
    int NUM_THREADS;
    int threadId;
    body_t *body;
};

Window win;
int width, height,
    x, y,
    border_width,
    display_width, display_height,
    screen;

GC gc;
XGCValues values;
Display *display;
XSizeHints size_hints;
Pixmap pm;
long valuemask=0;


void init_graph()
{
    display = XOpenDisplay (NULL);
    if (display == NULL)
    {
        fprintf (stderr, "Cannot connect to X server. \n");
        exit (-1);
    }

    screen = DefaultScreen (display);
    display_width = DisplayWidth (display, screen);
    display_height = DisplayHeight (display, screen);

    width = Xborder;
    height = Yborder;
    x = 10;
    y = 10;

    border_width = 1;
    win = XCreateSimpleWindow (display, RootWindow (display, screen),
                            x, y, width, height, border_width,
                            BlackPixel (display, screen), WhitePixel (display, screen));

    XSelectInput(display,win,ExposureMask|KeyPressMask);
    XMapWindow (display, win);

    gc = DefaultGC(display,screen);
    pm = XCreatePixmap(display,win,Xborder,Yborder,DefaultDepth(display,screen));
    XFillRectangle(display,pm,gc,0,0,Xborder,Yborder);
    XSetForeground(display,gc,WhitePixel(display,screen));

}


void init_body(body_t *body, int body_num)
{
    for (int i = 0; i < body_num; i++)
    {
        srand(time(0) + rand());
        body[i].mass = rand() % 1000 + 3000;
        body[i].x_v = 0;
        body[i].y_v = 0;
        body[i].x_f = 0;
        body[i].y_f = 0;
        body[i].x_pos = rand() % (maxX - minX) + minX;
        body[i].y_pos = rand() % (maxY - minY) + minY;

    }
}


void move(body_t *body, int localStart, int localEnd)
{
    double m, x_a, y_a, x_v0, y_v0, x_d, y_d;

    for (int i = localStart; i < localEnd; i++)
    {
        m = body[i].mass;
        x_a = body[i].x_f / m;
        y_a = body[i].y_f / m;
        x_v0 = body[i].x_v;
        y_v0= body[i].y_v;

        x_d = x_v0 * delta_t + 0.5 * x_a * pow(delta_t, 2);
        y_d = y_v0 * delta_t + 0.5 * y_a * pow(delta_t, 2);

        body[i].x_pos += x_d;
        body[i].y_pos += y_d;

        if (body[i].x_pos < 0 || body[i].x_pos > Xborder)
            body[i].x_v = -body[i].x_v;

        if (body[i].y_pos < 0 || body[i].y_pos > Yborder)
            body[i].y_v = -body[i].y_v;
    }
}


void *local_update(void *idp)
{
    Local *local = (Local*)idp;
    int localStart = local->threadId * local->local_body_num;
    int localEnd = (local->threadId + 1) * local->local_body_num;
    double x1, x2, y1, y2, m1, m2, r, f_x, f_y, xv1, xv2, yv1, yv2;

    if (local->threadId == local->NUM_THREADS)
        localEnd = local->body_num;

    for (int i = localStart; i < localEnd; i++)
    {
        local->body[i].x_f = 0;
        local->body[i].y_f = 0;
        for (int j = 0; j < local->body_num; j++)
        {
            x1 = local->body[i].x_pos;
            x2 = local->body[j].x_pos;
            y1 = local->body[i].y_pos;
            y2 = local->body[j].y_pos;
            m1 = local->body[i].mass;
            m2 = local->body[j].mass;

            if (x1 != x2 || y1 != y2) 
            {
                r = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
                f_x = G * m1 * m2 * (x2 - x1) / pow(r, 3);
                f_y = G * m1 * m2 * (y2 - y1) / pow(r, 3);

                if (r > 0.5)
                {
                    local->body[i].x_f += f_x;
                    local->body[i].y_f += f_y;
                }
                else 
                {
                    xv1 = local->body[i].x_v;
                    xv2 = local->body[j].x_v;
                    yv1 = local->body[i].y_v;
                    yv2 = local->body[j].y_v; 
                    local->body[i].x_v = ((m1 - m2) / (m1 + m2)) * xv1 + (2 * m2 / (m1 + m2)) * xv2;
                    local->body[j].x_v = (2 * m1 / (m1 + m2)) * xv1 + ((m2 - m1) / (m1 + m2)) * xv2;
                    local->body[i].y_v = ((m1 - m2) / (m1 + m2)) * yv1 + (2 * m2 / (m1 + m2)) * yv2;
                    local->body[j].y_v = (2 * m1 / (m1 + m2)) * yv1 + ((m2 - m1) / (m1 + m2)) * yv2;
                }
            }
        }
    }
    move(local->body,localStart,localEnd);
}


void NbodyPthread(int argc, char **argv)
{
    struct timeval timeStart, timeEnd;
    double runTime=0;
    int i, k, rc;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    int NUM_THREADS = atoi(argv[1]);

    pthread_t threads[NUM_THREADS];
    body_t *body = (body_t *)malloc(sizeof(body_t)*body_num);

    int local_body_num = body_num / NUM_THREADS;
    if (body_num % NUM_THREADS != 0)
        local_body_num++;

    printf("Name: Cui Yuncong\nStudent ID: 118010045\nAssignment 3, N-body simulation, Pthread Implementation\n");
    gettimeofday(&timeStart, NULL );

    if(PIC)
        init_graph();

    init_body(body, body_num);

    pthread_mutex_init(&mutex, NULL);

    Local *local = (Local *)malloc(sizeof(Local) * NUM_THREADS);

    for (i = 0; i < NUM_THREADS; i++)
    {
        local[i].body_num = body_num;
        local[i].local_body_num = local_body_num;
        local[i].NUM_THREADS = NUM_THREADS;
        local[i].threadId = i;
        local[i].body = body;
    }

    for (k = 0; k < Iteration; k++)
    {
        if (PIC)
        {
            XSetForeground(display, gc,0);
            XFillRectangle(display,pm,gc,0,0,Xborder,Yborder);
        }

        for (i = 0; i <NUM_THREADS; i++)
            rc = pthread_create(&threads[i], NULL,local_update, &local[i]);
        for (i = 0; i < NUM_THREADS; i++)
            pthread_join(threads[i], NULL);

        if (PIC)
        {   
            XFlush(display);
            XSetForeground(display,gc,WhitePixel (display, screen));
            for (i = 0; i < body_num; i++)
                XDrawPoint (display, pm, gc, (int)body[i].y_pos, (int)body[i].x_pos);
            XCopyArea(display,pm,win,gc,0,0,Xborder,Yborder,0,0);
        }

    }

    if (PIC)
    {
        XFreePixmap(display,pm);
        XCloseDisplay(display);
    }

    gettimeofday(&timeEnd, NULL);
    runTime = (timeEnd.tv_sec - timeStart.tv_sec ) + (double)(timeEnd.tv_usec -timeStart.tv_usec)/1000000;
    printf("runTime is  %lf\n", runTime);

    pthread_mutex_destroy(&mutex);
    pthread_exit(NULL);
}


int main(int argc, char *argv[])
{
    NbodyPthread(argc, argv);
    return 0;
}

                                                                


                                                                    






                                                       
