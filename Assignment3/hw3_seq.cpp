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

typedef struct body{
  double x_pos, y_pos;
  double x_v, y_v;
  double x_f, y_f; 
  double mass;
}body_t; 


void init_body(body_t *body, int body_num)
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
}


void calculation(body_t *body, int body_num)
{
    double x1, x2, y1, y2, m1, m2, xv1, xv2, yv1, yv2, r, f_x, f_y;

    for (int i = 0; i < body_num; i++)
    {
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


void move(body_t *body, int body_num)
{
    double m, x_a, y_a, x_v0, y_v0, x_d, y_d;

    for (int i = 0; i < body_num; i++)
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


void init_graph()
{
    display = XOpenDisplay (NULL);
    if (display == NULL)
    {
        fprintf (stderr, "Cannot connect to X server.\n");
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

	XSelectInput(display, win, ExposureMask|KeyPressMask);

    XMapWindow (display, win);

    gc = DefaultGC(display,screen);
	pm = XCreatePixmap(display,win,Xborder,Yborder,DefaultDepth(display,screen));
	XFillRectangle(display,pm,gc,0,0,Xborder,Yborder);
	XSetForeground(display,gc,WhitePixel(display,screen));
}


void NbodySequential(int argc, char **argv)
{
    struct timeval timeStart, timeEnd;
    double runTime;

    body_t *body = (body_t *)malloc(sizeof(body_t)*body_num);

    printf("Name: Cui Yuncong\nStudent ID: 118010045\nAssignment 3, N-body simulation, Sequential Implementation\n");
    gettimeofday(&timeStart, NULL );

    if (PIC)
        init_graph();

    init_body(body, body_num);

    for (int k = 0; k < Iteration; k++)
    {
        if (PIC)
        {
            XSetForeground(display,gc,0);
            XFillRectangle(display,pm,gc,0,0,Xborder,Yborder);
        }

        calculation(body, body_num);
        move(body, body_num);

        if (PIC && k == 0)
            XFlush(display);
        
        if (PIC)
        {
            XSetForeground(display,gc,WhitePixel (display, screen));
            for (int i = 0; i < body_num; i++)
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
}


int main(int argc, char *argv[])
{
    NbodySequential(argc, argv);
    return 0;
}
