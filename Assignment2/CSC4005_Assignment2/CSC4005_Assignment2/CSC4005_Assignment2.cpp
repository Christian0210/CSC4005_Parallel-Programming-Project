#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <mpi.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

#define  MASTER 0

double real_min;
double real_max;
double imag_min;
double imag_max;
int width;
int height;

int my_rank, comm_sz, len, local_width;
char hostname[MPI_MAX_PROCESSOR_NAME];
MPI_Status status;
MPI_Request request;
int tag;

typedef struct Compltype
{
    double real, imag;
}Compl;

struct Local
{
    int start;
    int end;
};

GC gc;
Display* display;
Window window;
int screen;

void initGraph(int width, int height)
{
    display = XOpenDisplay(NULL);
    if (display == NULL)
    {
        fprintf(stderr, "cannot open display\n");
        exit(1);
    }

    screen = DefaultScreen(display);

    int x = 0;
    int y = 0;

    int border_width = 0;

    window = XCreateSimpleWindow(display, RootWindow(display, screen), x, y, width, height, border_width, BlackPixel(display, screen), WhitePixel(display, screen));

    XGCValues values;
    long valuemask = 0;

    gc = XCreateGC(display, window, valuemask, &values);
    XSetForeground(display, gc, BlackPixel(display, screen));
    XSetBackground(display, gc, 0X0000FF00);
    XSetLineAttributes(display, gc, 1, LineSolid, CapRound, JoinRound);

    XMapWindow(display, window);
    XSync(display, 0);
}

void displayGraph(int x, int y, int repeats)
{
    if (isEnable)
    {
        XSetForeground(display, gc, 1024 * 1024 * (repeats % 256));
        XDrawPoint(display, window, gc, x, y);
    }
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

    do {                                         /* number of iterations */
        temp = z.real * z.real - z.imag * z.imag + c.real;
        z.imag = 2 * z.real * z.imag + c.imag;
        z.real = temp;
        lengthsq = z.real * z.real + z.imag * z.imag;
        count++;
    } while ((lengthsq < 4.0) && (count < max));

    return count;
}

void masterExecution(double real_min, double real_max, double imag_min, double imag_max, int width, int height)
{
    int receiveSize, index, receiveRank, realDataSize;
    Local* localX = (Local*)malloc((comm_sz - 1) * sizeof(Local));
    int* localInfo = (int*)malloc(sizeof(int) * 2);
    int* colors = (int*)malloc(width * height * sizeof(int));

    for (int i = 0; i < (comm_sz - 1); i++)
    {
        localX[i].start = i * local_width;
        if (((i * local_width) + local_width) > width)
            localX[i].end = width - 1;
        else
            localX[i].end = (i * local_width) + local_width - 1;
    }

    for (int i = 0; i < (comm_sz - 1); i++)
    {
        localInfo[0] = localX[i].start;
        localInfo[1] = localX[i].end;
        MPI_Isend(&localInfo[0], 2, MPI_INT, (i + 1), tag, MPI_COMM_WORLD, &request);
    }

    for (int i = 0; i < (comm_sz - 1); i++)
    {
        receiveSize = local_width * height;
        int* results = (int*)malloc(receiveSize * sizeof(int));

        MPI_Recv(results, receiveSize, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
        receiveRank = status.MPI_SOURCE;
        index = localX[receiveRank - 1].start * height;
        realDataSize = localX[receiveRank - 1].end - localX[receiveRank - 1].start + 1;
        memcpy(colors + index, results, (realDataSize * height * sizeof(int)));

        free(results);
        results = NULL;
        if (isEnable)
        {
            for (int k = localX[receiveRank - 1].start; k < localX[receiveRank - 1].end + 1; k++)
            {
                index = k * height;
                for (int j = 0; j < height; j++)
                    displayGraph(k, j, colors[index + j]);
            }
            XFlush(display);
        }
    }
    if (isEnable)
        sleep(10);
}


void  calculation(double real_min, double real_max, double imag_min, double imag_max,
    int start, int end, int width, int height, int* results)
{
    Compl z, c;
    int repeats;
    double scale_real = (real_max - real_min) / width;
    double scale_imag = (imag_max - imag_min) / height;
    double x, y;

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < height; j++)
        {
            x = real_min + i * scale_real;
            y = imag_min + j * scale_imag;
            z.real = 0.0;
            z.imag = 0.0;
            c.real = x;
            c.imag = y;
            repeats = cal_pixel(c);
            results[((i - start) * height) + j] = repeats;
        }
    }
}

void slaveExecution(double real_min, double real_max, double imag_min, double imag_max, int width, int height)
{
    int* localInfo = (int*)malloc(sizeof(int) * 2);
    MPI_Recv(&localInfo[0], 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    int matrixSize = local_width * height;
    int* colors = (int*)malloc(matrixSize * sizeof(int));
    calculation(real_min, real_max, imag_min, imag_max, localInfo[0], localInfo[1], width, height, colors);
    MPI_Isend(colors, matrixSize, MPI_INT, MASTER, tag, MPI_COMM_WORLD, &request);
}

int MandelbrotMPI(int argc, char** argv)
{
    double real_min = -2;
    double imag_min = -2;
    double real_max = 2;
    double imag_max = 2;
    int width = 200, height = 200;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Get_processor_name(hostname, &len);
    struct timeval timeStart, timeEnd, timeSystemStart;
    double runTime = 0, systemRunTime;

    local_width = width / (comm_sz - 1);

    if (width % (comm_sz - 1) != 0)
        local_width = local_width + 1;

    if (my_rank == MASTER)
    {
        printf("Name: Cui Yuncong\n");
        printf("Student ID: 118010045\n");
        printf("Assignment 2, Mandelbrot Set, MPI Implementation\n");
        gettimeofday(&timeStart, NULL);
        if (isEnable)
            initGraph(width, height);

        masterExecution(real_min, real_max, imag_min, imag_max, width, height);
        gettimeofday(&timeEnd, NULL);
        runTime = (timeEnd.tv_sec - timeStart.tv_sec) + (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000;
        printf("runTime is  %lf\n", runTime);
    }
    else
        slaveExecution(real_min, real_max, imag_min, imag_max, width, height);

    MPI_Finalize();
}

int main(int argc, char* argv[])
{
    MandelbrotMPI(argc, argv);
    return 0;
}

