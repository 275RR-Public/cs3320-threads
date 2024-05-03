/*
  Name: James Hofer
  ID: 1000199225
  Omega Compile: make all (added -std=c99 to make file)
*/

#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin,
					double ymax, int max, int num_threads );

struct thread_params
{
    struct bitmap *bm;
    double xmin, xmax, ymin, ymax;
    int max, start_height, end_height;
};


void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-n <max>    The maximum number of threads. (default=1)\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}


int main( int argc, char *argv[] )
{
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.
	const char *outfile = "mandelseries.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int    max = 1000;
	int    num_threads = 1;

	// For each command line argument given,
	// override the appropriate configuration value.
	while((c = getopt(argc,argv,"x:y:s:W:H:m:n:o:h"))!=-1) {
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'n':
				num_threads = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

	// Start stopwatch
	struct timeval begin_time;
  	struct timeval end_time;
	gettimeofday( &begin_time, NULL );
	
	// Display the configuration of the image.
	printf("\nmandelseries: x=%lf y=%lf scale=%lf max=%d threads=%d outfile=%s\n",
			xcenter,ycenter,scale,max,num_threads,outfile);

	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
	compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max,num_threads);

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandelseries: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}

	// Stop stopwatch and print time elapsed
	gettimeofday( &end_time, NULL );
	long time_to_execute = ( end_time.tv_sec * 1000000 + end_time.tv_usec ) -
                         ( begin_time.tv_sec * 1000000 + begin_time.tv_usec );
	printf("This code took *** %ld microseconds *** to execute\n\n", time_to_execute);

	return 0;
}


/*
Each thread computes a horizontal slice of the image
*/
void* compute_image_thread(void* arg)
{
    struct thread_params* data = (struct thread_params*)arg;

    int width = bitmap_width(data->bm);
    int height = bitmap_height(data->bm);

    // For every pixel in the image...
	for(int j = data->start_height; j < data->end_height; j++) {

        for(int i = 0; i < width; i++) {

			// Determine the point in x,y space for that pixel.
            double x = data->xmin + i * (data->xmax - data->xmin) / width;
            double y = data->ymin + j * (data->ymax - data->ymin) / height;

			// Compute the iterations at that point.
            int iters = iterations_at_point(x, y, data->max);

			// Set the pixel in the bitmap.
            bitmap_set(data->bm, i, j, iters);
        }
    }

    return NULL;
}


/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max".
*/
void compute_image(struct bitmap *bm, double xmin, double xmax, double ymin,
					double ymax, int max, int num_threads)
{
    // Determine the size of a slice based on the number of threads
	int height = bitmap_height(bm);
    int height_per_thread = height / num_threads;
    int remainder = height % num_threads;

    pthread_t tid[num_threads];
    struct thread_params data[num_threads];

    // Set the params for each thread
	for(int i = 0; i < num_threads; i++) {
        data[i].bm = bm;
        data[i].xmin = xmin;
        data[i].xmax = xmax;
        data[i].ymin = ymin;
        data[i].ymax = ymax;
        data[i].max = max;
        
		// Set the start and end height for each thread
		data[i].start_height = i * height_per_thread;
		
		// Adjust the end height for the last thread since it may have a remainder
		if(i == num_threads - 1) {
			data[i].end_height =  i * height_per_thread + height_per_thread + remainder;
		}
		else {
			data[i].end_height = (i + 1) * height_per_thread;
		}

        pthread_create(&tid[i], NULL, compute_image_thread, &data[i]);
    }

    for(int i = 0; i < num_threads; i++) {
        pthread_join(tid[i], NULL);
    }
}


/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/
int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}


/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/
int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	//int red = i % 255;
	//int blue = 255 - i % 255;
	return MAKE_RGBA(gray,gray,gray,0);
}
