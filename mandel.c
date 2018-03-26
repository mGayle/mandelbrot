/*Fixed segfault. seems to run well. :) 
 *
 */

#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>


struct Parameters{
	int num_threads; 
	int xmin; 
	int xmax; 
	int ymin; 
	int ymax; 
	float scale; 
}; 

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void *compute_image( void *paramsarg ); 
int get_num_lines( int image_height, int n ); 
struct Parameters* make_struct_array( struct Parameters* my_struct_ptr, int num_threads, 
			int image_width, int image_height, float scale, int num_lines );

struct Parameters my_struct[1000];
struct Parameters* my_struct_ptr = my_struct;  
struct bitmap *bm; 

float xcenter; 
float ycenter;
float scale; 
int image_width; 
int image_height;  
int max; 
int num_threads; 
const char *outfile;

	

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-n 	    Number of threads. (default=1)\n"); 
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");

	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

int main( int argc, char *argv[] )
{	
	//DON'T FORGET TO DELETE THIS WHEN YOU'RE READY TO TURN IT IN!!!!!!!
	FILE *f = fopen("mandel_vals.txt", "a"); 
	if (f == NULL)
	{
		printf("Error opening file!\n"); 
		exit(1); 
	}
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	outfile = "mandel.bmp";
	xcenter = 0;
	ycenter = 0;
	scale = 4;
	image_width = 500;
	image_height = 500;
	max = 1000;
	num_threads = 1; 

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

	int num_lines = get_num_lines(image_height, num_threads); 

	pthread_t tid[1000] = {0};  

	my_struct_ptr = make_struct_array(my_struct_ptr, num_threads, image_width, 
						image_height, scale, num_lines);   	
	struct timeval begin; 
	struct timeval end; 
	long int time_diff; 

		

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d threads=%d, outfile=%s\n"
		,xcenter,ycenter,scale,max,num_threads,outfile);

	// Create a bitmap of the appropriate size.
	bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));
			
	gettimeofday( &begin, NULL ); 
	// Compute the Mandelbrot image
	int m; 
	for( m = 0; m < num_threads; m++ )
	{
		if(pthread_create( &tid[m], NULL, compute_image,&my_struct[m]))
		{
			perror("Error creating thread: "); 
			exit( EXIT_FAILURE );  
		}
		if(pthread_join( tid[m], NULL))
		{
			perror("Problem with pthread_join: "); 
		}
	}

	gettimeofday( &end, NULL ); 
//	printf("\nBegin: %ld   End: %ld\n", (begin.tv_sec * 1000000 + begin.tv_usec),
//                                          (end.tv_sec * 1000000 + end.tv_usec)); 
	time_diff = ( end.tv_sec * 1000000 + end.tv_usec ) - ( begin.tv_sec * 1000000 + begin.tv_usec );
	printf("\nTime difference: %ld\n", time_diff); 	
	fprintf(f, "-x %f -y %f -s %f -W %d -H %d -m %d -n %d\n%ld\n\n", xcenter, ycenter, scale, 
						image_width, image_height, max, num_threads, time_diff); 
	fclose(f); 	
		
	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}

	return 0;
}

int get_num_lines(int image_height, int num_threads)
{
	return image_height/num_threads; 
}


struct Parameters* make_struct_array(struct Parameters* my_struct_ptr, int num_threads, int image_width, int image_height, float scale, int num_lines)
{

	my_struct_ptr[0].xmin = 0; 
	my_struct_ptr[0].xmax = image_width; 
	my_struct_ptr[0].ymin = 0; 
	my_struct_ptr[0].ymax = num_lines; 
	my_struct_ptr[0].scale = scale; 
	my_struct_ptr[0].num_threads = num_threads; 

	int remainder_lines; 
	remainder_lines = 500 % num_lines; 
	printf("REMAINDER LINES: %d\n", remainder_lines); 
	int i = 1; 
	while (i < num_threads)
	{
		my_struct_ptr[i].xmin = 0; 
		my_struct_ptr[i].xmax = image_width; 
		my_struct_ptr[i].ymin = my_struct_ptr[i-1].ymax; 
		my_struct_ptr[i].ymax = my_struct_ptr[i].ymin + num_lines;
		my_struct_ptr[i].scale = scale; 
		my_struct_ptr[i].scale = scale; 
		my_struct_ptr[i].num_threads = num_threads; 
		printf("in make struct function: struct[%d] y min: %d\n",i, my_struct_ptr[i].ymin);
		printf("in make struct function: struct[%d] y max: %d\n",i, my_struct_ptr[i].ymax);  
		i++; 
	}
	my_struct_ptr[i-1].ymax = my_struct_ptr[i-1].ymax + remainder_lines; 
	printf("in make struct function: struct[%d] y min: %d\n",i-1, my_struct_ptr[i-1].ymin);
	printf("in make struct function: struct[%d] y max: %d\n",i-1, my_struct_ptr[i-1].ymax);  
	return my_struct_ptr; 
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/


void *compute_image( void *paramsarg )
{
	struct Parameters *this_struct; 
	this_struct = (struct Parameters *) paramsarg; 	

	int i,j; 
	
	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

	double xmin = xcenter - scale; 
	double xmax = xcenter + scale; 
	double ymin = ycenter - scale; 
	double ymax = ycenter + scale; 

	// For every pixel in the image...
	

	for(j = this_struct->ymin; j < this_struct->ymax; j++) {
		for(i = this_struct->xmin; i < this_struct->xmax; i++) {
			// Determine the point in x,y space for that pixel.
			double x = xmin + i*(xmax-xmin)/width;
			double y = ymin + j*(ymax-ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,max);

			// Set the pixel in the bitmap.
			bitmap_set(bm,i,j,iters);
		}
	}
	return NULL; 
}


/*
void compute_image( struct Parameters params[] )
{
	int i,j;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

	double xmin = xcenter - scale; 
	double xmax = xcenter + scale; 
	double ymin = ycenter - scale; 
	double ymax = ycenter + scale; 

//	int num_threads = params[0].num_threads; 
	
//	float scale = params[0].scale; 
//	double xmin = params[0].xcenter - scale; 
//	double xmax = params[0].xcenter + scale; 
//	double ymin = params[0].ycenter - scale; 
//	double ymax = params[0].ycenter + scale; 
//	int max = params[0].max; 

//	double xmin = -4; 
//	double xmax = 4; 
//	double ymin = -4; 
//	double ymax = 4; 
//	int max = 1000; 

	int idx = 0; 
	// For every pixel in the image...
	while(idx < num_threads)
	{
		for(j=0;j<params[idx].ymax;j++) {

			for(i=0;i<params[idx].xmax;i++) {

				// Determine the point in x,y space for that pixel.
				double x = xmin + i*(xmax-xmin)/width;
				double y = ymin + j*(ymax-ymin)/height;

				// Compute the iterations at that point.
				int iters = iterations_at_point(x,y,max);

				// Set the pixel in the bitmap.
				bitmap_set(bm,i,j,iters);
			}
		}
		idx++; 
	}
}
*/
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
	return MAKE_RGBA(gray,0,gray,0);
}



