#include "bitmap.h"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h> 
#include <errno.h> 
#include <string.h> 
#include <sys/types.h> 
#include <sys/wait.h> 
#include <unistd.h> 

int iteration_to_color( int i, int max ); 
int iterations_at_point( double x, double y, int max ); 
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max ); 
void make_mandel(char*, float); 
void make_filename_array(void); 
void make_scale_array(void); 

int processes; 
int out_loop; 
pid_t pid; 
char filenames[50][100]; 
double scales[50]; 
int arr_index; 

double xcenter = -0.0001;
double ycenter = 0.637010;
double scale = .000002;
int    image_width = 900;
int    image_height = 900;
int    max = 2000;


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
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
} 

int main( int argc, char *argv[] )
{	
	int remainder;			//stores the remainder of 50/processes 
	if(argc <= 1)
	{
		printf("No arguments.\n"); 
		exit(1); 
	}
	processes = atoi(argv[1]);
	remainder = 50%processes;  
	arr_index = 0; 	
 	int i; 				//general variable for various loops	


	make_filename_array();		//create array of names for mandel1.bmp through mandel50.bmp 
	make_scale_array(); 		//make array of scales to be used within required range
	int j = 0;
	int is_remainder; 		//bool-like variable. Indicates whether there's a remainder needing to be dealt with.
	if(remainder > 0)		//If there's one, is_remainder = 1
	{
		is_remainder = 1; 
	}
	else
	{
		is_remainder = 0;       //If there's no remainder, is_remainder = 0
	}
	int remain = 0;			//Will set this when it's time to incorporate the remainder
	while(j < (50/processes + is_remainder))
	{
		if(j == 50/processes && ((j*processes + i + remain) < 50)) //if j has reached max and it's time to 
		{							   //incorporate the remainder in the loop, 	
			remain++; 					   //set the remainder to 1.	
		}
		i = 0; 
		while( i < processes && ((j*processes + i + remain) < 50)) //loop through inner loop the amount of times
		{							   //requested by the user. This loop creates
			pid = fork(); 					   //the processes. 	 
			if( pid == -1 )
			{
				perror("fork failed: "); 
				exit( EXIT_FAILURE ); 	
			}
			else if ( pid == 0 )
			{	//if we're in the child process, make_mandel with the current filenames[] index 
				//and current scales[] index.
				make_mandel(filenames[j*processes + i + remain], scales[j*processes + i + remain]); 
				exit( EXIT_SUCCESS );  
			}	
			else
			{	
				int status; 
				waitpid(pid, &status, 0);
			}
			i++;
 		}
 		j++;
	}
	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max )
{
	int i,j;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

	// For every pixel in the image...

	for(j=0;j<height;j++) {

		for(i=0;i<width;i++) {

			// Determine the point in x,y space for that pixel.
			double x = xmin + i*(xmax-xmin)/width;
			double y = ymin + j*(ymax-ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,max);

			// Set the pixel in the bitmap.
			bitmap_set(bm,i,j,iters);
		}
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
	int magenta = 255*i/max;
	return MAKE_RGBA(magenta,0,magenta,0);
}
/*
 *Function: make_filename_array
 *Parameters: void
 *Returns: void
 *Description: declares char name and assigns the appropriate 
 *mandelX.bmp filename to it. Then stores the filename in a 
 *global char* array. Increases X until we have reached mandel50.bmp
 */

void make_filename_array(void){
	int j = 0; 
	char name[100] = {0}; 

	while(j < 50)
	{
		sprintf(name, "mandel%d.bmp", j+1); 
		strcpy(filenames[j], name); 
		j++; 
	}
	
} 


void make_scale_array(void){
	double diff = 2.000000 - scale; 
	double interval = diff/49.0;
	int i = 1; 
	scales[0] = 2.00000; 
	while(i < 51)
	{
		scales[i] = scales[i-1] - interval; 
		i++; 
	}
}

void make_mandel(char* filename, float scale_val)
{
		// Create a bitmap of the appropriate size.
		struct bitmap *bm = bitmap_create(image_width,image_height);

		// Fill it with a dark blue, for debugging
		bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

		// Compute the Mandelbrot image
		compute_image(bm,xcenter-scale_val,xcenter+scale_val,ycenter-scale_val,ycenter+scale_val,max);

		// Save the image in the stated file.
		if(!bitmap_save(bm,filename)) {
			fprintf(stderr,"mandel: couldn't write to %s: %s\n",filename,strerror(errno));
		}
}
 

