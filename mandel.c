/// 
//  mandel.c
//  Based on example code found here:
//  https://users.cs.fiu.edu/~cpoellab/teaching/cop4610_fall22/project3.html
//
//  Converted to use jpg instead of BMP and other minor changes
//  
///
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>
#include "jpegrw.h"
#include <time.h>
#include <semaphore.h>
#include <pthread.h>
#include <fcntl.h>

#define NUM_FRAMES 50

// local routines
static int iteration_to_color( int i, int max );
static int iterations_at_point( double x, double y, int max );
static void compute_image( imgRawImage *img, double xmin, double xmax,
									double ymin, double ymax, int max );
static void show_help();

// Worker process to generate a portion of the image
void child_process(int id, double x_center, double y_center, double zoom, int num_frames, int image_width, int image_height, int max) {
    char filename[256];
    sprintf(filename, "mandel%d.jpg", (id+1));  // Assuming JPEG output
    printf("Child %d started processing\n", (id+1));

    // Create a raw image of the appropriate size.
    imgRawImage *img = initRawImage(image_width, image_height);

    // Fill it with a black color (assuming black is the default)
    setImageCOLOR(img, 0);

    // Compute the Mandelbrot image
    compute_image(img, x_center - zoom / 2, x_center + zoom / 2, y_center - zoom / 2, y_center + zoom / 2, max);

    // Save the image in the stated file
    storeJpegImageFile(img, filename);

    // Free the allocated memory
    freeRawImage(img);

    printf("Child %d finished processing\n", (id+1));
}

int main( int argc, char *argv[] )
{
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.
	const char *outfile = "mandel.jpg";
	double xcenter = 0;
	double ycenter = 0;
	double xscale = 4;
	double yscale = 0; // calc later
	int    image_width = 1000;
	int    image_height = 1000;
	int    max = 1000;
	int num_children = 1;

	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:o:n:h"))!=-1) {
		switch(c) 
		{
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				xscale = atof(optarg);
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
			case 'o':
				outfile = optarg;
				break;
			case 'n':
				num_children = atoi(optarg);
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

	// Calculate y scale based on x scale (settable) and image sizes in X and Y (settable)
	yscale = xscale / image_width * image_height;

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf xscale=%lf yscale=%1f max=%d outfile=%s\n",xcenter,ycenter,xscale,yscale,max,outfile);

	// Create a raw image of the appropriate size.
	imgRawImage* img = initRawImage(image_width,image_height);

	// Fill it with a black
	setImageCOLOR(img,0);

	// Compute the Mandelbrot image
	compute_image(img,xcenter-xscale/2,xcenter+xscale/2,ycenter-yscale/2,ycenter+yscale/2,max);

	// Save the image in the stated file.
	storeJpegImageFile(img,outfile);

	// free the mallocs
	freeRawImage(img);

	clock_t start_time = clock(); // Start the timer

    pid_t pid;

	int count = 0;

    // number of frames processed each round
    for (int frame = 0; frame < (NUM_FRAMES/num_children); frame++) {
        for (int i = 0; i < num_children; i++) {
            pid = fork();  // Fork a child process for each frame

            if (pid == 0) {
                // Child process: generate the specific frame assigned to it
                double frame_zoom = (1.0 + (count+i) * 0.025);  // Adjust the zoom per child
                child_process(i+count, xcenter, ycenter, frame_zoom, NUM_FRAMES, image_width, image_height, max);
                exit(0);  // Exit after completing the task for this frame
            }
        }
		count += num_children;
        // Wait for all child processes to finish before moving to the next frames
        for (int i = 0; i < num_children; i++) {
            wait(NULL);
        }
	}
	if(NUM_FRAMES % num_children != 0) {
		int extra = NUM_FRAMES % num_children;
		for (int i = 0; i < extra; i++) {
            pid = fork();  // Fork a child process for each frame

            if (pid == 0) {
                // Child process: generate the specific frame assigned to it
                double frame_zoom = (1.0 + (count+i) * 0.025);  // Adjust the zoom per child
                child_process(i+count, xcenter, ycenter, frame_zoom, NUM_FRAMES, image_width, image_height, max);
                exit(0);  // Exit after completing the task for this frame
            }
        }
        // Wait for all child processes to finish before moving to the next frames
        for (int i = 0; i < num_children; i++) {
            wait(NULL);
        }
    }

    clock_t end_time = clock(); // End the timer
    double runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Total runtime with %d children and %d frames: %.2f seconds\n", num_children, NUM_FRAMES, runtime);
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

	return iter;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void compute_image(imgRawImage* img, double xmin, double xmax, double ymin, double ymax, int max )
{
	int i,j;

	int width = img->width;
	int height = img->height;

	// For every pixel in the image...

	for(j=0;j<height;j++) {

		for(i=0;i<width;i++) {

			// Determine the point in x,y space for that pixel.
			double x = xmin + i*(xmax-xmin)/width;
			double y = ymin + j*(ymax-ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,max);

			// Set the pixel in the bitmap.
			setPixelCOLOR(img,i,j,iteration_to_color(iters,max));
		}
	}
}


/*
Convert a iteration number to a color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/
int iteration_to_color( int iters, int max )
{
	int color = 0xFFFFFF*iters/(double)max;
	return color;
}


// Show help message
void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates (X-axis). (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=1000)\n");
	printf("-H <pixels> Height of the image in pixels. (default=1000)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}
