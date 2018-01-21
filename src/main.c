#define _POSIX_C_SOURCE 2

#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "qdbmp.h"
#include "xorshift.h"

void print_usage()
{
	printf("Usage: buddhaBrot -o outputfile [-i iterations] [-c pointcount] [-x xsize] [-y ysize] [-p workerCount]\n");
	printf("Instead of using -i for iterations you can specify -n to start nebulabrot mode\n");
	printf("In Nebulabrot Mode, specify [-r iterations red channel] [-g iterations green channel] [-b iterations blue channel]\n");
}


/**
 * gets a random Point between -2-2i and 2+2i
 * Checks for the two bulbs of the mandelbrot set and does not return those points
 * @param r
 * @param i
 */
void randomPoint( double *r, double* i, uint64_t rand_state[2] )
{
	for( ;; ) {
		double x = ((double)xorshift128plus(rand_state) / XORSHIFT_MAX) * 4 - 2;
		double y = ((double)xorshift128plus(rand_state) / XORSHIFT_MAX) * 4 - 2;

		//check that the point is not within the primary or secondary bulb of the set
		if( (x+1)*(x+1) + y*y < 1.0/16)
			continue;

		double q = (x-0.25)*(x-0.25) + y*y;
		if( q*(q+(x-0.25)) < 0.25*y*y )
			continue;

		*r = x;
		*i = y;
		return;
	}
}

/**
 * Calculates the image coordinates for the complex coordinate
 * @param x
 * @param y
 * @param r real component
 * @param i imaginary component
 * @param xsize X dimension of the image
 * @param ysize Y dimension of the image
 */
void toCoordinate( unsigned int *x, unsigned int *y, double r, double i, unsigned int xsize, unsigned int ysize)
{
	*y = (unsigned int) ((r + 2) * ((ysize-1) / 4));
	*x = (unsigned int) ((i + 2) * ((xsize-1) / 4));
}

/**
 * Computes the next step in the mandelbrot sequence
 * Z_n = Z_{n-1}^2 + C
 *
 * r,i : Current number (Z_{n-1})
 *
 * cr, ci : C
 */
void nextPoint( double *r, double *i, double cr, double ci)
{
	double s = *r;
	double j = *i;

	*r = (s*s - j*j) + cr;
	*i = (2*s*j) + ci;
}

/**
 * Fills the given grid with a buddhabrot.
 * @param pointData
 * @param pointCount Number of points to trace
 * @param xsize
 * @param ysize
 * @param iterations Iteration Depth
 */
void buddhabrot(unsigned int ** pointData,unsigned int pointCount, unsigned int xsize, unsigned int ysize, unsigned int iterations)
{
	uint64_t rand_state[2] = {rand(),rand()};
	for( unsigned int i = 0; i < pointCount; ++i ) {
		double r = 0,i = 0;
		double cr,ci;
		randomPoint(&cr,&ci,rand_state);

		for(unsigned int iter = 0; iter < iterations; ++iter ) {
			nextPoint(&r, &i, cr, ci);

			if( r*r + i*i > 4 ) { // it escapes, so it will be traced
				r = i = 0;

				for( ;; ) {
					nextPoint(&r, &i, cr, ci);

					if( r*r + i*i > 4 )
						break;

					unsigned int x,y;
					toCoordinate(&x, &y, r, i, xsize, ysize);

					pointData[x][y]++;
				}
				break;
			}
		}
	}
}

struct buddhabrotArgument {
	unsigned int** pointData;
	unsigned int pointCount;
	unsigned int xsize;
	unsigned int ysize;
	unsigned int iterations;
};

void* buddhabrot_pthread( void* vp_arg )
{
	struct buddhabrotArgument *arg = vp_arg;
	buddhabrot(arg->pointData, arg->pointCount, arg->xsize, arg->ysize, arg->iterations);
	return 0;
}

unsigned char interpolateBrightness(unsigned int value, unsigned int max)
{
	//linear interpolation
//	return (unsigned char) (255 * ((double)value / max));

	//square interpolation
	double d = sqrt(value) / sqrt(max);

	//third root interpolation
//	double d = pow(value, 1.0/3) / pow(max,1.0/3);

	// square root of 4 interpolation
//	double d = sqrt(sqrt(value)) / sqrt(sqrt(max));

	// exponential interpolation log2
//	double d = log2(value) / log2(max);

	// exponential interpolation ln
//	double d = log(value) / log(max);

	return (unsigned char) (255 * d);
}

BMP* greyscaleImage( unsigned int pointCount, unsigned int xsize, unsigned int ysize, unsigned int iterations, unsigned int workercount )
{
	pthread_t workers[workercount];
	unsigned int ** data[workercount];
	struct buddhabrotArgument threadargs[workercount];

	unsigned int pointsPerWorker = pointCount / workercount;

	for( unsigned int i = 0; i < workercount; ++i ) {
		data[i] = malloc(sizeof(unsigned int*) * xsize);
		for(unsigned int x = 0; x < xsize; ++x)
			data[i][x] = calloc(ysize, sizeof(unsigned int));
		threadargs[i] = (struct buddhabrotArgument) {
				.pointData = data[i],
				.pointCount = pointsPerWorker,
				.xsize = xsize,
				.ysize = ysize,
				.iterations = iterations
		};
		pthread_create(workers+i, 0, buddhabrot_pthread, threadargs+i);
	}

	//join all data into canonicData and find the max in the grid
	unsigned int max = 0;
	pthread_join(workers[0], 0);
	unsigned int **canonicData = data[0];
	for( unsigned int i = 1; i < workercount; ++i ) {
		pthread_join(workers[i], 0);
		for( unsigned int x = 0; x < xsize; ++x ) {
			for(unsigned int y = 0; y < ysize; ++y ) {
				canonicData[x][y] += data[i][x][y];
				if( max < canonicData[x][y] )
					max = canonicData[x][y];
			}
			free(data[i][x]);
		}
		free(data[i]);
	}

	BMP* image = BMP_Create(xsize,ysize, 24);
	for(unsigned int x = 0; x < xsize; ++x ) {
		for(unsigned int y = 0; y < ysize; ++y ) {
			unsigned char val = interpolateBrightness(canonicData[x][y], max);
			BMP_SetPixelRGB(image, x,y, val, val, val);
		}
		free(canonicData[x]);
	}
	free(canonicData);

	return image;
}

BMP* colorImage( unsigned int pointCount, unsigned int xsize, unsigned int ysize, unsigned int iterationsR,
                 unsigned int iterationsG, unsigned iterationsB, unsigned int workercount )
{
	printf("rendering red\n");
	BMP* resultR = greyscaleImage(pointCount, xsize, ysize,iterationsR, workercount);
	printf("rendering green\n");
	BMP* resultG = greyscaleImage(pointCount, xsize, ysize,iterationsG, workercount);
	printf("rendering blue\n");
	BMP* resultB = greyscaleImage(pointCount, xsize, ysize,iterationsB, workercount);

	printf("combining\n");
	for(unsigned int x = 0; x < xsize; ++x )
		for(unsigned int y = 0; y < ysize; ++y ) {
			unsigned char r,g,b,t;
			BMP_GetPixelRGB(resultR, x, y, &r, &t, &t);
			BMP_GetPixelRGB(resultG, x, y, &t, &g, &t);
			BMP_GetPixelRGB(resultB, x, y, &t, &t, &b);
			BMP_SetPixelRGB(resultR, x, y, r, g, b);
		}
	BMP_Free(resultG);
	BMP_Free(resultB);

	return resultR;
}

int main(int argc,char *const *argv)
{
	unsigned int xsize = 2000;
	unsigned int ysize = 2000;
	unsigned int iterations = 2000;
	unsigned int points = 0xfffffff;
	unsigned int workercount = 2;
	char* outputfile = 0;
	char nebulaBrot = 0;
	unsigned int iterationsR = 0;
	unsigned int iterationsG = 0;
	unsigned int iterationsB = 0;


	//region Parse command line Arguments
	int c;
	while((c = getopt(argc, argv, "c:i:x:y:p:o:r:g:b:n")) != -1) {
		switch (c) {
			case 'c':
				points = atoi(optarg);
				break;
			case 'x':
				xsize = atoi(optarg);
				break;
			case 'y':
				ysize = atoi(optarg);
				break;
			case 'i':
				iterations = atoi(optarg);
				break;
			case 'p':
				workercount = atoi(optarg);
				break;
			case 'n':
				nebulaBrot = 1;
				break;
			case 'r':
				iterationsR = atoi(optarg);
				break;
			case 'g':
				iterationsG = atoi(optarg);
				break;
			case 'b':
				iterationsB =atoi(optarg);
				break;
			case 'o':
				outputfile = optarg;
				break;
			default:
				print_usage();
				exit(1);
		}
	}
	if( outputfile == 0 ) {
		print_usage();
		return 1;
	}
	//endregion

	BMP* image;
	if( nebulaBrot )
		image = colorImage(points, xsize, ysize, iterationsR, iterationsG, iterationsB, workercount);
	else
		image = greyscaleImage(points, xsize, ysize, iterations, workercount);

	if( strcmp(outputfile,"-") == 0 )
		BMP_WriteStdout(image);
	else
		BMP_WriteFile(image, outputfile);

	BMP_Free(image);

	return 0;
}
