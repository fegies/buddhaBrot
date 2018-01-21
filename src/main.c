#define _POSIX_C_SOURCE 2

#include "qdbmp.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <stdlib.h>

void fillPipe(unsigned int pointCount, int pipe);
void readPipe(unsigned int **points, unsigned int xsize, unsigned int ysize, int pipe);
void worker(int input, int output, unsigned int xsize, unsigned int ysize, unsigned int maxIterations);

void print_usage() {
	printf("Usage: buddhaBrot -o outputfile [-i iterations] [-c pointcount] [-x xsize] [-y ysize] [-p workerCount]\n");
	printf("Instead of using -i for iterations you can specify -n to start nebulabrot mode\n");
	printf("In Nebulabrot Mode, specify [-r iterations red channel] [-g iterations green channel] [-b iterations blue channel]\n");
}


void buddhabrot(unsigned int ** pointData,unsigned int pointCount, unsigned int xsize, unsigned int ysize, unsigned int iterations, unsigned int workercount) {
	//pipe between the writer and the worker processes
	int questionpipe[2];
	pipe(questionpipe);

	//start the writer
	int writeProcess = fork();
	if(writeProcess == 0) {
		fillPipe(pointCount,questionpipe[1]);
		close(questionpipe[1]);
		exit(0);
	}
	close(questionpipe[1]);

	//pipe between the worker process and this process
	int answerpipe[2];
	pipe(answerpipe);

	//spawn the worker processes
	for( int i = 0; i < workercount; ++i )
		if( fork() == 0 ) {
			worker(questionpipe[0], answerpipe[1], xsize, ysize, iterations);
			close(answerpipe[1]);
			exit(0);
		}
	close(answerpipe[1]);

	//get the image
	readPipe(pointData,xsize, ysize, answerpipe[0]);

	//collect child process zombies
	while( wait(0) > 0)
		;
}

int main(int argc,char *const *argv) {


	unsigned int xsize = 2000;
	unsigned int ysize = 2000;
	unsigned int iterations = 2000;
	unsigned int points = 20000000;
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

	if( !nebulaBrot ) {
		unsigned int **pointData = malloc(sizeof(unsigned int*) * xsize);
		for( unsigned int x = 0; x < xsize; ++x )
			pointData[x] = calloc(sizeof(unsigned int),ysize);

		printf("starting to render...\n");

		buddhabrot(pointData, points, xsize, ysize, iterations, workercount);

		printf("Processing Image...\n");

		unsigned int max = 0;
		for(unsigned int x = 0; x < xsize; ++x )
			for(unsigned int y = 0; y < ysize; ++y )
				if( pointData[x][y] > max )
					max = pointData[x][y];

		BMP *image = BMP_Create(xsize, ysize, 24);

		for(unsigned int x = 0; x < xsize; ++x ) {
			for(unsigned int y = 0; y < ysize; ++y ) {
				unsigned int v = pointData[x][y];
				unsigned char val = (unsigned char) (((float)v / max) * 255);
				BMP_SetPixelRGB(image, x,y, val,val,val);
			}
			free(pointData[x]);
		}
		free(pointData);

		BMP_WriteFile(image, outputfile);
		BMP_Free(image);
	}
	else {
		unsigned int **pointDataR = malloc(sizeof(unsigned int*) * xsize);
		unsigned int **pointDataG = malloc(sizeof(unsigned int*) * xsize);
		unsigned int **pointDataB = malloc(sizeof(unsigned int*) * xsize);
		for(unsigned int x = 0; x < xsize ; ++x ) {
			pointDataR[x] = calloc(sizeof(unsigned int),ysize);
			pointDataG[x] = calloc(sizeof(unsigned int),ysize);
			pointDataB[x] = calloc(sizeof(unsigned int),ysize);
		}

		printf("rendering red channel...\n");
		buddhabrot(pointDataR,points, xsize, ysize, iterationsR, workercount);
		printf("rendering green channel...\n");
		buddhabrot(pointDataG,points, xsize, ysize, iterationsG, workercount);
		printf("rendering blue channel...\n");
		buddhabrot(pointDataB,points,xsize,ysize,iterationsB,workercount);

		printf("gathering image data...\n");

		unsigned int maxR = 0, maxG = 0, maxB = 0;
		for( unsigned int x = 0; x < xsize; ++x )
			for(unsigned int y = 0; y < ysize; ++y ) {
				if( maxR < pointDataR[x][y])
					maxR = pointDataR[x][y];
				if( maxG < pointDataG[x][y])
					maxG = pointDataG[x][y];
				if( maxB < pointDataB[x][y])
					maxB = pointDataB[x][y];
			}

		BMP *image = BMP_Create(xsize, ysize, 24);

		for(unsigned int x = 0; x < xsize; ++x ) {
			for(unsigned int y = 0; y < ysize; ++y ) {
				unsigned char r = (unsigned char) (((float)(pointDataR[x][y]) / maxR) * 255);
				unsigned char g = (unsigned char) (((float)(pointDataG[x][y]) / maxG) * 255);
				unsigned char b = (unsigned char) (((float)(pointDataB[x][y]) / maxB) * 255);

				BMP_SetPixelRGB(image, x,y, r,g,b);
			}
			free(pointDataR[x]);
			free(pointDataG[x]);
			free(pointDataB[x]);
		}
		free(pointDataR);
		free(pointDataG);
		free(pointDataB);

		BMP_WriteFile(image, outputfile);
		BMP_Free(image);
	}

	return 0;
}

struct coordinate {
	unsigned int x;
	unsigned int y;
};
struct complex {
	double real;
	double imaginary;
};

#define BUFSIZE 255
void readPipe( unsigned int** points, unsigned int xsize, unsigned int ysize, int pipe ) {
	struct coordinate buffer[BUFSIZE];
	ssize_t readsize = 0;
	while ((readsize = read(pipe, buffer, BUFSIZE*sizeof(struct coordinate))) > 0) {
		if( readsize % sizeof(struct coordinate) != 0 ) {
			fprintf(stderr, "read a number of bytes that does not match answer struct size!");
			exit(1);
		}
		int itemsread = readsize / sizeof(struct coordinate);
		for(unsigned int i = 0; i < itemsread; ++i) {
			struct coordinate c = buffer[i];
			points[c.x][c.y]++;
		}
	}
}

void fillPipe(unsigned int pointCount, int pipe) {
	struct complex buffer[BUFSIZE];
	unsigned int bufpos = 0;
	srand(0); // because I want deterministic results
	unsigned int onePercent = pointCount / 100;

	printf("\n");
	for (unsigned int i = 0; i < pointCount; ++i) {
		if( i % onePercent == 0 ) {
			unsigned int percentage = i / onePercent;
			printf("\r\033[1A%d%%\n", percentage);
			printf("[%*c%*c",percentage+1,'#',100-percentage,']');
			fflush(stdout);
		}
		double x = ((double)rand() / RAND_MAX) * 4 - 2;
		double y = ((double)rand() / RAND_MAX) * 4 - 2;

		//check that the point is not within the primary or secondary bulb of the set
		// (if it is, it is in the Mandelbrot set, so it will not be rendered anyway)
		if( (x+1)*(x+1) + y*y < 1.0/16) {
			--i;
			continue;
		}
		double q = (x-0.25)*(x-0.25) + y*y;
		if( q*(q+(x-0.25)) < 0.25*y*y ) {
			--i;
			continue;
		}

		buffer[bufpos++] = (struct complex) { x,y };
		if( bufpos == BUFSIZE ) {
			write(pipe,buffer,bufpos*sizeof(struct complex));
			bufpos = 0;
		}
	}
	if( bufpos > 0)
		write(pipe,buffer,bufpos*sizeof(struct complex));
	printf("\r\033[1A100%%\n[%*c%*c\n", 100, '#', 0, ']');
}

void toCoordinate( struct coordinate* coord, double r, double i, unsigned int xsize, unsigned int ysize) {
	coord->y = (unsigned int) ((r + 2) * ((xsize-1) / 4));
	coord->x = (unsigned int) ((i + 2) * ((ysize-1) / 4));
}

void bufferadd( struct coordinate* buffer, int pipe, unsigned int* bufpos, struct coordinate value ) {
	buffer[(*bufpos)++] = value;
	if( *bufpos == BUFSIZE ) {
		write( pipe, buffer, BUFSIZE*sizeof(struct coordinate));
		*bufpos = 0;
	}
}

void worker(int input, int output, unsigned int xsize, unsigned int ysize, unsigned int maxIterations ) {
	struct complex inputBuffer[BUFSIZE];
	struct coordinate answerBuffer[BUFSIZE];
	ssize_t readsize = 0;
	unsigned int bufpos = 0;
	while ( (readsize = read(input,inputBuffer, sizeof(struct complex)*BUFSIZE)) > 0 ) {
		if( readsize % sizeof(struct complex) != 0) {
			printf("read a number of bytes that does not match struct size!\n");
			exit(1);
		}
		int itemsRead = readsize / sizeof(struct complex);
		for(int i = 0; i < itemsRead; ++i) {
			struct complex q = inputBuffer[i];
			double cr = q.real;
			double ci = q.imaginary;
			double r = 0;
			double i = 0;

			for( unsigned int iteration = 0; iteration < maxIterations; ++iteration ) {
				double s = (r*r - i*i) + cr; // Z_{n+1} = Z_n^2 + C
				double j = (2*r*i)   + ci;
				r = s;
				i = j;

				if((r*r + i*i) > 4.0) { // the point escapes so it has to be traced
					r = 0;
					i = 0;
					for( ;; ) {
						s = (r*r - i*i) + cr;
						j = (2*r*i)   + ci;
						r = s;
						i = j;

						double abs = r*r + i*i;
						if(abs > 4.0)
							break;
						struct coordinate val;
						toCoordinate(&val, r,i,xsize, ysize);
						bufferadd(answerBuffer, output, &bufpos, val);
					}
					break;
				}
			}
		}
	}
	if( bufpos > 0)
		write(output, answerBuffer, bufpos*sizeof(struct coordinate));
}