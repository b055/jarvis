/**
*****************************************************************************
**
**  File        : main.c
**
**  Abstract    : main and interrupt functions.
**
**  Functions   : main
**                Undefined_Handler
**                SWI_Handler
**                PAbort_Handler
**                DAbort_Handler
**                IRQ_Handler
**                FIQ_Handler
**
**  Environment : Atollic TrueSTUDIO(R)
**
**  Distribution: The file is distributed “as is,” without any warranty
**                of any kind.
**
**  (c)Copyright Atollic AB.
**  You may use this file as-is or modify it according to the needs of your
**  project. Distribution of this file (unmodified or modified) is not
**  permitted. Atollic AB permit registered Atollic TrueSTUDIO(R) users the
**  rights to distribute the assembled, compiled & linked contents of this
**  file as part of an application binary file, provided that it is built
**  using the Atollic TrueSTUDIO(R) toolchain.
**
**
*****************************************************************************
*/

/* Includes */
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>
/* Private typedef */

/* Private define  */

/* Private macro */

/* Private variables */


void init(void);
/* Private functions */

/**
**===========================================================================
**
**  Abstract: main program
**
**===========================================================================
*/

int windowSize,hop,numberSlices,xlength,hopOut;
long alpha;
int step;
int debug = 0;
/* macros */
#define TWO_PI (6.2831853071795864769252867665590057683943L)

	/* function prototypes */
	void fft(int N, double (*x)[2], double (*X)[2]);
	void fft_rec(int N, int offset, int delta,
	             double (*x)[2], double (*X)[2], double (*XX)[2]);
	void ifft(int N, double (*x)[2], double (*X)[2]);

	double** createFrames(double * );
	double * fusionFrames(double **);
	double * synthesis(double *,double * );
	void process(double * , double*);
	double** analyse(double * currentFrame);
	double * interpolate(double * outputScaled);


double *phaseCumulative;
double ** magPhase;
double (*currentFrameWindowed)[2];
double (*currentFrameFFT)[2];
static const double wn[256] = {0.000038,0.000339,0.000941,0.001844,0.003047,0.004549,0.006349,0.008447,0.010841,0.013530,0.016512,0.019785,0.023347,0.027196,0.031330,0.035747,0.040443,0.045416,0.050663,0.056180,0.061965,0.068014,0.074322,0.080888,0.087705,0.094771,0.102082,0.109631,0.117416,0.125432,0.133673,0.142135,0.150812,0.159700,0.168792,0.178084,0.187570,0.197244,0.207101,0.217134,0.227338,0.237705,0.248231,0.258908,0.269731,0.280692,0.291785,0.303004,0.314341,0.325791,0.337345,0.348997,0.360740,0.372567,0.384471,0.396444,0.408480,0.420571,0.432710,0.444889,0.457101,0.469340,0.481596,0.493864,0.506136,0.518404,0.530660,0.542899,0.555111,0.567290,0.579429,0.591520,0.603556,0.615529,0.627433,0.639260,0.651003,0.662655,0.674209,0.685659,0.696996,0.708215,0.719308,0.730269,0.741092,0.751769,0.762295,0.772662,0.782866,0.792899,0.802756,0.812430,0.821916,0.831208,0.840300,0.849188,0.857865,0.866327,0.874568,0.882584,0.890369,0.897918,0.905229,0.912295,0.919112,0.925678,0.931986,0.938035,0.943820,0.949337,0.954584,0.959557,0.964253,0.968670,0.972804,0.976653,0.980215,0.983488,0.986470,0.989159,0.991553,0.993651,0.995451,0.996953,0.998156,0.999059,0.999661,0.999962,0.999962,0.999661,0.999059,0.998156,0.996953,0.995451,0.993651,0.991553,0.989159,0.986470,0.983488,0.980215,0.976653,0.972804,0.968670,0.964253,0.959557,0.954584,0.949337,0.943820,0.938035,0.931986,0.925678,0.919112,0.912295,0.905229,0.897918,0.890369,0.882584,0.874568,0.866327,0.857865,0.849188,0.840300,0.831208,0.821916,0.812430,0.802756,0.792899,0.782866,0.772662,0.762295,0.751769,0.741092,0.730269,0.719308,0.708215,0.696996,0.685659,0.674209,0.662655,0.651003,0.639260,0.627433,0.615529,0.603556,0.591520,0.579429,0.567290,0.555111,0.542899,0.530660,0.518404,0.506136,0.493864,0.481596,0.469340,0.457101,0.444889,0.432710,0.420571,0.408480,0.396444,0.384471,0.372567,0.360740,0.348997,0.337345,0.325791,0.314341,0.303004,0.291785,0.280692,0.269731,0.258908,0.248231,0.237705,0.227338,0.217134,0.207101,0.197244,0.187570,0.178084,0.168792,0.159700,0.150812,0.142135,0.133673,0.125432,0.117416,0.109631,0.102082,0.094771,0.087705,0.080888,0.074322,0.068014,0.061965,0.056180,0.050663,0.045416,0.040443,0.035747,0.031330,0.027196,0.023347,0.019785,0.016512,0.013530,0.010841,0.008447,0.006349,0.004549,0.003047,0.001844,0.000941,0.000339,0.000038};
static double weird;
	int main(void)
	{
		int i = 0;
		double ** vector;
		double * x;
		double ** y;
		double ** outputy;
		double * timeStretched;
		double * previousPhaseFrame;

		hop = 32;
		windowSize = 256;
		xlength = 1024;
		step = 8;
		alpha = (pow(2,(step/12)));
		hopOut = round(alpha*hop);
		numberSlices = floor((double)(((xlength)-windowSize)/hop));
		phaseCumulative = malloc(windowSize * sizeof(double));
		x = malloc(xlength * sizeof(double));

		outputy = malloc(numberSlices * sizeof(double *));
		outputy[0] = malloc(numberSlices * windowSize * sizeof(double));

		previousPhaseFrame = malloc(windowSize * sizeof(double));

		init();

		//reads fromt the file
		FILE *fopen(), *fp;
		int c;
		fp = fopen("orig.txt","r");
		int index =-1;
		while (!feof(fp))
		{
			index++;
			if(index == xlength)
				break;
			fscanf(fp,"%lg",&x[index]);

		}
		fclose(fp);
		for(i = 1; i < numberSlices; i++)
			outputy[i] = outputy[0] + i * windowSize;

		y = createFrames(x);

		for(index = 0;index<windowSize;index++)
		{
			previousPhaseFrame[index] = 0;
		}
		for(index = 0;index<numberSlices;index++)
		{
			vector = analyse(y[index]);
			process(y[index],previousPhaseFrame);
			outputy[index] = synthesis(vector[0],vector[1]);
			//printf("%s","end of the major loop\n");
		}
		timeStretched = fusionFrames(outputy);


		for(i = 0;i<xlength;i++)
		{
			printf("%f",timeStretched[i]);
			printf("%s","\n");
		}
		// do the linear interpolation
		return 0;
	}





	/* creates the frames from the samples */
	double** createFrames(double * x)
	{
		int indexTimeStart,indexTimeEnd,i,j;
		double ** vectorFrames;

		vectorFrames = malloc(numberSlices * sizeof(double *));
		vectorFrames[0] = malloc(numberSlices * windowSize * sizeof(double));
		for(i = 1; i < numberSlices; i++)
		{
			vectorFrames[i] = vectorFrames[0] + i * windowSize;
		}
		for ( i = 0; i<numberSlices;i++)
		{
			indexTimeStart = (i)*hop;
			indexTimeEnd = indexTimeStart + windowSize;
			for(j = indexTimeStart;j<indexTimeEnd;j++)
			{
				vectorFrames[i][j-indexTimeStart] = x[j];
			}
		}
		return vectorFrames;
	}


	/* merges the frames back together */
	double * fusionFrames(double **framesMatrix)
	{
		double * vectorTime;
		int timeIndex,index,j;
		vectorTime = malloc (xlength * sizeof(double));
		timeIndex = 0;
		for(index = 0;index<xlength;index++)
			vectorTime[index]=0;

		for (index=0;index<numberSlices;index++)
		{
			for(j = timeIndex;j<timeIndex+windowSize;j++)
			{
				if(debug == 1)
				{
					printf("%d",j);printf("%s"," ");
					//printf("%d",timeIndex+windowSize);
					printf("%s","\n");
				}
				vectorTime[j] +=framesMatrix[index][j-timeIndex];
			}
			timeIndex += hopOut;
		}
		return vectorTime;
	}

	double * interpolate(double * outputScaled)
	{
		int index;
		double * inter;
		inter = malloc(xlength * sizeof(double));
		for(index = 0 ; index< xlength;index++)
		{
			//inter[index] = (outputScaled[index+1] -outputScaled[index])/;
		}
		return inter;
	}

	/*
	 *

	%% Analysis

		% Get current frame to be processed
		currentFrame = y(index,:);

		% Window the frame
		currentFrameWindowed = currentFrame .* wn' / sqrt(((winSize/hop)/2));

		% Get the FFT
		currentFrameWindowedFFT = fft(currentFrameWindowed);

		% Get the magnitude
		magFrame = abs(currentFrameWindowedFFT);

		% Get the angle
		phaseFrame = angle(currentFrameWindowedFFT);

	 */

	double** analyse(double * currentFrame)
	{
		int index;
		int curindex = 0;
		for(index = 0;index<windowSize;index++)
		{
				currentFrameWindowed[0][index] = currentFrame[index] * wn[index]/weird;
				currentFrameWindowed[1][index] = 0;
				if(debug ==1)
				{
					//printf("%f",currentFrame[index] );
					//printf("%s","  ");
					//printf("%d",curindex);
					//printf("%s","  ");
					//printf("%f",wn[curindex]);
					//printf("%s","  ");
					//printf("%f",currentFrame[0][index]);
					//printf("%s","  ");
					printf("%f",currentFrameWindowed[0][index]);
					printf("%s","  ");
					printf("%f",currentFrameWindowed[1][index]);
					printf("%s"," \n");
				}

		}
		if(debug ==1)
			printf("%s","listing from analyse \n");

		fft(windowSize<<1,currentFrameWindowed,currentFrameFFT);
		//magnitude32_32bIn(currentFrameFFt,windowSize*2);
		for(index = 0;index< windowSize;index++)
		{
			if(debug ==1)
			{
				printf("%f",currentFrameFFT[0][index]);
				printf("%s"," ");
				printf("%f",currentFrameFFT[1][index]);
				printf("%s","\n");
			}
			//considered getting these from the dude's library, but my head overheated
			magPhase[0][index] = sqrt(pow(currentFrameFFT[0][index],2) + pow(currentFrameFFT[1][index],2));
			magPhase[1][index]= atan(currentFrameFFT[1][index]/currentFrameFFT[0][index]);

		}
		return magPhase;
	}
	/*
	 *
	%% Processing

		% Get the phase difference
		deltaPhi = phaseFrame - previousPhase;
		previousPhase = phaseFrame;

		% Remove the expected phase difference
		deltaPhiPrime = deltaPhi - hop * 2*pi*(0:(winSize-1))/winSize;

		% Map to -pi/pi range
		deltaPhiPrimeMod = mod(deltaPhiPrime+pi, 2*pi) - pi;

		% Get the true frequency
		trueFreq = 2*pi*(0:(winSize-1))/winSize + deltaPhiPrimeMod/hop;

		% Get the final phase
			phaseCumulative = phaseCumulative + hopOut * trueFreq;

		% Remove the 60 Hz noise. This is not done for now but could be
		% achieved by setting some bins to zero.
	 *
	 */

	void process(double * phaseFrame, double* previousFrame)
	{
		int index;
		double * deltaPhi,* deltaPhiPrime,*deltaPhiPrimeMod,*trueFreq;
		deltaPhi = malloc(windowSize* sizeof(double));
		deltaPhiPrime = malloc(windowSize* sizeof(double));
		deltaPhiPrimeMod = malloc(windowSize* sizeof(double));
		trueFreq = malloc(windowSize* sizeof(double));
		//need to check whether these are vectors or not
		double precompute1 = hop * 2* M_PI/windowSize;
		double precompute2 = TWO_PI/windowSize;
		for(index = 0;index<windowSize;index++)
		{
			deltaPhi[index] = phaseFrame[index] - previousFrame[index];
			previousFrame[index] = phaseFrame[index];
			deltaPhiPrime[index] = deltaPhi[index] - precompute1 * index;

			deltaPhiPrimeMod[index] = fmod(deltaPhiPrime[index]+M_PI ,TWO_PI) -M_PI;
			trueFreq[index] = precompute2*index + deltaPhiPrimeMod[index]/hop;

			phaseCumulative[index] += hopOut+trueFreq[index];
			if(debug == 1)
			{
				printf("%f",phaseCumulative[index]);
				printf("%s","\n");
			}
		}
		if(debug == 1)
			printf("%s","finished the processing funcion\n");
	}

	/*
	 * %% Synthesis

		% Get the magnitude
		outputMag = magFrame;

		% Produce output frame
		outputFrame = real(ifft(outputMag .* exp(j*phaseCumulative)));

		% Save frame that has been processed
		outputy(index,:) = outputFrame .* wn' / sqrt(((winSize/hopOut)/2));
	 */

	double * synthesis(double *magFrame,double * phaseCumulative)
	{
		int index;
		double (*X)[2] = malloc(windowSize<<1*sizeof(double));
		for(index = 0; index <windowSize;index++)
		{
			X[0][index] = magFrame[index]*cos(phaseCumulative[index]);
			X[1][index] = magFrame[index]*sin(phaseCumulative[index]);
		}
		ifft(windowSize,X,X);

		for(index = 0;index<windowSize;index++)
		{
			magFrame[index] = X[0][index];
		}
		if(debug == 1)
			printf("%s","finishind the synthesis\n");
		return magFrame;
	}

	void init(void)
	{
		weird = sqrt(((windowSize/hop)>>1));
		magPhase = malloc(sizeof(double *)<<1);
		magPhase[0] = malloc(windowSize<<1 * sizeof(double));
		magPhase[1] = magPhase[0] + windowSize;
		currentFrameWindowed = malloc(windowSize<<1* sizeof(double));
		currentFrameFFT = malloc(windowSize<<1* sizeof(double));
	}



	/* FFT */
	void fft(int N, double (*x)[2], double (*X)[2])
	{
	  /* Declare a pointer to scratch space. */
	  double (*XX)[2] = malloc(2 * N * sizeof(double));

	  /* Calculate FFT by a recursion. */
	  fft_rec(N, 0, 1, x, X, XX);

	  /* Free memory. */
	  free(XX);
	}

	/* FFT recursion */
	void fft_rec(int N, int offset, int delta,
	             double (*x)[2], double (*X)[2], double (*XX)[2])
	{
	  int N2 = N/2;            /* half the number of points in FFT */
	  int k;                   /* generic index */
	  double cs, sn;           /* cosine and sine */
	  int k00, k01, k10, k11;  /* indices for butterflies */
	  double tmp0, tmp1;       /* temporary storage */

	  if(N != 2)  /* Perform recursive step. */
	    {
	      /* Calculate two (N/2)-point DFT's. */
	      fft_rec(N2, offset, 2*delta, x, XX, X);
	      fft_rec(N2, offset+delta, 2*delta, x, XX, X);

	      /* Combine the two (N/2)-point DFT's into one N-point DFT. */
	      for(k=0; k<N2; k++)
	        {
	          k00 = offset + k*delta;    k01 = k00 + N2*delta;
	          k10 = offset + 2*k*delta;  k11 = k10 + delta;
	          cs = cos(TWO_PI*k/(double)N); sn = sin(TWO_PI*k/(double)N);
	          tmp0 = cs * XX[k11][0] + sn * XX[k11][1];
	          tmp1 = cs * XX[k11][1] - sn * XX[k11][0];
	          X[k01][0] = XX[k10][0] - tmp0;
	          X[k01][1] = XX[k10][1] - tmp1;
	          X[k00][0] = XX[k10][0] + tmp0;
	          X[k00][1] = XX[k10][1] + tmp1;
	        }
	    }
	  else  /* Perform 2-point DFT. */
	    {
	      k00 = offset; k01 = k00 + delta;
	      X[k01][0] = x[k00][0] - x[k01][0];
	      X[k01][1] = x[k00][1] - x[k01][1];
	      X[k00][0] = x[k00][0] + x[k01][0];
	      X[k00][1] = x[k00][1] + x[k01][1];
	    }
	}

	/* IFFT */
	void ifft(int N, double (*x)[2], double (*X)[2])
	{
	  int N2 = N/2;       /* half the number of points in IFFT */
	  int i;              /* generic index */
	  double tmp0, tmp1;  /* temporary storage */

	  /* Calculate IFFT via reciprocity property of DFT. */
	  fft(N, X, x);
	  x[0][0] = x[0][0]/N;
	  x[0][1] = x[0][1]/N;
	  x[N2][0] = x[N2][0]/N;
	  x[N2][1] = x[N2][1]/N;
	  for(i=1; i<N2; i++)
	    {
	      tmp0 = x[i][0]/N;
	      tmp1 = x[i][1]/N;
	      x[i][0] = x[N-i][0]/N;
	      x[i][1] = x[N-i][1]/N;
	      x[N-i][0] = tmp0;
	      x[N-i][1] = tmp1;
	    }
	}
