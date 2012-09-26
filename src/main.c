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
double alpha;
int step;
int debug = 0;
int timeSize;
int interlength;
double *phaseCumulative;

//values only defined once...rather than repeatedly in loops
static double weird;
static double weird2;
static double precompute1;
static double precompute2;

/* macros */
#define TWO_PI (6.2831853071795864769252867665590057683943L)

/* function prototypes */
void fft(int N, double (*x)[2], double (*X)[2]);
void fft_rec(int N, int offset, int delta,
             double (*x)[2], double (*X)[2], double (*XX)[2]);
void ifft(int N, double (*x)[2], double (*X)[2]);

void createFrames(double * ,double(*)[1024]);
double * fusionFrames(double **);
double * synthesis(double(*)[2]);
void process(double(*)[2] , double*);
void analyse(double * currentFrame,double(*)[2]);
double * interpolate(double * outputScaled);

//1024 hanning wiindow
static const double wn[1024] ={0.0000,0.0000,0.0001,0.0002,0.0002,0.0003,0.0005,0.0006,0.0008,0.0009,0.0011,0.0014,0.0016,0.0018,0.0021,0.0024,0.0027,0.0030,0.0034,0.0038,0.0041,0.0045,0.0050,0.0054,0.0059,0.0063,0.0068,0.0073,0.0079,0.0084,0.0090,0.0096,0.0102,0.0108,0.0115,0.0121,0.0128,0.0135,0.0142,0.0150,0.0157,0.0165,0.0173,0.0181,0.0189,0.0197,0.0206,0.0215,0.0224,0.0233,0.0242,0.0252,0.0262,0.0271,0.0281,0.0292,0.0302,0.0313,0.0323,0.0334,0.0345,0.0357,0.0368,0.0380,0.0392,0.0404,0.0416,0.0428,0.0441,0.0453,0.0466,0.0479,0.0492,0.0506,0.0519,0.0533,0.0547,0.0561,0.0575,0.0589,0.0604,0.0618,0.0633,0.0648,0.0664,0.0679,0.0694,0.0710,0.0726,0.0742,0.0758,0.0774,0.0791,0.0807,0.0824,0.0841,0.0858,0.0875,0.0893,0.0910,0.0928,0.0946,0.0964,0.0982,0.1000,0.1019,0.1038,0.1056,0.1075,0.1094,0.1113,0.1133,0.1152,0.1172,0.1192,0.1212,0.1232,0.1252,0.1272,0.1293,0.1313,0.1334,0.1355,0.1376,0.1397,0.1419,0.1440,0.1462,0.1483,0.1505,0.1527,0.1549,0.1572,0.1594,0.1617,0.1639,0.1662,0.1685,0.1708,0.1731,0.1754,0.1778,0.1801,0.1825,0.1848,0.1872,0.1896,0.1920,0.1945,0.1969,0.1993,0.2018,0.2043,0.2067,0.2092,0.2117,0.2142,0.2167,0.2193,0.2218,0.2244,0.2269,0.2295,0.2321,0.2347,0.2373,0.2399,0.2425,0.2451,0.2478,0.2504,0.2531,0.2558,0.2585,0.2611,0.2638,0.2665,0.2693,0.2720,0.2747,0.2775,0.2802,0.2830,0.2857,0.2885,0.2913,0.2941,0.2969,0.2997,0.3025,0.3053,0.3081,0.3110,0.3138,0.3167,0.3195,0.3224,0.3252,0.3281,0.3310,0.3339,0.3368,0.3397,0.3426,0.3455,0.3484,0.3513,0.3543,0.3572,0.3601,0.3631,0.3660,0.3690,0.3719,0.3749,0.3779,0.3809,0.3838,0.3868,0.3898,0.3928,0.3958,0.3988,0.4018,0.4048,0.4078,0.4108,0.4138,0.4169,0.4199,0.4229,0.4259,0.4290,0.4320,0.4351,0.4381,0.4411,0.4442,0.4472,0.4503,0.4533,0.4564,0.4594,0.4625,0.4655,0.4686,0.4717,0.4747,0.4778,0.4808,0.4839,0.4870,0.4900,0.4931,0.4962,0.4992,0.5023,0.5054,0.5084,0.5115,0.5146,0.5176,0.5207,0.5237,0.5268,0.5299,0.5329,0.5360,0.5390,0.5421,0.5451,0.5482,0.5512,0.5543,0.5573,0.5604,0.5634,0.5665,0.5695,0.5725,0.5756,0.5786,0.5816,0.5846,0.5877,0.5907,0.5937,0.5967,0.5997,0.6027,0.6057,0.6087,0.6117,0.6147,0.6177,0.6206,0.6236,0.6266,0.6295,0.6325,0.6354,0.6384,0.6413,0.6443,0.6472,0.6501,0.6531,0.6560,0.6589,0.6618,0.6647,0.6676,0.6704,0.6733,0.6762,0.6791,0.6819,0.6848,0.6876,0.6905,0.6933,0.6961,0.6989,0.7017,0.7045,0.7073,0.7101,0.7129,0.7157,0.7184,0.7212,0.7239,0.7267,0.7294,0.7321,0.7348,0.7375,0.7402,0.7429,0.7456,0.7482,0.7509,0.7535,0.7562,0.7588,0.7614,0.7640,0.7666,0.7692,0.7718,0.7744,0.7769,0.7795,0.7820,0.7845,0.7870,0.7895,0.7920,0.7945,0.7970,0.7994,0.8019,0.8043,0.8068,0.8092,0.8116,0.8140,0.8163,0.8187,0.8211,0.8234,0.8257,0.8281,0.8304,0.8327,0.8349,0.8372,0.8395,0.8417,0.8439,0.8462,0.8484,0.8506,0.8527,0.8549,0.8571,0.8592,0.8613,0.8634,0.8655,0.8676,0.8697,0.8717,0.8738,0.8758,0.8778,0.8798,0.8818,0.8838,0.8857,0.8877,0.8896,0.8915,0.8934,0.8953,0.8972,0.8990,0.9009,0.9027,0.9045,0.9063,0.9081,0.9098,0.9116,0.9133,0.9150,0.9167,0.9184,0.9201,0.9218,0.9234,0.9250,0.9266,0.9282,0.9298,0.9313,0.9329,0.9344,0.9359,0.9374,0.9389,0.9403,0.9418,0.9432,0.9446,0.9460,0.9474,0.9488,0.9501,0.9514,0.9527,0.9540,0.9553,0.9566,0.9578,0.9590,0.9602,0.9614,0.9626,0.9638,0.9649,0.9660,0.9671,0.9682,0.9693,0.9703,0.9713,0.9724,0.9734,0.9743,0.9753,0.9762,0.9772,0.9781,0.9790,0.9798,0.9807,0.9815,0.9823,0.9831,0.9839,0.9847,0.9854,0.9861,0.9868,0.9875,0.9882,0.9889,0.9895,0.9901,0.9907,0.9913,0.9918,0.9924,0.9929,0.9934,0.9939,0.9944,0.9948,0.9953,0.9957,0.9961,0.9964,0.9968,0.9971,0.9974,0.9977,0.9980,0.9983,0.9985,0.9988,0.9990,0.9992,0.9993,0.9995,0.9996,0.9997,0.9998,0.9999,0.9999,1.0000,1.0000,1.0000,1.0000,0.9999,0.9999,0.9998,0.9997,0.9996,0.9995,0.9993,0.9992,0.9990,0.9988,0.9985,0.9983,0.9980,0.9977,0.9974,0.9971,0.9968,0.9964,0.9961,0.9957,0.9953,0.9948,0.9944,0.9939,0.9934,0.9929,0.9924,0.9918,0.9913,0.9907,0.9901,0.9895,0.9889,0.9882,0.9875,0.9868,0.9861,0.9854,0.9847,0.9839,0.9831,0.9823,0.9815,0.9807,0.9798,0.9790,0.9781,0.9772,0.9762,0.9753,0.9743,0.9734,0.9724,0.9713,0.9703,0.9693,0.9682,0.9671,0.9660,0.9649,0.9638,0.9626,0.9614,0.9602,0.9590,0.9578,0.9566,0.9553,0.9540,0.9527,0.9514,0.9501,0.9488,0.9474,0.9460,0.9446,0.9432,0.9418,0.9403,0.9389,0.9374,0.9359,0.9344,0.9329,0.9313,0.9298,0.9282,0.9266,0.9250,0.9234,0.9218,0.9201,0.9184,0.9167,0.9150,0.9133,0.9116,0.9098,0.9081,0.9063,0.9045,0.9027,0.9009,0.8990,0.8972,0.8953,0.8934,0.8915,0.8896,0.8877,0.8857,0.8838,0.8818,0.8798,0.8778,0.8758,0.8738,0.8717,0.8697,0.8676,0.8655,0.8634,0.8613,0.8592,0.8571,0.8549,0.8527,0.8506,0.8484,0.8462,0.8439,0.8417,0.8395,0.8372,0.8349,0.8327,0.8304,0.8281,0.8257,0.8234,0.8211,0.8187,0.8163,0.8140,0.8116,0.8092,0.8068,0.8043,0.8019,0.7994,0.7970,0.7945,0.7920,0.7895,0.7870,0.7845,0.7820,0.7795,0.7769,0.7744,0.7718,0.7692,0.7666,0.7640,0.7614,0.7588,0.7562,0.7535,0.7509,0.7482,0.7456,0.7429,0.7402,0.7375,0.7348,0.7321,0.7294,0.7267,0.7239,0.7212,0.7184,0.7157,0.7129,0.7101,0.7073,0.7045,0.7017,0.6989,0.6961,0.6933,0.6905,0.6876,0.6848,0.6819,0.6791,0.6762,0.6733,0.6704,0.6676,0.6647,0.6618,0.6589,0.6560,0.6531,0.6501,0.6472,0.6443,0.6413,0.6384,0.6354,0.6325,0.6295,0.6266,0.6236,0.6206,0.6177,0.6147,0.6117,0.6087,0.6057,0.6027,0.5997,0.5967,0.5937,0.5907,0.5877,0.5846,0.5816,0.5786,0.5756,0.5725,0.5695,0.5665,0.5634,0.5604,0.5573,0.5543,0.5512,0.5482,0.5451,0.5421,0.5390,0.5360,0.5329,0.5299,0.5268,0.5237,0.5207,0.5176,0.5146,0.5115,0.5084,0.5054,0.5023,0.4992,0.4962,0.4931,0.4900,0.4870,0.4839,0.4808,0.4778,0.4747,0.4717,0.4686,0.4655,0.4625,0.4594,0.4564,0.4533,0.4503,0.4472,0.4442,0.4411,0.4381,0.4351,0.4320,0.4290,0.4259,0.4229,0.4199,0.4169,0.4138,0.4108,0.4078,0.4048,0.4018,0.3988,0.3958,0.3928,0.3898,0.3868,0.3838,0.3809,0.3779,0.3749,0.3719,0.3690,0.3660,0.3631,0.3601,0.3572,0.3543,0.3513,0.3484,0.3455,0.3426,0.3397,0.3368,0.3339,0.3310,0.3281,0.3252,0.3224,0.3195,0.3167,0.3138,0.3110,0.3081,0.3053,0.3025,0.2997,0.2969,0.2941,0.2913,0.2885,0.2857,0.2830,0.2802,0.2775,0.2747,0.2720,0.2693,0.2665,0.2638,0.2611,0.2585,0.2558,0.2531,0.2504,0.2478,0.2451,0.2425,0.2399,0.2373,0.2347,0.2321,0.2295,0.2269,0.2244,0.2218,0.2193,0.2167,0.2142,0.2117,0.2092,0.2067,0.2043,0.2018,0.1993,0.1969,0.1945,0.1920,0.1896,0.1872,0.1848,0.1825,0.1801,0.1778,0.1754,0.1731,0.1708,0.1685,0.1662,0.1639,0.1617,0.1594,0.1572,0.1549,0.1527,0.1505,0.1483,0.1462,0.1440,0.1419,0.1397,0.1376,0.1355,0.1334,0.1313,0.1293,0.1272,0.1252,0.1232,0.1212,0.1192,0.1172,0.1152,0.1133,0.1113,0.1094,0.1075,0.1056,0.1038,0.1019,0.1000,0.0982,0.0964,0.0946,0.0928,0.0910,0.0893,0.0875,0.0858,0.0841,0.0824,0.0807,0.0791,0.0774,0.0758,0.0742,0.0726,0.0710,0.0694,0.0679,0.0664,0.0648,0.0633,0.0618,0.0604,0.0589,0.0575,0.0561,0.0547,0.0533,0.0519,0.0506,0.0492,0.0479,0.0466,0.0453,0.0441,0.0428,0.0416,0.0404,0.0392,0.0380,0.0368,0.0357,0.0345,0.0334,0.0323,0.0313,0.0302,0.0292,0.0281,0.0271,0.0262,0.0252,0.0242,0.0233,0.0224,0.0215,0.0206,0.0197,0.0189,0.0181,0.0173,0.0165,0.0157,0.0150,0.0142,0.0135,0.0128,0.0121,0.0115,0.0108,0.0102,0.0096,0.0090,0.0084,0.0079,0.0073,0.0068,0.0063,0.0059,0.0054,0.0050,0.0045,0.0041,0.0038,0.0034,0.0030,0.0027,0.0024,0.0021,0.0018,0.0016,0.0014,0.0011,0.0009,0.0008,0.0006,0.0005,0.0003,0.0002,0.0002,0.0001,0.0000,0.0000};
//static const double wn [20] = {  22.2135971069296e-003,    86.8806128420025e-003,   188.255099070633e-003,    317.329487816802e-003,    462.634953206788e-003,   611.260466978157e-003,    750.000000000000e-003,    866.525935914913e-003, 950.484433951210e-003,    994.415413112564e-003,    994.415413112564e-003,    950.484433951210e-003, 866.525935914913e-003,    750.000000000000e-003,    611.260466978157e-003,    462.634953206788e-003,317.329487816802e-003, 188.255099070633e-003,    86.8806128420025e-003,    22.2135971069296e-003};

	int main(void)
	{
		int i = 0;
		double (*magPhase)[2];
		double * x;
		
		double (*y)[1024];
		double **outputy;
		double * timeStretched;
		double * previousPhaseFrame;
		double * synthesizedOutput;


		//hard coded the test values ie hop, windowSize and the input length
		hop = 128;
		windowSize =1024;
		xlength =  7411968 + hop*3; // adds zeros to the input... think this is for time scaling sort of
		step = 10;
		alpha = pow(2,(step/(float)12));
		hopOut = round(alpha*hop);
		numberSlices = (((xlength-windowSize)/hop));
		
		phaseCumulative = malloc(windowSize * sizeof(double));
		x = malloc(xlength * sizeof(double));

		outputy = malloc(numberSlices * sizeof(double *));
		outputy[0] = malloc(numberSlices * windowSize * sizeof(double));

		previousPhaseFrame = malloc(windowSize * sizeof(double));
		magPhase = malloc(2*windowSize*sizeof(double));
		init();
	
		//reads fromt the file
		FILE *fopen(), *fp;

		fp = fopen("../pitchshifter/test.txt","r"); //hard coded the test input source
		int index =0;
		for(index = 0;index<hop*3;index++)
		{
			x[index] = 0;
		}
		while (!feof(fp) )
		{
			index++;
			fscanf(fp,"%lg",&x[index]);
			if(debug == 1)
			{
				printf("%lg",x[index]);
				printf("%s","\n");
			}
		}
		fclose(fp);
		
		for(i = 1; i < numberSlices; i++)
			outputy[i] = outputy[0] + i * windowSize;
		
		if(xlength %hop != 0)
			xlength = numberSlices*hop + windowSize;
		y = malloc((xlength/hop)*windowSize* sizeof(double *)); //note the number of column is different from the number of slicerSlices
		
		createFrames(x,y);
		
		for(index = 0;index<windowSize;index++)
		{
			previousPhaseFrame[index] = 0;
		}
		for(index = 0;index<numberSlices;index++)
		{
			analyse(y[index],magPhase);
			process(magPhase,previousPhaseFrame);
			outputy[index] = synthesis(magPhase);			
		}
		timeStretched = fusionFrames(outputy);//should check this timeStretched value before linear interpolation, to see what comes out

		
		// do the linear interpolation
		synthesizedOutput = interpolate(timeStretched);

		for(i = 0;i<interlength;i++)
		{
			printf("%f",synthesizedOutput[i]);
			printf("%s","\n");
		}
		return 0;
	}





	/* creates the frames from the samples */
	void createFrames(double * x, double (*vectorFrames)[windowSize])
	{
		int indexTimeStart,indexTimeEnd,i,j;
		
		//truncate if not a multiple of hop
		

		for ( i = 0; i<numberSlices;i++)
		{
			indexTimeStart = (i)*hop+1;
			indexTimeEnd = indexTimeStart-1 + windowSize;

			for(j = indexTimeStart;j<=indexTimeEnd+ windowSize;j++)
			{
				vectorFrames[i][j-indexTimeStart] = x[j];
			}
		}

		if(debug ==1)
			for(i = 0;i<xlength/hop;i++)
			{
				for(j = 0;j<windowSize;j++)
				{
					printf("%f",vectorFrames[i][j]);
					printf("%s"," ");
				}
				printf("%s","\n");
			}
	}


	/* merges the frames back together */
	// think there's a bug in here somewhere because i'm losing data at the end
	double * fusionFrames(double **framesMatrix)
	{
		double * vectorTime;
		int timeIndex,index,j;
		timeSize = numberSlices * hopOut - hopOut+windowSize ;
		vectorTime = malloc (timeSize* sizeof(double));
		timeIndex = 0;
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
		
		
		if(debug ==1)
			for(index = 0;index<timeSize;index++)
			{
				printf("%f",vectorTime[index]);
				printf("%s"," ");
				printf("%s","\n");
			}
		return vectorTime;
	}


	//does the linear interpolation
	double * interpolate(double * outputScaled)
	{
		int index;
		int aindex;
		double * inter;
		double temp = 0;
		interlength = floor(timeSize/alpha);
		inter = malloc(interlength * sizeof(double));

		int up,down;
		for(index = 0 ; index< interlength;index++,temp+=alpha)
		{
			down = floor(temp);
			up = ceil(temp);
			if(up == down)
			{
			//exact point
				inter[index] = outputScaled[up];
			}
			else if(up>=interlength)
			{
			//uses linear extrapolation if the value is out of bounds
				inter[index] = outputScaled[interlength] + (up-interlength)*(outputScaled[interlength]-outputScaled[interlength-1])/(interlength - interlength-1);
			}
			else
			{
				//linear interpolation
				inter[index] = outputScaled[down]+((temp - down)*(outputScaled[up]-outputScaled[down]))/(up-down);
			}
			if(debug == 1)
			{
				printf("%f",inter[index]);
				printf("%s","\n");
			}
		}
		return inter;
	}

	/*
	 *

	%% Analysis

		Analyses the input frame

	 */

	void analyse(double * currentFrame,double(*magPhase)[2] )
	{
		int index;
		int curindex = 0;
		
		double (*currentFrameWindowed) [2];
		double (*currentFrameFFT) [2];
		currentFrameWindowed = malloc(windowSize * 2* sizeof(double));
		currentFrameFFT = malloc(windowSize* 2* sizeof(double));
		for(index = 0;index<windowSize;index++)
		{
				currentFrameWindowed[index][0] = (currentFrame[index] * wn[index])/weird;
				currentFrameWindowed[index][1] = 0;
				if(debug ==1)
				{
					printf("%lg",currentFrameFFT[index][0]);
					printf("%s","\t\t");
					printf("%lg",currentFrameWindowed[index][1]);
					printf("%s","\t\t");
					printf("%s"," \n");
				}
		}
		
		if(debug ==1)
			printf("%s","listing from analyse \n");

		fft(windowSize,currentFrameWindowed,currentFrameFFT);//creates the ffts
		for(index = 0;index< windowSize;index++)
		{
			if(debug ==1)
			{
				printf("%lg",currentFrameWindowed[index][0]);
				printf("%s","\t\t");
				printf("%lg",currentFrameWindowed[index][1]);
				printf("%s","\t\t");
				printf("%lg",currentFrameFFT[index][0]);
				printf("%s","\t\t");
				printf("%lg",currentFrameFFT[index][1]);
				printf("%s","i\n");
			}
			//considered getting these from the dude's library, but my head overheated
			magPhase[index][0] = sqrt(pow(currentFrameFFT[index][0],2) + pow(currentFrameFFT[index][1],2));
			if(currentFrameFFT[index][1] == 0)
			{
				magPhase[index][1] = 0;
			}
			else
				magPhase[index][1]= atan((currentFrameFFT[index][1]/currentFrameFFT[index][0]));

		}
		free (currentFrameWindowed);
		free (currentFrameFFT);
	}
	/*
	 *
	%% Processing

		% Get the phase difference
		

		% Remove the 60 Hz noise. This is not done for now but could be
		% achieved by setting some bins to zero.
	 *
	 */

	//this does the processing, just elementary maths in here
	void process(double(*magPhase)[2] , double* previousFrame)
	{
		int index;
		double * deltaPhi,* deltaPhiPrime,*deltaPhiPrimeMod,*trueFreq;
		//don't need these in the loop
		
		
		deltaPhi = malloc(windowSize* sizeof(double));
		deltaPhiPrime = malloc(windowSize* sizeof(double));
		deltaPhiPrimeMod = malloc(windowSize* sizeof(double));
		trueFreq = malloc(windowSize* sizeof(double));
		
		for(index = 0;index<windowSize;index++)
		{
			deltaPhi[index] = magPhase[index][1] - previousFrame[index];
			
			previousFrame[index] = magPhase[index][1];
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
		free(deltaPhi);
		free(deltaPhiPrime);
		free(deltaPhiPrimeMod);
		free (trueFreq);
		if(debug == 1)
			printf("%s","finished the processing funcion\n");
	}

	/*
	 * %% Synthesis

		
	 */

	double * synthesis(double (*magPhase)[2])
	{
		int index;
		double (*x)[2] = malloc(windowSize*2*sizeof(double)); //frequency domain
		double (*X)[2] = malloc(windowSize*2*sizeof(double)); //time domain
		double * magFrame = malloc(windowSize*sizeof(double));
		for(index = 0; index <windowSize;index++)
		{
			x[index][0] = magPhase[index][0]*cos(phaseCumulative[index]);
			x[index][1] = magPhase[index][0]*sin(phaseCumulative[index]);
			
			if(debug == 1)
			{
				printf("%lg",x[index][0]);
				printf("%s","\t");
				printf("%lg",x[index][1]);
				printf("%s","\n");
			}
		};
		
		ifft(windowSize,X,x); //inverse FFT
		
		for(index = 0;index<windowSize;index++)
		{
			magFrame[index] = X[index][0]*(wn[index]/weird2);
			if(debug ==1)
			{
				printf("%lg",magFrame[index]);
				printf("%s","\n");
			}
		}
		free(X);
		free(x);
		if(debug == 1)
			printf("%s","finishind the synthesis\n");
		return magFrame;
	}

//think i can stack more stuff into here for optimisation
	void init(void)
	{
		//values use over and over in the loop, should probably make them global
		weird = sqrt(((windowSize/(float)hop)/2));
		weird2 = sqrt(((windowSize/(float)hopOut)/2));
		precompute2 = TWO_PI/windowSize;
		precompute1 = hop * precompute2;
	}


/* land of the FFT's from some other dude's library, will reference when the sun comes up*/


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
