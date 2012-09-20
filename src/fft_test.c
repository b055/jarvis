/*----------------------------------------------
   fft_test.c - demonstration program for fft.c
  ----------------------------------------------*/

/******************************************************************************
 * This program demonstrates how to use the file fft.c to calculate an FFT    *
 * of given time-domain samples, as well as to calculate an inverse FFT       *
 * (IFFT) of given frequency-domain samples.  First, N complex-valued time-   *
 * domain samples x, in rectangular form (Re x, Im x), are read from a        *
 * specified file; the 2N values are assumed to be separated by whitespace.   *
 * Then, an N-point FFT of these samples is found by calling the function     *
 * fft, thereby yielding N complex-valued frequency-domain samples X in       *
 * rectangular form (Re X, Im X).  Next, an N-point IFFT of these samples is  *
 * is found by calling the function ifft, thereby recovering the original     *
 * samples x.  Finally, the calculated samples X are saved to a specified     *
 * file, if desired.                                                          *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "fft.c"

int main()
{
  int i;                    /* generic index */
  char file[FILENAME_MAX];  /* name of data file */
  int N;                    /* number of points in FFT */
  double (*x)[2];           /* pointer to time-domain samples */
  double (*X)[2];           /* pointer to frequency-domain samples */
  double dummy;             /* scratch variable */
  FILE *fp;                 /* file pointer */

  /* Get name of input file of time-domain samples x. */
  printf("Input file for time-domain samples x(n)? ");
  scanf("%s", file);

  /* Read through entire file to get number N of points in FFT. */
  if(!(fp = fopen(file, "r")))
    {
      printf("   File \'%s\' could not be opened!", file);
      exit(EXIT_FAILURE);
    }
  N=0;
  while(fscanf(fp, "%lg%lg", &dummy, &dummy) == 2) N++;
  printf("N = %d", N);

  /* Check that N = 2^n for some integer n >= 1. */
  if(N >= 2)
    {
      i = N;
      while(i==2*(i/2)) i = i/2;  /* While i is even, factor out a 2. */
    }  /* For N >=2, we now have N = 2^n iff i = 1. */
  if(N < 2 || i != 1)
    {
      printf(", which does not equal 2^n for an integer n >= 1.");
      exit(EXIT_FAILURE);
    }

  /* Allocate time- and frequency-domain memory. */
  x = malloc(2 * N * sizeof(double));
  X = malloc(2 * N * sizeof(double));

  /* Get time-domain samples. */
  rewind(fp);
  for(i=0; i<N; i++) fscanf(fp, "%lg%lg", &x[i][0], &x[i][1]);
  fclose(fp);

  /* Calculate FFT. */
  fft(N, x, X);

  /* Print time-domain samples and resulting frequency-domain samples. */
  printf("\nx(n):");
  for(i=0; i<N; i++) printf("\n   n=%d: %12f %12f", i, x[i][0], x[i][1]);
  printf("\nX(k):");
  for(i=0; i<N; i++) printf("\n   k=%d: %12f %12f", i, X[i][0], X[i][1]);

  /* Clear time-domain samples and calculate IFFT. */
  for(i=0; i<N; i++) x[i][0] = x[i][1] = 0;
  ifft(N, x, X);

  /* Print recovered time-domain samples. */
  printf("\nx(n):");
  for(i=0; i<N; i++) printf("\n   n=%d: %12f %12f", i, x[i][0], x[i][1]);

  /* Write frequency-domain samples X to a file, if desired. */
  printf("\nOutput file for frequency-domain samples X(k)?"
         "\n   (if none, abort program): ");
  scanf("%s", file);
  if(!(fp = fopen(file, "w")))
    {
      printf("   File \'%s\' could not be opened!", file);
      exit(EXIT_FAILURE);
    }
  for(i=0; i<N; i++) fprintf(fp, "%23.15e  %23.15e\n", X[i][0], X[i][1]);
  fclose(fp);
  printf("Samples X(k) were written to file %s.", file);

  /* Free memory. */
  free(x);
  free(X);

  return 0;
}

/*==============================================================================
 * Program output (example)
 *==============================================================================
 * Input file for time-domain samples x(n)? data.txt
 * N = 8
 * x(n):
 *    n=0:     3.600000     2.600000
 *    n=1:     2.900000     6.300000
 *    n=2:     5.600000     4.000000
 *    n=3:     4.800000     9.100000
 *    n=4:     3.300000     0.400000
 *    n=5:     5.900000     4.800000
 *    n=6:     5.000000     2.600000
 *    n=7:     4.300000     4.100000
 * X(k):
 *    k=0:    35.400000    33.900000
 *    k=1:     3.821320     0.892893
 *    k=2:    -5.800000    -3.300000
 *    k=3:     5.971068     7.042641
 *    k=4:    -0.400000   -14.700000
 *    k=5:    -0.421320     2.307107
 *    k=6:    -1.600000    -3.900000
 *    k=7:    -8.171068    -1.442641
 * x(n):
 *    n=0:     3.600000     2.600000
 *    n=1:     2.900000     6.300000
 *    n=2:     5.600000     4.000000
 *    n=3:     4.800000     9.100000
 *    n=4:     3.300000     0.400000
 *    n=5:     5.900000     4.800000
 *    n=6:     5.000000     2.600000
 *    n=7:     4.300000     4.100000
 * Output file for frequency-domain samples X(k)?
 *    (if none, abort program): X.txt
 * Samples X(k) were written to file X.txt.
 */
