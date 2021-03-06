/* http://www.engr.colostate.edu/ECE423/lab06/Lab_Notes/Lab_Notes_6.pdf */
#include <stdio.h>
#include <math.h>

#define N 8
#define M 3
#define N2 4
#define PI 3.14159265358979 // fixed point approx of PI

typedef struct {float real, imag; } COMPLEX;

void FFT_func(COMPLEX *X, COMPLEX *W) {

  COMPLEX temp;

  short i,j,k;
  short i_lower;
  short step;
  short stage;
  short DFTpts;
  short numBF;

  for(i=0;i<N2;i++)
  {
    W[i].real = cos(2.0*PI*i/N);
    W[i].imag = sin(2.0*PI*i/N);
  }

  j = 0;

  for (i=1;i<(N-1);i++)
  {
    k = N2;
    while(k<=j)
    {
      j = j-k;
      k = k>>1;
    }

    j = j + k;

    if (i<j)
    {
      temp.real = X[j].real;
      temp.imag = X[j].imag;
      X[j].real = X[i].real;
      X[j].imag = X[i].imag;
      X[i].real = temp.real;
      X[i].imag = temp.imag;
    }

  }

  step = N2;
  for (stage=1;stage <= M; stage++)
  {
    DFTpts = 1 << stage;
    numBF = DFTpts/2;
    k = 0;

    for (j=0;j<numBF;j++)
    {
      for (i=j;i<N;i+= DFTpts)
      {

        i_lower = i + numBF;
        temp.real = X[i_lower].real*W[k].real + X[i_lower].imag*W[k].imag;
        temp.imag = X[i_lower].imag*W[k].real - X[i_lower].real*W[k].imag;
  
        X[i_lower].real = X[i].real - temp.real;
        X[i_lower].imag = X[i].imag - temp.imag;
  
        X[i].real = X[i].real + temp.real;
        X[i].imag = X[i].imag + temp.imag;
      }
      k += step;
    }
    step = step/2;
  }

}

int main()
{
  short i;
  COMPLEX X[N];
  COMPLEX W[N2];

  // calculate twiddle factors
  for(i=0;i<N2;i++)
  {
    W[i].real = cos(2.0*PI*i/N);
    W[i].imag = sin(2.0*PI*i/N);
  }

  //initialize input array
  for(i=0;i<N;i++)
  {
    X[i].imag = 0.0;
  }

  X[0].real = 0.3535;
  X[1].real = 0.3535;
  X[2].real = 0.6464;
  X[3].real = 1.0607;
  X[4].real = 0.3535;
  X[5].real = -1.0607;
  X[6].real = -1.3535;
  X[7].real = -0.3535;

  FFT_func(X,W);

  for(i=0;i<N;i++)
  {
    printf("X[%d] = \t%10.8f + j %10.8f\n", i, X[i].real,X[i].imag);
  }
  return 0;
}