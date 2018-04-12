/*
The MIT License (MIT)

Copyright (c) <2017> tekito2

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/
/*
  --Example--

#include <iostream>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"

void canny(int rows, int cols, const unsigned char *src, unsigned char *dst, int th1, int th2, int grad);

int main() {
    // Open the image of Lena
    cv::Mat image = cv::imread("lena.png");

    cv::Mat x;
    cv::cvtColor(image, x, CV_RGB2GRAY);

    // create clone
    cv::Mat z = x.clone();
    // Just borrow the data member of cv::Mat.
    canny(x.rows, x.cols, x.data, z.data, 60, 200, 1);

    cv::imshow("Lena and tiny Canny", z);
    // Waiting to a key
    cv::waitKey(0);

    return 0;
}

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void judge(int l, 
       double* a1,
       double* a2,
       double* a3, 
       double* b2,
       double* b22,
       unsigned char* c2);

void  determine(int l,
        double* a2, 
        unsigned char* c1,
        unsigned char* c2,
        unsigned char* c3, 
        unsigned char* dst, 
        int max,
        int min);

void canny(int rows,
       int cols,
       const unsigned char *src,
       unsigned char *dst,
       int th1,
       int th2,
       int grad)
{
  int i, j;

  //Sobel kernel
  int sky[9] = { -1, -2, -1, 0, 0, 0, 1, 2, 1 };
  int skx[9] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };

  int * sky1 = sky;
  int * sky2 = sky + 3;
  int * sky3 = sky + 6;

  int * skx1 = skx;
  int * skx2 = skx + 3;
  int * skx3 = skx + 6;

  double fx;
  double fy;
  double mag;

  const int maxth = (th1 > th2) ? th1 : th2;
  const int minth = (th1 < th2) ? th1 : th2;

  // for Amplitude
  double * a1 = (double *) malloc(sizeof(double) * cols);
  double * a2 = (double *) malloc(sizeof(double) * cols);
  double * a3 = (double *) malloc(sizeof(double) * cols);

  // for Fx
  double * b1 = (double *) malloc(sizeof(double) * cols);
  double * b2 = (double *) malloc(sizeof(double) * cols);
  double * b3 = (double *) malloc(sizeof(double) * cols);

  // for Fy
  double * b21 = (double *) malloc(sizeof(double) * cols);
  double * b22 = (double *) malloc(sizeof(double) * cols);
  double * b23 = (double *) malloc(sizeof(double) * cols);
  double * dtmp;

  unsigned char * c1 = (unsigned char *) malloc(sizeof(unsigned char) * cols);
  unsigned char * c2 = (unsigned char *) malloc(sizeof(unsigned char) * cols);
  unsigned char * c3 = (unsigned char *) malloc(sizeof(unsigned char) * cols);
  unsigned char * ctmp;

  const unsigned char * p1 = src;
  const unsigned char * p2 = src + cols;
  const unsigned char * p3 = src + 2*cols;

  for(j = 0; j < cols; ++j ){
    *(a1 + j) = *(a2 + j) = *(a3 + j) = 0.0;
    *(b1 + j) = *(b2 + j) = *(b3 + j) = 0.0;
    *(b21 + j) = *(b22 + j) = *(b23 + j) = 0.0;
    *(c1 + j) = *(c2 + j) = *(c3 + j) = 0;
  }

  //zero clear
  memset(dst, 0, cols*rows);

  //1st line
  for(i = 1; i < cols - 1; ++i ){

    fx = (*(p2 + i)) * (*skx1 ) +  (*((p2 + i) + 1)) * (*(skx1 + 1)) +  (*((p2 + i) + 2)) * (*(skx1 + 2)) 
      + (*(p2 + i)) * (*skx2 ) +  (*((p2 + i) + 1)) * (*(skx2 + 1)) +  (*((p2 + i) + 2)) * (*(skx2 + 2)) 
      + (*(p3 + i)) * (*skx3 ) +  (*((p3 + i) + 1)) * (*(skx3 + 1)) +  (*((p3 + i) + 2)) * (*(skx3 + 2));

    fy = (*(p2 + i)) * (*sky1 ) +  (*((p2 + i) + 1)) * (*(sky1 + 1)) +  (*((p2 + i) + 2)) * (*(sky1 + 2)) 
      + (*(p2 + i)) * (*sky2 ) +  (*((p2 + i) + 1)) * (*(sky2 + 1)) +  (*((p2 + i) + 2)) * (*(sky2 + 2)) 
      + (*(p3 + i)) * (*sky3 ) +  (*((p3 + i) + 1)) * (*(sky3 + 1)) +  (*((p3 + i) + 2)) * (*(sky3 + 2));

    if( grad == 0 ){
      fx = (fx < 0) ? -fx : fx;
      fy = (fy < 0) ? -fy : fy;
      mag = fx + fy;
    }else{
      mag = sqrt( 1.0*fx*fx + 1.0*fy*fy );
    }

    *(a1 + i) = mag;
    *(b1 + i) = fx;
    *(b21 + i) = fy;
  }

  //2nd line
  for(i = 1; i < cols - 1; ++i ){

    fx = (*(p1 + i)) * (*skx1 ) +  (*((p1 + i) + 1)) * (*(skx1 + 1)) +  (*((p1 + i) + 2)) * (*(skx1 + 2)) 
      + (*(p2 + i)) * (*skx2 ) +  (*((p2 + i) + 1)) * (*(skx2 + 1)) +  (*((p2 + i) + 2)) * (*(skx2 + 2)) 
      + (*(p3 + i)) * (*skx3 ) +  (*((p3 + i) + 1)) * (*(skx3 + 1)) +  (*((p3 + i) + 2)) * (*(skx3 + 2));

    fy = (*(p1 + i)) * (*sky1 ) +  (*((p1 + i) + 1)) * (*(sky1 + 1)) +  (*((p1 + i) + 2)) * (*(sky1 + 2)) 
      + (*(p2 + i)) * (*sky2 ) +  (*((p2 + i) + 1)) * (*(sky2 + 1)) +  (*((p2 + i) + 2)) * (*(sky2 + 2)) 
      + (*(p3 + i)) * (*sky3 ) +  (*((p3 + i) + 1)) * (*(sky3 + 1)) +  (*((p3 + i) + 2)) * (*(sky3 + 2));

    if( grad == 0 ){
      fx = (fx < 0) ? -fx : fx;
      fy = (fy < 0) ? -fy : fy;
      mag = fx + fy;
    }else{
      mag = sqrt( 1.0*fx*fx + 1.0*fy*fy );
    }

    *(a2 + i) = mag;
    *(b2 + i) = fx;
    *(b22 + i) = fy;
  }

  dst += cols;

  //3rd line and the blow
  for( j = 2; j < rows -1 ; ++j ){
    for( i = 1; i < cols - 1; ++i ){

      fx = (*(p1 + i)) * (*skx1 ) +  (*((p1 + i) + 1)) * (*(skx1 + 1)) +  (*((p1 + i) + 2)) * (*(skx1 + 2)) 
    + (*(p2 + i)) * (*skx2 ) +  (*((p2 + i) + 1)) * (*(skx2 + 1)) +  (*((p2 + i) + 2)) * (*(skx2 + 2)) 
    + (*(p3 + i)) * (*skx3 ) +  (*((p3 + i) + 1)) * (*(skx3 + 1)) +  (*((p3 + i) + 2)) * (*(skx3 + 2));

      fy = (*(p1 + i)) * (*sky1 ) +  (*((p1 + i) + 1)) * (*(sky1 + 1)) +  (*((p1 + i) + 2)) * (*(sky1 + 2)) 
    + (*(p2 + i)) * (*sky2 ) +  (*((p2 + i) + 1)) * (*(sky2 + 1)) +  (*((p2 + i) + 2)) * (*(sky2 + 2)) 
    + (*(p3 + i)) * (*sky3 ) +  (*((p3 + i) + 1)) * (*(sky3 + 1)) +  (*((p3 + i) + 2)) * (*(sky3 + 2));

      if( grad == 0 ){
    fx = (fx < 0) ? -fx : fx;
    fy = (fy < 0) ? -fy : fy;
    mag = fx + fy;
      }else{
    mag = sqrt( 1.0*fx*fx + 1.0*fy*fy );
      }

      *(a3 + i) = mag;
      *(b3 + i) = fx;
      *(b23 + i) = fy;

    }// for( i = 1; i < cols - 1; ++i ){

    //judge whether the magnitude is maximum candidate or not
    judge(cols, a1, a2, a3, b2, b22, c2);

    //determine if the candidate is real or not.
    determine(cols, a2, c1, c2, c3, dst, maxth, minth);

    dst += cols;
    p1  += cols;
    p2  += cols;
    p3  += cols;

    //ring buffer
    dtmp = a1;
    a1 = a2;
    a2 = a3;
    a3 = dtmp;

    dtmp = b1;
    b1 = b2;
    b2 = b3;
    b3 = dtmp;

    dtmp = b21;
    b21 = b22;
    b22 = b23;
    b23 = dtmp;

    ctmp = c1;
    c1 = c2;
    c2 = c3;
    c3 = ctmp;

  }// for( j = 2; j < rows -1 ; ++j ){

  //most bottom line
  for(i = 1; i < cols - 1; ++i ){

    fx = (*(p1 + i)) * (*skx1 ) +  (*((p1 + i) + 1)) * (*(skx1 + 1)) +  (*((p1 + i) + 2)) * (*(skx1 + 2)) 
      + (*(p2 + i)) * (*skx2 ) +  (*((p2 + i) + 1)) * (*(skx2 + 1)) +  (*((p2 + i) + 2)) * (*(skx2 + 2)) 
      + (*(p2 + i)) * (*skx3 ) +  (*((p2 + i) + 1)) * (*(skx3 + 1)) +  (*((p2 + i) + 2)) * (*(skx3 + 2));

    fy = (*(p1 + i)) * (*sky1 ) +  (*((p1 + i) + 1)) * (*(sky1 + 1)) +  (*((p1 + i) + 2)) * (*(sky1 + 2)) 
      + (*(p2 + i)) * (*sky2 ) +  (*((p2 + i) + 1)) * (*(sky2 + 1)) +  (*((p2 + i) + 2)) * (*(sky2 + 2)) 
      + (*(p2 + i)) * (*sky3 ) +  (*((p2 + i) + 1)) * (*(sky3 + 1)) +  (*((p2 + i) + 2)) * (*(sky3 + 2));

    if( grad == 0 ){
      fx = (fx < 0) ? -fx : fx;
      fy = (fy < 0) ? -fy : fy;
      mag = fx + fy;
    }else{
      mag = sqrt( 1.0*fx*fx + 1.0*fy*fy );
    }

    *(a3 + i) = mag;
    *(b3 + i) = fx;
    *(b23 + i) = fy;
  }

  //judge whether maximum candidate or not
  judge(cols, a1, a2, a3, b2, b22, c2);

  //determine if the candidate is real or not.
  determine(cols, a2, c1, c2, c3, dst, maxth, minth);

  //release memories
  free( a1 );
  free( a2 );
  free( a3 );

  free( b1 );
  free( b2 );
  free( b3 );

  free( b21 );
  free( b22 );
  free( b23 );

  free( c1 );
  free( c2 );
  free( c3 );

}// void canny(...

void judge(int l, 
       double* a1, double* a2, double* a3, 
       double* b2, double* b22, unsigned char* c2)
{
  int i;
  double mid1, mid2;
  double tg;

  for( i = 1; i < l - 1; ++i ){

    if( b2[i] != 0.0 && b22[i] != 0.0 ){
      tg  =  b22[i] / b2[i];

      if ( tg > 0.0 && tg <= 1.0){
    mid1 = *(a2 + i + 1) + tg * (*(a1 + i + 1) - *(a2 + i + 1));
    mid2 = *(a2 + i - 1) + tg * (*(a3 + i - 1) - *(a2 + i - 1));
    if( *(a2 + i) > mid1 && *(a2 + i) > mid2 )
      *(c2 + i) = 1;
    else
      *(c2 + i) = 0;
      }
      if ( tg > 1.0 ){ 
    tg =  b2[i] / b22[i];
    mid1 = *( a1 + i ) + tg * (*(a1 + i + 1) - *( a1 + i ));
    mid2 = *( a3 + i ) + tg * (*(a3 + i - 1) - *( a3 + i ));
    if( *(a2 + i) > mid1 && *(a2 + i) > mid2 )
      *(c2 + i) = 1;
    else
      *(c2 + i) = 0;
      }
      else if( tg >= -1.0 && tg < 0.0){
    tg =  -1.0 * b22[i] / b2[i];
    mid1 = *( a2 + i - 1 ) + tg * (*(a1 + i - 1) - *( a2 + i - 1 ));
    mid2 = *( a2 + i + 1 ) + tg * (*(a3 + i + 1) - *( a2 + i + 1 ));
    if( *(a2 + i) > mid1 && *(a2 + i) > mid2 )
      *(c2 + i) = 1;
    else
      *(c2 + i) = 0;
      }
      else if ( tg < -1.0 ){ 
    tg =  -1.0 * b2[i] / b22[i];
    mid1 = *( a1 + i ) + tg * (*(a1 + i - 1) - *( a1 + i ));
    mid2 = *( a3 + i ) + tg * (*(a3 + i + 1) - *( a3 + i ));
    if( *(a2 + i) > mid1 && *(a2 + i) > mid2 )
      *(c2 + i) = 1;
    else
      *(c2 + i) = 0;
      }
    }
    else if(b2[i] == 0.0 && b22[i] != 0.0){
      if( *(a2 + i) > *(a1 + i) && *(a2 + i) > *(a3 + i) )
      *(c2 + i) = 1;
      else
    *(c2 + i) = 0;
    }
    else if(b2[i] != 0.0 && b22[i] == 0.0){
      if( *(a2 + i) > *(a2 + i - 1) && *(a2 + i) > *(a2 + i + 1) )
    *(c2 + i) = 1;
      else
    *(c2 + i) = 0;
    }
    else
    *(c2 + i) = 0;
  }

}

void  determine(int l,
        double* a2, 
        unsigned char* c1,
        unsigned char* c2,
        unsigned char* c3, 
        unsigned char* dst, 
        int max,
        int min)
{
  int i;
  int s;

  for( i = 1; i < l - 1; ++i ){
    if ( *(c2 + i) == 1) {
      if ( *(a2 + i) >= max )
    *(dst + i) = 255;
      else if ( *(a2 + i) >= min && *(a2 + i) <= max ){
    s = *(c1 + i - 1) + *(c1 + i) + *(c1 + i + 1)
      + *(c2 + i - 1) + *(c2 + i + 1) +
      + *(c3 + i - 1) + *(c3 + i) + *(c3 + i + 1);
    if ( s != 0 )
      *(dst + i) = 255;
      }
    }
  }
}

