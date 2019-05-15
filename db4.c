
#include "arm_math.h"
#include <stdlib.h>
#include "math.h"

#include "db4.h"



static int upsampling(const float32_t*,float32_t*,int,int);

/*
@Brief: This function calculates detail coefficients of Daubechies 4 discrete
 wavelet decomposition of input at d_level  using arm DSP library(size of input shoud be 2^n)
 @Param: input is array of data samples, Coefficients is array of detail coefficients,
 d_level is the target octave and input_size is the size of input array
 @output: returns -2 if inputs size are incorrect, -1 if memory allocation
 fails and 0 if result is valid
*/

int db4( float32_t *Input , float32_t * Coefficients, int d_level,int input_size)
{  
  if(pow(2,(d_level+1))>input_size)
	  return -2;
	  
  int i,j;    
  
  float32_t Hi_D [ FILTER_LENGTH ]={-0.230380,0.71484,-0.63088,-0.02798,0.18703,0.03084,-0.03288,-0.01060};// Hi pass filter Coefficients of Daubechies 4
  
  float32_t Lo_D [ FILTER_LENGTH ]={-0.01060,0.03288,0.03084,-0.18703,-0.02798,0.63088,0.71484,0.23038}; //Low pass filter Coefficients of Daubechies 4
  
  int filter_taps = FILTER_LENGTH;// Number of filter_taps in the filter
  
  int u_filter_size = FILTER_LENGTH ;// Size of upsamplingd filter
  
  float32_t *U_Lo_D,*U_Hi_D;// upsamplingd Low pass and highpass filters

  int num_e; // number of neglected elements in convolution result to keep each octave size constant 
  
  float32_t *lo , *hi , *buf ;// Pointers to filters , convolution and coefficients arrays
   
  int level_size;
   
  level_size = input_size + (pow(2,(d_level-1)))*FILTER_LENGTH - 1; 
  
  // Computes DWT for each level
  buf = ( float32_t *) malloc ( sizeof ( float32_t )* level_size);
  
  if(buf==NULL){
     return(-1);
    }
  for (i = 0; i < (d_level-1) ; i++)
  {
    // Takes input for the first level decomposition
    if(i == 0)// first step, don't upsampling the filter 
    {
     lo = & Lo_D [0];
     arm_conv_f32 (Input , input_size , lo , filter_taps , buf);         
    }    
    else
    {
      //upsampling low pass filters
      u_filter_size=2*u_filter_size;
      if(i>1){
      free (U_Lo_D);
      }
      U_Lo_D = ( float32_t *) malloc ( sizeof ( float32_t )* u_filter_size);
	  if(U_Lo_D==NULL){
		  free (buf); buf = NULL ;
          return(-1);
      }
	  
      upsampling(&(Lo_D [0]),&(U_Lo_D[0]),8,(u_filter_size/8));
      lo=&(U_Lo_D[0]);
      filter_taps=u_filter_size;
      arm_conv_f32 (Coefficients , input_size , lo , filter_taps , buf);     
     }
      // Store values of buf in Coefficients array,
      num_e=(pow(2,i)*FILTER_LENGTH*0.5);                          // repeated code!!!!!!!!!!!!!!!!!!!!!
      for (j =num_e; j <(input_size+num_e); j++){
       Coefficients[j-num_e]=buf[j];         
    }
  }
    // Calculate the last step detail coefficients
      
    //upsampling high pass filters
      u_filter_size=2*u_filter_size;
      if(i>1){
      free (U_Lo_D);
      }
      U_Hi_D = ( float32_t *) malloc ( sizeof ( float32_t )* u_filter_size);
	  if(U_Hi_D==NULL){
		  free (buf); buf = NULL ;
          return(-1);
      }
      upsampling(&(Hi_D [0]),&(U_Hi_D[0]),8,(u_filter_size/8));
      hi=&(U_Hi_D[0]);
      filter_taps=u_filter_size;
      arm_conv_f32 (Coefficients , input_size , hi , filter_taps , buf);
      num_e=(pow(2,i)*FILTER_LENGTH*0.5);
     for (j =num_e; j <(input_size+num_e); j++){
       Coefficients[j-num_e]=buf[j];    
     } 
    // Memory release
    free (buf); buf = NULL ;
    free (U_Hi_D); U_Hi_D = NULL ;    
  return(0);
}

/*@Brief: This function upsamplings the wavelet filter coefficients
  @Param: input is array of filter coefficients, output is array of upsamplingd 
          filter coefficients and level is size of output to the size of input that can be 2,4,8,16,...
  @output: returns -1 if input size isn't valid, zero if result is valid 
*/
int upsampling(const float32_t* input,float32_t* output,int length,int level)
{
 int i,j;
 if ((level<2)|(length<1)){
 return(-1);
 }
 for(i=0;i<length;i++){
   output[level*i]=input[i];
   for (j=1;j<level;j++){   
   output[level*i+j]=0;   
   }  
 }
 return(0);
}



