
#include "arm_math.h"
#include <stdlib.h>
#include "math.h"

#include "db4.h"

static float32_t filterHiCoeff [] = {-0.230380, 0.71484, -0.63088, -0.02798, 0.18703, 0.03084, -0.03288, -0.01060}; // High pass filter coefficients of Daubechies 4

static float32_t filterLoCoeff [] = {-0.01060, 0.03288, 0.03084, -0.18703, -0.02798, 0.63088, 0.71484, 0.23038}; //Low pass filter coefficients of Daubechies 4

static int upsampling(const float32_t *, float32_t *, int length, int level);

/*
@Brief: This function calculates detail coefficients of Daubechies 4 discrete
 wavelet decomposition of input at dLevel (decomposition level)  using arm DSP library
@Param:
   input : array of data samples
   coefficients : array of detail coefficients of target octave
   dLevel: the target octave
   inputSize: the size of input array ( must be power of 2 and more than 2^(dLevel + 1) )
@Output: returns -2 if inputs size is short or filter sizes are different, -1 if dynamic memory allocation fails and 0 if result is valid
*/

int db4(float32_t * input, float32_t * coefficients, int dLevel, int inputSize)
{
	if ((pow(2, (dLevel + 1)) > inputSize)||(sizeof(filterHiCoeff)!=sizeof(filterLoCoeff)))
	{
		return -2;
	}	

	int loopCounter = 0;
	
	int filterLength = (int)(sizeof(filterHiCoeff)/sizeof(filterHiCoeff[0])); // filter length

	int filterTaps = filterLength;// Number of filterTaps in the filter

	int UpsampledFilterLength = filterLength ;// Size of upsamplingd filter

	float32_t * pbufLoCoeff = NULL;
	
	float32_t *pbufHiCoeff = NULL; // upsamplingd Low pass and highpass filters

	int neglectedSize=0; // number of neglected elements in convolution result to keep each octave size constant		

	int dataBufSize = inputSize + (pow(2, (dLevel - 1))) * filterLength - 1;
	
	float32_t * pdataBuf = (float32_t *) malloc(sizeof(float32_t) * dataBufSize);	

	if (pdataBuf == NULL)
	{
		return (-1);
	}
	// Compute DWT for each level
	for (loopCounter = 0; loopCounter < dLevel ; loopCounter++)
	{
		// Takes input for the first level decomposition
		if (loopCounter == 0) // first step, don't upsampling the filter
		{			
			arm_conv_f32(input, inputSize, &filterLoCoeff[0], filterTaps, pdataBuf);
		}
		else
		{
			//upsampling low pass filters
			UpsampledFilterLength = 2 * UpsampledFilterLength;

			if (loopCounter == (dLevel - 1)) // Calculate the last step detail coefficients - upsampling high pass filters
			{
				pbufHiCoeff = (float32_t *) malloc(sizeof(float32_t) * UpsampledFilterLength);

				if (pbufHiCoeff == NULL)
				{
					free(pdataBuf);
					pdataBuf = NULL ;
					return (-1);
				}

				upsampling(&(filterHiCoeff [0]), &(pbufHiCoeff[0]), 8, (UpsampledFilterLength / 8));				
				filterTaps = UpsampledFilterLength;
				arm_conv_f32(coefficients, inputSize, &(pbufHiCoeff[0]), filterTaps, pdataBuf);
				free(pbufHiCoeff);
	            pbufHiCoeff = NULL ;
			}
			else
			{
				pbufLoCoeff = (float32_t *) malloc(sizeof(float32_t) * UpsampledFilterLength);

				if (pbufLoCoeff == NULL)
				{
					free(pdataBuf);
					pdataBuf = NULL ;
					return (-1);
				}

				upsampling(&(filterLoCoeff [0]), &(pbufLoCoeff[0]), 8, (UpsampledFilterLength / 8));				
				filterTaps = UpsampledFilterLength;
				arm_conv_f32(coefficients, inputSize, &(pbufLoCoeff[0]), filterTaps, pdataBuf);
				free(pbufLoCoeff);
				pbufLoCoeff = NULL;
			}
		}
		// Store values of dataBuf in coefficients array and keep the size the same as input data
		neglectedSize = (pow(2, loopCounter) * filterLength * 0.5);
		memcpy(&coefficients[0],&pdataBuf[neglectedSize],inputSize);
	}
	// Memory release
	free(pdataBuf);
	pdataBuf = NULL ;
	return (0);
}

/*
@Brief: This function upsamples the wavelet filter coefficients
@Param:
    arr_in:  array of input  (filter coefficients)
    arr_out: array of output (upsampled filter coefficients)
    length:  size of arr_in
    level:   size of arr_out to the size of arr_in that can be 2,4,8,16,...
@Output: returns -1 if arr_in size or level isn't valid, zero if result is valid
*/
int upsampling(const float32_t * arr_in, float32_t * arr_out, int length, int level)
{
	int elementCounter = 0;
	int paddingCounter = 0;
	if ((level < 2) | (length < 1))
	{
		return (-1);
	}
	for (elementCounter = 0; elementCounter < length; elementCounter++)
	{
		arr_out[level * elementCounter] = arr_in[elementCounter];
		for (paddingCounter = 1; paddingCounter < level; paddingCounter++)
		{
			arr_out[level * elementCounter + paddingCounter] = 0;
		}
	}
	return (0);
}
