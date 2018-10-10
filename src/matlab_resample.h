//RESAMPLE  Change the sampling rate of a signal.
//   Y = RESAMPLE(UpFactor, DownFactor, InputSignal, OutputSignal) resamples the sequence in 
//   vector InputSignal at UpFactor/DownFactor times and stores the resampled data to OutputSignal.
//   OutputSignal is UpFactor/DownFactor times the length of InputSignal. UpFactor and DownFactor must be 
//   positive integers.

//This function is translated from Matlab's Resample funtion. 

//Author: Haoqi Bai

//modified by Miha Mohorčič (move to namespace, use std17 instead of boost,..)

#ifndef VISECG2_MATLAB_RESAMPLE_H
#define VISECG2_MATLAB_RESAMPLE_H

#include <stdint.h>
#include <vector>

namespace ve::dsp::resample {
  void matlab_resample(uint32_t upFactor, uint32_t downFactor,
              std::vector<double> &inputSignal, std::vector<double> &outputSignal);
}
#endif