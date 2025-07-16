#pragma once

#include <math.h>
#include <zephyr/kernel.h>
#include "arm_math.h"
#include "arm_sorting.h"
#include "filter.h"
#include "merge.h"

// Smallest possible value (prevent divide by zero errors)
const float epsilon = 1/(2^52);

// MFCC consts, these should be function parameters but const for us
const float lowfreq = 133.33;
const float linc = 200/3;
const float logsc = 1.0711703;
const int num_lin_filt = 13;
const int num_log_filt = 27;
const int num_filt_total = num_lin_filt + num_log_filt;
const int freq_size = num_filt_total + 2;

// Delta tracking
float* current_window;
int first_window = 1;

// Main feature extraction function
void feature_extraction(const int* window, const int window_size, const int sampling_rate, const float** mfcc_fbank, float** output);

// Individual feature functions (called internally)
void zero_crossing_rate(const float frame[], const int size, float* zcr);
void rms(const float frame[], int size, float* rms);
void spectral_centroid_spread(const float fft_magnitude[], const float f[], const int size, float* centroid, float* spread);
void mfcc_filter_banks(const int sampling_rate, const int num_fft, float** frequencies, float*** fbank);
void mfcc(const float fft_magnitude[], const float** fbank, const int fft_mag_size, const int fbank_size, const int num_mfcc_feats, float** output);
void fft_peaks(const float fft_magnitude[], const float fft_freq[], const int num_fft, float** peaks, float** peak_freqs);

// Helper functions
float fpow(float val, const int pow);
void dct_type_2_ortho(const float input[], float** output, int n);