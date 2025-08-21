#pragma once

#include <math.h>
#include <zephyr/kernel.h>
#include "arm_math.h"
#include "model.h"

// Smallest possible value (prevent divide by zero errors)
const float epsilon = 0.0000000000000000000000000000000000000000001;

// Tracking of prev window data
float32_t prev_audio_window[2048];
float32_t prev_pressure_window[32];
int first_window = 1;

float32_t prev_pressure_rms = 0;

// Main feature extraction function
void predict_actuation(const int16_t* audio_data, const int audio_window_size, const float* pressure_data);