#include <zephyr/kernel.h>

// Filled automatically by model_exporter.py
// Model ID: 20250630_7

const int num_filters = 4;
const float filters_sos[6] = {};

// Applies all filtering
void apply_lp_filter(float* signal, int nsamp){

}

// Adapted from https://github.com/scipy/scipy/blob/main/scipy/signal/_sosfilt.pyx
void sosfilter(float* signal, int nsamp, float* sos, int nsos){
    // Initialize delay memory bank to 0 (rest)
    // This can be set to specific values, so it could be a function argument, but in our case it'll always start as 0
    float zi[nsos*2] = {0};

    for (int n = 0; n < nsamp; ++n){
        float xn = signal[n];

        for (int s = 0; s < nsos; ++s){
            // Direct II transposed structure
            // See https://www.dsprelated.com/freebooks/filters/Four_Direct_Forms.html for a breakdown of forms
            
            // Note that zi[s*6+3] = a0, and is not typically shown in block diagrams and not used in the calculations below
            // The output node, y[n] would be divided by a0, but a0 is typically 1 (and will always be 1 in our case)
            // because all other coefficients are divided by a0 (rather than applying it to the output node)

            signal[n] = sos[s*6]*xn + zi[s*2];

            zi[s*2] = sos[s*6+1]*xn - sos[s*6+4]*signal[n]+zi[s*2+1];
            zi[s*2+1] = sos[s*6+2]*xn - sos[s*6+5]*signal[n];
        }
    }
}