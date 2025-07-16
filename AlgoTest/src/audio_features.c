#include "audio_features.h"

void zero_crossing_rate(const float frame[], const int size, float* zcr){
    // TODO
}

void rms(const float frame[], int size, float* rms){
    float tmp = 0;
    for(int i = 0; i < size; ++i){
        tmp += frame[i] * frame[i];
    }
    tmp /= size;
    tmp = sqrtf(tmp);

    *rms = tmp;
}

void spectral_centroid_spread(const float fft_magnitude[], const float f[], const int size, float* centroid, float* spread){
    float numerator = 0;
    float denominator = epsilon;
    for(int i = 0; i < size; ++i){
        numerator += f[i]*fft_magnitude[i];
        denominator += fft_magnitude[i];
    }

    *centroid = numerator / denominator;
    float sum = 0;
    for(int i = 0; i < size; ++i){
        sum += fft_magnitude[i]*fpow((f[i] - *centroid), 2);
    }

    *spread = sqrtf(sum / denominator);
}

float fpow(float val, const int pow){
    // Should use an external library for this
    // No special case consideration, no overflow protection
    for(int i = 0; i < pow; ++i){
        val *= val;
    }
    
    return val;
}

void mfcc_filter_banks(const int sampling_rate, const int num_fft, float** frequencies, float*** fbank){
    /* 
    MFCC Calculations
    Parameters:
        sampling_rate: 
        num_fft: the number of points in the fft
    Return variables (by address):
        frequencies: the centre frequencies in the filter bank
        fbank: the filter bank, 2D array of shape (filters, num_fft)
            - each filter in the MFCC filter bank is triangular with a response of 1 
            at its center frequency (stored in frequencies) and decreasing linearly towards 
            0 until the next frequency
    Adapted from https://pypi.python.org/pypi/scikits.talkbox 
    */

    float heights[num_filt_total] = {0};
    float nfreqs[num_fft] = {0};
    for(int i = 0; i < num_fft; ++i){
        nfreqs[num_fft] = ((float)i / (float)num_fft) * (float)sampling_rate;
    }
    
    // Allocate return vars
    // I hate this
    /*NOTE: I don't currently use the frequencies array for anything outside of this function, so we could
    just return the filter bank*/
    *frequencies = (float*)malloc(freq_size*sizeof(float));
    *fbank = (float**)malloc(num_filt_total*sizeof(float));
    for (int i = 0; i < num_filt_total; ++i) {
        *fbank[i] = (float *)malloc(num_fft * sizeof(float));
    }

    // Compute centre frequencies
    for(int i = 0; i < num_lin_filt; ++i){
        *frequencies[i] = lowfreq + (float)i*linc;
    }
    for(int i = num_lin_filt; i < freq_size; ++i){
        *frequencies[i] = *frequencies[num_lin_filt-1] * fpow(logsc, i-num_lin_filt+1);
    }

    // Compute the y part of the slopes
    for(int i = 0; i < num_filt_total; ++i){
        heights[i] = 2 / (*frequencies[i+2] - *frequencies[i]);
    }

    // Compute filterbank coefficients
    for(int i = 0; i < num_filt_total; ++i){
        float low_freqs  = *frequencies[i];
        float cent_freqs = *frequencies[i + 1]; // The centre frequency for the ith filter bank
        float high_freqs = *frequencies[i + 2];

        int low_id = (int)(low_freqs*num_fft/sampling_rate) + 1;    // The idx (ie the fbank column) where the triangle begins
        int cent_id = (int)(cent_freqs*num_fft/sampling_rate) + 1;  // The idx where the triangle is at its peak
        int high_id = (int)(high_freqs*num_fft/sampling_rate) + 1;  // The idx where the triangle ends

        float lslope = heights[i] / (cent_freqs - low_freqs);
        float rslope = heights[i] / (high_freqs - cent_freqs);

        // NOTE: it would be more efficient to ONLY loop through the part of fbank that is non-zero,
        // however, malloc does not necessary initialize the fbank arrays to 0
        // Should see if it is possible to initialize arrays to 0
        for(int j = 0; j < num_fft; ++j){
            if(low_id <= j && j < cent_id){
                *fbank[i][j] = lslope * (nfreqs[j] - low_freqs);
            }
            else if(cent_id <= j && j < high_id){
                *fbank[i][j] = rslope * (high_freqs - nfreqs[j]);
            }
            else{
                *fbank[i][j] = 0;
            }
        }
    }
}

void dct_type_2_ortho(const float input[], float** output, int n) {
    int fft_size = n * 2;

    float x_ext[fft_size];
    float fft_out[fft_size];

    *output = (float*)malloc(n*sizeof(float));
    
    // Create even-symmetric extension of the input
    for (int i = 0; i < n; i++) {
        x_ext[i] = input[i];
        x_ext[fft_size - 1 - i] = input[i];  // symmetric
    }

    // Compute real FFT of x_ext (length 2N)
    arm_rfft_fast_instance_f32 rfft;
    arm_rfft_fast_init_f32(&rfft, fft_size);
    arm_rfft_fast_f32(&rfft, x_ext, fft_out, 0);

    // DCT-II from FFT real part with orthonormal scaling
    // Scaling defined in https://docs.scipy.org/doc/scipy/tutorial/fft.html#type-ii-dct
    const float scale = sqrtf(1/(2 * (float)n));
    for (int k = 0; k < n; k++) {
        float re = fft_out[2 * k];
        float im = fft_out[2 * k + 1];

        // Compute twiddle factor: exp(-j * π * k / (2N))
        float theta = (PI * k) / (2.0f * N);
        float cos_theta = cosf(theta);
        float sin_theta = sinf(theta);

        // Real part of FFT[k] * twiddle[k]
        // NOTE: you can prove this by doing a binomial expansion between FFT[k] and the twiddle factor
        float dct_val = re * cos_theta + im * sin_theta;

        // Apply orthonormal normalization
        if (k == 0) {
            output[k] = dct_val * sqrtf(1.0f / (4.0f * N));
        } else {
            output[k] = dct_val * sqrtf(1.0f / (2.0f * N));
        }
    }
}

void mfcc(const float fft_magnitude[], const float** fbank, const int fft_mag_size, const int fbank_size, const int num_mfcc_feats, float** output){
    // Get correlated cepstral coefficients
    float result[fbank_size] = {0};
    for(int i = 0; i < fbank_size; ++i){
        for(int j = 0; j < fft_mag_size; ++j){
            result[i] += fft_magnitude[j]*fbank[i][j];
        }
    }

    // Decorrelate the cepstral coefficients using DCT
    float* dct;
    dct_type_2_ortho(result, dct, fbank_size);

    // Get only the first num_mfcc_feats coefficients, typically 13
    *output = (float*)malloc(num_mfcc_feats*sizeof(float));
    for(int i = 0; i < num_mfcc_feats; ++i){
        output[i] = dct[i];
    }
}

void fft_peaks(const float fft_magnitude[], const float fft_freq[], const int num_fft, const int num_peaks, float** peaks, float** peak_freqs){
    // Sort fft magnitude
    // Create an Element array to store fft mag values and the indices of each mag
    Element arr[num_fft];
    for(int i = 0; i < num_fft; ++i){
        arr[num_fft].value = fft_magnitude[i];
        arr[num_fft].original_index = i;
    }

    // Returns sorted array (of type Element) in descending order
    merge_sort(arr, 0, num_fft-1);
    *peaks = (float*)malloc(num_peaks*sizeof(float));
    *peak_freqs = (float*)malloc(num_peaks*sizeof(float));

    // Get fft peaks and associated frequencies
    for(int i = 0; i < num_peaks; ++i){
        peaks[i] = arr[i].value;
        peak_freqs[i] = fft_freq[arr[i].original_index];
    }
}

void feature_extraction(const int* window, const int window_size, const int sampling_rate, const float** mfcc_fbank, float** output){
    // Convert to float
    float window_flt[window_size] = {0};
    for(int i = 0; i < window_size; ++i){
        window_flt[i] = (float)window[i];
    }

    // Apply low-pass filter
    apply_lp_filter(window_flt, window_size);
    uint8_t numStages = 1;
    float32_t pState[1] = {0};
    float32_t pCoeffs[1] = {0};
    arm_biquad_cascade_df2T_instance_f64 S1 = {numStages, pState, pCoeffs};

    // Down sample
    float window2[window_size/2] = {0};
    for (int i = 0; i < window_size; i++) {
      window2[i] = window_flt[i*2];
    }

    // Apply hanning window
    float32_t window[window_size];
    float32_t signal[window_size];
    arm_window_hanning_f32(window, window_size);
    arm_mult_f32(window_flt, window, signal, window_size);

    // FFT calculation
    int num_fft = (int)(window_size/2);
    float32_t fft_out[num_fft];
    arm_rfft_fast_instance_f32 rfft;
    arm_rfft_fast_init_f32(&rfft, num_fft);
    arm_rfft_fast_f32(&rfft, signal, fft_out, 0);

    // Get FFT magnitude
    float32_t fft_mag[num_fft];
    // arm_cmplx_mag_f32(fft_out, )

    // Extract features
    int num_mfccs = 13;
    float* mfccs = NULL;
    mfcc(fft_out, mfcc_fbank, num_fft, )

}