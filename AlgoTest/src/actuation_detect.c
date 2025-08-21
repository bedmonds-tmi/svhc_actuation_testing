#include "actuation_detect.h"

void predict_actuation(const int16_t* audio_data, const int audio_data_cnt, const float* pressure_data){
    printk("Test\n");
    printk("%d\n", audio_data_cnt);
    // Error check
    if(audio_data_cnt != audio_window_length){
        return;
    }

    // When the first half-window of data is received store and do nothing
    if(first_window == 1){
        printk("First window received\n");
        first_window = 0;
        for(int i = 0; i < audio_data_cnt; ++i){
            prev_audio_window[i] = audio_data[i];
        }
        for(int i = 0; i < pressure_window_length; ++i){
            prev_pressure_window[i] = (float32_t)pressure_data[i];
        }

        return;
    }

    printk("Merging windows; ");
    int64_t us1 = k_uptime_ticks() * 1000000LL / CONFIG_SYS_CLOCK_TICKS_PER_SEC;
    // Merge previous half-window and new half window
    float32_t audio_window[audio_window_length*2];  // audio_window_length is num samples AFTER down sampling, mult by 2 to get current window size
    for(int i = 0; i < audio_window_length; ++i){
        audio_window[i] = prev_audio_window[i];
        audio_window[i+audio_window_length-1] = (float32_t)audio_data[i];
        prev_audio_window[i] = (float32_t)audio_data[i];
    }
    float32_t pressure_window[pressure_window_length];
    for(int i = 0; i < (int)(pressure_window_length/2); ++i){
        pressure_window[i] = prev_pressure_window[i];
        pressure_window[i+(int)(pressure_window_length/2)-1] = (float32_t)pressure_data[i];
    }
    int64_t us2 = k_uptime_ticks() * 1000000LL / CONFIG_SYS_CLOCK_TICKS_PER_SEC;
    printk("Took %lld us\n", us2-us1);

    // This should go somewhere else
    int audio_sample_rate = 8000;

    // Apply low-pass filter to audio data
    /*
    *   void arm_biquad_cascade_df2T_f32 	( 	const arm_biquad_cascade_df2T_instance_f32 *  	S,
    *       const float32_t *  	pSrc,
    *       float32_t *  	pDst,
    *       uint32_t  	blockSize 
    *   ) 	
    */
    printk("Audio filtering\n");
    float32_t audio_data_filt[audio_window_length*2];
    int audio_filter_stages = 1;
    float32_t p_state[audio_filter_stages];
    arm_biquad_cascade_df2T_instance_f32 S1 = {audio_filter_stages, p_state, audio_filter_sos};
    arm_biquad_cascade_df2T_f32(&S1, audio_window, audio_data_filt, audio_window_length*2);

    printk("Audio downsampling\n");
    // Down sample audio data
    float audio_window_dwnsmpl[audio_window_length];
    for (int i = 0; i < audio_window_length*2; i++) {
      audio_window_dwnsmpl[i] = audio_data_filt[i*2];
    }

    // Apply hanning window to audio and pressure data
    float32_t audio_signal[audio_window_length];
    arm_mult_f32(audio_window_dwnsmpl, audio_hanning_window, audio_signal, audio_window_length);

    float32_t pressure_signal[pressure_window_length];
    arm_mult_f32(pressure_window, pressure_hanning_window, pressure_signal, pressure_window_length);

    printk("FFT\n");
    // Audio FFT calculation
    int num_fft = (int)(audio_data_cnt/2);
    float32_t fft_raw[2*num_fft];
    arm_rfft_fast_instance_f32 rfft;
    arm_rfft_fast_init_f32(&rfft, num_fft);
    arm_rfft_fast_f32(&rfft, audio_signal, fft_raw, 0);

    // Get audio FFT magnitude
    float32_t fft_mag[num_fft];
    arm_cmplx_mag_f32(fft_raw, fft_mag, num_fft);

    // Initialize audio features
    float power = 0;
    float centroid = 0;

    float power_band1 = 0;
    float power_band2 = 0;

    // Numerator and denominator for spectral centroid
    float numerator = 0;
    float total_audio_power = epsilon;
    
    int min_idx = (int)(800*num_fft/audio_sample_rate);   // hacky high pass filter

    for (int i = min_idx; i < num_fft; ++i){
        numerator += audio_fft_freqs[i]*fft_mag[i];
        total_audio_power += fft_mag[i];

        power_band1 += fft_mag[i]*audio_trifils[2][i];
        power_band2 += fft_mag[i]*audio_trifils[3][i];
    }

    centroid = numerator / total_audio_power;
    float pow_centroid = centroid*power;

    // Pressure data filter
    // Apply low-pass filter to audio data
    /*
    *   void arm_biquad_cascade_df2T_f32 	( 	const arm_biquad_cascade_df2T_instance_f32 *  	S,
    *       const float32_t *  	pSrc,
    *       float32_t *  	pDst,
    *       uint32_t  	blockSize 
    *   ) 	
    */
    float32_t pressure_data_filt[pressure_window_length];
    int pressure_filter_stages = 2;
    float32_t p_state1[pressure_filter_stages];
    arm_biquad_cascade_df2T_instance_f32 S2 = {pressure_filter_stages, p_state1, pressure_filter_sos};
    arm_biquad_cascade_df2T_f32(&S2, pressure_signal, pressure_data_filt, pressure_window_length);

    // Pressure RMS
    float32_t pressure_rms;
    arm_rms_f32(pressure_data_filt, (uint32_t)pressure_window_length, &pressure_rms);

    float interaction1 = (float)(pressure_rms*total_audio_power);
    float interaction2 = (float)(pressure_rms*centroid);
    float interaction3 = (float)(pressure_rms*interaction1);

    float delta_pressure_rms = (float)(pressure_rms - prev_pressure_rms);
    float feature_vector[9] = {(float)pressure_rms, delta_pressure_rms, 
                                total_audio_power, power_band1, power_band2, 
                                pow_centroid, 
                                interaction1, interaction2, interaction3};

    printk("%d\n", model_predict(feature_vector));
    return;
}