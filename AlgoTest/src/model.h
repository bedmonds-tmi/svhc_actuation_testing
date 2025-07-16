#include <zephyr/kernel.h>

// Filled automatically by model_exporter.py
// Model ID: 20250630_7

#define WINDOW_SIZE 8000
#define STEP_SIZE 2000

const char features[15][30] = {"spectral_centroid*rms", "rms", "delta spectral_centroid*rms", "delta rms", "mfcc_4", "mfcc_1", "delta fft_peak3", "delta mfcc_avg", "spectral_spread", "delta fft_peak1", "fft_peak3", "delta mfcc_1", "delta fft_peak2", "fft_peak2", "spectral_centroid"};
const int dimensions = 15;

const float weight[15] = {0.8253315129401085, 0.0, 0.22979798817980895, 0.0, 0.6215227378189154, 0.0, 0.0, 0.0, 0.5582961542879631, 0.0, 0.0, 0.003473601308559734, 0.0, -0.018791246666432837, 0.09918869162966702};
const float intercept = -1.6850998917180784;

const float means[15] = {271663.0498156863, 193.84623917136557, 41682.912324814046, 26.89486044602858, -0.6284134519704544, -11.781066533480727, 3.9666141199505294, 0.011705385663154738, 1358.2522773602589, 4.806058309263547, 39.08466752872953, 0.1485788200470457, 4.270694942487929, 43.22207904966284, 2907.6153267588093};
const float stds[15] = {268354.34381315, 178.4971628971753, 142379.48607752298, 92.99485437344127, 0.7049277102837309, 4.336151147296494, 14.666780793550243, 0.05104512850593677, 111.47925472898721, 19.800468486863736, 33.183246676522835, 0.7389947658972233, 16.44520534965847, 36.77710153543518, 618.2460558786829};

int predict_linearSVM(float *vector){
    /*
    Produces a prediction using a linear SVM based on a single window feature-vector.
    
    Parameters:
        *vector: an array of non-standardized features for a single time-window, length {dimensions}
    Returns:
        classification: 1 (actuation) or -1 (not actuation)
    */

    int value = 0;

    for(int i=0; i < 15; i++){
        vector[i] = (vector[i] - means[i]) / stds[i];
        value += weight[i] * vector[i];
    }

    value = value - intercept;

    if (value>0){
        return 1;
    }
    else{
        return -1;
    }
}