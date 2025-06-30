#include <array>
#include <chrono>
#include <cmath>

#include <sndfile.h>
#include <fftw3.h>

#include "constants.cpp"

using namespace std;


constexpr array<double, FFT_LENGTH> makeHannWindow() {
    array<double, FFT_LENGTH> window{};
    for (int n = 0; n < window.size(); n++){
        window[n] = 0.5*(1-cos(2*M_PI*n/FFT_LENGTH));
    }
    return window;
}
static const array<double, FFT_LENGTH> hannWindow = makeHannWindow();


void printFileInfo(SF_INFO* info){
    printf("Sample Rate = %d Hz\n", info->samplerate);
    printf("Channels = %d\n", info->channels);
    printf("Format = 0x%x\n", info->format);
    printf("Sections = %d\n", info->sections);
    printf("Seekable = %d\n", info->seekable);
    printf("Frames = %ld\n", info->frames);
}


double powerTodB(double A){
    return 10.0 * log10(A + 1e-10);
}


double clampdB(double db) {
    return min(0.0, max(-DB_LOW, db));
}


double magnitude(fftw_complex c){
    return sqrt(pow(c[0], 2) + pow(c[1], 2));
}


void scale_complex(fftw_complex c, fftw_complex cNew, double scale){
    cNew[0] = c[0]*scale;
    cNew[1] = c[1]*scale;
}


int64_t getTime(){
    auto time = chrono::high_resolution_clock::now();
    auto dur = time.time_since_epoch();
    auto ms = chrono::duration_cast<chrono::milliseconds>(dur).count();
    return ms;
}


vector<int> linspace(double start, double stop, int num) {
    vector<int> res(num);
    double step = (stop - start)/num;
    for (int i = 0; i < num; i++){
        res[i] = start + (i * step);
    }
    return res;
}


double aWeightCurve(double f){
    double f2 = f * f;
    double num = 12194 * 12194 * f2 * f2;
    double den = (f2 + 20.6 * 20.6) * sqrt((f2 + 107.7 * 107.7) * (f2 + 737.9 * 737.9)) * (f2 + 12194 * 12194);
    return num/den;
}

double aWeightingdB(double f){
    return 20 * log10(aWeightCurve(f)) + 2.0;
}


vector<pair<int, int>> createOctaveBands(int sampleRate) {
    double fMax = sampleRate / 2.0;

    // Create logarithmically spaced frequency edges
    vector<double> freqEdges(REQUESTED_NUMBER_OF_POINTS + 1);
    for (int i = 0; i <= REQUESTED_NUMBER_OF_POINTS; ++i) {
        double fraction = static_cast<double>(i) / REQUESTED_NUMBER_OF_POINTS;
        freqEdges[i] = 20.0 * pow(fMax / 20.0, fraction);
    }

    vector<pair<int, int>> bands;
    for (int i = 0; i < REQUESTED_NUMBER_OF_POINTS; ++i) {
        double fStart = freqEdges[i];
        double fEnd = freqEdges[i + 1];

        int binStart = static_cast<int>(ceil(fStart * FFT_LENGTH / sampleRate));
        int binEnd = static_cast<int>(ceil(fEnd * REQUESTED_NUMBER_OF_POINTS / sampleRate));

        // Clamp to FFT range
        binStart = clamp(binStart, 0, FFT_OUT_LENGTH - 1);
        binEnd = clamp(binEnd, binStart + 1, FFT_OUT_LENGTH);

        bands.emplace_back(binStart, binEnd);
    }

    return bands;
}

// vector<pair<int,int>> createOctaveBands(int samplerate){
//     int fftBins = FFT_LENGTH / 2;
//     double f_min = 20.0;
//     double f_max = samplerate / 2;
//     double ratio = pow(f_max / f_min, 1.0 / REQUESTED_NUMBER_OF_POINTS);

//     vector<pair<int, int>> bands;
//     double prev = f_min;
//     for (int i = 0; i < REQUESTED_NUMBER_OF_POINTS; ++i){
//         double next = prev * ratio;
//         int startBin = ceil(prev * FFT_LENGTH / samplerate);
//         int endBin = floor(next * FFT_LENGTH / samplerate);
//         bands.push_back({startBin, endBin});
//         prev = next;
//     }
//     return bands;
// }




typedef struct SongData {
    double* leftChannel;
    double* rightChannel;
    long int totalFrames;
    long int lastFrame;
    int samplerate;

    double* leftSignalBuffer;
    double* rightSignalBuffer;

    fftw_complex* leftOutBuffer;
    fftw_complex* rightOutBuffer;

    fftw_plan leftPlan;
    fftw_plan rightPlan;

    double soundLevel;
    std::vector<double> levels;

    SongData(double* l, double* r, long int t, int fs)
        : leftChannel(l), rightChannel(r), totalFrames(t), lastFrame(0), soundLevel(0.0), samplerate(fs), levels(FFT_LENGTH){      
        leftSignalBuffer = fftw_alloc_real(sizeof(double) * FFT_LENGTH);
        rightSignalBuffer = fftw_alloc_real(sizeof(double) * FFT_LENGTH);
        
        leftOutBuffer = fftw_alloc_complex(sizeof(fftw_complex) * FFT_OUT_LENGTH);
        rightOutBuffer = fftw_alloc_complex(sizeof(fftw_complex) * FFT_OUT_LENGTH);

        leftPlan = fftw_plan_dft_r2c_1d(
            FFT_LENGTH,
            leftSignalBuffer,
            leftOutBuffer,
            FFTW_ESTIMATE
        );

        rightPlan = fftw_plan_dft_r2c_1d(
            FFT_LENGTH,
            rightSignalBuffer,
            rightOutBuffer,
            FFTW_ESTIMATE
        );
    }

    ~SongData() {
        fftw_free(leftOutBuffer);
        fftw_free(rightOutBuffer);
        fftw_free(leftSignalBuffer);
        fftw_free(rightSignalBuffer);
        fftw_destroy_plan(rightPlan);
        fftw_destroy_plan(leftPlan);
    }
} SongData;

