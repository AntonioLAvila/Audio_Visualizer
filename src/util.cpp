#include <array>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <sndfile.h>
#include <fftw3.h>

#include "constants.h"
#include "util.h"

using namespace std;


constexpr array<double, FFT_LENGTH> makeHannWindow() {
    array<double, FFT_LENGTH> window{};
    for (int n = 0; n < window.size(); n++){
        window[n] = 0.5*(1-cos(2*M_PI*n/FFT_LENGTH));
    }
    return window;
}
const array<double, FFT_LENGTH> hannWindow = makeHannWindow();
const double hannWindowSum = 0.5 * FFT_LENGTH;


const array<double, 7> bhCoefs = {0.27105140069342, -0.43329793923448, 0.21812299954311, -0.06592544638803, 0.01081174209837, -0.00077658482522, 0.00001388721735};
constexpr array<double, FFT_LENGTH> makeBlackmanHarrisWindow() {
    array<double, FFT_LENGTH> window = {0};
    for (int n = 0; n < FFT_LENGTH; n++){
        for (int k = 0; k < bhCoefs.size(); k++){
            window[n] += bhCoefs[k] * (2*M_PI*k*n / FFT_LENGTH);
        }
    }
    return window;
}
const array<double, FFT_LENGTH> blackmanHarrisWindow = makeBlackmanHarrisWindow();
const double blackmanHarrisWindowSum = reduce(blackmanHarrisWindow.begin(), blackmanHarrisWindow.end(), 0);


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
    return std::min(0.0, std::max(-DB_LOW, db));
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
    double den = (f2 + 20.6 * 20.6) * (f2 + 12194 * 12194) * sqrt((f2 + 107.7 * 107.7) * (f2 + 737.9 * 737.9));
    return num/den;
}

double aWeightdB(double f){
    return 20 * log10(aWeightCurve(f)) + 2.0;
}


const array<double, 32> nominalOneThirdOctaveFrequencies = \
    {16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000};


vector<pair<int, int>> createOctaveBands(int fs){
    vector<pair<int, int>> bounds;
    static double fd = pow(10, 0.05);
    for (const double &fCenter : nominalOneThirdOctaveFrequencies){
        double fLower = fCenter / fd;
        double fUpper = fCenter * fd;

        int kLower = fLower * FFT_LENGTH / fs;
        int kUpper = fUpper * FFT_LENGTH / fs;

        bounds.emplace_back(kLower, kUpper);
    }
    return bounds;
}

