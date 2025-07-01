#pragma once

#include <fftw3.h>
#include <vector>

#include "constants.h"

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
    std::array<double, FFT_OUT_LENGTH> PSD;

    SongData(double* l, double* r, long int t, int fs)
        : leftChannel(l), rightChannel(r), totalFrames(t), lastFrame(0), soundLevel(0.0), samplerate(fs), PSD({0}){      
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