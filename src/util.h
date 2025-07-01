#pragma once

#include <vector>
#include <array>
#include <sndfile.h>
#include <fftw3.h>
#include <chrono>
#include <numeric>

#include "constants.h"

using namespace std;

extern const array<double, FFT_LENGTH> hannWindow;
extern const double hannWindowSum;

extern const array<double, FFT_LENGTH> blackmanHarrisWindow;
extern const double blackmanHarrisWindowSum;

void printFileInfo(SF_INFO* info);
double powerTodB(double A);
double clampdB(double db);
void scale_complex(fftw_complex c, fftw_complex cNew, double scale);
int64_t getTime();
vector<int> linspace(double start, double stop, int num);
double aWeightCurve(double f);
double aWeightdB(double f);

vector<pair<int, int>> createOctaveBands(int sampleRate);
