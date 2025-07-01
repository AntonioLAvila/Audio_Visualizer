#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <vector>
#include <math.h>
#include <ranges>

#include <portaudio.h>
#include <sndfile.h>
#include <SFML/Graphics.hpp>
#include <fftw3.h>

#include "constants.h"
#include "util.h"
#include "songdata.h"

using namespace std;

static int audioCallback(
    const void* inputBuffer,
    void* outputBuffer,
    unsigned long framesPerBuffer,
    const PaStreamCallbackTimeInfo* timeInfo,
    PaStreamCallbackFlags statusFlags,
    void* userData
){
    SongData* song = (SongData*)userData;
    float* out = (float*)outputBuffer;
    
    long int fftIdx;
    if (song->lastFrame + FFT_LENGTH >= song->totalFrames) {
        fftIdx = song->totalFrames - FFT_LENGTH;
    } else {
        fftIdx = song->lastFrame;
    }

    double sum = 0;
    int frameWriteCount = 0;
    for (int n = 0; n < FFT_LENGTH; n++){
        // write out song
        if (n < framesPerBuffer){
            *out++ = song->leftChannel[song->lastFrame];
            *out++ = song->rightChannel[song->lastFrame];
            double powerLeft = song->leftChannel[song->lastFrame] * song->leftChannel[song->lastFrame];
            double powerRight = song->rightChannel[song->lastFrame] * song->rightChannel[song->lastFrame];
            sum += (powerLeft + powerRight)/2;
            frameWriteCount += 1;
        }
        // window
        song->leftSignalBuffer[n] = song->leftChannel[fftIdx+n] * blackmanHarrisWindow[n];
        song->rightSignalBuffer[n] = song->rightChannel[fftIdx+n] * blackmanHarrisWindow[n];
    }
    // write out rest of song
    while (frameWriteCount < framesPerBuffer){
        *out++ = song->leftChannel[song->lastFrame];
        *out++ = song->rightChannel[song->lastFrame];
        double powerLeft = song->leftChannel[song->lastFrame] * song->leftChannel[song->lastFrame];
        double powerRight = song->rightChannel[song->lastFrame] * song->rightChannel[song->lastFrame];
        sum += (powerLeft + powerRight)/2;
        frameWriteCount += 1;
    }
    song->lastFrame += frameWriteCount;
    song->soundLevel = sum/framesPerBuffer;

    // fft
    fftw_execute(song->leftPlan); fftw_execute(song->rightPlan);

    // calc power
    for (int i = 0; i < FFT_OUT_LENGTH; i++){
        fftw_complex leftW, rightW;
        leftW[0] = song->leftOutBuffer[i][0] / blackmanHarrisWindowSum;
        leftW[1] = song->leftOutBuffer[i][1] / blackmanHarrisWindowSum;
        rightW[0] = song->rightOutBuffer[i][0] / blackmanHarrisWindowSum;
        rightW[1] = song->rightOutBuffer[i][1] / blackmanHarrisWindowSum;

        double powerLeft = leftW[0]*leftW[0] + leftW[1]*leftW[1];
        double powerRight = rightW[0]*rightW[0] + rightW[1]*rightW[1];

        if (0 < i && i < FFT_OUT_LENGTH - 1) {
            powerLeft *= 2;  // Account for symmetric frequencies
            powerRight *= 2;
        }

        double f = i * song->samplerate / FFT_LENGTH;
        double aWeight = aWeightCurve(f);
        powerLeft *= aWeight;
        powerRight *= aWeight;

        double power = (powerLeft + powerRight)/2;
        song->levels[i] = power;
    }
    return 0;
}


class BarVizualizer{
    private:
        vector<sf::RectangleShape> bars;
        vector<double> heights;
        double barWidth;
        double numBars;
        sf::RenderWindow* window;
        static constexpr float alpha = 0.8;
    
    public:
        BarVizualizer(sf::RenderWindow* win, int n){
            barWidth = SCREEN_WIDTH/n;
            numBars = n;
            window = win;
            bars = vector<sf::RectangleShape>(numBars);
            heights = vector<double>(numBars, 0.f);
            for (int i = 0; i < numBars; i++){
                bars[i] = sf::RectangleShape(sf::Vector2f(barWidth, 0.f));
                bars[i].setPosition(barWidth*i, SCREEN_HEIGHT);
                bars[i].setOrigin(0,1);
                bars[i].setFillColor(sf::Color::Cyan);
            }
        }

        void setHeights(vector<pair<int,int>> bands, vector<double> levels){
            for (int i = 0; i < bands.size(); i++){
                auto [start, end] = bands[i];
                double sum = 0;
                int count = 0;
                for (int j = start; j < end; j++){
                    sum += levels[j];
                    count++;
                }
                double avgPower = count > 0 ? sum / count : 0;

                double db = clampdB(powerTodB(avgPower));
                float norm = (db + DB_LOW)/DB_LOW;

                norm = pow(norm, 1.5);

                float height = norm * SCREEN_HEIGHT / 2;
                float curWidth = bars[i].getSize().x;

                heights[i] = alpha*height + (1-alpha)*heights[i]; // EMA

                bars[i].setSize(sf::Vector2f(curWidth, -heights[i]));
            }
        }

        void draw(){
            for (int i = 0; i < numBars; i++){
                window->draw(bars[i]);
            }
        }
};


int main(void){
    SF_INFO info;
    SNDFILE* file  = sf_open(SONG_DIR, SFM_READ, &info);
    printFileInfo(&info);
    const int samplerate = info.samplerate;
    const float duration = info.frames/samplerate;

    // Load entire song samples and close the file
    double* samples = static_cast<double*>(fftw_malloc((sizeof(double) * info.frames * info.channels)));
    double* leftChannel = static_cast<double*>(fftw_malloc((sizeof(double) * info.frames)));
    double* rightChannel = static_cast<double*>(fftw_malloc((sizeof(double) * info.frames)));
    int readcount = (int)sf_readf_double(file, samples, info.frames);
    for (int i = 0; i < readcount; i++){
        leftChannel[i] = samples[2*i];
        rightChannel[i] = samples[2*i+1];
    }
    fftw_free(samples);
    SongData song(leftChannel, rightChannel, info.frames, info.samplerate);
    sf_close(file);

    // Start PA
    PaError err = Pa_Initialize();
    if(err != paNoError) Pa_Terminate();
    // Create PA stream
    PaStream* stream;
    err = Pa_OpenDefaultStream(
        &stream,
        0,
        2,
        paFloat32,
        samplerate,
        FRAMES_PER_BUFFER,
        audioCallback,
        &song
    );
    if(err != paNoError) Pa_Terminate();
    // Start PA stream
    err = Pa_StartStream(stream);
    if(err != paNoError) Pa_Terminate();
    int64_t startTime = getTime();

    // Graphic design is my passion
    sf::RenderWindow window(sf::VideoMode(SCREEN_WIDTH, SCREEN_HEIGHT), "Vizualizer");
    vector<pair<int, int>> bands = createOctaveBands(info.samplerate);

    int nBars = bands.size();

    // for (const auto &[low, high] : bands){
    //     printf("%d, %d\n", low, high);
    // }



    BarVizualizer viz = BarVizualizer(&window, nBars);
    // circle
    double r = 300;
    sf::CircleShape shape(300.f);
    shape.setOrigin(-SCREEN_WIDTH/2+r, -SCREEN_HEIGHT/2+r);
    shape.setFillColor(sf::Color::Green);
    float scale = 0;
    const float alpha = 0.5;
    while (getTime() - startTime < (duration-0.5)*1000){
        sf::Event event;
        while (window.pollEvent(event)){
            if (event.type == sf::Event::Closed){
                window.close();
                goto close;
            }
        }
        viz.setHeights(bands, song.levels);
        window.clear();
        viz.draw();
        scale = (1-alpha)*scale + alpha*song.soundLevel;
        shape.setOrigin(-SCREEN_WIDTH/2+r*scale, -SCREEN_HEIGHT/2+r*scale);
        shape.setRadius(r*scale);
        window.draw(shape);
        window.display();
    }

close:
    fftw_free(rightChannel);
    fftw_free(leftChannel);
    // Stop and close stream
    err = Pa_StopStream(stream);
    if(err != paNoError) Pa_Terminate();
    err = Pa_CloseStream(stream);
    if(err != paNoError) Pa_Terminate();
    // Terminate
    Pa_Terminate();
    return err;
}