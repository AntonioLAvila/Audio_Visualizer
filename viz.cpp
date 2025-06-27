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

#define SCREEN_WIDTH 1920
#define SCREEN_HEIGHT 1080
#define FRAMES_PER_BUFFER 2048
#define FFT_LENGTH 2048
#define FFT_OUT_LENGTH 2048/2 + 1
#define REQUESTED_NUMBER_OF_POINTS 128

double* fftWindow = (double*)malloc(FFT_LENGTH * sizeof(double));

void fillHannWindow(){
    for (int n = 0; n < FFT_LENGTH; n++){
        fftWindow[n] = 0.5*(1-cos(2*M_PI*n/(FFT_LENGTH - 1)));
    }
}


void printFileInfo(SF_INFO* info){
    printf("Sample Rate = %d Hz\n", info->samplerate);
    printf("Channels = %d\n", info->channels);
    printf("Format = 0x%x\n", info->format);
    printf("Sections = %d\n", info->sections);
    printf("Seekable = %d\n", info->seekable);
    printf("Frames = %ld\n", info->frames);
}


double powerTodB(double A){
    return 10 * log10(A + 1e-10);
}


double clampdB(double db) {
    return std::min(0.0, std::max(-120.0, db));
}


double magnitude(fftw_complex c){
    return sqrt(pow(c[0],2) + pow(c[1],2));
}


void scale(fftw_complex c, double scale){
    c[0] = c[0]*scale;
    c[1] = c[1]*scale;
}


int64_t getTime(){
    auto time = std::chrono::high_resolution_clock::now();
    auto dur = time.time_since_epoch();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    return ms;
}


std::vector<int> linspace(double start, double stop, int num) {
    std::vector<int> res(num);
    double step = (stop - start)/num;
    for (int i = 0; i < num; i++){
        res[i] = start + (i * step);
    }
    return res;
}


std::vector<std::pair<int, int>> getBands(std::vector<int>& indices) {
    std::vector<std::pair<int, int>> bands;
    for (size_t i = 0; i < indices.size() - 1; i++) {
        bands.push_back({indices[i], indices[i + 1]});
    }
    return bands;
}


double aWeighting(double freq){
    double f2 = freq * freq;
    double num = 12200 * 12200 * f2 * f2;
    double den = (f2 + 20.6 * 20.6) * sqrt((f2 + 107.7 * 107.7) * (f2 + 737.9 * 737.9)) * (f2 + 12200 * 12200);
    return 20 * log10(num / den) + 2.0; // +2 dB calibration
}


typedef struct SongData {
    double* leftChannel;
    double* rightChannel;
    long int totalFrames;
    long int lastFrame;
    int samplerate;

    double soundLevel;
    std::vector<double> levels;

    fftw_complex* leftOutBuffer;
    fftw_complex* rightOutBuffer;

    SongData(double* l, double* r, long int t, int fs)
        : leftChannel(l), rightChannel(r), totalFrames(t), lastFrame(0), soundLevel(0.0), samplerate(fs), levels(FFT_LENGTH)
    {
        leftOutBuffer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_OUT_LENGTH);
        rightOutBuffer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFT_OUT_LENGTH);
    }

    ~SongData() {
        fftw_free(leftOutBuffer);
        fftw_free(rightOutBuffer);
    }
} SongData;


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
    double sum = 0;
    
    long int fftIdx;
    if (song->lastFrame + FFT_LENGTH >= song->totalFrames) {
        fftIdx = song->totalFrames - FFT_LENGTH;
    } else {
        fftIdx = song->lastFrame;
    }
    // copy and window
    double leftSignal[FFT_LENGTH];
    double rightSignal[FFT_LENGTH];
    for (int n = 0; n < FFT_LENGTH; n++){
        leftSignal[n] = song->leftChannel[fftIdx+n] * fftWindow[n];
        rightSignal[n] = song->rightChannel[fftIdx+n] * fftWindow[n];
    }
    // fft
    fftw_plan leftPlan = fftw_plan_dft_r2c_1d(
        FFT_LENGTH,
        leftSignal,
        song->leftOutBuffer,
        FFTW_ESTIMATE
    );
    fftw_plan rightPlan = fftw_plan_dft_r2c_1d(
        FFT_LENGTH,
        rightSignal,
        song->rightOutBuffer,
        FFTW_ESTIMATE
    );
    fftw_execute(leftPlan); fftw_execute(rightPlan);

    for (int i = 0; i < FFT_LENGTH; i++){
        fftw_complex leftW, rightW;
        if (i < FFT_OUT_LENGTH){
            leftW[0] = song->leftOutBuffer[i][0];
            leftW[1] = song->leftOutBuffer[i][1];

            rightW[0] = song->rightOutBuffer[i][0];
            rightW[1] = song->rightOutBuffer[i][1];
        }else {
            int mirror_i = FFT_LENGTH - i;

            leftW[0] = song->leftOutBuffer[mirror_i][0];
            leftW[1] = -song->leftOutBuffer[mirror_i][1];

            rightW[0] = song->rightOutBuffer[mirror_i][0];
            rightW[1] = -song->rightOutBuffer[mirror_i][1];
        }
        double magLeft = magnitude(leftW);
        double magRight = magnitude(rightW);

        // Normalize
        double scale = 1.0 / FFT_LENGTH;
        if (i == 0 || i == FFT_LENGTH / 2){
            magLeft *= scale;
            magRight *= scale;
        }else{
            magLeft *= 2.0 * scale;
            magRight *= 2.0 * scale;
        }

        double power = (magLeft*magLeft + magRight*magRight)/2;
        double freq = i * song->samplerate / FFT_LENGTH;
        double db = powerTodB(sqrt(power)) + aWeighting(freq);
        song->levels[i] = clampdB(db);
    }

    for (int i = 0; i < framesPerBuffer; i++){
        *out++ = song->leftChannel[song->lastFrame];
        *out++ = song->rightChannel[song->lastFrame];
        double powerLeft = song->leftChannel[song->lastFrame] * song->leftChannel[song->lastFrame];
        double powerRight = song->rightChannel[song->lastFrame] * song->rightChannel[song->lastFrame];
        sum += (powerLeft + powerRight)/2;
        song->lastFrame += 1;
    }
    song->soundLevel = sum/framesPerBuffer;
    fftw_destroy_plan(leftPlan); fftw_destroy_plan(rightPlan);
    return 0;
}


class BarVizualizer{
    private:
        std::vector<sf::RectangleShape> bars;
        std::vector<double> heights;
        double barWidth;
        double numBars;
        sf::RenderWindow* window;
    
    public:
        BarVizualizer(sf::RenderWindow* win, int n, double width){
            barWidth = width;
            numBars = n;
            window = win;
            bars = std::vector<sf::RectangleShape>(numBars);
            heights = std::vector<double>(numBars, 0.f);
            for (int i = 0; i < numBars; i++){
                bars[i] = sf::RectangleShape(sf::Vector2f(barWidth, 0.f));
                bars[i].setPosition(barWidth*i, SCREEN_HEIGHT);
                bars[i].setOrigin(0,1);
                bars[i].setFillColor(sf::Color::Cyan);
            }
        }
        
        void setHeights(std::vector<int> indices, std::vector<double> levels){
            for (int i = 0; i < indices.size(); i++){
                float norm = (levels[indices[i]] + 120)/120; // [0, 1]
                norm = pow(norm, 1.5); // gamma compression
                float height = norm * SCREEN_HEIGHT/2;
                float curWidth = bars[i].getSize().x;
                bars[i].setSize(sf::Vector2f(curWidth, -height));
            }
        }

        void setHeightsFromBands(std::vector<std::pair<int,int>> bands, std::vector<double> levels){
            for (int i = 0; i < bands.size(); i++){
                auto [start, end] = bands[i];
                double sum = 0.0;
                int count = 0;
                for (int j = start; j < end; j++){
                    sum += levels[j];
                    count++;
                }
                double avg = count > 0 ? sum / count : -120;
                float norm = (avg + 120)/120;
                norm = pow(norm, 1.5); // gamma correction
                float height = norm * SCREEN_HEIGHT / 2;
                float curWidth = bars[i].getSize().x;
                bars[i].setSize(sf::Vector2f(curWidth, -height));
            }
        }

        void draw(){
            for (int i = 0; i < numBars; i++){
                window->draw(bars[i]);
            }
        }
};


int main(void){
    fillHannWindow();
    SF_INFO info;
    SNDFILE* file  = sf_open("./music/snow.wav", SFM_READ, &info);
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
    std::vector<int> logLinearIndices = linspace(0, song.levels.size(), REQUESTED_NUMBER_OF_POINTS);

    for (int i = 0; i < logLinearIndices.size(); i++){
        printf("%d, %d\n", logLinearIndices[i], logLinearIndices[i]*info.samplerate/FFT_LENGTH);
    }
    printf("FFT Resolution: %dHz\n", info.samplerate/FFT_LENGTH);

    int numBars = logLinearIndices.size();
    double barWidth = SCREEN_WIDTH/numBars;
    BarVizualizer viz = BarVizualizer(&window, numBars, barWidth);
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
        viz.setHeights(logLinearIndices, song.levels);
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
    free(fftWindow);
    // Stop and close stream
    err = Pa_StopStream(stream);
    if(err != paNoError) Pa_Terminate();
    err = Pa_CloseStream(stream);
    if(err != paNoError) Pa_Terminate();
    // Terminate
    Pa_Terminate();
    return err;
}