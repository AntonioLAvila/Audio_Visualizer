#pragma once

#define SCREEN_WIDTH 1920
#define SCREEN_HEIGHT 1080
#define FRAMES_PER_BUFFER 2048
#define FFT_LENGTH 4096
#define FFT_OUT_LENGTH (2048 / 2 + 1)
#define REQUESTED_NUMBER_OF_POINTS 64
#define DB_LOW 120.0

#ifndef SONG_DIR
#define SONG_DIR "music/snow.wav"
#endif