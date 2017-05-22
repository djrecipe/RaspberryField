#ifndef SPECTROMETER_H
#define SPECTROMETER_H

#include <alsa/asoundlib.h>
#include <sndfile.h>

#include <led-matrix.h>
#include <signal.h>

#include <cstdint>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "Config.h"
#include "glcdfont.h"
#include "gpu_fft.h"
#include "GridTransformer.h"
#include "mailbox.h"

#define PCM_DEVICE "plughw:1,0"
#define FFT_JOBS 1
#define FFT_LOOPS 1
#define FFT_LOG 9
#define BIN_COUNT 16
#define BIN_DEPTH 8
#define SAMP_RATE 11025

class Spectrometer
{
    public:
    Spectrometer(char* config_path);
    ~Spectrometer();
    
    void Start();
    void Stop();
    private:
    
    int panelWidth = 0;
    int panelHeight = 0;
    int mailbox;
    struct GPU_FFT *fft;
    
    snd_pcm_t *pcm_handle;
    
    bool running;
    Config* config;    
    rgb_matrix::GPIO io;
    GridTransformer grid;
    rgb_matrix::RGBMatrix* canvas;
    bool** pixels;
	
	unsigned char* lib_logo;
    
    void GetBins(short* buffer, int* bins, bool logarithmic);
    void InitializeAudioDevice();
    void InitializeLEDMatrix(char* config_path);
    void InitializeFFT();
    void PrintBars(int bins[][BIN_COUNT], bool** pixels, bool print_black);
    void PrintBitmap(int bins[][BIN_COUNT], bool** pixels, unsigned char* data);
    void PrintBlack(bool** pixels);
	void PrintRadial(int bins[][BIN_COUNT], bool** pixels);
    void PrintText(int x, int y, const std::string& message, int r = 255, int g = 255, int b = 255);
	void ReadBitmap(char* filename, unsigned char* data);
      
};

#endif