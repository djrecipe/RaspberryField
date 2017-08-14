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
#include "RandomSequenceOfUnique.h"

// capture device
#define PCM_DEVICE "plughw:1,0"
// number of fft jobs (always 1?)
#define FFT_JOBS 1
// fft iterations
#define FFT_LOOPS 1
// determines fft length & buffer size
#define FFT_LOG 9
#define FULL_SCALE 100.0
// # of frequency bins
#define BIN_COUNT 16
// history count for each frequency bin (for normalization)
#define TOTAL_BIN_DEPTH 64
// history count for each frequency bin (for display)
#define BIN_DEPTH 8
// number of "fan blades" per bin value
#define RADIAL_FAN_COUNT 4
// angle between each radial "fan blade" (angle = pi/value)
#define RADIAL_FAN_SPACING 32.0
// capture sample rate
#define SAMP_RATE 11025
// sigmoid numerator value
// with logarithmic also enabled, decreasing this value allows you to stretch the sigmoid shape along the x-axis
#define SIGMOID_NUMERATOR 10.0
// sigmoid X offset (higher = more low frequency attenuation, more amplitude required to hit gain threshold)
// with logarithmic also enabled, increasing this number will exponentially increase the amount of amplitude required to hit full scale db
#define SIGMOID_OFFSET 25.0
// sigmoid sloep (higher = slopes slower, less effect of attenuation/gain)
// with logarithmic also enabled, increasing this number will result in a sharper corner and a closer resemblance to the 20log10 function
// *** this parameter is of great interest
#define SIGMOID_SLOPE 12.0

enum FFTOptions { None = 0, Logarithmic = 1, Sigmoid = 2, Autoscale = 4};

inline FFTOptions operator|(FFTOptions a, FFTOptions b){return static_cast<FFTOptions>(static_cast<int>(a)|static_cast<int>(b));}

class Spectrometer
{
	enum DisplayMode { Bars=0, Bitmap=1, Radial=2 };
	
	
    public:
    Spectrometer(char* config_path);
    ~Spectrometer();
    
    void Start();
    void Stop();
	
	
    private:
    
    int panelWidth = 0;
    int panelHeight = 0;
	float seconds = 0.0;
	float animationDuration = 1.0;
    int mailbox;
    struct GPU_FFT *fft;
    
    snd_pcm_t *pcm_handle;
    
    bool running;
    Config* config;    
    rgb_matrix::GPIO io;
    GridTransformer* grid;
    rgb_matrix::RGBMatrix* canvas;
	
	DisplayMode displayMode;
	
	std::vector<unsigned char*> logos;
    
    void GetBins(short* buffer, int* bins);
    unsigned int GetBitmapIndex(float seconds);
	int GetRandomNumber(int min, int max);
    void InitializeAudioDevice();
    void InitializeLEDMatrix(char* config_path);
    void InitializeFFT();
    void NormalizeBins(int bins[][BIN_COUNT], int normalized_bins[][BIN_COUNT], FFTOptions options);
    void PrintBars(int bins[][BIN_COUNT]);
    void PrintBitmap(int bins[][BIN_COUNT], unsigned char* data);
	void PrintRadial(int bins[][BIN_COUNT], float seconds);
    void PrintText(int x, int y, const std::string& message, int r = 255, int g = 255, int b = 255);
	void ReadBitmap(const char* filename, unsigned char* data);
    double SigmoidFunction(double value);
      
};

#endif