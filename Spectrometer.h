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
// # of frequency bins
#define BIN_COUNT 16
// history count for each frequency bin (for normalization)
#define TOTAL_BIN_DEPTH 64
// history count for each frequency bin (for display)
#define BIN_DEPTH 8
// number of "fan blades" per bin value
#define RADIAL_FAN_COUNT 8
// angle between each radial "fan blade" (angle = pi/value)
#define RADIAL_FAN_SPACING 32.0
// capture sample rate
#define SAMP_RATE 11025


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
    int mailbox;
    struct GPU_FFT *fft;
    
    snd_pcm_t *pcm_handle;
    
    bool running;
    Config* config;    
    rgb_matrix::GPIO io;
    GridTransformer grid;
    rgb_matrix::RGBMatrix* canvas;
	
	DisplayMode displayMode;
	
	unsigned char* lib_logo;
    
    void GetBins(short* buffer, int* bins, bool logarithmic);
	void GetExcludedPixels(bool** exclude);
	int GetRandomNumber(int min, int max);
    void InitializeAudioDevice();
    void InitializeLEDMatrix(char* config_path);
    void InitializeFFT();
    void NormalizeBins(int bins[][BIN_COUNT], int normalized_bins[][BIN_COUNT]);
    void PrintBars(int bins[][BIN_COUNT], bool** pixels);
    void PrintBitmap(int bins[][BIN_COUNT], bool** pixels, unsigned char* data);
    void PrintBlack(bool** pixels, bool** exclude);
	void PrintRadial(int bins[][BIN_COUNT], bool** pixels, float seconds);
    void PrintText(int x, int y, const std::string& message, int r = 255, int g = 255, int b = 255);
	void ReadBitmap(char* filename, unsigned char* data);
	//void RemoveExclusions(bool** exclude);
      
};

#endif