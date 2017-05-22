#include "Spectrometer.h"

using namespace std;
using namespace rgb_matrix;

Spectrometer::Spectrometer(char* config_path)
{
    this->running = false;
    this->InitializeAudioDevice();
    this->InitializeFFT();
    this->InitializeLEDMatrix(config_path);   
    this->pixels = new bool*[this->panelWidth];
    for(int i=0; i<this->panelWidth; i++)
    {
        this->pixels[i] = new bool[this->panelHeight];
    }
    this->lib_logo = new unsigned char[this->panelWidth * this->panelHeight * 3];
    this->ReadBitmap("lib-logo.bmp", this->lib_logo);
    return;
}

void Spectrometer::GetBins(short* buffer, int* bins, bool logarithmic)
{
    int full_count = 1<<FFT_LOG;
    float frequencies[BIN_COUNT + 1]= {20.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,600.0,750.0,1000.0,2000.0,3000.0,5000.0,7500.0};
    float avgs[BIN_COUNT];
    int counts[BIN_COUNT];
    float ratio = 1.0;
    
    int i=0,j=0;
    
    for(i=0; i<BIN_COUNT; i++)
    {
        avgs[i] = 0.0;
        counts[i] = 0;
    }
    
    // input buffer
    for (i=0; i<full_count; i++)
    {
        fft->in[i].re = (float)buffer[i];
        fft->in[i].im =  0.0;
    }
    
    // execute
    gpu_fft_execute(this->fft); // call one or many times

    for(i=0; i<full_count/2; i++)
    {
        float frequency = (float)i * ((float)(SAMP_RATE)/(float)(full_count));
        for(j=0; j<BIN_COUNT; j++)
        {
            if(frequency >= frequencies[j] && frequency < frequencies[j+1])
            {
                avgs[j] += sqrt(pow(this->fft->out[i].re, 2) + pow(this->fft->out[i].re, 2));
                counts[j]++;
                break;
            }
        }
    }
    for(i=0; i<BIN_COUNT; i++)
    {
		if(counts[i] ==0)
			avgs[i] = 0.0;
		else
			avgs[i] = avgs[i] / (float)counts[i];
        // convert to db
        if(logarithmic)
            avgs[i] = 20.0 * log10(avgs[i]);
        // normalize
        ratio = pow((float)i/(float)BIN_COUNT + 0.5, 2.0);
        bins[i] = fmin(fmax(avgs[i] - 100.0 * ratio, 0.0) * 3.0, 100.0);
    }
    return;
}

void Spectrometer::InitializeAudioDevice()
{
    fprintf(stderr, "Initializing Audio Device\n");
    int err;
    /* Open the PCM device in playback mode */
    if((err = snd_pcm_open(&pcm_handle, PCM_DEVICE, SND_PCM_STREAM_CAPTURE , 0)) < 0)
	{
        fprintf (stderr, "cannot open audio device %s (%s)\n", PCM_DEVICE, snd_strerror (err));
        exit (1);
	}

    /* Allocate parameters object and fill it with default values*/
    snd_pcm_hw_params_t *params;
    if((err = snd_pcm_hw_params_malloc(&params)) < 0)
    {
        fprintf (stderr, "cannot allocate hardware parameter structure (%s)\n", snd_strerror (err));
        exit (1);
    }
    if((err = snd_pcm_hw_params_any(pcm_handle, params)) < 0)
    {
        fprintf (stderr, "cannot initialize hardware parameter structure (%s)\n",
        snd_strerror (err));
    }
    /* Set parameters */
    if((err = snd_pcm_hw_params_set_access(pcm_handle, params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0)
    {
        fprintf (stderr, "cannot set access type (%s)\n", snd_strerror (err));
        exit (1);
    }
    
    if((err = snd_pcm_hw_params_set_format(pcm_handle, params, SND_PCM_FORMAT_S16_LE)) < 0)
    {
        fprintf (stderr, "cannot set sample format (%s)\n", snd_strerror (err));
        exit (1);
    }
    if((err = snd_pcm_hw_params_set_channels(pcm_handle, params, 1)) < 0)
    {
        fprintf (stderr, "cannot set channel count (%s)\n", snd_strerror (err));
        exit (1);
    }
    
    if((err = snd_pcm_hw_params_set_rate(pcm_handle, params, SAMP_RATE, 0)) < 0)
    {
        fprintf (stderr, "cannot set sample rate (%s)\n", snd_strerror (err));
        exit (1);
    }
    if((err = snd_pcm_hw_params(pcm_handle, params)) < 0)
    {
        fprintf (stderr, "cannot set parameters (%s)\n", snd_strerror (err));
        exit (1);
    }
    
    snd_pcm_hw_params_free (params);

    if ((err = snd_pcm_prepare (pcm_handle)) < 0)
    {
        fprintf (stderr, "cannot prepare audio interface for use (%s)\n", snd_strerror (err));
        exit (1);
    }
    
    fprintf(stderr,"Successfully Initialized Audio Device\n");
    return;
}
/*
SNDFILE* InitializeAudioFile(char* path)
{
    SNDFILE *infile = NULL;
    fprintf(stderr, "\tFile: %s\n", path);
    infile = sf_open(path, SFM_READ, &sfinfo);
    fprintf(stderr,"\tChannels: %d\n", sfinfo.channels);
    fprintf(stderr,"\tSample rate: %d\n", sfinfo.samplerate);
    fprintf(stderr,"\tSections: %d\n", sfinfo.sections);
    fprintf(stderr,"\tFormat: %d\n", sfinfo.format);
    return infile;
}
*/
void Spectrometer::InitializeFFT()
{
    fprintf(stderr, "Initializing FFT\n");
    this->mailbox = mbox_open();
    int ret = gpu_fft_prepare(this->mailbox, FFT_LOG, GPU_FFT_REV, FFT_JOBS, &(this->fft));
      
    switch(ret)
    {
        case -1: throw std::runtime_error("Unable to enable V3D. Please check your firmware is up to date.\n");
        case -2: throw std::runtime_error("log2_N=%d not supported.  Try between 8 and 22.\n");
        case -3: throw std::runtime_error("Out of memory.  Try a smaller batch or increase GPU memory.\n");
        case -4: throw std::runtime_error("Unable to map Videocore peripherals into ARM memory space.\n");
        case -5: throw std::runtime_error("Can't open libbcm_host.\n");
    }
    int buffer_size = 1<<FFT_LOG;
    fprintf(stderr, "\tBuffer Size: %d\n", buffer_size);
    fprintf(stderr, "Successfully Initialized FFT\n");   
    return;
}
void Spectrometer::InitializeLEDMatrix(char* config_path)
{
    fprintf(stderr, "Initializing LED Matrix\n");
    fprintf(stderr, "\tFile: %s\n", config_path);
    
    // initialize parameters
    this->config = new Config(config_path);
    int width = this->config->getPanelWidth();
    int height = this->config->getPanelHeight();
    int chain_length = this->config->getChainLength();
    int parallel_count = this->config->getParallelCount(); 
    this->panelWidth = width * chain_length;
    this->panelHeight = height;
    fprintf(stderr, "\tSize: %d x %d\n", this->panelWidth, this->panelHeight);
    
    // Initialize GPIO
    if (!this->io.Init())
    {
      throw runtime_error("Error while initializing rpi-led-matrix library");
    }
    // Initialize Canvas
    this->canvas = new RGBMatrix(&io, height, chain_length, parallel_count); 
    this->grid = this->config->getGridTransformer();
    // set max pixel value (aka 'brightness')
    this->grid->SetMaxBrightness(225);
    // turn on horizontal mirroring (due to incorrect wiring)
    this->grid->SetMirrorX(true);
    // set grid transformer (for virtual calls)
    this->canvas->SetTransformer(&grid);
    // clear canvas
    this->canvas->Fill(0, 0, 0);
    
    fprintf(stderr, "Successfully Initialized LED Matrix\n");   
    return;
}
void Spectrometer::PrintBars(int bins[][BIN_COUNT], bool** pixels, bool print_black)
{
    // intialize parameters
	int col_width = (float)this->panelWidth/(float)BIN_COUNT;
    int i=0,j=0,r=0,g=0,b=0,x=0,x_index=0,y=0,bar_height=0;
    float ratio=0.0,gain=0.0,r_gain=0.0,g_gain=0.0,b_gain=0.0;
    for(i=0; i<BIN_DEPTH; i++)
    {
        // calculate base gain (based on age)
        gain = (float)(BIN_DEPTH - i - 1)/BIN_DEPTH;
        // iterate through bins
        for(j=0; j<BIN_COUNT; j++)
        {
            // calculate amplitude
            ratio = (float)bins[i][j] / CENTER_AMPLITUDE;
            bar_height = (int)(ratio * (float)this->panelHeight);
            // iterate across
            for(x=0; x<col_width; x++)
            {
                // calculate x index
                x_index = j * col_width + x;
                // iterate upwards
                for(y=0; y<=bar_height; y++)
                {
                    // check if pixel is already occupied
                    if(pixels[x_index][y])
                        continue;
                    // calculate colors
                    r_gain = 1.0*gain * ((float)j/(float)(BIN_COUNT/2)); // increases with frequency
                    g_gain = 1.0*gain * ratio;                     // increases with amplitude
                    b_gain = 1.0*gain * (1.0/ratio);               // decreases with amplitude
                    r = 110.0*r_gain;
                    g = 110.0*g_gain;
                    b = 110.0*b_gain;
                    // draw
                    this->canvas->SetPixel(x_index, y, r, g, b);
                    pixels[x_index][y] = true;
                }
            }
        }
    }
	return;
}
void Spectrometer::PrintBlack(bool** pixels)
{
    int x=0,y=0;
    // iterate across
    for(x=0; x<this->panelWidth; x++)
    {
        // iterate upwards
        for(y=0; y<this->panelHeight; y++)
        {
            // check if pixel is already occupied
            if(pixels[x][y])
                continue;
            // draw
            this->canvas->SetPixel(x, y, 0, 0, 0);
            pixels[x][y] = true;
        }
    }
    return;
}
void Spectrometer::PrintBitmap(int bins[][BIN_COUNT], bool** pixels, unsigned char * data)
{
    // initialize parameters
    int x = 0, y = 0, r = 0, g = 0, b = 0, i = 0, j=0, max_index = 0;
    float decay = 0.0, bin_gain = 1.0, red_gain = 1.0, green_gain = 1.0, blue_gain = 1.0;
    // draw
    // calculate gain (based on amplitude)
    for(x =0; x<this->panelWidth; x++)
    {
        // iterate upwards
        for(y=0; y<this->panelHeight; y++)
        {
            // check if pixel is already occupied
            if(pixels[x][y])
                continue;
            // iterate through bin depthwise
            for(i=0; i<BIN_DEPTH; i++)
            {
                // calculate decay (based on age)
                decay = (float)(BIN_DEPTH - i - 1)/BIN_DEPTH;
                // iterate through bins and calculate gains
                for(j=0; j<BIN_COUNT; j++)
                {
                    // calculate base gain (based on bin amplitude)
                    bin_gain = (float)bins[i][j] / CENTER_AMPLITUDE;
                    // calculate red gain (increases with bin frequency)
                    red_gain += bin_gain * (float)(j+1)/(float)BIN_COUNT;
                    // calculate green gain (increases towards center frequency)
                    green_gain += bin_gain * (float)(BIN_COUNT/2-abs(j-BIN_COUNT/2))/(float)(BIN_COUNT/2);
                    // calculate blue gain (decreases with bin frequency)
                    blue_gain += bin_gain * (float)(BIN_COUNT - j)/(float)BIN_COUNT; 
                }
                red_gain = red_gain/(float)BIN_COUNT + 0.25;
                green_gain = green_gain/(float)BIN_COUNT + 0.25;
                blue_gain = blue_gain/(float)BIN_COUNT + 0.25;
                // calculate data index
                int index = y * this->panelWidth * 3 + x*3;
                // calculate color
                r = (float)data[index]*red_gain*decay;
                g = (float)data[index+1]*green_gain*decay;
                b = (float)data[index+2]*blue_gain*decay;
                if(r > 50 || g > 50 || b > 50)
                {
                    this->canvas->SetPixel(x, y, r, g, b);
                    pixels[x][y] = true;
                }
            }  
        }
    }
    return;
}

void Spectrometer::PrintText(int x, int y, const string& message, int r, int g, int b)
{
  // Loop through all the characters and print them starting at the provided
  // coordinates.
  for (auto c: message) {
    // Loop through each column of the character.
    for (int i=0; i<5; ++i) {
      unsigned char col = glcdfont[c*5+i];
      x += 1;
      // Loop through each row of the column.
      for (int j=0; j<8; ++j)
      {
        // Put a pixel for each 1 in the column byte.
        if ((col >> j) & 0x01)
        {
          this->canvas->SetPixel(x, y+j, r, g, b);
        }
      }
    }
    // Add a column of padding between characters.
    x += 1;
  }
}
void Spectrometer::PrintRadial(int bins[][BIN_COUNT], bool** pixels)
{
	// intialize parameters
	int col_width = (float)this->panelWidth/(float)BIN_COUNT;
    int i=0,j=0,k=0,r=0,g=0,b=0,x=0,y=0;
    float freq_ratio=0.0,amplitude_ratio=0.0,gain=0.0,base_angle=0.0,angle=0.0;
	// iterate through bins depthwise
    for(i=0; i<BIN_DEPTH; i++)
    {
		// calculate gain (1 -> 0) 
		gain = (float)(BIN_DEPTH-i)/(float)BIN_DEPTH;
		// iterate through bins
        for(j=0; j<BIN_COUNT; j++)
        {
			// calculate ratios
			freq_ratio = (float)(j+1)/(float)BIN_COUNT;
            amplitude_ratio = (float)bins[i][j] / CENTER_AMPLITUDE;
			// calculate color(s) (based on amplitude)
			g = 150.0*amplitude_ratio*gain;
            base_angle = freq_ratio * M_PI * 2.0 + M_PI/2.0;
			// iterate through set of close angles
			for(k =0; k<4; k++)
			{
				// calculate angle (based on frequency)
				angle = base_angle + (M_PI*(k-2))/16.0;
				// calculate output coordinates
				x = 31 + 32.0 * sin(angle) * amplitude_ratio;
				y = 15 + 16.0 * cos(angle) * amplitude_ratio;
				// check if already written
				if(!pixels[x][y])
				{
					// calculate color(s) (based on frequency)
					r = (100.0+x*3)*gain;
					b = (100.0+y*3)*gain;
					// draw
					this->canvas->SetPixel(x, y, r, b, g);
					pixels[x][y] = true;
				}
			}
        }
    }
	return;
}
void Spectrometer::ReadBitmap(char* filename, unsigned char* data)
{
    fprintf(stderr, "Reading bitmap (%s)\n", filename);
    
    // open file
    FILE* f = fopen(filename, "rb");
    if(f == NULL)
        throw "Argument Exception";

    // read header
    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f);
    int width = *(int*)&info[18];
    int height = *(int*)&info[22];
    // assert correct dimensions
    assert(width == this->config->getPanelWidth() * this->config->getChainLength());
    assert(height == this->config->getPanelHeight());

    //
    int row_padded = (width*3 + 3) & (~3);
    unsigned char tmp;

    int x = 0, y = 0, index = 0;
    for(y = 0; y < height; y++)
    {
        fread(&data[y*row_padded], sizeof(unsigned char), row_padded, f);
        for(x = 0; x < width; x++)
        {
            // Convert (B, G, R) to (R, G, B)
            index = y*row_padded+(x*3);
            tmp = data[index];
            data[index] = data[index+2];
            data[index+2] = tmp;
        }
    }
    fclose(f);
    fprintf(stderr, "Successfully read bitmap\n");
    return;  
}
void Spectrometer::Start()
{    
    fprintf(stderr, "Initializing display loop\n");
    // initialize
    this->running = true;
    this->canvas->Fill(0, 0, 0);
    this->canvas->SetPixel(1, 1, 50, 50, 50);
    int bins[BIN_DEPTH][BIN_COUNT];
    int i=0, j=0;
    for(i=0; i<BIN_DEPTH; i++)
    {
        for(j=0; j<BIN_COUNT; j++)
        {
            bins[i][j] = 0;   
        }
    }
    // create buffer
    int err;
    int buffer_size = 1<<FFT_LOG;
    short buf[buffer_size];
    // loop
    fprintf(stderr, "Successfully initialized display loop\n");
    fprintf(stderr, "Starting display loop\n");
    while((err = snd_pcm_readi (pcm_handle, buf, buffer_size)) == buffer_size)
    {
        // intialize data
        for(i=0; i<this->panelWidth; i++)
        {
            for(j=0; j<this->panelHeight; j++)
            {
                pixels[i][j] = false;
            }
        }
        if(!this->running)
        {
            fprintf(stderr, "\nAborting\n");
            break;
        }
        // move bins
        for(i=BIN_DEPTH-1; i>0; i--)
        {
            for(j=0; j<BIN_COUNT; j++)
            {
                bins[i][j] = bins[i-1][j];  
            }
        }
        // retrieve new bins
        this->GetBins(buf, &bins[0][0], true);

        // get time
        struct timespec time;
        clock_gettime(CLOCK_REALTIME, &time);
    
        // display
        if(time.tv_sec%60<10)
        {
            this->PrintBitmap(bins, pixels, this->lib_logo);
        }
		else if(time.tv_sec%60<30)
		{
            this->PrintRadial(bins, pixels);
		}
        else
        {  
            this->PrintBars(bins, pixels, true);
        }
        this->PrintBlack(pixels);
        // sleep
        usleep(10000);
    }
    fprintf(stderr, "Exiting display loop\n");
    return;
}

void Spectrometer::Stop()
{
    this->running = false;
    return;
}

Spectrometer::~Spectrometer()
{
    delete this->lib_logo;
    delete this->canvas;
    delete this->config;
    this->Stop();
    this->canvas->Clear();
    snd_pcm_drain(this->pcm_handle);
    snd_pcm_close(this->pcm_handle);
    gpu_fft_release(this->fft);
    for(int i=0; i<this->panelWidth; i++)
    {
        delete &(this->pixels[i]);
    }
    delete this->pixels;
    return;
}