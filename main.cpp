#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>

#include "Spectrometer.h"

Spectrometer* spectrometer;

void sigintHandler(int s)
{
    spectrometer->Stop();
    return;
}

int main(int argc, char** argv)
{
    // Expect two command line parameter with display configuration filename.
    if (argc != 2)
    {
      throw std::runtime_error("USAGE: spectrometer [config path]");
    }
    
    spectrometer = new Spectrometer(argv[1]);
    
    signal(SIGINT, sigintHandler);
    
    sleep(1);
    
    spectrometer->Start();
    delete spectrometer;
    fprintf(stderr, "Exiting Program\n");
    return 0;
}
