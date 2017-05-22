# Configure the rpi-rgb-led-matrix library here:
# The -DADAFRUIT_RGBMATRIX_HAT value configures the library to use the Adafruit
# LED matrix HAT wiring, and the -DRGB_SLOWDOWN_GPIO=1 value configures the
# library to work with a Raspberry Pi 2.  For a Pi 1 (or perhaps even on a Pi 2,
# but I found it necessary in my testing) you can remove the -DRGB_SLOWDOWN_GPIO=1
# option.  You can also add any other rpi-rgb-led-matrix library defines here
# to configure the library for more special needs.  See the library's docs for
# details on options:
#   https://github.com/hzeller/rpi-rgb-led-matrix/blob/master/lib/Makefile
export DEFINES = -DADAFRUIT_RGBMATRIX_HAT -DRGB_SLOWDOWN_GPIO=1

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O3 -I. -I./rpi-rgb-led-matrix/include -I/opt/vc/include -I/opt/vc/include/interface/vcos/pthreads -I/opt/vc/include/interface/vmcs_host -I/opt/vc/include/interface/vmcs_host/linux -L./rpi-rgb-led-matrix/lib -L/opt/vc/lib
LIBS = -lrgbmatrix -lrt -lm -lpthread -lbcm_host -lconfig++ -ldl -lasound -lsndfile

# Makefile rules:
all: hex spectrometer

spectrometer: main.o mailbox.o gpu_fft.o gpu_fft_base.o gpu_fft_shaders.o gpu_fft_twiddles.o Spectrometer.o GridTransformer.o Config.o glcdfont.o ./rpi-rgb-led-matrix/lib/librgbmatrix.a
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

./rpi-rgb-led-matrix/lib/librgbmatrix.a:
	$(MAKE) -C ./rpi-rgb-led-matrix/lib

.PHONY: clean

clean:
	rm -f *.o rpi-fb-matrix spectrometer
	$(MAKE) -C ./rpi-rgb-led-matrix/lib clean
