// Matrix configuration parsing class declaration.
// Author: Tony DiCola
#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>

#include "GridTransformer.h"

class Config {
public:
  Config(const std::string& filename);
  ~Config();

  // Attribute accessors:
  float getAnimationDuration(int set_index) const {
	return _animation_durations[set_index];
  }
  
  int getDisplayWidth() const {
    return _display_width;
  }
  int getDisplayHeight() const {
    return _display_height;
  }
  int getLEDCutoff() const {
	return _led_cutoff;
  }
  int getLEDMaxBrightness() const {
	return _led_max_brightness;
  }
  int getPanelWidth() const {
    return _panel_width;
  }
  int getPanelHeight() const {
    return _panel_height;
  }
  int getChainLength() const {
    return _chain_length;
  }
  const char* getImage(int set_index, int image_index)
  {
	  return (*_image_sets[set_index])[image_index].c_str();
  }
  int getImageCount(int index) const {
	return _image_sets[index]->size();
  }
  int getImageSetCount() const {
	return _image_sets.size();
  }
  int getImageSetDuration() const {
	  return _image_set_duration;
  } 
  int getParallelCount() const {
    return _parallel_count;
  }
  GridTransformer* getGridTransformer() const
  {
		return new GridTransformer(_display_width, _display_height, _panel_width, _panel_height, _chain_length, _panels);
  }
  bool hasCropOrigin() const {
    return (_crop_x > -1) && (_crop_y > -1);
  }
  int getCropX() const {
    return _crop_x;
  }
  int getCropY() const {
    return _crop_y;
  }

private:
  int _display_width,
      _display_height,
      _panel_width,
      _panel_height,
      _chain_length,
      _parallel_count,
      _crop_x,
      _crop_y;
  int _led_cutoff,
	  _led_max_brightness;
  int _image_set_duration;
  std::vector<float> _animation_durations;
  std::vector<GridTransformer::Panel> _panels;
  std::vector<std::vector<std::string>*> _image_sets;
};

#endif
