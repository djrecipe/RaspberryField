// LED matrix library transformer to map a rectangular canvas onto a complex
// chain of matrices.
// Author: Tony DiCola
#ifndef GRIDTRANSFORMER_H
#define GRIDTRANSFORMER_H

#include <cassert>
#include <vector>
#include <math.h>

#include "led-matrix.h"


class GridTransformer: public rgb_matrix::Canvas, public rgb_matrix::CanvasTransformer {
public:
  struct Panel {
    int order;
    int rotate;
    int parallel;
  };
  GridTransformer();
  GridTransformer(int width, int height, int panel_width, int panel_height,
                  int chain_length, const std::vector<Panel>& panels);
  virtual ~GridTransformer();

  // Canvas interface implementation:
  virtual int width() const {
    return _width;
  }
  virtual int height() const {
    return _height;
  }
  virtual void Clear() {
    assert(_source != NULL);
    _source->Clear();
  }
  virtual void Fill(uint8_t red, uint8_t green, uint8_t blue) {
    assert(_source != NULL);
    _source->Fill(red, green, blue);
  }
  
	// GridTransformer methods
	void FillRemaining(uint8_t red, uint8_t green, uint8_t blue);
	bool GetPixelState(int x, int y);
	void SetCutoff(int value);
	void SetMaxBrightness(int value);
	void SetMirrorX(bool value);
	void SetMirrorY(bool value);
	
	// Canvas methods
	virtual void SetPixel(int x, int y, uint8_t red, uint8_t green, uint8_t blue);
	virtual rgb_matrix::Canvas* Transform(rgb_matrix::Canvas* source);

	// Other attribute accessors.
	int getRows() const
	{
		return _rows;
	}
	int getColumns() const
	{
		return _cols;
	}

private:
  int _width,
      _height,
      _panel_width,
      _panel_height,
      _chain_length,
      _rows,
      _cols;
	int cutoff;
	int maxBrightness;
	bool mirrorX, mirrorY, overrideCutoff;
	bool** pixels;
	rgb_matrix::Canvas* _source;
	std::vector<Panel> _panels;
	
	void ResetPixels();
};

#endif
