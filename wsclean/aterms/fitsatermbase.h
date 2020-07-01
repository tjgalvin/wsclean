#ifndef FITS_ATERM_BASE_H
#define FITS_ATERM_BASE_H

#include "atermbase.h"
#include "cache.h"

#include "../fitsreader.h"
#include "../windowfunction.h"

#include <aocommon/uvector.h>

class FitsATermBase : public ATermBase
{
public:
	FitsATermBase(size_t nAntenna, const CoordinateSystem& coordinateSystem);
	~FitsATermBase();
	
	virtual double AverageUpdateTime() const override;
	
	void SetTukeyWindow(double padding)
	{
		_window = WindowFunction::Tukey;
		_padding = padding;
	}
	
	void SetWindow(WindowFunction::Type window)
	{
		_window = window;
	}
	
	void SetDownSample(bool downsample)
	{
		_downsample = downsample;
	}
	
protected:
	void initializeFromFiles(std::vector<FitsReader>& readers);
	
	bool findFilePosition(std::complex<float>* buffer, double time, double frequency, size_t& timeindex, bool& requiresRecalculation);
	
	void storeInCache(double frequency, const std::complex<float>* buffer);
	
	void readAndResample(FitsReader& reader, size_t fileIndex, aocommon::UVector<double>& scratch, aocommon::UVector<double>& output);
	
	struct Timestep {
		double time;
		size_t readerIndex;
		size_t imgIndex;
	};
	std::vector<Timestep> _timesteps;
	
	size_t Width() const { return _width; }
	size_t Height() const { return _height; }
	size_t AllocatedWidth() const { return _allocatedWidth; }
	size_t AllocatedHeight() const { return _allocatedHeight; }
	size_t NAntenna() const { return _nAntenna; }
	size_t NFrequencies() const { return _nFrequencies; }
	double DL() const { return _dl; }
	double DM() const { return _dm; }
	double PhaseCentreDL() const { return _phaseCentreDL; }
	double PhaseCentreDM() const { return _phaseCentreDM; }
	
private:
	void regrid(const FitsReader& reader, double* dest, const double* source);
	
	Cache _cache;
	size_t _curTimeindex;
	double _curFrequency;
	size_t _nFrequencies, _nAntenna, _width, _height;
	double _ra, _dec, _dl, _dm, _phaseCentreDL, _phaseCentreDM;
	size_t _allocatedWidth, _allocatedHeight;
	bool _downsample;
	WindowFunction::Type _window;
	double _padding;
	std::unique_ptr<class FFTResampler> _resampler;
};

#endif
