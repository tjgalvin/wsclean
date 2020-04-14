#include "fitsatermbase.h"

#include "../wsclean/logger.h"

#include "../units/imagecoordinates.h"

#include "../fftresampler.h"

#include <algorithm>
#include <limits>

FitsATermBase::FitsATermBase(size_t nAntenna, size_t width, size_t height, double ra, double dec, double dl, double dm, double phaseCentreDL, double phaseCentreDM, size_t atermSize) :
	_cache(nAntenna * 4 * width * height),
	_curTimeindex(0),
	_curFrequency(0),
	_nFrequencies(0),
	_nAntenna(nAntenna),
	_width(width),
	_height(height),
	_ra(ra), _dec(dec),
	_dl(dl), _dm(dm),
	_phaseCentreDL(phaseCentreDL),
	_phaseCentreDM(phaseCentreDM),
	_allocatedWidth(atermSize),
	_allocatedHeight(atermSize),
	_window(WindowFunction::Rectangular),
	_padding(1.0)
{
}

FitsATermBase::~FitsATermBase()
{ }

void FitsATermBase::initializeFromFiles(std::vector<FitsReader>& readers)
{
	// Sort the readers on observation time
	std::sort(readers.begin(), readers.end(), [](const FitsReader& a, const FitsReader& b)->bool
		{ return a.TimeDimensionStart() < b.TimeDimensionStart(); }
	);
	_nFrequencies = readers.front().NFrequencies();
	for(size_t readerIndex=0; readerIndex!=readers.size(); ++readerIndex)
	{
		const FitsReader& reader = readers[readerIndex];
		if(_nFrequencies != reader.NFrequencies())
			throw std::runtime_error("A-term FITS files have inconsistent number of frequencies");
		if(reader.NAntennas() == 1 && _nAntenna != 1)
		{
			Logger::Debug << "A-term fits file has length of one in antenna direction:\n"
				"Using single item for all " << _nAntenna << " antennas.\n";
		}
		else if(reader.NAntennas() != _nAntenna)
		{
			std::ostringstream str;
			str << "FITS file for A-terms has incorrect number of antennas. Measurement set has "
			<< _nAntenna << " antennas, a-term FITS file has " << reader.NAntennas() << " antennas.";
			throw std::runtime_error(str.str());
		}
		double time0 = reader.TimeDimensionStart();
		if(!_timesteps.empty() && time0 < _timesteps.back().time)
			throw std::runtime_error("Time axis of FITS files seem to overlap (start of fitsfile " + std::to_string(readerIndex) +
				" (t=" + std::to_string(time0) + " was before end of previous fitsfile)");
		for(size_t i=0; i!=reader.NTimesteps(); ++i)
			_timesteps.emplace_back(Timestep{time0 + i*reader.TimeDimensionIncr(), readerIndex, i});
	}
	_curTimeindex = std::numeric_limits<size_t>::max();
	_curFrequency = std::numeric_limits<double>::max();
}

double FitsATermBase::AverageUpdateTime() const
{
	if(_timesteps.size() < 2)
		return 60.0*30.0;
	else
		return (_timesteps.back().time - _timesteps.front().time) / (_timesteps.size() - 1);
}

bool FitsATermBase::findFilePosition(std::complex<float>* buffer, double time, double frequency, size_t& timeindex, bool& requiresRecalculation)
{
	requiresRecalculation = false;
	if(_curTimeindex == std::numeric_limits<size_t>::max())
	{
		requiresRecalculation = true;
		_cache.Reset();
		_curTimeindex = 0;
	}
	
	bool finishedSearch = false;
	while(_curTimeindex+1 < _timesteps.size() && !finishedSearch)
	{
		// Do we need to calculate a next timestep?
		double curTime = _timesteps[_curTimeindex].time;
		double nextTime = _timesteps[_curTimeindex+1].time;
		// If we are closer to the next timestep, use the next.
		if(std::fabs(nextTime - time) < std::fabs(curTime - time))
		{
			++_curTimeindex;
			requiresRecalculation = true;
			_cache.Reset();
			finishedSearch = false;
		}
		else {
			finishedSearch = true;
		}
	}
	timeindex = _curTimeindex;
	
	if(!requiresRecalculation)
	{
		// If we are here, it means that the timestep didn't
		// change. So if the frequency also didn't change, we're done...
		if(_curFrequency == frequency)
			return false;
		// If it did change: do we have this frequency in the cache?
		size_t cacheIndex = _cache.Find(frequency);
		if(cacheIndex == Cache::NOT_FOUND)
			requiresRecalculation = true;
		else {
			_cache.Get(cacheIndex, buffer);
			_curFrequency = frequency;
			return true;
		}
	}
	
	return true;
}

void FitsATermBase::storeInCache(double frequency, const std::complex<float>* buffer)
{
	_curFrequency = frequency;
	_cache.Store(frequency, buffer);
}

void FitsATermBase::readAndResample(FitsReader& reader, size_t fileIndex, ao::uvector<double>& scratch, ao::uvector<double>& output)
{
	if(_resampler == nullptr)
	{
		_resampler.reset(new FFTResampler(_allocatedWidth, _allocatedHeight, _width, _height, 1, false));
		if(_window == WindowFunction::Tukey)
			_resampler->SetTukeyWindow(double(_allocatedWidth) / _padding, false);
		else
			_resampler->SetWindowFunction(_window, true);
	}
	
	reader.ReadIndex(output.data(), fileIndex);
	
	// First, the image is regridded on a smaller image that fits in the kernel support allocated for the aterms
	regrid(reader, scratch.data(), output.data());
	
	// Now, the small image is enlarged so that it matches the kernel size
	_resampler->Resample(scratch.data(), output.data());
}

void FitsATermBase::regrid(const FitsReader& reader, double* dest, const double* source)
{
	size_t inWidth = reader.ImageWidth(), inHeight = reader.ImageHeight();
	double inPixelSizeX = reader.PixelSizeX(), inPixelSizeY = reader.PixelSizeY();
	double inPhaseCentreRA = reader.PhaseCentreRA(), inPhaseCentreDec = reader.PhaseCentreDec();
	double inPhaseCentreDL = reader.PhaseCentreDL(), inPhaseCentreDM = reader.PhaseCentreDM();
	
	// The full size is regridded onto the 'Nyquist-sampled' image to remove high-frequency
	// components. atermDL/DM are the pixelsizes of the Nyquist-sampled image.
	double atermDL = _dl * _width / _allocatedWidth;
	double atermDM = _dm * _height / _allocatedHeight;
	/**
	 * If phase centra of input and output are the same, i.e. they have the same
	 * tangential plane, a few calculations can be saved.
	 */
	bool samePlane = inPhaseCentreRA == _ra && inPhaseCentreDec == _dec;
	
	size_t index = 0;
	for(size_t y=0; y!=_allocatedHeight; ++y)
	{
		for(size_t x=0; x!=_allocatedWidth; ++x)
		{
			double l, m;
			ImageCoordinates::XYToLM(x, y, atermDL, atermDM, _allocatedWidth, _allocatedHeight, l, m);
			l += _phaseCentreDL;
			m += _phaseCentreDM;
			if(!samePlane)
			{
				double pixra, pixdec;
				ImageCoordinates::LMToRaDec(l, m, _ra, _dec, pixra, pixdec);
				ImageCoordinates::RaDecToLM(pixra, pixdec, inPhaseCentreRA, inPhaseCentreDec, l, m);
			}
			l -= inPhaseCentreDL;
			m -= inPhaseCentreDM;
			int inX, inY;
			ImageCoordinates::LMToXY(l, m, inPixelSizeX, inPixelSizeY, inWidth, inHeight, inX, inY);
			if(inX < 0 || inY < 0 || inX >= int(inWidth) || inY >= int(inHeight))
				dest[index] = 0;
			else {
				dest[index] = source[inX + inY * inWidth];
			}
			++index;
		}
	}
}

