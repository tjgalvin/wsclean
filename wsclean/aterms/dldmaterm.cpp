#include "dldmaterm.h"

#include "../units/imagecoordinates.h"

#include "../wsclean/logger.h"

#include "../banddata.h"

DLDMATerm::DLDMATerm(size_t nAntenna, size_t width, size_t height, double ra, double dec, double dl, double dm, double phaseCentreDL, double phaseCentreDM, size_t atermSize) :
	FitsATermBase(nAntenna, width, height, ra, dec, dl, dm, phaseCentreDL, phaseCentreDM, atermSize),
	_updateInterval(60),
	_previousTime(0)
{
}

void DLDMATerm::Open(const std::vector<std::string>& filenames)
{
	_readers.reserve(filenames.size());
	for(const std::string& filename : filenames)
	{
		_readers.emplace_back(filename, true, true);
		if(_readers.back().NMatrixElements() != 2)
			throw std::runtime_error("FITS file for dl,dm offsets did not have 2 matrix elements in it");
	}
	initializeFromFiles(_readers);
}

bool DLDMATerm::Calculate(std::complex<float>* buffer, double time, double frequency, const double* uvwInM)
{
	size_t timeIndex;
	bool requiresRecalculation;
	bool positionChanged = findFilePosition(buffer, time, frequency, timeIndex, requiresRecalculation);
	bool outdated = std::fabs(time - _previousTime) > _updateInterval;
	if(!positionChanged && !outdated)
		return false;
	else {
		if(requiresRecalculation || outdated)
		{
			_previousTime = time;
			readImages(buffer, timeIndex, frequency, uvwInM);
			storeInCache(frequency, buffer);
		}
		return true;
	}
}

void DLDMATerm::readImages(std::complex<float>* buffer, size_t timeIndex, double frequency, const double* uvwInM)
{
	Logger::Debug << "DLDMATerm::readImages(buffer, " << timeIndex << ", " << frequency << ", uvw={... ";
	if(NAntenna()>1)
	{
		Logger::Debug << uvwInM[3] << ',' << uvwInM[4] << ',' << uvwInM[5];
	}
	Logger::Debug << " ... })\n";
	const size_t freqIndex = round((frequency - _readers.front().FrequencyDimensionStart())
		/ _readers.front().FrequencyDimensionIncr());
	const size_t imgIndex = _timesteps[timeIndex].imgIndex * NFrequencies() + freqIndex;
	FitsReader& reader = _readers[_timesteps[timeIndex].readerIndex];
	_scratch.resize(AllocatedWidth() * AllocatedHeight());
	_dlImage.resize(std::max(Width()*Height(), reader.ImageWidth() * reader.ImageHeight()));
	_dmImage.resize(std::max(Width()*Height(), reader.ImageWidth() * reader.ImageHeight()));
	readAndResample(reader, imgIndex * 2, _scratch, _dlImage);
	readAndResample(reader, imgIndex * 2 + 1, _scratch, _dmImage);

	double wavel = BandData::FrequencyToLambda(frequency);
	for(size_t antennaIndex = 0; antennaIndex != NAntenna(); ++antennaIndex)
	{
		double uvw[3] = {
			uvwInM[antennaIndex*3] / wavel,
			uvwInM[antennaIndex*3+1] / wavel,
			uvwInM[antennaIndex*3+2] / wavel
		};
		
		std::complex<float>* antennaBuffer = buffer + antennaIndex * Width()*Height()*4;
		evaluateDLDM(antennaBuffer, _dlImage.data(), _dmImage.data(), uvw);
	}
}

void DLDMATerm::evaluateDLDM(std::complex<float>* dest, const double* dl, const double* dm, const double* uvwInL)
{
	// For a single source at l,m, we have:
	//   dl = oldl - newl
	//      oldl: measured l
	//      newl: "real" (model) l
	//   V(u,v,w) = I(l, m) exp 2pi i ( lu + mv + nw ), with dn = sqrt(1-l^2-m^2) - sqrt(1-(l+dl)^2-(m+dm)^2) ~ 0;
	// dl,dm are the shifts in l,m. Given a0 as reference antenna, with u0, v0, w0 as coords,
	//   dphase = phase[ I(l, m) exp -2pi i ( dl(u-u0) + dm(v-v0) + dn(w-w0) ) ]
	//          = -2pi i ( l(u-u0) + m(v-v0) + dn(w-w0) )
	// The baselines uvw are already referenced to the first antenna (i.e. uvw=0,0,0 for antenna 0), so
	// uvwInL[0] is (u-u0).
	const double
		u = uvwInL[0],
		v = uvwInL[1],
		w = uvwInL[2];
	for(size_t y=0; y!=Height(); ++y)
	{
		for(size_t x=0; x!=Width(); ++x)
		{
			double l, m;
			ImageCoordinates::XYToLM(x, y, DL(), DM(), Width(), Height(), l, m);
			l += PhaseCentreDL();
			m += PhaseCentreDM();
			double lproj = l+(*dl), mproj = m+(*dm);
			double lmSq = l*l + m*m, lmprojSq = lproj*lproj + mproj*mproj;
			double dn;
			if(lmSq >= 1.0 || lmprojSq >= 1.0)
				dn = 0.0;
			else
				dn = std::sqrt(1.0 - lmprojSq) - std::sqrt(1.0 - lmSq);
			dest[0] = std::polar(1.0, 2.0*M_PI*(u*(*dl) + v*(*dm) + w*dn));
			dest[1] = 0.0;
			dest[2] = 0.0;
			dest[3] = dest[0];
			
			++dl; ++dm; dest += 4;
		}
	}
}
