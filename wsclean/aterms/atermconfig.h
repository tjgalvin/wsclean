#ifndef ATERM_CONFIG_H
#define ATERM_CONFIG_H

#include <string>

#include "atermbase.h"
#include "atermbeam.h"
#include "atermstub.h"
#include "dishaterm.h"
#include "dldmaterm.h"
#include "fitsaterm.h"
#include "lofarbeamterm.h"
#include "mwabeamterm.h"
#include "telescope.h"

#include "../matrix2x2.h"
#include "../parsetreader.h"

#include "../wsclean/wscleansettings.h"

class ATermConfig : public ATermBase
{
public:
	ATermConfig(casacore::MeasurementSet& ms, size_t nAntenna, const CoordinateSystem& coordinateSystem,
		const WSCleanSettings& settings
	) :
		_ms(ms),
		_nAntenna(nAntenna),
		_coordinateSystem(coordinateSystem),
		_settings(settings)
	{ }
	
	void Read(const std::string& parset)
	{
		ParsetReader reader(parset);
		std::vector<std::string> aterms = reader.GetStringList("aterms");
		if(aterms.empty())
			throw std::runtime_error("No a-term correction given in parset (aterms key is an empty list)");
		
		for(const std::string atermName : aterms)
		{
			std::string atermType = reader.GetStringOr(atermName + ".type", atermName);
			if(atermType == "tec")
			{
				std::vector<std::string> tecFiles = reader.GetStringList(atermName + ".images");
				std::unique_ptr<FitsATerm> f(new FitsATerm(_nAntenna, _coordinateSystem));
				f->OpenTECFiles(tecFiles);
				std::string windowStr = reader.GetStringOr(atermName + ".window", "raised-hann");
				WindowFunction::Type window = WindowFunction::GetType(windowStr);
				if(window == WindowFunction::Tukey)
					f->SetTukeyWindow(double(_settings.paddedImageWidth) / _settings.trimmedImageWidth);
				else
					f->SetWindow(window);
				f->SetDownSample( reader.GetBoolOr(atermName + ".downsample", true) );
				_aterms.emplace_back(std::move(f));
			}
			else if(atermType == "diagonal")
			{
				std::vector<std::string> diagFiles = reader.GetStringList(atermName + ".images");
				std::unique_ptr<FitsATerm> f(new FitsATerm(_nAntenna, _coordinateSystem));
				f->OpenDiagGainFiles(diagFiles);
				std::string windowStr = reader.GetStringOr(atermName + ".window", "raised-hann");
				WindowFunction::Type window = WindowFunction::GetType(windowStr);
				if(window == WindowFunction::Tukey)
					f->SetTukeyWindow(double(_settings.paddedImageWidth) / _settings.trimmedImageWidth);
				else
					f->SetWindow(window);
				f->SetDownSample( reader.GetBoolOr(atermName + ".downsample", true) );
				_aterms.emplace_back(std::move(f));
			}
			else if(atermType == "dldm")
			{
				std::vector<std::string> dldmFiles = reader.GetStringList(atermName + ".images");
				std::unique_ptr<DLDMATerm> f(new DLDMATerm(_nAntenna, _coordinateSystem));
				f->Open(dldmFiles);
				f->SetUpdateInterval( reader.GetDoubleOr("dldm.update_interval", 5.0*60.0) );
				std::string windowStr = reader.GetStringOr(atermName + ".window", "raised-hann");
				WindowFunction::Type window = WindowFunction::GetType(windowStr);
				if(window == WindowFunction::Tukey)
					f->SetTukeyWindow(double(_settings.paddedImageWidth) / _settings.trimmedImageWidth);
				else
					f->SetWindow(window);
				f->SetDownSample( reader.GetBoolOr(atermName + ".downsample", true) );
				_aterms.emplace_back(std::move(f));
			}
			else if(atermType == "beam")
			{
				std::unique_ptr<ATermBeam> beam;
				switch(Telescope::GetType(_ms))
				{
					case Telescope::AARTFAAC:
					case Telescope::LOFAR: {
						bool differential = reader.GetBoolOr("beam.differential", false);
						bool useChannelFrequency = reader.GetBoolOr("beam.usechannelfreq", true);
						std::unique_ptr<LofarBeamTerm> lofarBeam(new LofarBeamTerm(_ms, _coordinateSystem, _settings.dataColumnName));
						lofarBeam->SetUseDifferentialBeam(differential);
						lofarBeam->SetUseChannelFrequency(useChannelFrequency);
						beam = std::move(lofarBeam);
						break;
					}
					case Telescope::MWA: {
						std::unique_ptr<MWABeamTerm> mwaTerm(new MWABeamTerm(_ms, _coordinateSystem));
						mwaTerm->SetSearchPath(_settings.mwaPath);
						beam = std::move(mwaTerm);
						break;
					}
					case Telescope::VLA: {
						beam.reset(new DishATerm(_ms, _coordinateSystem));
						break;
					}
					default: {
						// This is here to make sure ATermStub compiles. This call should be the
						// same as the call for LofarBeamTerm(..)
						beam.reset(new ATermStub(_ms, _coordinateSystem, _settings.dataColumnName));
						throw std::runtime_error("Can't make beam for this telescope");
					}
				}
				double updateInterval = reader.GetDoubleOr("beam.update_interval", _settings.beamAtermUpdateTime);
				beam->SetUpdateInterval(updateInterval);	
				_aterms.emplace_back(std::move(beam));
			}
			_aterms.back()->SetSaveATerms(false, _settings.prefixName);  // done by config after combining
		}
		Logger::Debug << "Constructed an a-term configuration with " << _aterms.size() << " terms.\n";
		if(_aterms.empty())
		{
			throw std::runtime_error("The specified a-term configuration does not define any terms to apply");
		}
		if(_aterms.size() > 1)
		{
			_previousAterms.resize(_aterms.size());
			for(ao::uvector<std::complex<float>>& buf : _previousAterms)
				buf.resize(_coordinateSystem.width * _coordinateSystem.height * _nAntenna * 4);
		}
	}
	
	virtual bool Calculate(std::complex<float>* buffer, double time, double frequency, size_t fieldId, const double* uvwInM) final override
	{
		if(_aterms.size() == 1)
		{
			bool result = _aterms.front()->Calculate(buffer, time, frequency, fieldId, uvwInM);
			if(result)
				saveATermsIfNecessary(buffer, _nAntenna, _coordinateSystem.width, _coordinateSystem.height);
			return result;
		}
		else {
			bool isUpdated = false;
			for(size_t i=0; i!=_aterms.size(); ++i)
			{
				bool atermUpdated = _aterms[i]->Calculate(_previousAterms[i].data(), time, frequency, fieldId, uvwInM);
				isUpdated = isUpdated || atermUpdated;
			}
			
			if(isUpdated)
			{
				std::copy(_previousAterms[0].begin(), _previousAterms[0].end(), buffer);
				for(size_t i=1; i!=_aterms.size(); ++i)
				{
					for(size_t j=0; j!=_coordinateSystem.width*_coordinateSystem.height*_nAntenna*4; j+=4)
					{
						std::complex<float> scratch[4];
						Matrix2x2::ATimesB(scratch, &_previousAterms[i][j], &buffer[j]);
						Matrix2x2::Assign(&buffer[j], scratch);
					}
				}
				saveATermsIfNecessary(buffer, _nAntenna, _coordinateSystem.width, _coordinateSystem.height);
			}
			
			return isUpdated;
		}
	}
	
	virtual double AverageUpdateTime() const override
	{
		double avgTime = _aterms.front()->AverageUpdateTime();
		for(size_t i=1; i<_aterms.size(); ++i)
			avgTime = std::min(avgTime, _aterms[i]->AverageUpdateTime());
		return avgTime;
	}
private:
	casacore::MeasurementSet& _ms;
	size_t _nAntenna;
	CoordinateSystem _coordinateSystem;
	std::vector<std::unique_ptr<ATermBase>> _aterms;
	std::vector<ao::uvector<std::complex<float>>> _previousAterms;
	const WSCleanSettings& _settings;
};

#endif
