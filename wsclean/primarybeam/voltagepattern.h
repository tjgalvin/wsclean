#ifndef VOLTAGE_PATTERN_H
#define VOLTAGE_PATTERN_H

#include "../uvector.h"

#include <complex>

// Holds the information for a symmetric voltage pattern
struct VoltagePattern
{
	// These are the radial (one-dimensional) values of the beam
	// It is array of size nsamples x nfrequencies, where the sample index is least significant (fastest changing)
	ao::uvector<double> values;
	
	// Array of since nfrequencies
	ao::uvector<double> frequencies;
	
	double inverseIncrementRadius;
	double maximumRadiusArcMin;
	
	size_t NSamples() const { return values.size() / frequencies.size(); }
	
	const double* FreqIndexValues(size_t freqIndex) const { return &values[freqIndex * NSamples()]; }
	
	void EvaluatePolynomial(const ao::uvector<double>& coefficients, bool invert);
	
	void Render(class PrimaryBeamImageSet& beamImages,
		double pixelScaleX, double pixelScaleY, 
		double phaseCentreRA, double phaseCentreDec,
		double pointingRA, double pointingDec,
		double phaseCentreDL, double phaseCentreDM,
		double frequencyHz) const;
	
	void Render(std::complex<float>* aterm,
		size_t width, size_t height,
		double pixelScaleX, double pixelScaleY, 
		double phaseCentreRA, double phaseCentreDec,
		double pointingRA, double pointingDec,
		double phaseCentreDL, double phaseCentreDM,
		double frequencyHz) const;
	
private:
	// Only works when frequencies.size() > 1
	ao::uvector<double> interpolateValues(double freq) const;
	// Works for any frequencies.size(), including when 1
	const double* interpolateValues(double frequencyHz, ao::uvector<double>& interpolatedVals) const;
	
	double lmMaxSquared(double frequencyHz) const;
};

#endif
