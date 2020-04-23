#ifndef VOLTAGE_PATTERN_H
#define VOLTAGE_PATTERN_H

#include "../uvector.h"

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
	
	ao::uvector<double> interpolateValues(double freq) const;
	
	void EvaluateInversePolynomial(const ao::uvector<double>& coefficients);
	
	void Render(class PrimaryBeamImageSet& beamImages,
		double pixelScaleX, double pixelScaleY, 
		double phaseCentreRA, double phaseCentreDec,
		double phaseCentreDL, double phaseCentreDM,
		double frequencyHz) const;
};

#endif
