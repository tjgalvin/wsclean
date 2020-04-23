#ifndef VLA_BEAM_H
#define VLA_BEAM_H

class VLABeam
{
public:
	void GetCoefficients();
	
private:
	void getCoefficients(double freq);
	static char determineFeed(double freq, double freqCenter = 0.0);
	static void limitFreqForBand(char band, double& freq)
};

#endif
