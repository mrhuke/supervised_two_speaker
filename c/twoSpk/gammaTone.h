#ifndef MYGAMMATONE_H
#define MYGAMMATONE_H

#include "tool.h"

#define BW_CORRECTION       1.019      			/* ERB bandwidth correction 4th order */
#define DB 60									/* default average loudness of the utterance */
#define MIDDLE_EAR_SIZE 29

struct gammaTonePara
{
	float lCF, rCF;
	int nChan;
	int sf;
};

class CGammaToneFilter
{
	double rPart[4], iPart[4];

	double dR, dI, ddR, ddI;

	double gain, twoPiT, z, f1, f2, cf, bw;

	int delay;

	void initFilter(void);

	void oneStep(double input);

public:

	CGammaToneFilter(double cf, double bc, double sf, double midEarCoeff=1);

	~CGammaToneFilter(){;}

	double loudnessLevelInPhons(double dB, double freq);

	void filtering(double *input, unsigned int sigLength, double *response);
};

class CGammaToneFilterBank
{
	CGammaToneFilter **gf;	// each filter impulse responses

	double lowerCF, upperCF;		// cutoff frequencies

	double HzToERBRate(double Hz){ return( 21.4*log10(Hz*0.00437 + 1.0) ); }

	double ERBRateToHz(double ERBRate){ return( (pow(10, ERBRate/21.4) - 1) / 0.00437 ); }

	// equal-loudness contours
	double f[MIDDLE_EAR_SIZE];
	double af[MIDDLE_EAR_SIZE];
	double bf[MIDDLE_EAR_SIZE];
	double tf[MIDDLE_EAR_SIZE];
	void initMidEar();

public:

	unsigned int nChan, nSample;

	double sf;
	double **response;
	double *cf;

	CGammaToneFilterBank(unsigned int numChan=128, double lCF=50, double uCF=8000, double sampleFreq=16000, bool midEar=true);

	~CGammaToneFilterBank();

	void filtering(double *input, unsigned int sigLength);

	double* synthesis(double *input, vector< vector<double> > mask, unsigned int fRate=160, unsigned int winLength=320, double sf=16000, double lf=50, double hf=8000, bool bMidEar=true);

	double loudnessLevelInPhons(double dB, double freq);	// equal-loudness contours
};


#endif
