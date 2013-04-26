#ifndef FILTER_H
#define FILTER_H

#include "tool.h"

class CFilter
{
	int fs;
	int fLength, nOrder, nFFT;

	double beta;

	double *rValue;
	double *iValue;

	double *inputR;
	double *inputI;

	void hammingPara(double delta, double transBw, int &fLength);

	void hammingBandPass(double wc, double wn);
	void hammingDiffOneBandPass(double wc, double wn);
	void hammingDiffTwoBandPass(double wc, double wn);

	void kaiserPara(double delta, double transBw);

	void kaiserLowPass(double wn);
	void diffKaiserLowPass(double wn);

public:
	CFilter(double ws, double wn, double ts, int fs, int type);

	~CFilter();

	void filtering(double *input, double *output, int sigLength);
};

#endif