#ifndef FRONTEND_H
#define FRONTEND_H

#include "gammaTone.h"
#include "filter.h"

#include <iostream>
#include <vector>
using namespace std;

#define MAX_DELAY   0.0125					// maxmum time-lag = 12.5 ms
#define MIN_D 0.002
#define MAX_D 0.015

struct acfCompute
{
        double sum[662];
        double sumS[662];
        double inputXR[2048];
        double inputXI[2048];
        double inputXWR[2048];
        double inputXWI[2048];
};

class CFrontEnd{

	double **response, **envelope;  // filtered signals
	double sf, *cf, acfLen, acfOrder;
	unsigned int fRate, winLen, nSample;				// fRate and winLen are in samples

	acfCompute data;
        void fftACF(double *, int);
        void getEnvelope();

	double min_delay, max_delay;

public:
	vector<double>::size_type nFrame, nChan;
	
	vector< vector<double> > clg;                                                   // cochleagram
        vector< vector<vector<double> > > acf;                                  // correlogram, #frame x #channel x #delay
        vector< vector<vector<double> > > evAcf;                                // correlogram, #frame x #channel x #delay
        vector< vector<double> > zc;                                                    // zero-crossing rates
        vector< vector<double> > evZc;                                          // envelope zero-crossing rates
        vector< vector<double> > cc;                                                    // cross-channel correlation
        vector< vector<double> > evCc;                                          // envelope cross-channel correlation
        vector< vector< vector<double> > > f6;                                  // 6-D features
	
	CFrontEnd(CGammaToneFilterBank *gBank, double r=0.01, double w=0.02);

	void getCochleagram_MATLAB();						// cochleagram (MatLab version)
	void getCrossChanCorr();                                                        // normalized cross-channel correlation
        void getCorrelogram();
        void getZC();
        void get6f(vector<double> &pitch);                                      // Guoning's 6-dimensional feature
};

#endif
