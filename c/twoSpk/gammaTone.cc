#include "gammaTone.h"

CGammaToneFilterBank::CGammaToneFilterBank(unsigned int numChan, double lCF, double uCF, double sampleFreq, bool bMidEar)
{
	nChan = numChan;
	lowerCF = lCF;
	upperCF = uCF;
	sf = sampleFreq;

	// allocate memory
	response = new double*[nChan];
	cf = new double[nChan];

	double lowerERB = HzToERBRate(lowerCF);
	double upperERB = HzToERBRate(upperCF);

	double spaceERB = (nChan > 1) ? (upperERB-lowerERB)/(nChan-1) : 0;
	
	// midear adjustment?
	if (bMidEar) {initMidEar();}
	
	gf = new CGammaToneFilter*[nChan];
	for (int chan=0; chan<nChan; chan++)
	{
		cf[chan] = ERBRateToHz(lowerERB + chan*spaceERB);
		if (bMidEar){			
			// use mid-ear coefficients
			gf[chan] = new CGammaToneFilter(cf[chan], BW_CORRECTION, sf, pow(10,(loudnessLevelInPhons(60,cf[chan])-60)/20));
		}
		else
			// not use mid-ear coefficients
			gf[chan] = new CGammaToneFilter(cf[chan], BW_CORRECTION, sf);
	}
}

CGammaToneFilter::CGammaToneFilter(double centerF, double bC, double sf, double midEarCoeff)
{
	cf = centerF;
	bw = 24.7 * ( cf*0.00437 + 1.0) * bC;
	
	// phase compensation
//  double d = 1.5/PI/bw * double(sf);	
//	if (fmod(d, 1)>=0.5) delay=floor(d)+1;
//	else delay=floor(d);

	delay = 0;
	
	twoPiT = 2*PI/sf;	
	gain = midEarCoeff * pow(twoPiT*bw, 4) / 3.0; // mid-ear adjustment
	
	z = exp(-twoPiT * bw);

	f1 = cos(cf * twoPiT) * z;
	f2 = sin(cf * twoPiT) * z;
}

void CGammaToneFilter::initFilter(void)
{
	for (int n=0; n<4; n++)
	{
		rPart[n] = 0;
		iPart[n] = 0;
	}
}

void CGammaToneFilter::oneStep(double input)
{
	double x[4], y[4];
	
	for (int i=0; i<4; i++)
	{
		x[i] = f1*rPart[i] - f2*iPart[i];
		y[i] = f2*rPart[i] + f1*iPart[i];
	}

	rPart[0] = input * f1 + x[0];
	iPart[0] = input * f2 + y[0];
		
	rPart[1] = rPart[0] + x[1];
	iPart[1] = iPart[0] + y[1];
		
	rPart[2] = rPart[1] + x[1] + x[2];
	iPart[2] = iPart[1] + y[1] + y[2];
		
	rPart[3] = rPart[2] + x[1] + 2*x[2] + x[3];
	iPart[3] = iPart[2] + y[1] + 2*y[2] + y[3];

	dR = 3 * rPart[2] - twoPiT*bw * rPart[3] - twoPiT*cf* iPart[3];
	dI = 3 * iPart[2] - twoPiT*bw * iPart[3] + twoPiT*cf* rPart[3];
}

void CGammaToneFilter::filtering(double *input, unsigned int sigLength, double *response)
{
	int k;

	initFilter();

	for (unsigned int n=0; n<sigLength; n++)
	{
		k = n-delay;
		if (k>=0) response[k] = gain * rPart[3];

		oneStep(input[n]);
	}	

	for (unsigned int n=0; n<delay; n++)
	{
		k = n-delay+sigLength;
		if (k>=0) response[k] = gain * rPart[3];

		oneStep(0);
	}
}

CGammaToneFilterBank::~CGammaToneFilterBank()
{
	for(unsigned int chan=0; chan<nChan; chan++)
	{
		delete [] response[chan];
		delete gf[chan];
	}
}

void CGammaToneFilterBank::filtering(double *input, unsigned int sigLength)
{
	cout<<"Gammatone filtering...";
	nSample = sigLength;

	for (unsigned int c=0; c<nChan; c++)
	{
		processBar(c,nChan-1);

		response[c] = new double[sigLength];
		gf[c]->filtering(input, sigLength, response[c]);
	}
}


double* CGammaToneFilterBank::synthesis(double *mixture, vector< vector<double> > mask, unsigned int fRate, unsigned int winLength, double sf, double lf, double hf, bool bMidEar)
{
	// first pass of filtering
	filtering(mixture, nSample);	

	// prepare for second pass
	double* coswin = new double[winLength];
	for (int t=0; t<winLength; t++)
		coswin[t] = (1 + cos(2*PI*t/winLength - PI))/2;
	double lowerERB = HzToERBRate(lowerCF);
	double upperERB = HzToERBRate(upperCF);
	double spaceERB = (nChan > 1) ? (upperERB-lowerERB)/(nChan-1) : 0;
	
	// midear adjustment?
	if (bMidEar) {initMidEar();}
	double *midEarCoeff = new double[nChan];
	for (int chan=0; chan<nChan; chan++)
	{
		if (bMidEar)
		{
			double cf = ERBRateToHz(lowerERB + chan*spaceERB);
			midEarCoeff[chan] = pow(10,(loudnessLevelInPhons(60,cf)-60)/20);
		}
		else
		{
			midEarCoeff[chan] = 1;
		}
	}
	
	// second pass filtering
	double **temp1 = new double*[nChan];
	double **temp2 = new double*[nChan];
	double *r = new double[nSample];
	for (unsigned int t=0; t<nSample; t++)
		r[t] = 0;
	for (int c=0; c<nChan; ++c)
	{
		// time reverse filter output & normalize out mid-ear coefficients
		temp1[c] = new double[nSample];
		for (unsigned int t=0; t<nSample; t++)
			temp1[c][t] = response[c][nSample-1-t]/midEarCoeff[c];
		
		// second pass filtering
		temp2[c] = new double[nSample];
		gf[c]->filtering(temp1[c], nSample, temp2[c]);

		// time reverse again & normalize out mid-ear coefficients
		for (unsigned int t=0; t<nSample; t++)
			temp1[c][t] = temp2[c][nSample-1-t]/midEarCoeff[c];

		// mask value can be binary or rational
		double *weight = new double[nSample];
		for (unsigned int t=0; t<nSample; t++)
			weight[t] = 0;
		int increment=winLength/fRate;
		for (unsigned int m=0; m<mask[1].size(); m++){            
			if (m<increment-1)
			{
				for (unsigned int t=0; t<(m+1)*fRate; t++)
					weight[t] = weight[t] + mask[c][m]*coswin[(increment-1-m)*fRate+t];
			}
			else
			{
				for (unsigned int t=(m-1)*fRate; t<(m-1)*fRate+winLength; t++)
					if(t<nSample)
						weight[t] = weight[t] + mask[c][m]*coswin[t-(m-1)*fRate];
			}
		}

		for (unsigned int t=0; t<nSample; t++)
			r[t] = r[t] + temp1[c][t]*weight[t];
	}
	
	return r;
}


void CGammaToneFilterBank::initMidEar()
{
	f[0]=20.0;     af[0]=2.347;  bf[0]=0.00561;   tf[0]=74.3;
	f[1]=25.0;     af[1]=2.190;  bf[1]=0.00527;   tf[1]=65.0;
	f[2]=31.5;     af[2]=2.050;  bf[2]=0.00481;   tf[2]=56.3;
	f[3]=40.0;     af[3]=1.879;  bf[3]=0.00404;   tf[3]=48.4;
	f[4]=50.0;     af[4]=1.724;  bf[4]=0.00383;   tf[4]=41.7;
	f[5]=63.0;     af[5]=1.579;  bf[5]=0.00286;   tf[5]=35.5;
	f[6]=80.0;     af[6]=1.512;  bf[6]=0.00259;   tf[6]=29.8;
	f[7]=100.0;    af[7]=1.466;  bf[7]=0.00257;   tf[7]=25.1;
	f[8]=125.0;    af[8]=1.426;  bf[8]=0.00256;   tf[8]=20.7;
	f[9]=160.0;    af[9]=1.394;  bf[9]=0.00255;   tf[9]=16.8;
	f[10]=200.0;   af[10]=1.372; bf[10]=0.00254;  tf[10]=13.8;
	f[11]=250.0;   af[11]=1.344; bf[11]=0.00248;  tf[11]=11.2;
	f[12]=315.0;   af[12]=1.304; bf[12]=0.00229;  tf[12]=8.9;
	f[13]=400.0;   af[13]=1.256; bf[13]=0.00201;  tf[13]=7.2;
	f[14]=500.0;   af[14]=1.203; bf[14]=0.00162;  tf[14]=6.0;
	f[15]=630.0;   af[15]=1.135; bf[15]=0.00111;  tf[15]=5.0;
	f[16]=800.0;   af[16]=1.062; bf[16]=0.00052;  tf[16]=4.4;
	f[17]=1000.0;  af[17]=1.000; bf[17]=0.00000;  tf[17]=4.2;
	f[18]=1250.0;  af[18]=0.967; bf[18]=-0.00039; tf[18]=3.7;
	f[19]=1600.0;  af[19]=0.943; bf[19]=-0.00067; tf[19]=2.6;
	f[20]=2000.0;  af[20]=0.932; bf[20]=-0.00092; tf[20]=1.0;
	f[21]=2500.0;  af[21]=0.933; bf[21]=-0.00105; tf[21]=-1.2;
	f[22]=3150.0;  af[22]=0.937; bf[22]=-0.00104; tf[22]=-3.6;
	f[23]=4000.0;  af[23]=0.952; bf[23]=-0.00088; tf[23]=-3.9;
	f[24]=5000.0;  af[24]=0.974; bf[24]=-0.00055; tf[24]=-1.1;
	f[25]=6300.0;  af[25]=1.027; bf[25]=0.00000;  tf[25]=6.6;
	f[26]=8000.0;  af[26]=1.135; bf[26]=0.00089;  tf[26]=15.3;
	f[27]=10000.0; af[27]=1.266; bf[27]=0.00211;  tf[27]=16.4;
	f[28]=12500.0; af[28]=1.501; bf[28]=0.00488;  tf[28]=11.6;
}


double CGammaToneFilterBank::loudnessLevelInPhons(double dB, double freq)
{
	int i=0;
	double afy, bfy, tfy;

	if ((freq<20.0) | (freq>12500.0)) {
		cerr<<"Cannot compute a outer/middle ear gain for that frequency."<<endl;
		return 0;
	}
	while (f[i] < freq) i++;
	afy=af[i-1]+(freq-f[i-1])*(af[i]-af[i-1])/(f[i]-f[i-1]);
	bfy=bf[i-1]+(freq-f[i-1])*(bf[i]-bf[i-1])/(f[i]-f[i-1]);
	tfy=tf[i-1]+(freq-f[i-1])*(tf[i]-tf[i-1])/(f[i]-f[i-1]);
	return 4.2+afy*(dB-tfy)/(1.0+bfy*(dB-tfy));
}
