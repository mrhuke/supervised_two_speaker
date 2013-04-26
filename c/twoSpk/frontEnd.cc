#include "frontEnd.h"

CFrontEnd::CFrontEnd(CGammaToneFilterBank *gBank, double r, double w)
{	
	sf = gBank->sf;
	cf = gBank->cf;
	nChan = gBank->nChan;		
	fRate = static_cast<unsigned int>(sf*r);
        winLen = static_cast<unsigned int>(sf*w);
	nFrame = static_cast<unsigned int>(gBank->nSample/fRate);	
	nSample = gBank->nSample;	

	response = gBank->response;
        envelope = 0;
}

void CFrontEnd::getCochleagram_MATLAB()
{
	cout<<"Getting Cochleagram...";

	// initialize Cochleagram
	clg.resize(nChan);
	for(vector<double>::size_type i=0; i!=clg.size(); ++i){
		clg[i].resize(nFrame);
	}

	int increment = static_cast<int>(winLen/fRate);
	for (vector<double>::size_type c=0; c!=clg.size(); ++c){
		processBar(c,clg.size());
		for (vector<double>::size_type m=0; m!=clg[c].size(); ++m){
			if (m<increment-1)
			{
				for (vector<double>::size_type t=0; t<fRate; ++t)
					if (t>=0 && t<nSample)
						clg[c][m] += response[c][t]*response[c][t];
			}
			else
			{
				for (vector<double>::size_type t=(m-increment+1)*fRate; t<(m-increment+1)*fRate+winLen; ++t)
				{	
					if (t>=0 && t<nSample)
						clg[c][m] += response[c][t]*response[c][t];
				}
			}			
		}
	}
	cout<<"Done!"<<endl;
}

void CFrontEnd::fftACF(double *response, int frame)
{
        int step, delay;
        int window = winLen;
        int sigLength = nSample;

        for(step=0; step<acfLen; step++)
                data.inputXR[step] = data.inputXI[step] = data.inputXWR[step] = data.inputXWI[step] = 0;

        for(step=0; step<(window+max_delay); step++)
        {
                int tim = (frame+2) * window/2 - (step+1);
                if ( (tim>= 0)  && (tim<sigLength) )
                {
                        data.inputXR[step] = response[tim];

                        if(step<window) data.inputXWR[step] = data.inputXR[step];
                }
        }

        data.sum[0]=data.sumS[0]=0;
        for(step=0; step<window; step++)
        {
                data.sumS[0] += data.inputXR[step]*data.inputXR[step]/double(window);
        }

        for(step=1; step<max_delay; step++)
        {
                double f1 = data.inputXR[step+window-1] - data.inputXR[step-1];
                double f2 = data.inputXR[step+window-1] + data.inputXR[step-1];

                data.sumS[step] = data.sumS[step-1] + f1*f2/double(window);

                if (data.sumS[step]<0) data.sumS[step]=0;
        }

        fft(data.inputXR, data.inputXI, acfOrder, 1);
        fft(data.inputXWR, data.inputXWI, acfOrder, 1);

        for(step=0; step<acfLen; step++)
        {
                double f1 = data.inputXR[step]*data.inputXWR[step] + data.inputXI[step]*data.inputXWI[step];
                double f2 = -data.inputXR[step]*data.inputXWI[step] + data.inputXI[step]*data.inputXWR[step];

                data.inputXR[step]=f1;
                data.inputXI[step]=f2;
        }

        fft(data.inputXR, data.inputXI, acfOrder, -1);

        for(delay=0; delay<max_delay; delay++)
        {
                data.inputXR[delay] /= double(acfLen*window);
                data.inputXR[delay] /= double(sqrt(data.sumS[0]*data.sumS[delay]+1e-5*data.sumS[0])) / data.sumS[0]; //+1e-100));
        }
}

void CFrontEnd::getEnvelope()
{
        envelope = new double*[nChan];
        for(int chan=0; chan<nChan; chan++)
                envelope[chan]=new double[nSample];

        double bP1 = 50, bP2 = 450, bPTs = 20;
        CFilter *bandPass = new CFilter( (bP1+bP2)/sf, (bP2-bP1)/sf, bPTs*2/sf, sf, 0);
        for(int chan=0; chan<nChan; chan++)
        {
                double *temp = new double[nSample];
                for(int n=0; n<nSample; n++)
                        temp[n] = (response[chan][n]>0) ? response[chan][n]:0;
                bandPass->filtering(temp, envelope[chan], nSample);
                delete[] temp;
        }
        delete bandPass;
}

void CFrontEnd::getCorrelogram()
{
        min_delay = sf*MIN_D;   // min_delay = 2 ms
        max_delay = sf*MAX_D;   // max_delay = 15 ms

        acfLen = max_delay + winLen;
        acfOrder = ceil(log(acfLen)/log(2.0));
        acfLen = pow(2.0, acfOrder);

        if (envelope==0)
                getEnvelope();

        evAcf.resize(nFrame);
        acf.resize(nFrame);
        for(int n=0; n<nFrame; n++)
        {
                evAcf[n].resize(nChan);
                acf[n].resize(nChan);
                processBar(n,nFrame-1);
                for(int chan=0; chan<nChan; chan++)
                {
                        fftACF(envelope[chan], n);
                        for(int m=0; m<max_delay; m++)
                                evAcf[n][chan].push_back(data.inputXR[m]);

                        fftACF(response[chan], n);
                        for(int m=0; m<max_delay; m++)
                                acf[n][chan].push_back(data.inputXR[m]);
                }
        }
}

void CFrontEnd::getZC()
{
        if (envelope==0)
                getEnvelope();

        double max_delay = sf*MAX_D;
        evZc.resize(nFrame);
        zc.resize(nFrame);
        for(int n=0; n<nFrame; n++)
        {
                for(int chan=0; chan<nChan; chan++)
                {
                        evZc[n].push_back(zeroCross(evAcf[n][chan], max_delay));
                        zc[n].push_back(zeroCross(acf[n][chan], max_delay));
                }
        }
}  

void CFrontEnd::get6f(vector<double> &pitch)
{

	cout<<"."<<flush;
        getEnvelope(); 

	 cout<<"."<<flush;
        getCorrelogram();

	cout<<"."<<flush;
        getZC();

	f6.resize(nFrame);
        for (vector<double>::size_type m=0; m!=f6.size(); m++){
                f6[m].resize(nChan);
                for (vector<double>::size_type c=0; c!=f6[m].size(); c++)
                {
                        f6[m][c].resize(6);

                        f6[m][c][0] = acf[m][c][pitch[m]]/acf[m][c][0];
                        f6[m][c][1] = evAcf[m][c][pitch[m]]/evAcf[m][c][0];
                        f6[m][c][2] = (pitch[m]+1)/(zc[m][c]+1e-10)/2;
                        f6[m][c][3] = (pitch[m]+1)/(evZc[m][c]+1e-10)/2;
                        f6[m][c][4] = static_cast<int>(f6[m][c][2]); if (f6[m][c][2]-f6[m][c][4]>=.5) f6[m][c][4]++;
                        f6[m][c][5] = static_cast<int>(f6[m][c][3]); if (f6[m][c][3]-f6[m][c][5]>=.5) f6[m][c][5]++;
                        f6[m][c][2] = fabs(f6[m][c][2]-f6[m][c][4]);
                        f6[m][c][3] = fabs(f6[m][c][3]-f6[m][c][5]);
                        f6[m][c][4] *= 0.05; if (f6[m][c][4]>1) f6[m][c][4]=1; if (f6[m][c][4]==0) f6[m][c][4]=2;
                        f6[m][c][5] *= 0.05; if (f6[m][c][5]>1) f6[m][c][5]=1; if (f6[m][c][5]==0) f6[m][c][5]=2;
                }
        }
}
