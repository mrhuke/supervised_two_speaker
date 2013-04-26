#include "modelAdapt.h"

// normalize mixture to 60 dB and adapt both gmm1 and gmm2
void modelAdapt::adapt(double* sig, size_t nSample, GMM &gmm1, GMM &gmm2)
{
	normWavSig(sig, nSample, 60);

	double alpha = pow(10,isnr/20);
        double I2 = 60-10*log10(1+alpha*alpha);
        double I1 = I2 + isnr;

 	double gain1 = log(pow(10,(I1-60)/10));  //note the original model is trained using 60-dB mixtures
        for (vector<double>::size_type k=0; k!=gmm1.nGau; ++k){
		for (vector<double>::size_type c=0; c!=gmm1.gauss[k].vecSize; ++c)
                	gmm1.gauss[k].me[c] = gain1 + gmm1.gauss[k].me[c];
        }

	double gain2 = log(pow(10,(I2-60)/10));
        for (vector<double>::size_type k=0; k!=gmm2.nGau; ++k){
		for (vector<double>::size_type c=0; c!=gmm2.gauss[k].vecSize; ++c){
                	gmm2.gauss[k].me[c] = gain2 + gmm2.gauss[k].me[c];
                }
        }

	#ifdef DEBUG
        cout<<"SNR="<<isnr<<",I1="<<I1<<",I2="<<I2<<endl;
        cout<<"gain1="<<gain1<<",gain2="<<gain2<<endl<<endl;
        #endif
}


// keep gmm1 at 60-dB, change gmm2 and mixture
void modelAdapt::adapt(double *sig, size_t nSample, GMM &gmm2)
{
	double I2 = 60 - isnr;	// gmm1 is assumed to 60-dB
	double gain2 = log(pow(10,(I2-60)/10));  // again, original models are trained using 60-dB mixtures
	for (vector<double>::size_type k=0; k!=gmm2.nGau; ++k){
        	for (vector<double>::size_type c=0; c!=gmm2.gauss[k].vecSize; ++c){
                	gmm2.gauss[k].me[c] = gain2 + gmm2.gauss[k].me[c];
                }
        }

	double I = 60 + 10*log10(1+pow(10,-isnr/10));
        normWavSig(sig, nSample, I);

	#ifdef DEBUG
        cout<<"SNR="<<isnr<<",I="<<I<<endl;
        cout<<"gain2="<<gain2<<endl<<endl;
        #endif
}
