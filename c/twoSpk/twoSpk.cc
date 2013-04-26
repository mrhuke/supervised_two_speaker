#include "gammaTone.h"
#include "frontEnd.h"
#include "acoust.h"
#include "acoustDynamic.h"
#include "modelAdapt.h"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <fstream>

int main(int argc, const char **argv)
{
	bool arg = (argc==9);
	arg = arg || (argc==10 && atoi(argv[3])==4);

	if ( !arg )
	{
		cout << "Usage : twoSpk iSNR nChan Ver Sid1 Sid2 nGau Input (bW) Output" << endl << endl;
		cout << "iSNR  : Input SNR of the mixture (in dB)" << endl;
		cout << "nChan : Number of gammatone filters (128)" << endl;
		cout << "Ver   :-1 - State estimation based on premixed clean utterances" << endl;
		cout << "        0 - Reddy & Raj'07 system (matched train and test conditions)" << endl;
		cout << "        1 - MMSE estimation based on (Eq. 25) in Reddy & Raj'07" << endl;
		cout << "        2 - MAP estimation" << endl;
		cout << "        3 - Acoustic dynamics + MAP" << endl;
		cout << "        4 - Acoustic dynamics (BS) + MAP" << endl;
		cout << "Sid1  : Speaker identity (1-34) of the target" << endl;
		cout << "Sid2  : Speaker identity (1-34) of the interferer" << endl;
		cout << "nGau  : # of Gaussians (8,16,32,64,128 or 256) in a GMM" << endl;
		cout << "Input : Name of the input signal file (in ASCii format)" << endl;
		cout << "bW    : Beam width for viterbi decoding (optional)" << endl;
		cout << "Output: Name of the output soft mask" << endl << endl;
		cout << "Note  : - This implementation is in the cochleagram domain using log-cochleagram as features" << endl;
		cout << "        - GMMs are trained for 34 speakers in SSC corpus" << endl;
		cout << "        - Models are trained 16 kHz signals; 128-channel gammatone filterbank with range [50 8000] Hz" << endl;
		cout << "        - Gammatone filtering and cochleagram results match those of the MATLAB code"<<endl;		
		exit(1);
	}

	double isnr = atof(argv[1]);	
	unsigned int nChan = atoi(argv[2]);
	double ver = atof(argv[3]);
	int sid[2];
	sid[0] = atoi(argv[4]);
	sid[1] = atoi(argv[5]);
	unsigned int nGau = atoi(argv[6]);
	string fn(argv[7]);
	string outfn, pcfn1, pcfn2;
	int bW;
	// handling input files
	switch ( static_cast<int>(ver) )
	{
		case -1:
		case 0:
		case 1:
		case 2:
		case 3:
			outfn.assign(argv[8]);
			break;		
		case 4:
			outfn.assign(argv[9]);
			bW = atoi(argv[8]);
			break;
	}
        
	// read GMMs	
	stringstream tmp;
	tmp << "mdl/60dB." << nChan << ".hmmdefs." << nGau;
	string fn_gmm(tmp.str());
	CGMM *mdl = new CGMM(fn_gmm, 34);

	// read input signal
	unsigned int nSample;
	double *sig = readFileAsciiDouble(fn, nSample);

	// spkear models	
	GMM gmm1 = mdl->allGMM[sid[0]-1];
	GMM gmm2 = mdl->allGMM[sid[1]-1];
	
	// model adaptation based on estimated input SNR
	if (ver){
		modelAdapt mdl_adaptor(isnr); // construct with estimated input SNR
		//mdl_adaptor.adapt(sig, nSample, gmm2);
		mdl_adaptor.adapt(sig, nSample, gmm1, gmm2);
	}

	// gammatone filtering
	CGammaToneFilterBank *gBank = new CGammaToneFilterBank(nChan,50,8000,16000,1);
	gBank->filtering(sig,nSample);

	// feature extraction
	CFrontEnd *fe = new CFrontEnd(gBank);
	fe->getCochleagram_MATLAB();

	// log-amplitude
	vector< vector<double> > f;
	f.resize(fe->clg.size());
	for (vector<double>::size_type c=0; c!=fe->clg.size(); ++c)
	{
		f[c].resize(fe->clg[c].size());
		for (vector<double>::size_type m=0; m!=fe->clg[c].size(); ++m)
			f[c][m] = log(fe->clg[c][m]+EPS);
	}


	///////////////////////////////////////////////
	// Separation

	// construct a acoustic object with basic pdf & cdf values	
	AcoustDym acoust(f,gmm1,gmm2);

	// state sequence detection
	cout<<"State detection by "<<flush;
	switch (static_cast<int>(ver))
	{
		case -1: 		     // prior states from clean utterances
			cout<<"premixed utterances..."<<endl;
			acoust.acousticProbs();
			acoust.pri_gauss(fn,1);
			acoust.pri_gauss(fn,2);
			break;
		case 0:
		case 1:		             // MMSE = no hidden states
			cout<<"MMSE..."<<endl;
			break;	
		case 2:
			cout<<"MAP..."<<endl;
			acoust.acoust_gauss(); // MAP (indepednent frames)
			break;
		case 3:
			cout<<"Viterbi with acoust dynamics..."<<endl;
			acoust.viterbi_gauss();	// temporal + MAP
			break;
		case 4:
			cout<<"Viterbi (bW="<<bW<<") with acoust dynamics..."<<endl;
			acoust.viterbi_gauss(bW);  // temporal + beam search (width=bW)
                        break;
	}

	// mask estimation	
	cout<<"Mask estimation..."<<endl;
	switch (static_cast<int>(ver))
        {
		case 0:	
		case 1:
			acoust.mmse_est(); // only for MMSE
			break;
		case -1:
		case 2:
		case 3:
		case 4:
			acoust.map_est();	
			break;
	}
	cout<<endl;

	// output mask
	ofstream of(outfn.c_str());
	if (of.is_open())  of<<acoust<<flush;
	else { cerr<<"Opening file "<<outfn<<" failed."<<endl; return 1; }
	of.close();
	of.clear();
	
	// output joint search likelihood (only for the "searchSNR" method)
	if (ver==4.1){
		string fit = outfn+".fitness";
		of.open(fit.c_str());
		of<<acoust.fitness()<<flush;
		of.close();
	}	
	
	return 0;
}  
