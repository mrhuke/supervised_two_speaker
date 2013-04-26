#ifndef ACOUSTDYNAMIC_H
#define ACOUSTDYNAMIC_H

#include "acoust.h"
#include "gmm.h"

#include <utility>
#include <vector>
using namespace std;

class AcoustDym: public Acoust
{
	typedef vector< vector<gState> > gState2D;  // a 2-D likelihood matrix storing frame-level observation function
	
	vector3D<double>::type pb;  // Viterbi probabilities
	vector<gState2D> gPath;  // 2-D Viterbi path (only detecting Gaussian states)

	vector<double> pi1, pi2; // Gaussian priors
	vector2D<double>::type gt1, gt2; // Gaussian trans. prob.
	
	int bW;
	vector<gState> gbeam;
	
	void initPriors(GMM &, GMM &); // load acoustic dynamics and other pre-trained info

public:
	typedef pair<double,gState> p_n_gState;
	
	// constructor
	AcoustDym(const vector2D<double>::type &f, GMM &gmm1, GMM &gmm2) : Acoust(f, gmm1, gmm2) 
		{initPriors(gmm1,gmm2);}
	
	// Viterbi decoding	
	void viterbi_gauss();	// exhaust Viterbi search for Gaussian states
	void viterbi_gauss(int bW);	// beam search

	double fitness();
};

// comparators for beam search
bool icompare_p_n_gstate(const AcoustDym::p_n_gState &, const AcoustDym::p_n_gState &);

#endif
