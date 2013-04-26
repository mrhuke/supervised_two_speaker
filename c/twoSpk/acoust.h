#ifndef ACOUST_H
#define ACOUST_H

#include "gmm.h"
#include "tool.h"

#include <utility>
#include <vector>
using namespace std;

#define DBLMAX 1.7976931348623158e+308 /* max double value */
#define EPS 2.2250738585072014e-308 /* min positive double value */

class Acoust
{
protected:
	typedef pair<int, int> gState;

	GMM gmm1, gmm2;
	const vector2D<double>::type f;  // log-cocheagram
	vector<gState> state;  // detected Gaussian sequence

	vector3D<double>::type px, py, cx, cy;
	vector3D<double>::type p_ij_z;
	vector2D<double>::type p;

public:
	
	int nFrame, nChan, nGau;

	Acoust(const vector2D<double>::type&, GMM&, GMM&); // constructor

	void acousticProbs();
	void acoust_gauss();	// Maximum a posterior estimate of Gaussian states
	void pri_gauss(string file, int spkInd);	// Detect states using clean utterances (ideal)

	void mmse_est();	// MMSE for mask estimation
	void map_est();		// MAP for mask estimation

	friend ostream& operator<<(ostream&, const Acoust&); //output soft mask (the 2D-vector 'p')
};

#endif
