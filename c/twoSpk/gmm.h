#ifndef MYGMM_H
#define MYGMM_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
using namespace std;

#define PI (3.1415926535897932384626433832795)

class Gaussian
{
public:
	vector<double> me;
	vector<double> var;
	double weight;
	double vecSize;
	double gConst;

	double calpdf(double,unsigned int);	// log-liklihood	
	double calcdf(double,unsigned int);

	double calpdf_vec(vector<double>);	// log-likelihood
};

struct GMM
{
	vector<Gaussian> gauss;
	string name;
	unsigned int nGau, dim;
};



class CGMM
{	

public:
	unsigned int spkrNum;
	vector<GMM> allGMM;
	
	CGMM(string, unsigned int); // construct with all HTK GMMs from a single file
};

#endif
