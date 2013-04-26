#include "acoust.h"
#include "tool.h"
#include "gammaTone.h"
#include "frontEnd.h"
#include "gmm.h"

#include <algorithm>
#include <iterator>
#include <numeric>

ostream& operator<<(ostream& os, const Acoust& a)
{
        for (vector< vector<double> >::const_iterator it1=a.p.begin(); it1!=a.p.end(); ++it1){
                for(vector<double>::const_iterator it2=it1->begin(); it2!=it1->end(); ++it2)
                        os<<*it2<<" ";
                os<<endl;
        }
        return os;
}

Acoust::Acoust(const vector2D<double>::type &fea, GMM &g1, GMM &g2) : 
		f(fea), gmm1(g1), gmm2(g2)
{
	nFrame = f[0].size();
	nChan = f.size();
	nGau = gmm1.nGau;  // assume gmm1.nGau = gmm2.nGau

	state.resize(nFrame);
	px.resize(nFrame);
	py.resize(nFrame);
	cx.resize(nFrame);
	cy.resize(nFrame);
	p.resize(nFrame);
	p_ij_z.resize(nFrame);	
}

void Acoust::acousticProbs()
{
	cout<<"Compute acoustic probabilities (pdfs & cdfs)..."<<flush;

	for (unsigned int m=0; m!=nFrame; ++m)
        {
	        processBar(m,nFrame-1);

		p[m].resize(nChan);

		// pdf & cdf
        	px[m].resize(gmm1.nGau);
	        py[m].resize(gmm2.nGau);
	        cx[m].resize(gmm1.nGau);
        	cy[m].resize(gmm2.nGau);
		for (unsigned int c=0; c!=nChan; ++c)
	        {
                        for (unsigned int i=0; i<gmm1.nGau; i++)
                        {
                                double tmp = gmm1.gauss[i].calpdf(f[c][m],c);
                                px[m][i].push_back(exp(tmp));
                                tmp = gmm1.gauss[i].calcdf(f[c][m],c);
                                cx[m][i].push_back(tmp);
                        }
                        for (unsigned int j=0; j<gmm2.nGau; j++)
                        {
                                double tmp = gmm2.gauss[j].calpdf(f[c][m],c);
                                py[m][j].push_back(exp(tmp));
                                tmp = gmm2.gauss[j].calcdf(f[c][m],c);
                                cy[m][j].push_back(tmp);
                        }
                }

		// posterior prob.
		double pz0=0;
		vector3D<double>::type pz;
		pz.resize(gmm1.nGau);
		p_ij_z[m].resize(gmm1.nGau);
		for (unsigned int i=0; i<gmm1.nGau; i++)
	        {
			pz[i].resize(gmm2.nGau);
                        p_ij_z[m][i].resize(gmm2.nGau);
                        for (unsigned int j=0; j<gmm2.nGau; j++)
                        {
				pz[i][j].resize(nChan);
				double tmpProd = 1;
                                for (unsigned int d=0; d<nChan; d++)
                                {
                		        pz[i][j][d] = px[m][i][d]*cy[m][j][d] + py[m][j][d]*cx[m][i][d];
                                        tmpProd *=pz[i][j][d];
                                }
                                p_ij_z[m][i][j] = gmm1.gauss[i].weight * gmm2.gauss[j].weight * tmpProd;
                                pz0 += p_ij_z[m][i][j];
                         }
                }
		for (unsigned int i=0; i<gmm1.nGau; i++)
                {
                        for (unsigned int j=0; j<gmm2.nGau; j++)
			{
				p_ij_z[m][i][j] /= (pz0+EPS);
			}
		}
	}
}

void Acoust::mmse_est()
{
	acousticProbs();

	for (vector<double>::size_type m=0; m!=p_ij_z.size(); ++m)
        {
		processBar(m,p_ij_z.size()-1);
		for (unsigned int i=0; i<gmm1.nGau; i++)
	        {
			for (unsigned int j=0; j<gmm2.nGau; j++)
			{
                	        for (unsigned int d=0; d<nChan; d++)
				{
					double upper = px[m][i][d]*cy[m][j][d];
					double lower = px[m][i][d]*cy[m][j][d] + py[m][j][d]*cx[m][i][d];
                        		p[m][d] += p_ij_z[m][i][j] * ( upper/(lower+EPS) );                     
				}
			}
                }
       }
}

void Acoust::acoust_gauss()
{
	acousticProbs();

	for (vector<double>::size_type m=0; m!=p_ij_z.size(); ++m)
	{
		double max_value = 0;
		vector<double>::iterator max_ind;
		for (vector< vector<double> >::iterator it=p_ij_z[m].begin(); it!=p_ij_z[m].end(); ++it)
		{
			max_ind = max_element(it->begin(), it->end());
	                if (*max_ind>max_value){                                                
				max_value = *max_ind;
                	        state[m].first = static_cast<int>(it-p_ij_z[m].begin());
                        	state[m].second = static_cast<int>(max_ind-it->begin());
	                }
	        }
	}
	cout<<endl;
}

void Acoust::pri_gauss(string fileRoot, int spkrInd)
{
	unsigned int nSample;

	// hand-coded file names
	stringstream ss;
	ss<<fileRoot.substr(0,fileRoot.size()-5);
	ss<<"."<<spkrInd<<".val2";

        double *sig = readFileAsciiDouble(ss.str(), nSample);
        normWavSig(sig, nSample, 60);

	cout<<endl;
	CGammaToneFilterBank *gBank = new CGammaToneFilterBank(nChan,50,8000,16000,1);
        gBank->filtering(sig,nSample);

	CFrontEnd *fe = new CFrontEnd(gBank);
        fe->getCochleagram_MATLAB();

	vector< vector<double> > f;
	f.resize(fe->clg.size());
        for (vector<double>::size_type c=0; c!=fe->clg.size(); ++c)
        {
        	f[c].resize(fe->clg[c].size());
                for (vector<double>::size_type m=0; m!=fe->clg[c].size(); ++m)
                	f[c][m] = log(fe->clg[c][m]+EPS);
        }

	GMM gmm;
	if (spkrInd==1)
		gmm = gmm1;
	else if (spkrInd==2)
		gmm = gmm2;

	for (unsigned int m=0; m!=fe->nFrame; ++m)
        {		
 	        vector< vector<double> > p;
                vector<double> upper;
                p.resize(gmm.nGau);
		upper.resize(gmm.nGau);
	        for (unsigned int c=0; c!=fe->nChan; ++c)
                {
			for (unsigned int i=0; i<gmm.nGau; i++)
                        {
                	        double tmp = gmm.gauss[i].calpdf(f[c][m],c);
                                p[i].push_back(exp(tmp));
                        }
                }
		for (unsigned int i=0; i<gmm.nGau; ++i)
                {
                	double tmpProd = 1;
                        for (unsigned int c=0; c<fe->nChan; c++)
                        	tmpProd *= p[i][c];
                        upper[i] = gmm.gauss[i].weight * tmpProd;
                }
                vector<double>::iterator tmp = max_element(upper.begin(), upper.end());

		if (spkrInd==1)
			state[m].first = static_cast<int>(tmp-upper.begin());
		else if (spkrInd==2)
			state[m].second = static_cast<int>(tmp-upper.begin());
	}

	delete sig;
	delete gBank;
}

void Acoust::map_est()
{
	for (vector<double>::size_type m=0; m!=p_ij_z.size(); ++m)
        {
		p[m].resize(nChan);
		int i=state[m].first;
		int j=state[m].second;
	        for (unsigned int d=0; d<nChan; d++)
			p[m][d] = px[m][i][d] * cy[m][j][d] / (px[m][i][d]*cy[m][j][d]+py[m][j][d]*cx[m][i][d]+EPS);
	}
}
