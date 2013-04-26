#include "acoustDynamic.h"
#include "tool.h"
#include <string>
#include <numeric>
#include <algorithm>

void AcoustDym::initPriors(GMM &gmm1, GMM &gmm2)
{
	// Gaussian priors
	pi1.resize(gmm1.gauss.size());
	pi2.resize(gmm2.gauss.size());
	for (vector<Gaussian>::size_type k=0; k!=gmm1.gauss.size(); ++k)
	{
		pi1[k] = gmm1.gauss[k].weight;
		pi2[k] = gmm2.gauss[k].weight;
	}

	// transition probabilities
	string fn1(gmm1.name.substr(4)), fn2(gmm2.name.substr(4));	
	fn1 = "mdl/trans.256-"+fn1;  // hand-code filenames
	fn2 = "mdl/trans.256-"+fn2;
	cout<<"Loading transitions, "<<flush;
	ifstream in_str1(fn1.c_str()), in_str2(fn2.c_str());
	if (in_str1.is_open() && in_str2.is_open())	
	{
		in_str1>>gt1;
		in_str2>>gt2;
	}
	else
	{
		cerr<<"Open file "<<fn1<<","<<fn2<<" error!"<<endl;
		return;
	}
	// normalize
	for (vector< vector<double> >::size_type i=0; i!=gt1.size(); ++i)
	{	
		double sum1=accumulate(gt1[i].begin(),gt1[i].end(),0);
		double sum2=accumulate(gt2[i].begin(),gt2[i].end(),0);
		for (vector<double>::size_type j=0; j!=gt1[i].size(); ++j)
		{
			gt1[i][j] /= sum1;
			gt2[i][j] /= sum2;
		}
	}
	cout<<endl;
}

// Exhaustive search via 2-D Viterbi 
void AcoustDym::viterbi_gauss()
{
	acousticProbs();

	#ifdef DEBUG
	vector<double> pctPrune;
	#endif

	for(int m=0; m<nFrame; ++m)
	{
		processBar(m,nFrame-1);

		vector2D<double>::type pt;	// observation prob. until current frame
		pt.resize(pi1.size());

		gState2D curPath;	// current 2-D viterbi state path
		curPath.resize(pi1.size());
		
		// (u,v) -> (i,j)
		for (vector<double>::size_type i=0; i!=pi1.size(); ++i)
		{
			pt[i].resize(pi2.size());
			curPath[i].resize(pi2.size());
			for (vector<double>::size_type j=0; j!=pi2.size(); ++j)
			{
				if (m==0)
				{
					pt[i][j] = log(p_ij_z[m][i][j]); // forward viterbi log-liklihood 
					curPath[i][j].first = -1; // boundary values
					curPath[i][j].second = -1;
				}
				else
				{
					double maxP=-DBLMAX;
					gState maxS;
					vector< vector<double> > tmp;				
					tmp.resize(pi1.size());	
					
					// prune (u,v) in (u,v)->(i,j)
					vector< vector2D<double>::type::size_type > ind_us;
					vector< vector2D<double>::type::size_type > ind_vs;
					for (vector2D<double>::type::size_type k=0; k!=gt1.size(); ++k)
					{
						if (gt1[k][i]>0) // purne trans.
							ind_us.push_back(k);
                                                if (gt2[k][j]>0)
                                                        ind_vs.push_back(k);
                                        }
					#ifdef DEBUG		
					double r = static_cast<double>(ind_us.size())*ind_vs.size()/(pi1.size()*pi2.size());
				        pctPrune.push_back( 100*( 1 - r ) );
				        #endif

					// viterbi search
					for (int ind_u=0; ind_u!=ind_us.size(); ++ind_u)
					{
						vector< vector<double> >::size_type u = ind_us[ind_u];
						tmp[u].resize(pi2.size());
						for (int ind_v=0; ind_v!=ind_vs.size(); ++ind_v)
						{
							int v = ind_vs[ind_v];
							double t_part1 = log( gt1[u][i]+EPS );
							double t_part2 = log( gt2[v][j]+EPS );
							double o_part = log( p_ij_z[m][i][j]/(pi1[i]*pi2[j]+EPS) + EPS );
							tmp[u][v] = pb[pb.size()-1][u][v] + o_part + t_part1 + t_part2;
							if (tmp[u][v]>maxP)
							{
								maxP = tmp[u][v];
								maxS.first = u;
								maxS.second = v;
							}
						}
					}
					pt[i][j]=maxP;
					curPath[i][j]=maxS;
				}
			}
		}

		pb.push_back(pt);
		gPath.push_back(curPath);
	}	
	#ifdef DEBUG
	cout<<"Average state pruned percent: "<< accumulate(pctPrune.begin(),pctPrune.end(),0)/pctPrune.size() << "\%" << endl;
        #endif

	// back trace
	double maxP = -DBLMAX;
	for (vector<double>::size_type i=0; i!=pi1.size(); ++i)
        {
                for (vector<double>::size_type j=0; j!=pi2.size(); ++j)
                {
                	if (pb[nFrame-1][i][j]>maxP)
                        {
                        	maxP = pb[nFrame-1][i][j];
                                state[nFrame-1].first = i;
                                state[nFrame-1].second = j;
                        }
                }
        }
	for (int m=nFrame-1; m>0; m--)
		state[m-1] = gPath[m][state[m].first][state[m].second];
}

// Gaussian state detection using Viterbi search with a beam width
void AcoustDym::viterbi_gauss(int bW)
{
	acousticProbs();

	gbeam.resize(bW);

	for(int m=0; m<nFrame; ++m)
	{
		processBar(m,nFrame-1);

		gState2D curPath; // Viterbi path at current frame (store the most likely states of previous frame)
		curPath.resize(pi1.size());
	
		vector<p_n_gState> state2s; // store likelihoods-and-asscoiated-states for pruning
		
		vector2D<double>::type pt; // Viterbi probability, a 2-D matrix storing the accumulated likelihoods until the current frame
		pt.resize(pi1.size());
	
		// (u,v) -> (i,j)
		for (vector<double>::size_type i=0; i!=pi1.size(); ++i)
		{
			pt[i].resize(pi2.size());
			curPath[i].resize(pi2.size());
			for (vector<double>::size_type j=0; j!=pi2.size(); ++j)
			{
				if (m==0)
				{
					pt[i][j] = log(p_ij_z[m][i][j]); // forward viterbi log-liklihood 
					curPath[i][j].first = -1; // boundary values
					curPath[i][j].second = -1;
				}
				else
				{
					double maxP=-DBLMAX;
					gState maxS;
					maxS.first = -1;
					maxS.second = -1;
					double tmp;				
					
					// viterbi search (only search in beam)
				        for (vector<gState>::iterator it=gbeam.begin(); it!=gbeam.end(); ++it)
					{	
						int u = it->first;
						int v = it->second;
						// cond. trans. prob.
						double t_part1 = log( gt1[u][i]+EPS );
						double t_part2 = log( gt2[v][j]+EPS );
						double o_part = log( p_ij_z[m][i][j]/(pi1[i]*pi2[j]+EPS) + EPS );
						tmp = pb[pb.size()-1][u][v] + t_part1 + t_part2 + o_part;
						if (tmp>maxP)
						{
							maxP = tmp;
							maxS.first = u;
							maxS.second = v;
						}
					}
					pt[i][j]=maxP;
					curPath[i][j]=maxS;
				}
				// store states for pruning
				gState tmp;
                                tmp.first = i;
                                tmp.second = j;

                                p_n_gState curState2(pt[i][j],tmp);
                                state2s.push_back(curState2);
			}
		}

		pb.push_back(pt);
		gPath.push_back(curPath);

		// state pruning
		sort(state2s.begin(), state2s.end(), icompare_p_n_gstate);
                for (vector<gState>::size_type kk=0; kk!=gbeam.size(); ++kk)
                        gbeam[kk] = state2s[kk].second;
	}	

	// back trace
	double maxP = -DBLMAX;
	for (vector<double>::size_type i=0; i!=pi1.size(); ++i)
        {
                for (vector<double>::size_type j=0; j!=pi2.size(); ++j)
                {
                	if (pb[nFrame-1][i][j]>maxP)
                        {
                        	maxP = pb[nFrame-1][i][j];
                                state[nFrame-1].first = i;
                                state[nFrame-1].second = j;
                        }
                }
        }
	for (int m=nFrame-1; m>0; m--)
		state[m-1] = gPath[m][state[m].first][state[m].second];
}

double AcoustDym::fitness()
{
	return pb[nFrame-1][state[nFrame-1].first][state[nFrame-1].second];
}

bool icompare_p_n_gstate(const AcoustDym::p_n_gState &l, const AcoustDym::p_n_gState &r)
{
        return l.first > r.first;
}
