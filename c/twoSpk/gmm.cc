#include "gmm.h"

double Gaussian::calpdf(double data, unsigned int d)
{
        double out = -.5*(data-me[d])*(data-me[d])/var[d] - .5*log(2*PI*var[d]);
        return out;
}

double Gaussian::calpdf_vec(vector<double> data)
{
        double s1,s2,out;
        s1=0; s2=1;
        for (vector<double>::size_type k=0; k!=data.size(); ++k){
                s1 += (data[k]-me[k])*(data[k]-me[k])/var[k];
                s2 *= var[k];
        }

        out = -.5*log(s1)-data.size()/2*log(2*PI)-.5*log(s2);
        return out;
}

double Gaussian::calcdf(double data, unsigned d)
{
        data = (data-me[d])/sqrt(var[d]);

        double L, K, w ;

        double const a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937;
        double const a4 = -1.821255978, a5 = 1.330274429;

        L = fabs(data);
        K = 1.0 / (1.0 + 0.2316419 * L);
        w = 1.0 - 1.0 / sqrt(2 * PI) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K,3) + a4 * pow(K,4) + a5 * pow(K,5));

        if (data < 0 ){
                w= 1.0 - w;
        }
        return w;
}

CGMM::CGMM(string fn, unsigned int spkrNum)
{
	cout<<"Loading GMMs from "<<"\""<<fn<<"\""<<"..."<<endl;
	
	string tline, tmpstr;
	istringstream tline_stream;

	ifstream ifile(fn.c_str());
	if (!ifile.is_open())
	{
		cerr<<"Open file error: "<<fn<<endl;
		return;
	}
	else
	{	
		this->spkrNum = spkrNum;
		allGMM.resize(spkrNum);
		for (vector<GMM>::size_type k=0; k!=allGMM.size(); k++)
		{
			if (k==0){
				//skip the first 3 lines
				getline(ifile,tline); // ~o
				getline(ifile,tline); // <STREAMINFO> %d %d
				getline(ifile,tline); // <VecSize> %d<NULLD><USER><DIAGC>
			}

			getline(ifile,tline); //~h "%s"
			tline_stream.str(tline);
			tline_stream>>tmpstr;
			tline_stream>>tmpstr;
			allGMM[k].name = tmpstr.substr(1,tmpstr.length()-2);
			tline_stream.clear();
									
			getline(ifile,tline); //<BeginHMM>
			getline(ifile,tline); //<NumStates> 3
			getline(ifile,tline); //<STATE> 2

			getline(ifile,tline); //'<NUMMIXES> %d\n'
			tline_stream.str(tline);
			tline_stream>>tmpstr;
			tline_stream>>allGMM[k].nGau;
			tline_stream.clear();
			allGMM[k].gauss.resize(allGMM[k].nGau);
									
			for (vector<Gaussian>::size_type gInd=0; gInd<allGMM[k].gauss.size(); gInd++)
			{				
				getline(ifile,tline); //<MIXTURE> %d %f
				tline_stream.str(tline);
				tline_stream>>tmpstr;
				tline_stream>>tmpstr;
				tline_stream>>allGMM[k].gauss[gInd].weight;
				tline_stream.clear();
				
				getline(ifile,tline); //<MEAN> %d
				tline_stream.str(tline);
				tline_stream>>tmpstr;
				tline_stream>>allGMM[k].gauss[gInd].vecSize;
				allGMM[k].gauss[gInd].me.resize(allGMM[k].gauss[gInd].vecSize);	
				tline_stream.clear();
				getline(ifile,tline); 
				tline_stream.str(tline);
				for (int i=0; i<allGMM[k].gauss[gInd].me.size(); i++)
					tline_stream >> allGMM[k].gauss[gInd].me[i];
				tline_stream.clear();

				getline(ifile,tline); //<VARIANCE> %d
				getline(ifile,tline);
				tline_stream.str(tline);
				allGMM[k].gauss[gInd].var.resize(allGMM[k].gauss[gInd].vecSize);	
				for (int i=0; i<allGMM[k].gauss[gInd].var.size(); i++)
					tline_stream >> allGMM[k].gauss[gInd].var[i];
				tline_stream.clear();

				getline(ifile,tline); //<GCONST> %f
				tline_stream.str(tline);
				tline_stream>>tmpstr; 
				tline_stream>>allGMM[k].gauss[gInd].gConst;
				tline_stream.clear();
			}

			allGMM[k].dim = allGMM[k].gauss[0].vecSize;
			//cout<<"Model "<< allGMM[k].name << " loaded, "<< allGMM[k].nGau <<" Gaussians."<< endl;

			//skip the last 5 lines
			getline(ifile,tline);
			getline(ifile,tline);
			getline(ifile,tline);
			getline(ifile,tline);
			getline(ifile,tline);
		}
	}

	ifile.close();
}
