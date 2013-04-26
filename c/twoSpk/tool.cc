#include "tool.h"
#include <dirent.h>
#include <errno.h>
#include <numeric>
#include <iterator>
using namespace std;

void sum(const vector< vector<double> > &in, int type, vector<double> &out)
{
	if (type==1)	
	{
		for (vector<double>::size_type j=0; j!=in[0].size(); ++j)
                {
			double tmp=0;
			for (vector<double>::size_type i=0; i!=in.size(); ++i)
			{
	                        tmp += in[i][j];
			}
        	        out.push_back(tmp);
                }
	}
	else if (type==2)
	{
		for (vector< vector<double> >::size_type i=0; i!=in.size(); ++i)
		{
			double tmp = accumulate(in[i].begin(),in[i].end(),0);
			out.push_back(tmp);
		}		
	}
}

void operator+=(vector<double> &x, const vector<double> &y)
{
        for (vector<double>::size_type k=0; k<x.size(); ++k)
                x[k] += y[k];
}

vector<double> operator+(const vector<double> &x, const vector<double> &y)
{
        vector<double> z(x);
        z += y;
        return z;
}

void operator*=(vector<double> &x, const double &a)
{
        for (vector<double>::size_type k=0; k<x.size(); ++k)
                x[k] *= a;
}

vector<double> operator*(const vector<double> &x, const double &a)
{
        vector<double> z(x);
        z *= a;
        return z;
}

istream& operator>>(istream &in, vector<string> &out)
{
        string line;
        while(getline(in,line))
                out.push_back(line);
        return in;
}

istream& operator>>(istream &in, vector< vector<double> > &s)
{
	string line;
	while(getline(in,line))
	{
		stringstream ss(line);
		vector<double> vc;
		double tmp;
		while(ss>>tmp)
			vc.push_back(tmp);
		s.push_back(vc);
	}
	return in;	
}

istream& operator>>(istream &in, vector<double> &s)
{
        double val;
	while (in>>val)        
		s.push_back(val);
        return in;
}

bool getDir(string fdFN, string ext, vector<string> &files)
{	
	DIR *dp;
 	struct dirent *dirp;
	if((dp  = opendir(fdFN.c_str())) == NULL) {
		cout << "Error(" << errno << ") opening " << fdFN << endl;
        	return errno;
   	}	

	while ( ((dirp = readdir(dp)) != NULL)) {
		string ent = string(dirp->d_name);
		int pos = ent.find_first_of(".");
		if (ent.substr(pos+1)==ext)
	        	files.push_back(fdFN+"/"+ent);
	}
	closedir(dp);
	return 0;
}

istream& dlmread(istream &in, vector<double> &s, const char delim)
{
	string line;
	while (getline(in,line))
	{
		stringstream ss(line);
		string tmp;
		while (getline(ss,tmp,delim))
			s.push_back(atof(tmp.c_str()));
	}
	return in;
}

ostream& operator<<(ostream &out, const vector<double> &s)
{
        for (vector<double>::const_iterator it=s.begin(); it!=s.end(); ++it)
        	out<<*it<<" "<<flush;
        return out;
}

ostream& operator<<(ostream &out, const vector< vector<double> > &s)
{
        ostream_iterator<double> out_iter(out," ");
        for (vector< vector<double> >::const_iterator it=s.begin(); it!=s.end(); ++it)
        {
                for(vector<double>::const_iterator it0=it->begin(); it0!=it->end(); ++it0)
                        *out_iter++=*it0;
                out<<endl;
        }
	return out;
}

ostream& operator<<(ostream &out, const vector< vector< vector<double> > > &s)
{
	ostream_iterator<double> out_iter(out," ");
	for (vector< vector< vector<double> > >::size_type i=0; i!=s.size(); ++i)
	{
		for (vector< vector<double> >::size_type j=0; j!=s[i].size(); ++j)
	        {	
        	        for(vector<double>::size_type k=0; k!=s[i][j].size(); ++k)
                	        *out_iter++=s[i][j][k];
	                out<<endl;
        	}
	}
	return out;
}

ostream& operator<<(ostream &out, const vector< pair<int,int> > &s)
{
	ostream_iterator<double> out_iter(out," ");
	for (vector< pair<int,int> >::const_iterator it=s.begin(); it!=s.end(); ++it)
        {
		*out_iter++ = it->first;
		*out_iter++ = it->second;
		out<<endl;
        }
	return out;
}

void loadFile(string filename, vector< vector<double> > &data)
{
        ifstream f_str(filename.c_str());
        string line;
        if (f_str.is_open())
        {
                vector<double> row;
                while(getline(f_str,line))
                {
                        istringstream iss(line);
                        double val;
                        while(iss>>val)
				row.push_back(val);
                        data.push_back(row);
                        row.clear();
                }
        }
        else
        {
                cerr<<"Open file "<<filename<<" error!"<<endl;
        }
}

double* readFileAsciiDouble(string filename, unsigned int &nSample)
{
	//cout<<"Reading file "<<filename<<" ";

	double* in; 	
	ifstream ifile(filename.c_str());

	if (!ifile.is_open())
	{
		cerr<<"Open file error: "<<filename<<endl;
		return NULL;
	}
	else
	{				
		double tmp[MAX_SAMPLE];
		nSample = 0;
		while (!ifile.eof()){
			ifile>>tmp[nSample++];			
		}
		nSample--;
		
		in=new double[nSample];		
		for (int k=0; k<nSample; k++){
			in[k] = tmp[k];
		}
		ifile.close();		
//		cout<<nSample<<" samples read."<<endl;
	}

	return in;
}

int* readFileAsciiInt(string filename, unsigned int &nSample)
{
        cout<<"Reading file "<<filename<<" ";

        int* in;
        ifstream ifile(filename.c_str());

        if (!ifile.is_open())
        {
                cerr<<"Open file error: "<<filename<<endl;
                return NULL;
        }
        else
        {
                int tmp[MAX_SAMPLE];
                nSample = 0;
                while (!ifile.eof()){
                        ifile>>tmp[nSample++];
                }
                nSample--;

                in=new int[nSample];
                for (int k=0; k<nSample; k++){
                        in[k] = tmp[k];
                }
                ifile.close();
//                cout<<nSample<<" samples read."<<endl;
        }

        return in;
}

void writeFileAsciiDouble(string fn, double *mat, unsigned int n)
{
	ofstream ofile(fn.c_str());
	if (!ofile.is_open()){
		cerr<<"Open file error: "<<fn<<endl;
	}
	else
	{		
		for (int r=0; r<n; r++)
			ofile<<mat[r]<<" ";
	}
	ofile.close();	
}

void writeFileAsciiVecDouble(string fn, vector<double> mat)
{
	ofstream ofile(fn.c_str());
	if (!ofile.is_open()){
		cerr<<"Open file error: "<<fn<<endl;
	}
	else
	{		
		for (vector<double>::size_type k=0; k!=mat.size(); ++k)
			ofile<<mat[k]<<" ";
	}
	ofile.close();	
}

void writeMatAsciiDouble(string fn, double **mat, unsigned int row, unsigned int col)
{
	ofstream ofile(fn.c_str());
	if (!ofile.is_open()){
		cerr<<"Open file error: "<<fn<<endl;
	}
	else
	{		
		for (int r=0; r<row; r++){
			for (int c=0; c<col; c++)
				ofile<<mat[r][c]<<" ";
			ofile<<endl;
		}
	}
	ofile.close();	
}

void writeMatAsciiVecDouble(string fn, vector< vector<double> > mat)
{
	ofstream ofile(fn.c_str());
	if (!ofile.is_open()){
		cerr<<"Open file error: "<<fn<<endl;
	}
	else
	{		
		for (int r=0; r!=mat.size(); ++r){
			for (int c=0; c!=mat[r].size(); ++c)
				ofile<<mat[r][c]<<" ";
			ofile<<endl;
		}
	}
	ofile.close();	
}

double sinc(double tim, double wn)
{
	if(tim == 0) return(wn);
	else return( sin(wn * PI * tim) / (PI * tim) );
}

double sincDiffOne(double tim, double wn)
{
	if(tim == 0) return(0);
	else
	{
		return( cos(wn * PI * tim) * wn / tim - sin(wn * PI *tim) / (PI * tim * tim) );
	}
}

double sincDiffTwo(double tim, double wn)
{
	if(tim == 0) return( -wn*wn*wn/PI/PI);
	else
	{
		return( -wn*wn*PI*sin(wn*PI*tim)/tim - 2*wn*cos(wn*PI*tim)/tim/tim + 2*sin(wn*PI*tim)/PI/tim/tim/tim );
	}
}

void fft(double *inputR, double *inputI, int N, double direct)
{
	long sigL, i, j, k, n, period, twoPeriod;
	double tmpR, tmpI, uR, uI, wR, wI;

	sigL = long(pow((double)2, N));

	j = 1;
	for(i=1; i<sigL; i++)
	{
		if(i < j)
		{
			tmpR = inputR[j-1];
			tmpI = inputI[j-1];

			inputR[j-1] = inputR[i-1];
			inputI[j-1] = inputI[i-1];

			inputR[i-1] = tmpR;
			inputI[i-1] = tmpI;
		}

		k = sigL/2;
		while (k < j){ j -=  k;	k /= 2;	}
		j += k;
	}

	for(n=1; n<=N; n++ )
    {  
		twoPeriod = long(pow(2.0, n));
        period = twoPeriod/2;
        uR = 1.0; 
        uI = 0.0; 
        wR = double( cos( PI/period ) ); 
        wI = double( -1.0 * sin( PI/period * direct) );

        for(j=0; j<period; j++ ) 
        {  
			for(i=j; i<sigL; i+=twoPeriod)
			{
				tmpR = inputR[i+period]*uR - inputI[i+period]*uI;
                tmpI = inputR[i+period]*uI + inputI[i+period]*uR;
				
				inputR[i+period] = inputR[i] - tmpR; 
				inputI[i+period] = inputI[i] - tmpI; 
				inputR[i] += tmpR ; 
				inputI[i] += tmpI; 
			}
			tmpR = uR*wR - uI*wI; 
			tmpI = uR*wI + uI*wR; 
			uR = tmpR; 
			uI = tmpI; 
		} 
	} 
}

int maxPos(double *input, int d1, int d2)
{
	int judge=d1;
	double mv=input[d1];
	for(int step=d1+1; step<=d2; step++)
	{	
		if (input[step]>mv)
		{
			judge=step;
			mv=input[step];
		}
	}

	return(judge);
}

int minPos(double *input, int d1, int d2)
{
	int judge=d1;
	double mv=input[d1];
	for(int step=d1+1; step<=d2; step++)
	{	
		if (input[step]<mv)
		{
			judge=step;
			mv=input[step];
		}
	}

	return(judge);
}

double bessi0(double x)
{
	double ax,ans;
	double y;

	ax = fabs(x);
	if (ax < 3.75)
	{
		y = x/3.75;
		y *= y;
		ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y*0.45813e-2)))));
	}
	else
	{
		y = 3.75 / ax;
		ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
	}
	return ans;
}

double bessi1(double x)
{
	double ax, y, ans;
	if (fabs(x)<3.75)
	{
		y = x/3.75;
		y *= y;
		ans = x * ( 0.5 + y*( 0.87890594 + y*( 0.51498869 + y*( 0.15084934 + y*( 0.02658733 + y* (0.00301532 + y*0.00032411 ) ) ) ) ) );
	}
	else
	{	
		ax = fabs(x);
		y = 3.75/ax;
		ans = (exp(ax)/sqrt(ax)) * (0.39894228 + y*( -0.03988024 + y*( -0.00362018 + y*( 0.00163801 + y*( -0.01031555 + y*( 0.02282967 + y*( -0.02895312 + y*( 0.01787654 + y*(-0.00420059) ) ) ) ) ) ) ) );

		if (x<0) ans *= -1;
	}
   
    return(ans);
}

double zeroCross(double *input, int sLen)
{
	int nCross=-1;
	double currentSign=1, sp, ep;
	double zC=sLen*2;

	if (input[0]<0) currentSign=-1;
	
	for(int n=1; n<sLen; n++)
	{
		if ((input[n]*currentSign)<0)
		{
			nCross++;
			if (nCross==0) sp=n;
			else ep=n;

			currentSign *= -1;

			if ( ((nCross%2)==0) && (nCross>0))
				zC=(ep-sp)/double(nCross);
		}
	}

	if(nCross==1) zC=ep*2.0/3.0;

	return(zC);
}

double zeroCross(vector<double> &input, int sLen)
{
        int nCross=-1;
        double currentSign=1, sp, ep;
        double zC=sLen*2;

        if (input[0]<0) currentSign=-1;

        for(int n=1; n<sLen; n++)
        {
                if ((input[n]*currentSign)<0)
                {
                        nCross++;
                        if (nCross==0) sp=n;
                        else ep=n;

                        currentSign *= -1;

                        if ( ((nCross%2)==0) && (nCross>0))
                                zC=(ep-sp)/double(nCross);
                }
        }

        if(nCross==1) zC=ep*2.0/3.0;

        return(zC);
}

void processBar(double cur, double tot)
{
	//backspace
	if (cur!=tot)
	{
		printf("%04.1f%%", 100*cur/tot);
		cout << "\b\b\b\b\b" << flush;
	}
	else
		cout<<"Done!"<<endl;
}

unsigned int* readUnsignedInt(string filename, unsigned int &nSample)
{
	cout<<"Reading file "<<filename<<" ";

	unsigned int *in;
	ifstream ifile(filename.c_str());

	if (!ifile.is_open())
	{
		cerr<<"Open file error: "<<filename<<endl;
		return NULL;
	}
	else
	{				
		unsigned int tmp[MAX_SAMPLE];
		nSample = 0;
		while (!ifile.eof()){
			ifile>>tmp[nSample++];			
		}
		
		in=new unsigned int[nSample];		
		for (int k=0; k<nSample; k++){
			in[k] = tmp[k];
		}
		ifile.close();		
		cout<<nSample<<" samples read."<<endl;
	}

	return in;
}

void normWavSig(double* sig, unsigned int length, double dB)
{	
	double eng=0, scale;
	for (unsigned int i=0; i<length; i++)
		eng += sig[i]*sig[i];	
	eng /=length;
	
	scale = pow(10.0, (dB-10*log10(eng))/20 );
    
	eng = 0;
	for (unsigned int i=0; i<length; i++)
	{
		sig[i] = sig[i]*scale; // 16-bit wav signals			
		eng +=sig[i]*sig[i];
	}
	eng /= length;
	
//	cout<<"Normalized to "<<10*log10(eng)<<" dB average intensity."<<endl;
}

double normpdf(double data, double me, double var)
{	
	double out = -.5*(data-me)*(data-me)/var - .5*log(2*PI*var);
	return out;
}

double normcdf(double data, double me, double var)
{
	data = (data-me)/sqrt(var);

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
