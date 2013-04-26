#ifndef TOOL_H
#define TOOL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

#define PI (3.1415926535897932384626433832795)
#define MAX_SAMPLE 80000					// 5s-long utterance at 16kHz

template<typename T> struct vector2D {
     typedef vector< vector<T> > type;
};
template<typename T> struct vector3D {
     typedef vector< vector< vector<T> > > type;
};
template<typename T> struct vector4D {
     typedef vector< vector< vector< vector<T> > > > type;
};
template<typename T> struct vector5D {
     typedef vector< vector< vector< vector< vector<T> > > > > type;
};

void processBar(double cur, double tot);
void normWavSig(double* sig, unsigned int length, double dB);  // normalize the signal to have dB (dB) energy on average

ostream& operator<<(ostream &, const vector<double> &);
ostream& operator<<(ostream &, const vector< vector<double> > &);
ostream& operator<<(ostream &, const vector< vector< vector<double> > > &);
ostream& operator<<(ostream &, const vector< pair<int,int> > &);

istream& operator>>(istream &, vector< vector<double> > &);
istream& operator>>(istream &, vector<double> &);
istream& operator>>(istream &, vector<string> &);
istream& dlmread(istream &, vector<double> &, char);
bool getDir(string,string,vector<string> &);

vector<double> operator+(const vector<double> &, const vector<double> &);
void operator+=(vector<double> &, const vector<double> &);
vector<double> operator*(const vector<double> &, const double &);
void operator*=(vector<double> &, const double &);

void sum(const vector< vector<double> > &, int, vector<double> &);

void loadFile(string, vector< vector<double> > &);
double* readFileAsciiDouble(string filename, unsigned int &nSample);
int* readFileAsciiInt(string filename, unsigned int &nSample);
void writeFileAsciiDouble(string fn, double *mat, unsigned int n);
void writeFileAsciiVecDouble(string fn, vector<double> mat);
void writeMatAsciiDouble(string fn, double **mat, unsigned int row, unsigned int col);
void writeMatAsciiVecDouble(string fn, vector< vector<double> > mat);
unsigned int* readUnsignedInt(string filename, unsigned int &nSample);

inline double sigmoid(double x){ return( 2/(1+exp(-2*x)) -1 ); }
inline double hamming(double x){ return( 0.54 + 0.46*cos(x*PI) ); }
inline double hammingDiffOne(double x){ return( -0.46*PI*sin(x*PI) ); }
inline double hammingDiffTwo(double x){ return( -0.46*PI*PI*cos(x*PI) ); }

double sinc(double tim, double wn);
double sincDiffOne(double tim, double wn);
double sincDiffTwo(double tim, double wn);
double bessi0(double x);
double bessi1(double x);

inline bool isNotZero(double x) { return x!=0; }
int maxPos(double *input, int d1, int d2);
int minPos(double *input, int d1, int d2);
double normpdf(double data, double me, double var);
double normcdf(double data, double me, double var);

void fft(double *inputR, double *inputI, int N, double direct);
double zeroCross(double *input, int sLen);
double zeroCross(vector<double> &, int sLen);



#endif
