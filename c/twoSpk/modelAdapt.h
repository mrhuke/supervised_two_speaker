#ifndef MODEL_ADAPT_H
#define MODEL_ADAPT_H

#include "gmm.h"
#include "tool.h"

class modelAdapt{
	double isnr;  // estimated input SNR
public:
	modelAdapt(): isnr(0) {}
	modelAdapt(double s): isnr(s) {}

	void adapt(double*, size_t, GMM&, GMM&);  //normlize mixture to 60-dB, and modify two models
	void adapt(double*, size_t, GMM&);  // adapt mixture and only the interfering model
};

#endif
