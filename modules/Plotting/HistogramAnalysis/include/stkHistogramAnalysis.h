/*
 * stkHistogramAnalysis.h
 *
 *  Created on: 14 Aug 2015
 *      Author: lewish
 */

#ifndef __stkHistogramAnalysis_h
#define __stkHistogramAnalysis_h

#include <string>
#include <iostream>
#include <memory>
#include <cmath>
#include "stkImage.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TFile.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "tbb/tbb.h"

namespace stk{


class HistogramAnalysis{

public:
	HistogramAnalysis();
	virtual ~HistogramAnalysis();


	void LoadRootFile(std::string filename);
	void GenerateXSlice();
	void TwoDFit(std::shared_ptr< stk::Image< float > >buImage, std::shared_ptr< stk::Image< float > >noBuImage);
	void TwoDFit(std::shared_ptr< stk::Image< float > >image);
	void MaxMin();
	void OneDFit(std::shared_ptr< stk::Image< float > >image);
	std::shared_ptr< stk::Image<float> > kFunc();
	void ErrorFunction();
	void PeakValueFind();
	void SlopeAnalysis(int x, int y);
	void OutsideBeam();




private:

	double r();
	TH2 *inputHistogram;
	TFile *inputFile;
	struct kernelParameters
	{
		double x, y, xs, ys, z,alpha, beta, a, omega, F, Az, mu, nu, count;

	};

	kernelParameters m_kernelParam;

};


}

#include "stkHistogramAnalysis.hxx"

#endif /* __stkImageHistogramAnalysis_h */
