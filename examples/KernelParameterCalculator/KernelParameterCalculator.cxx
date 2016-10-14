/*
 * KernelParameterCalculator.cxx
 *
 *  Created on: 26 Aug 2015
 *      Author: lewish
 */




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
#include "TH1D.h"
#include "TGraph.h"
#include <iostream>
#include <sstream>
#include "RooGenericPdf.h"

#include <chrono>
#include <ctime>
#include <ratio>
#include <cmath>


#include "stkImageHistogram.h"//write to histogram
#include "stkRawImageIO.h"//write pedestal in the frames
#include "stkRawImageStackIO.h"//load in the frames
#include "stkImageStack.h"//hold the files to be loaded
#include "stkImage.h"//will hold the result
#include "stkImageSum.h"//used to sum the stack
#include "stkImageDivision.h"//used to divide through to get final result
#include "stkImageVariance.h"//used for variance
#include "stkImageMask.h"
#include "stkImageBadPixelAverage.h"
#include "stkImageMinus.h"
#include "stkImageResize.h"
//#include "stkImageFlatFieldCorrection.h"
#include "stkITKAdapter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "stkITKTools.h"
#include "itkImageFileWriter.h"
//#include "stkDarkFrameDetection.h"
#include "tbb/tbb.h"
#include "itkConvolutionImageFilter.h"


int main(int argc, char*argv[])
{
	TFile* file=new TFile("example.root", "RECREATE");
	TGraph graph(10);

	graph.SetPoint(0, 0.5, 1.005);
	graph.SetPoint(1, 1, 1.010);
	graph.SetPoint(2, 1.5, 1.015);
	graph.SetPoint(3, 2, 1.021);
	graph.SetPoint(4, 2.5, 1.026);
	graph.SetPoint(5, 3, 1.031);
	graph.SetPoint(6, 3.5, 1.037);
	graph.SetPoint(7, 4, 1.042);
	graph.SetPoint(8, 4.5, 1.044);

	graph.Draw("*");

	RooRealVar x("x", "x", 0, 20);//Define Roo variable
	RooRealVar y("y", "y", 0.9, 1.1);//Define Roo variable

	TH1F* graphHist = graph.GetHistogram();
	graphHist->Draw();
	graphHist->Write();
	graph.Draw();
	graph.Write();
	//RooDataHist  dataGraph("No Bu Hist", "No Bu Hist",x, graphHist);
	RooPlot* frame = x.frame( );
	//dataGraph.plotOn(frame);



	TH1F* signal_histBu = new TH1F("OAR", "OAR", 38, 1.5, 20);
	signal_histBu->SetBinContent(1,100);
	signal_histBu->SetBinContent(2,98.8);
	signal_histBu->SetBinContent(3,96.6);
	signal_histBu->SetBinContent(4,94.4);
	signal_histBu->SetBinContent(5,92.0);
	signal_histBu->SetBinContent(6, 89.6);
	signal_histBu->SetBinContent(7, 87.2);
	signal_histBu->SetBinContent(8, 84.8);
	signal_histBu->SetBinContent(9, 82.5);
	signal_histBu->SetBinContent(10, 80.2);
	signal_histBu->SetBinContent(11,77.95);
	signal_histBu->SetBinContent(12, 75.7);
	signal_histBu->SetBinContent(13, 73.65);
	signal_histBu->SetBinContent(14, 71.6);
	signal_histBu->SetBinContent(15, 69.5);
	signal_histBu->SetBinContent(16, 67.4);
	signal_histBu->SetBinContent(17, 65.5);
	signal_histBu->SetBinContent(18, 63.6);
	signal_histBu->SetBinContent(19, 61.75);
	signal_histBu->SetBinContent(20, 59.9);
	signal_histBu->SetBinContent(21, 58.15);
	signal_histBu->SetBinContent(22, 56.4);
	signal_histBu->SetBinContent(23, 54.75);
	signal_histBu->SetBinContent(24, 53.1);
	signal_histBu->SetBinContent(25, 51.6);
	signal_histBu->SetBinContent(26, 50.1);
	signal_histBu->SetBinContent(27, 48.65);
	signal_histBu->SetBinContent(28, 47.2);
	signal_histBu->SetBinContent(29,45.85);
	signal_histBu->SetBinContent(30, 44.5);
	signal_histBu->SetBinContent(31, 43.25);
	signal_histBu->SetBinContent(32, 42);
	signal_histBu->SetBinContent(33, 40.8);
	signal_histBu->SetBinContent(34, 39.6);
	signal_histBu->SetBinContent(35,38.45);
	signal_histBu->SetBinContent(36, 37.3);
	signal_histBu->SetBinContent(37,36.2);
	signal_histBu->SetBinContent(38,  35.1);

	signal_histBu->Draw();
	signal_histBu->Write();

	RooRealVar alpha("alpha", "alpha",2.8, 0, 10);//Define Roo variable
	RooRealVar beta("beta", "beta",9.08 ,0, 10);//Define Roo variable
	RooRealVar omega("omega", "omega",1.335 , 0, 10);//Define Roo variable
	RooRealVar a("a", "a", 0.08 ,0, 10);//Define Roo variable
	RooRealVar F("F", "F", 0.704, 0, 20);//Define Roo variable
	RooRealVar A("A", "A", 100, 0, 200);//Define Roo variable
	RooRealVar mu("mu", "F", 0.0404, 0.02, 20);//Define Roo variable
	RooRealVar nu("nu", "F", 0.002, 0.0002, 20);//Define Roo variable
	RooDataHist  dataGraph("No Bu Hist", "No Bu Hist",x, signal_histBu);

	RooGenericPdf kernelFunction("kernelFunction","kernelFunction","A*exp(-mu*x*(1-nu*x))",RooArgSet(x,A, mu, nu));//Defining top hat PDF]


	//(1/(2*3.14159265*.0001))*((1-exp(-alpha*x))*(alpha*(1-F)*exp(-alpha*.0001)+beta*F*exp(-beta*.0001))+((a*pow(x,2))/pow(omega*.0001+x,2)))
	kernelFunction.fitTo(dataGraph);


	dataGraph.plotOn(frame);
	kernelFunction.plotOn(frame);
	frame->Draw();
	frame->Write();
	file->Close();

	return 1;

}
