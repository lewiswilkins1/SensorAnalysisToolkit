/*
 * stkHistogramAnalysis.hxx
 *
 *  Created on: 14 Aug 2015
 *      Author: lewish
 */

#ifndef __stkHistogramAnalysis_hxx
#define __stkHistogramAnalysis_hxx



#include "stkHistogramAnalysis.h"
#include "RooRealVar.h"
#include "RooGlobalFunc.h"
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
#include "RooChebychev.h"
#include "RooHistPdf.h"
#include "RooGenericPdf.h"
#include "RooMinuit.h"
#include "RooExtendPdf.h"
#include <cmath>
#include "tbb/tbb.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooExponential.h"
#include "RooEffProd.h"
#include "RooFormulaVar.h"
#include <memory>
#include "TMath.h"
#include "TH1.h"


using namespace RooFit;

Double_t LogisticFunction(Double_t *x, Double_t *par) {
	//***********************************************
	Double_t fitval;
	

	Double_t max = par[0];
	Double_t steepness = par[1];
	Double_t midpoint = par[2];
	Double_t xval = x[0];


	


	if (xval>650&&xval<940)
	{
		fitval = max / (1 + std::exp(-steepness*(xval-midpoint))) + par[3];
		return fitval;
	}
	else
	{
		TF1::RejectPoint();
		return 0;
	}	
}
Double_t ErfGauss(Double_t *x, Double_t *par) {
	//***********************************************
	Double_t fitval;
	Double_t var = (x[0]-par[2])/par[3];
	Double_t varc = (x[0]-par[4])/par[5];

	Double_t mean = par[0];
	Double_t sigma = par[1];
	Double_t xval = x[0];


	if(xval<343){


		fitval = par[6]*TMath::Gaus(xval, mean, sigma)+par[7]*(TMath::Erf(var)+1);
	}

	else if(xval>660){

		fitval =par[6]* TMath::Gaus(xval, mean, sigma)+par[8]*TMath::Erfc(varc);
	}


	return fitval;
}

Double_t Gauss(Double_t *x, Double_t *par) {
	//***********************************************
	Double_t fitval;
	Double_t xval = x[0];
	Double_t mean = par[0];
	Double_t sigma = par[1];
	fitval = par[2]*TMath::Gaus(xval,mean,sigma);

	return fitval;
}

Bool_t reject;
Double_t GaussRangeBig(Double_t *x, Double_t *par) {
	//***********************************************
	Double_t fitval;
	Double_t xval = x[0];
	Double_t mean = par[0];
	Double_t sigma = par[1];

	if (reject && x[0] > 100 && x[0] < 920) {
		TF1::RejectPoint();
		return 0;
	}



	fitval = par[2]*TMath::Gaus(xval,mean,sigma);


	return fitval;
}
Double_t GaussRangeSmall(Double_t *x, Double_t *par) {
	//***********************************************
	Double_t fitval;
	Double_t xval = x[0];
	Double_t mean = par[0];
	Double_t sigma = par[1];

	if ((reject && x[0] > 290 && x[0] < 730 )||( reject &&x[0]<100 && x[0]>920)) {
		TF1::RejectPoint();
		return 0;
	}



	fitval = par[2]*TMath::Gaus(xval,mean,sigma);


	return fitval;
}

Double_t Erfunc(Double_t *x, Double_t *par) {
	//***********************************************
	Double_t fitval;
	Double_t var = (x[0]-par[1])/par[2];



	Double_t xval = x[0];





	fitval = par[0]*(TMath::Erf(var)+1)+par[3];



	return fitval;
}


Double_t Erfn(Double_t *x, Double_t *par) {
	//***********************************************
	Double_t fitval;
	Double_t var = (x[0]-par[1])/par[2];
	Double_t varc = (x[0]-par[4])/par[5];


	Double_t xval = x[0];


	if(xval>60&&xval<180){


		fitval = par[0]*(TMath::Erf(var)+1)+par[6];
	}

	else if(xval>750&&xval<830){

		fitval =par[3]*TMath::Erfc(varc)+par[7];
	}
	else
	{
		TF1::RejectPoint();
	}


	return fitval;
}




Double_t expFunc(Double_t *x, Double_t *par) {
	//***********************************************
	/*3 170-270, 610-710
	 *4 40-140, 730-830
	 *6 25-75, 855-955
	 */
	Double_t fitval;
	Double_t var = (x[0]/par[0]);

	if(x[0]>=35 && x[0] <=145 )
	{
		fitval = par[1] * TMath::Exp(var);
	}
	// else if(x[0] > 730 && x[0] < 830)
	// {
	// 	fitval = par[2] - par[3]*TMath::Exp(var);
	// }

	else
	{
		TF1::RejectPoint();
		fitval = 0;
	}


	return fitval;
}

namespace stk{


HistogramAnalysis::HistogramAnalysis(){


	m_kernelParam.Az = 111.56;
	m_kernelParam.F = 0.704;
	m_kernelParam.a = 0.0809;
	m_kernelParam.alpha = 2.8;
	m_kernelParam.beta = 9.08;
	m_kernelParam.count = 0;
	m_kernelParam.mu = 5.91158e-02;
	m_kernelParam.nu = 1.32742e-08;
	m_kernelParam.omega = 1.335;
	m_kernelParam.x = 0;
	m_kernelParam.y = 0;
	m_kernelParam.xs = 0;
	m_kernelParam.ys = 0;
	m_kernelParam.z = 5;


}

HistogramAnalysis::~HistogramAnalysis(){}


void HistogramAnalysis::LoadRootFile(std::string filename){
	filename+=".root";
	const char * charName = filename.c_str();
	inputFile = new TFile(charName, "UPDATE");
	std::cout<<"ROOT files loaded."<<std::endl;



}

void HistogramAnalysis::GenerateXSlice(){




	TH2* noBuildUp = (TH2*)inputFile->Get("NoBu");//Load in the two histogram
	TH2* buildUp = (TH2*)inputFile->Get("Bu");


	RooRealVar x("x", "x", -0.5, 1023.5);//Define Roo variable




	TH1* inputNoBu = noBuildUp->ProjectionX("Y Slice NoBu", 512, 512);//Taking projection of first histogram
	TH1D* sliceNoBu = new TH1D("Y Slice w NoBu", "Y Slice", 1024, -0.5,1023.5);//Putting into a new TH1 histogram. It removes the error bar, not pointless
	for(int ibin=1; ibin <= inputNoBu->GetNbinsX(); ibin++ ){
		double value = inputNoBu->GetBinContent(ibin);
		sliceNoBu->SetBinContent(ibin,value);

	}




	TH1* inputBu = buildUp->ProjectionX("Y Slice Bu", 512, 512);//Taking projection of second histogram
	TH1D* sliceBu = new TH1D("Y Slice w Bu", "Y Slice", 1024, -0.5,1023.5);//Same as above
	for(int ibin=1; ibin <= inputBu->GetNbinsX(); ibin++ ){
		double value = inputBu->GetBinContent(ibin);
		sliceBu->SetBinContent(ibin,value);

	}


	RooDataHist  SliceRooNoBu("Y slice", "Y Slice",RooArgList(x), sliceNoBu);//Creating two roo hists for slices
	RooDataHist  SliceRooBu("Y slice", "Y Slice",RooArgList(x), sliceBu);

	RooPlot* frame = x.frame(Title("5x5 NoBU 200 MU") );//Creating two frames to plot slices on
	RooPlot* frame2 = x.frame(Title("5x5 BU 200 MU") );
	SliceRooNoBu.plotOn( frame );//Plotting slices on frames
	SliceRooBu.plotOn( frame2 );




	RooRealVar x_mean("x_Mean","mean of gaussian",512,0,1024);//Setting up gaussian parameters
	RooRealVar x_sigma("x_Sigma","width of gaussian",256, 0, 1024);
	RooGaussian gx("gx","gx",x,x_mean,x_sigma);//Setting up gaussian


	RooRealVar sharpness("sharpness","sharpness",50, 5, 100);//Setting up tophat parameters
	RooRealVar topwidth("topwidth","topwidth", 150, 100,250);
	RooRealVar offset("offset","offset", 512, 200,600);
	RooRealVar up("up","up", 50, 0,200);
	RooRealVar fsig("fsig","signal fraction", 0.5,0,1);//Variabel for PDF fraction
	RooGenericPdf g("g","g","sqrt((1+pow(sharpness,2))/(1+pow(sharpness,2)*pow(cos((x-offset)/topwidth),2)))*cos((x-offset)/topwidth)+1",RooArgSet(x,sharpness,topwidth,offset));//Defining top hat PDF

	//RooRealVar fsig("fsig","signal fraction", 0.5,0,1);//Variabel for PDF fraction
	RooAddPdf model("model", "model", RooArgList(gx,g), fsig);//Combines the gaussian and tophat PDF
	x.setRange("range", 100,950);//Sets range for buildup option

	model.fitTo(SliceRooNoBu);//Fits combined model to first hist
	model.plotOn(frame);//Plots fit

	TH1* bkg_histNoBu = gx.createHistogram("Background HistNo Bu", x, Binning(1024));//Creates histogram from gaussian bkg fit
	RooPlot* frame3 = x.frame(Title("PDF Components"));//Creates new frame

	model.plotOn(frame3, Components(g),LineColor(kRed));//Plots the two components of the combined function on a seperate fram
	model.plotOn(frame3, Components(gx), LineColor(kBlue));



	model.fitTo(SliceRooBu, Range("range"));//Fits to other histogram but over range to exclude build up anomaly at edge
	model.plotOn(frame2);

	TH1* bkg_histBu = gx.createHistogram("Background Hist", x, Binning(1024));//Generates a histogrm for gaussian bkg fit


	TH1D* signal_histBu = new TH1D("Signal Hist Bu", "Signal Hist", 1024, -0.5, 1023.5);//creating two new histograms for signal plots
	TH1D* signal_histNoBu = new TH1D("Signal Hist No Bu", "Signal Hist", 1024, -0.5, 1023.5);


	double integralNoBu(0);//Integrating the slice
	for(int ibin=1; ibin <= bkg_histNoBu->GetNbinsX(); ibin++ ){
		integralNoBu += sliceNoBu->GetBinContent(ibin);
	}


	double maximum = sliceNoBu->GetMaximum();//Finds maximum
	for(int ibin=1; ibin <= bkg_histNoBu->GetNbinsX(); ibin++ ){//Removes the background from the signal
		double bkg_value = 6.51667e-01*integralNoBu*bkg_histNoBu->GetBinContent(ibin);
		double signalAndBkgValue = sliceNoBu->GetBinContent(ibin);
		double signal = signalAndBkgValue - bkg_value;
		signal_histNoBu->SetBinContent(ibin,signal);
	}


	double sigIntegralNoBu(0);//Integrates signal
	for(int ibin=1; ibin <= bkg_histNoBu->GetNbinsX(); ibin++ ){
		sigIntegralNoBu+= signal_histNoBu->GetBinContent(ibin);
	}

	for(int ibin=1; ibin <= bkg_histNoBu->GetNbinsX(); ibin++ ){//Normalises the signal
		double normsig  = signal_histNoBu->GetBinContent(ibin)/sigIntegralNoBu;
		signal_histNoBu->SetBinContent(ibin, normsig);
	}



	double integralBu(0);//Same as above for other histogram
	double maximum2 = sliceBu->GetMaximum();
	for(int ibin=1; ibin <= bkg_histBu->GetNbinsX(); ibin++ ){
		integralBu += sliceBu->GetBinContent(ibin);

	}

	for(int ibin=1; ibin <= bkg_histBu->GetNbinsX(); ibin++ ){
		double bkg_value = 2.49357e-01*integralBu*bkg_histBu->GetBinContent(ibin);
		double signalAndBkgValue = sliceBu->GetBinContent(ibin);
		double signal = signalAndBkgValue - bkg_value;
		signal_histBu->SetBinContent(ibin,signal);
	}



	double sigIntegralBu(0);
	for(int ibin=1; ibin <= bkg_histBu->GetNbinsX(); ibin++ ){
		sigIntegralBu+= signal_histBu->GetBinContent(ibin);
	}
	for(int ibin=1; ibin <= bkg_histBu->GetNbinsX(); ibin++ ){
		double normsig  = signal_histBu->GetBinContent(ibin)/sigIntegralBu;
		signal_histBu->SetBinContent(ibin, normsig);
	}




	//Drawing stuff
	TCanvas* signals = new TCanvas("signals", "signals", 1200,800);
	signals->cd(1);

	signal_histNoBu->Draw("hist");
	signal_histBu->Draw("hist same");
	signal_histNoBu->Write();
	signal_histBu->Write();
	bkg_histBu->Draw();
	bkg_histBu->Write();
	bkg_histNoBu->Draw();
	bkg_histNoBu->Write();


	TCanvas* c = new TCanvas("plot", "plot", 1200,800);
	c->cd(1);

	frame->Draw();
	frame->Write();
	c->cd(2);
	frame2->Draw();

	frame2->Write();

	frame3->Draw();

	frame3->Write();

	inputFile->Close();
}

void HistogramAnalysis::TwoDFit(std::shared_ptr< stk::Image< float > >buImage, std::shared_ptr< stk::Image< float > >noBuImage){




	RooRealVar x("x", "x", -0.5, 1023.5);//Define Roo variable
	RooRealVar y("y", "y", -0.5, 1023.5);//Define Roo variable
	TH2* noBuildUp = (TH2*)inputFile->Get("NoBu");//Load in the two histogram
	TH2* buildUp = (TH2*)inputFile->Get("Bu");

	RooDataHist  SliceRooNoBu("No Bu Hist", "No Bu Hist",RooArgList(x,y), noBuildUp);//Creating two roo hists for slices
	RooDataHist  SliceRooBu("BU Hist", "BU Hist",RooArgList(x,y), buildUp);

	RooRealVar x_mean("x_Mean","mean of gaussian x",512,0,1024);//Setting up gaussian parameters
	RooRealVar x_sigma("x_Sigma","width of gaussian x",256, 100, 600);
	RooGaussian gx("gx","gx",x,x_mean,x_sigma);//Setting up gaussian

	RooRealVar y_mean("y_Mean","mean of gaussian y",512,0,1024);//Setting up gaussian parameters
	RooRealVar y_sigma("y_Sigma","width of gaussian y",256, 100, 600);
	RooGaussian gy("gy","gy",y,y_mean,y_sigma);//Setting up gaussian

	std::cout<<"PDFs Created"<<std::endl;
	RooRealVar sharpness_x("sharpness_x","sharpness",50, 5, 100);//Setting up tophat parameters
	RooRealVar topwidth_x("topwidth_x","topwidth", 150, 100,250);
	RooRealVar offset_x("offset_x","offset", 512, 200,600);
	//RooRealVar up("up","up", 50, 0,200);
	RooRealVar fsig_x("fsig_x","signal fraction", 0.5,0,1);//Variabel for PDF fraction
	RooGenericPdf tophat_x("tophat_x","tophat_x","sqrt((1+pow(sharpness_x,2))/(1+pow(sharpness_x,2)*pow(cos((x-offset_x)/topwidth_x),2)))*cos((x-offset_x)/topwidth_x)+1",RooArgSet(x,sharpness_x,topwidth_x,offset_x));//Defining top hat PDF


	RooRealVar sharpness_y("sharpness_y","sharpness",50, 5, 100);//Setting up tophat parameters
	RooRealVar topwidth_y("topwidth_y","topwidth", 150, 100,250);
	RooRealVar offset_y("offset_y","offset", 512, 200,600);
	//RooRealVar up("up","up", 50, 0,200);
	RooRealVar fsig_y("fsig_y","signal fraction", 0.5,0,1);//Variabel for PDF fraction
	RooGenericPdf tophat_y("tophat_y","tophat_y","sqrt((1+pow(sharpness_y,2))/(1+pow(sharpness_y,2)*pow(cos((y-offset_y)/topwidth_y),2)))*cos((y-offset_y)/topwidth_y)+1",RooArgSet(y,sharpness_y,topwidth_y,offset_y));//Defining top hat PDF

	std::cout<<"Tophats done"<<std::endl;


	RooAddPdf model_x("model_x", "model", RooArgList(gx,tophat_x), fsig_x);//Combines the gaussian and tophat PDF
	RooAddPdf model_y("model_y", "model", RooArgList(gy,tophat_y), fsig_y);//Combines the gaussian and tophat PDF

	RooProdPdf fullModel("fullModel","bkgModel",model_x,model_y);



	RooProdPdf backgroundModel("bkgModel","bkgModel",gx,gy);
	//	RooProdPdf signalModel("signalModel", "signalModel",tophat_x, tophat_y);
	//


	x.setRange("LFT",100.5,210.5);
	y.setRange("LFT",100.5,950.5);

	x.setRange("RGT",800.5,950.5);
	y.setRange("RGT",100.5,950.5);

	x.setRange("T",210.5,800.5);
	y.setRange("T",800.5,950.5);

	x.setRange("BM",210.5,800.5);
	y.setRange("BM",100.5,250.5);

	x.setRange("LFT_BU",-0.5,150.5);
	y.setRange("LFT_BU",-0.5,1023.5);

	x.setRange("RGT_BU",800.5,1023.5);
	y.setRange("RGT_BU",-0.5,1023.5);

	x.setRange("T_BU",150.5,800.5);
	y.setRange("T_BU",800.5,1023.5);

	x.setRange("BM_BU",150.5,800.5);
	y.setRange("BM_BU",-0.5,150);




	std::cout<<"Models combined"<<std::endl;

	backgroundModel.fitTo(SliceRooNoBu,Range("LFT_BU, RGT_BU, T_BU, BM_BU"));
	TH1* bkg_histNoBu = backgroundModel.createHistogram("Background Hist No Bu",x,Binning(1024),YVar(y,Binning(1024)));

	backgroundModel.fitTo(SliceRooBu,Range("LFT, RGT, T, BM"));
	TH1* bkg_histBu = backgroundModel.createHistogram("Background Hist Bu",x,Binning(1024),YVar(y,Binning(1024)));
	//
	//

	TH2D* signal_histBu = new TH2D("Signal Hist Bu", "Signal Hist", 1024, -0.5, 1023.5, 1024, -0.5, 1023.5);//creating two new histograms for signal plots
	TH2D* signal_histNoBu = new TH2D("Signal Hist No Bu", "Signal Hist", 1024, -0.5, 1023.5, 1024, -0.5, 1023.5);

	double integralNoBu(0);//Integrating the slice
	for(int ibin=1; ibin <= 1024*1024; ibin++ ){
		integralNoBu += noBuildUp->GetBinContent(ibin);
	}
	std::cout<<bkg_histNoBu->GetNcells()<<std::endl;


	double n(0);
	double sideband(0);
	for( int ibin=1; ibin<= 1024*1024; ibin++){
		if(ibin<150+n*1024&&ibin>n*1024){
			sideband+= noBuildUp->GetBinContent(ibin);

		}
		else if(ibin>800+n*1024&&ibin<(n+1)*1024){
			sideband+= noBuildUp->GetBinContent(ibin);

		}
		else if(ibin>150+n*1024&&ibin<800+n*1024&&n<150)
		{
			sideband+= noBuildUp->GetBinContent(ibin);
		}
		else if(ibin>150+n*1024&&ibin<800+n*1024&&n>=800)
		{
			sideband+= noBuildUp->GetBinContent(ibin);
		}
		if(ibin%1024==0){
			n++;
		}

	}



	double maximum = noBuildUp->GetMaximum();//Finds maximum
	for(int ibin=1; ibin <= 1024*1024; ibin++ ){//Removes the background from the signal
		double bkg_value = sideband*bkg_histNoBu->GetBinContent(ibin);
		double signalAndBkgValue = noBuildUp->GetBinContent(ibin);
		double signal = signalAndBkgValue - bkg_value;
		signal_histNoBu->SetBinContent(ibin,signal);
	}


	double sigIntegralNoBu(0);//Integrates signal
	for(int ibin=1; ibin <= 1024*1024; ibin++ ){
		sigIntegralNoBu+= signal_histNoBu->GetBinContent(ibin);
	}

	for(int ibin=1; ibin <= 1024*1024; ibin++ ){//Normalises the signal
		double normsig  = signal_histNoBu->GetBinContent(ibin)/sigIntegralNoBu;
		signal_histNoBu->SetBinContent(ibin, normsig);
	}



	double integralBu(0);//Integrating the slice
	for(int ibin=1; ibin <= 1024*1024; ibin++ ){
		integralNoBu += buildUp->GetBinContent(ibin);
	}
	std::cout<<bkg_histNoBu->GetNcells()<<std::endl;


	n=0;
	sideband=0;



	for( int ibin=1; ibin<= 1024*1024; ibin++){
		if(ibin<210+n*1024&&ibin>n*1024+100){
			sideband+= buildUp->GetBinContent(ibin);

		}
		else if(ibin>800+n*1024&&ibin<n*1024+950){
			sideband+= buildUp->GetBinContent(ibin);

		}
		else if(ibin>250+n*1024&&ibin<800+n*1024&&n<210&&n>100)
		{
			sideband+= buildUp->GetBinContent(ibin);
		}
		else if(ibin>210+n*1024&&ibin<800+n*1024&&n>=800&&n<950)
		{
			sideband+= buildUp->GetBinContent(ibin);
		}
		if(ibin%1024==0){
			n++;
		}

	}



	maximum = buildUp->GetMaximum();//Finds maximum
	for(int ibin=1; ibin <= 1024*1024; ibin++ ){//Removes the background from the signal
		double bkg_value = sideband*bkg_histBu->GetBinContent(ibin);
		double signalAndBkgValue = buildUp->GetBinContent(ibin);
		double signal = signalAndBkgValue - bkg_value;
		signal_histBu->SetBinContent(ibin,signal);
	}


	double sigIntegralBu(0);//Integrates signal
	for(int ibin=1; ibin <= 1024*1024; ibin++ ){
		sigIntegralBu+= signal_histBu->GetBinContent(ibin);
	}

	for(int ibin=1; ibin <= 1024*1024; ibin++ ){//Normalises the signal
		double normsig  = signal_histBu->GetBinContent(ibin)/sigIntegralBu;
		signal_histBu->SetBinContent(ibin, normsig);
	}







	std::cout<<"Hist made"<<std::endl;
	signal_histNoBu->Draw();
	signal_histNoBu->Write();
	signal_histBu->Draw();
	signal_histBu->Write();
	bkg_histNoBu->Draw();
	bkg_histNoBu->Write();
	bkg_histBu->Draw();
	bkg_histBu->Write();
	std::cout<<"written to file"<<std::endl;

	std::cout<<signal_histNoBu->GetNbinsX()<<" "<<signal_histNoBu->GetNbinsY()<<std::endl;
	for( int i =1; i<=1024; i++){
		for( int  j=1;j<=1024;j++){
			double noBu = signal_histNoBu->GetBinContent(i,j)*1000000;
			double bu = signal_histBu->GetBinContent(i,j)*1000000;
			//std::cout<<bu<<std::endl;
			buImage->SetPixelAt((i-1)*1024+j-1, bu);
			noBuImage->SetPixelAt((i-1)*1024+j-1,noBu);
		}
	}

	inputFile->Close();


}


void HistogramAnalysis::TwoDFit(std::shared_ptr< stk::Image< float > >image){





	RooRealVar x("x", "x", -0.5, 1023.5);//Define Roo variable
	RooRealVar y("y", "y", -0.5, 1023.5);//Define Roo variable
	TH2* noBuildUp = (TH2*)inputFile->Get("NoBu");//Load in the two histogram


	RooDataHist  SliceRooNoBu("No Bu Hist", "No Bu Hist",RooArgList(x,y), noBuildUp);//Creating two roo hists for slices


	RooRealVar x_mean("x_Mean","mean of gaussian x",512,0,1024);//Setting up gaussian parameters
	RooRealVar x_sigma("x_Sigma","width of gaussian x",256, 100, 600);
	RooGaussian gx("gx","gx",x,x_mean,x_sigma);//Setting up gaussian

	RooRealVar y_mean("y_Mean","mean of gaussian y",512,0,1024);//Setting up gaussian parameters
	RooRealVar y_sigma("y_Sigma","width of gaussian y",256, 100, 600);
	RooGaussian gy("gy","gy",y,y_mean,y_sigma);//Setting up gaussian

	std::cout<<"PDFs Created"<<std::endl;





	RooProdPdf backgroundModel("bkgModel","bkgModel",gx,gy);
	//	RooProdPdf signalModel("signalModel", "signalModel",tophat_x, tophat_y);


	x.setRange("LFT_BU",-0.5,220.5);
	y.setRange("LFT_BU",-0.5,1023.5);

	x.setRange("RGT_BU",800.5,1023.5);
	y.setRange("RGT_BU",-0.5,1023.5);

	x.setRange("T_BU",220.5,800.5);
	y.setRange("T_BU",800.5,1023.5);

	x.setRange("BM_BU",220.5,800.5);
	y.setRange("BM_BU",-0.5,200);




	std::cout<<"Models combined"<<std::endl;

	backgroundModel.fitTo(SliceRooNoBu,Range("LFT_BU, RGT_BU, T_BU, BM_BU"));
	TH1* bkg_histNoBu = backgroundModel.createHistogram("Background Hist No Bu",x,Binning(1024),YVar(y,Binning(1024)));



	TH2D* signal_histNoBu = new TH2D("Signal Hist No Bu", "Signal Hist", 1024, -0.5, 1023.5, 1024, -0.5, 1023.5);

	double integralNoBu(0);//Integrating the slice
	for(int ibin=1; ibin <= 1024*1024; ibin++ ){
		integralNoBu += noBuildUp->GetBinContent(ibin);
	}
	std::cout<<bkg_histNoBu->GetNcells()<<std::endl;


	double n(0);
	double sideband(0);
	for( int ibin=1; ibin<= 1024*1024; ibin++){
		if(ibin<220+n*1024&&ibin>n*1024){
			sideband+= noBuildUp->GetBinContent(ibin);

		}
		else if(ibin>800+n*1024&&ibin<(n+1)*1024){
			sideband+= noBuildUp->GetBinContent(ibin);

		}
		else if(ibin>220+n*1024&&ibin<800+n*1024&&n<220)
		{
			sideband+= noBuildUp->GetBinContent(ibin);
		}
		else if(ibin>220+n*1024&&ibin<800+n*1024&&n>=800)
		{
			sideband+= noBuildUp->GetBinContent(ibin);
		}
		if(ibin%1024==0){
			n++;
		}

	}



	double maximum = noBuildUp->GetMaximum();//Finds maximum
	for(int ibin=1; ibin <= 1024*1024; ibin++ ){//Removes the background from the signal
		double bkg_value = sideband*bkg_histNoBu->GetBinContent(ibin);
		double signalAndBkgValue = noBuildUp->GetBinContent(ibin);
		double signal = signalAndBkgValue - bkg_value;
		signal_histNoBu->SetBinContent(ibin,signal);
	}


	double sigIntegralNoBu(0);//Integrates signal
	for(int ibin=1; ibin <= 1024*1024; ibin++ ){
		sigIntegralNoBu+= signal_histNoBu->GetBinContent(ibin);
	}

	//	for(int ibin=1; ibin <= 1024*1024; ibin++ ){//Normalises the signal
	//		double normsig  = signal_histNoBu->GetBinContent(ibin)/sigIntegralNoBu;
	//		signal_histNoBu->SetBinContent(ibin, normsig);
	//	}
	//









	std::cout<<"Hist made"<<std::endl;
	signal_histNoBu->Draw();
	signal_histNoBu->Write();

	bkg_histNoBu->Draw();
	bkg_histNoBu->Write();

	std::cout<<"written to file"<<std::endl;

	std::cout<<signal_histNoBu->GetNbinsX()<<" "<<signal_histNoBu->GetNbinsY()<<std::endl;
	for( int i =1; i<=1024; i++){
		for( int  j=1;j<=1024;j++){
			double noBu = signal_histNoBu->GetBinContent(i,j);
			image->SetPixelAt((i-1)*1024+j-1,noBu);
		}
	}

	inputFile->Close();


}

void HistogramAnalysis::MaxMin(){

	TH1* slice = (TH1*)inputFile->Get("0");


	double m1(0),m2=1e12, m3(0);


	for(int i=1; i<1024;i++){

		if(i>50&&i<200){
			if(slice->GetBinContent(i)>m1){
				m1=slice->GetBinContent(i);
			}
		}

		if(i>800&&i<1020){
			if(slice->GetBinContent(i)>m3){
				m3=slice->GetBinContent(i);
			}

		}

		if(i>200&&i<800){
			if(slice->GetBinContent(i)<m2){
				m2=slice->GetBinContent(i);
			}
		}


	}
	double max(0);

	if(m1>m3){
		max=m1;
	}
	else
	{
		max=m3;
	}

	std::cout<<"Symmetry = "<<((m1/m3)-1)*100<<"% Flatness = "<<((max/m2)-1)*100<<"%"<<std::endl;

}


void HistogramAnalysis::OneDFit(std::shared_ptr< stk::Image< float > >image){


	TH2* noBuildUp = (TH2*)inputFile->Get("NoBu");//Load in the two histogram

	TH2F* output_signal = new TH2F("Output Signal", "Output Signal", 1024, -0.5,1023.5,  1024, -0.5,1023.5 );
	int slice(0);

	TH1F* mean = new TH1F("Mean", "Mean",1024, -0.5,1023.5 );
	TH1F* mean_dist = new TH1F("Mean dist", "Mean dist",70, 400,500 );
	TH1F* sigma = new TH1F("sigma", "sigma",1024, -0.5,1023.5 );
	TH1F* sigma_dist = new TH1F("sigma dist", "sigma dist",500, 280,1000 );

	gErrorIgnoreLevel = kWarning;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
	double mean_param(0), sig_param(0), sharp_param(0), topwidth_param(0), fsig_param(0), offset_param(0);
	for(slice=1; slice<=1024; slice++){

		RooRealVar x("x", "x", -0.5, 1023.5);//Define Roo variable
		TH1* inputNoBu = noBuildUp->ProjectionX("Y Slice NoBu", slice, slice);//Taking projection of first histogram
		TH1D* sliceNoBu = new TH1D("Y Slice w NoBu", "Y Slice", 1024, -0.5,1023.5);//Putting into a new TH1 histogram. It removes the error bar, not pointless
		for(int ibin=1; ibin <= inputNoBu->GetNbinsX(); ibin++ ){
			double value = inputNoBu->GetBinContent(ibin);
			sliceNoBu->SetBinContent(ibin,value);

		}





		RooDataHist  SliceRooNoBu("Y slice", "Y Slice",RooArgList(x), sliceNoBu);//Creating two roo hists for slices
		std::stringstream ss;
		ss <<slice;
		const char* num = ss.str().c_str();

		RooPlot* frame = x.frame(Title(num), Bins(1024) );//Creating two frames to plot slices on
		SliceRooNoBu.plotOn( frame );//Plotting slices on frames

		RooRealVar x_mean("x_Mean","mean of gaussian",480,320,630);//Setting up gaussian parameters
		RooRealVar x_sigma("x_Sigma","width of gaussian",391,200,800);
		RooGaussian gx("gx","gx",x,x_mean,x_sigma);//Setting up gaussian
		RooRealVar norm("norm", "norm",2.6e7,0,10000000000);
		RooExtendPdf gE("gE","gE", gx,norm);



		x.setRange("BU", 100.5,950.5);

		RooMinuit m(gE);
		m.setPrintLevel(-1000);
		x.setRange("LFT", -0.5,250.5);//Sets range for buildup option
		x.setRange("RGT", 700.5, 1023.5);
		x.setRange("FULL", -0.5, 1023.5);

		gE.fitTo(SliceRooNoBu,Range("LFT,RGT"),Verbose(kFALSE), PrintLevel(-1));



		RooAbsReal* lftInt = gE.createIntegral(x, Range("LFT")) ;
		RooAbsReal *fullInt = gE.createIntegral(x, Range("FULL")) ;
		std::cout << "gx_Int[x] = " << lftInt->getVal() <<std::endl ;
		std::cout<< fullInt<<std::endl;
		double left = lftInt->getVal();
		double full = fullInt->getVal();
		double ratio = full/left;
		std::cout<<"lft/full = "<<ratio<<std::endl;


		//gx.plotOn(frame, Range(-0.5, 1023.5),NormRange("LFT,RGT") );

		//model.fitTo(SliceRooNoBu,Range("BU"),PrintLevel(-1), Verbose(kFALSE));//Fits combined model to first hist
		if(slice==270||slice==340||slice==440||slice==512||slice==640||slice==720||slice==120){
			mean_param=x_mean.getValV(), sig_param=x_sigma.getValV();

			gE.plotOn(frame,Range("FULL"));
			//model.plotOn(frame);
			frame->Draw();
			frame->Write();

		}

		if(slice==490){
			mean_param=x_mean.getValV(), sig_param=x_sigma.getValV();

		}






		//	RooPlot* frame3 = x.frame(Title("PDF Components"));//Creates new frame

		//	model.plotOn(frame3, Components(g),LineColor(kRed));//Plots the two components of the combined function on a seperate fram
		//	model.plotOn(frame3, Components(gx), LineColor(kBlue));



		TH1D* signal_histNoBu = new TH1D("Signal Hist No Bu", "Signal Hist", 1024, -0.5, 1023.5);
		TH1* bkg_histNoBu = gE.createHistogram("Background HistNo Bu", x, Binning(1024));//Creates histogram from gaussian bkg fit
		if(slice==512){

			bkg_histNoBu->Draw();
			bkg_histNoBu->Write();
		}

		double sideBand(0);//Integrating the slice
		for(int ibin=1; ibin <= bkg_histNoBu->GetNbinsX(); ibin++ ){

			if( (ibin<=120) && (ibin>=1)){
				sideBand += sliceNoBu->GetBinContent(ibin);
			}
			//						if((ibin>=830) && (ibin<=1024)){
			//							sideBand += sliceNoBu->GetBinContent(ibin);
			//						}
		}

		double integral(0);//Integrating the slice
		for(int ibin=1; ibin <= bkg_histNoBu->GetNbinsX(); ibin++ ){


			integral += sliceNoBu->GetBinContent(ibin);

		}
		//		double bkgIntegral(0);
		//		for( int ibin=1; ibin<=bkg_histNoBu->GetNbinsX(); ibin++){
		//
		//			double temp = bkg_histNoBu->GetBinContent(ibin);
		//			bkg_histNoBu->SetBinContent(ibin, temp*integral);
		//			bkgIntegral+=temp*integral;
		//
		//		}
		//		double ratio = bkgIntegral/integral;

		std::cout<<"side band "<<sideBand<<"norm "<<norm.getValV()<<std::endl;
		double intVal = bkg_histNoBu->GetBinContent(1);
		std::cout<<"intval "<<integral<<std::endl;
		std::cout<<"Gx "<<gE.getVal()<<std::endl;
		RooArgSet nset(x);
		std::cout<<"gx norm "<<gE.getVal(&nset)<<std::endl;

		double maximum = sliceNoBu->GetMaximum();//Finds maximum
		for(int ibin=1; ibin <= bkg_histNoBu->GetNbinsX(); ibin++ ){//Removes the background from the signal
			double bkg_value =ratio*bkg_histNoBu->GetBinContent(ibin);
			double signalAndBkgValue = sliceNoBu->GetBinContent(ibin);
			double signal = signalAndBkgValue - bkg_value;
			if(ibin ==1){
				std::cout<<"signal "<<signal<<"bkg "<<bkg_value<<"sigbkg "<<signalAndBkgValue<<std::endl;
			}
			signal_histNoBu->SetBinContent(ibin,signal);
			output_signal->Fill(ibin, slice, signal);
		}


		double sigIntegralNoBu(0);//Integrates signal
		for(int ibin=1; ibin <= bkg_histNoBu->GetNbinsX(); ibin++ ){
			sigIntegralNoBu+= signal_histNoBu->GetBinContent(ibin);
		}

		for(int ibin=1; ibin <= bkg_histNoBu->GetNbinsX(); ibin++ ){//Normalises the signal
			double normsig  = signal_histNoBu->GetBinContent(ibin);
			signal_histNoBu->SetBinContent(ibin, normsig);
		}


		std::cout<<"Slice : "<<slice<<std::endl;
		std::cout<<"\t Mean = "<<x_mean.getValV()<<std::endl;
		std::cout<<"\t Sigma = "<<x_sigma.getValV()<<std::endl;

		mean->Fill(slice,x_mean.getValV());
		sigma->Fill(slice,x_sigma.getValV());




		delete inputNoBu;
		delete bkg_histNoBu;
		delete sliceNoBu;
		delete signal_histNoBu;



	}

	for( int i =1; i<=1024; i++){
		for( int  j=1;j<=1024;j++){
			double noBu = output_signal->GetBinContent(i,j);
			image->SetPixelAt((i-1)*1024+j-1,noBu);
		}
	}
	std::cout<<"Middle Slice"<<std::endl;
	std::cout<<"\t Mean \t "<<mean_param<<std::endl;
	std::cout<<"\t Sigma \t "<<sig_param<<std::endl;


	mean->Draw("hist");
	sigma->Draw("hist");

	mean_dist->Draw("hist");
	sigma_dist->Draw("hist");

	mean->Write();
	sigma->Write();

	mean_dist->Write();
	sigma_dist->Write();
	output_signal->Draw();
	output_signal->Write();
	inputFile->Close();




}


std::shared_ptr< stk::Image<float> > HistogramAnalysis::kFunc(){

	std::shared_ptr< stk::Image<float> > k( new stk::Image<float>(6,6));
	double temp(0);
	for( m_kernelParam.xs=-2.5; m_kernelParam.xs<3.5; m_kernelParam.xs++){
		for( m_kernelParam.ys=-2.5; m_kernelParam.ys<3.5;m_kernelParam.ys++){
			double kval(0);

			if(r()==0){
				kval = 1;
			}
			else
			{
				kval = 111.56*std::exp(-5.91158e-02*m_kernelParam.z*(1-1.32742e-08*m_kernelParam.z))*(1/(2*M_PI*r()))*((1-std::exp(-m_kernelParam.alpha*m_kernelParam.z))*(m_kernelParam.alpha*(1-m_kernelParam.F)*std::exp(-m_kernelParam.alpha*r())+m_kernelParam.beta*m_kernelParam.F*std::exp(-m_kernelParam.beta*r()))+((m_kernelParam.a*std::pow(m_kernelParam.z,2))/std::pow(m_kernelParam.omega*r()+m_kernelParam.z,2)));
			}

			k->SetPixelAt((m_kernelParam.xs+2.5)*6+m_kernelParam.ys+2.5,kval);
			temp+=kval;
		}
	}

	return k;
}

void  HistogramAnalysis::ErrorFunction(){

	TH2* noBuildUp = (TH2*)inputFile->Get("NoBu");//Load in the two histogram

	TH1* inputNoBu = noBuildUp->ProjectionX("Y Slice NoBu", 512, 512);//Taking projection of first histogram
	TH1D* sliceNoBu = new TH1D("Y Slice w NoBu", "Y Slice", 1024, -0.5,1023.5);//Putting into a new TH1 histogram. It removes the error bar, not pointless
	TH1D* normalised = new TH1D("normalised", "normalised", 1024, -0.5,1023.5);//Putting into a new TH1 histogram. It removes the error bar, not pointless
	TH1D* normalised_max = new TH1D("normalised_max", "normalised_max", 1024, -0.5,1023.5);//Putting into a new TH1 histogram. It removes the error bar, not pointless
	double integral(0);
	double max(0);
	for(int ibin=1; ibin <= inputNoBu->GetNbinsX(); ibin++ ){
		double value = inputNoBu->GetBinContent(ibin);
		integral+=value;
		if(value>max){
			max=value;
		}
		sliceNoBu->SetBinContent(ibin,value);
	}

	for(int i=1; i<=inputNoBu->GetNbinsX(); i++){
		double value = inputNoBu->GetBinContent(i);
		normalised->SetBinContent(i,value/integral );
		normalised_max->SetBinContent(i,value/max);
	}
	normalised_max->Draw();
	normalised_max->Write();
	normalised->Draw();
	normalised->Write();
	//	TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);

	//	TF1 *f2 = new TF1("f2", "[0]*ROOT::Math::erfc((x-[1])/[3])+[2]	", 680, 820);
	//
	//	f2->SetParameters(25000000, 750, 12500000, 0.5);
	//	sliceNoBu->Fit("f1","r","m");
	//	sliceNoBu->Fit("f2","r","m");
	//	sliceNoBu->Draw("hist");
	//	f1->Draw("sames");
	//	f2->Draw("sames");
	//	c1->Update();
	//	c1->Draw();
	//	c1->Write();
	//	f1->Write();
	//	f2->Write();
	//	sliceNoBu->Write();

	RooRealVar x_lft("x_lft", "x_lft", 200.5, 310.5);
	RooRealVar a_lft("a_lft", "a_lft",350,0, 500);
	RooRealVar b_lft("b_lft", "b_lft",27,0,40);
	RooRealVar c_lft("c_lft", "c_lft",1.3e6,1e6 , 1e7);
	RooRealVar d_lft("d_lft", "d_lft",9.2e5, 1e5, 1e7);

	RooRealVar x_mean_lft("x_Mean_lft","mean of gaussian",340,0,900);//Setting up gaussian parameters
	RooRealVar x_sigma_lft("x_Sigma_lft","width of gaussian",300, 0, 1024);
	RooGaussian gx_lft("gx","gx",x_lft,x_mean_lft,x_sigma_lft);//Setting up gaussian
	RooGenericPdf erf("erf", "erf", "(ROOT::Math::erf((x_lft-a_lft)/b_lft))", RooArgList(x_lft,a_lft,b_lft));

	RooRealVar fsig_lft("fsig_lft","signal fraction", 0.5,0,1);//Variabel for PDF fraction
	RooAddPdf model("model", "model", RooArgList(gx_lft,erf), fsig_lft);//Combines the gaussian and tophat PDF
	RooPlot* frame_lft = x_lft.frame(Title("_lft"), Bins(1024) );
	RooDataHist  SliceRooNoBu_lft("Y slice", "Y Slice",RooArgList(x_lft), sliceNoBu);
	model.fitTo(SliceRooNoBu_lft);

	SliceRooNoBu_lft.plotOn(frame_lft);
	model.plotOn(frame_lft);
	frame_lft->Write();




	RooRealVar x_rgt("x_rgt", "x_rgt", 780.5, 880.5);
	RooRealVar a_rgt("a_rgt", "a_rgt",750,500, 1023);
	RooRealVar b_rgt("b_rgt", "b_rgt",27,0,40);
	RooRealVar c_rgt("c_rgt", "c_rgt",4.2e5,1e5 , 1e7);
	RooRealVar d_rgt("d_rgt", "d_rgt",9.2e5, 1e5, 1e7);

	RooRealVar x_mean_rgt("x_Mean_rgt","mean of gaussian",340,0,900);//Setting up gaussian parameters
	RooRealVar x_sigma_rgt("x_Sigma_rgt","width of gaussian",300, 0, 1024);
	RooGaussian gx_rgt("gx_rgt","gx_rgt",x_rgt,x_mean_rgt,x_sigma_rgt);//Setting up gaussian
	RooGenericPdf erfc("erfc", "erfc", "(ROOT::Math::erfc((x_rgt-a_rgt)/b_rgt))", RooArgList(x_rgt,a_rgt,b_rgt));

	RooRealVar fsig_rgt("fsig_rgt","signal fraction", 0.5,0,1);//Variabel for PDF fraction
	RooAddPdf model_rgt("model_rgt", "model_rgt", RooArgList(gx_rgt,erfc), fsig_rgt);//Combines the gaussian and tophat PDF
	RooPlot* frame_rgt = x_rgt.frame(Title("_rgt"), Bins(1024) );
	RooDataHist  SliceRooNoBu_rgt("Y slice_rgt", "Y Slice",RooArgList(x_rgt), sliceNoBu);
	model_rgt.fitTo(SliceRooNoBu_rgt);


	RooPlot* frameComp_lft = x_lft.frame(Title("PDF Components left"));//Creates new frame

	model.plotOn(frameComp_lft, Components(gx_lft),LineColor(kRed));//Plots the two components of the combined function on a seperate fram
	model.plotOn(frameComp_lft, Components(erf), LineColor(kBlue));

	RooPlot* frameComp_rgt = x_rgt.frame(Title("PDF Components right"));//Creates new frame

	model_rgt.plotOn(frameComp_rgt, Components(gx_rgt),LineColor(kRed));//Plots the two components of the combined function on a seperate fram
	model_rgt.plotOn(frameComp_rgt, Components(erfc), LineColor(kBlue));

	SliceRooNoBu_rgt.plotOn(frame_rgt);
	model_rgt.plotOn(frame_rgt);
	frame_rgt->Write();
	frameComp_rgt->Write();
	frameComp_lft->Write();
	std::cout<<"Name \t Value \t  Error :"<<std::endl;
	std::cout<<"Err_sig_left \t"<<b_lft.getValV()<<" \t"<<b_lft.getError()<<std::endl;
	std::cout<<"Gaus_sig_left \t"<<x_sigma_lft.getValV()<<" \t"<<x_sigma_lft.getError()<<std::endl;
	std::cout<<"Gaus_mean_left \t"<<x_mean_lft.getValV()<<" \t"<<x_mean_lft.getError()<<std::endl;
	std::cout<<"fsig_left \t"<<fsig_lft.getValV()<<" \t"<<fsig_lft.getError()<<std::endl;
	std::cout<<"Err_sig_right \t"<<b_rgt.getValV()<<" \t"<<b_rgt.getError()<<std::endl;
	std::cout<<"Gaus_sig_right \t"<<x_sigma_rgt.getValV()<<" \t"<<x_sigma_rgt.getError()<<std::endl;
	std::cout<<"Gaus_mean_right \t"<<x_mean_rgt.getValV()<<" \t"<<x_mean_rgt.getError()<<std::endl;
	std::cout<<"fsig_right \t"<<fsig_rgt.getValV()<<" \t"<<fsig_rgt.getError()<<std::endl;




	TF1 *funca=new TF1("fit",GaussRangeBig,-0.5,1023.5,3);
	//	funca->SetParameters(512,256,314, 40, 750, 40,100e3,1500e3,9.1500e3);
	//	funca->SetParLimits(0, 0, 1000);
	//	funca->SetParNames("Gaussian_Mean","Gaussian_Sigma", "Erf_Offset", "Erf_Sigma", "Erfc_Offset", "Erfc_Sigma","Gaussian_Normalisation", "Efr_Normalisation", "Erfc_Normalisation");
	//	funca->SetParLimits(1, 0, 1000);
	//	funca->SetParLimits(2, 0, 1000);
	//	funca->SetParLimits(3, 30, 60);
	//	funca->SetParLimits(4, 0, 1000);
	//	funca->SetParLimits(5, 30, 60);
	//	funca->SetParLimits(6, 0, 5e6);
	//	funca->SetParLimits(7, 1e4, 1e7);
	//	funca->SetParLimits(8, 1e4, 1e7);
	//	funca->SetLineColor(3);
	//	funca->SetLineWidth(2);
	funca->SetParameters(640, 500, 1e5);
	funca->SetParLimits(0, 300, 700);
	funca->SetParLimits(1, 300, 1500);
	funca->SetParLimits(2, 0, 5e6);
	reject = kTRUE;
	sliceNoBu->Fit("fit");
	sliceNoBu->Draw("hist");
	funca->Draw("sames");
	funca->Draw();
	funca->Write();
	sliceNoBu->Write();
	TF1* bkg_hist = new TF1("bkg_hist", Gauss,-0.5, 1023.5,3);
	bkg_hist->SetParameter(0, funca->GetParameter(0));
	bkg_hist->SetParameter(1, funca->GetParameter(1));
	bkg_hist->SetParameter(2, funca->GetParameter(2));
	bkg_hist->Draw();
	bkg_hist->Write();
	std::cout<<funca->GetParameter(0)<<std::endl;

	TH1F* bkg = new TH1F("bkg", "bkg", 1024, -0.5 , 1023.5);

	bkg->Eval(bkg_hist);


	bkg->Draw();
	bkg->Write();
	TH1D* signal_hist = new TH1D("signal_hist", "signal_hist", 1024, -0.5, 1023.5);
	for (int i=1; i<=1024; i++){
		double temp = bkg->GetBinContent(i);
		double signal = sliceNoBu->GetBinContent(i)- temp;
		signal_hist->SetBinContent(i, signal);


	}





	for( int slice =1; slice<=1024;slice++){
		std::stringstream ss;
		ss <<slice;
		const char* num = ss.str().c_str();
		if(slice==775||slice==760||slice==740||slice==720||slice==520||slice==310||slice==295||slice==280){
			TH1D* temp = new TH1D(num, num, 1024,-0.5,1023.5);
			temp = noBuildUp->ProjectionX("Y Slice NoBu", slice, slice);
			temp->SetTitle(num);
			temp->Write();
		}

	}


	TF1 *f1 = new TF1("f1",GaussRangeSmall, -0.5, 1023.5, 3);
	f1->SetParameters(512, 400	, 1e5);
	f1->SetParLimits(0, 0,1000 );
	f1->SetParLimits(1, 0, 1000);
	signal_hist->Fit(f1);
	signal_hist->Draw("hist");
	f1->Draw("sames");
	f1->Draw();
	f1->Write();


	signal_hist->Write();

	TF1* bkg2_hist = new TF1("bkg2_hist", Gauss,-0.5, 1023.5,3);
	bkg2_hist->SetParameter(0,f1->GetParameter(0));
	bkg2_hist->SetParameter(1,f1->GetParameter(1));
	bkg2_hist->SetParameter(2,f1->GetParameter(2));

	TH1F* bkg2 = new TH1F("bkg2", "bkg", 1024, -0.5 , 1023.5);

	bkg2->Eval(bkg2_hist);


	bkg2->Draw();
	bkg2->Write();
	TH1D* signal_hist2 = new TH1D("signal2_hist", "signal_hist", 1024, -0.5, 1023.5);
	for (int i=1; i<=1024; i++){
		double temp = bkg2->GetBinContent(i);
		double signal = signal_hist->GetBinContent(i)- temp;
		signal_hist2->SetBinContent(i, signal);


	}



	TF1 *erfn = new TF1("erfn", Erfn, -0.5, 1023.5, 6);
	erfn->SetParameters(1e5, 400, 20, 8e5, 660, 20);
	erfn->SetParNames("Left Normalisation", "Left Offset", "Left Width","Right Normalisation", "Right Offset", "Right Width");
	erfn->SetParLimits(0, 1e5, 5e6);
	erfn->SetParLimits(1, 0, 512);
	erfn->SetParLimits(2, 1, 50);
	erfn->SetParLimits(3, 1e5, 5e6);
	erfn->SetParLimits(4, 512, 1024);
	erfn->SetParLimits(5, 1, 50);
	signal_hist2->Fit(erfn);
	signal_hist2->Draw("hist");
	erfn->Draw("sames");
	erfn->Write();

	signal_hist2->Write();


	std::cout<<"Name \t Value \t Error"<<std::endl;
	std::cout<<erfn->GetParName(0)<<" \t"<<erfn->GetParameter(0)<<" \t"<<erfn->GetParError(0)<<std::endl;
	std::cout<<erfn->GetParName(1)<<" \t"<<erfn->GetParameter(1)<<" \t"<<erfn->GetParError(1)<<std::endl;
	std::cout<<erfn->GetParName(2)<<" \t"<<erfn->GetParameter(2)<<" \t"<<erfn->GetParError(2)<<std::endl;
	std::cout<<erfn->GetParName(3)<<" \t"<<erfn->GetParameter(3)<<" \t"<<erfn->GetParError(3)<<std::endl;
	std::cout<<erfn->GetParName(4)<<" \t"<<erfn->GetParameter(4)<<" \t"<<erfn->GetParError(4)<<std::endl;
	std::cout<<erfn->GetParName(5)<<" \t"<<erfn->GetParameter(5)<<" \t"<<erfn->GetParError(5)<<std::endl;
	std::cout<<"Chi2"<<" \t"<<erfn->GetChisquare()<<std::endl;

	std::cout<<noBuildUp->GetBinContent(110,528)<<" \t"<<noBuildUp->GetBinContent(504,528)<<" \t"<<noBuildUp->GetBinContent(896,528)<<std::endl;
}

void HistogramAnalysis::PeakValueFind(){


	TH2* noBuildUp = (TH2*)inputFile->Get("NoBu");//Load in the two histogram

	TH1* inputNoBu = noBuildUp->ProjectionX("Y Slice NoBu", 490, 490);//Taki
	inputNoBu->Draw("hist");
	inputNoBu->Write();
	TH1D* sliceNoBu = new TH1D("Y Slice w NoBu", "Y Slice", 1024, -0.5,1023.5);
	for(int ibin=1; ibin <= inputNoBu->GetNbinsX(); ibin++ ){
		double value = inputNoBu->GetBinContent(ibin);
		sliceNoBu->SetBinContent(ibin,value);
	}

	double lastMax(0);
	double lastMin(0);
	double check(0), maximum(0), min(0), avg(0), count(0);
	bool max = false;
	for (int i = 350; i<580; i++){
		double prev = sliceNoBu->GetBinContent(i-1);
		double cur = sliceNoBu->GetBinContent(i);
		double nxt = sliceNoBu->GetBinContent(i+1);
		if(cur>prev&&cur>nxt&&(i-lastMax)>3&&max!=TRUE)
		{
			std::cout<<"Max \t"<<cur<<" \t"<<i<<std::endl;
			lastMax =i;
			maximum=cur;
			max=TRUE;
			check++;

		}
		else if(cur<prev&&cur<nxt&&(i-lastMin)>3&&max==TRUE )
		{
			std::cout<<"Min \t"<<cur<<" \t"<<i<<std::endl;
			lastMin=i;
			max=FALSE;
			check++;
			min=cur;

		}
		if(check==2){
			std::cout<<"Ratio \t"<<maximum/min<<std::endl;
			avg+=(maximum/min);
			count++;
			check=0;
		}
	}

	std::cout<<"Avg: "<<avg/count<<" "<<count<<std::endl;


	TF1 *erfn = new TF1("erfn", Erfunc, 340, 360, 4);
	erfn->SetParameters(1e6, 355, 120,5e5);
	erfn->SetParNames("Left Normalisation", "Left Offset", "Left Width","Right Normalisation", "Right Offset", "Right Width");
	erfn->SetParLimits(0, 1e5, 5e7);
	erfn->SetParLimits(1, 300, 400);
	erfn->SetParLimits(2, 1, 200);
	erfn->SetParLimits(3, 0, 5e6);
	erfn->SetParLimits(4, 0, 50e6);

	sliceNoBu->Fit(erfn);
	sliceNoBu->Draw("hist");
	erfn->Draw("sames");
	erfn->Write();

}
double HistogramAnalysis::r(){

	return std::sqrt(std::pow((m_kernelParam.xs-m_kernelParam.x),2)+std::pow((m_kernelParam.ys-m_kernelParam.y),2));
}

void HistogramAnalysis::SlopeAnalysis(int x, int y){

	TH2* mainImage = (TH2*)inputFile->Get("mainImage");
	TH2* varianceImage = (TH2*)inputFile->Get("varianceImage");

	TH1* mainImageSliceX = mainImage->ProjectionX("mainImageSliceX", x, x);//Taking projection of first histogram
	TH1* varianceImageSliceX = varianceImage->ProjectionX("varianceImageSliceX", x, x );
	TH1* mainImageSliceY = mainImage->ProjectionY("mainImageSliceY", y, y);//Taking projection of first histogram
	TH1* varianceImageSliceY = varianceImage->ProjectionY("varianceImageSliceY", y, y );

	TH1F * XsliceWErrors  = new TH1F("XsliceWErrors","XsliceWErrors", 1024, -0.5, 1023.5 );
	TH1F * YsliceWErrors  = new TH1F("YsliceWErrors","YsliceWErrors", 1024, -0.5, 1023.5 );
	for(int i = 1; i<=1024; i++)
	{
		XsliceWErrors->SetBinContent(i, mainImageSliceX->GetBinContent(i));
		XsliceWErrors->SetBinError(i, varianceImageSliceX->GetBinContent(i));
		YsliceWErrors->SetBinContent(i, mainImageSliceY->GetBinContent(i));
		YsliceWErrors->SetBinError(i, varianceImageSliceY->GetBinContent(i));
	}
	double avgADC = 0;
	double avgADCErr = 0;
	for (int i = 0; i < 100; ++i)
	{
		avgADC+=XsliceWErrors->GetBinContent(450+i);
		avgADCErr+=XsliceWErrors->GetBinError(450+i);
	}
	avgADCErr=avgADCErr/100;
	avgADC=avgADC/100;
	std::cout<<"Mid point: "<<avgADC<<"±"<<avgADCErr<<std::endl;
	TF1 * fittingFunc = new TF1("fittingFunc", expFunc, 0, 1023 ,4);
	TF1 * fittingFuncrgt = new TF1("fittingFunclrgt", expFunc, 0, 1023 ,4);
	TF1 * errorFuncFit = new TF1("errrorfn",Erfn, 0 , 1023, 8 );
	errorFuncFit->SetParameters(2e5, 200, 20, 2e5, 700, 73);
	errorFuncFit->SetParNames("Left Normalisation", "Left Offset", "Left Width","Right Normalisation", "Right Offset", "Right Width");
	//	errorFuncFit->SetParLimits(0, 1e5, 5e6);
	//	errorFuncFit->SetParLimits(1, 0, 512);
	//	errorFuncFit->SetParLimits(2, 1, 50);
	//	errorFuncFit->SetParLimits(3, 1e5, 5e6);
	//	errorFuncFit->SetParLimits(4, 512, 1024);
	//	errorFuncFit->SetParLimits(5, 1, 50);

	// fittingFunc->SetParameter(0,9.72062e+01 );
	// //fittingFunc->SetParLimits(0, 9.72062e+01 , 9.72062e+01);
	// fittingFunc->SetParameter(1,1e4);
	// fittingFunc->SetParameter(2,2e5);
	// fittingFunc->SetParameter(3,9e1);



	TF1 * logistic  = new TF1("logistic", LogisticFunction, 0, 1023, 4);
	logistic->SetParameter(0, 1e6);
	logistic->SetParameter(2, 790);
	XsliceWErrors->Fit(logistic);
	//YsliceWErrors->Fit(errorFuncFit);

	XsliceWErrors->Write();
	YsliceWErrors->Write();
	std::cout<<logistic->GetChisquare()/logistic->GetNDF()<<std::endl;
	std::cout<<"Chi2 "<<logistic->GetChisquare()<<" NDF "<<logistic->GetNDF()<<std::endl;
	std::cout<<errorFuncFit->GetChisquare()/errorFuncFit->GetNDF()<<std::endl;
	std::cout<<"Chi2 "<<errorFuncFit->GetChisquare()<<" NDF "<<errorFuncFit->GetNDF()<<std::endl;
}



void HistogramAnalysis::OutsideBeam()
{
	TH2* mainImage = (TH2*)inputFile->Get("mainImage");
	TH1F *signalHist = new TH1F("signal", "signal", 100, 0, 10000);
	for(int x = 900; x <= 1024; x++)
	{
		for(int y = 0; y<=200; y++)
		{
			signalHist->Fill(mainImage->GetBinContent(x,y));
		}
	}

	signalHist->Write();
}
}


#endif /* MODULES_PLOTTING_HISTOGRAMANALYSIS_INCLUDE_STKHISTOGRAMANALYSIS_HXX_ */















