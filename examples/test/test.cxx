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
#include "TMultiGraph.h"


int main(int argc, char*argv[])
{

	typedef stk::Image<float> stkImageType;
	stk::IO<float> myIO;
	std::shared_ptr< stkImageType > myImage( new stkImageType(4096, 4096) );
	myImage->Resize();
	myImage->Delete();
	TFile* file = new TFile("Graphs.root", "RECREATE");


	TGraph* noBU_dose = new TGraph(4);
	TGraph* dose = new TGraph(5);

	TMultiGraph *mg = new TMultiGraph();

	noBU_dose->SetPoint(0,3, 267360/330602);
	noBU_dose->SetPoint(1,4, 292944/330602);
	noBU_dose->SetPoint(2,6,304688/330602);
	noBU_dose->SetPoint(3,8, 1);
	noBU_dose->SetMarkerStyle(5);
	noBU_dose->SetMarkerColor(kRed);
	noBU_dose->SetLineColor(0);
	noBU_dose->SetTitle("No Build Up Dose");
	noBU_dose->SetDrawOption("P");
	noBU_dose->GetYaxis()->SetTitle("Output Factor");
	noBU_dose->GetXaxis()->SetTitle("Field Size (cm)");


	noBU_dose->Draw("p");
	noBU_dose->Write();


	dose->SetPoint(0,3,412922/515083);
	dose->SetPoint(1,4,438956/515083);
	dose->SetPoint(2,5,459216/515083);
	dose->SetPoint(3,6,486479/515083);
	dose->SetPoint(4,8,1);
	
	dose->SetMarkerStyle(7);
	//noBU_dose->SetMarkerColor(kRed);
	dose->SetLineColor(kRed);
	dose->GetYaxis()->SetTitle("Output Factor");
	dose->GetXaxis()->SetTitle("Field Size (cm)");
	dose->SetTitle("Documentation Ouput Factor");
	dose->SetDrawOption("P");
	dose->Draw("p");

	dose->Write();
	TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
	c1->SetGrid();
	c1->SetTitle("No Build Up Dose");
	mg->Add(noBU_dose);
	mg->Add(dose);

	mg->Draw("ap");
	mg->GetYaxis()->SetTitle("Output Factor");
	mg->GetXaxis()->SetTitle("Field Size (cm)");
	mg->SetTitle("No Build Up Dose");
	mg->Write();

	c1->BuildLegend();
	c1->Update();
	c1->Write();
	TCanvas *c2 = new TCanvas("c2", "c2", 1200, 800);

	TGraph* BU_dose = new TGraph(6);
	TMultiGraph *mg2 = new TMultiGraph();

	c2->SetGrid();
	BU_dose->SetPoint(0,2, 1.059238417);
	BU_dose->SetPoint(1,3, 1.088970822);
	BU_dose->SetPoint(2,4, 1.069642644);
	BU_dose->SetPoint(3,5, 1);
	BU_dose->SetPoint(4,6, 0.890914836);
	BU_dose->SetPoint(5,7, 0.737202452);

	BU_dose->SetMarkerStyle(5);
	BU_dose->SetMarkerColor(kRed);
	BU_dose->SetLineColor(0);
	BU_dose->SetTitle("Build Up Dose");
	BU_dose->GetYaxis()->SetTitle("Output Factor");
	BU_dose->GetXaxis()->SetTitle("Field Size (cm)");
	BU_dose->Draw("p");
	BU_dose->Write();

	mg2->Add(BU_dose);
	mg2->Add(dose);

	mg2->Draw("ap");
	mg2->GetYaxis()->SetTitle("Output Factor");
	mg2->GetXaxis()->SetTitle("Field Size (cm)");
	mg2->SetTitle("Build Up Dose");
	mg2->Write();
	c2->BuildLegend();
	c2->Write();

	TCanvas *c3 = new TCanvas("c3", "c3", 1200, 800);
	c3->SetGrid();

	TGraph* BU_sigma = new TGraph(6);

	BU_sigma->SetPoint(0,2, 385);
	BU_sigma->SetPoint(1,3, 354);
	BU_sigma->SetPoint(2,4, 327);
	BU_sigma->SetPoint(3,5, 316);
	BU_sigma->SetPoint(4,6, 331);
	BU_sigma->SetPoint(5,7, 333);

	BU_sigma->SetMarkerStyle(2);
	BU_sigma->SetLineColor(0);
	BU_sigma->GetXaxis()->SetTitle("Field Size (cm)");
	BU_sigma->GetYaxis()->SetTitle("Sigma (px)");
	BU_sigma->SetTitle("Sigma Values With Build Up");
	BU_sigma->Draw("AP");
	BU_sigma->Write();

	c3->Update();
	TCanvas *c4 = new TCanvas("c4", "c4", 1200, 800);
	c3->Write();
	TGraph* noBU_sigma = new TGraph(5);
	c4->SetGrid();
	noBU_sigma->SetPoint(0,2, 436);
	noBU_sigma->SetPoint(1,3, 429);
	noBU_sigma->SetPoint(2,4, 413);
	noBU_sigma->SetPoint(3,5, 401);
	noBU_sigma->SetPoint(4,6, 412);
	noBU_sigma->SetMarkerStyle(2);
	noBU_sigma->SetLineColor(0);
	noBU_sigma->GetXaxis()->SetTitle("Field Size (cm)");
	noBU_sigma->GetYaxis()->SetTitle("Sigma (px)");
	noBU_sigma->SetTitle("Sigma Values With No Build Up");
	noBU_sigma->Draw("AP");
	noBU_sigma->Write();
	c4->Write();



	TCanvas *c5 = new TCanvas("c5", "c5", 1200, 800);

	TMultiGraph *mg3 = new TMultiGraph();
	TGraph* leftSigma = new TGraph(8);
	TGraph* rightSigma = new TGraph(8);
	c5->SetGrid();
	leftSigma->SetPoint(0, 2, 21.245);
	leftSigma->SetPoint(1, 3, 24.1642);
	leftSigma->SetPoint(2, 4, 24.1031);
	leftSigma->SetPoint(3, 5, 23.4030);
	leftSigma->SetPoint(4, 6, 23.1145);
	leftSigma->SetPoint(5, 7, 25.253);
	leftSigma->SetPoint(6, 8, 25.1320);
	leftSigma->SetPoint(7, 9, 22.0187);
	leftSigma->SetLineColor(0);
	leftSigma->SetMarkerStyle(2);
	leftSigma->SetMarkerColor(kRed);
	leftSigma->SetTitle("Left Sigma Values");
	leftSigma->GetXaxis()->SetTitle("Field Size (cm)");
	leftSigma->GetYaxis()->SetTitle("Sigma (px)");
	leftSigma->Draw("ap");

	rightSigma->SetPoint(0, 2, 22.1498);
	rightSigma->SetPoint(1, 3, 21.3888);
	rightSigma->SetPoint(2, 4, 22.8794);
	rightSigma->SetPoint(3, 5, 20.5747);
	rightSigma->SetPoint(4, 6, 23.7880);
	rightSigma->SetPoint(5, 7, 22.4044);
	rightSigma->SetPoint(6, 8, 22.7710);
	rightSigma->SetPoint(7, 9, 23.2504);
	rightSigma->SetLineColor(0);
	rightSigma->SetMarkerStyle(5);
	rightSigma->SetMarkerColor(kBlue);
	rightSigma->SetTitle("Right Sigma Values");
	rightSigma->GetXaxis()->SetTitle("Field Size (cm)");
	rightSigma->GetYaxis()->SetTitle("Sigma (px)");
	rightSigma->Draw("ap");


	mg3->Add(leftSigma);
	mg3->Add(rightSigma);
	mg3->Draw("ap");
	mg3->GetYaxis()->SetTitle("Sigma");
	mg3->GetXaxis()->SetTitle("Field Size (cm)");
	mg3->Write();
	c5->Update();
	c5->Write();


	TCanvas *c6 = new TCanvas("c6", "c6", 1200, 800);

	TGraph* ratio_1000 = new TGraph(8);
	c6->SetGrid();
	ratio_1000->SetPoint(0, 2, 1.142);
	ratio_1000->SetPoint(1, 3, 1.118);
	ratio_1000->SetPoint(2, 4, 1.147);
	ratio_1000->SetPoint(3, 5, 1.150);
	ratio_1000->SetPoint(4, 6, 1.151);
	ratio_1000->SetPoint(6, 8, 1.157);
	ratio_1000->SetPoint(7, 9, 1.153);

	ratio_1000->SetLineColor(0);
	ratio_1000->SetMarkerStyle(5);
	ratio_1000->SetMarkerColor(kBlue);
	ratio_1000->SetTitle("Peak Ratios -Slice 1000");
	ratio_1000->GetXaxis()->SetTitle("Field Size (cm)");
	ratio_1000->GetYaxis()->SetTitle("Ratio Max/Min");

	ratio_1000->Draw("AP");
	c6->Update();
	c6->Write();

	TCanvas *c7 = new TCanvas("c7", "c7", 1200, 800);

	TGraph* ratio_512 = new TGraph(8);
	c7->SetGrid();
	ratio_512->SetPoint(0, 2, 1.08885);
	ratio_512->SetPoint(1, 3, 1.07962);
	ratio_512->SetPoint(2, 4, 1.06203);
	ratio_512->SetPoint(3, 5, 1.05016);
	ratio_512->SetPoint(4, 6, 1.03595);
	ratio_512->SetPoint(6, 8, 1.01575);
	ratio_512->SetPoint(7, 9, 1.0122);


	ratio_512->SetLineColor(0);
	ratio_512->SetMarkerStyle(5);
	ratio_512->SetMarkerColor(kBlue);
	ratio_512->SetTitle("Peak Ratios - Slice 490(Centre)");
	ratio_512->GetXaxis()->SetTitle("Field Size (cm)");
	ratio_512->GetYaxis()->SetTitle("Ratio Max/Min");

	ratio_512->Draw("AP");
	c7->Update();
	c7->Write();



	TCanvas *c8 = new TCanvas("c8", "c8", 1200, 800);

	TGraph* outputFactors = new TGraph(8);
	c8->SetGrid();
	outputFactors->SetPoint(0, 2, 0.1471/0.4020);
	outputFactors->SetPoint(1, 3, 0.1716/0.4020);
	outputFactors->SetPoint(2, 4, 0.3012/0.4020);
	outputFactors->SetPoint(3, 5, 0.4020/0.4020);
	outputFactors->SetPoint(4, 6, 0.6324/0.4020);
	outputFactors->SetPoint(6, 8, 0.6403/0.4020);
	outputFactors->SetPoint(7, 9, 0.6892/0.4020);


	outputFactors->SetLineColor(0);
	outputFactors->SetMarkerStyle(5);
	outputFactors->SetMarkerColor(kBlue);
	outputFactors->SetTitle("Ouput Factors - Normalised to 5x5 Field");
	outputFactors->GetXaxis()->SetTitle("Field Size (cm)");
	outputFactors->GetYaxis()->SetTitle("Output factor");

	TMultiGraph *mg4 = new TMultiGraph();
	mg4->Add(dose);
	mg4->Add(outputFactors);
	outputFactors->Draw("AP");

	mg4->Draw("ap");
		mg4->GetYaxis()->SetTitle("Output Factor");
		mg4->GetXaxis()->SetTitle("Field Size (cm)");
		mg4->SetTitle("");
		mg4->Write();
	c8->Update();
	c8->Write();





	TCanvas *c9 = new TCanvas("c9", "c9", 1200, 800);

	TGraph* sigma_noBU_convoluted = new TGraph(3);
	c9->SetGrid();
	sigma_noBU_convoluted->SetPoint(0, 4, 354.27);
	sigma_noBU_convoluted->SetPoint(1, 5, 382.673);
	sigma_noBU_convoluted->SetPoint(2, 6, 417.437);



	sigma_noBU_convoluted->SetLineColor(0);
	sigma_noBU_convoluted->SetMarkerStyle(5);
	sigma_noBU_convoluted->SetMarkerColor(kBlue);
	sigma_noBU_convoluted->SetTitle("Gaussian widths - NoBU - tophat/gaussian");
	sigma_noBU_convoluted->GetXaxis()->SetTitle("Field Size (cm)");
	sigma_noBU_convoluted->GetYaxis()->SetTitle("Width (px)");

	sigma_noBU_convoluted->Draw("AP");
	c9->Update();
	c9->Write();




	TCanvas *c10 = new TCanvas("c10", "c10", 1200, 800);

	TGraph* sigma_BU_convoluted = new TGraph(3);
	c10->SetGrid();
	sigma_BU_convoluted->SetPoint(0, 4,293.365);
	sigma_BU_convoluted->SetPoint(1, 5, 345.679);
	sigma_BU_convoluted->SetPoint(2, 6, 404.346);



	sigma_BU_convoluted->SetLineColor(0);
	sigma_BU_convoluted->SetMarkerStyle(5);
	sigma_BU_convoluted->SetMarkerColor(kBlue);
	sigma_BU_convoluted->SetTitle("Gaussian widths - BU - tophat/gaussian");
	sigma_BU_convoluted->GetXaxis()->SetTitle("Field Size (cm)");
	sigma_BU_convoluted->GetYaxis()->SetTitle("Width (px)");

	sigma_BU_convoluted->Draw("AP");
	c10->Update();
	c10->Write();



	TCanvas *c11 = new TCanvas("c11", "c11", 1200, 800);

	TGraph* sigma_BU_erfc = new TGraph(7);
	c11->SetGrid();
	sigma_BU_erfc->SetPoint(0, 2, 265.479);
	sigma_BU_erfc->SetPoint(1, 3, 267.348);
	sigma_BU_erfc->SetPoint(2, 4, 288.046);
	sigma_BU_erfc->SetPoint(3, 5, 300.984);
	sigma_BU_erfc->SetPoint(4, 6, 334.529);
	sigma_BU_erfc->SetPoint(5, 7, 379.4);
	sigma_BU_erfc->SetPoint(6, 8, 396.039);


	sigma_BU_erfc->SetLineColor(0);
	sigma_BU_erfc->SetMarkerStyle(5);
	sigma_BU_erfc->SetMarkerColor(kBlue);
	sigma_BU_erfc->SetTitle("Gaussian widths - BU - erfn");
	sigma_BU_erfc->GetXaxis()->SetTitle("Field Size (cm)");
	sigma_BU_erfc->GetYaxis()->SetTitle("Width (px)");

	sigma_BU_erfc->Draw("AP");
	c11->Update();
	c11->Write();


	file->Close();
	return 1;
}
