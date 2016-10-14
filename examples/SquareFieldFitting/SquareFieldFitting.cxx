/*
 * SquareFieldFitting.cxx
 *
 *  Created on: 28 Aug 2015
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
#include "TFile.h"
#include "TH3.h"

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
#include "stkHistogramAnalysis.h"


int main(int argc, char*argv[])
{


	typedef itk::Image< float, 2 >  ImageType;
	typedef stk::Image<float> stkImageType;
	std::string filePathandfilename;
	int rows, cols, framesize;
	std::stringstream inputs;
	for(int iArg=1; iArg < argc; iArg++)
	{
		inputs << argv[iArg] << ' ';//get the arguments
	}

	inputs >>  filePathandfilename >>   rows >> cols;//write the inputs
	framesize = cols*rows;

	stk::IO<float> myIO;
	std::shared_ptr< stkImageType > myImage( new stkImageType );
	myImage->Initialise( myIO.ReadImage( filePathandfilename ), rows, cols );
	myImage->Resize();

	stk::ImageHistogram<TH2F, float> myImageHistogramX;

	myImageHistogramX.SetYAxisTitle("Row");
	myImageHistogramX.SetYAxisLog(false);
	myImageHistogramX.SetNumberOfYBins(1024);
	myImageHistogramX.SetYLimits( -0.5, 1023.5 );
	myImageHistogramX.SetXAxisTitle( "Col" );
	myImageHistogramX.SetNumberOfXBins( 1024);
	myImageHistogramX.SetXLimits(  -0.5, 1023.5  );

	myImageHistogramX.SetGridY(false);
	myImageHistogramX.SetGridX(false);

	myImageHistogramX.SetStatBoxOptions(0);
	myImageHistogramX.SetOutputFileNameAndPath("X");
	myImageHistogramX.SetOutputFileType( stk::FileType::PNG );

	myImageHistogramX.Generate2DHistogram( myImage, "NoBu");
	myImageHistogramX.SaveHistogram();
	myImageHistogramX.ExportRoot();

	stk::HistogramAnalysis analysis;

	analysis.LoadRootFile();
	analysis.OneDFit(myImage);





	ImageType::Pointer itkImage = ImageType::New();

	stk::ITKAdapter adapter;
	adapter.STKtoITK(itkImage, myImage, 1024, 1024);
	myImage->Delete();
	std::shared_ptr< stkImageType > kernel;
	kernel = analysis.kFunc();

	ImageType::Pointer itkKernel = ImageType::New();

	ImageType::Pointer itkConvolutionOutput = ImageType::New();
	std::shared_ptr< stkImageType > convolvedImage (new stkImageType(1024,1024));

	adapter.STKtoITK(itkKernel, kernel, 5, 5);


	typedef itk::ConvolutionImageFilter<ImageType> FilterType;
	FilterType::Pointer convolutionFilter = FilterType::New();


	convolutionFilter->SetInput(itkImage);
	convolutionFilter->SetKernelImage(itkKernel);
	convolutionFilter->Update();
	itkConvolutionOutput = convolutionFilter->GetOutput();
	itkConvolutionOutput->Update();
	itkConvolutionOutput->DisconnectPipeline();



	adapter.ITKtoSTK(itkConvolutionOutput,convolvedImage );

	TFile* file = new TFile("output.root", "RECREATE");
	TH2F* tempMap = new TH2F("TempMap", "Map",1024, -0.5, 1023.5,1024, -0.5, 1023.5 );
	TH2F* map = new TH2F("Map", "Map",1024, -0.5, 1023.5,1024, -0.5, 1023.5 );

	for(int i = 1; i<=1024; i++){
		for(int j=1; j<=1024; j++){
			double temp = convolvedImage->GetPixelAt(((i-1)*1024)+(j-1));
			tempMap->Fill(i,j, temp);
		}
	}
	double avg;
	for(int j=486; j<=494; j++){
		for(int i=491; i<=499; i++){

			avg+=tempMap->GetBinContent(i,j);
		}
	}
	avg=avg/81;

	tempMap->Draw();
	tempMap->Write();
	double binVal=tempMap->GetBinContent(504,528);
	std::cout<<"Central Pixel Value = "<<binVal<<std::endl;
	double sf =5.1314e-08;
	std::cout<<"Scale Factor = "<<sf<<std::endl;
	for(int i = 1; i<=1024; i++){
			for(int j=1; j<=1024; j++){
				double temp = tempMap->GetBinContent(i,j)*sf;
				map->SetBinContent(i,j,temp);
			}
		}
	double outputDose(0);
	for(int j=486; j<=494; j++){
			for(int i=491; i<=499; i++){
				outputDose+=map->GetBinContent(i,j);
			}
		}
	outputDose=outputDose/81;
	std::cout<<"Dose at centre = "<<outputDose<<std::endl;
	map->Draw();
	map->Write();
	return 1;

}


