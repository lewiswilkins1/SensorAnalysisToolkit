/*
 * PeakValue.cxx
 *
 *  Created on: 21 Sep 2015
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


	analysis.PeakValueFind();


	std::cout<<"done"<<std::endl;

	return 1;
}


