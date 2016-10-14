/*
 * ROOTAnalysis.cxx
 *
 *  Created on: 22 Apr 2016
 *      Author: lewish
 */




/*
 * ErrorFunctionFit.cxx
 *
 *  Created on: 14 Sep 2015
 *      Author: lewish
 */





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

#include "stkHistogramAnalysis.h"




int main(int argc, char*argv[])
{



	typedef stk::Image<float> stkImageType;

	std::string filePathandfilename, varianceFilePathandfilename;
	int rows=4096, cols=4096, framesize, x,y;
	std::stringstream inputs;
	for(int iArg=1; iArg < argc; iArg++)
	{
		inputs << argv[iArg] << ' ';//get the arguments
	}

	inputs >>  filePathandfilename >> x >> y;//write the inputs
	framesize = cols*rows;
	varianceFilePathandfilename = filePathandfilename + "_Variance.raw";
	std::string rootName = filePathandfilename;
	std::cout<<rootName<<std::endl;
	filePathandfilename.append(".raw");

	stk::IO<float> myIO;
	std::shared_ptr< stkImageType > mainImage( new stkImageType );
	mainImage->Initialise( myIO.ReadImage( filePathandfilename ), rows, cols );
	mainImage->Resize();

	std::shared_ptr< stkImageType > varianceImage( new stkImageType );
	varianceImage->Initialise( myIO.ReadImage( varianceFilePathandfilename ), rows, cols );
	varianceImage->Resize();



	stk::ImageHistogram<TH2F, float> mainImageHistogram;
	mainImageHistogram.SetYAxisTitle("Row");
	mainImageHistogram.SetYAxisLog(false);
	mainImageHistogram.SetNumberOfYBins(1024);
	mainImageHistogram.SetYLimits( -0.5, 1023.5 );
	mainImageHistogram.SetXAxisTitle( "Col" );
	mainImageHistogram.SetNumberOfXBins( 1024);
	mainImageHistogram.SetXLimits(  -0.5, 1023.5  );
	mainImageHistogram.SetGridY(false);
	mainImageHistogram.SetGridX(false);
	mainImageHistogram.SetStatBoxOptions(0);
	mainImageHistogram.Generate2DHistogram( mainImage, "mainImage");
	mainImageHistogram.Generate2DHistogram( varianceImage, "varianceImage");
	mainImageHistogram.ExportRoot(rootName);









	stk::HistogramAnalysis analysis;

	analysis.LoadRootFile(rootName);


	analysis.SlopeAnalysis(x,y);
	//analysis.OutsideBeam();

	std::cout<<"done"<<std::endl;




	return 1;
}
