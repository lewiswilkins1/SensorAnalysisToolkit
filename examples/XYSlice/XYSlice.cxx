/*
 * XYSlice.cxx
 *
 *  Created on: 4 Aug 2015
 *      Author: lewish
 */



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
#include "stkImageFlatFieldCorrection.h"
#include "stkHistogramAnalysis.h"


int main(int argc, const char** argv){

	std::string filePathandfilename, otherfilename;
	int startingframe, numOfFrames, rows, cols, framesize;

	std::stringstream inputs;
	for(int iArg=1; iArg < argc; iArg++)
	{
		inputs << argv[iArg] << ' ';//get the arguments
	}

	inputs >> filePathandfilename>> otherfilename  >> rows >> cols ;//write the inputs
	framesize = cols*rows;


	stk::IO<float> myIO;
	std::shared_ptr< stk::Image<float> > myImage( new stk::Image<float> );
	myImage->Initialise( myIO.ReadImage( filePathandfilename ), rows, cols );


	stk::IO<float> myIO2;
		std::shared_ptr< stk::Image<float> > myImage2( new stk::Image<float> );
		myImage2->Initialise( myIO2.ReadImage( otherfilename ), rows, cols );


	//Image to hold resized result
	std::shared_ptr<stk::Image<float> > myResultResize (new stk::Image<float>(1024,1024));
	std::shared_ptr<stk::Image<float> > myResultResize2 (new stk::Image<float>(1024,1024));

	stk::ImageResize mySize;

	//Image resizing
	mySize.ResizeImage(myImage , myResultResize);
	mySize.ResizeImage(myImage2 , myResultResize2);

		typename std::vector<float>::iterator imageStart;
		typename std::vector<float>::iterator imageEnd;

		imageStart=myResultResize->StartImage();
		imageEnd=myResultResize->EndImage();


		float min(0), max(0), min2(0), max2(0);
		min = *std::min_element(imageStart,imageEnd);
		max = *std::max_element(imageStart,imageEnd);


		std::cout<<"Max: "<<max<<" Min: "<<min<<std::endl;




	for(int iElements=0; iElements<1024*1024; iElements++){
		float Temp, temp2;
		Temp = (myResultResize->GetPixelAt(iElements));

		myResultResize->SetPixelAt(iElements, Temp);


	}





	stk::ImageHistogram<TH2F, float> myImageHistogramX;

	myImageHistogramX.SetYAxisTitle("Normalised Counts");
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

	myImageHistogramX.Generate2DHistogram( myResultResize, "NoBu");
	myImageHistogramX.Generate2DHistogram( myResultResize2, "Bu");
	//myImageHistogramX.GenerateXSlice(512);
	myImageHistogramX.SaveHistogram();
	myImageHistogramX.ExportRoot();

stk::HistogramAnalysis analysis;

analysis.LoadRootFile();
analysis.TwoDFit();




}

