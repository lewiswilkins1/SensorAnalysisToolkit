/*
 * BasicImageProcessing.cxx
 *
 *  Created on: 28 Jul 2015
 *      Author: lewish
 */


/*
 * BasicImageProcessing.cxx
 *
 *  Created on: 28 Jul 2015
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
#include "stkDarkFrameDetection.h"


int main(int argc, const char** argv){



	if(argc != 18 )
	{
		std::cout<<"Usage: BasicImageProcessing [PEDFILEPATH] [PEDFILENAMEANDFORMAT] [LIGHTFILEPATH] [LIGHTFILENAMEANDFORMAT] [RAWFILEPATH] [RAWFILENAMEANDFORMAT] [OUTPUTFILENAMEANDPATH]  [ROWS] [COLS] [PEDSTARTINGFRAME] [PEDNUMBEROFFILES] [LIGHTSTARTINGFRAME] [LIGHTNUMBEROFFILES] [RAWSTARTINGFRAME] [RAWNUMBEROFFILES]"<<std::endl;
		std::cout<<"Only 16 bit is supported in this example"<<std::endl;
		return 0;
	}



	//Variables to hold the input parameters

	std::string pedfilePath, pedfileNameAndFormat, lightfilePath, lightfileNameAndFormat, filePath, outFileNameAndPath, fileNameAndFormat, outFilePath, outFileName, outFileExtension;
	int startingframe, numOfFrames, rows, cols, framesize, pedStartingframe, pedNumOfFrames, lightStartingframe, lightNumOfFrames;


	std::stringstream inputs;
	for(int iArg=1; iArg < argc; iArg++)
	{
		inputs << argv[iArg] << ' ';//get the arguments
	}

	inputs >> pedfilePath >> pedfileNameAndFormat  >> lightfilePath >> lightfileNameAndFormat >> filePath >> fileNameAndFormat >> outFilePath >> outFileName >> outFileExtension >>  rows >> cols >> pedStartingframe >> pedNumOfFrames >> lightStartingframe >> lightNumOfFrames >> startingframe >> numOfFrames;//write the inputs
	framesize = cols*rows;
	//outFileNameAndPath = outFilePath+outFileName+outFileExtension;

	std::cout<<"Input Works"<<std::endl;
	outFileNameAndPath = outFilePath+outFileName+outFileExtension;
	std::cout<<outFileNameAndPath<<std::endl;
	//Load ped image stack
	stk::IOImageStack<unsigned short> mypedIO;
	std::shared_ptr< stk::ImageStack<unsigned short> > myPedImageStack( new stk::ImageStack<unsigned short> );
	myPedImageStack->Initialise( mypedIO.ReadImageStack( pedfilePath, pedfileNameAndFormat, pedStartingframe, pedNumOfFrames, framesize ), pedNumOfFrames, rows, cols );
	std::cout<<"ped Works"<<std::endl;
	//Load light image stack
	stk::IOImageStack<unsigned short> mylightIO;
	std::shared_ptr< stk::ImageStack<unsigned short> > mylightStack( new stk::ImageStack<unsigned short> );
	mylightStack->Initialise( mylightIO.ReadImageStack( lightfilePath, lightfileNameAndFormat, lightStartingframe, lightNumOfFrames, framesize ), lightNumOfFrames, rows, cols );
	std::cout<<"light Works"<<std::endl;
	//Load Image Stack
	//	stk::IOImageStack<unsigned short> myIO;
	//	std::shared_ptr< stk::ImageStack<unsigned short> > myImageStack( new stk::ImageStack<unsigned short> );
	//	myImageStack->Initialise( myIO.ReadImageStack( filePath, fileNameAndFormat, startingframe, numOfFrames, framesize ), numOfFrames, rows, cols );


	std::cout<<"Images loaded."<<std::endl;



	/*
	 * Declaring things
	 *
	 */

	//Sum
	stk::ImageSum mySummer;

	//divide
	stk::ImageDivision myDivider;

	//variance
	stk::ImageVariance myVarianceCalc;

	//mask
	stk::ImageMask myMask;

	//Pixel average
	stk::ImageBadPixelAverage myAverage;

	//Ped Remover
	stk::ImageMinus myMinus;

	//resizer
	stk::ImageResize mySize;

	//FieldCorrecter
	stk::ImageFlatFieldCorrection myField;

	stk::DarkFrameDetection darkFrameDetect;


	/*
	 * Image initialising
	 */

	//Image to hold pedestal
	std::shared_ptr< stk::Image<float> > myPedestal ( new stk::Image<float>(4096,4096) );

	//Image to hold variance
	std::shared_ptr<stk::Image<float> > myVarianceImage ( new stk::Image<float>(4096,4096) );

	//Image to hold mask
	std::shared_ptr<stk::Image<bool> > myMaskImage ( new stk::Image<bool>(4096,4096) );

	//Image to hold result
	std::shared_ptr<stk::Image<float> > myResult (new stk::Image<float>(4096,4096));

	//Image to hold gain
	std::shared_ptr<stk::Image<float> > myGain (new stk::Image<float>(4096,4096));

	//Image to hold resized result
	std::shared_ptr<stk::Image<float> > myResultResize (new stk::Image<float>(1024,1024));

	//Image to hold resized result
	std::shared_ptr<stk::Image<float> > myLightImage (new stk::Image<float>(1024,1024));






	/*
	 * Begin arithmetic
	 */

	//Creating pedestal image.
	mySummer.SumImageStack( myPedImageStack, myPedestal );
	myDivider.DivideImage(myPedestal, static_cast<float>(myPedImageStack->NumberOfImageInStack()) );

	std::cout<<"Pedestal Calculated."<<std::endl;
	float pedAvg(0);
	for (int iElements = 0; iElements<framesize; iElements++)
	{

		pedAvg += myPedestal->GetPixelAt(iElements);
	}
	pedAvg=pedAvg/framesize;


	float pedVar(0);
	for (int iElements=0; iElements<framesize; iElements++)
	{
		pedVar+= (myPedestal->GetPixelAt(iElements)-pedAvg)*(myPedestal->GetPixelAt(iElements)-pedAvg);
	}
	pedVar=std::sqrt(pedVar/framesize);



	darkFrameDetect.DetectDarkFrames(myResult ,myPedestal, pedAvg,filePath, fileNameAndFormat, rows, cols, numOfFrames, startingframe  );

	//std::cout<<myResult->GetPixelAt(9454)<<std::endl;
	//Variance calculation

	std::cout<<"Ped of Pixel 0: "<<myPedestal->GetPixelAt(0)<<std::endl;
	std::cout<<"Var of Pixel 0: "<<myVarianceImage->GetPixelAt(0)<<std::endl;
	std::cout<<"Ped subtracted stack pixel 0: "<<myPedImageStack->GetPixelAt(0)<<std::endl;

	myVarianceCalc.VarianceImageStack(myPedImageStack, myPedestal, myVarianceImage);

	std::cout<<"Var post calc of Pixel 0: "<<myVarianceImage->GetPixelAt(0)<<std::endl;


	float varAvg(0), varVar(0);
	for (int iElements = 0; iElements<framesize; iElements++)
	{
		varAvg += myVarianceImage->GetPixelAt(iElements);
	}
	varAvg=(varAvg/(float)framesize);
	std::cout<<"Ave Var: "<<varAvg<<std::endl;

	for (int iElements=0; iElements<framesize; iElements++)
	{
		varVar+= (myVarianceImage->GetPixelAt(iElements)-varAvg)*(myVarianceImage->GetPixelAt(iElements)-varAvg);
	}
	varVar=std::sqrt(varVar/(float)framesize);

	std::cout<<"Variance calculated."<<std::endl;

	std::cout<<"Images integrated."<<std::endl;



	//	int framesUsed(0);
	//	darkFrameDetect.DetectDarkFrames(myLightImage, pedAvg, lightfilePath, lightfileNameAndFormat, rows, cols, numOfFrames, lightStartingframe, framesUsed );
	//	myDivider.DivideImage(myLightImage, static_cast<float>(framesUsed));
	//
	//
	//
	//	//Flat field correction
	//	myField.CorrectImage(myLightImage, myPedestal, myResult, myGain);

	std::cout<<"Cut values Ped: "<<pedAvg<<" Var: "<<varVar<<std::endl;

	//Mask generation
	myMask.MaskImage(myResult, myVarianceImage, myMaskImage , myPedestal,myGain, varAvg+3*varVar, varAvg-2*varVar, pedAvg+2*pedVar, pedAvg-2*pedVar);

	//Averaging the bad pixels
	myAverage.BadPixelAverage(myResult, myMaskImage);

	std::cout<<"Bad pixels removed."<<std::endl;



	//Image resizing
	mySize.ResizeImage(myResult , myResultResize);

	stk::IO<float> imageIO;
	imageIO.WriteImage( myResult, outFileNameAndPath );






	stk::ImageHistogram<TH2F, float> myImageHistogram;
	myImageHistogram.SetTitle(fileNameAndFormat);
	myImageHistogram.SetYAxisTitle("Row");
	myImageHistogram.SetYAxisLog(false);
	myImageHistogram.SetNumberOfYBins(1024);
	myImageHistogram.SetYLimits( -0.5, 1023.5 );
	myImageHistogram.SetXAxisTitle( "Col" );
	myImageHistogram.SetNumberOfXBins( 1024);
	myImageHistogram.SetXLimits(  -0.5, 1023.5  );

	myImageHistogram.SetGridY(false);
	myImageHistogram.SetGridX(false);

	myImageHistogram.SetStatBoxOptions(0);
	myImageHistogram.SetOutputFileNameAndPath(outFilePath+outFileName);
	myImageHistogram.SetOutputFileType( stk::FileType::PNG );

	myImageHistogram.Generate2DHistogram( myResultResize, "Map");
	myImageHistogram.SaveHistogram();




	stk::ImageHistogram<TH2F, float> myImageHistogramX;

	myImageHistogramX.SetYAxisTitle("ADC");
	myImageHistogramX.SetYAxisLog(false);
	myImageHistogramX.SetNumberOfYBins(1024);
	myImageHistogramX.SetYLimits( -0.5, 1023.5 );
	myImageHistogramX.SetXAxisTitle( "Row" );
	myImageHistogramX.SetNumberOfXBins( 1024);
	myImageHistogramX.SetXLimits(  -0.5, 1023.5  );

	myImageHistogramX.SetGridY(false);
	myImageHistogramX.SetGridX(false);

	myImageHistogramX.SetStatBoxOptions(0);



	std::string xEnd = "X";
	std::string histTitleX = outFileName+"_"+xEnd;
	outFilePath.append(outFileName+"_");
	outFilePath.append(xEnd);

	myImageHistogramX.SetOutputFileNameAndPath(outFilePath);
	myImageHistogramX.SetTitle(histTitleX);
	outFilePath.pop_back();
	myImageHistogramX.SetOutputFileType( stk::FileType::PNG );

	myImageHistogramX.Generate2DHistogram( myResultResize, "X Slice");
	myImageHistogramX.GenerateXSlice(512);
	myImageHistogramX.SaveHistogram();



	stk::ImageHistogram<TH2F, float> myImageHistogramY;

	myImageHistogramY.SetYAxisTitle("ADC");
	myImageHistogramY.SetYAxisLog(false);
	myImageHistogramY.SetNumberOfYBins(1024);
	myImageHistogramY.SetYLimits( -0.5, 1023.5 );
	myImageHistogramY.SetXAxisTitle( "Col" );
	myImageHistogramY.SetNumberOfXBins( 1024);
	myImageHistogramY.SetXLimits(  -0.5, 1023.5  );

	myImageHistogramY.SetGridY(false);
	myImageHistogramY.SetGridX(false);

	std::string yEnd = "Y";
	myImageHistogramY.SetStatBoxOptions(0);
	std::string histTitleY = outFileName+"_"+yEnd;


	myImageHistogramY.SetTitle(histTitleY);
	myImageHistogramY.SetOutputFileNameAndPath(outFilePath.append(yEnd));
	myImageHistogramY.SetOutputFileType( stk::FileType::PNG );

	myImageHistogramY.Generate2DHistogram( myResultResize, "Y Slice");
	myImageHistogramY.GenerateYSlice(512);
	myImageHistogramY.SaveHistogram();


	return 1;
}
