/*
 * GammaMap.cxx
 *
 *  Created on: 24 Aug 2015
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
double r(double &x, double &y, double &xs, double &ys);
std::shared_ptr< stk::Image<float> > kFunc(double &x, double &y, double &z, double &alpha, double &beta,double &a, double &omega, double &F, std::shared_ptr< stk::Image<float> > &image);

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
	std::shared_ptr<stk::Image<float> > convolvedNoBu (new stk::Image<float>(1024,1024));
	std::shared_ptr<stk::Image<float> > convolvedBu (new stk::Image<float>(1024,1024));
	stk::ImageResize mySize;

	//Image resizing
	mySize.ResizeImage(myImage , myResultResize);
	mySize.ResizeImage(myImage2 , myResultResize2);



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

	myImageHistogramX.Generate2DHistogram( myResultResize, "NoBu");
	myImageHistogramX.Generate2DHistogram( myResultResize2, "Bu");
	//myImageHistogramX.GenerateXSlice(512);
	myImageHistogramX.SaveHistogram();
	myImageHistogramX.ExportRoot();

	stk::HistogramAnalysis analysis;

	analysis.LoadRootFile();
	analysis.TwoDFit(myResultResize, myResultResize2);

	stk::ImageHistogram<TH2F, float> myImageHistogramY;
	myImageHistogramY.SetOutputFileNameAndPath("X");
	myImageHistogramY.SetOutputFileType( stk::FileType::PNG );
	myImageHistogramY.Generate2DHistogram( myResultResize2, "Bu");
	myImageHistogramY.SaveHistogram();



	typedef itk::Image< float, 2 >  ImageType;

	ImageType::Pointer itkImage = ImageType::New();
	ImageType::Pointer itkImage2 = ImageType::New();
	ImageType::Pointer itkBlurredImage = ImageType::New();
	ImageType::Pointer itkBlurredImage2 = ImageType::New();
	ImageType::Pointer itkTranslated = ImageType::New();

	ImageType::RegionType region;
	ImageType::IndexType start;


	ImageType::SizeType size;
	start[0] = 0;
	start[1] = 0;
	size[0] = rows;
	size[1] = cols;
	region.SetSize(size);
	region.SetIndex(start);
	itkBlurredImage->SetRegions(region);
	itkBlurredImage->Allocate();
	itkBlurredImage2->SetRegions(region);
	itkBlurredImage2->Allocate();
	itkTranslated->SetRegions(region);
	itkTranslated->Allocate();



	stk::ITKAdapter adapter;
	adapter.STKtoITK(itkImage, myResultResize, 1024, 1024);
	adapter.STKtoITK(itkImage2, myResultResize2, 1024, 1024);
	stk::ITKTools gaussian;

	itkTranslated = gaussian.ImageRegistration(itkImage2, itkImage);

	itkTranslated->Update();



	double x(0), y(0), xs(0), ys(0), z(0),alpha(0), beta(0), a(0), omega(0), F(0);
	alpha = 2.8;
	beta = 9.08;
	omega = 1.335;
	F = 0.704;
	z=1;
	x=0;
	y=0;
	a=0.0809;
	std::shared_ptr< stk::Image<float> > k( new stk::Image<float>(5,5));
	std::shared_ptr< stk::Image<float> > kernel;
	kernel = kFunc(x,y,z,alpha,beta,a,omega,F,k);

	ImageType::Pointer itkKernel = ImageType::New();

	ImageType::Pointer itkConvolutionOutput = ImageType::New();
	ImageType::Pointer itkConvolutionOutput2 = ImageType::New();

	adapter.STKtoITK(itkKernel, kernel, 5, 5);

	std::cout<<"Convolving"<<std::endl;
	typedef itk::ConvolutionImageFilter<ImageType> FilterType;
	FilterType::Pointer convolutionFilter = FilterType::New();
	FilterType::Pointer convolutionFilter2 = FilterType::New();

	convolutionFilter->SetInput(itkImage);
	convolutionFilter->SetKernelImage(itkKernel);
	convolutionFilter->Update();
	itkConvolutionOutput = convolutionFilter->GetOutput();
	itkConvolutionOutput->Update();


	convolutionFilter2->SetInput(itkTranslated);
	convolutionFilter2->SetKernelImage(itkKernel);
	convolutionFilter2->Update();
	itkConvolutionOutput2 = convolutionFilter2->GetOutput();
	itkConvolutionOutput2->Update();

	adapter.ITKtoSTK(itkConvolutionOutput,convolvedNoBu );
	adapter.ITKtoSTK(itkConvolutionOutput2,convolvedBu );


	std::shared_ptr< stk::Image< float > > gammaMap( new stk::Image<float>(1024,1024));

std::cout<<"starting gamma map"<<std::endl;
	double disM(0), dosM(0), pass(0);
	disM= 30;
	dosM=0.03;
	for(int x=20; x<1014; x++){
		for( int y=20; y<1014; y++){

			double gamma=10;
			for(int xs=-20; xs<20; xs++){
				for( int ys=-20; ys<20;ys++){


//
					double temp = std::sqrt( (std::pow(xs,2)+std::pow(ys,2))/std::pow(disM,2) + std::pow( convolvedBu->GetPixelAt((x+xs)*1024+(y+ys)) - convolvedNoBu->GetPixelAt(x*1024+y), 2 )/std::pow(dosM, 2));
					if(temp==0){

						std::cout<<"x: "<<x<<" y: "<<y<<" xs: "<<xs<<" ys: "<<ys<<" Dc: "<<convolvedBu->GetPixelAt((x+xs)*1024+(y+ys))<<" Dm: "<<convolvedNoBu->GetPixelAt(x*1024+y)<<std::endl;

					}
					if(temp<gamma){
						gamma=temp;
					}
				}
			}
			//std::cout<<gamma<<std::endl;
			gammaMap->SetPixelAt(x*1024+y, gamma);
			if(gamma<=1){
				pass++;
			}


		}
	}


std::cout<<"Percentage pass = "<<(pass/1024*1024)*100<<"%"<<std::endl;

	stk::ImageHistogram<TH2F, float> myImageHistogram;

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
		myImageHistogram.SetOutputFileNameAndPath("gammaMap");
		myImageHistogram.SetOutputFileType( stk::FileType::PNG );

		myImageHistogram.Generate2DHistogram( gammaMap, "GammaMap");
//		myImageHistogram.Generate2DHistogram( convolvedNoBu, "convolvedNoBu");
//		myImageHistogram.Generate2DHistogram( convolvedBu, "convolvedNoBuv");
		myImageHistogram.SaveHistogram();
		//myImageHistogramX.ExportRoot();
	return 1;
}

double r(double &x, double &y, double &xs, double &ys){

	double r;

	r=std::sqrt(std::pow((xs-x),2)+std::pow((ys-y),2));
	return r;
}


std::shared_ptr< stk::Image<float> > kFunc(double &x, double &y, double &z, double &alpha, double &beta, double &a, double &omega, double &F,std::shared_ptr< stk::Image<float> > &image ){

	std::shared_ptr< stk::Image<float> > k( new stk::Image<float>(6,6));
	double temp(0);
	for(double xs=-2.5; xs<3.5; xs++){
		for(double ys=-2.5; ys<3.5;ys++){
			double kval;

			if(r(x,y,xs,ys)==0){
				kval = 1;
			}
			else
			{
				kval = (1/(2*M_PI*r(x,y,xs,ys)))*((1-std::exp(-alpha*z))*(alpha*(1-F)*std::exp(-alpha*r(x,y,xs,ys))+beta*F*std::exp(-beta*r(x,y,xs,ys)))+((a*std::pow(z,2))/std::pow(omega*r(x,y,xs,ys)+z,2)));

			}
			//std::cout<<"here"<<std::endl;
			k->SetPixelAt((xs+2.5)*6+ys+2.5,kval);
			temp+=kval;
		}
	}
	std::cout<<temp<<std::endl;
	return k;
}
