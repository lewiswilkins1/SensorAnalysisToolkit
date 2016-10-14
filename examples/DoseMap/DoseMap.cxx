/*
 * DoseMap.cxx
 *
 *  Created on: 27 Aug 2015
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

double r(double &x, double &y, double &xs, double &ys);
std::shared_ptr< stk::Image<float> > kFunc(double &x, double &y, double &z, double &alpha, double &beta,double &a, double &omega, double &F,double &A, double &mu, double &nu, std::shared_ptr< stk::Image<float> > &image);



int main(int argc, char*argv[])
{



	std::string pedfilePath,fileNameAndFormat, lightfilePath, lightfileNameAndFormat, filePathandfilename;
	int startingframe, numOfFrames, rows, cols, framesize;
	std::stringstream inputs;
	for(int iArg=1; iArg < argc; iArg++)
	{
		inputs << argv[iArg] << ' ';//get the arguments
	}

	inputs >>  filePathandfilename >>   rows >> cols;//write the inputs
	framesize = cols*rows;




	stk::IO<float> myIO;
	std::shared_ptr< stk::Image<float> > myImage( new stk::Image<float> );
	myImage->Initialise( myIO.ReadImage( filePathandfilename ), rows, cols );


	std::shared_ptr<stk::Image<float> > myResultResize (new stk::Image<float>(1024,1024));

	stk::ImageResize mySize;
	mySize.ResizeImage(myImage , myResultResize);


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
	myImageHistogramX.SaveHistogram();
	myImageHistogramX.ExportRoot();

	stk::HistogramAnalysis analysis;

	analysis.LoadRootFile();
	analysis.TwoDFit(myResultResize);





	typedef itk::Image< float, 2 >  ImageType;
	ImageType::Pointer itkImage = ImageType::New();

	stk::ITKAdapter adapter;
	adapter.STKtoITK(itkImage, myResultResize, 1024, 1024);

	//TFile* file=new TFile("example.root", "RECREATE");

	TFile* file = new TFile("output.root", "RECREATE");
	TH2F* doseMap = new TH2F("Dose Map", "Dose Map",1024, -0.5, 1023.5, 114, 1.5, 30  );
	doseMap->SetDirectory(file);
	double x(0), y(0), xs(0), ys(0), z(0),alpha(0), beta(0), a(0), omega(0), F(0), Az(0), mu(0), nu(0), count(1);
	alpha = 2.8;
	beta = 9.08;
	omega = 1.335;
	F = 0.704;
	//z=1;
	x=0;
	y=0;
	a=0.0809;
	Az= 111.56;
	mu = 5.91158e-02;
	nu = 1.32742e-08;
	for(z=1.5; z<30; z+=0.25){



		std::shared_ptr< stk::Image<float> > k( new stk::Image<float>(5,5));
		std::shared_ptr< stk::Image<float> > kernel;
		kernel = kFunc(x,y,z,alpha,beta,a,Az, mu, nu, omega,F,k);

		ImageType::Pointer itkKernel = ImageType::New();

		ImageType::Pointer itkConvolutionOutput = ImageType::New();
		std::shared_ptr<stk::Image<float> > convolvedImage (new stk::Image<float>(1024,1024));

		adapter.STKtoITK(itkKernel, kernel, 5, 5);

		//std::cout<<"Convolving"<<std::endl;
		typedef itk::ConvolutionImageFilter<ImageType> FilterType;
		FilterType::Pointer convolutionFilter = FilterType::New();


		convolutionFilter->SetInput(itkImage);
		convolutionFilter->SetKernelImage(itkKernel);
		convolutionFilter->Update();
		itkConvolutionOutput = convolutionFilter->GetOutput();
		itkConvolutionOutput->Update();
		itkConvolutionOutput->DisconnectPipeline();



		adapter.ITKtoSTK(itkConvolutionOutput,convolvedImage );
		if(z==2){
			TH2F* map = new TH2F("Map", "Map",1024, -0.5, 1023.5,1024, -0.5, 1023.5 );

			for(int i = 1; i<=1024; i++){
						for(int j=1; j<=1024; j++){
							map->Fill(i,j, convolvedImage->GetPixelAt(((i-1)*1024)+(j-1)));

						}
						}

			map->Draw();
			map->Write();
		}
		TH2F* temp = new TH2F("temp", "temp", 1024, -0.5, 1023.5,1024, -0.5, 1023.5);

		for(int i = 1; i<=1024; i++){
			for(int j=1; j<=1024; j++){

				//std::cout<<i<<" "<<j<<" "<<count<<std::endl;
				temp->Fill(i,j, convolvedImage->GetPixelAt((i-1)*1024+(j-1)));

			}
		}
		std::cout<<"val "<<convolvedImage->GetPixelAt((512)*1024+(512))<<std::endl;
		TH1* slice = temp->ProjectionY("proj", 512);
		for( int i=1; i<=1024;i++){
			doseMap->Fill(i,z, slice->GetBinContent(i));
		}
		std::cout<<slice->GetBinContent(512)<<" z = "<<z<<std::endl;
		delete temp;
		delete slice;
		count++;
		//std::cout<<"Z = "<<count<<" done."<<std::endl;




	}
	doseMap->GetXaxis()->SetTitle( "Row" ) ;
	doseMap->GetYaxis()->SetTitle( "Depth (cm)") ;
	doseMap->SetTitle("Dose Map 5x5 No BU");
	doseMap->Draw("COLZ");
	std::cout<<"Drawn"<<std::endl;
	doseMap->Write();
	std::cout<<"Written to file."<<std::endl;
	file->Close();
	return 1;
}


std::shared_ptr< stk::Image<float> > kFunc(double &x, double &y, double &z, double &alpha, double &beta, double &a, double &omega, double &F,double &A, double &mu, double &nu, std::shared_ptr< stk::Image<float> > &image ){

	std::shared_ptr< stk::Image<float> > k( new stk::Image<float>(6,6));
	double temp(0);
	for(double xs=-2.5; xs<3.5; xs++){
		for(double ys=-2.5; ys<3.5;ys++){
			double kval(0);

			if(r(x,y,xs,ys)==0){
				kval = 1;
			}
			else
			{
				kval = 111.56*std::exp(-5.91158e-02*z*(1-1.32742e-08*z))*(1/(2*M_PI*r(x,y,xs,ys)))*((1-std::exp(-alpha*z))*(alpha*(1-F)*std::exp(-alpha*r(x,y,xs,ys))+beta*F*std::exp(-beta*r(x,y,xs,ys)))+((a*std::pow(z,2))/std::pow(omega*r(x,y,xs,ys)+z,2)));

			}
			//std::cout<<"k = "<<kval<<" z = "<<z<<std::endl;
			k->SetPixelAt((xs+2.5)*6+ys+2.5,kval);
			temp+=kval;
		}
	}
	//std::cout<<temp<<std::endl;
	return k;
}

double r(double &x, double &y, double &xs, double &ys){

	double r;

	r=std::sqrt(std::pow((xs-x),2)+std::pow((ys-y),2));
	return r;
}
