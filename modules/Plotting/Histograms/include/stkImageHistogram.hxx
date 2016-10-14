/*
 * stkImageHistogram.hxx
 *
 *  Created on: 10 Jun 2015
 *      Author: phrfp
 */

#ifndef __stkImageHistogram_hxx
#define __stkImageHistogram_hxx

#include "stkImageHistogram.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TFile.h"
namespace stk{

template<class T_Hist, typename T_Pixeltype>
ImageHistogram<T_Hist, T_Pixeltype>::ImageHistogram() {

	m_HistSetup.mainTitle = "";
	m_HistSetup.ouputfileNameAndPath ="";
	m_HistSetup.histFileType = FileType::PNG;
	m_HistSetup.statOptions= 1111;
	m_HistSetup.x_axisTitle="";
	m_HistSetup.numberOfXBins=0;
	m_HistSetup.xstart=0;
	m_HistSetup.xend=0;
	m_HistSetup.gridX =false;
	m_HistSetup.numberOfYBins=0;
	m_HistSetup.ystart=0;
	m_HistSetup.yend=0;
	m_HistSetup.y_axisTitle ="";
	m_HistSetup.logY=false;
	m_HistSetup.gridY=false;
	slice = false;








}



template<class T_Hist, typename T_Pixeltype>
ImageHistogram<T_Hist, T_Pixeltype>::~ImageHistogram(){


	m_histogram.reset();



}

template<class T_Hist, typename T_Pixeltype>
void ImageHistogram<T_Hist, T_Pixeltype>::GenerateHistogram( const std::shared_ptr< stk::Image<T_Pixeltype> > imageToHistogram,const std::string &histName){

	m_histogram.reset( new T_Hist( m_HistSetup.mainTitle.c_str() , m_HistSetup.mainTitle.c_str(),
			m_HistSetup.numberOfXBins, m_HistSetup.xstart, m_HistSetup.xend) );

	typename std::vector<T_Pixeltype>::iterator imageIt;
	for( imageIt = imageToHistogram->StartImage(); imageIt != imageToHistogram->EndImage(); imageIt++ )
		m_histogram->Fill( (*imageIt) );

}

template< class T_Hist, typename T_Pixeltype >
void ImageHistogram<T_Hist, T_Pixeltype>::Generate2DHistogram( const std::shared_ptr< stk::Image<T_Pixeltype> >  imageToHistogram,const std::string histName ){

	std::cout<<imageToHistogram->NumberOfRows()<<std::endl;
	m_histogram.reset( new T_Hist( m_HistSetup.mainTitle.c_str() , m_HistSetup.mainTitle.c_str(),
			m_HistSetup.numberOfXBins, m_HistSetup.xstart, m_HistSetup.xend,
			m_HistSetup.numberOfYBins, m_HistSetup.ystart, m_HistSetup.yend) );

	for( int irows = 0; irows < imageToHistogram->NumberOfRows(); irows++ ){
		for( int icols = 0; icols < imageToHistogram->NumberOfColumns(); icols++ ){
			m_histogram->Fill( icols, irows, imageToHistogram->GetPixelAt(irows,icols) );
			m_histogram->SetName(histName.c_str());


		}

	}
	TH2Container.push_back(m_histogram);


}








template<class T_Hist, typename T_Pixeltype>
void ImageHistogram<T_Hist, T_Pixeltype>::GenerateXSlice(const int &bin){

	for(int iHist=0; iHist<TH2Container.size();iHist++){


		TH1Container.push_back(TH2Container.at(iHist)->ProjectionX(std::to_string(iHist).c_str(), bin,bin));


	}
	slice = true;
	std::string title = "X Slice";
	m_HistSetup.mainTitle = title;





}

template<class T_Hist, typename T_Pixeltype>
void ImageHistogram<T_Hist, T_Pixeltype>::GenerateYSlice(const int &bin){


	for(int iHist=0; iHist<TH2Container.size();iHist++){


		TH1Container.push_back(TH2Container.at(iHist)->ProjectionY(std::to_string(iHist).c_str(), bin,bin));


	}
	slice =true;
	m_HistSetup.mainTitle = "Y Slice";


}


template<class T_Hist, typename T_Pixeltype>
void ImageHistogram<T_Hist, T_Pixeltype>::GenerateXSlice(const int &bin, const int &end){


	for(int iHist=0; iHist<TH2Container.size();iHist++){


		TH1Container.push_back(TH2Container.at(iHist)->ProjectionX(std::to_string(iHist).c_str(), bin,end));


	}
	slice = true;
	std::string title = "X Slice";
	m_HistSetup.mainTitle = title;
	std::cout<<m_HistSetup.mainTitle<<std::endl;
	m_HistSetup.statOptions= 1111111;



}

template<class T_Hist, typename T_Pixeltype>
void ImageHistogram<T_Hist, T_Pixeltype>::GenerateYSlice(const int &bin, const int &end){


	for(int iHist=0; iHist<TH2Container.size();iHist++){


		TH1Container.push_back(TH2Container.at(iHist)->ProjectionY(std::to_string(iHist).c_str(), bin,end));


	}
	slice =true;
	m_HistSetup.mainTitle = "aaa";


}


template<class T_Hist, typename T_Pixeltype>
void ImageHistogram<T_Hist, T_Pixeltype>::SaveHistogram(){


	std::shared_ptr<TCanvas> histCanvas( new TCanvas() );
	histCanvas->SetGridx( m_HistSetup.gridX );
	histCanvas->SetGridy( m_HistSetup.gridY );
	histCanvas->SetLogy( m_HistSetup.logY );
	gStyle->SetOptStat( m_HistSetup.statOptions );



	if(slice==true){
		TH1Container.at(0)->Draw();

	}

	for(int i=0; i<TH2Container.size();i++){
	TH2Container.at(i)->GetXaxis()->SetTitle( m_HistSetup.x_axisTitle.c_str() ) ;
	TH2Container.at(i)->GetYaxis()->SetTitle( m_HistSetup.y_axisTitle.c_str() ) ;

	TH2Container.at(i)->Draw("COLZ");
	//TH2Container.at(0)->Write();
	}

	std::shared_ptr<TImage> img(TImage::Create());
	img->FromPad( histCanvas.get() );





	std::string tempName=m_HistSetup.ouputfileNameAndPath;

	switch(  m_HistSetup.histFileType()  ){
	case FileType::JPG:
		tempName += ".jpg";
		img->SetImageQuality (TAttImage::EImageQuality::kImgBest);
		img->WriteImage( tempName.c_str());
		break;
	case FileType::PNG:
		tempName += ".png";
		img->WriteImage( tempName.c_str());
		break;
	default:
		break;
	}
	histCanvas.reset();

	//histCanvas->Clear();
	img.reset();



}

template<class T_Hist, typename T_Pixeltype>
void ImageHistogram<T_Hist, T_Pixeltype>::ExportRoot(std::string outputName){
	
	outputName+=".root";
	const char * filename = outputName.c_str();
	TFile *myFile =  new TFile(filename,"RECREATE");

	for(int i=0; i<TH2Container.size(); i++){
		TH2Container.at(i)->Write();
		std::cout<<"Writing"<<std::endl;
	}
	for(int i=0; i<TH1Container.size(); i++){
			TH1Container.at(i)->Write();
			std::cout<<"Writing"<<std::endl;
		}

	myFile->Close();

}


}


#endif /* MODULES_PLOTTING_HISTOGRAMS_INLCUDE_STKIMAGEHISTOGRAM_HXX_ */
