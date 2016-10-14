/*
 * stkITKAdapter.hxx
 *
 *  Created on: 4 Aug 2015
 *      Author: lewish
 */

#ifndef __stkITKAdapter_hxx
#define __stkITKAdapter_hxx

#include "stkITKAdapter.h"
#include "itkImage.h"


namespace stk{


ITKAdapter::ITKAdapter(){}
ITKAdapter::~ITKAdapter(){}



template<typename T_Pixeltype>
void ITKAdapter::STKtoITK(itk::Image< float, 2 >::Pointer itkOutput, const std::shared_ptr< stk::Image<T_Pixeltype> > inputImage, const int &rows, const int &cols){


	typedef itk::Image<float,2> ImageType;
	ImageType::RegionType region;
	ImageType::IndexType start;


	ImageType::SizeType size;
	start[0] = 0;
	start[1] = 0;
	size[0] = rows;
	size[1] = cols;
	region.SetSize(size);
	region.SetIndex(start);
	itkOutput->SetRegions(region);
	itkOutput->Allocate();



	typedef itk::ImageRegionIterator< ImageType > IteratorType;
	IteratorType out( itkOutput, region );
	typename std::vector<T_Pixeltype>::iterator inputStart;
	typename std::vector<T_Pixeltype>::iterator inputEnd;


	for( inputStart=inputImage->StartImage(), out.GoToBegin(); !out.IsAtEnd(); ++inputStart, ++out){
		out.Set(*inputStart);

	}

}


template<typename T_Pixeltype>
void ITKAdapter::ITKtoSTK(itk::Image< float, 2 >::Pointer itkInput, const std::shared_ptr< stk::Image<T_Pixeltype> > outputImage){

	typedef itk::Image<float,2> ImageType;

	typedef itk::ImageRegionIterator<ImageType > IteratorType;
		IteratorType in( itkInput, itkInput->GetLargestPossibleRegion() );
		typename std::vector<T_Pixeltype>::iterator outputStart;
		typename std::vector<T_Pixeltype>::iterator outputEnd;


		for( outputStart=outputImage->StartImage(), in.GoToBegin(); !in.IsAtEnd(); ++outputStart, ++in){

			*outputStart = in.Get();

		}



}

}

#endif /* __stkITKAdapter_hxx */
