/*
 * stkImage.hxx
 *
 *  Created on: 2 Jun 2015
 *      Author: phrfp
 */

#ifndef __stkImage_hxx
#define __stkImage_hxx
#include "stkImage.h"


namespace stk{

template<typename T_PixelType>//Default constructor
Image<T_PixelType>::Image() : m_row(0), m_col(0){//init values to 0

}

template<typename T_PixelType>
Image<T_PixelType>::Image( const int &row, const int &col ) : m_row(row), m_col(col){
	m_buffer.reset( new std::vector<T_PixelType>( (m_row*m_col), 0) );//create image with element initialised to 0
}

template<typename T_PixelType>
Image<T_PixelType>::~Image(){
	m_buffer.reset();//delete managed object
}

template<typename T_PixelType>
void Image<T_PixelType>::Initialise( const std::shared_ptr< std::vector<T_PixelType> > buffer, const int &row, const int &col ){

	m_buffer = buffer; //take ownership of buffer
	m_row = row;
	m_col = col;
}

template<typename T_PixelType>
void Image<T_PixelType>::Resize(){

	int pix=0;
	for(int iElements=0; iElements<(m_row*m_col-4);iElements+=4)
	{
		T_PixelType tempPixel=0;
		T_PixelType count=0;


		for(int iKernel=iElements; iKernel<(iElements+(4*m_row));iKernel++){


			tempPixel += m_buffer->at(iKernel);
			count++;

			if(count==4)
			{
				iKernel=iKernel+(m_row-4);
				count=0;

			}

		}


		m_buffer->at(pix)=tempPixel/16;
		pix++;
		if((iElements+4)%m_row==0)
		{
			iElements+=(3*m_row);

		}
	}
	m_row = 1024;
	m_col = 1024;

	m_buffer->resize(m_row*m_col);
	m_buffer->shrink_to_fit();

}
template<typename T_PixelType>
void Image<T_PixelType>::Delete(){
	std::cout<<m_buffer->capacity()<<std::endl;
	m_buffer->clear();
	m_buffer->shrink_to_fit();
	std::cout<<m_buffer->capacity()<<std::endl;
}
}


#endif /* __stkImage_hxx */
