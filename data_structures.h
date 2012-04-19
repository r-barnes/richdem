#ifndef _data_structures_included
#define _data_structures_included

#define NDEBUG
#define BOOST_UBLAS_NDEBUG	//TODO: Not sure if this does anything
#include <boost/numeric/ublas/matrix.hpp>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>

template <class T>
class array2d : public boost::numeric::ublas::matrix<T>{
	public:
		int cellsize;
//		std::string xllcorner,yllcorner;	//TODO: Should be string
		double xllcorner,yllcorner;
		long data_cells;
		T no_data;

		long width() const;
		long height() const ;
		array2d ();
		template<class U> array2d (const array2d<U> &copyfrom, bool do_resize=false);
		long estimated_output_size();
		void init(T val);
		void print_block(std:: ostream& out, int minx, int maxx, int miny, int maxy, int precision=0, std::streamsize swidth=2); //TODO
		bool operator==(const array2d<T> &other) const;
		template<class U> friend std::ostream& operator<<(std::ostream &out, const array2d<U> &arr);
		T max() const;
		T min() const;
//		double avg() const; //This should use the Kaham summation algorithm with a running-average code, or something like that.
};

template <class T>
long array2d<T>::width() const {
	return boost::numeric::ublas::matrix<T>::size1();
}

template <class T>
long array2d<T>::height() const {
	return boost::numeric::ublas::matrix<T>::size2();
}

template <class T>
array2d<T>::array2d(){
	cellsize=-1;
	xllcorner=-1;
	yllcorner=-1;
	data_cells=-1;
	no_data=-1;
}

template <class T>
template <class U>
array2d<T>::array2d(const array2d<U> &copyfrom, bool do_resize){
	cellsize=copyfrom.cellsize;
	xllcorner=copyfrom.xllcorner;
	yllcorner=copyfrom.yllcorner;
	data_cells=copyfrom.data_cells;
	no_data=copyfrom.no_data;
	if(do_resize)
		boost::numeric::ublas::matrix<T>::resize(copyfrom.width(),copyfrom.height());
}

template <class T>
void array2d<T>::init(T val){
	#pragma omp parallel for
	for(int x=0;x<width();x++)
	for(int y=0;y<height();y++)
		boost::numeric::ublas::matrix<T>::operator()(x,y)=val;
}

template <> inline long array2d<float>::estimated_output_size(){return 9*this->width()*this->height();}
template <> inline long array2d<char>::estimated_output_size(){return 4*this->width()*this->height();}
template <> inline long array2d<bool>::estimated_output_size(){return 2*this->width()*this->height();}
template <> inline long array2d<unsigned int>::estimated_output_size(){return 9*this->width()*this->height();}

template <class T>
std::ostream& operator<< (std::ostream &out, const array2d<T> &arr){
	std::streamsize width=out.width();
	for(int y=0;y<arr.height();y++){
		for(int x=0;x<arr.width();x++)
			if(sizeof(T)==1)	//TODO: An ugly way of detecting chars
				out<<std::setw(width)<<(int)arr(x,y)<<" ";
			else
				out<<std::fixed<<std::setw(width)<<arr(x,y)<<" ";
		out<<std::endl;
	}
	return out;
}

template <class T>
bool array2d<T>::operator==(const array2d<T> &other) const {
	if(width()!=other.width() || height()!=other.height())
		return false; 
	for(int x=0;x<width();x++)
	for(int y=0;y<height();y++)
		if(boost::numeric::ublas::matrix<T>::operator()(x,y)!=other(x,y))
			return false;
	return true;
}

template <class T>
T array2d<T>::max() const {
	bool init=false;
	T maxval=no_data;
	T temp;
	//TODO: OpenMP 3.1 min/max reduction operators can speed this up
	for(int x=0;x<width();x++)
	for(int y=0;y<height();y++){
		temp=boost::numeric::ublas::matrix<T>::operator()(x,y);
		if(temp==no_data) continue;
		if(init && temp>maxval)
			maxval=temp;
		else if (!init){
			maxval=temp;
			init=true;
		}
	}
	return maxval;
}

template <class T>
T array2d<T>::min() const {
	bool init=false;
	T minval=no_data;
	T temp;
	//TODO: OpenMP 3.1 min/max reduction operators can speed this up
	for(int x=0;x<width();x++)
	for(int y=0;y<height();y++){
		temp=boost::numeric::ublas::matrix<T>::operator()(x,y);
		if(temp==no_data) continue;
		if(init && temp<minval)
			minval=temp;
		else if (!init){
			minval=temp;
			init=true;
		}
	}
	return minval;
}

typedef array2d<double> double_2d;
typedef array2d<float> float_2d;
typedef array2d<signed char> char_2d;
typedef array2d<bool> bool_2d;
typedef array2d<unsigned int> uint_2d;
typedef array2d<int> int_2d;

typedef struct grid_cell_typez {
	int x;
	int y;
	float z;
	grid_cell_typez(int x0, int y0, float z0):x(x0),y(y0),z(z0){}
	grid_cell_typez(){}
} grid_cellz;

typedef struct grid_cell_type {
	int x;
	int y;
	grid_cell_type(int x0, int y0){
		x=x0;
		y=y0;
	}
	grid_cell_type(){}
} grid_cell;

#endif
