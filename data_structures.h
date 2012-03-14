#ifndef _data_structures_included
#define _data_structures_included

#include <boost/numeric/ublas/matrix.hpp>
#include <cstdio>

template <class T>
class array2d : public boost::numeric::ublas::matrix<T>{
	public:
		int cellsize;
		double xllcorner,yllcorner;
		long data_cells;
		T no_data;

		long width() const;
		long height() const ;
		array2d ();
		template<class U> array2d (const array2d<U> &copyfrom, bool do_resize=false);
		long estimated_output_size();
		int print(FILE *fout, int x, int y);
		void init(T val);
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

template <> inline int array2d<float>::print(FILE *fout, int x, int y){
	return fprintf(fout, "%.3f ",boost::numeric::ublas::matrix<float>::operator()(x,y));
}
template <> inline int array2d<char>::print(FILE *fout, int x, int y){
	return fprintf(fout, "%d ",boost::numeric::ublas::matrix<char>::operator()(x,y));
}
template <> inline int array2d<bool>::print(FILE *fout, int x, int y){
	return fprintf(fout, "%d ",boost::numeric::ublas::matrix<bool>::operator()(x,y));
}
template <> inline int array2d<unsigned int>::print(FILE *fout, int x, int y){
	return fprintf(fout, "%d ",boost::numeric::ublas::matrix<unsigned int>::operator()(x,y));
}




typedef array2d<float> float_2d;
typedef array2d<signed char> char_2d;
typedef array2d<bool> bool_2d;
typedef array2d<unsigned int> uint_2d;
typedef array2d<int> int_2d;




typedef struct grid_cell_typez {
	int x;
	int y;
	float z;
	grid_cell_typez(int x0, int y0, float z0){
		x=x0;
		y=y0;
		z=z0;
	}
} grid_cellz;

typedef struct grid_cell_type {
	int x;
	int y;
	grid_cell_type(int x0, int y0){
		x=x0;
		y=y0;
	}
} grid_cell;

#endif
