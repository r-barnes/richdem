#ifndef _data_structures_included
#define _data_structures_included

#include <boost/numeric/ublas/matrix.hpp>

template <class T>
class array2d : public boost::numeric::ublas::matrix<T>{
	public:
		int cellsize;
		double xllcorner,yllcorner;
		long data_cells;
		T no_data;

		long width();
		long height();
};

template <class T>
long array2d<T>::width(){
	return boost::numeric::ublas::matrix<T>::size1();
}

template <class T>
long array2d<T>::height(){
	return boost::numeric::ublas::matrix<T>::size2();
}

typedef array2d<float> float_2d;
typedef array2d<char> char_2d;
typedef array2d<bool> bool_2d;
typedef array2d<unsigned int> uint_2d;

#endif
