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
		array2d ();
		template<class U> 
		array2d (array2d<U> &copyfrom);
};

template <class T>
long array2d<T>::width(){
	return boost::numeric::ublas::matrix<T>::size1();
}

template <class T>
long array2d<T>::height(){
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
array2d<T>::array2d(array2d<U> &copyfrom){
	cellsize=copyfrom.cellsize;
	xllcorner=copyfrom.xllcorner;
	yllcorner=copyfrom.yllcorner;
	data_cells=copyfrom.data_cells;
	no_data=copyfrom.no_data;
}

typedef array2d<float> float_2d;
typedef array2d<char> char_2d;
typedef array2d<bool> bool_2d;
typedef array2d<unsigned int> uint_2d;

#endif
