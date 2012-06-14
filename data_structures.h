#ifndef _data_structures_included
#define _data_structures_included

#define NDEBUG
#define BOOST_UBLAS_NDEBUG	//TODO: Not sure if this does anything
#include <boost/numeric/ublas/matrix.hpp>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>

//Use this D8 method //TODO
//321
//4 0
//567
const int dinf_dx[9]={1,1,0,-1,-1,-1,0,1,1};
const int dinf_dy[9]={0,-1,-1,-1,0,1,1,1,0};
const int dinf_d8_inverse[9]={4,5,6,7,0,1,2,3,4};

template <class T>
class array2d : public boost::numeric::ublas::matrix<T>{
	public:
		double cellsize;
//		std::string xllcorner,yllcorner;	//TODO: Should be string
		double xllcorner,yllcorner;
		long data_cells;
		T no_data;

		array2d ();
		template<class U> array2d (const array2d<U> &copyfrom);
		long width() const
			{return boost::numeric::ublas::matrix<T>::size1();}
		long height() const
			{return boost::numeric::ublas::matrix<T>::size2();}
		template<class U> void copyprops (const array2d<U> &copyfrom);
		long estimated_output_size();
		void init(T val);
		void print_block(std:: ostream& out, int minx, int maxx, int miny, int maxy, int precision=0, std::streamsize swidth=2); //TODO
		void surroundings(int x0, int y0, int precision=3) const; //TODO
		bool operator==(const array2d<T> &other) const;
		template<class U> friend std::ostream& operator<<(std::ostream &out, const array2d<U> &arr);
		T max() const;
		T min() const;
		inline bool in_grid(int x, int y) const
			{return (x>=0 && y>=0 && x<width() && y<height());}
		inline bool interior_grid(int x, int y) const
			{return (x>=1 && y>=1 && x<width()-1 && y<height()-1);}
		inline bool edge_grid(int x, int y) const
			{return (x==0 || y==0 || x==width()-1 || y==height()-1);}
		T& operator()(int x, int y)
			{return boost::numeric::ublas::matrix<T>::operator()(x,y);}
		const T& operator()(int x, int y) const
			{return boost::numeric::ublas::matrix<T>::operator()(x,y);}
		void resize(int width, int height, bool preserve=false)
			{boost::numeric::ublas::matrix<T>::resize(width,height,preserve);}
		void low_pass_filter();
		void high_pass_filter();
};

template <class T>
array2d<T>::array2d(){
	cellsize=-1;
	xllcorner=-1;
	yllcorner=-1;
	data_cells=-1;
	no_data=-1;
}

template<class T>
template<class U>
void array2d<T>::copyprops (const array2d<U> &copyfrom){
	cellsize=copyfrom.cellsize;
	xllcorner=copyfrom.xllcorner;
	yllcorner=copyfrom.yllcorner;
	data_cells=copyfrom.data_cells;
	no_data=copyfrom.no_data;
	resize(copyfrom.width(),copyfrom.height());
}

template <class T>
template <class U>
array2d<T>::array2d(const array2d<U> &copyfrom){
	array2d();
	copyprops(copyfrom);
}

template <class T>
void array2d<T>::init(T val){
	#pragma omp parallel for
	for(int x=0;x<width();x++)
	for(int y=0;y<height();y++)
		operator()(x,y)=val;
}

template <> inline long array2d<float>::estimated_output_size(){return 9*this->width()*this->height();}
template <> inline long array2d<char>::estimated_output_size(){return 4*this->width()*this->height();}
template <> inline long array2d<bool>::estimated_output_size(){return 2*this->width()*this->height();}
template <> inline long array2d<unsigned int>::estimated_output_size(){return 9*this->width()*this->height();}

//TODO: This is probably only useful for testing, since the file_io thing uses its own output
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
		if(operator()(x,y)!=other(x,y))
			return false;
	return true;
}

template <class T>
T array2d<T>::max() const {
	T maxval=no_data;
	#pragma omp parallel for collapse(2) reduction(max:maxval)
	for(int x=0;x<width();x++)
	for(int y=0;y<height();y++){
		T temp=operator()(x,y);
		if(temp==no_data)
			continue;
		else if(temp>maxval || maxval==no_data)
			maxval=temp;
	}
	return maxval;
}

template <class T>
T array2d<T>::min() const {
	T minval=no_data;
	#pragma omp parallel for collapse(2) reduction(min:minval)
	for(int x=0;x<width();x++)
	for(int y=0;y<height();y++){
		T temp=operator()(x,y);
		if(temp==no_data)
			continue;
		else if (temp<minval || minval==no_data)
			minval=temp;
	}
	return minval;
}



//A low pass filter smooths the data by reducing local variation and removing noise. The low pass filter calculates the average (mean) value for each 3 x 3 neighborhood. The effect is that the high and low values within each neighborhood will be averaged out, reducing the extreme values in the data.
//	http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Neighborhood%20filters
template <class T>
void array2d<T>::low_pass_filter(){	//TODO: Should provide diagnostics, I think
	array2d<T> filtered(*this);
	#pragma omp parallel for collapse(2)
	for(int x=0;x<width();x++)
	for(int y=0;y<height();y++){
		filtered(x,y)=operator()(x,y);
		if(operator()(x,y)==no_data)
			continue;

		int ncount=1;	//Middle cell is defined
		for(int n=0;n<8;n++){
			int nx=x+dinf_dx[n];
			int ny=y+dinf_dy[n];
			if(operator()(nx,ny)!=no_data){
				ncount++;
				filtered(x,y)+=operator()(nx,ny);
			}
		}
		filtered(x,y)/=ncount;
	}
	#pragma omp parallel for collapse(2)
	for(int x=0;x<width();x++)
	for(int y=0;y<height();y++)
		operator()(x,y)=filtered(x,y);
}


//The high pass filter accentuates the comparative difference in the values with its neighbors. A high pass filter calculates the focal sum statistic for each cell of the input using a 3 x 3 weighted kernel neighborhood. It brings out the boundaries between features (for example, where a water body meets the forest), thus sharpening edges between objects. The high pass filter is referred to as an edge enhancement filter.
//	http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Neighborhood%20filters
//Dinf Dirs			Weights
//321		-0.7	-1.0	-0.7
//4 0		-1.0	6.8		-1.0
//567		-0.7	-1.0	-0.7
//Weights sum to zero because they are normalized, according to ArcGIS
template <class T>
void array2d<T>::high_pass_filter(){	//TODO: Should provide diagnostics, I think
	const float weights[8]={-1.0,-0.7,-1.0,-0.7,1.0,-0.7,-1.0,-0.7};
	array2d<float> filtered(*this);
	#pragma omp parallel for collapse(2)
	for(int x=0;x<width();x++)
	for(int y=0;y<height();y++){
		if(operator()(x,y)==no_data){
			filtered(x,y)=no_data;
			continue;
		}

		filtered(x,y)=operator()(x,y)*6.8;
		for(int n=0;n<8;n++){
			int nx=x+dinf_dx[n];
			int ny=y+dinf_dy[n];
			if(operator()(nx,ny)!=no_data){
				filtered(x,y)+=operator()(nx,ny)*weights[n];
			}
		}
	}
	#pragma omp parallel for collapse(2)
	for(int x=0;x<width();x++)
	for(int y=0;y<height();y++)
		operator()(x,y)=(T)filtered(x,y);
}



typedef array2d<double> double_2d;
typedef array2d<float> float_2d;
typedef array2d<signed char> char_2d;
typedef array2d<bool> bool_2d;
typedef array2d<unsigned int> uint_2d;
typedef array2d<int> int_2d;



typedef struct grid_cell_type {
	int x;
	int y;
	grid_cell_type(int x0, int y0):x(x0),y(y0){}
	grid_cell_type(){}
} grid_cell;



typedef struct grid_cell_typez {
	int x;
	int y;
	float z;
	grid_cell_typez(int x0, int y0, float z0):x(x0),y(y0),z(z0){}
	grid_cell_typez(){}
} grid_cellz;

class grid_cellz_compare{
	bool reverse;
	public:
		grid_cellz_compare(const bool& revparam=false){reverse=revparam;}
		bool operator() (const grid_cellz &lhs, const grid_cellz &rhs) const{
			if (reverse) return (lhs.z<rhs.z);
			else return (lhs.z>rhs.z);
		}
};



typedef struct grid_cell_typezk {
	int x;
	int y;
	float z;
	int k;
	grid_cell_typezk(int x0, int y0, float z0, int k0):x(x0),y(y0),z(z0),k(k0){}
	grid_cell_typezk(){}
} grid_cellzk;

class grid_cellzk_compare{
	bool reverse;
	public:
		grid_cellzk_compare(const bool& revparam=false){reverse=revparam;}
		bool operator() (const grid_cellzk &lhs, const grid_cellzk &rhs) const{
			if (reverse) return (lhs.z<rhs.z || (lhs.z==rhs.z && lhs.k<rhs.k));
			else return (lhs.z>rhs.z || (lhs.z==rhs.z && lhs.k>rhs.k));
		}
};

#endif
