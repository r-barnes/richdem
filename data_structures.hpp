/**
  @file
  Template code defining data structures used throughout the package.
  Defines 2D arrays, grid cells, and priority queues for grid cells,
  among other.

  Richard Barnes (rbarnes@umn.edu), 2012
*/
#ifndef _data_structures_included
#define _data_structures_included

#define NDEBUG
#define BOOST_UBLAS_NDEBUG  //TODO: Not sure if this does anything
#include <boost/numeric/ublas/matrix.hpp>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <queue>

//Use this D8 method //TODO
//321
//4 0
//567
///Definition of x offsets of D-inf neighbours
const int dinf_dx[9]={1,1,0,-1,-1,-1,0,1,1};
///Definition of y offsets of D-inf neighbours
const int dinf_dy[9]={0,-1,-1,-1,0,1,1,1,0};
const int dinf_d8_inverse[9]={4,5,6,7,0,1,2,3,4};

template <class T>
class array2d : protected boost::numeric::ublas::matrix<T>{
  public:
    ///Value corresponding the length of one edge of a square DEM cell
    double cellsize;
//    std::string xllcorner,yllcorner;  //TODO: Should be string
    ///Global grid location of lower left x-coordinate
    double xllcorner;
    ///Global grid location of lower left y-coordinate
    double yllcorner;
    ///Number of cells containing data (excludes NO_DATA cells)
    long data_cells;
    ///NO_DATA value. The cell does not contain data, should not be processed.
    T no_data;

    array2d ();
    ///Creates array2d with no data, but invokes copyprops()
    template<class U> array2d (const array2d<U> &copyfrom);
    long width() const
      {return boost::numeric::ublas::matrix<T>::size1();}
    long height() const
      {return boost::numeric::ublas::matrix<T>::size2();}
    ///Copys everything but the data from another array2d
    template<class U> void copyprops (const array2d<U> &copyfrom);
    ///Prints an estimate of the file size were the array printed in ASCII
    long estimated_output_size();
    ///Sets all the cells of an array2d to "val"
    void init(T val);
    void print_block(
      std:: ostream& out, int minx, int maxx, int miny, int maxy,
      int precision=0, std::streamsize swidth=2
    ); //TODO
    void surroundings(int x0, int y0, int precision=3) const; //TODO
    ///Returns true if each cell of the array2d equals its counterpart in "other"
    bool operator==(const array2d<T> &other) const;
    template<class U> friend std::ostream& operator<<(
      std::ostream &out, const array2d<U> &arr
    );
    T max() const;
    T min() const;
    ///Returns true if (x,y) is within the bounds of the array2d
    inline bool in_grid(int x, int y) const
      {return (x>=0 && y>=0 && x<width() && y<height());}
    ///Returns true if (x,y) is not an edge cell
    inline bool interior_grid(int x, int y) const
      {return (x>=1 && y>=1 && x<width()-1 && y<height()-1);}
    ///Returns true if (x,y) lies on the edge of the array2d
    inline bool edge_grid(int x, int y) const
      {return (x==0 || y==0 || x==width()-1 || y==height()-1);}
    ///Returns a reference to (x,y)
    T& operator()(int x, int y)
      {return boost::numeric::ublas::matrix<T>::operator()(x,y);}
    ///Returns a const reference to (x,y)
    const T& operator()(int x, int y) const
      {return boost::numeric::ublas::matrix<T>::operator()(x,y);}
    ///Resizes the array2d. May or may not be destructive to existing data.
    void resize(int width, int height, bool preserve);
    ///Destroys all data in the array2d.
    void clear() {boost::numeric::ublas::matrix<T>::clear();}
    void low_pass_filter();
    void high_pass_filter();
    void print_random_sample(int n=1, int seed=1) const;

    template <class U> array2d<T>& operator*(U scalar);
    template <class U> array2d<T>& operator*(const array2d<U> B);
    template <class U> array2d<T>& operator+=(const array2d<U> B);
};

template <class T>
array2d<T>::array2d(){
  cellsize=-1;
  xllcorner=-1;
  yllcorner=-1;
  data_cells=-1;
  no_data=-1;
}

template <class T>
void array2d<T>::resize(int width, int height, bool preserve=false){
  fprintf(
    stderr,
    "\n\tApprox RAM requirement: %.2fMB\n",
    (float)width/1024. * (float)height/1024. * (float)sizeof(T)
  );
  boost::numeric::ublas::matrix<T>::resize(width,height,preserve);
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
  *this=array2d();
  copyprops(copyfrom);
}

template <class T>
void array2d<T>::init(T val){
  #pragma omp parallel for
  for(int x=0;x<width();x++)
  for(int y=0;y<height();y++)
    operator()(x,y)=val;
}

template <> inline long array2d<float>::estimated_output_size(){
  return 9*this->width()*this->height();
}
template <> inline long array2d<char>::estimated_output_size(){
  return 4*this->width()*this->height();
}
template <> inline long array2d<bool>::estimated_output_size(){
  return 2*this->width()*this->height();
}
template <> inline long array2d<unsigned int>::estimated_output_size(){
  return 9*this->width()*this->height();
}

//TODO: Probably only useful for testing, 
//  since the file_io thing uses its own output
template <class T>
std::ostream& operator<< (std::ostream &out, const array2d<T> &arr){
  std::streamsize width=out.width();
  for(int y=0;y<arr.height();y++){
    for(int x=0;x<arr.width();x++)
      if(sizeof(T)==1)  //TODO: An ugly way of detecting chars
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
template <class U>
array2d<T>& array2d<T>::operator*(U scalar) {
  #pragma omp parallel for collapse(2)
  for(int x=0;x<width();x++)
  for(int y=0;y<height();y++)
    if(operator()(x,y)!=no_data)
      operator()(x,y)*=(T)scalar;
  return *this;
}

template <class T>
template <class U>
array2d<T>& array2d<T>::operator*(array2d<U> B) {
  #pragma omp parallel for collapse(2)
  for(int x=0;x<width();x++)
  for(int y=0;y<height();y++)
    if(operator()(x,y)==no_data || B(x,y)==B.no_data)
      operator()(x,y)=no_data;
    else
      operator()(x,y)*=(T)B(x,y);
  return *this;
}

template <class T>
template <class U> array2d<T>& array2d<T>::operator+=(const array2d<U> B){
  #pragma omp parallel for collapse(2)
  for(int x=0;x<width();x++)
  for(int y=0;y<height();y++)
    if(operator()(x,y)==no_data || B(x,y)==B.no_data)
      operator()(x,y)=no_data;
    else
      operator()(x,y)+=(T)B(x,y);
  return *this;
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


//array2d.print_random_sample
/**
  @brief  Prints one or more random data_cells from the grid
  @author Richard Barnes (rbarnes@umn.edu)

  Prints one or more random data_cells from the grid. Note that if the grid is
  mostly no_data cells, the function may take a long time to complete as it
  will have to make many false attempts at finding data cells.

  @post The grid must have data_cells>0, or the function will throw an error

  @param[in]    n
    The number of data_cells to print. Default is 1.
  @param[in]    seed
    Seed to use. The default value is 1.
*/
template <class T>
void array2d<T>::print_random_sample(int n, int seed) const {
  if(data_cells==0)
    throw "Called print_random_sample() on a grid with no data_cells";

  srand(seed);

  for(int i=0;i<n;i++)
    while (true){
      int x=rand()%width();
      int y=rand()%height();
      if(operator()(x,y)!=no_data){
        std::cout<<"("<<x<<","<<y<<") = "<<operator()(x,y)<<std::endl;
        break;
      }
    }
}


/**
  @brief Smooths data by reducing local variation and removing noise with a neighbourhood average
  @author Richard Barnes (rbarnes@umn.edu)

  A low pass filter smooths the data by reducing local variation and removing
  noise. The low pass filter calculates the average (mean) value for each
  3 x 3 neighborhood. The effect is that the high and low values within each
  neighborhood will be averaged out, reducing the extreme values in the data.

  The filtered grid is initialized to the unfiltered grid and used to store
  the responses. The filtered grid is then copied into the unfilitered grid.

  See: http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Neighborhood%20filters

  @todo What if this overflows?
  @todo Progress bar? Other diagnostics?
*/
template <class T>
void array2d<T>::low_pass_filter(){
  array2d<T> filtered;
  filtered=*this;
  #pragma omp parallel for collapse(2)
  for(int x=0;x<width();x++)
  for(int y=0;y<height();y++){
    if(operator()(x,y)==no_data)
      continue;

    int ncount=1;  //Middle cell is defined
    for(int n=0;n<8;n++){
      int nx=x+dinf_dx[n];
      int ny=y+dinf_dy[n];
      if(in_grid(nx,ny) && operator()(nx,ny)!=no_data){
        ncount++;
        filtered(x,y)+=operator()(nx,ny);
      }
    }
    filtered(x,y)/=ncount;
  }

  *this=filtered;
}

/**
  @brief Accentuates differences between the cell and its neighbours.
  @author Richard Barnes (rbarnes@umn.edu)

  The high pass filter accentuates the comparative difference in the values
  with its neighbors. A high pass filter calculates the focal sum statistic
  for each cell of the input using a 3 x 3 weighted kernel neighborhood. It
  brings out the boundaries between features (for example, where a water body
  meets the forest), thus sharpening edges between objects. The high pass
  filter is referred to as an edge enhancement filter.

  See: http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Neighborhood%20filters

  \verbatim
  Dinf Dirs           Weights
  321            -0.7  -1.0  -0.7
  4 0            -1.0   6.8  -1.0
  567            -0.7  -1.0  -0.7
  \endverbatim

  Weights sum to zero because they are normalized, according to ArcGIS

  The filtered grid is initialized to the unfiltered grid and used to store
  the responses. The filtered grid is then copied into the unfilitered grid.

  @todo What if this overflows?
  @todo Progress bar? Other diagnostics?
*/
template <class T>
void array2d<T>::high_pass_filter(){
  const float weights[8]={-1.0,-0.7,-1.0,-0.7,1.0,-0.7,-1.0,-0.7};
  array2d<float> filtered;
  filtered=*this;
  #pragma omp parallel for collapse(2)
  for(int x=0;x<width();x++)
  for(int y=0;y<height();y++){
    if(operator()(x,y)==no_data)
      continue;

    filtered(x,y)*=6.8;
    for(int n=0;n<8;n++){
      int nx=x+dinf_dx[n];
      int ny=y+dinf_dy[n];
      if(in_grid(nx,ny) && operator()(nx,ny)!=no_data){
        filtered(x,y)+=operator()(nx,ny)*weights[n];
      }
    }
  }

  filtered=*this;
}



typedef array2d<double>       double_2d;
typedef array2d<float>        float_2d;
typedef array2d<signed char>  char_2d;
typedef array2d<bool>         bool_2d;
typedef array2d<unsigned int> uint_2d;
typedef array2d<int>          int_2d;


/// Stores the (x,y) coordinates of a grid cell
class grid_cell {
  public:
    int x; ///< Grid cell's x-coordinate
    int y; ///< Grid cell's y-coordinate
    /// Initiate the grid cell without coordinates; should generally be avoided
    grid_cell(){}
    /// Initiate the grid cell to the coordinates (x0,y0)
    grid_cell(int x, int y):x(x),y(y){}
};


/// Stores the (x,y,z) coordinates of a grid cell; useful for priority sorting
/// with \ref grid_cellz_compare
/// @todo z-coordinate should be templated
class grid_cellz : public grid_cell {
  public:
    float z;         ///< Grid cell's z-coordinate
    grid_cellz(int x, int y, float z): grid_cell(x,y), z(z) {}
    grid_cellz(){}
    bool operator< (const grid_cellz& a) const { return z< a.z; }
    bool operator> (const grid_cellz& a) const { return z> a.z; }
    bool operator>=(const grid_cellz& a) const { return z>=a.z; }
    bool operator<=(const grid_cellz& a) const { return z<=a.z; }
    bool operator==(const grid_cellz& a) const { return z==a.z; }
    bool operator!=(const grid_cellz& a) const { return !operator==(a); }
};



/// Stores the (x,y,z) coordinates of a grid cell and a priority indicator k;
/// used by grid_cellz_pq
/// @todo z-coordinate should be templated
class grid_cellzk : public grid_cellz {
  public:
    int k;           ///< Used to store an integer to make sorting stable
    grid_cellzk(int x, int y, float z, int k): grid_cellz(x,y,z), k(k) {}
    grid_cellzk(){}
    bool operator< (const grid_cellzk& a) const { return z< a.z || (z==a.z && k<a.k); }
    bool operator> (const grid_cellzk& a) const { return z> a.z || (z==a.z && k>a.k); }
};

///A priority queue of grid_cells, sorted by ascending height
class grid_cellz_pq : public std::priority_queue<grid_cellz, std::vector<grid_cellz>, std::greater<grid_cellz> > {
  public:
    void push_cell(int x, int y, float z){
      std::priority_queue<grid_cellz, std::vector<grid_cellz>, std::greater<grid_cellz> >::push(grid_cellz(x,y,z));
    }
};

///A priority queue of grid_cells, sorted by ascending height or, if heights
///if heights are equal, by the order of insertion
class grid_cellzk_pq : public std::priority_queue<grid_cellzk, std::vector<grid_cellzk>, std::greater<grid_cellzk> > {
  private:
    int count;
  public:
    grid_cellzk_pq() : count(0) {}
    void push(const grid_cellz &a){
      std::priority_queue<grid_cellzk, std::vector<grid_cellzk>, std::greater<grid_cellzk> >::push(grid_cellzk(a.x,a.y,a.z,count++));
    }
    void push_cell(int x, int y, float z){
      std::priority_queue<grid_cellzk, std::vector<grid_cellzk>, std::greater<grid_cellzk> >::push(grid_cellzk(x,y,z,count++));
    }
};

#endif
