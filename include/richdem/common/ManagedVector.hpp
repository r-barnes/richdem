#ifndef _richdem_managed_vector_hpp_
#define _richdem_managed_vector_hpp_

#include <memory>
//#include <iostream>

namespace richdem {

///ManagedVector works like a regular vector, but can wrap external memory
template<class T>
class ManagedVector {
 private:
  template<typename> friend class ManagedVector;

  std::unique_ptr<T[]> _data;
  bool   _owned = true;        ///< If this is true, we are responsible for clean-up of the data
  std::size_t _size  = 0;           ///< Number of elements being managed

 public:
  ///Creates an empty ManagedVector
  ManagedVector() = default;

  ///Creates a ManagedVector with \p size members each set to \p default_val
  ///
  ///@param[in] size         Number of elements to be created in the vector
  ///@param[in] default_val  Initial value of the elements
  ManagedVector(std::size_t size, T default_val=T()){
    _size = size;
    _data.reset(new T[size]);   //Create unique pointer to newly allocated memory
    for(std::size_t i=0;i<size;i++)  //Initialize all members to default value
      _data[i] = default_val;
    //std::cerr<<"ManagedVector construct with owned data at "<<((void*)_data.get())<<std::endl;
  }

  ///Creates a ManagedVector which wraps \p data0 of length \p size0
  ///
  ///@param[in] data         Memory to wrap
  ///@param[in] size         Number of elements to wrap
  ManagedVector(T* data, std::size_t size){
    //std::cerr<<"ManagedVector construct with unowned data at "<<((void*)data0)<<std::endl;
    _data.reset(data);
    _size  = size;
    _owned = false;
  }

  //Copy from other
  template<class U>
  ManagedVector(const ManagedVector<U> &other) {
    //std::cerr<<"Copying ManagedVector from U"<<std::endl;
    _size = other.size();
     _data.reset(new T[other.size()]);
    for(std::size_t i=0;i<other.size();i++)
      _data[i] = other._data[i];    
  }

  //TODO: Should this clear old data?
  //Copy Constructor
  ManagedVector(const ManagedVector<T> &other) {
    //std::cerr<<"Copying ManagedVector from T"<<std::endl;
    _size = other.size();
    _data.reset(new T[other.size()]);
    for(std::size_t i=0;i<other.size();i++)
      _data[i] = other._data[i];    
  }

  //TODO: Should this clear old data?
  //Move Constructor
  template<class U>
  ManagedVector(ManagedVector<U> &&other) noexcept {
    //std::cerr<<"Move constructor ManagedVector"<<std::endl;
    _size       = other._size;
    _data       = std::move(other._data);
    _owned       = other._owned;
    other._owned = true;
    other._size = 0;
  }

  ~ManagedVector(){
    if(!_owned){
      //std::cerr<<"ManagedVector releasing data at "<<((void*)_data.get())<<std::endl;
      _data.release();
    }
    // if(_data.get()!=nullptr){
    //   std::cerr<<"ManagedVector freeing its data at "<<((void*)_data.get())<<std::endl;
    // } else {  
    //   std::cerr<<"ManagedVector destructing with no data..."<<std::endl;
    // }
  }

  //Copy assignment operator
  template<class U>
  ManagedVector<T>& operator=(const ManagedVector<U>& other){
    //std::cerr<<"ManagedVector copy assignment of "<<((void*)other._data.get())<<std::endl;
    ManagedVector<T> tmp(other);  // re-use copy-constructor
    *this = std::move(tmp);       // re-use move-assignment
    return *this;
  }

  /** Move assignment operator */
  template<class U>
  ManagedVector<T>& operator=(ManagedVector<U>&& other) noexcept {
    //std::cerr<<"ManagedVector move assignment of "<<((void*)other._data.get())<<std::endl;
    _size        = other.size();
    _data        = std::move(other._data);
    _owned       = other._owned;
    other._owned = true;
    other._size  = 0;
    return *this;
  }

  ///Get a raw pointer to the managed data
  ///
  ///@return A raw pointer to the managed data
  inline T* data() {
    return _data.get();
  }

  ///Get a raw constant pointer to the managed data
  ///
  ///@return A raw constant pointer to the managed data
  inline const T* data() const {
    return _data.get();
  }

  ///Are there more than zero elements being managed?
  ///
  ///@return True, if zero elements are managed; otherwise, false
  inline bool empty() const {
    return _size==0;
  }

  ///Get the number of elements being managed
  ///
  ///@return The number of elements being managed
  inline std::size_t size() const {
    return _size;
  }

  ///Determine whether the ManagedVector owns the memory it is managing
  ///
  ///@return True, if this ManagedVector owns its memory; otherwise, false
  inline bool owned() const {
    return _owned;
  }

  //TODO: Keep old memory?
  void resize(std::size_t new_size) {
    if(_size==new_size)
      return;
    if(!_owned)
      throw std::runtime_error("Cannot resize unowned memory!");

    //auto oldmem = _data.get();

    _data.reset(); //Clear old memory before allocating new memory
    _data.reset(new T[new_size]);

    //std::cerr<<"Resizing ManagedVector from "<<((void*)oldmem)<<" to "<<((void*)_data.get())<<std::endl;

    _size = new_size;
  }

  inline T& operator[](std::size_t i){
    return _data[i];
  }

  inline const T& operator[](std::size_t i) const {
    return _data[i];
  }  
};

}

#endif
