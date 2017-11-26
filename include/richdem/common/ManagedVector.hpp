#ifndef _richdem_managed_vector_hpp_
#define _richdem_managed_vector_hpp_

#include <memory>

namespace richdem {

template<class T>
class ManagedVector {
 private:
  template<typename> friend class ManagedVector;

  std::unique_ptr<T[]> _data;
  bool   owned = true;
  size_t _size = 0;

 public:
  ManagedVector() = default;

  ManagedVector(size_t size, T default_val=T()){
    _size = size;
    _data.reset(new T[size]);
    for(size_t i=0;i<size;i++)
      _data[i] = default_val;
  }

  ManagedVector(T* data0, size_t size0){
    _data.reset(data0);
    _size = size0;
    owned = false;
  }

  //Copy from other
  template<class U>
  ManagedVector(const ManagedVector<U> &other) {
    _size = other.size();
     _data.reset(new T[other.size()]);
    for(size_t i=0;i<other.size();i++)
      _data[i] = other._data[i];    
  }

  //Copy Constructor
  ManagedVector(const ManagedVector<T> &other) {
    _size = other.size();
    _data.reset(new T[other.size()]);
    for(size_t i=0;i<other.size();i++)
      _data[i] = other._data[i];    
  }

  //Move Constructor
  template<class U>
  ManagedVector(ManagedVector<U> &&other) noexcept {
    _size       = other._size;
    _data       = std::move(other._data);
    owned       = other.owned;
    other.owned = true;
    other._size = 0;
  }

  ~ManagedVector(){
    if(!owned)
      _data.release();
  }

  //Copy assignment operator
  template<class U>
  ManagedVector<T>& operator=(const ManagedVector<U>& other){
    ManagedVector<T> tmp(other);  // re-use copy-constructor
    *this = std::move(tmp);       // re-use move-assignment
    return *this;
  }

  /** Move assignment operator */
  template<class U>
  ManagedVector<T>& operator=(ManagedVector<U>&& other) noexcept {
    _size       = other.size();
    _data       = std::move(other._data);
    owned       = other.owned;
    other.owned = true;
    other._size = 0;
    return *this;
  }

  T* data() {
    return _data.get();
  }

  bool empty() const {
    return _size==0;
  }

  size_t size() const {
    return _size;
  }

  //TODO: Keep old memory?
  void resize(size_t new_size) {
    if(_size!=new_size && !owned)
      throw std::runtime_error("Cannot resize unowned memory!");

    _data.reset(new T[new_size]);
    _size = new_size;
  }

  T& operator[](size_t i){
    return _data[i];
  }

  const T& operator[](size_t i) const {
    return _data[i];
  }  
};

}

#endif
