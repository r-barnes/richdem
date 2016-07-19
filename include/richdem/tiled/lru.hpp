#ifndef _lru_hpp_
#define _lru_hpp_

#include <list>
#include <unordered_map>

template<class T>
class LRU {
 private:
  typedef std::list<T> cachetype;
  int len;
  int maxlen;
  T   last;
  std::unordered_map<T, typename std::list<T>::iterator> visited;

 public:
  cachetype cache;

  LRU(){
    len    = 0;
    maxlen = -1;
  }

  void insert(const T &entry){
    if(entry==last)
      return;

    if(visited.count(entry)){
      auto existing_entry = visited[entry];
      cache.splice(cache.begin(),cache,existing_entry);
    } else {
      cache.push_front(entry);
      visited[entry]=cache.begin();
      len++;
    }

    last = entry;
  }

  int size() const {
    return len;
  }

  bool full() const {
    return len==maxlen;
  }

  void setCapacity(int n){
    maxlen = n;
  }

  int getCapacity() const {
    return maxlen;
  }

  T back() const {
    return cache.back();
  }

  void pop_back() {    
    visited.erase(cache.back());
    cache.pop_back();
    len--;
  }

  void prune(){
    if(maxlen == -1) return;

    while(size()>maxlen)
      pop_back();
  }
};

#endif