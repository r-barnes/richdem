/**
  @file
  @brief Defines a Least-Recently Used (LRU) cache class.

  Richard Barnes (rbarnes@umn.edu), 2016
*/
#ifndef _lru_hpp_
#define _lru_hpp_

#include <list>
#include <unordered_map>

namespace richdem {

///@brief A Least-Recently Used (LRU) cache
template<class T>
class LRU {
 private:
  typedef std::list<T> cachetype; ///< Container used for storage by the cache
  int len;                        ///< Number of items in the cache
  int maxlen;                     ///< Maximum size the cache is allowed to be
  T   last;                       ///< Copy of the last inserted item. Speeds up insertion. TODO: This'd be better as a pointer
  std::unordered_map<T, typename std::list<T>::iterator> visited; ///< Used for O(1) access to members

 public:
  cachetype cache;                ///< The cache

  ///@brief Construct a new LRU
  LRU(){
    len    = 0;
    maxlen = -1;
  }

  ///@brief Insert an item into the LRU.
  ///
  ///The item is either added to the queue or its entry is moved to the top of
  ///the queue. If the item is new and the length of the queue is greater than
  ///maxlen, then the least recently seen item is evicted from the queue.
  ///
  ///@param entry The item to add to the queue.
  void insert(const T &entry){
    //Used to quickly move on if we've seen this item before (TODO: Should use a
    //pointer in case T is not a POD.)
    if(entry==last)
      return;

    if(visited.count(entry)){                     //Item is already in the queue
      auto existing_entry = visited[entry];       //Get item location
      cache.splice(cache.begin(),cache,existing_entry); //Move item to the front
    } else {                                      //Item is not in the queue
      cache.push_front(entry);                    //Add item to the front
      visited[entry] = cache.begin();             //Make a note of its location
      len++;                                      //Queue just got bigger
    }

    last = entry;                                 //Make a note of the new item
  }

  ///@brief Returns the number of itmes in the LRU cache
  ///@return Number of items in the LRU cache
  int size() const {
    return len;
  }

  ///@brief Is the LRU cache full?
  ///@return True if the LRU cache is full; otherwise, false.
  bool full() const {
    return len==maxlen;
  }

  ///@brief Set the maximum capacity of the LRU cache
  void setCapacity(int n){
    maxlen = n;
  }

  ///@brief Returns the capacity of the LRU cache
  ///@return The capacity of the LRU cache
  int getCapacity() const {
    return maxlen;
  }

  ///@brief Return the least-recently used item in the LRU cache.
  ///@return The least-recently used item in the LRU cache.
  T back() const {
    return cache.back();
  }

  ///@brief Evict the least-recently used item out of the LRU cache.
  void pop_back() {    
    visited.erase(cache.back());
    cache.pop_back();
    len--;
  }

  ///@brief Evict items from the LRU cache until it is within its capacity.
  void prune(){
    if(maxlen == -1) return;

    while(size()>maxlen)
      pop_back();
  }
};

}

#endif
