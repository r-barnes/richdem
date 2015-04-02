#ifndef _disjoint_dense_int_sets_
#define _disjoint_dense_int_sets_

#include <vector>
#include <cassert>

class DisjointDenseIntSet {
 private:
  std::vector<unsigned int> rank;
  std::vector<int> parent;
  void checkSize(unsigned int n){
    if(n+1<=rank.size())
      return;

    rank.resize(n+1,0);
    parent.resize(n+1,-1);
  }
 public:
  DisjointDenseIntSet(){}
  DisjointDenseIntSet(int max_element){
    assert(max_element>=0);
    rank.resize(max_element+1,0);
    parent.resize(max_element+1);
    for(int i=0;i<=max_element;i++)
      parent[i] = i;
  }
  void makeSet(int n){
    assert(n>=0);
    checkSize(n);
    parent[n] = n;
    rank[n]   = 0;
  }
  unsigned int maxElement() const {
    return rank.size()-1;
  }
  int findSet(int n){
    assert(n>=0);
    if(parent[n]==-1 || parent[n]==n)
      return n;
    else
      return parent[n] = findSet(parent[n]);
  }
  void unionSet(int a, int b){
    assert(a>=0);
    assert(b>=0);
    auto roota = findSet(a);
    auto rootb = findSet(b);
    if(roota==rootb || roota==-1 || rootb==-1)
      return;

    if(rank[roota]<rank[rootb]){
      parent[roota] = rootb;
    } else if(rank[roota]>rank[rootb]) {
      parent[rootb] = roota;
    } else {
      parent[rootb] = roota;
      rank[roota]++;
    }
  }
  bool sameSet(int a, int b){
    assert(a>=0);
    assert(b>=0);
    return findSet(a)==findSet(b);
  }
};

#endif