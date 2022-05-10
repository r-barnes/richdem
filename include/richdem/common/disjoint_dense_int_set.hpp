//DisjointDenseIntSet
//Author: Richard Barnes (rijard.barnes@gmail.com)
#ifndef _disjoint_dense_int_sets_
#define _disjoint_dense_int_sets_

#include <vector>
#include <cassert>


//A disjoint-set/union-find class. Starting from a collection of sets, this data
//structure efficiently keeps track of which sets have been merged. It is
//assumed that every integer between 0 and some maximum value N is a set, so the
//data structure takes O(N) space. If only `findSet()` and `unionSet()` are
//used, then all accesses are in O(a(N)) time, where `a()` is the inverse
//Ackermann function. For all practical purposes, this is `(1). If
//`mergeAintoB()` is used then `findSet()` can have a worst-case of `O(N)`.
class DisjointDenseIntSet {
 private:
  //How many sets are children of this set. Initially 0.
  std::vector<unsigned int> rank;
  //Which set is this set's parent. May be the set itself.
  std::vector<unsigned int> parent;

  //When a set value X is passed to the data structure, this checks to see if
  //the set exists. If not, the set is created, along with all the sets between
  //the old maximum set id and X.
  void checkSize(unsigned int newN){
    if(newN<rank.size())               //Does the set exist?
      return;                          //Yup.
                                       //Nope.
    const auto old_size = rank.size(); //Get old maximum set value

    //Resize so that `newN` is a valid value. None of the new sets have
    //children.
    rank.resize(newN+1,0);

    //Resize so that `newN` is a valid value
    parent.resize(newN+1);

    //Ensure that each new set is its own parent since they have not yet been
    //merged.
    for(auto i=old_size;i<newN+1;i++)
      parent[i] = i;
  }

 public:
  //Construct a DisjointDenseIntSet without any sets. Sets will be dynamically
  //created as the data structure is used.
  DisjointDenseIntSet(){}

  //Create a DisjointDenseIntSet with `N` initial sets preallocated. More sets
  //can be dynamically allocated as the data structure is used.
  DisjointDenseIntSet(unsigned int N){
    rank.resize(N,0);
    parent.resize(N);
    for(unsigned int i=0;i<N;i++)
      parent[i] = i;
  }

  //Explicitly creates a set. May incidentally create several intermediate sets
  //of `n` is more than one larger than the maximum set id previously seen.
  void makeSet(unsigned int n){
    checkSize(n);
  }

  //Returns the highest set id.
  unsigned int maxElement() const {
    return rank.size()-1;
  }

  //Follows a set's chain of parents until a set which is its own parent is
  //reached. This ultimate parent's id is returned as the representative id of
  //the set in question. Note that this collapses the chain of parents so that
  //after this method has run every set between the one in question and the
  //ultimate parent points to the ultimate parent. This means that while the
  //first call to this function may take `O(N)` lookups in the worst-case (less
  //due to the use of ranks, as explained below), subsequent calls to any set in
  //the chain will take `O(1)` time. This technique is known as "path
  //compression".
  unsigned int findSet(unsigned int n){
    if(n>=parent.size()){
      throw std::runtime_error("DisjointDenseIntSet::findSet(" + std::to_string(n) + ") is looking for a set outside the valid range, which is [0," + std::to_string(parent.size()) + ")!");
    }
    if(parent[n]==n){                  //Am I my own parent?
      return n;                        //Yes: I represent the set in question.
    } else {                           //No.
      //Who is my parent's ultimate parent? Make them my parent.
      return parent[n] = findSet(parent[n]);
    }
  }

  //Join two sets into a single set. Note that we "cannot" predict the `id` of
  //the resulting set ahead of time.
  void unionSet(unsigned int a, unsigned int b){
    const auto roota = findSet(a); //Find the ultimate parent of A
    const auto rootb = findSet(b); //Find the ultimate parent of B
    //Note that the foregoing collapses any chain of parents so that each set in
    //the chain points to the ultimate parent. Therefore, any subsequent call to
    //`findSet` involving any set in the chain will take `O(1)` time.

    //If A and B already share a parent, then they do not need merging.
    if(roota==rootb)
      return;

    //If we always naively tacked A onto B then we could develop a worst-case
    //scenario in which each set pointed to just one other set in a long, linear
    //chain. If this happened then calls to `findSet()` would take `O(N)` time.
    //Instead, we keep track of how many child sets each set has and ensure that
    //the shorter tree of sets becomes part of the taller tree of sets. This
    //ensures that the tree does not grow taller unless the two trees were of
    //equal height in which case the resultant tree is taller by 1. In essence,
    //this bounds the depth of any query to being `log_2(N)`. However, due to
    //the use of path compression above, the query path is actually less than
    //this.

    if(rank[roota]<rank[rootb]){          //Is A shorter?
      parent[roota] = rootb;              //Make A a child of B.
    } else if(rank[roota]>rank[rootb]) {  //Is B shorter?
      parent[rootb] = roota;              //Make B a child of A.
    } else {                              //A and B are the same height
      parent[rootb] = roota;              //Arbitrarily make B a child of A
      rank[roota]++;                      //Increase A's height.
    }
  }

  //Using `unionSet` merges two sets in a way which does not allow us to decide
  //which set is the parent; however, `unionSet` helps guarantee fast queries.
  //`mergeAintoB` sacrifices speed but preserves parenthood by always making A a
  //child of B, regardless of the height of `B`.
  void mergeAintoB(unsigned int a, unsigned int b){
    checkSize(a);
    checkSize(b);
    parent[a] = b;
    if(rank[a]==rank[b]){
      rank[b]++;
    } else if(rank[a]>rank[b]){
      rank[b] = rank[a]+1;
    } else {
      //If `rank[b]>rank[a]` then making A a child of B does not increase B's
      //height.
    }
  }

  //Returns true if A and B belong to the same set.
  bool sameSet(unsigned int a, unsigned int b){
    return findSet(a)==findSet(b);
  }
};

#endif
