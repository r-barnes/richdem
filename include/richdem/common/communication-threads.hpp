#ifndef _communication_hpp_
#define _communication_hpp_

#error This test module is not yet ready for compilation.

#ifdef COMM_MPI
  #include <mpi.h>
#elif COMM_THREAD
  #include <thread>
#else
  #error Must delcare COMM_MPI or COMM_THREAD
#endif

#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/archives/binary.hpp>
#include <sstream>
#include <vector>
#include <queue>
#include <map>
#include <atomic>
#include <iterator>
#include <sstream>
#include <cassert>
#include <iostream>
#include <atomic> //Used for barrier
#include <list>   //Used for message passing
#include <chrono> //Used for message passing

#define _unused(x) ((void)x) //Used for asserts

#ifdef COMM_MPI
  static long bytes_sent = 0;
  static long bytes_recv = 0;
#elif COMM_THREAD
  std::vector<long> bytes_sent;
  std::vector<long> bytes_recv;
#endif

#ifdef COMM_THREAD
  std::vector<std::thread> threadpool;
  std::map<int,std::thread::id> id_to_tid;
  std::map<std::thread::id,int> tid_to_id;
#endif

#ifdef COMM_THREAD
  struct msg {
    int tag;
    int from;
    std::vector<char> msg;
  };

  std::vector< std::list< struct msg > > > msgqueue;

  static std::mutex comm_msg_mutex;
  static void CommPostMessage(int post_to, int tag, const std::vector<char> &msg){
    std::lock_guard<std::mutex> lock(comm_msg_mutex);
    msgqueue[post_to].emplace_back(tag,CommRank(),msg);
  }

  static void CommRecvMessage(int tag, int from, const std::vector<char> msg){
    int my_id = CommRank();
    while(true){
      {
        std::lock_guard<std::mutex> lock(comm_msg_mutex);
        if(!msgqueue[my_id].empty()){
          for(auto i=msgqueue[my_id].begin();i!=msgqueue[my_id].end();++i){
            if(i->tag!=tag || i->from!=from)
              continue;
            msg = i->msg;
            msgqueue[my_id].erase(i);
            return;
          }
        }
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
  }
#endif

template<class Fn>
void CommInit(int n, Fn&& fn, int *argc, char ***argv){
  #ifdef COMM_MPI
    MPI_Init(argc,argv);
    fn(argc,argv);
  #elif COMM_THREAD
    if(n==-1)
      n = std::thread::hardware_concurrency();

    //Do this here to be sure that CommSize returns the correct value
    threadpool.resize(n);
    bytes_sent.resize(n,0);
    bytes_recv.resize(n,0);

    for(int i=0;i<n;i++){
      threadpool[n]=std::thread(fn,argc,argv);
      id_to_tid[i]=threadpool.back().get_id();
      tid_to_id[threadpool.back().get_id()]=i;
    }
  #endif
}

template<class T, class U>
void CommSend(const T* a, const U* b, int dest, int tag){
  std::vector<char> omsg;
  {
    std::stringstream ss(std::stringstream::in|std::stringstream::out|std::stringstream::binary);
    ss.unsetf(std::ios_base::skipws);
    cereal::BinaryOutputArchive archive(ss);
    archive(*a);
    if(b!=nullptr)
      archive(*b);

    std::copy(std::istream_iterator<char>(ss), std::istream_iterator<char>(), std::back_inserter(omsg));
  }

  #ifdef COMM_MPI
    bytes_sent += omsg.size();
  #elif COMM_THREAD
    bytes_sent[CommRank()] += omsg.size();
  #endif

  #ifdef COMM_MPI
    int ret = MPI_Send(omsg.data(), omsg.size(), MPI_BYTE, dest, tag, MPI_COMM_WORLD);
    assert(ret==MPI_SUCCESS);
    _unused(ret);
  #elif COMM_THREAD
    CommPostMessage(dest,tag,omsg);
  #endif
}

template<class T>
void CommSend(const T* a, nullptr_t, int dest, int tag){
  CommSend(a, (int*)nullptr, dest, tag);
}

int CommGetTag(int from){
  #ifdef COMM_MPI
    MPI_Status status;
    MPI_Probe(from, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  #endif

  return status.MPI_TAG;
}

int CommGetSource(){
  #ifdef COMM_MPI
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    return status.MPI_SOURCE;
  #endif
}

int CommRank(){
  #ifdef COMM_MPI
    int rank;
    int ret = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    assert(ret==MPI_SUCCESS);
    _unused(ret);
    return rank;
  #elif COMM_THREAD
    return tid_to_id.at(std::this_thread::get_id());
  #endif
}

int CommSize(){
  #ifdef COMM_MPI
    int size;
    int ret = MPI_Comm_size (MPI_COMM_WORLD, &size);
    assert(ret==MPI_SUCCESS);
    _unused(ret);
    return size;
  #elif COMM_THREAD
    return threadpool.size();
  #endif
}

void CommAbort(int errorcode){
  #ifdef COMM_MPI
    MPI_Abort(MPI_COMM_WORLD, errorcode);
  #endif
}

template<class T, class U>
void CommRecv(T* a, U* b, int from){
  #ifdef COMM_MPI
    MPI_Status status;
    MPI_Probe(from, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    int msg_size;
    MPI_Get_count(&status, MPI_BYTE, &msg_size);

    // Allocate a buffer to hold the incoming data
    char* buf = (char*)malloc(msg_size);
    assert(buf!=NULL);

    MPI_Recv(buf, msg_size, MPI_BYTE, from, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    #ifdef COMM_MPI
      bytes_recv += msg_size;
    #elif COMM_THREAD
      bytes_recv[CommRank()] += msg_size;
    #endif

    std::stringstream ss(std::stringstream::in|std::stringstream::out|std::stringstream::binary);
    ss.unsetf(std::ios_base::skipws);
    ss.write(buf,msg_size);
    free(buf);
  #endif

  {
    cereal::BinaryInputArchive archive(ss);
    archive(*a);
    if(b!=nullptr)
      archive(*b);
  }
}

template<class T>
void CommRecv(T* a, nullptr_t, int from){
  CommRecv(a, (int*)nullptr, from);
}

template<class T>
void CommBroadcast(T *datum, int root){
  #ifdef COMM_MPI
    MPI_Bcast(datum, 1, MPI_INT, root, MPI_COMM_WORLD);
  #endif
}

void CommFinalize(){
  #ifdef COMM_MPI
    MPI_Finalize();
  #endif
}

int CommBytesSent(){
  #ifdef COMM_MPI
    return bytes_sent;
  #elif COMM_THREAD
    return bytes_sent[CommRank()];
  #endif
}

int CommBytesRecv(){
  #ifdef COMM_MPI
    return bytes_recv;
  #elif COMM_THREAD
    return bytes_recv[CommRank()];
  #endif
}

void CommBytesReset(){
  #ifdef COMM_MPI
    bytes_recv = 0;
    bytes_sent = 0;
  #elif COMM_THREAD
    bytes_sent[CommRank()] = 0;
    bytes_recv[CommRank()] = 0;
  #endif
}


void CommBarrier(){
  #ifdef COMM_THREAD
    static std::atomic<int> tcount = 0;
    tcount++;
    while(tcount<CommSize()-1){} //Busy loop until all threads present
  #endif
}

#endif
