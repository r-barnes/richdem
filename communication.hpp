#ifndef _communication_hpp_
#define _communication_hpp_

#include <mpi.h>
#include <sstream>

#define _unused(x) ((void)x) //Used for asserts

template<class T, class U>
void CommSend(const T* a, const U* b, int dest, int tag){
  std::vector<uint8_t> omsg;
  {
    std::ostringstream oss;
    cereal::BinaryOutputArchive archive(oss);
    archive(*a);
    if(b!=nullptr)
      archive(*b);
    std::copy(omsg.begin(), omsg.end(), std::ostream_iterator<char>(oss));
  }

  int ret = MPI_Send(omsg.data(), omsg.size(), MPI_BYTE, dest, tag, MPI_COMM_WORLD);
  assert(ret==MPI_SUCCESS);
  _unused(ret);
}

template<class T>
void CommSend(const T* a, nullptr_t, int dest, int tag){
  CommSend(a, (int*)nullptr, dest, tag);
}

int CommGetTag(int from){
  MPI_Status status;
  MPI_Probe(from, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  return status.MPI_TAG;
}

int CommGetSource(){
  MPI_Status status;
  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  return status.MPI_SOURCE;
}

int CommRank(){
  int rank;
  int ret = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  assert(ret==MPI_SUCCESS);
  _unused(ret);
  return rank;
}

int CommSize(){
  int size;
  int ret = MPI_Comm_size (MPI_COMM_WORLD, &size);
  assert(ret==MPI_SUCCESS);
  _unused(ret);
  return size;
}

void CommAbort(int errorcode){
  MPI_Abort(MPI_COMM_WORLD, errorcode);
}

template<class T, class U>
void CommRecv(T* a, U* b, int from){
  MPI_Status status;
  MPI_Probe(from, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  int msg_size;
  MPI_Get_count(&status, MPI_BYTE, &msg_size);

  std::istringstream iss;

  // Allocate a buffer to hold the incoming data
  char* buf = (char*)malloc(msg_size);
  assert(buf!=NULL);

  MPI_Recv(buf, msg_size, MPI_BYTE, from, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
  iss.read(buf,msg_size);
  free(buf);

  {
    cereal::BinaryInputArchive archive(iss);
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
  MPI_Bcast(&datum, 1, MPI_INT, root, MPI_COMM_WORLD);
}

#endif