#ifndef _communication_hpp_
#define _communication_hpp_

#include <mpi.h>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/archives/binary.hpp>
#include <sstream>
#include <vector>
#include <iterator>
#include <sstream>
#include <cassert>
#include <iostream>

//For non-busy looping
#include <thread>
#include <chrono>

#define _unused(x) ((void)x) //Used for asserts

typedef unsigned long long comm_count_type;

static comm_count_type bytes_sent = 0;
static comm_count_type bytes_recv = 0;

//Consumers call CommGetTag() when they are waiting for work from the Producer.
//In many MPI implementations this leads to busy looping that can prevent
//processors from doing other useful work, such as apears here in the @retainall
//strategy. To avoid this, I have the processors sleep and periodically check
//for additional work. This variable sets the delay between checks.
static auto MPI_BUSY_LOOP_SLEEP = std::chrono::milliseconds(50);

void CommInit(int *argc, char ***argv){
  MPI_Init(argc,argv);
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

  bytes_sent += omsg.size();

  int ret = MPI_Send(omsg.data(), omsg.size(), MPI_BYTE, dest, tag, MPI_COMM_WORLD);
  assert(ret==MPI_SUCCESS);
  _unused(ret);
}

template<class T>
void CommSend(const T* a, std::nullptr_t, int dest, int tag){
  CommSend(a, (int*)nullptr, dest, tag);
}

int CommGetTag(int from){
  int flag;
  MPI_Status status;

  do {
    MPI_Iprobe(from, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    std::this_thread::sleep_for(MPI_BUSY_LOOP_SLEEP);
  } while (!flag);

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

  std::stringstream ss(std::stringstream::in|std::stringstream::out|std::stringstream::binary);
  ss.unsetf(std::ios_base::skipws);

  // Allocate a buffer to hold the incoming data
  char* buf = (char*)malloc(msg_size);
  assert(buf!=NULL);

  MPI_Recv(buf, msg_size, MPI_BYTE, from, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  bytes_recv += msg_size;

  ss.write(buf,msg_size);
  free(buf);

  {
    cereal::BinaryInputArchive archive(ss);
    archive(*a);
    if(b!=nullptr)
      archive(*b);
  }
}

template<class T>
void CommRecv(T* a, std::nullptr_t, int from){
  CommRecv(a, (int*)nullptr, from);
}

template<class T>
void CommBroadcast(T *datum, int root){
  MPI_Bcast(datum, 1, MPI_INT, root, MPI_COMM_WORLD);
}

void CommFinalize(){
  MPI_Finalize();
}

comm_count_type CommBytesSent(){
  return bytes_sent;
}

comm_count_type CommBytesRecv(){
  return bytes_recv;
}

void CommBytesReset(){
  bytes_recv = 0;
  bytes_sent = 0;
}

#endif