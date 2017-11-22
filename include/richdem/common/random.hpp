//This file contains a number of functions for getting seeding random number
//generators and pulling numbers from them in a thread-safe manner.
#ifndef _richdem_random_hpp_
#define _richdem_random_hpp_

#include <random>
#include <string>

namespace richdem {

///Maximum number of threads this class should deal with
#define PRNG_THREAD_MAX 32

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num()  0
  #define omp_get_num_threads() 1
  #define omp_get_max_threads() 1
#endif

typedef std::string RandomEngineState;

typedef std::mt19937 our_random_engine;

//Returns a PRNG engine specific to the calling thread
our_random_engine& rand_engine();

//Seeds the PRNG engines using entropy from the computer's random device
void seed_rand(unsigned long seed);

//Returns an integer value on the closed interval [from,thru]
//Thread-safe
int uniform_rand_int(int from, int thru);

//Returns an floating-point value on the interval [from,thru)
//Thread-safe
double uniform_rand_real(double from, double thru);

//Returns a Gaussian-distributed value with specified mean and standard
//deviation. Thread-safe
double normal_rand(double mean, double stddev);

template<class T>
T uniform_bits(){
  std::uniform_int_distribution<T> 
    dist(std::numeric_limits<T>::lowest(),std::numeric_limits<T>::max());
  return dist( rand_engine() );
}

RandomEngineState SaveRandomState();
void SetRandomState(const RandomEngineState &res);

}

#endif
