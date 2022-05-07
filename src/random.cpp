#include <richdem/common/random.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>

namespace richdem {

static int rd_get_thread_num(){
  #ifdef _OPENMP
    return omp_get_thread_num();
  #else
    return 0;
  #endif
}

our_random_engine& rand_engine(){
  static std::array<our_random_engine, PRNG_THREAD_MAX> e;
  return e[rd_get_thread_num()];
}


//Be sure to read: http://www.pcg-random.org/posts/cpp-seeding-surprises.html
//and http://www.pcg-random.org/posts/cpps-random_device.html
void seed_rand(uint64_t seed){
  #pragma omp parallel default(none) shared(seed) //All threads must come here
  {
    #pragma omp critical    //But only one at a time
    if(seed==0){
      std::array<std::uint_least32_t, std::mt19937::state_size> seed_data;
      std::random_device r;
      std::generate_n(seed_data.data(), std::mt19937::state_size, std::ref(r));
      std::seed_seq q(std::begin(seed_data), std::end(seed_data));
      rand_engine().seed(q);
    } else {
      rand_engine().seed( seed*rd_get_thread_num() );
    }
  }
}


int uniform_rand_int(int from, int thru){
  static std::array<std::uniform_int_distribution<>, PRNG_THREAD_MAX> d;
  using parm_t = std::uniform_int_distribution<>::param_type;
  return d[rd_get_thread_num()]( rand_engine(), parm_t{from, thru} );
}


double uniform_rand_real(double from, double thru){
  static std::array<std::uniform_real_distribution<>, PRNG_THREAD_MAX> d;
  using parm_t = std::uniform_real_distribution<>::param_type;
  return d[rd_get_thread_num()]( rand_engine(), parm_t{from, thru} );
}


double normal_rand(double mean, double stddev){
  static std::array<std::normal_distribution<double>, PRNG_THREAD_MAX> d;
  using parm_t = std::normal_distribution<double>::param_type;
  return d[rd_get_thread_num()]( rand_engine(), parm_t{mean, stddev} );
}

RandomEngineState SaveRandomState(){
  std::ostringstream oss;
  oss<<rand_engine();
  return oss.str();
}

void SetRandomState(const RandomEngineState &res){
  std::istringstream iss(res);
  iss>>rand_engine();
}

}
