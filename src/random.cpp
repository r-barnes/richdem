#include <richdem/common/random.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>

namespace richdem {

our_random_engine& rand_engine(){
  static our_random_engine e[PRNG_THREAD_MAX];
  return e[omp_get_thread_num()];
}


//Be sure to read: http://www.pcg-random.org/posts/cpp-seeding-surprises.html
//and http://www.pcg-random.org/posts/cpps-random_device.html
void seed_rand(unsigned long seed){
  #pragma omp parallel      //All threads must come here
  {
    #pragma omp critical    //But only one at a time
    if(seed==0){
      std::uint_least32_t seed_data[std::mt19937::state_size];
      std::random_device r;
      std::generate_n(seed_data, std::mt19937::state_size, std::ref(r));
      std::seed_seq q(std::begin(seed_data), std::end(seed_data));
      rand_engine().seed(q);
    } else
      rand_engine().seed( seed*omp_get_thread_num() );
  }
}


int uniform_rand_int(int from, int thru){
  static std::uniform_int_distribution<> d[PRNG_THREAD_MAX];
  using parm_t = std::uniform_int_distribution<>::param_type;
  return d[omp_get_thread_num()]( rand_engine(), parm_t{from, thru} );
}


double uniform_rand_real(double from, double thru){
  static std::uniform_real_distribution<> d[PRNG_THREAD_MAX];
  using parm_t = std::uniform_real_distribution<>::param_type;
  return d[omp_get_thread_num()]( rand_engine(), parm_t{from, thru} );
}


double normal_rand(double mean, double stddev){
  static std::normal_distribution<double> d[PRNG_THREAD_MAX];
  using parm_t = std::normal_distribution<double>::param_type;
  return d[omp_get_thread_num()]( rand_engine(), parm_t{mean, stddev} );
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
