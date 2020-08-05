#include <iostream>
#include <string>
#include <cstdlib>
#include <richdem/common/version.hpp>
#include <richdem/depressions/Lindsay2016.hpp>
#include <richdem/common/Array2D.hpp>

using namespace richdem;

template<class T>
int PerformAlgorithm(
  std::string outputname,
  LindsayMode mode,
  bool     eps_gradients,
  bool     fill_depressions,
  uint32_t max_path_len,
  double   max_depth,
  std::string analysis,
  Array2D<T> elevation
){
  elevation.loadData();

  Lindsay2016(elevation, mode, eps_gradients, fill_depressions, max_path_len, (T)max_depth);

  elevation.saveGDAL(outputname,analysis);

  return 0;
}

#include "router.hpp"

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=8){
    std::cerr<<"Eliminate all depressions via flooding."<<std::endl;
    std::cerr<<argv[0]<<" <Input> <Output name> <Mode> <Epsilon> <Fill Depressions> <Max Path Length> <Max Depth>"<<std::endl;
    std::cerr<<"\t<Mode>             - COMPLETE, SELECTIVE, CONSTRAINED"<<std::endl;
    std::cerr<<"\t<Epsilon>          - EPS/NOEPS: Whether to fill with epsilon-gradients"<<std::endl;
    std::cerr<<"\t<Fill Depressions> - FILL/NOFILL: Whether to fill depressions after breaching"<<std::endl;
    std::cerr<<"\t<Max Path Length>  - Integer [0, Inf): How many cells long a breach path can be"<<std::endl;
    std::cerr<<"\t<Max Depth>        - Float [0, Inf): How deep a breach path can be"<<std::endl;
    return -1;
  }

  LindsayMode mode;
  if(argv[3]==std::string("COMPLETE"))
    mode = LindsayMode::COMPLETE_BREACHING;
  else if(argv[3]==std::string("SELECTIVE"))
    mode = LindsayMode::SELECTIVE_BREACHING;
  else if(argv[3]==std::string("CONSTRAINED"))
    mode = LindsayMode::CONSTRAINED_BREACHING;
  else
    throw std::runtime_error("Unknown breaching mode!");

  bool eps_gradients;
  if(argv[4]==std::string("EPS"))
    eps_gradients = true;
  else if(argv[4]==std::string("NOEPS"))
    eps_gradients = false;
  else
    throw std::runtime_error("Unknown filling mode!");

  bool fill_depressions;
  if(argv[5]==std::string("FILL"))
    fill_depressions = true;
  else if(argv[5]==std::string("NOFILL"))
    fill_depressions = false;
  else
    throw std::runtime_error("Unknown filling mode!");

  uint32_t max_path_len = std::stoul(argv[6]);
  double   max_depth    = std::stod(argv[7]);

  return PerformAlgorithm(argv[1],argv[2],mode,eps_gradients,fill_depressions,max_path_len,max_depth,analysis);
}