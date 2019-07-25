#include <iostream>
#include <string>
#include <cstdlib>
#include <richdem/common/version.hpp>
#include <richdem/methods/flow_accumulation.hpp>
#include <richdem/common/Array2D.hpp>
using namespace richdem;

template<class T>
int PerformAlgorithm(std::string output, int algorithm, float param, std::string analysis, Array2D<T> dem){
  dem.loadData();

  Array2D<double> accum(dem,1);

  switch(algorithm){
    case 1: //D8     - O'Callaghan/Marks (1984)
      FA_D8(dem,accum);                                 break;
    case 2: //Rho8   - Fairfield & Leymarie (1991)
      FA_Rho8(dem,accum);                               break;
    case 3: //MD8    - Quinn (1991)
      FA_Quinn(dem,accum);                              break;
    case 4: //MD8    - Holmgren (1994)
      FA_Holmgren(dem,accum,param);                     break;
    case 5: //MD8    - Freeman (1991) 
      FA_Freeman(dem,accum,param);                      break;
    case 6: //D∞     - Tarboton (1997)
      FA_Tarboton(dem,accum);                           break;
    case 7: //MD∞    - Seibert & McGlynn (2007)
      std::cerr<<"This FA is temporarily disabled."; break;
      //FA_SeibertMcGlynn(dem,accum,param);               break;
    case 8: //D8-LTD - Orlandini et al. (2003)
      std::cerr<<"This FA is temporarily disabled."; break;
      //FA_Orlandini(dem,accum,OrlandiniMode::LTD,param); break;
    case 9: //D8-LAD - Orlandini et al. (2003)
      std::cerr<<"This FA is temporarily disabled."; break;
      //FA_Orlandini(dem,accum,OrlandiniMode::LAD,param); break;
  }

  accum.scale(accum.getCellArea());

  accum.saveGDAL(output,analysis);

  return 0;
}

#include "router.hpp"



int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);
  
  int   algorithm = 0;
  float param     = 0;

  if(argc<4 || argc>5){
    std::cerr<<"Calculate flow accumulation in terms of upstream area"<<std::endl;
    std::cerr<<argv[0]<<" <DEM file> <Output File> <Algorithm #> [Parameter]"<<std::endl;
    std::cerr<<"Algorithms:"<<std::endl;
    std::cerr<<" 1: D8     - O'Callaghan/Marks (1984)"<<std::endl;
    std::cerr<<" 2: Rho8   - Fairfield & Leymarie (1991)"<<std::endl;
    std::cerr<<" 3: MD8    - Quinn (1991)"<<std::endl;
    std::cerr<<" 4: MD8    - Holmgren (1994). Requires the parameter x. Suggested value: 4.0-6.0"<<std::endl;
    std::cerr<<" 5: MD8    - Freeman (1991). Requires the parameter p. Suggested value: 1.1"<<std::endl;
    std::cerr<<" 6: D∞     - Tarboton (1997)"<<std::endl;
    std::cerr<<" 7: MD∞    - Seibert & McGlynn (2007). Requires the parameter x. Suggested value: 1.0"<<std::endl;
    std::cerr<<" 8: D8-LTD - Orlandini et al. (2003). Requires the parameter x. Suggested value: 1.0"<<std::endl;
    std::cerr<<" 9: D8-LAD - Orlandini et al. (2003). Requires the parameter x. Suggested value: 1.0"<<std::endl;
    return -1;
  }

  algorithm = std::stoi(argv[3]);

  switch(algorithm){
    case 4:
    case 5:
    case 7:
    case 8:
    case 9:
      if(argc!=5){
        std::cerr<<"Parameter value is required!"<<std::endl;
        return -1;
      } else {
        param = std::stof(argv[4]);
      }
  }

  return PerformAlgorithm(std::string(argv[1]),std::string(argv[2]),algorithm,param,analysis);
}