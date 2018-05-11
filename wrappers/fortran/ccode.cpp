#include <richdem/common/Array2D.hpp>
#include <richdem/depressions/depressions.hpp>
#include <richdem/methods/flow_accumulation.hpp>
using namespace std;

namespace rd = richdem;

extern "C" {

void FloodDepressions(float *const dem, const int width, const int height, const float no_data){
   rd::Array2D<float> rd_dem(dem, width, height);
   rd_dem.setNoData(no_data);

   rd::FillDepressions(rd_dem);
}

}
