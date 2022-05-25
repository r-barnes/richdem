#ifndef _richdem_upslope_cells_hpp_
#define _richdem_upslope_cells_hpp_

#include <richdem/common/Array2D.hpp>
#include <richdem/common/Array3D.hpp>
#include <richdem/methods/upslope_cells_functions.hpp>

namespace richdem {

template<class T, class U>
void init(
    Array2D<T> &to_init,
    const Array2D<U> &example,
    const int &val
){
    to_init.resize(example);
    to_init.setAll(val);
    to_init.setNoData(val);
}


// TODO: do we have to resize  and initialize the results array?
template<class mask_t,               class result_t> void UC_mask_props(const Array2D<mask_t> &mask,    const Array3D<float> &props, Array2D<result_t> &upslope){init(upslope,mask,0);  std::queue<GridCell> expansion = queue_from_mask(mask,upslope)             ; upslope_cells_props<Topology::D8>(expansion, props, upslope);}
template<class mask_t, class elev_t, class result_t> void UC_mask_mflow(const Array2D<mask_t> &mask,    const Array2D<elev_t> &elev, Array2D<result_t> &upslope){init(upslope,elev,0);  std::queue<GridCell> expansion = queue_from_mask(mask,upslope)             ; upslope_cells_mflow<Topology::D8>(expansion, elev,  upslope);}
template<                            class result_t> void UC_line_props(int x0, int y0, int x1, int y1, const Array3D<float> &props, Array2D<result_t> &upslope){init(upslope,props,0); std::queue<GridCell> expansion = queue_from_linepoints(x0,y0,x1,y1,upslope); upslope_cells_props<Topology::D8>(expansion, props, upslope);}
template<              class elev_t, class result_t> void UC_line_mflow(int x0, int y0, int x1, int y1, const Array2D<elev_t> &elev, Array2D<result_t> &upslope){init(upslope,elev,0);  std::queue<GridCell> expansion = queue_from_linepoints(x0,y0,x1,y1,upslope); upslope_cells_mflow<Topology::D8>(expansion, elev,  upslope);}

}
#endif //_richdem_upslope_cells_