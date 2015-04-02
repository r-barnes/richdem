#ifndef _fill_extent_
#define _fill_extent_

#include <boost/serialization/map.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

template<class elev_t>
class FillExtent {
 private:
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version){
    ar & xmin;
    ar & xmax;
    ar & ymin;
    ar & ymax;
    ar & elevation;
    ar & init;
  }
 public:
  int xmin, xmax, ymin, ymax;
  elev_t elevation;
  bool init;
  FillExtent(){
    init = false;
  }
  FillExtent(elev_t elevation0){
    elevation = elevation0;
    init      = false;
  }
  void expand(int x, int y){
    if(!init){
      xmin = xmax = x;
      ymin = ymax = y;
      init = true;
      return;
    }

    if(x<xmin)
      xmin = x;
    else if(x>xmax)
      xmax = x;

    if(y<ymin)
      ymin = y;
    else if(y>ymax)
      ymax = y;
  }
  void expand(const FillExtent &o){
    if(!o.init) return;
    if(!init){
      xmin      = o.xmin;
      xmax      = o.xmax;
      ymin      = o.ymin;
      ymax      = o.ymax;
      elevation = o.elevation;
      init      = true;
    } else {
      if(elevation!=o.elevation){
        std::cerr<<"Problem! Depression fill elevations do not match! "; //TODO
        std::cerr<<elevation<<" "<<o.elevation<<std::endl;
      }
      xmin = std::min(xmin,o.xmin);
      xmax = std::max(xmax,o.xmax);
      ymin = std::min(ymin,o.ymin);
      ymax = std::max(ymax,o.ymax);
    }
  }
};

#endif