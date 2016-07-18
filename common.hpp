#ifndef __distpf_common_hpp__
#define __distpf_common_hpp__
#include <sys/time.h>
#include <queue>

//Used for indicating whether a block is on the edge of the larger DEM and which
//edges it is adjacent to
const uint8_t GRID_LEFT   = 1;
const uint8_t GRID_TOP    = 2;
const uint8_t GRID_RIGHT  = 4;
const uint8_t GRID_BOTTOM = 8;

/// Stores the (x,y) coordinates of a grid cell
class GridCell {
 public:
  int x; ///< Grid cell's x-coordinate
  int y; ///< Grid cell's y-coordinate
  /// Initiate the grid cell without coordinates; should generally be avoided
  GridCell(){}
  /// Initiate the grid cell to the coordinates (x0,y0)
  GridCell(int x, int y):x(x),y(y){}
};

template<class elev_t>
class GridCellZ : public GridCell {
 public:
  elev_t z;         ///< Grid cell's z-coordinate
  GridCellZ(int x, int y, elev_t z): GridCell(x,y), z(z) {}
  GridCellZ(){}
  bool operator> (const GridCellZ<elev_t>& a) const { return z>a.z; }
};

template<typename elev_t>
using GridCellZ_pq = std::priority_queue<GridCellZ<elev_t>, std::vector<GridCellZ<elev_t> >, std::greater<GridCellZ<elev_t> > >;


class Timer{
  private:
    timeval start_time;      ///<Last time the timer was started
    double accumulated_time; ///<Accumulated running time since creation
    bool running;            ///<True when the timer is running
    ///Number of seconds between two time objects
    double timediff(timeval beginning, timeval end){
      long seconds, useconds;
      seconds  = end.tv_sec  - beginning.tv_sec;
      useconds = end.tv_usec - beginning.tv_usec;
      return seconds + useconds/1000000.0;
    }
  public:
    Timer(){
      accumulated_time=0;
      running=false;
    }
    void start(){
      if(running)
        throw "Timer was already started!";
      running=true;
      gettimeofday(&start_time, NULL);
    }
    void stop(){
      if(!running)
        throw "Timer was already stopped!";
      running=false;
      timeval end_time;
      gettimeofday(&end_time, NULL);

      accumulated_time+=timediff(start_time,end_time);
    }
    double accumulated(){
      if(running)
        throw "Timer is still running!";
      return accumulated_time;
    }
    double lap(){
      if(!running)
        throw "Timer was not started!";
      timeval lap_time;
      gettimeofday(&lap_time, NULL);
      return timediff(start_time,lap_time);
    }
    void reset(){
      accumulated_time = 0;
      running          = false;
    }
};

#endif