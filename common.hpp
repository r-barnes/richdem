#define FLOWDIR_NO_DATA ((uint8_t)255)

//D8 Directions
///x offsets of D8 neighbours, from a central cell
const int dx[9]={0,-1,-1,0,1,1,1,0,-1};  //TODO: These should be merged with my new dinf_d8 to reflect a uniform and intelligent directional system
///y offsets of D8 neighbours, from a central cell
const int dy[9]={0,0,-1,-1,-1,0,1,1,1};
//234
//105
//876
const int d8_inverse[9] = {0,5,6,7,8,1,2,3,4};

//ArcGIS uses:
//32 64 128
//16  0   1
// 8  4   2
const uint8_t d8_arcgis[9] = {0,16,32,64,128,1,2,4,8};

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
  bool operator> (const GridCellZ& a) const { return z>a.z; }
};

typedef int     label_t;
typedef uint8_t flowdir_t;