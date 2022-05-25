#include "doctest.h"

#include <richdem/common/constants.hpp>
#include <richdem/common/Array2D.hpp>
#include <richdem/common/loaders.hpp>
#include <richdem/misc/misc_methods.hpp>
#include <richdem/richdem.hpp>
#include <richdem/methods/upslope_cells_functions.hpp>
#include <richdem/methods/upslope_cells.hpp>

#include <filesystem>
#include <queue>

namespace fs = std::filesystem;
using namespace richdem;

TEST_CASE("ManagedVector Construction"){
  auto TestFunc = [](){ManagedVector<int> vec(30,5); return vec;};

  ManagedVector<int>   vec(30,4);
  ManagedVector<float> vecfl(30,4.5);

  auto vec2 = vec;
  auto vec3(vec);
  auto vec4 = TestFunc();

  vec2 = vecfl;
}

TEST_CASE("ManagedVector Resizing"){
  std::vector<int> datavec(30,4);

  ManagedVector<int> mvec(datavec.data(), 30);

  REQUIRE_THROWS(mvec.resize(40));
  REQUIRE_THROWS(mvec.resize(10));
  REQUIRE_NOTHROW(mvec.resize(30));
}

#ifdef USEGDAL
TEST_CASE("Test padding on load") {
  Array2D<int> temp;
  LoadGDAL("ones_block.dem", temp, 1, 1);
  for(int y=0;y<temp.height();y++)
  for(int x=0;x<temp.width();x++){
    if((x==0 || y==0 || x==temp.width()-1 || y==temp.height()-1)){
      REQUIRE(temp(x,y)==0);
    } else {
      REQUIRE(temp(x,y)==1);
    }
  }
}
#endif

TEST_CASE( "Array2D works" ) {

  SUBCASE( "A 7x11 Array2D<float>" ) {
    Array2D<float> arr(7,11);

    REQUIRE( arr.width()  == 7  );
    REQUIRE( arr.height() == 11 );
    REQUIRE( arr.size()   == 77 );

    /*
      0   1  2  3 4   5  6
      7   8  9 10 11 12 13
      14 15 16 17 18 19 20
    */

    int x,y;

    arr.iToxy(0,x,y);
    REQUIRE( x==0 );
    REQUIRE( y==0 );

    arr.iToxy(arr.size()-1,x,y);
    REQUIRE( x==arr.width()-1  );
    REQUIRE( y==arr.height()-1 );

    REQUIRE( arr.xyToI(2,1)==9  );
    REQUIRE( arr.xyToI(6,2)==20 );

    int val = 0;
    for(int y=0;y<arr.height();y++)
    for(int x=0;x<arr.width();x++)
        arr(x,y) = val++;

    REQUIRE(arr(2,1)==9);
    REQUIRE(arr(6,2)==20);
    REQUIRE(arr(0,0)==0);
            }

  SUBCASE("Resize test"){
    Array2D<float> arr;
    arr.resize(3,5,4);
    REQUIRE( arr.width()  == 3  );
    REQUIRE( arr.height() == 5  );
    REQUIRE( arr.size()   == 15 );

    REQUIRE( arr(2,1)==4 );
    REQUIRE( arr(2,2)==4 );
    REQUIRE( arr(0,0)==4 );
            }

  SUBCASE("Copy test"){
    Array2D<float> arr0(3,5,4);
    auto arr1 = arr0;
    auto arr2 = arr0;
    auto arr3 = arr0;
    auto arr4 = arr0;
    auto arr5 = arr0;
            }

  SUBCASE("Clear"){
    Array2D<float> arr0(3,5,2);
    auto arr1 = arr0;
    arr0.clear();
  }

  SUBCASE("Initializer list construction"){
    const Array2D<float> arr = {
      { 1, 2, 3},
      { 4, 5, 6},
      { 7, 8, 9},
      {10,11,12}
    };
    CHECK(arr.width()==3);
    CHECK(arr.height()==4);
    CHECK(arr(0,0)==1);
    CHECK(arr(2,0)==3);
    CHECK(arr(0,2)==7);
    CHECK(arr(2,3)==12);
  }
}


TEST_CASE("Checking flow accumulation") {
  for(auto p: fs::directory_iterator("flow_accum")){
    fs::path this_path = p.path();
    if(this_path.extension()==".d8"){
        Array2D<d8_flowdir_t> fds(this_path.string(), false);
        Array2D<int32_t>   correct_ans(this_path.replace_extension("out").string(), false);
        Array2D<int32_t>   my_ans;
        d8_flow_accum(fds,my_ans);
        REQUIRE( correct_ans == my_ans );
      }
    }
  }





TEST_CASE("Checking GridCellZk_pq") {
  GridCellZk_low_pq<int> pq;

  SUBCASE("Testing Elevation Ordering"){
    pq.emplace(0,0,0);
    pq.emplace(1,0,1);
    pq.emplace(2,0,2);
    pq.emplace(3,0,3);
    REQUIRE(pq.top().x==0); pq.pop();
    REQUIRE(pq.top().x==1); pq.pop();
    REQUIRE(pq.top().x==2); pq.pop();
    REQUIRE(pq.top().x==3); pq.pop();
    REQUIRE(pq.empty()==true);
  }

  SUBCASE("Testing Insertion Ordering"){
    pq.emplace(0,0,0);
    pq.emplace(1,0,0);
    pq.emplace(2,0,0);
    pq.emplace(3,0,0);
    REQUIRE(pq.top().x==0); pq.pop();
    REQUIRE(pq.top().x==1); pq.pop();
    REQUIRE(pq.top().x==2); pq.pop();
    REQUIRE(pq.top().x==3); pq.pop();
    REQUIRE(pq.empty()==true);
  }

  SUBCASE("Testing Mixed Ordering"){
    pq.emplace(0,0,0);
    pq.emplace(1,0,1);
    pq.emplace(2,0,1);
    pq.emplace(3,0,2);
    REQUIRE(pq.top().x==0); pq.pop();
    REQUIRE(pq.top().x==1); pq.pop();
    REQUIRE(pq.top().x==2); pq.pop();
    REQUIRE(pq.top().x==3); pq.pop();
    REQUIRE(pq.empty()==true);
  }
}

TEST_CASE("Checking GridCellZk_pq") {
  GridCellZk_high_pq<int> pq;

  SUBCASE("Testing Elevation Ordering"){
    pq.emplace(0,0,0);
    pq.emplace(1,0,1);
    pq.emplace(2,0,2);
    pq.emplace(3,0,3);
    REQUIRE(pq.top().x==0); pq.pop();
    REQUIRE(pq.top().x==1); pq.pop();
    REQUIRE(pq.top().x==2); pq.pop();
    REQUIRE(pq.top().x==3); pq.pop();
    REQUIRE(pq.empty()==true);
  }

  SUBCASE("Testing Insertion Ordering"){
    pq.emplace(0,0,0);
    pq.emplace(1,0,0);
    pq.emplace(2,0,0);
    pq.emplace(3,0,0);
    REQUIRE(pq.top().x==3); pq.pop();
    REQUIRE(pq.top().x==2); pq.pop();
    REQUIRE(pq.top().x==1); pq.pop();
    REQUIRE(pq.top().x==0); pq.pop();
    REQUIRE(pq.empty()==true);
  }

  SUBCASE("Testing Mixed Ordering"){
    pq.emplace(0,0,0);
    pq.emplace(1,0,1);
    pq.emplace(2,0,1);
    pq.emplace(3,0,2);
    REQUIRE(pq.top().x==0); pq.pop();
    REQUIRE(pq.top().x==2); pq.pop();
    REQUIRE(pq.top().x==1); pq.pop();
    REQUIRE(pq.top().x==3); pq.pop();
    REQUIRE(pq.empty()==true);
  }
}


TEST_CASE("Checking depression filling") {
  RDLOG_DEBUG<<"About to load depressions/testdem1.dem";
  Array2D<int> elevation_orig("depressions/testdem1.dem");
  RDLOG_DEBUG<<"Loaded depressions/testdem1.dem";

  SUBCASE("PriorityFlood_Original"){
    auto elevation = elevation_orig;
    PriorityFlood_Original<Topology::D8>(elevation);
    Array2D<int> manually_flooded("depressions/testdem1.all.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("PriorityFlood_Barnes2014"){
    auto elevation = elevation_orig;
    PriorityFlood_Barnes2014<Topology::D8>(elevation);
    Array2D<int> manually_flooded("depressions/testdem1.all.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("Zhou2016"){
    auto elevation = elevation_orig;
    PriorityFlood_Zhou2016(elevation);
    Array2D<int> manually_flooded("depressions/testdem1.all.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("Wei2018"){
    auto elevation = elevation_orig;
    PriorityFlood_Wei2018(elevation);
    Array2D<int> manually_flooded("depressions/testdem1.all.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("FillDepressions"){
    auto elevation = elevation_orig;
    FillDepressions<Topology::D8>(elevation);
    Array2D<int> manually_flooded("depressions/testdem1.all.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("PriorityFlood_Barnes2014_max_dep"){
    auto elevation = elevation_orig;
    PriorityFlood_Barnes2014_max_dep<Topology::D8>(elevation,1);
    elevation.printAll();
    Array2D<int> manually_flooded("depressions/testdem1.1.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("PriorityFlood_Barnes2014_max_dep"){
    auto elevation = elevation_orig;
    PriorityFlood_Barnes2014_max_dep<Topology::D8>(elevation,2);
    elevation.printAll();
    Array2D<int> manually_flooded("depressions/testdem1.2.out");
    REQUIRE(elevation==manually_flooded);
  }

}



TEST_CASE("Checking depression breaching") {
  RDLOG_DEBUG<<"About to load depressions/testdem1.dem";
  Array2D<int> elevation_orig("breaching/testdem1.dem");
  RDLOG_DEBUG<<"Loaded breaching/testdem1.dem";

  SUBCASE("Lindsay2016 Complete Breaching"){
    auto elevation = elevation_orig;
    CompleteBreaching_Lindsay2016<Topology::D8>(elevation);
    elevation.printAll();
    Array2D<int> manually_flooded("breaching/testdem1.complete.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("Lindsay2016 Selective Breaching (Length=2, Depth=9e99)"){
    auto elevation = elevation_orig;
    Lindsay2016(elevation,LindsayMode::SELECTIVE_BREACHING,false,false,2,9999);
    elevation.printAll();
    Array2D<int> manually_flooded("breaching/testdem1.selective-len2-depth9999.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("Lindsay2016 Selective Breaching (Length=4, Depth=9e99)"){
    auto elevation = elevation_orig;
    Lindsay2016(elevation,LindsayMode::SELECTIVE_BREACHING,false,false,4,9999);
    elevation.printAll();
    Array2D<int> manually_flooded("breaching/testdem1.selective-len4-depth9999.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("Lindsay2016 Selective Breaching (Length=4, Depth=2)"){
    auto elevation = elevation_orig;
    Lindsay2016(elevation,LindsayMode::SELECTIVE_BREACHING,false,false,4,2);
    elevation.printAll();
    Array2D<int> manually_flooded("breaching/testdem1.selective-len4-depth2.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("Lindsay2016 Selective Breaching (Length=4, Depth=2, Fill Depressions)"){
    auto elevation = elevation_orig;
    Lindsay2016(elevation,LindsayMode::SELECTIVE_BREACHING,false,true,4,2);
    elevation.printAll();
    Array2D<int> manually_flooded("breaching/testdem1.selective-len4-depth2-filldep.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("Lindsay2016 Selective Breaching (Length=4, Depth=8)"){
    auto elevation = elevation_orig;
    Lindsay2016(elevation,LindsayMode::SELECTIVE_BREACHING,false,false,4,8);
    elevation.printAll();
    Array2D<int> manually_flooded("breaching/testdem1.selective-len4-depth8.out");
    REQUIRE(elevation==manually_flooded);
  }

  SUBCASE("Lindsay2016 Constrained Breaching (Length=4, Depth=3)"){
    auto elevation = elevation_orig;
    Lindsay2016(elevation,LindsayMode::CONSTRAINED_BREACHING,false,false,4,3);
    elevation.printAll();
    Array2D<int> manually_flooded("breaching/testdem1.constrained-len4-depth3.out");
    REQUIRE(elevation==manually_flooded);
  }

}


TEST_CASE("Checking flow accumulation") {
  Array2D<float> beauford("beauford/beauford.tif");
  PriorityFlood_Wei2018(beauford);

  SUBCASE("FA_Tarboton")            {Array2D<double> accum(beauford); FA_Tarboton            (beauford, accum); }
  SUBCASE("FA_Dinfinity")           {Array2D<double> accum(beauford); FA_Dinfinity           (beauford, accum); }
  SUBCASE("FA_Holmgren")            {Array2D<double> accum(beauford); FA_Holmgren            (beauford, accum, 0.8); }
  SUBCASE("FA_Quinn")               {Array2D<double> accum(beauford); FA_Quinn               (beauford, accum); }
  SUBCASE("FA_Freeman")             {Array2D<double> accum(beauford); FA_Freeman             (beauford, accum, 0.8); }
  SUBCASE("FA_FairfieldLeymarieD8") {Array2D<double> accum(beauford); FA_FairfieldLeymarieD8 (beauford, accum); }
  SUBCASE("FA_FairfieldLeymarieD4") {Array2D<double> accum(beauford); FA_FairfieldLeymarieD4 (beauford, accum); }
  SUBCASE("FA_Rho8")                {Array2D<double> accum(beauford); FA_Rho8                (beauford, accum); }
  SUBCASE("FA_Rho4")                {Array2D<double> accum(beauford); FA_Rho4                (beauford, accum); }
  SUBCASE("FA_D8")                  {Array2D<double> accum(beauford); FA_D8                  (beauford, accum); }
  SUBCASE("FA_D4")                  {Array2D<double> accum(beauford); FA_D4                  (beauford, accum); }

}



TEST_CASE("BucketFill"){
  SUBCASE("Middle"){
    Array2D<int> check_raster = {
      {1, 1, 1, 1, 1, 1, 1},
      {1, 1, 2, 1, 1, 2, 1},
      {1, 1, 2, 1, 2, 1, 1},
      {1, 1, 2, 2, 1, 1, 1},
      {1, 1, 1, 1, 1, 1, 1},
    };

    Array2D<int> set_raster(check_raster.width(), check_raster.height(), 0);

    std::vector<size_t> q;
    q.emplace_back(set_raster.xyToI(2,1));

    BucketFill<Topology::D8>(check_raster, set_raster, 2, 4, q);

    const Array2D<int> good_raster = {
      {0, 0, 0, 0, 0, 0, 0},
      {0, 0, 4, 0, 0, 4, 0},
      {0, 0, 4, 0, 4, 0, 0},
      {0, 0, 4, 4, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0},
    };

    CHECK(set_raster==good_raster);
  }

  SUBCASE("Edges"){
    Array2D<int> check_raster = {
      {1, 1, 1, 1, 1, 2, 1},
      {1, 1, 2, 1, 2, 1, 1},
      {1, 1, 2, 1, 1, 2, 1},
      {1, 1, 2, 2, 1, 1, 1},
      {1, 2, 1, 1, 1, 1, 1},
    };

    Array2D<int> set_raster(check_raster.width(), check_raster.height(), 0);

    BucketFillFromEdges<Topology::D8>(check_raster, set_raster, 2, 4);

    const Array2D<int> good_raster = {
      {0, 0, 0, 0, 0, 4, 0},
      {0, 0, 4, 0, 4, 0, 0},
      {0, 0, 4, 0, 0, 4, 0},
      {0, 0, 4, 4, 0, 0, 0},
      {0, 4, 0, 0, 0, 0, 0},
    };

    CHECK(set_raster==good_raster);
  }
}

TEST_CASE("DirectionsMatchExpectations"){
  SUBCASE("D8"){
    const auto dx_func = get_dx_for_topology<Topology::D8>();
    const auto dx_good = {0, -1, -1,  0,  1, 1, 1, 0, -1};
    CHECK(std::equal(dx_func.begin(), dx_func.end(), dx_good.begin()));

    const auto dy_func = get_dy_for_topology<Topology::D8>();
    const auto dy_good = {0,  0, -1, -1, -1, 0, 1, 1,  1};
    CHECK(std::equal(dy_func.begin(), dy_func.end(), dy_good.begin()));

    CHECK(!std::equal(dx_func.begin(), dx_func.end(), dy_func.begin()));

    CHECK(get_nmax_for_topology<Topology::D8>() == 8);
  }

  SUBCASE("D4"){
    const auto dx_func = get_dx_for_topology<Topology::D4>();
    const auto dx_good =  {0, -1,  0, 1, 0};
    CHECK(std::equal(dx_func.begin(), dx_func.end(), dx_good.begin()));

    const auto dy_func = get_dy_for_topology<Topology::D4>();
    const auto dy_good = {0,  0, -1, 0, 1};
    CHECK(std::equal(dy_func.begin(), dy_func.end(), dy_good.begin()));

    CHECK(!std::equal(dx_func.begin(), dx_func.end(), dy_func.begin()));

    CHECK(get_nmax_for_topology<Topology::D4>() == 4);
  }
}
TEST_CASE("Checking Catchment Delineation Generic"){
  SUBCASE("mask"){
    SUBCASE("single"){
    Array2D<uint8_t> mask = {
      {0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0},
      {0, 0, 1, 0, 0, 0},
      {0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0},
    };
    mask.setNoData(0);
    Array2D<uint8_t> upslope(mask,0);

    std::queue<GridCell> q_expected;
    q_expected.emplace(2,2);
    const Array2D<uint8_t> upslope_e = {
      {0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0},
      {0, 0, 2, 0, 0, 0},
      {0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0},
    };

    std::queue<GridCell> q = queue_from_mask(mask, upslope);
    std::cout<< mask.noData() << std::endl;
    CHECK(q.front().x == q_expected.front().x);
    CHECK(q.front().y == q_expected.front().y);
    CHECK(q.size() == 1);
    CHECK(upslope(2,2) == 2);
    }
    SUBCASE("multi"){
    Array2D<uint8_t> mask = {
      {0, 0, 0, 0, 0, 0},
      {0, 1, 0, 0, 0, 0},
      {0, 0, 1, 0, 0, 0},
      {0, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0},
    };
    mask.setNoData(0);
    Array2D<uint8_t> upslope(mask,0);

    std::queue<GridCell> q_expected;
    q_expected.emplace(1,1);
    q_expected.emplace(1,3);
    q_expected.emplace(2,2);

    const Array2D<uint8_t> upslope_e = {
      {0, 0, 0, 0, 0, 0},
      {0, 2, 0, 0, 0, 0},
      {0, 0, 2, 0, 0, 0},
      {0, 2, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0},
    };

    std::queue q = queue_from_mask(mask, upslope);
    CHECK(q.size() == 3);
    CHECK(q == q_expected);
    CHECK(upslope == upslope_e);
    }
    // TODO: make more tests?
  }

  SUBCASE("line"){
    SUBCASE("horizontal"){
      int x0=0;
      int y0=0;
      int x1=2;
      int y1=0;

      std::queue<GridCell> q_expected;
      q_expected.emplace(0,0);
      q_expected.emplace(1,0);
      q_expected.emplace(2,0);

      Array2D<uint8_t> upslope_e = {
        {2, 2, 2},
        {0, 0, 0},
        {0, 0, 0},
      };
      upslope_e.setNoData(0);
      Array2D<uint8_t> upslope(upslope_e,0);
      upslope.setNoData(0);

      std::queue<GridCell> q = queue_from_linepoints(x0,y0,x1,y1,upslope);

      CHECK(upslope == upslope_e);
      CHECK(q == q_expected);
    }
    //TODO: add more tests. definately for slope>1, so that swapping is
    //routinely tested. (I now did it by inspection, and it works)
  }

  Array2D<int> dem = {
    {1,2,3,4,5},
    {1,2,3,4,5},
    {1,2,3,4,5},
    {1,2,3,4,5},
    {1,2,3,4,5},
  };
  dem.setNoData(0);

  SUBCASE("mflow"){
    
    Array2D<uint8_t> result(dem,0);
    std::queue<GridCell> expansion;
    expansion.emplace(0,2);

    Array2D<uint8_t> expected = {
      {0, 0, 1, 1, 1},
      {0, 1, 1, 1, 1},
      {0, 1, 1, 1, 1},
      {0, 1, 1, 1, 1},
      {0, 0, 1, 1, 1}
    };
    result.setNoData(0);
    expected.setNoData(0); //the function sets nodata to 0

    upslope_cells_mflow<Topology::D8>(expansion, dem, result);
    CHECK(result == expected);
  }

  SUBCASE("props (take quinn)"){
    //more methods is not necessary
    Array2D<uint8_t> result(dem,0);
    result.setNoData(0);
    std::queue<GridCell> expansion;
    expansion.emplace(0,2);

    //construct the proportions array
    Array3D<float> props(dem);
    FM_Quinn(dem,props);

    Array2D<uint8_t> expected = {
      {0, 0, 0, 0, 0},
      {0, 1, 1, 1, 0},
      {0, 1, 1, 1, 0},
      {0, 1, 1, 1, 0},
      {0, 0, 0, 0, 0}      // TODO: maybe add edge cells if they are higher? is slightly complex though, cuz Array3D does not support that...
    };                     // so flow always goes off the edge of the grid, hence the 0s.
    expected.setNoData(0); //the function sets nodata to 0

    upslope_cells_props<Topology::D8>(expansion, props, result);
    CHECK(result == expected);
  }
}

TEST_CASE("Checking upslope cells") {
  for(auto p: fs::directory_iterator("upslope_cells")){
    fs::path this_path = p.path();
    if(this_path.extension()==".dem"){
        Array2D<int>     dem      (this_path,                            false);
        Array2D<uint8_t> mask     (this_path.replace_extension("mask"),  false);
        Array2D<int32_t> ans_mflow(this_path.replace_extension("mflow"), false);
        Array2D<int32_t> ans_quinn(this_path.replace_extension("quinn"), false);
        Array2D<int32_t> ans_d8   (this_path.replace_extension("d8"),    false);
        Array2D<int32_t> result;
        Array3D<float>   props(dem);

        UC_mask_mflow(mask,dem,result);
        CHECK(result == ans_mflow);
        result.printAll("mflow result");
        ans_mflow.printAll("expected result");

        FM_Quinn(dem, props);
        UC_mask_props(mask,props,result);
        CHECK(result == ans_quinn);
        result.printAll("Quinn upslope result");
        ans_quinn.printAll("expected result");

        FM_D8(dem, props);
        UC_mask_props(mask,props,result);
        CHECK(result == ans_d8);
        result.printAll("D8 upslope result");
        ans_d8.printAll("expected result");


      }
    }
  }
