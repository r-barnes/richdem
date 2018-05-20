//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
//#include "catch/catch.hpp"
#include "doctest.h"
#include "richdem/common/Array2D.hpp"

#include <richdem/richdem.hpp>
using namespace richdem;

#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

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
}


TEST_CASE("Checking flow accumulation") {
  for(auto p: fs::directory_iterator("flow_accum")){
    fs::path this_path = p.path();
    if(this_path.extension()==".d8"){
        Array2D<d8_flowdir_t> fds(this_path, false);
        Array2D<int32_t>   correct_ans(this_path.replace_extension("out"), false);
        Array2D<int32_t>   my_ans;
        d8_flow_accum(fds,my_ans);
        REQUIRE( correct_ans == my_ans );
      }
    }
  }  





TEST_CASE("Checking GridCellZk_pq") {
  GridCellZk_pq<int> pq;

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
