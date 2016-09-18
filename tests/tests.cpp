#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch/catch.hpp"
#include "richdem/common/Array2D.hpp"

#include "richdem/methods/d8_methods.hpp"

#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

SCENARIO( "Array2D works", "[Array2D]" ) {

  GIVEN( "A 7x11 Array2D<float>" ) {
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

/*

        WHEN( "the size is increased" ) {
            v.resize( 10 );

            THEN( "the size and capacity change" ) {
                REQUIRE( v.size() == 10 );
                REQUIRE( v.capacity() >= 10 );
            }
        }
        WHEN( "the size is reduced" ) {
            v.resize( 0 );

            THEN( "the size changes but not capacity" ) {
                REQUIRE( v.size() == 0 );
                REQUIRE( v.capacity() >= 5 );
            }
        }
        WHEN( "more capacity is reserved" ) {
            v.reserve( 10 );

            THEN( "the capacity changes but not the size" ) {
                REQUIRE( v.size() == 5 );
                REQUIRE( v.capacity() >= 10 );
            }
        }
        WHEN( "less capacity is reserved" ) {
            v.reserve( 0 );

            THEN( "neither size nor capacity are changed" ) {
                REQUIRE( v.size() == 5 );
                REQUIRE( v.capacity() >= 5 );
            }
        }*/
  }
}


TEST_CASE("Checking flow accumulation", "[FlowAcc]") {
  for(auto p: fs::directory_iterator("flow_accum")){
    fs::path this_path = p.path();
    if(this_path.extension()==".d8"){
      SECTION(p.path().c_str()){
        Array2D<d8_flowdir_t> fds(this_path, false);
        Array2D<int32_t>   correct_ans(this_path.replace_extension("out"), false);
        Array2D<int32_t>   my_ans;
        d8_upslope_area(fds,my_ans);
        REQUIRE( correct_ans == my_ans );
      }
    }
  }  
}

SCENARIO( "Checking flow accumulation", "[FlowAcc]" ) {


    // {
    //     Array2D<flowdir_t> arr()
    //     Array2D<int32_t>   garr
    //   d8_upslope_area
    // }
}