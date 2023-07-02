#include "doctest.h"

#include <richdem/depressions/fill_spill_merge.hpp>
#include <richdem/terrain_generation.hpp>

#include <random>
#include <sstream>

using namespace richdem;
using namespace richdem::dephier;

#ifdef CODE_COVERAGE
  #pragma message "FSM is using a small number of test cases for code coverage estimation. Disable code coverage to enable more extensive tests."
  const int number_of_small_tests = 50;
  const int number_of_large_tests = 5;
#else
  #pragma message "FSM is using a large number of test cases to judge correctness. Enabling code coverage will reduce the number of test cases used."
  const int number_of_small_tests = 6000;
  const int number_of_large_tests = 500;
#endif

template<class T>
double MaxArrayDiff(const Array2D<T> &a, const Array2D<T> &b){
  double max_diff = 0;
  for(auto i=a.i0();i<a.size();i++){
    max_diff = std::max(max_diff, (double)std::abs(a(i)-b(i)));
  }
  return max_diff;
}

template<class T>
bool ArrayValuesEqual(const Array2D<T> &a, const Array2D<T> &b){
  for(auto i=a.i0();i<a.size();i++){
    if(a(i)!=b(i))
      return false;
  }
  return true;
}

template<class T>
bool ArrayValuesAllEqual(const Array2D<T> &a, const T val){
  for(auto i=a.i0();i<a.size();i++){
    if(a(i)!=val)
      return false;
  }
  return true;
}



std::mt19937_64 gen;

Array2D<double> random_terrain(std::mt19937_64 &gen, const int min_size, const int max_size){
  static std::uniform_int_distribution<uint32_t> seed_dist;

  std::uniform_int_distribution<int> size_dist(min_size, max_size);

  const auto size = size_dist(gen);
  Array2D<double> terrain(size, size);
  generate_perlin_terrain(terrain, seed_dist(gen));
  return terrain;
}

Array2D<double> random_integer_terrain(std::mt19937_64 &gen, const int min_size, const int max_size){
  static std::uniform_int_distribution<uint32_t> seed_dist;

  std::uniform_int_distribution<int> size_dist(min_size, max_size);

  const auto size = size_dist(gen);
  Array2D<double> dem(size, size);
  generate_perlin_terrain(dem, seed_dist(gen));
  for(auto i=dem.i0();i<dem.size();i++){
    dem(i) *= 100;
    dem(i) = static_cast<int>(dem(i));
  }

  return dem;
}


TEST_CASE("Depression volume"){
  CHECK_EQ(DepressionVolume(2, 5, 10), 0);
  CHECK_EQ(DepressionVolume(3, 5, 10), 5);
  CHECK_EQ(DepressionVolume(4, 5, 10), 10);
}

TEST_CASE("Determine water level"){
  SUBCASE("Depression volume exactly equals water volume"){
    double sill_wtd = -2;
    const auto water_level = DetermineWaterLevel(sill_wtd, 10, 4, 5, 10);
    CHECK_EQ(sill_wtd, -2);
    CHECK_EQ(water_level, 4); //Water elevation equals the sill elevation
  }

  SUBCASE("Water volume is less than the depression volume"){
    double sill_wtd = -2;
    const auto water_level = DetermineWaterLevel(sill_wtd, 8, 4, 5, 10);
    CHECK_EQ(sill_wtd, -2);
    CHECK_EQ(water_level, 18/5.0); //Water elevation equals the sill elevation
  }

  SUBCASE("Water volume is greater than the depression volume"){
    double sill_wtd = -2;
    const auto water_level = DetermineWaterLevel(sill_wtd, 12, 4, 5, 10);
    CHECK_EQ(sill_wtd, 0);
    //Water elevation equals the sill elevation since the sill absorbs all
    //excess
    CHECK_EQ(water_level, 4);
  }
}



TEST_CASE("MoveWaterIntoPits 1"){
  const Array2D<double> topo = {
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
      {-9,  9,  9,  9,  9,  9,  9,  9,  9, -9},
      {-9,  9,  8,  8,  8,  8,  7,  6,  9, -9},
      {-9,  9,  8,  7,  8,  7,  6,  5,  9, -9},
      {-9,  9,  8,  7,  8,  6,  5,  4,  9, -9},
      {-9,  9,  8,  8,  8,  5,  4,  3,  9, -9},
      {-9,  9,  7,  6,  5,  4,  3,  2,  9, -9},
      {-9,  9,  7,  6,  5,  4,  3,  1,  9, -9},
      {-9,  9,  9,  9,  9,  9,  9,  9,  9, -9},
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
  };

  Array2D<dh_label_t> label(topo.width(), topo.height(), NO_DEP);
  label.setEdges(OCEAN);

  Array2D<flowdir_t> flowdirs(topo.width(), topo.height(), NO_FLOW);
  Array2D<double> wtd(topo.width(), topo.height(), 0);

  auto DH = GetDepressionHierarchy<double,Topology::D8>(topo, label, flowdirs);

  const Array2D<dh_label_t> label_good = {
    {0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0},
    {0,0,2,2,1,1,1,1,0,0},
    {0,0,2,2,1,1,1,1,0,0},
    {0,0,2,2,1,1,1,1,0,0},
    {0,0,1,1,1,1,1,1,0,0},
    {0,0,1,1,1,1,1,1,0,0},
    {0,0,1,1,1,1,1,1,0,0},
    {0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0}
  };

  const Array2D<flowdir_t> flowdirs_good = {
    {0,0,0,0,0,0,0,0,0,0},
    {0,8,4,4,4,4,4,4,6,0},
    {0,8,6,7,6,6,6,7,6,0},
    {0,8,6,7,6,6,6,7,6,0},
    {0,8,5,0,6,6,6,7,6,0},
    {0,8,6,6,6,6,6,7,6,0},
    {0,8,6,5,6,5,6,7,6,0},
    {0,8,5,4,5,4,5,0,6,0},
    {0,6,6,6,6,6,6,6,6,0},
    {0,0,0,0,0,0,0,0,0,0}
  };

  CHECK(ArrayValuesEqual(label,label_good));
  CHECK(ArrayValuesEqual(flowdirs,flowdirs_good));

  wtd.setAll(1);

  MoveWaterIntoPits<double, double>(topo, label, flowdirs, DH, wtd);

  CHECK(ArrayValuesAllEqual(wtd,0.0));

  CHECK_EQ(DH.at(0).water_vol, 64);
  CHECK_EQ(DH.at(1).water_vol, 30);
  CHECK_EQ(DH.at(2).water_vol,  6);

  CHECK_EQ(DH.at(0).parent, NO_VALUE);
  CHECK_EQ(DH.at(1).parent, 3);
  CHECK_EQ(DH.at(2).parent, 3);
  CHECK_EQ(DH.at(3).parent, 0);

  CHECK(std::isnan(DH.at(0).dep_vol));
  CHECK_EQ(DH.at(1).dep_vol,  73);
  CHECK_EQ(DH.at(2).dep_vol,   2);
  CHECK_EQ(DH.at(3).dep_vol, 111);
}
//TODO: Add second test case with more tests and clearer outlets


void MoveWaterIntoPitsRepeatedly(const int count, const int min_size, const int max_size){
  #pragma omp parallel for
  for(int i=0;i<count;i++){
    std::stringstream oss;

    Array2D<double> dem;

    #pragma omp critical
    {
      oss<<gen;
      dem = random_integer_terrain(gen, min_size, max_size);
      std::cerr<<"MoveWaterIntoPits Repeatedly #"<<i<<std::endl;
    }

    Array2D<dh_label_t> labels  (dem.width(), dem.height(), NO_DEP );
    Array2D<flowdir_t>  flowdirs(dem.width(), dem.height(), NO_FLOW);

    dem.setEdges(-1);
    labels.setEdges(OCEAN);

    const auto labels_orig = labels;

    auto deps1 = GetDepressionHierarchy<double,Topology::D8>(dem, labels, flowdirs);

    Array2D<double> wtd(dem.width(), dem.height(), 1);

    MoveWaterIntoPits<double,double>(dem, labels, flowdirs, deps1, wtd);

    labels = labels_orig;
    auto deps2 = GetDepressionHierarchy<double,Topology::D8>(dem, labels, flowdirs);

    wtd.setAll(0);

    for(const auto &dep: deps1){
      if(dep.water_vol>0 && dep.pit_cell!=NO_VALUE){
        wtd(dep.pit_cell) = dep.water_vol;
      }
    }

    MoveWaterIntoPits<double,double>(dem, labels, flowdirs, deps2, wtd);

    CHECK_EQ(deps1.size(), deps2.size());
    for(size_t i=1;i<deps1.size();i++){
      CHECK_EQ(deps1.at(i).water_vol, deps2.at(i).water_vol);
    }
  }
}

TEST_CASE("MoveWaterIntoPits Repeatedly"){
  MoveWaterIntoPitsRepeatedly(number_of_small_tests,  10,  30);
  MoveWaterIntoPitsRepeatedly(number_of_large_tests, 100, 300);
}



TEST_CASE("Backfill Depression"){
  const Array2D<double> topo = {
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
      {-9,  6,  6,  6,  6,  6,  6,  6,  6, -9},
      {-9,  6,  1,  6,  1,  1,  6,  1,  6, -9},
      {-9,  6,  1,  6,  1,  3,  6,  1,  6, -9},
      {-9,  6,  1,  6,  2,  1,  4,  1,  6, -9},
      {-9,  6,  1,  6,  1,  1,  6,  1,  6, -9},
      {-9,  6,  1,  6,  6,  6,  6,  1,  6, -9},
      {-9,  6,  1,  1,  1,  1,  1,  1,  6, -9},
      {-9,  6,  6,  6,  6,  6,  6,  6,  6, -9},
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
  };

  Array2D<double> wtd(topo.width(), topo.height());

  std::vector<flat_c_idx> cells_affected = {
    topo.xyToI(4,2), topo.xyToI(5,2),
    topo.xyToI(4,3), topo.xyToI(5,3),
    topo.xyToI(4,4), topo.xyToI(5,4),
    topo.xyToI(4,5), topo.xyToI(5,5),
  };

  BackfillDepression(4.0, topo, wtd, cells_affected);

  const Array2D<double> wtd_good = {
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
      { 0,  0,  0,  0,  3,  1,  0,  0,  0,  0},
      { 0,  0,  0,  0,  2,  3,  0,  0,  0,  0},
      { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  };

  CHECK(ArrayValuesEqual(wtd,wtd_good));
}



TEST_CASE("FillDepressions"){
  const Array2D<double> topo = {
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
      {-9,  6,  6,  6,  6,  6,  6,  6,  6, -9},
      {-9,  6,  1,  6,  1,  1,  6,  1,  6, -9},
      {-9,  6,  1,  6,  3,  3,  6,  1,  6, -9},
      {-9,  6,  1,  6,  2,  1,  4,  1,  6, -9},
      {-9,  6,  1,  6,  1,  1,  6,  1,  6, -9},
      {-9,  6,  1,  6,  6,  6,  6,  1,  6, -9},
      {-9,  6,  1,  1,  1,  1,  1,  1,  6, -9},
      {-9,  6,  6,  6,  6,  6,  6,  6,  6, -9},
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
  };

  const Array2D<dh_label_t> label = {
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  2,  2,  2,  2,  1,  1,  0},
      { 0,  1,  1,  1,  1,  1,  1,  1,  1,  0},
      { 0,  1,  1,  1,  1,  1,  1,  1,  1,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  };

  Array2D<double> wtd(topo.width(), topo.height(), 0.0);

  const auto pit_cell = topo.xyToI(4,2);
  const auto out_cell = topo.xyToI(4,3);

  REQUIRE(topo(pit_cell)==1);

  const std::unordered_set<dh_label_t> dep_labels = {2};

  SUBCASE("No water to add"){
    wtd(4,3) = -0.5;
    const auto wtd_good = wtd;
    const double water_vol = 0;
    FillDepressions(pit_cell, out_cell, dep_labels, water_vol, topo, label, wtd);
    CHECK_EQ(wtd, wtd_good);
  }

  SUBCASE("Standard Case"){
    wtd.setAll(0);
    const double water_vol = 3.0;
    FillDepressions(pit_cell, out_cell, dep_labels, water_vol, topo, label, wtd);

    const auto W = 1.5;
    const Array2D<double> wtd_good = {
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  W,  W,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    };

    CHECK(MaxArrayDiff(wtd,wtd_good)<1e-6);
  }

  SUBCASE("Sill Cell Aborbs some water"){
    wtd(4,3) = -1;
    const double water_vol = 5.0;
    FillDepressions(pit_cell, out_cell, dep_labels, water_vol, topo, label, wtd);

    const auto W = 2.0;
    const Array2D<double> wtd_good = {
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  W,  W,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    };

    CHECK(MaxArrayDiff(wtd,wtd_good)<1e-6);
  }

  SUBCASE("Passes over a saddle"){
    const double water_vol = 19.0;
    const auto out_cell = topo.xyToI(6,4);
    FillDepressions(pit_cell, out_cell, dep_labels, water_vol, topo, label, wtd);

    const Array2D<double> wtd_good = {
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
        { 0,  0,  0,  0,  1,  1,  0,  0,  0,  0},
        { 0,  0,  0,  0,  2,  3,  0,  0,  0,  0},
        { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    };

    CHECK(MaxArrayDiff(wtd,wtd_good)<1e-6);
  }

  SUBCASE("Passes over a saddle and sill absorbs some"){
    const double water_vol = 19.5;
    const auto out_cell = topo.xyToI(6,4);
    wtd(6,4) = -1;
    FillDepressions(pit_cell, out_cell, dep_labels, water_vol, topo, label, wtd);

    const Array2D<double> wtd_good = {
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
        { 0,  0,  0,  0,  1,  1,  0,  0,  0,  0},
        { 0,  0,  0,  0,  2,  3,-.5,  0,  0,  0},
        { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    };

    CHECK(MaxArrayDiff(wtd,wtd_good)<1e-6);
  }
}



void RandomizedHeavyFloodingVsPriorityFlood(const int count, const int min_size, const int max_size){
  #pragma omp parallel for
  for(int i=0;i<count;i++){
    std::stringstream oss;

    Array2D<double> dem;

    #pragma omp critical
    {
      oss<<gen;
      dem = random_terrain(gen, min_size, max_size);
      std::cerr<<"Randomized Heavy Flooding vs Priority-Flood #"<<i<<std::endl;
    }

    Array2D<dh_label_t> labels  (dem.width(), dem.height(), NO_DEP );
    Array2D<flowdir_t>  flowdirs(dem.width(), dem.height(), NO_FLOW);

    dem.setEdges(-1);
    labels.setEdges(OCEAN);

    auto deps = GetDepressionHierarchy<double,Topology::D8>(dem, labels, flowdirs);

    //wtd with a *lot* of initial surface water
    Array2D<double> wtd(dem.width(), dem.height(), 100);

    try {
      FillSpillMerge(dem, labels, flowdirs, deps, wtd);
    } catch (const std::exception &e) {
      std::cerr<<"FillSpillMerge failed because of \""<<e.what()<<"\" with width = "<<dem.width()<<" height = "<<dem.height()<<" state = "<<oss.str()<<std::endl;
      throw e;
    }

    for(auto i=dem.i0(); i<dem.size(); i++){
      if(!dem.isNoData(i))
        dem(i) += wtd(i);
    }

    auto comparison_dem = dem;
    PriorityFlood_Zhou2016(comparison_dem);

    CHECK_MESSAGE(
      MaxArrayDiff(comparison_dem,dem)<1e-6,
      "Randomized Heavy Flooding vs Priority-Flood failed with width = ", dem.width(), " height = ", dem.height(), " state = ", oss.str()
    );
  }
}

TEST_CASE("Randomized Heavy Flooding vs Priority-Flood"){
  RandomizedHeavyFloodingVsPriorityFlood(number_of_small_tests,  10,  30);
  RandomizedHeavyFloodingVsPriorityFlood(number_of_large_tests, 100, 300);
}



void RandomizedTestingOfRepeatedFSM(const int count, const int min_size, const int max_size){
  #pragma omp parallel for
  for(int i=0;i<count;i++){
    std::stringstream oss;

    Array2D<double> dem;
    #pragma omp critical
    {
      oss<<gen;
      dem = random_integer_terrain(gen, min_size, max_size);
      std::cerr<<"Randomized Testing of Repeated FSM #"<<i<<std::endl;
    }

    Array2D<dh_label_t> label   (dem.width(), dem.height(), NO_DEP);
    Array2D<flowdir_t>  flowdirs(dem.width(), dem.height(), NO_FLOW);
    Array2D<double>     wtd     (dem.width(), dem.height(), 0);

    auto do_fsm = [&](){
      //Make sure the edges are identifiable as an ocean
      label.setAll(NO_DEP);
      dem.setEdges(-1);
      label.setEdges(OCEAN);

      auto DH = GetDepressionHierarchy<double,Topology::D8>(dem, label, flowdirs);
      try {
        FillSpillMerge(dem, label, flowdirs, DH, wtd);
      } catch (const std::exception &e) {
        std::cerr<<"FillSpillMerge failed because of \""<<e.what()<<"\" with width = "<<dem.width()<<" height = "<<dem.height()<<" state = "<<oss.str()<<std::endl;
        throw e;
      }
    };

    //Initially distribute the water
    wtd.setAll(1);
    do_fsm();

    const auto first_wtd = wtd;

    //Distribute it a second time
    do_fsm();

    CHECK_MESSAGE(
      MaxArrayDiff(first_wtd,wtd)<1e-6,
      "Randomized Testing of Repeated FSM failed with width = ", dem.width(), " height = ", dem.height(), " state = ", oss.str()
    );
  }
}

TEST_CASE("Randomized Testing of Repeated FSM"){
  RandomizedTestingOfRepeatedFSM(number_of_small_tests,  10,  30);
  RandomizedTestingOfRepeatedFSM(number_of_large_tests, 100, 300);
}



TEST_CASE("PQ Issue"){
  const Array2D<double> topo = {
      {-9, -9, -9, -9, -9, -9, -9, 9, -9},
      {-9,  5,  5,  5,  5,  5,  5, 5, -9},
      {-9,  5,  1,  1,  5,  1,  1, 5, -9},
      {-9,  5,  1,  1,  5,  1,  1, 5, -9},
      {-9,  2,  1,  1,  3,  1,  1, 5, -9},
      {-9,  5,  1,  1,  5,  1,  1, 5, -9},
      {-9,  5,  5,  5,  5,  5,  5, 5, -9},
      {-9, -9, -9, -9, -9, -9, -9, 9, -9}
  };

  Array2D<dh_label_t> labels  (topo.width(), topo.height(), NO_DEP );
  Array2D<flowdir_t>  flowdirs(topo.width(), topo.height(), NO_FLOW);

  labels.setEdges(OCEAN);

  auto deps = GetDepressionHierarchy<double, Topology::D8>(topo, labels, flowdirs);

  flowdirs.printAll("Flowdirs", 1, 0);

  Array2D<double> wtd(topo.width(), topo.height(), 0);

  wtd(5,5) = 17;

  FillSpillMerge(topo, labels, flowdirs, deps, wtd);
}



TEST_CASE("PQ Issue 2"){
  const Array2D<double> topo = {
      {-9, -9,  -9, -9, -9},
      {-9,  0,   0,  0, -9},
      {-9,  0,  -8,  0, -9},
      {-9,  0,  -8,  0, -9},
      {-9,  0,   0,  0, -9},
      {-9, -9,  -9, -9, -9}
  };

  const Array2D<dh_label_t> labels = {
      {0,  0,   0,  0,  0},
      {0,  0,   0,  0,  0},
      {0,  0,   1,  0,  0},
      {0,  0,   1,  0,  0},
      {0,  0,   0,  0,  0},
      {0,  0,   0,  0,  0}
  };

  const auto pit_cell = topo.xyToI(2,2);
  const auto out_cell = topo.xyToI(3,3);

  const std::unordered_set<dh_label_t> dep_labels = {1};

  double water_vol = 2.0;

  Array2D<double> wtd(topo.width(), topo.height(), 0.0);

  FillDepressions(pit_cell, out_cell, dep_labels, water_vol, topo, labels, wtd);
}



void RandomizedIncrementalVsBigDump(const int count, const int min_size, const int max_size){
  #pragma omp parallel for
  for(int i=0;i<count;i++){
    std::stringstream oss;

    Array2D<double> dem;
    #pragma omp critical
    {
      oss<<gen;
      dem = random_integer_terrain(gen, min_size, max_size);
      std::cerr<<"Randomized Incremental vs Big Dump #"<<i<<std::endl;
    }

    Array2D<dh_label_t> label   (dem.width(), dem.height(), NO_DEP);
    Array2D<flowdir_t>  flowdirs(dem.width(), dem.height(), NO_FLOW);
    Array2D<double>     wtd     (dem.width(), dem.height(), 0);

    //Make sure the edges are identifiable as an ocean
    label.setAll(NO_DEP);
    dem.setEdges(-1);
    label.setEdges(OCEAN);

    auto DH = GetDepressionHierarchy<double,Topology::D8>(dem, label, flowdirs);

    //Initially distribute the water
    wtd.setAll(1);
    FillSpillMerge(dem, label, flowdirs, DH, wtd);

    const auto big_dump_wtd = wtd;

    wtd.setAll(0);
    for(int i=0;i<10;i++){
      for(auto i=wtd.i0();i<wtd.size();i++){
        wtd(i) += 0.1;
      }
      FillSpillMerge(dem, label, flowdirs, DH, wtd);
    }


    REQUIRE_MESSAGE(
      MaxArrayDiff(big_dump_wtd,wtd)<1e-6,
      "Randomized Testing of Repeated FSM failed with width = ", dem.width(), " height = ", dem.height(), " state = ", oss.str()
    );
  }
}

TEST_CASE("Randomized Testing of Incremental FSM vs Big Dump"){
  RandomizedIncrementalVsBigDump(number_of_small_tests,  10,  30);
  RandomizedIncrementalVsBigDump(number_of_large_tests, 100, 300);
}



void RandomizedMassConservation(const int count, const int min_size, const int max_size){
  #pragma omp parallel for
  for(int i=0;i<count;i++){
    std::stringstream oss;
    std::uniform_real_distribution<double> surface_amount_dist(0,1);

    Array2D<double> dem;
    double surface_water_amount;
    #pragma omp critical
    {
      oss<<gen;
      dem = random_integer_terrain(gen, min_size, max_size);
      surface_water_amount = surface_amount_dist(gen);
      std::cerr<<"RandomizedMassConservation #"<<i<<std::endl;
    }

    Array2D<dh_label_t> label   (dem.width(), dem.height(), NO_DEP);
    Array2D<flowdir_t>  flowdirs(dem.width(), dem.height(), NO_FLOW);
    Array2D<double>     wtd     (dem.width(), dem.height(), surface_water_amount);

    //Make sure the edges are identifiable as an ocean
    dem.setEdges(-1);
    label.setEdges(OCEAN);

    auto DH = GetDepressionHierarchy<double,Topology::D8>(dem, label, flowdirs);

    FillSpillMerge(dem, label, flowdirs, DH, wtd);

    double sum = 0;
    for(auto i=wtd.i0();i<wtd.size();i++){
      sum += wtd(i);
    }

    sum += DH.at(OCEAN).water_vol;

    REQUIRE_MESSAGE(
      sum==doctest::Approx(surface_water_amount*wtd.size()),
      "Randomized Testing of Repeated FSM failed with width = ", dem.width(), " height = ", dem.height(), " state = ", oss.str()
    );
  }
}

TEST_CASE("RandomizedMassConservation"){
  RandomizedMassConservation(number_of_small_tests,  10,  30);
  RandomizedMassConservation(number_of_large_tests, 100, 300);
}