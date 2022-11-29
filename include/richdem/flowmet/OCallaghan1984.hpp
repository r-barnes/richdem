#ifndef _richdem_OCallaghan1984_hpp_
#define _richdem_OCallaghan1984_hpp_

#include <richdem/common/Array2D.hpp>
#include <richdem/common/Array3D.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <richdem/common/constants.hpp>
#include <richdem/common/logger.hpp>

namespace richdem
{

  // TODO: Add Marks et al (1984)
  template <Topology topo, class elev_t>
  void FM_OCallaghan(
      const Array2D<elev_t> &elevations,
      Array3D<float> &props)
  {
    RDLOG_ALG_NAME << "O'Callaghan (1984)/Marks (1984) D8/D4 Flow Accumulation";
    RDLOG_CITATION << "O'Callaghan, J.F., Mark, D.M., 1984. The Extraction of Drainage Networks from Digital Elevation Data. Computer vision, graphics, and image processing 28, 323--344.";
    RDLOG_CONFIG << "topology = " << TopologyName(topo);

    constexpr auto dx = get_dx_for_topology<topo>();
    constexpr auto dy = get_dy_for_topology<topo>();
    constexpr auto nmax = get_nmax_for_topology<topo>();

    props.setAll(NO_FLOW_GEN);
    props.setNoData(NO_DATA_GEN);

    ProgressBar progress;
    progress.start(elevations.size());

#pragma omp parallel for collapse(2)
    for (int y = 0; y < elevations.height(); y++)
      for (int x = 0; x < elevations.width(); x++)
      {
        ++progress;

        if (elevations.isNoData(x, y))
        {
          props(x, y, 0) = NO_DATA_GEN;
          continue;
        }

        if (elevations.isEdgeCell(x, y))
          continue;

        const elev_t e = elevations(x, y);

        int lowest_n = 0;
        elev_t lowest_n_elev = std::numeric_limits<elev_t>::max();
        for (int n = 1; n <= nmax; n++)
        {
          const int nx = x + dx[n];
          const int ny = y + dy[n];

          if (elevations.isNoData(nx, ny)) // TODO: Don't I want water to drain this way?
            continue;

          const elev_t ne = elevations(nx, ny);

          if (ne >= e)
            continue;

          if (ne < lowest_n_elev)
          {
            lowest_n_elev = ne;
            lowest_n = n;
          }
        }

        if (lowest_n == 0)
          continue;

        props(x, y, 0) = HAS_FLOW_GEN;

        assert(e >= elevations(e + elevations.nshift(lowest_n))); // Ensure flow goes downhill

        props(x, y, lowest_n) = 1;
      }
    progress.stop();
  }

  template <class elev_t>
  void FM_D8(const Array2D<elev_t> &elevations, Array3D<float> &props)
  {
    FM_OCallaghan<Topology::D8>(elevations, props);
  }

  template <class elev_t>
  void FM_D4(const Array2D<elev_t> &elevations, Array3D<float> &props)
  {
    FM_OCallaghan<Topology::D4>(elevations, props);
  }

}

#endif
