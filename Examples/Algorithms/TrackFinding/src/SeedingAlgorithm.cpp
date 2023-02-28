// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/FpeMonitor.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <csignal>
#include <limits>
#include <stdexcept>
#include <signal.h>
#include <unistd.h>
#include <iomanip>


#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"

ActsExamples::SeedingAlgorithm::SeedingAlgorithm(
    ActsExamples::SeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  m_cfg.seedFinderConfig =
      m_cfg.seedFinderConfig.toInternalUnits().calculateDerivedQuantities();
  m_cfg.seedFinderOptions =
      m_cfg.seedFinderOptions.toInternalUnits().calculateDerivedQuantities(
          m_cfg.seedFinderConfig);
  m_cfg.seedFilterConfig = m_cfg.seedFilterConfig.toInternalUnits();
  m_cfg.gridConfig = m_cfg.gridConfig.toInternalUnits();
  m_cfg.gridOptions = m_cfg.gridOptions.toInternalUnits();
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collections");
  }
  for (const auto& i : m_cfg.inputSpacePoints) {
    if (i.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks output collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collection");
  }

  if (m_cfg.gridConfig.rMax != m_cfg.seedFinderConfig.rMax and
      m_cfg.allowSeparateRMax == false) {
    throw std::invalid_argument(
        "Inconsistent config rMax: using different values in gridConfig and "
        "seedFinderConfig. If values are intentional set allowSeparateRMax to "
        "true");
  }

  if (m_cfg.seedFilterConfig.deltaRMin != m_cfg.seedFinderConfig.deltaRMin) {
    throw std::invalid_argument("Inconsistent config deltaRMin");
  }

  if (m_cfg.gridConfig.deltaRMax != m_cfg.seedFinderConfig.deltaRMax) {
    throw std::invalid_argument("Inconsistent config deltaRMax");
  }

  static_assert(std::numeric_limits<decltype(
                    m_cfg.seedFinderConfig.deltaRMaxTopSP)>::has_quiet_NaN,
                "Value of deltaRMaxTopSP must support NaN values");

  static_assert(std::numeric_limits<decltype(
                    m_cfg.seedFinderConfig.deltaRMinTopSP)>::has_quiet_NaN,
                "Value of deltaRMinTopSP must support NaN values");

  static_assert(std::numeric_limits<decltype(
                    m_cfg.seedFinderConfig.deltaRMaxBottomSP)>::has_quiet_NaN,
                "Value of deltaRMaxBottomSP must support NaN values");

  static_assert(std::numeric_limits<decltype(
                    m_cfg.seedFinderConfig.deltaRMinBottomSP)>::has_quiet_NaN,
                "Value of deltaRMinBottomSP must support NaN values");

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMaxTopSP)) {
    m_cfg.seedFinderConfig.deltaRMaxTopSP = m_cfg.seedFinderConfig.deltaRMax;
  }

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMinTopSP)) {
    m_cfg.seedFinderConfig.deltaRMinTopSP = m_cfg.seedFinderConfig.deltaRMin;
  }

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMaxBottomSP)) {
    m_cfg.seedFinderConfig.deltaRMaxBottomSP = m_cfg.seedFinderConfig.deltaRMax;
  }

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMinBottomSP)) {
    m_cfg.seedFinderConfig.deltaRMinBottomSP = m_cfg.seedFinderConfig.deltaRMin;
  }

  if (m_cfg.gridConfig.zMin != m_cfg.seedFinderConfig.zMin) {
    throw std::invalid_argument("Inconsistent config zMin");
  }

  if (m_cfg.gridConfig.zMax != m_cfg.seedFinderConfig.zMax) {
    throw std::invalid_argument("Inconsistent config zMax");
  }

  if (m_cfg.seedFilterConfig.maxSeedsPerSpM !=
      m_cfg.seedFinderConfig.maxSeedsPerSpM) {
    throw std::invalid_argument("Inconsistent config maxSeedsPerSpM");
  }

  if (m_cfg.gridConfig.cotThetaMax != m_cfg.seedFinderConfig.cotThetaMax) {
    throw std::invalid_argument("Inconsistent config cotThetaMax");
  }

  if (m_cfg.gridConfig.minPt != m_cfg.seedFinderConfig.minPt) {
    throw std::invalid_argument("Inconsistent config minPt");
  }

  if (m_cfg.gridOptions.bFieldInZ != m_cfg.seedFinderOptions.bFieldInZ) {
    throw std::invalid_argument("Inconsistent config bFieldInZ");
  }

  if (m_cfg.gridConfig.zBinEdges.size() - 1 != m_cfg.zBinNeighborsTop.size() &&
      m_cfg.zBinNeighborsTop.empty() == false) {
    throw std::invalid_argument("Inconsistent config zBinNeighborsTop");
  }

  if (m_cfg.gridConfig.zBinEdges.size() - 1 !=
          m_cfg.zBinNeighborsBottom.size() &&
      m_cfg.zBinNeighborsBottom.empty() == false) {
    throw std::invalid_argument("Inconsistent config zBinNeighborsBottom");
  }

  if (!m_cfg.seedFinderConfig.zBinsCustomLooping.empty()) {
    // check if zBinsCustomLooping contains numbers from 1 to the total number
    // of bin in zBinEdges
    for (size_t i = 1; i != m_cfg.gridConfig.zBinEdges.size(); i++) {
      if (std::find(m_cfg.seedFinderConfig.zBinsCustomLooping.begin(),
                    m_cfg.seedFinderConfig.zBinsCustomLooping.end(),
                    i) == m_cfg.seedFinderConfig.zBinsCustomLooping.end()) {
        throw std::invalid_argument(
            "Inconsistent config zBinsCustomLooping does not contain the same "
            "bins as zBinEdges");
      }
    }
  }

  if (m_cfg.seedFinderConfig.useDetailedDoubleMeasurementInfo) {
    m_cfg.seedFinderConfig.getTopHalfStripLength.connect(
        [](const void*, const SimSpacePoint& sp) -> float {
          return sp.topHalfStripLength();
        });

    m_cfg.seedFinderConfig.getBottomHalfStripLength.connect(
        [](const void*, const SimSpacePoint& sp) -> float {
          return sp.bottomHalfStripLength();
        });

    m_cfg.seedFinderConfig.getTopStripDirection.connect(
        [](const void*, const SimSpacePoint& sp) -> Acts::Vector3 {
          return sp.topStripDirection();
        });

    m_cfg.seedFinderConfig.getBottomStripDirection.connect(
        [](const void*, const SimSpacePoint& sp) -> Acts::Vector3 {
          return sp.bottomStripDirection();
        });

    m_cfg.seedFinderConfig.getStripCenterDistance.connect(
        [](const void*, const SimSpacePoint& sp) -> Acts::Vector3 {
          return sp.stripCenterDistance();
        });

    m_cfg.seedFinderConfig.getTopStripCenterPosition.connect(
        [](const void*, const SimSpacePoint& sp) -> Acts::Vector3 {
          return sp.topStripCenterPosition();
        });
  }

  m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SimSpacePoint>>(m_cfg.seedFilterConfig);
  m_seedFinder = Acts::SeedFinder<SimSpacePoint>(m_cfg.seedFinderConfig);
}


void dumpNVec(const Acts::NeighborhoodVector &a_vec, Acts::SpacePointGrid<ActsExamples::SimSpacePoint> &grid) {
   for (const auto &val : a_vec) {
      std::cout << " " << val;
      try {
         std::cout << "[#" << grid.at(val).size() << "]" << std::endl;
      }
      catch (...) {
         std::cout << "[E]" << std::endl;
      }
   }
   std::cout << std::endl;
}
void dumpGrid(const Acts::SpacePointGrid<ActsExamples::SimSpacePoint> &grid) {
   std::size_t nbins = grid.size();
   for (unsigned int bin_i=0u; bin_i<nbins; ++bin_i) {
      try {
         if (!grid.at(bin_i).empty()) {
            auto idx = grid.localBinsFromGlobalBin(bin_i);
            auto ll = grid.lowerLeftBinEdge(idx);
            auto ur = grid.upperRightBinEdge(idx);
            std::cout << std::setw(5) << bin_i;
            for (unsigned int axis_i=0u; axis_i < Acts::SpacePointGrid<ActsExamples::SimSpacePoint>::DIM; ++axis_i) {
               std::cout << " " << std::setw(3) << idx[axis_i]
                         << ":" << std::setw(14) << ll.at(axis_i) << ".." << std::setw(14) << ur.at(axis_i);
            }
            std::cout << " #" << std::setw(4) << grid.at(bin_i).size();
            for (const auto &elm : grid.at(bin_i)) {
               std::cout << " {"
                  // << std::setw(14) << sqrt( sqr( elm->x() ) + sqr(elm->y()) ) << ", "
                         << std::setw(14) << elm->radius() << ", " << std::setw(14)<< elm->phi() << ", "
                         << std::setw(14) << elm->z() << "}";
            }
            std::cout << std::endl;
         }
      }
      catch (...) {
         std::cout << "[" << bin_i << " E]" << std::endl;
      }
   }
   std::cout << std::endl;
}


void dumpSPs(const Acts::SpacePointGrid<ActsExamples::SimSpacePoint> &grid,
             const Acts::GeometryContext &geo_ctx,
             const Acts::TrackingGeometry &trackingGeometry,
             float rmin, float rmax,
             const ActsExamples::AlgorithmContext *ctx=nullptr, const std::string &prefix="") {
   std::size_t nbins = grid.size();

   std::vector<std::pair<unsigned int,unsigned int> > hits;
   for (unsigned int bin_i=0u; bin_i<nbins; ++bin_i) {
      try {
         if (!grid.at(bin_i).empty()) {
            for (unsigned int hit_i=0; hit_i< grid.at(bin_i).size(); ++hit_i) {
               hits.push_back( std::make_pair(bin_i, hit_i));
            }
         }
      }
      catch (...) {
         std::cerr << "[" << bin_i << " E]" << std::endl;
      }
   }
   std::sort( hits.begin(), hits.end(), [&grid](const std::pair<unsigned int, unsigned int> &a,
                                                const std::pair<unsigned int, unsigned int> &b) {
                 const double r_a = grid.at(a.first).at(a.second)->radius();
                 const double r_b = grid.at(b.first).at(b.second)->radius();
                 return r_a < r_b ;
              });

   const ActsExamples::MeasurementContainer *measurements=nullptr;
   if (ctx) {
      measurements= &ctx->eventStore.get<ActsExamples::MeasurementContainer>("measurements");
   }

   // dummy truth
   std::cout << prefix
             << 1  // PDG ID
             << " " << 1  // barcode
             << " " << 1. // pt
             << " " << 0. << " " << 0. << " " << 0.; // direction

   Acts::Vector3 dummy_momentum{};

   unsigned int n_hits=0;
   for (const std::pair<unsigned int, unsigned int> &hit : hits ) {
      unsigned int bin_i=hit.first;
      unsigned int elm_i=hit.second;
      if (!grid.at(bin_i).empty()) {
         try {

         // measurements
         const auto &sp = grid.at(bin_i).at(elm_i)->sp();
         ++n_hits;

         for(const auto &sl : sp.sourceLinks() ) {
            const auto geoId = sl.geometryId();
            (void) sl.get<ActsExamples::IndexSourceLink>();
            const Acts::Surface* a_surface=trackingGeometry.findSurface(geoId);
            if (a_surface) {
               ++n_hits;
            }
         }
         }
         catch (...) {}
      }
   }

   //   std::cout << " " << hits.size();
   std::cout << " " << n_hits;
   for (const std::pair<unsigned int, unsigned int> &hit : hits ) {
      unsigned int bin_i=hit.first;
      unsigned int elm_i=hit.second;
      if (!grid.at(bin_i).empty()) {
         try {
         // space point
         std::cout << " " << -1; // number of contributions
         std::cout << " " << grid.at(bin_i).at(elm_i)->x()
                   << " " << grid.at(bin_i).at(elm_i)->y()
                   << " " << grid.at(bin_i).at(elm_i)->z();
         std::cout << " " << -2; // number of surface courners. -2 -> box edges

         auto idx = grid.localBinsFromGlobalBin(bin_i);
         auto ll = grid.lowerLeftBinEdge(idx);
         auto ur = grid.upperRightBinEdge(idx);

         std::cout << " " << rmin*cos(ll[0]) << " " << rmin*sin(ll[0]) << " " << ll[1];
         std::cout << " " << rmax*cos(ur[0]) << " " << rmax*sin(ur[0]) << " " << ur[1];

         // measurements
         const auto &sp = grid.at(bin_i).at(elm_i)->sp();

         for(const auto &sl : sp.sourceLinks() ) {
            const auto geoId = sl.geometryId();
            const ActsExamples::IndexSourceLink& islink = sl.get<ActsExamples::IndexSourceLink>();
            const Acts::Surface* a_surface=trackingGeometry.findSurface(geoId);
            if (a_surface) {
               Acts::Vector3 glob;
               if (measurements) {
                  std::visit([&geo_ctx, &dummy_momentum, a_surface, &glob](const auto& a_measurement) {
                        auto expander = a_measurement.expander();
                        Acts::BoundVector par = expander * a_measurement.parameters();
                        // extract local position
                        Acts::Vector2 lpar(par[Acts::eBoundLoc0], par[Acts::eBoundLoc1]);
                        glob = a_surface->localToGlobal(geo_ctx,
                                                        lpar,
                                                        dummy_momentum);
                     },measurements->at(islink.index()) );
               }
               else {
                  glob=Acts::Vector3(sp.x(),sp.y(),sp.z());
               }
               std::cout << " " << -1; // number of contributions
               std::cout << " " << glob[0]
                         << " " << glob[1]
                         << " " << glob[2];
               const Acts::AnnulusBounds *annulus_bounds = dynamic_cast<const Acts::AnnulusBounds *>( &a_surface->bounds() );
               std::vector<Acts::Vector2> corners;
               if (annulus_bounds) {
                  corners = annulus_bounds->corners();
               }
               else {
                  const Acts::PlanarBounds *bounds = dynamic_cast<const Acts::PlanarBounds *>( &a_surface->bounds() );
                  if (bounds) {
                     corners = bounds->vertices(1);
                  }
               }
               std::cout << " " << corners.size();
               for (const Acts::Vector2 &a_corner  : corners) {
                  assert( a_surface != nullptr);
                  Acts::Vector3 glob_corner = a_surface->localToGlobal(geo_ctx,
                                                                a_corner,
                                                                dummy_momentum);
                  std::cout << " " << glob_corner[0] << " " << glob_corner[1] << " " << glob_corner[2];
               }
            }
         }
         }
         catch (...) {
            std::cerr << "[" << bin_i << " E]" << std::endl;
         }
      }
   }
   std::cout << std::endl;
}

void dumpSPs(const Acts::SpacePointGrid<ActsExamples::SimSpacePoint> &grid,
             const Acts::TrackingGeometry &trackingGeometry,
             float rmin, float rmax,
             const ActsExamples::AlgorithmContext *ctx=nullptr) {
   Acts::GeometryContext geo_ctx;
   dumpSPs(grid,geo_ctx, trackingGeometry, rmin, rmax,ctx);
}

namespace {
   bool first_occ=true;
}
ActsExamples::ProcessCode ActsExamples::SeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // construct the combined input container of space point pointers from all
  // configured input sources.
  // pre-compute the total size required so we only need to allocate once
   if (m_cfg.stop && first_occ) {
   std::cout << "DEBUG " << getpid() << " wait for signal CONT " << std::endl;
   kill(getpid(), SIGSTOP);
   first_occ=false;
   }

  const Acts::TrackingGeometry *tg=m_cfg.trackingGeometry.get();
  const Acts::GeometryContext& geo_ctx= ctx.geoContext;

  size_t nSpacePoints = 0;
  for (const auto& isp : m_cfg.inputSpacePoints) {
    nSpacePoints += ctx.eventStore.get<SimSpacePointContainer>(isp).size();
  }

  std::vector<const SimSpacePoint*> spacePointPtrs;
  spacePointPtrs.reserve(nSpacePoints);
  for (const auto& isp : m_cfg.inputSpacePoints) {
    for (const auto& spacePoint :
         ctx.eventStore.get<SimSpacePointContainer>(isp)) {
      // since the event store owns the space points, their pointers should be
      // stable and we do not need to create local copies.
      spacePointPtrs.push_back(&spacePoint);
    }
  }

  // construct the seeding tools
  // covariance tool, extracts covariances per spacepoint as required
  auto extractGlobalQuantities =
      [=](const SimSpacePoint& sp, float, float,
          float) -> std::pair<Acts::Vector3, Acts::Vector2> {
    Acts::Vector3 position{sp.x(), sp.y(), sp.z()};
    Acts::Vector2 covariance{sp.varianceR(), sp.varianceZ()};
    return std::make_pair(position, covariance);
  };

  // extent used to store r range for middle spacepoint
  Acts::Extent rRangeSPExtent;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
      Acts::BinFinder<SimSpacePoint>(m_cfg.zBinNeighborsBottom,
                                     m_cfg.numPhiNeighbors));
  auto topBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
      Acts::BinFinder<SimSpacePoint>(m_cfg.zBinNeighborsTop,
                                     m_cfg.numPhiNeighbors));
  auto grid = Acts::SpacePointGridCreator::createGrid<SimSpacePoint>(
      m_cfg.gridConfig, m_cfg.gridOptions);
  auto spacePointsGrouping = Acts::BinnedSPGroup<SimSpacePoint>(
      spacePointPtrs.begin(), spacePointPtrs.end(), extractGlobalQuantities,
      bottomBinFinder, topBinFinder, std::move(grid), rRangeSPExtent,
      m_cfg.seedFinderConfig, m_cfg.seedFinderOptions);

  // safely clamp double to float
  float up = Acts::clampValue<float>(
      std::floor(rRangeSPExtent.max(Acts::binR) / 2) * 2);

  /// variable middle SP radial region of interest
  const Acts::Range1D<float> rMiddleSPRange(
      std::floor(rRangeSPExtent.min(Acts::binR) / 2) * 2 +
          m_cfg.seedFinderConfig.deltaRMiddleMinSPRange,
      up);

  dumpSPs(spacePointsGrouping.grid(),
          geo_ctx,
          *tg,
          rRangeSPExtent.min(Acts::binR), rRangeSPExtent.max(Acts::binR),&ctx, "PARTSEED");


  // run the seeding
  static thread_local SimSeedContainer seeds;
  seeds.clear();
  static thread_local decltype(m_seedFinder)::SeedingState state;

  auto group = spacePointsGrouping.begin();
  auto groupEnd = spacePointsGrouping.end();
  for (; !(group == groupEnd); ++group) {
    m_seedFinder.createSeedsForGroup(
        m_cfg.seedFinderOptions, state, std::back_inserter(seeds),
        group.bottom(), group.middle(), group.top(), rMiddleSPRange);
  }

  // extract proto tracks, i.e. groups of measurement indices, from tracks seeds
  size_t nSeeds = seeds.size();
  static thread_local ProtoTrackContainer protoTracks;
  protoTracks.clear();

  protoTracks.reserve(nSeeds);
  for (const auto& seed : seeds) {
    ProtoTrack& protoTrack = protoTracks.emplace_back();
    protoTrack.reserve(seed.sp().size());
    for (auto spacePointPtr : seed.sp()) {
      for (const auto& slink : spacePointPtr->sourceLinks()) {
        const IndexSourceLink& islink = slink.get<IndexSourceLink>();
        protoTrack.emplace_back(islink.index());
      }
    }
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePointPtrs.size() << " space points");

  ctx.eventStore.add(m_cfg.outputSeeds, SimSeedContainer{seeds});
  ctx.eventStore.add(m_cfg.outputProtoTracks, ProtoTrackContainer{protoTracks});
  return ActsExamples::ProcessCode::SUCCESS;
}
