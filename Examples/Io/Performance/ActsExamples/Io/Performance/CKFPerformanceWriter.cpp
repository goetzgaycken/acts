// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "Acts/Definitions/Units.hpp"

#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/EventDataView3D.hpp"

#include <numeric>
#include <stdexcept>

#include <unistd.h>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TVectorF.h>

namespace ActsExamples {
   std::ostream &operator<<(std::ostream &out, const Stat &a ) {
      out <<           std::setw(14) << a.min() << " < "
          <<           std::setw(14) << a.mean()
          << " +- " << std::setw(14) << a.rms()
          << " < "  << std::setw(14) << a.max()
          << " / "  << std::setw(9)  << a.n();
      if (!a.bins().empty()) {
         out << std::endl;
         out << a.lowerEdge() << ", " << a.lowerEdge()+a.binWidth() << ".." << a.upperEdge() << ": ";
         out << std::setw(5) << a.bins()[0] << " |";
         for (unsigned int bin_i=1; bin_i<a.bins().size()-1; ++bin_i) {
            out << " " << std::setw(5) << a.bins()[bin_i];
         }
         out << " | " << std::setw(5) << a.bins().back() ;
      }
      return out;
   }
}


ActsExamples::CKFPerformanceWriter::CKFPerformanceWriter(
    ActsExamples::CKFPerformanceWriter::Config cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputTrajectories, "CKFPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_effPlotTool(m_cfg.effPlotToolConfig, lvl),
      m_fakeRatePlotTool(m_cfg.fakeRatePlotToolConfig, lvl),
      m_duplicationPlotTool(m_cfg.duplicationPlotToolConfig, lvl),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig, lvl) {
  // trajectories collection name is already checked by base ctor
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::invalid_argument("Could not open '" + m_cfg.filePath + "'");
  }

  // initialize the plot tools
  m_effPlotTool.book(m_effPlotCache);
  m_fakeRatePlotTool.book(m_fakeRatePlotCache);
  m_duplicationPlotTool.book(m_duplicationPlotCache);
  m_trackSummaryPlotTool.book(m_trackSummaryPlotCache);

  std::cout << "DEBUG ActsExamples::CKFPerformanceWriter pid " << getpid() << " sleep..."<< std::endl;
  sleep(15);
  std::cout << "DEBUG ActsExamples::CKFPerformanceWriter ... continue." << std::endl;
}

ActsExamples::CKFPerformanceWriter::~CKFPerformanceWriter() {
  m_effPlotTool.clear(m_effPlotCache);
  m_fakeRatePlotTool.clear(m_fakeRatePlotCache);
  m_duplicationPlotTool.clear(m_duplicationPlotCache);
  m_trackSummaryPlotTool.clear(m_trackSummaryPlotCache);
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::CKFPerformanceWriter::endRun() {
  float eff_tracks = float(m_nTotalMatchedTracks) / m_nTotalTracks;
  float fakeRate_tracks = float(m_nTotalFakeTracks) / m_nTotalTracks;
  float duplicationRate_tracks =
      float(m_nTotalDuplicateTracks) / m_nTotalTracks;

  float eff_particle = float(m_nTotalMatchedParticles) / m_nTotalParticles;
  float fakeRate_particle = float(m_nTotalFakeParticles) / m_nTotalParticles;
  float duplicationRate_particle =
      float(m_nTotalDuplicateParticles) / m_nTotalParticles;

  ACTS_DEBUG("nTotalTracks                = " << m_nTotalTracks);
  ACTS_DEBUG("nTotalMatchedTracks         = " << m_nTotalMatchedTracks);
  ACTS_DEBUG("nTotalDuplicateTracks       = " << m_nTotalDuplicateTracks);
  ACTS_DEBUG("nTotalFakeTracks            = " << m_nTotalFakeTracks);

  ACTS_INFO(
      "Efficiency with tracks (nMatchedTracks/ nAllTracks) = " << eff_tracks);
  ACTS_INFO(
      "Fake rate with tracks (nFakeTracks/nAllTracks) = " << fakeRate_tracks);
  ACTS_INFO("Duplicate rate with tracks (nDuplicateTracks/nAllTracks) = "
            << duplicationRate_tracks);
  ACTS_INFO("Efficiency with particles (nMatchedParticles/nTrueParticles) = "
            << eff_particle);
  ACTS_INFO("Fake rate with particles (nFakeParticles/nTrueParticles) = "
            << fakeRate_particle);
  ACTS_INFO(
      "Duplicate rate with particles (nDuplicateParticles/nTrueParticles) = "
      << duplicationRate_particle);

  dumpStat();

  auto write_float = [&](float f, const char* name) {
    TVectorF v(1);
    v[0] = f;
    m_outputFile->WriteObject(&v, name);
  };

  if (m_outputFile != nullptr) {
    m_outputFile->cd();
    m_effPlotTool.write(m_effPlotCache);
    m_fakeRatePlotTool.write(m_fakeRatePlotCache);
    m_duplicationPlotTool.write(m_duplicationPlotCache);
    m_trackSummaryPlotTool.write(m_trackSummaryPlotCache);
    write_float(eff_tracks, "eff_tracks");
    write_float(fakeRate_tracks, "fakerate_tracks");
    write_float(duplicationRate_tracks, "duplicaterate_tracks");
    write_float(eff_particle, "eff_particles");
    write_float(fakeRate_particle, "fakerate_particles");
    write_float(duplicationRate_particle, "duplicaterate_particles");
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

class Counter {
public:
   Counter &operator++() { ++m_counter; return *this;}
   unsigned int counts() const { return m_counter; }
   unsigned int m_counter=0;
};

namespace {
   using ActsExamples::sqr;
   enum EStat {
      kMeasPerTraj,
      kPartPerTraj,
      kMeasPt,
      kPartPt,
      kNSharedMeasPerTraj,
      kTrajPerPart,
      kTrajPerPartWith2Meas,
      kTrajPerPartWith3Meas,
      kTrajPerPartWith50Min3Meas,
      kTrajPerPartWith90Min3Meas,
      kTrajPerPartWith100Min3Meas,
      kNParticles,
      kNStat
   };
   constexpr std::array<const char *,kNStat> statNames() {
      return std::array<const char *,kNStat> {
            "Hits per track",
            "Particles per track",
            "pt at hit",
            "pt of particle",
            "shared hits per track",
            "tracks per particle",
            "tracks per particle >=2 part. hits",
            "tracks per particle >=3 part. hits",
            "tracks per particle >=3 and >=50% part hits",
            "tracks per particle >=3 and >=90% part hits",
            "tracks per particle >=3 and >=100% part hits",
            "particles per hit",
            };
   }
   template <class T1, class T2>
   void dumpStatPtBins(const T1 &stat_pt_bins, const T2 &stat) {
     static const std::array<const char *,kNStat> stat_names( statNames());
     std::size_t max_length=0;
     for (const char * a_name : stat_names) {
        max_length=std::max(max_length, strlen(a_name));
     }
     for (unsigned int pt_i=0; pt_i < stat_pt_bins.size(); ++pt_i) {
        std::cout << "-- Statistics pt bin " << stat_pt_bins[pt_i] << " ... "<< std::endl;
        for (unsigned int stat_i=0; stat_i<stat.at(pt_i).size(); ++stat_i ) {
           std::cout << stat_names.at(stat_i) << std::setw( 31 - std::min(static_cast<std::size_t>(30),strlen(stat_names.at(stat_i))) )
                     << " " << stat.at(pt_i).at(stat_i)
                     << std::endl;
        }
        std::cout << std::endl;
     }
   }


   template <class T>
   double transverseMomentum(const T &param) {
      return param[Acts::eBoundQOverP] != 0. ? std::sin(param[Acts::eBoundTheta]) / std::abs(param[Acts::eBoundQOverP]) : 0.;
   }

   Acts::Vector3 globalCoords(const Acts::TrackingGeometry &tracking_geometetry,
                              const Acts::GeometryContext& gctx,
                              const ActsExamples::Measurement& meas) {
      const auto& slink =
         std::visit([](const auto& x) { return &x.sourceLink(); }, meas);

      const auto geoId = slink->geometryId();

      const Acts::Surface* surface = tracking_geometetry.findSurface(geoId);
      Acts::Vector2 localPos = std::visit(
                                             [](const auto& measurement) {
                                                auto expander = measurement.expander();
                                                Acts::BoundVector par = expander * measurement.parameters();
                                                Acts::Vector2 lpar(par[Acts::eBoundLoc0], par[Acts::eBoundLoc1]);
                                                return lpar;
                                             },
                                             meas);
      Acts::Vector3 globalPos = surface->localToGlobal(gctx, localPos, Acts::Vector3());
      return globalPos;
   }
   double computeChi2ForHit(const Acts::Vector2 &hitPos,
                            const Acts::Vector2 &fitPos,
                            const Acts::Vector2 &fitPosCovDiag) {
      double chi2=0.;
      for (unsigned int idx=0; idx < hitPos.size(); ++idx) {
         if (fitPosCovDiag[idx]>1e-6) {
            chi2+=sqr(hitPos[idx]-fitPos[idx])/fitPosCovDiag[idx];
         }
      }
      return chi2;
   }
   std::pair<double, unsigned int> modifiedChi2(const std::tuple<double, unsigned int, std::array<float, 4> > &chi2) {
      static unsigned int debug_counter=100;
      std::pair<double, unsigned int>  ret = std::make_pair(std::get<0>(chi2),std::get<1>(chi2) ) ;
      if (ret.second>1 && ret.first>0.) {
         double old_chi2=ret.first;// / ret.second;
         for (unsigned int arr_i=0; arr_i< std::get<2>(chi2).size();  ++arr_i) {
            double mod_chi2 = (ret.first - std::get<2>(chi2).at(arr_i));
               // / (ret.second - 1);
            if (++debug_counter<100) {
               std::cerr << "DEBUG modifiedChi2 " << arr_i << " : " << old_chi2
                         << " -> " << ret.first
                         << " - " << std::get<2>(chi2).at(arr_i) << " / " << (ret.second - 1) 
                         << " = " << mod_chi2
                         << " ( " << (mod_chi2/old_chi2) << " )"
                         << std::endl;
            }
            if ((old_chi2 - mod_chi2) * ret.second < 3 * old_chi2 ) break;
            old_chi2 = mod_chi2;
            ret.first -= std::get<2>(chi2).at(arr_i);
            --ret.second;
            if (ret.second<=1 || ret.second <= 3 * ret.first) break;
         }
      }
      return ret;
   }
      
   void dumpParameters(const std::string &header,
                       const std::vector<std::array<double,7> > &trajectory_params,
                       const unsigned int floats_per_column);

   void dumpParameters(const std::string &header,
                       const ActsExamples::AlgorithmContext& ctx,
                       const std::map<unsigned int,Counter> &trajectory_counts,
                       const std::map<unsigned int,const Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy> &hit_states,
                       std::vector<double> &pt_out,
                       const unsigned int floats_per_column,
                       bool really_dump_params) {
      std::vector<std::array<double,7> > trajectory_params;
      trajectory_params.reserve( trajectory_counts.size() );
      bool has_param=false;
      pt_out.resize(trajectory_counts.size(),-std::numeric_limits<double>::max());
      unsigned idx=0;
      for (const std::pair<const unsigned int,Counter> &trajectory : trajectory_counts) {
         std::map<unsigned int,const Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy>::const_iterator
            traj_iter = hit_states.find(trajectory.first);
         if (traj_iter != hit_states.end()) {
            const auto &boundParams = traj_iter->second.parameters();
            Acts::Vector3 direction = Acts::makeDirectionUnitFromPhiTheta(boundParams[Acts::eBoundPhi],
                                                                          boundParams[Acts::eBoundTheta]);

            Acts::Vector2 local(boundParams[Acts::eBoundLoc0], boundParams[Acts::eBoundLoc1]);
            Acts::Vector3 position(traj_iter->second.referenceSurface().localToGlobal(ctx.geoContext, local, direction));
            trajectory_params.emplace_back(std::array<double,7>{ Acts::VectorHelpers::phi(position),
                                                                Acts::VectorHelpers::theta(position),
                                                                Acts::VectorHelpers::perp(position),
                                                                position[Acts::ePos2],
                                                                boundParams[Acts::eBoundPhi],
                                                                boundParams[Acts::eBoundTheta],
                                                                std::copysign ( transverseMomentum(boundParams)/Acts::UnitConstants::GeV,
                                                                                boundParams[Acts::eBoundQOverP])
               });
            pt_out.at(idx)= std::abs(trajectory_params.back()[6]);
            has_param=true;
         }
         else {
            trajectory_params.emplace_back(std::array<double,7>{});
         }
         ++idx;
      }
      if (has_param && really_dump_params) {
         dumpParameters(header, trajectory_params, floats_per_column);
      }
   }
   void dumpPerigeeParameters(const std::string &header,
                              const ActsExamples::AlgorithmContext& ctx,
                              const std::map<unsigned int,Counter> &trajectory_counts,
                              const std::map<unsigned int, std::pair<size_t, size_t> > &trajectoryIdMap,
                              const ActsExamples::TrajectoriesContainer& trajectories,
                              std::vector<double> &pt_out,
                              const unsigned int floats_per_column) {
      std::vector<std::array<double,7> > trajectory_params;
      trajectory_params.reserve( trajectory_counts.size() );
      bool has_param=false;
      pt_out.resize(trajectory_counts.size(),-std::numeric_limits<double>::max());
      unsigned idx=0;
      for (const std::pair<const unsigned int,Counter> &trajectory : trajectory_counts) {
         std::map<unsigned int, std::pair<size_t, size_t> >::const_iterator traj_iter = trajectoryIdMap.find(trajectory.first);
         if (traj_iter != trajectoryIdMap.end()) {
            const auto& traj = trajectories[traj_iter->second.first];
            const auto& fittedParameters = traj.trackParameters(traj_iter->second.second);
            const auto& momentum = fittedParameters.momentum();
            Acts::Vector3 position(fittedParameters.position(ctx.geoContext));
            trajectory_params.emplace_back(std::array<double,7>{ Acts::VectorHelpers::phi(position),
                                                                Acts::VectorHelpers::theta(position),
                                                                Acts::VectorHelpers::perp(position),
                                                                position[Acts::ePos2],
                                                                Acts::VectorHelpers::phi(momentum),
                                                                Acts::VectorHelpers::theta(momentum),
                                                                std::copysign ( Acts::VectorHelpers::perp(momentum)/Acts::UnitConstants::GeV,
                                                                                fittedParameters.charge())
               });
            pt_out.at(idx)= std::abs(trajectory_params.back()[6]);
            has_param=true;
         }
         else {
            trajectory_params.emplace_back(std::array<double,7>{});
         }
         ++idx;
      }
      if (has_param) {
         dumpParameters(header, trajectory_params, floats_per_column);
      }
   }

   void dumpParameters(const std::string &header,
                       const std::vector<std::array<double,7> > &trajectory_params,
                       unsigned int floats_per_column) {
      static std::array<std::string,7> param_names{"phi","theta","r","z","phi","theta","(+-)pt"};
      for (unsigned int offset_i=0;offset_i<trajectory_params.size();) {
         unsigned int traj_end  = std::min(static_cast<unsigned int>(trajectory_params.size()), offset_i+floats_per_column);
         for (unsigned int param_i=0; param_i<7; ++param_i) {
            if (param_i==0 && offset_i==0) {
               std::cout << std::setw(4) << header << ") ";
            }
            else {
               std::cout << std::setw(4) << " " << "  ";
            }
            std::cout << std::setw(7) << param_names.at(param_i) << " ";
            for (unsigned int traj_i=offset_i; traj_i<traj_end;++traj_i) {
               const std::array<double, 7> &param = trajectory_params[traj_i];
               if (std::abs(param[5])>0) {
                  std::cout << std::setw(14) << param[param_i];
               }
               else {
                  std::cout << std::setw(14) << " ";
               }
            }
            std::cout << std::endl;
         }
         offset_i = traj_end;
         std::cout << std::endl;
      }
   }
}

void  ActsExamples::CKFPerformanceWriter::dumpStat() {
  dumpStatPtBins( m_ptBins, m_stat);
}

ActsExamples::ProcessCode ActsExamples::CKFPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoriesContainer& trajectories) {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
  // The number of majority particle hits and fitted track parameters
  using RecoTrackInfo = std::pair<size_t, Acts::BoundTrackParameters>;
  using Acts::VectorHelpers::perp;

  // Read truth input collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);
  const auto& seedIdx =
     ctx.eventStore.get<std::vector<std::size_t> >("TrackSeedIdx");
  const auto& protoTracks =
      ctx.eventStore.get<ProtoTrackContainer>("extended_proto_tracks");
  const auto& measurements =
     ctx.eventStore.get<ActsExamples::MeasurementContainer>("measurements");


  // Counter of truth-matched reco tracks
  std::map<ActsFatras::Barcode, std::vector<RecoTrackInfo>> matched;
  // Counter of truth-unmatched reco tracks
  std::map<ActsFatras::Barcode, size_t> unmatched;
  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Vector of input features for neural network classification
  std::vector<float> inputFeatures(3);

  if (m_cfg.dumpDuplicates) {
     std::cout << "particles " << hitParticlesMap.size() << std::endl;
  }
  std::map<std::pair<size_t, size_t>, unsigned int >              trajectoryId;
  std::map<unsigned int, std::pair<size_t, size_t> >              trajectoryIdMap;
  std::map<ActsFatras::Barcode, unsigned int >                    shortParticleId;
  std::map<ActsFatras::Barcode, size_t >                          barcodeToParticle;
  std::map<ActsFatras::Barcode, std::map<unsigned int, Counter> > trajectoriesPerParticle;
  std::map<Index, std::set< unsigned int > >                      hitTrajectoryMap;
  std::map<unsigned int, std::vector<Index> >                     trajectoryToHits;
  std::map<unsigned int, std::map<ActsFatras::Barcode, Counter> > trajectoryBarcodeMap;
  std::map<unsigned int, std::map<unsigned int,
                                  const Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy> > hitToTrackState;

  static constexpr std::array<double,5 > stat_pt_bins {
        -1. * Acts::UnitConstants::GeV,
        0.  * Acts::UnitConstants::GeV,
        1.  * Acts::UnitConstants::GeV,
        2.5 * Acts::UnitConstants::GeV,
        5   * Acts::UnitConstants::GeV};
  constexpr unsigned int floats_per_col =10;
  std::array<std::vector<Stat>,stat_pt_bins.size() > stat;
  for (std::vector<Stat> &a_stat : stat ) {
     std::vector<Stat> tmp{
         Stat(15,-0.5,30-0.5), // MeasPerTraj
         Stat(10,-0.5,10-0.5), // PartPerTraj
         Stat(10,0,10.), //  kMeasPt,
         Stat(10,0,10.), //  kPartPt,
         Stat(15,-0.5,30-0.5), //kNSharedMeasPerTraj
         Stat(15,0-0.5, 30-0.5), // kTrajPerPart,
         Stat(15,0-0.5, 30-0.5),// kTrajPerPartWith2Meas,
         Stat(15,0-0.5, 30-0.5),// kTrajPerPartWith3Meas,
         Stat(15,0-0.5, 30-0.5),// kTrajPerPartWith50Min3Meas,
         Stat(15,0-0.5, 30-0.5),// kTrajPerPartWith90in3Meas,
         Stat(15,0-0.5, 30-0.5),// kTrajPerPartWith100Min3Meas,
         Stat(10,0-0.5, 10-0.5)// kNParticles
         };
     a_stat = std::move(tmp);
     if (a_stat.size() != kNStat) { throw std::logic_error("Stat array has wrong dimension."); }
  }
  for (const SimParticle &a_particle : particles) {
     double particle_pt = a_particle.transverseMomentum();
     for (unsigned int pt_i=0; pt_i < stat_pt_bins.size(); ++pt_i) {
        if (particle_pt < stat_pt_bins.at(pt_i) ) break;
        stat.at(pt_i)[kPartPt].add(std::abs(particle_pt));
     }
  }

  for (unsigned int meas_i=0; meas_i<measurements.size(); ++meas_i) {
     std::array<unsigned int, stat_pt_bins.size()> hit_share_count{};
     for (auto hitParticle : makeRange(hitParticlesMap.equal_range(meas_i))) {
        auto particle_iter = particles.find(hitParticle.second);
        if (particle_iter != particles.end()) {
           double particle_pt = particle_iter->transverseMomentum();
           for (unsigned int pt_i=0; pt_i < stat_pt_bins.size(); ++pt_i) {
              if (particle_pt < stat_pt_bins.at(pt_i) ) break;
              ++hit_share_count[pt_i];
           }
        }
     }
     for (unsigned int pt_i=0; pt_i < hit_share_count.size(); ++pt_i) {
        stat[pt_i]   [kNParticles].add( hit_share_count[pt_i] );
     }
  }

  std::map<unsigned int, unsigned int> trajToSeed;
  std::vector<unsigned int> processed;
  std::vector<double>       hit_r;
  std::vector<double>       traj_assoc_pt;
  std::vector<std::map<ActsFatras::Barcode,std::tuple<double, unsigned int, std::array<float, 4> > > >      hit_chi2;
  std::map<ActsFatras::Barcode, unsigned int >                                      particle_n_hits;

  unsigned int max_hit_index=0;
  unsigned int plot_i=0;
  std::set<unsigned int> visualised;

  static unsigned int debug_counter=0;
  // Loop over all trajectories
  for (auto [itraj, trackTip] : trackTips) {
    const auto& traj = trajectories[itraj];
    const auto& mj = traj.multiTrajectory();
    auto trajState =
       Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
    unsigned int seed_i=(itraj<=seedIdx.size() ? seedIdx[itraj] : std::numeric_limits<unsigned int>::max() );
    if (itraj>=seedIdx.size()) {
       std::cout << "WARNING " << __FUNCTION__ << " index mismatch " << itraj << " !< " << seedIdx.size() << std::endl;
    }

    if (traj.hasTrajectory(trackTip)) {
       std::pair<std::map<std::pair<size_t, size_t>, unsigned int >::iterator, bool>
          ret = trajectoryId.insert( std::make_pair( std::make_pair(itraj,trackTip), trajectoryId.size() ));
       if (ret.second) {
          trajectoryIdMap.insert( std::make_pair( ret.first->second, ret.first->first) );
       }
       unsigned int traj_id=ret.first->second;
       if (traj_id>=hit_chi2.size()) { hit_chi2.resize( traj_id + 1); }
       trajToSeed.insert( std::make_pair(traj_id, seed_i));
       bool visualize=false;
       bool dump_all=false;
       int state_i=0;
       traj.multiTrajectory().visitBackwards(trackTip, [&](const auto& state) {
          // no truth info with non-measurement state
          if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
             return true;
          }

          // register all particles that generated this hit
          const auto& sl = static_cast<const IndexSourceLink&>(state.uncalibrated());
          auto hitIndex = sl.index();
          max_hit_index=std::max(max_hit_index,hitIndex);
          hitToTrackState[hitIndex].insert( std::make_pair(traj_id, state) );
          hitTrajectoryMap[hitIndex].insert(traj_id);
          trajectoryToHits[traj_id].push_back(hitIndex);
          for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hitIndex))) {
             /*std::pair<std::map<ActsFatras::Barcode, unsigned int >::const_iterator, bool >
               particle_insert =*/ shortParticleId.insert( std::make_pair(hitParticle.second, shortParticleId.size()));
             ++(trajectoryBarcodeMap[ traj_id ][hitParticle.second]);
             ++(trajectoriesPerParticle[hitParticle.second][ traj_id ]);

             {
                const auto &boundParams = state.parameters();
                const auto &cov = state.covariance();

                Acts::Vector2(boundParams[Acts::eBoundLoc0], boundParams[Acts::eBoundLoc1]);
                std::size_t n_measurements = 0;
                Acts::Vector2 localPos = std::visit(
                                                    [&n_measurements](const auto& measurement) {
                                                       auto expander = measurement.expander();
                                                       Acts::BoundVector par = expander * measurement.parameters();
                                                       Acts::Vector2 lpar(par[Acts::eBoundLoc0], par[Acts::eBoundLoc1]);
                                                       n_measurements=measurement.size();
                                                       return lpar;
                                                    },
                                                    measurements.at(hitIndex));
                if (cov(Acts::eBoundLoc1,Acts::eBoundLoc1)<1e-6 || dump_all) {
                   if (debug_counter < 1000) {
                   visualize=true;
                   Acts::Vector3 pos( globalCoords(*(m_cfg.trackingGeometry), ctx.geoContext,measurements.at(hitIndex)) );
                   const auto &meas_calibrated = state.effectiveCalibrated();
                   const auto &calib_cov = state.effectiveCalibratedCovariance();
                   if (debug_counter<2) {
                      dump_all=true;
                   }

                   if (m_cfg.dumpDuplicates) {
                   std::cout << "DEBUG " 
                             << (cov(Acts::eBoundLoc1,Acts::eBoundLoc1)<1e-6 ? "suspicious loc1 cov in " : "" )
                             << " trajectory " << itraj << ", " << trackTip
                             << " state " << --state_i << " measurement : (r,phi,z) = ( " << Acts::VectorHelpers::perp(pos)
                             << ", " << Acts::VectorHelpers::phi(pos)
                             << ", " << pos[2] << " )"
                             << " : has "
                             << ( state.hasPredicted() ?  "PREDICTED " : "")
                             << ( state.hasFiltered()  ?  "FILTERED "  : "")
                             << ( state.hasSmoothed()  ?  "SMOOTHED "  : "") << std::endl
                             << " measurement " <<  meas_calibrated.transpose() << std::endl
                             << calib_cov << std::endl
                             << " parameters: " << boundParams.transpose() << std::endl
                             << " cov predicted: " << std::endl << ( state.hasPredicted() ?  state.predictedCovariance() : cov) << std::endl
                             << " cov filtered: " << std::endl << ( state.hasFiltered()   ?  state.filteredCovariance() : cov) << std::endl
                             << " cov smoothed: " << std::endl<< ( state.hasSmoothed()   ?  state.smoothedCovariance() : cov) << std::endl
                             << std::endl;
                   }
                   }
                }
                double a_hit_chi2 = computeChi2ForHit(/*ctx.geoContext,*/
                                                      localPos,
                                                      Acts::Vector2(boundParams[Acts::eBoundLoc0], boundParams[Acts::eBoundLoc1]),
                                                      Acts::Vector2(cov(Acts::eBoundLoc0,Acts::eBoundLoc0), cov(Acts::eBoundLoc1,Acts::eBoundLoc1))
                                                      );

                std::pair<std::map<ActsFatras::Barcode,std::tuple<double, unsigned int, std::array<float, 4> > >::iterator,bool>
                   insert_result  = hit_chi2.at(traj_id).insert(std::make_pair(hitParticle.second, std::make_tuple(a_hit_chi2,1,std::array<float,4>{})));
                if (!insert_result.second && insert_result.first != hit_chi2.at(traj_id).end()) {
                   std::get<0>(insert_result.first->second) += a_hit_chi2;
                   ++std::get<1>(insert_result.first->second);
                   if (std::get<2>(insert_result.first->second).back() < a_hit_chi2) {
                      unsigned int arr_i=std::get<2>(insert_result.first->second).size();
                      for (;arr_i-->1 && std::get<2>(insert_result.first->second).at(arr_i-1)<a_hit_chi2;) {
                         std::get<2>(insert_result.first->second).at(arr_i) = std::get<2>(insert_result.first->second).at(arr_i-1);
                      }
                      std::get<2>(insert_result.first->second).at(arr_i)=a_hit_chi2;
                   }
                }

             if (debug_counter++ < 1000) {
                const auto &meas_calibrated = state.effectiveCalibrated();
                const auto &calib_cov = state.effectiveCalibratedCovariance();
                if (m_cfg.dumpDuplicates) {
                std::cout << "DEBUG hit " << hitIndex << " measurement " <<  localPos
                          << "(" << meas_calibrated << ")"
                          <<  " fit " << Acts::Vector2(boundParams[Acts::eBoundLoc0], boundParams[Acts::eBoundLoc1])
                          << " +- " << Acts::Vector2(cov(Acts::eBoundLoc0,Acts::eBoundLoc0), cov(Acts::eBoundLoc1,Acts::eBoundLoc1))
                          << " (+- " << Acts::Vector2(calib_cov(Acts::eBoundLoc0,Acts::eBoundLoc0), calib_cov(Acts::eBoundLoc1,Acts::eBoundLoc1))
                          << " )"
                          << " -> " << a_hit_chi2
                          << " (" << state.chi2() << ") "
                          << " -> " << ( insert_result.first != hit_chi2.at(traj_id).end() ? std::get<0>(insert_result.first->second) : 0.)
                          << " /  " << ( insert_result.first != hit_chi2.at(traj_id).end() ? std::get<1>(insert_result.first->second) : static_cast<unsigned int>(0))
                          << " | " << n_measurements
                          << std::endl;
                }
             }
             }

          }
          return true;
       });
       if (visualize && m_cfg.dumpDuplicates) {
          Acts::ObjVisualization3D visualizer;
          Acts::EventDataView3D::drawMultiTrajectory(visualizer,  traj.multiTrajectory(), trackTip, ctx.geoContext);
          std::stringstream file_name;
          file_name << "/tmp/traj_" << plot_i << "_" << itraj<< "_" << trackTip << ".obj";
          std::ofstream out(file_name.str().c_str());
          visualizer.write(out);
          visualizer.clear();
          ++plot_i;
       }
       unsigned int        max_count=0;
       ActsFatras::Barcode max_barcode=0;
       for (const std::pair<const ActsFatras::Barcode, Counter> &part : trajectoryBarcodeMap[ traj_id ]) {
          if (part.second.counts() > max_count) {
             max_barcode=part.first;
             max_count=part.second.counts();
          }
       }
       auto particle_iter = particles.find(max_barcode);
       double particle_pt= particle_iter != particles.end()
          ? particle_iter->transverseMomentum()
          : -std::numeric_limits<double>::max();

       if (traj_id>=traj_assoc_pt.size()) {
          traj_assoc_pt.resize(traj_id+1);
       }
       traj_assoc_pt.at(traj_id) = particle_pt;
       for (unsigned int pt_i=0; pt_i < stat_pt_bins.size(); ++pt_i) {
          if (particle_pt < stat_pt_bins.at(pt_i) ) break;
          stat.at(pt_i)[kPartPerTraj].add( trajectoryBarcodeMap[ traj_id ].size() );
          stat[pt_i]   [kMeasPerTraj].add( trajectoryToHits[traj_id].size() );
       }
    }
  }

  for (const std::pair<const ActsFatras::Barcode, std::map<unsigned int, Counter> > &particle_to_trajectories : trajectoriesPerParticle) {
     auto particle_iter = particles.find(particle_to_trajectories.first);
     double particle_pt= particle_iter != particles.end()
        ? particle_iter->transverseMomentum()
        : -std::numeric_limits<double>::max();

     std::array<unsigned int,6> n_trajectories{};
     for (const  std::pair<const unsigned int, Counter> &trajectory_hits : particle_to_trajectories.second) {
        if (trajectory_hits.second.counts()>0) { ++n_trajectories[0]; }
        if (trajectory_hits.second.counts()>=2) { ++n_trajectories[1]; }
        if (trajectory_hits.second.counts()>=3) {
           static const std::array<double,4> fractions{0., 0.5,0.9, .999};
           unsigned int dest_i=2;
           unsigned int n_trajectory_hits = trajectoryToHits.at( trajectory_hits.first ).size();
           for (const double &fraction : fractions) {
              if (trajectory_hits.second.counts() < fraction * n_trajectory_hits || dest_i >= n_trajectories.size() ) break;
              ++n_trajectories[ dest_i ];
              ++dest_i;
           }
        }
     }

     for (unsigned int pt_i=0; pt_i < stat_pt_bins.size(); ++pt_i) {
        if (particle_pt < stat_pt_bins.at(pt_i) ) break;
        stat.at(pt_i)[kTrajPerPart].add( n_trajectories[0] );
        stat.at(pt_i)[kTrajPerPartWith2Meas].add( n_trajectories[1] );
        stat.at(pt_i)[kTrajPerPartWith3Meas].add( n_trajectories[2] );
        stat.at(pt_i)[kTrajPerPartWith50Min3Meas].add( n_trajectories[3] );
        stat.at(pt_i)[kTrajPerPartWith90Min3Meas].add( n_trajectories[4] );
        stat.at(pt_i)[kTrajPerPartWith100Min3Meas].add( n_trajectories[5] );
     }
  }


  processed.reserve(trajectoryIdMap.size());
  hit_r.resize( max_hit_index+1, 0.);

  for (auto [itraj, trackTip] : trackTips) {
    const auto& traj = trajectories[itraj];
    const auto& mj = traj.multiTrajectory();
    auto trajState =
       Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);

    if (traj.hasTrajectory(trackTip)) {
       unsigned int missing_hit=0;
       std::map<ActsFatras::Barcode, Counter> particle_counts;
       std::map<unsigned int,Counter> trajectory_counts;

       std::map<std::pair<size_t, size_t>, unsigned int >::iterator
          id_iter = trajectoryId.find( std::make_pair(itraj,trackTip) );
       unsigned int traj_id=(id_iter != trajectoryId.end() ? id_iter->second : std::numeric_limits<unsigned int>::max());
       {
          std::vector<unsigned int>::const_iterator processed_iter = std::lower_bound(processed.begin(), processed.end(), traj_id);
          if (processed_iter == processed.end() || *processed_iter != traj_id) {
             processed.insert( processed_iter, traj_id);
          }
          else {
             continue;
          }
       }
       const std::vector<Index> &hits = trajectoryToHits.at(traj_id);
       std::vector<Index> all_hits(hits);
       for (Index hitIndex: all_hits) {
          Acts::Vector3 pos( globalCoords(*(m_cfg.trackingGeometry), ctx.geoContext,measurements.at(hitIndex)) );
          hit_r.at(hitIndex) = Acts::VectorHelpers::perp(pos);
       }
       std::sort(all_hits.begin(),all_hits.end());
       unsigned int n_shared_hits=0;
       for (const Index &hitIndex : hits) {

          std::map<Index, std::set< unsigned int > >::const_iterator
             hit_iter = hitTrajectoryMap.find(hitIndex);
          if (hit_iter != hitTrajectoryMap.end()) {
             if (hit_iter->second.size()>1) {
                ++n_shared_hits;
             }
             for (unsigned int sibling_traj : hit_iter->second) {
                std::vector<unsigned int>::const_iterator processed_iter = std::lower_bound(processed.begin(), processed.end(), sibling_traj);
                if (processed_iter == processed.end() || *processed_iter !=sibling_traj) {
                   processed.insert( processed_iter, sibling_traj);
                }
                const std::vector<Index> &sibling_hits = trajectoryToHits.at(sibling_traj);
                for (const Index &sibling_hit : sibling_hits) {
                   std::vector<Index>::iterator new_hit_iter = std::lower_bound(all_hits.begin(),all_hits.end(),sibling_hit);
                   if (new_hit_iter == all_hits.end() || *new_hit_iter != sibling_hit) {
                      all_hits.insert( new_hit_iter, sibling_hit);
                      Acts::Vector3 pos( globalCoords(*(m_cfg.trackingGeometry), ctx.geoContext,measurements.at(sibling_hit)) );
                      hit_r.at(sibling_hit) = Acts::VectorHelpers::perp(pos);
                   }
                }
             }
          }
       }
       const double particle_pt = traj_assoc_pt.at(traj_id);
       for (unsigned int pt_i=0; pt_i < stat_pt_bins.size(); ++pt_i) {
          if (particle_pt < stat_pt_bins.at(pt_i) ) break;
          stat.at(pt_i)[kNSharedMeasPerTraj].add(n_shared_hits);
       }

       std::sort(all_hits.begin(),all_hits.end(),
                 [&hit_r](Index hit_idx_a, Index hit_idx_b) {
                    return hit_r[hit_idx_a] < hit_r[hit_idx_b]; });

       for (const Index &hitIndex : all_hits ) {
          for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hitIndex))) {
             ++(particle_counts[hitParticle.second]);

             std::pair<std::map<ActsFatras::Barcode, unsigned int >::iterator, bool>
                ret = particle_n_hits.insert( std::make_pair(hitParticle.second, 1 ));
             if (!ret.second && ret.first != particle_n_hits.end()) {
                ++ret.first->second;
             }
          }

          std::map<Index, std::set< unsigned int > >::const_iterator
             hit_iter = hitTrajectoryMap.find(hitIndex);
          if (hit_iter != hitTrajectoryMap.end()) {
             for(unsigned int a_traj_id : hit_iter->second) {
                ++(trajectory_counts[a_traj_id]);
             }
          }
          else {
             ++missing_hit;
          }
       }

       if (particle_counts.size()>1 || trajectory_counts.size()>1) {
          if (m_cfg.dumpDuplicates) {
          std::cout << "trajectory " << std::setw(9) << itraj << " tip: " << std::setw(9) << trackTip
                    << " id " << traj_id
                    << std::endl;
          std::cout << "----------------------------- " << std::endl;
          std::cout << std::setw(9) << " " << " "
                    << std::setw(14) << " "
                    << std::setw(14) << " "
                    << std::setw(14) << " " << " ";

          for (const std::pair<const ActsFatras::Barcode, Counter> &a_particle : particle_counts) {
             // std::map<ActsFatras::Barcode, unsigned int >::const_iterator
             //    particle_iter = shortParticleId.find( a_particle.first );
             // if (particle_iter != shortParticleId.end()) {
             //    std::cout << std::setw(12) << particle_iter->second;
             // }
             // else {
                std::stringstream msg;
                msg << a_particle.first;
                std::cout << std::setw(16) << msg.str();
                //             }
          }
          std::cout << " | ";
          for (const std::pair<const unsigned int,Counter> a_trajectory : trajectory_counts) {
             std::cout << std::setw(6) << a_trajectory.first;
          }
          std::cout << std::endl;
          std::cout << std::setw(9) << " " << " "
                    << std::setw(14) << " "
                    << std::setw(14) << " "
                    << std::setw(14) << " " << " ";

          for (const std::pair<const ActsFatras::Barcode, Counter> &a_particle : particle_counts) {
             std::map<ActsFatras::Barcode, unsigned int >::const_iterator
                particle_iter = shortParticleId.find( a_particle.first );
             if (particle_iter != shortParticleId.end()) {
                 std::cout << std::setw(16) << particle_iter->second;
             }
             else {
                std::cout << std::setw(16) << " ";
             }
          }
          std::cout << " | ";
          for (const std::pair<const unsigned int,Counter> a_trajectory : trajectory_counts) {
             std::map<unsigned int, std::pair<size_t, size_t> >::const_iterator
                iter = trajectoryIdMap.find(a_trajectory.first);
             if (iter != trajectoryIdMap.end()) {
                std::cout << std::setw(6) << iter->second.first;
             }
             else {
                std::cout << std::setw(6) << " ";
             }
          }
          std::cout << std::endl;
          std::cout << std::setw(9) << " " << " "
                    << std::setw(14) << " "
                    << std::setw(14) << " "
                    << std::setw(14) << " " << " ";

          for (const std::pair<const ActsFatras::Barcode, Counter> &a_particle : particle_counts) {
             (void) a_particle; // avoid unused variable warning
             std::cout << std::setw(16) << " ";
          }
          std::cout << " | ";
          for (const std::pair<const unsigned int,Counter> a_trajectory : trajectory_counts) {
             std::map<unsigned int, std::pair<size_t, size_t> >::const_iterator
                iter = trajectoryIdMap.find(a_trajectory.first);
             if (iter != trajectoryIdMap.end()) {
                std::cout << std::setw(6) << iter->second.second;
             }
             else {
                std::cout << std::setw(6) << " ";
             }
          }
          std::cout << std::endl;

          std::cout << std::setw(9) << " " << " "
                    << std::setw(14) << " "
                    << std::setw(14) << " "
                    << std::setw(14) << " " << " " ;
          for (const std::pair<const ActsFatras::Barcode, Counter> &a_particle : particle_counts) {
             std::cout << std::setw(16) << a_particle.second.counts();
          }
          std::cout << " | ";
          for (const std::pair<const unsigned int,Counter> a_trajectory : trajectory_counts) {
             std::cout << std::setw(6) << a_trajectory.second.counts();
          }
          std::cout << std::endl << std::endl;

          // -- particles and trajectories
          for (const Index &hitIndex : all_hits ) {
             std::cout << std::setw(9) << hitIndex << " ";
             Acts::Vector3 position( globalCoords(*(m_cfg.trackingGeometry), ctx.geoContext,measurements.at(hitIndex)) );
             std::cout << std::setw(14) << Acts::VectorHelpers::perp(position)
                       << std::setw(14) << position[2]
                       << std::setw(14) << Acts::VectorHelpers::phi(position)
                       << " ";

             for (const std::pair<const ActsFatras::Barcode, Counter> &a_particle : particle_counts) {
                bool has=false;
                for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hitIndex))) {
                   if (hitParticle.second == a_particle.first) {
                      has=true;
                      break;
                   }
                }
                if (has) {
                   std::map<ActsFatras::Barcode, unsigned int >::const_iterator
                      particle_iter = shortParticleId.find( a_particle.first );
                   if (particle_iter != shortParticleId.end()) {
                      std::cout << std::setw(16) << particle_iter->second;
                   }
                   else {
                      std::stringstream msg;
                      msg << a_particle.first;
                      std::cout << std::setw(16) << msg.str();
                   }
                 }
                else {
                   std::cout << std::setw(16) << " ";
                }
             }

             std::cout << " | ";

             std::map<Index, std::set< unsigned int > >::const_iterator
                hit_iter = hitTrajectoryMap.find(hitIndex);
             if (hit_iter != hitTrajectoryMap.end()) {
                for (const std::pair<const unsigned int,Counter> a_trajectory : trajectory_counts) {
                   bool has=hit_iter->second.find(a_trajectory.first) != hit_iter->second.end();
                   if (has) {
                      unsigned int a_traj_id = a_trajectory.first;
                      bool outlier=false;
                      for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hitIndex))) {
                         std::map<ActsFatras::Barcode,std::tuple<double, unsigned int, std::array<float, 4> > >::const_iterator
                            matched_n_hits_iter = hit_chi2.at(a_traj_id).find(hitParticle.second);
                         if (matched_n_hits_iter != hit_chi2.at(a_traj_id).end()) {
                            std::map<const unsigned int, std::map<unsigned int,
                                                                  const Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy> >::const_iterator
                               hit_state_iter = hitToTrackState.find(hitIndex);
                            if (hit_state_iter != hitToTrackState.end()) {
                               std::map<unsigned int,
                                        const Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy>::const_iterator
                                  state_iter = hit_state_iter->second.find(a_trajectory.first);
                               if (state_iter != hit_state_iter->second.end()) {
                                  const auto &boundParams = state_iter->second.parameters();
                                  const auto &cov = state_iter->second.covariance();
                                  Acts::Vector2(boundParams[Acts::eBoundLoc0], boundParams[Acts::eBoundLoc1]);
                                  Acts::Vector2 localPos = std::visit(
                                                                      [](const auto& measurement) {
                                                                         auto expander = measurement.expander();
                                                                         Acts::BoundVector par = expander * measurement.parameters();
                                                                         Acts::Vector2 lpar(par[Acts::eBoundLoc0], par[Acts::eBoundLoc1]);
                                                                         return lpar;
                                                                      },
                                                                      measurements.at(hitIndex));
                                  double a_hit_chi2 = computeChi2ForHit(/*ctx.geoContext,*/
                                                                        localPos,
                                                                        Acts::Vector2(boundParams[Acts::eBoundLoc0], boundParams[Acts::eBoundLoc1]),
                                                                        Acts::Vector2(cov(Acts::eBoundLoc0,Acts::eBoundLoc0), cov(Acts::eBoundLoc1,Acts::eBoundLoc1))
                                                                        );
                                  // double traj_chi2 = std::get<0>(matched_n_hits_iter->second);
                                  // unsigned int ndf = std::get<1>(matched_n_hits_iter->second);
                                  //              if (a_hit_chi2 * ndf > 3 * traj_chi2 ) {
                                  if (a_hit_chi2 > 6  ) {
                                     outlier=true;
                                  }
                               }
                            }
                         }
                         
                      }

                      bool is_seed=false;
                      std::map< unsigned int, unsigned int>::const_iterator
                         iter = trajToSeed.find(a_trajectory.first);
                      if (iter != trajToSeed.end()) {
                         if (iter->second < protoTracks.size()) {
                            for (const auto &idx : protoTracks[iter->second]) {
                               if (idx ==  hitIndex) {
                                  is_seed=true;
                                  break;
                               }
                            }
                         }
                      }
                      std::cout << std::setw(5) << a_trajectory.first << (is_seed ? (outlier ? "x" : "*") : (outlier ? ">" : " "));
                   }
                   else {
                      std::cout << std::setw(5) << " " << " ";
                   }
                }
             }
             std::cout << std::endl;

          }

          static const std::array<const char *,9> pass_labels{"xi2phi","xi2theta","xi2q/p","xi2","match", "xi2hits", "ndf", "xi2hits'", "ndf'"};
          for (const std::pair<const ActsFatras::Barcode, Counter> &a_particle : particle_counts) {
             for (unsigned int offset_i=0; offset_i<trajectory_counts.size(); ) {
                unsigned int offset_end = std::min(static_cast<std::size_t>(offset_i+floats_per_col),trajectory_counts.size());
                for (unsigned int pass_i=0; pass_i<pass_labels.size(); ++pass_i) {
                   auto particle_iter = particles.find(a_particle.first);
                   if (particle_iter != particles.end()) {

                      if (pass_i==0) {
                         std::stringstream tmp;
                         tmp << a_particle.first;
                         std::cout << std::setw(16) << tmp.str();
                      }
                      else {
                         std::cout << std::setw(16) << " ";
                      }
                      std::cout << " " << std::setw(8) << pass_labels[pass_i] << ") ";
                      // std::cout << std::setw(14) << " "
                      //           << std::setw(14) << " "
                      //           << std::setw(14) << " "
                      //           << " ";

                      // for (const std::pair<const ActsFatras::Barcode, Counter> &a_particle_dummy : particle_counts) {
                      //    if (a_particle_dummy.second.counts()> 0) {
                      //       std::cout << std::setw(14) << " " ;
                      //    }
                      // }
                      // std::cout << " | ";
                      unsigned int traj_i=0;
                      for (const std::pair<const unsigned int,Counter> a_trajectory : trajectory_counts) {
                         if (traj_i++>= offset_i) {
                            if (traj_i > offset_end) break;

                            if (pass_i<4) {
                               Acts::Vector3 direction( particle_iter->unitDirection() );
                               std::array<double,3> ref{Acts::VectorHelpers::phi(direction),
                                  Acts::VectorHelpers::theta(direction),
                                  1./copysign(particle_iter->absoluteMomentum(),particle_iter->charge())};

                               std::map<unsigned int, std::pair<size_t, size_t> >::const_iterator traj_iter = trajectoryIdMap.find(a_trajectory.first);
                               if (traj_iter != trajectoryIdMap.end()) {
                                  const auto& a_traj = trajectories[traj_iter->second.first];
                                  const auto& fittedParameters = a_traj.trackParameters(traj_iter->second.second);

                                  if (fittedParameters.covariance()) {
                                     const auto& params = fittedParameters.parameters();
                                     const auto& cov = fittedParameters.covariance().value();
                                     std::array<Acts::BoundIndices,3> indices{
                                        Acts::eBoundPhi,
                                        Acts::eBoundTheta,
                                        Acts::eBoundQOverP};
                                     double chi2=0.;
                                     if (pass_i<indices.size()) {
                                        unsigned int idx=pass_i;
                                        chi2 += sqr( params[indices[idx] ] - ref[idx] ) / cov(indices[idx],indices[idx]);
                                        //  std::cout << params[indices[idx] ] << " - " << ref[idx] << " / " << sqrt(cov(indices[idx],indices[idx])) << " ";
                                     }
                                     else {
                                        for (unsigned int idx=0; idx < indices.size(); ++idx) {
                                           chi2 += sqr( params[indices[idx] ] - ref[idx] ) / cov(indices[idx],indices[idx]);
                                        }
                                     }
                                     std::cout << std::setw(14) << chi2;
                                  }
                                  else {
                                     std::cout << std::setw(14) << " ";
                                  }
                               }
                               else {
                                  std::cout << std::setw(14) << " ";
                               }
                            }
                            else {
                               switch (pass_i) {
                               case 4 : {
                                  unsigned int a_traj_id = a_trajectory.first;
                                  std::map<ActsFatras::Barcode, unsigned int >::const_iterator
                                     n_particle_hits_iter =  particle_n_hits.find(a_particle.first);
                                  std::map<ActsFatras::Barcode,std::tuple<double, unsigned int, std::array<float, 4> > >::const_iterator
                                     matched_n_hits_iter = hit_chi2.at(a_traj_id).find(a_particle.first);
                                  if (n_particle_hits_iter != particle_n_hits.end() && matched_n_hits_iter != hit_chi2.at(a_traj_id).end()) {
                                     std::cout << std::setw(14) << ( n_particle_hits_iter->second > 0
                                                                     ? std::get<1>(matched_n_hits_iter->second) / (1. * n_particle_hits_iter->second)
                                                                     : 0.);
                                  }
                                  else {
                                     std::cout << std::setw(14) << " ";
                                  }
                                  
                                  break;
                               }
                               case 5 : {
                                  unsigned int a_traj_id = a_trajectory.first;
                                  std::map<ActsFatras::Barcode,std::tuple<double, unsigned int, std::array<float, 4> > >::const_iterator
                                     matched_n_hits_iter = hit_chi2.at(a_traj_id).find(a_particle.first);
                                  if (matched_n_hits_iter != hit_chi2.at(a_traj_id).end()) {
                                     std::cout << std::setw(14) << ( std::get<1>(matched_n_hits_iter->second) > 0
                                                                     ? std::get<0>(matched_n_hits_iter->second) / std::get<1>(matched_n_hits_iter->second)
                                                                     : 0. );
                                  }
                                  else {
                                     std::cout << std::setw(14) << " ";
                                  }
                                  break;
                               }
                               case 6 : {
                                  unsigned int a_traj_id = a_trajectory.first;
                                  std::map<ActsFatras::Barcode,std::tuple<double, unsigned int, std::array<float, 4> > >::const_iterator
                                     matched_n_hits_iter = hit_chi2.at(a_traj_id).find(a_particle.first);
                                  if (matched_n_hits_iter != hit_chi2.at(a_traj_id).end()) {
                                     std::cout << std::setw(14) << std::get<1>(matched_n_hits_iter->second);
                                  }
                                  else {
                                     std::cout << std::setw(14) << " ";
                                  }
                                  break;
                               }
                               case 7 : {
                                  unsigned int a_traj_id = a_trajectory.first;
                                  std::map<ActsFatras::Barcode,std::tuple<double, unsigned int, std::array<float, 4> > >::const_iterator
                                     matched_n_hits_iter = hit_chi2.at(a_traj_id).find(a_particle.first);
                                  if (matched_n_hits_iter != hit_chi2.at(a_traj_id).end()) {
                                     std::pair<double, unsigned int> mod_chi2 = modifiedChi2(matched_n_hits_iter->second);
                                     std::cout << std::setw(14) << ( mod_chi2.second > 1
                                                                     ? (mod_chi2.first / mod_chi2.second)
                                                                     : mod_chi2.first );
                                  }
                                  else {
                                     std::cout << std::setw(14) << " ";
                                  }
                                  break;
                               }
                               case 8 : {
                                  unsigned int a_traj_id = a_trajectory.first;
                                  std::map<ActsFatras::Barcode,std::tuple<double, unsigned int, std::array<float, 4> > >::const_iterator
                                     matched_n_hits_iter = hit_chi2.at(a_traj_id).find(a_particle.first);
                                  if (matched_n_hits_iter != hit_chi2.at(a_traj_id).end()) {
                                     std::pair<double, unsigned int> mod_chi2 = modifiedChi2(matched_n_hits_iter->second);
                                     std::cout << std::setw(14) << mod_chi2.second ;
                                  }
                                  else {
                                     std::cout << std::setw(14) << " ";
                                  }
                                  break;
                               }
                               }
                                  
                            }
                         }
                      }
                   }
                   std::cout << std::endl;
                }
                offset_i = offset_end;
             }
          }
          }


          // -- seed
          if (m_cfg.dumpDuplicates) {
          std::vector<std::vector<Index> > seed_cluster_idx;
          for (unsigned int offset_i=0;offset_i<trajectory_counts.size();) {
             unsigned int traj_end  = std::min(static_cast<unsigned int>(trajectory_counts.size()), offset_i+floats_per_col);
             seed_cluster_idx.clear();
             seed_cluster_idx.reserve(floats_per_col);
             bool have_idx=true;
             for (unsigned pass_i=0; have_idx; ++pass_i) {
                unsigned int counter_i=0;
                unsigned int elm_i=0;
                have_idx=(pass_i==0);
                std::cout << std::setw(16) << " " << " " << std::setw(8) << "seedhits" << "  "; 
                for (const std::pair<const unsigned int,Counter> trajectory : trajectory_counts) {
                   if (counter_i++>=offset_i) {
                      if (pass_i==0) {
                         seed_cluster_idx.emplace_back();
                         std::map< unsigned int, unsigned int>::const_iterator
                            iter = trajToSeed.find(trajectory.first);
                         std::cout << std::setw(14);
                         if (iter != trajToSeed.end()) {
                            std::cout << iter->second;
                            if (iter->second < protoTracks.size()) {
                               for (const auto &idx : protoTracks[iter->second]) {
                                  seed_cluster_idx.back().push_back(idx);
                               }
                            }
                         }
                         else {
                            std::cout << " ";
                         }
                         if (counter_i>=traj_end) { break; }
                      }
                      else {
                         std::cout << std::setw(14);
                         if (pass_i <= seed_cluster_idx.at(elm_i).size()) {
                            std::cout << seed_cluster_idx[elm_i][pass_i-1];
                            have_idx=true;
                         }
                         else {
                            std::cout << " ";
                         }
                         if (counter_i>=traj_end) { break; }
                      }
                      ++elm_i;
                   }
                }
                std::cout << std::endl;
             }
             offset_i = traj_end;
             std::cout << std::endl;
          }
          std::cout << std::endl;
          }

          // -- particles
          if (m_cfg.dumpDuplicates) {
             std::cout << std::setw(4) << " ";
             for (const std::pair<const ActsFatras::Barcode, Counter> &a_particle : particle_counts) {
                if (a_particle.second.counts()> 0) {
                   auto particle_iter = particles.find(a_particle.first);
                   if (particle_iter != particles.end()) {
                      std::stringstream barcode_str;
                      barcode_str << particle_iter->particleId();
                      std::cout << std::setw(14) << barcode_str.str();
                   }
                   else {
                      std::stringstream barcode_str;
                      barcode_str << a_particle.first;
                      std::cout << std::setw(14) << barcode_str.str();
                   }
                }
             }
             std::cout << std::endl;
             std::vector< std::array<double,7> > particle_param;
             particle_param.reserve( particle_counts.size() );
             for (const std::pair<const ActsFatras::Barcode, Counter> &a_particle : particle_counts) {
                if (a_particle.second.counts()> 0) {
                   auto particle_iter = particles.find(a_particle.first);
                   if (particle_iter != particles.end()) {
                      Acts::Vector3 position( particle_iter->position() );
                      Acts::Vector3 direction( particle_iter->unitDirection() );
                      particle_param.emplace_back(std::array<double,7>({Acts::VectorHelpers::phi(position),
                                                  Acts::VectorHelpers::theta(position),
                                                  Acts::VectorHelpers::perp(position),
                                                  position[2],
                                                  Acts::VectorHelpers::phi(direction),
                                                  Acts::VectorHelpers::theta(direction),
                                                  std::copysign(particle_iter->transverseMomentum(), particle_iter->charge())} ));
                   }
                }
             }
             dumpParameters("PART", particle_param, floats_per_col);
          }
          // -- parameters
          if (m_cfg.dumpDuplicates) {
             std::vector<double> traj_pt;
             dumpPerigeeParameters("ORIG", ctx, trajectory_counts, trajectoryIdMap, trajectories,traj_pt, floats_per_col);
          }
          for (const Index &hitIndex : all_hits ) {
             std::map<const unsigned int, std::map<unsigned int,
                                                   const Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy> >::const_iterator
                hit_state_iter = hitToTrackState.find(hitIndex);
             if (hit_state_iter != hitToTrackState.end()) {
                std::vector<double> traj_pt;
                std::stringstream header;
                header << hitIndex;
                dumpParameters(header.str(), ctx,trajectory_counts, hit_state_iter->second, traj_pt, floats_per_col, m_cfg.dumpDuplicates);
                unsigned idx=0;
                for (const std::pair<const unsigned int,Counter> &trajectory : trajectory_counts) {
                   double a_traj_pt=traj_pt.at(idx);
                   if (a_traj_pt >= 0.) {
                      const double a_particle_pt = traj_assoc_pt.at(trajectory.first);
                      for (unsigned int pt_i=0; pt_i < stat_pt_bins.size(); ++pt_i) {
                         if (a_particle_pt < stat_pt_bins.at(pt_i) ) break;
                         stat.at(pt_i)[kMeasPt].add(a_traj_pt);
                      }
                   }
                   ++idx;
                }
             }
          }
          
       }
    }
  }

  if (m_stat.empty()) {
     std::lock_guard<std::mutex> stat_lock(m_mutex);
     m_stat.resize(stat_pt_bins.size());
     m_ptBins.clear();
     m_ptBins.reserve(stat_pt_bins.size());
     std::copy(stat_pt_bins.begin(),stat_pt_bins.end(),std::back_inserter(m_ptBins));
     std::copy(stat.begin(), stat.end(), std::back_inserter(m_stat));
     m_stat.resize(m_stat.size());
     for (unsigned int pt_i=0; pt_i < stat.size(); ++pt_i) {
        m_stat.at(pt_i).resize(stat[pt_i].size());
        for (unsigned int stat_i=0; stat_i < stat[pt_i].size(); ++stat_i) {
           m_stat.at(pt_i).at(stat_i) = stat.at(pt_i).at(stat_i);
        }
     }
  }
  else {
     std::lock_guard<std::mutex> stat_lock(m_mutex);
     for (unsigned int pt_i=0; pt_i < stat.size(); ++pt_i) {
        for (unsigned int stat_i=0; stat_i < stat[pt_i].size(); ++stat_i) {
           m_stat.at(pt_i).at(stat_i) += stat.at(pt_i).at(stat_i);
        }
     }
  }
  dumpStatPtBins( stat_pt_bins, stat);

  // Loop over all trajectories
  for (std::size_t iTraj = 0; iTraj < trajectories.size(); ++iTraj) {
    const auto& traj = trajectories[iTraj];
    const auto& mj = traj.multiTrajectory();
    for (auto trackTip : traj.tips()) {
      auto trajState =
          Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);

      // Reco track selection
      //@TODO: add interface for applying others cuts on reco tracks:
      // -> pT, d0, z0, detector-specific hits/holes number cut
      if (trajState.nMeasurements < m_cfg.nMeasurementsMin) {
        continue;
      }
      // Check if the reco track has fitted track parameters
      if (not traj.hasTrackParameters(trackTip)) {
        ACTS_WARNING(
            "No fitted track parameters for trajectory with entry index = "
            << trackTip);
        continue;
      }
      const auto& fittedParameters = traj.trackParameters(trackTip);
      // Requirement on the pT of the track
      const auto& momentum = fittedParameters.momentum();
      const auto pT = perp(momentum);
      if (pT < m_cfg.ptMin) {
        continue;
      }
      // Fill the trajectory summary info
      m_trackSummaryPlotTool.fill(m_trackSummaryPlotCache, fittedParameters,
                                  trajState.nStates, trajState.nMeasurements,
                                  trajState.nOutliers, trajState.nHoles,
                                  trajState.nSharedHits);

      // Get the majority truth particle to this track
      identifyContributingParticles(hitParticlesMap, traj, trackTip,
                                    particleHitCounts);
      if (particleHitCounts.empty()) {
        ACTS_WARNING(
            "No truth particle associated with this trajectory with entry "
            "index = "
            << trackTip);
        continue;
      }
      // Get the majority particleId and majority particle counts
      // Note that the majority particle might be not in the truth seeds
      // collection
      ActsFatras::Barcode majorityParticleId =
          particleHitCounts.front().particleId;
      size_t nMajorityHits = particleHitCounts.front().hitCount;

      // Check if the trajectory is matched with truth.
      // If not, it will be classified as 'fake'
      bool isFake = false;
      if (nMajorityHits * 1. / trajState.nMeasurements >=
          m_cfg.truthMatchProbMin) {
        matched[majorityParticleId].push_back(
            {nMajorityHits, fittedParameters});
      } else {
        isFake = true;
        unmatched[majorityParticleId]++;
      }
      // Fill fake rate plots
      m_fakeRatePlotTool.fill(m_fakeRatePlotCache, fittedParameters, isFake);

      // Use neural network classification for duplication rate plots
      // Currently, the network used for this example can only handle
      // good/duplicate classification, so need to manually exclude fake tracks
      if (m_cfg.duplicatedPredictor && !isFake) {
        inputFeatures[0] = trajState.nMeasurements;
        inputFeatures[1] = trajState.nOutliers;
        inputFeatures[2] = trajState.chi2Sum * 1.0 / trajState.NDF;
        // predict if current trajectory is 'duplicate'
        bool isDuplicated = m_cfg.duplicatedPredictor(inputFeatures);
        // Add to number of duplicated particles
        if (isDuplicated) {
          m_nTotalDuplicateTracks++;
        }
        // Fill the duplication rate
        m_duplicationPlotTool.fill(m_duplicationPlotCache, fittedParameters,
                                   isDuplicated);
      }
      // Counting number of total trajectories
      m_nTotalTracks++;
    }
  }

  // Use truth-based classification for duplication rate plots
  if (!m_cfg.duplicatedPredictor) {
    // Loop over all truth-matched reco tracks for duplication rate plots
    for (auto& [particleId, matchedTracks] : matched) {
      // Sort the reco tracks matched to this particle by the number of majority
      // hits
      std::sort(matchedTracks.begin(), matchedTracks.end(),
                [](const RecoTrackInfo& lhs, const RecoTrackInfo& rhs) {
                  return lhs.first > rhs.first;
                });
      for (size_t itrack = 0; itrack < matchedTracks.size(); itrack++) {
        const auto& [nMajorityHits, fittedParameters] =
            matchedTracks.at(itrack);
        // The tracks with maximum number of majority hits is taken as the
        // 'real' track; others are as 'duplicated'
        bool isDuplicated = (itrack != 0);
        // the track is associated to the same particle
        if (isDuplicated) {
          m_nTotalDuplicateTracks++;
        }
        // Fill the duplication rate
        m_duplicationPlotTool.fill(m_duplicationPlotCache, fittedParameters,
                                   isDuplicated);
      }
    }
  }

  // Loop over all truth particle seeds for efficiency plots and reco details.
  // These are filled w.r.t. truth particle seed info
  for (const auto& particle : particles) {
    if (particle.transverseMomentum() < m_cfg.ptMin) {
      continue;
    }
    auto particleId = particle.particleId();
    // Investigate the truth-matched tracks
    size_t nMatchedTracks = 0;
    bool isReconstructed = false;
    auto imatched = matched.find(particleId);
    if (imatched != matched.end()) {
      nMatchedTracks = imatched->second.size();
      // Add number for total matched tracks here
      m_nTotalMatchedTracks += nMatchedTracks;
      m_nTotalMatchedParticles += 1;
      // Check if the particle has more than one matched track for the duplicate
      // rate
      if (nMatchedTracks > 1) {
        m_nTotalDuplicateParticles += 1;
      }
      isReconstructed = true;
    }
    // Fill efficiency plots
    m_effPlotTool.fill(m_effPlotCache, particle, isReconstructed);
    // Fill number of duplicated tracks for this particle
    m_duplicationPlotTool.fill(m_duplicationPlotCache, particle,
                               nMatchedTracks - 1);

    // Investigate the fake (i.e. truth-unmatched) tracks
    size_t nFakeTracks = 0;
    auto ifake = unmatched.find(particleId);
    if (ifake != unmatched.end()) {
      nFakeTracks = ifake->second;
      m_nTotalFakeTracks += nFakeTracks;
      // unmatched is a map of majority particle id to # of tracks associated
      // with that particle
      m_nTotalFakeParticles += 1;
    }
    // Fill number of reconstructed/truth-matched/fake tracks for this particle
    m_fakeRatePlotTool.fill(m_fakeRatePlotCache, particle, nMatchedTracks,
                            nFakeTracks);
    m_nTotalParticles += 1;
  }  // end all truth particles

  return ProcessCode::SUCCESS;
}
