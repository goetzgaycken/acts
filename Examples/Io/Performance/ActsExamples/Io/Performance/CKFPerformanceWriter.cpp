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

#include <numeric>
#include <stdexcept>

#include <unistd.h>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TVectorF.h>

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

  std::cout << "particles " << hitParticlesMap.size() << std::endl;
  std::map<std::pair<size_t, size_t>, unsigned int >              trajectoryId;
  std::map<unsigned int, std::pair<size_t, size_t> >              trajectoryIdMap;
  std::map<ActsFatras::Barcode, unsigned int >                    shortParticleId;
  std::map<Index, std::set< unsigned int > >                      hitTrajectoryMap;
  std::map<unsigned int, std::vector<Index> >                     trajectoryToHits;
  std::map<unsigned int, std::map<ActsFatras::Barcode, Counter> > trajectoryBarcodeMap;
  std::map<unsigned int, std::map<unsigned int,
                                  const Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy> > hitToTrackState;
  std::map<unsigned int, unsigned int> trajToSeed;
  std::vector<unsigned int> processed;

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
       trajToSeed.insert( std::make_pair(traj_id, seed_i));
       traj.multiTrajectory().visitBackwards(trackTip, [&](const auto& state) {
          // no truth info with non-measurement state
          if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
             return true;
          }

          // register all particles that generated this hit
          const auto& sl = static_cast<const IndexSourceLink&>(state.uncalibrated());
          auto hitIndex = sl.index();
          hitToTrackState[hitIndex].insert( std::make_pair(traj_id, state) );
          hitTrajectoryMap[hitIndex].insert(traj_id);
          trajectoryToHits[traj_id].push_back(hitIndex);
          for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hitIndex))) {
             /*std::pair<std::map<ActsFatras::Barcode, unsigned int >::const_iterator, bool >
               particle_insert =*/ shortParticleId.insert( std::make_pair(hitParticle.second, shortParticleId.size()));
             ++(trajectoryBarcodeMap[ traj_id ][hitParticle.second]);
          }
          return true;
       });
    }
  }
  processed.reserve(trajectoryIdMap.size());
  
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
       std::sort(all_hits.begin(),all_hits.end());
       for (const Index &hitIndex : hits) {
          std::map<Index, std::set< unsigned int > >::const_iterator
             hit_iter = hitTrajectoryMap.find(hitIndex);
          if (hit_iter != hitTrajectoryMap.end()) {
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
                   }
                }
             }
          }
       }
       

       for (const Index &hitIndex : all_hits ) {

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
          for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hitIndex))) {
             ++(particle_counts[hitParticle.second]);
          }
       }

       if (particle_counts.size()>1 || trajectory_counts.size()>1) {
          
          std::cout << "trajectory " << std::setw(9) << itraj << " tip: " << std::setw(9) << trackTip
                    << " id " << traj_id
                    << std::endl;
          std::cout << "----------------------------- " << std::endl;
          std::cout << std::setw(9) << " " << " ";
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
          std::cout << std::setw(9) << " " << " ";
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
                std::cout << std::setw(4) << iter->second.first;
                
             }
             else {
                std::cout << std::setw(4) << " ";
             }
          }
          std::cout << std::endl;
          std::cout << std::setw(9) << " " << " ";
          for (const std::pair<const ActsFatras::Barcode, Counter> &a_particle : particle_counts) {
             (void) a_particle; // avoid unused variable warning
             std::cout << std::setw(16) << " ";
          }
          std::cout << " | ";
          for (const std::pair<const unsigned int,Counter> a_trajectory : trajectory_counts) {
             std::map<unsigned int, std::pair<size_t, size_t> >::const_iterator
                iter = trajectoryIdMap.find(a_trajectory.first);
             if (iter != trajectoryIdMap.end()) {
                std::cout << std::setw(4) << iter->second.second;
             }
             else {
                std::cout << std::setw(4) << " ";
             }
          }
          std::cout << std::endl;
          
          std::cout << std::setw(9) << " " << " ";
          for (const std::pair<const ActsFatras::Barcode, Counter> &a_particle : particle_counts) {
             std::cout << std::setw(16) << a_particle.second.counts();
          }
          std::cout << " | ";
          for (const std::pair<const unsigned int,Counter> a_trajectory : trajectory_counts) {
             std::cout << std::setw(4) << a_trajectory.second.counts();
          }
          std::cout << std::endl << std::endl;

          // -- particles and trajectories
          for (const Index &hitIndex : all_hits ) {
             std::cout << std::setw(9) << hitIndex << " ";
             
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
                      std::cout << std::setw(4) << a_trajectory.first << (is_seed ? "*" : " ");
                   }
                   else {
                      std::cout << std::setw(4) << " " << " ";
                   }
                }
             }
             std::cout << std::endl;
 
          }
          
          // -- seed 
          std::vector<std::vector<Index> > seed_cluster_idx;
          for (unsigned int offset_i=0;offset_i<trajectory_counts.size();) {
             unsigned int traj_end  = std::min(static_cast<unsigned int>(trajectory_counts.size()), offset_i+10);
             seed_cluster_idx.clear();
             seed_cluster_idx.reserve(10);
             bool have_idx=true;
             for (unsigned pass_i=0; have_idx; ++pass_i) {
                unsigned int counter_i=0;
                unsigned int elm_i=0;
                have_idx=(pass_i==0);
                std::cout << std::setw(4) << " " << "  ";
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

          // -- parameters
          for( const std::pair<const unsigned int, std::map<unsigned int,
                                                           const Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy> >
                  &hit_states : hitToTrackState ) {
             std::vector<std::array<double,6> > trajectory_params;
             trajectory_params.reserve( trajectory_counts.size() );
             bool has_param=false;
             for (const std::pair<const unsigned int,Counter> trajectory : trajectory_counts) {
                std::map<unsigned int,const Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy>::const_iterator
                   traj_iter = hit_states.second.find(trajectory.first);
                if (traj_iter != hit_states.second.end()) {
                   const auto &boundParams = traj_iter->second.parameters();
                   Acts::Vector3 direction = Acts::makeDirectionUnitFromPhiTheta(boundParams[Acts::eBoundPhi],
                                                                                 boundParams[Acts::eBoundTheta]);

                   Acts::Vector2 local(boundParams[Acts::eBoundLoc0], boundParams[Acts::eBoundLoc1]);
                   Acts::Vector3 position(traj_iter->second.referenceSurface().localToGlobal(ctx.geoContext, local, direction));
                   trajectory_params.emplace_back(std::array<double,6>{ position[Acts::ePos0],
                         position[Acts::ePos1],
                         position[Acts::ePos2],
                         boundParams[Acts::eBoundPhi],
                         boundParams[Acts::eBoundTheta],
                         boundParams[Acts::eBoundQOverP]});
                   has_param=true;
                }
                else {
                   trajectory_params.emplace_back(std::array<double,6>{0.,0.,0.,0.,0.,0.});
                }
             }
             if (has_param) {
                for (unsigned int offset_i=0;offset_i<trajectory_params.size();) {
                   unsigned int traj_end  = std::min(static_cast<unsigned int>(trajectory_params.size()), offset_i+10);
                   for (unsigned int param_i=0; param_i<6; ++param_i) {
                      if (param_i==0 && offset_i==0) {
                         std::cout << std::setw(4) << hit_states.first << ") ";
                      }
                      else {
                         std::cout << std::setw(4) << " " << "  ";
                      }
                      for (unsigned int traj_i=offset_i; traj_i<traj_end;++traj_i) {
                         const std::array<double, 6> &param = trajectory_params[traj_i];
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
       }
    }
  }

  
  // Loop over all trajectories
  for (std::size_t iTraj = 0; iTraj < trajectories.size(); ++iTraj) {
    const auto& traj = trajectories[iTraj];
    const auto& mj = traj.multiTrajectory();
    for (auto trackTip : traj.tips()) {
      // @TODO: Switch to using this directly from the track
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
