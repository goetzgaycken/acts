#pragma once

#include <limits>
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"

namespace ActsExamples {

/// "concept" for a seed filter
/// The seed filter allows to filter seeds using information from previously
/// found trajectories.
template <typename traj_t>
class SeedPassThroughFilter {
public:
   using TrajectoryIDType = std::size_t;

   bool filterSeed([[maybe_unused]] std::size_t iseed) { return false; }
   /// update the seed filter using the information from the latest trajectory
   /// @return the next trajector ID, where a unique id is only assigned to trajectories with e.g. measurements
   TrajectoryIDType update(TrajectoryIDType traj_id,
                           [[maybe_unused]] const traj_t &multi_trajectory,
                           const std::vector<Acts::MultiTrajectoryTraits::IndexType> &tips) {
      return traj_id + tips.size();
   }
};

class SeedFilter {
public:
   using TrajectoryIDType = unsigned short;
   constexpr static unsigned int NTrajectoriesPerHit = 16;

   template <typename T_Index, unsigned int N>
   constexpr std::array<T_Index, N> initialTrajectoriesPerHit(T_Index def_val) {
      std::array<T_Index, N> tmp;
      std::fill(tmp.begin(),tmp.end(), def_val);
      return tmp;
   }


   SeedFilter(const ActsExamples::ProtoTrackContainer &proto_tracks,
              std::size_t maxMeasurementIndex)
      : protoTracks(proto_tracks)
   {
      trajectories_per_hit.resize(maxMeasurementIndex,
                                  initialTrajectoriesPerHit<TrajectoryIDType,
                                                            NTrajectoriesPerHit>(std::numeric_limits<TrajectoryIDType>::max()));
   }

   bool filterSeed(std::size_t iseed) {
      unsigned int n_traj=0;
      last_seed=iseed;
      std::array<unsigned short,NTrajectoriesPerHit> shared_traj;
      for (const auto &hitIndex : protoTracks[iseed]) {
         if (n_traj==0) {
            for (auto a_traj_id : trajectories_per_hit.at(hitIndex) ) {
               if (a_traj_id == std::numeric_limits<TrajectoryIDType>::max()) break;
               if (n_traj>=shared_traj.size()) {
                  break;
               }
               shared_traj[n_traj++]=a_traj_id;
            }
            if (n_traj==0) return false;
         }
         else {
            // remove all trjectories from the list which do not contain this seed hit;
            for (unsigned int a_traj_i=0; a_traj_i < n_traj; ++a_traj_i) {
               // @TODO for the likely small number of trajectories per hit, lineare search
               //       is likely faster. lower_bound is just used for convenience.
               auto traj_iter = std::lower_bound(trajectories_per_hit.at(hitIndex).begin(),
                                                 trajectories_per_hit.at(hitIndex).end(),
                                                 shared_traj[a_traj_i]);
               if (traj_iter == trajectories_per_hit.at(hitIndex).end()) {
                  --n_traj;
                  if (n_traj==0) return false;
                  for (;a_traj_i < n_traj; ++a_traj_i) {
                     shared_traj[a_traj_i]=shared_traj[a_traj_i+1];
                  }
                  break;
               }
            }
         }
      }
      if (n_traj>0) {
         std::cout << "DEBUG SeedFilter seed " <<  iseed << " (future traj id: " << last_traj <<  ")" << " filtered out because its hits are on " << n_traj << " trajectories :";
         for (unsigned int i=0; i<n_traj; ++i) {
            std::cout << " " << shared_traj[i];
         }
         for (const auto &hitIndex : protoTracks[iseed]) {
            std::cout << ", " << hitIndex << " :" ;
            for(auto a_traj_id : trajectories_per_hit.at(hitIndex) ) {
               if (a_traj_id == std::numeric_limits<TrajectoryIDType>::max()) break;
               std::cout <<  " " << a_traj_id;
            }
         }
         std::cout << std::endl;
      }
      return  n_traj>0;
   }
   TrajectoryIDType update(TrajectoryIDType traj_id,
               const Acts::MultiTrajectory<Acts::VectorMultiTrajectory> &multi_trajectory,
               const std::vector<Acts::MultiTrajectoryTraits::IndexType> &tips) {
      // mark all hits of the trajectory as being used
      for (auto trackTip : tips) {
         bool has_measurements=false;
         multi_trajectory.visitBackwards(trackTip, [this,
                                                    traj_id,
                                                    &has_measurements](const auto& state) {
            // no truth info with non-measurement state
            if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
               return true;
            }

            // register all particles that generated this hit
            const auto& sl = state.uncalibratedSourceLink().template get<ActsExamples::IndexSourceLink>();
            auto hitIndex = sl.index();
            if (hitIndex>=this->trajectories_per_hit.size()) {
               this->trajectories_per_hit.resize(hitIndex+1);
            }
            auto &hit_trajectories = this->trajectories_per_hit.at(hitIndex);

            // trajectories are added in increasing order, thus once the elment value
            // is larger than the ID of the current trajectory, an unused element is reached.
            unsigned int counter=0;
            for (auto &elm : hit_trajectories) {
               ++counter;
               if (elm>traj_id) {
                  elm = traj_id;
                  break;
               }
            }
            this->max_counter=std::max(this->max_counter,counter);
            if (hit_trajectories.back() != std::numeric_limits<unsigned short>::max()) {
               ++this->n_errors;
            }
            has_measurements=true;
            return true;
         });

         if (has_measurements) {
            std::cout << "DEBUG SeedFilter register new trajectory " <<  traj_id << " for seed " << last_seed << std::endl;
            ++traj_id;
         }
      }
      last_seed=std::numeric_limits<unsigned int>::max();
      last_traj=traj_id;
      return traj_id;
   }

   const ActsExamples::ProtoTrackContainer &protoTracks;
   std::vector< std::array<TrajectoryIDType, NTrajectoriesPerHit> > trajectories_per_hit;
   unsigned int n_errors = 0;
   unsigned int max_counter=0;
   unsigned int last_traj=0;
   unsigned int last_seed=std::numeric_limits<unsigned int>::max();
};
}
