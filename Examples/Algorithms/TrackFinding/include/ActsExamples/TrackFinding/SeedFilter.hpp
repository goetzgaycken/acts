#pragma once

#include <limits>
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"

namespace ActsExample {

class SeedFilter {
public:
   TrajectoryIDType = using unsigned short;
   constexpr unsigned int NTrajectoriesPerHit = 16;

   template <typename T_Index, unsigned int N>
   constexpr std::array<T_Index, N> initialTrajectoriesPerHit(T_Index def_val) {
      std::array<T_Index, N> tmp;
      std::fill(tmp.begin(),tmp.end(), def_val);
      return tmp;
   }


   SeedFilter(const ProtoTrackContainer &proto_track,
              std::size_t maxMeasurementIndex)
      : protoTrack(proto_track)
   {
      trajectories_per_hit.resize(maxMeasurementIndex,
                                  initialTrajectoriesPerHit<TrajectoryIDType,
                                                            NTrajectoriesPerHit>(std::numeric_limits<TrajectoryIDType>::max()));
   }

   bool filterSeed(std::size_t iseed) {
      unsigned int n_traj=0;
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
                  for (;a_traj_i < n_traj; ++a_traj_i) {
                     shared_traj[a_traj_i]=shared_traj[a_traj_i+1];
                  }
                  break;
               }
            }
         }
      }
      return  n_traj>0;
   }
   bool update(TrajectoryIDType traj_id,
               Acts::MultiTrajectory<Acts::VectorMultiTrajectory> &multi_trajectory,
               std::vector<Acts::MultiTrajectoryTraits::IndexType> &tips) {
      // mark all hits of the trajectory as being used
      for (auto trackTip : tips) {
         bool has_measurements=false;
         multi_trajectory.visitBackwards(trackTip, [&this,
                                                    traj_id,
                                                    &has_measurements](const auto& state) {
            // no truth info with non-measurement state
            if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
               return true;
            }

            // register all particles that generated this hit
            const auto& sl = static_cast<const ActsExamples::IndexSourceLink&>(state.uncalibrated());
            auto hitIndex = sl.index();
            if (hitIndex>=trajectories_per_hit.size()) {
               this->trajectories_per_hit.resize(hitIndex+1);
            }
            auto &hit_trajectories = this->trajectories_per_hit.at(hitIndex);

            // trajectories are added in increasing order, thus once the elment value
            // is larger than the ID of the current trajectory, an unused element is reached.
            unsigned int counter=0;
            for (auto &elm : hit_trajectories) {
               if (elm>traj_id) {
                  elm = traj_id;
                  break;
               }
               ++counter;
            }
            this->max_counter=std::max(this->max_counter,counter);
            if (hit_trajectories.beck() != std::numeric_limits<unsigned short>::max()) {
               ++this->n_errors;
            }
            has_measurements=true;
         });
         
         if (has_measurements) {
            ++traj_id;
         }
      }
      return traj_id;
   }

   const ProtoTrackContainer &protoTrack;
   std::vector< std::array<TrajectoryIDType, NTrajectoriesPerHit> > trajectories_per_hit;
   unsigned int n_errors = 0;
   unsigned int max_counter=0;
};
}
