// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Result.hpp"

#include <vector>

namespace Acts {

/// Return type of the `BranchStopper` delegate for the
/// CombinatorialKalmanFilter
enum class CombinatorialKalmanFilterBranchStopperResult {
  Continue,
  StopAndDrop,
  StopAndKeep,
};

/// Extension struct which holds the delegates to customize the CKF behavior
template <typename track_container_t>
struct CombinatorialKalmanFilterExtensions {
  using traj_t = typename track_container_t::TrackStateContainerBackend;
  using TrackProxy = typename track_container_t::TrackProxy;
  using TrackStateProxy = typename track_container_t::TrackStateProxy;

  using BranchStopperResult = CombinatorialKalmanFilterBranchStopperResult;

  using Calibrator = typename KalmanFitterExtensions<traj_t>::Calibrator;
  using Updater = typename KalmanFitterExtensions<traj_t>::Updater;
  using BranchStopper =
      Delegate<BranchStopperResult(const TrackProxy&, const TrackStateProxy&)>;

  /// The Calibrator is a dedicated calibration algorithm that allows to
  /// calibrate measurements using track information, this could be e.g. sagging
  /// for wires, module deformations, etc.
  Calibrator calibrator{
      DelegateFuncTag<detail::voidFitterCalibrator<traj_t>>{}};

  /// The updater incorporates measurement information into the track parameters
  Updater updater{DelegateFuncTag<detail::voidFitterUpdater<traj_t>>{}};

  /// The branch stopper is called during the filtering by the Actor.
  BranchStopper branchStopper{DelegateFuncTag<voidBranchStopper>{}};

 private:

  /// Default branch stopper which will never stop
  /// @return false
  static BranchStopperResult voidBranchStopper(
      const TrackProxy& /*track*/, const TrackStateProxy& /*trackState*/) {
    return BranchStopperResult::Continue;
  }
};
}
