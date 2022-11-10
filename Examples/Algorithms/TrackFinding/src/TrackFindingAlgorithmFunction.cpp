// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedFilter.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <random>
#include <stdexcept>

namespace {

using Updater = Acts::GainMatrixUpdater;
using Smoother = Acts::GainMatrixSmoother;

using Stepper = Acts::EigenStepper<>;
using Navigator = Acts::Navigator;
using Propagator = Acts::Propagator<Stepper, Navigator>;
using CKF =
    Acts::CombinatorialKalmanFilter<Propagator, Acts::VectorMultiTrajectory>;

struct TrackFinderFunctionImpl
    : public ActsExamples::TrackFindingAlgorithm::TrackFinderFunction {
  CKF trackFinder;

  TrackFinderFunctionImpl(CKF&& f) : trackFinder(std::move(f)) {}

  ActsExamples::TrackFindingAlgorithm::TrackFinderResult operator()(
      const ActsExamples::TrackParametersContainer& initialParameters,
      const ActsExamples::TrackFindingAlgorithm::TrackFinderOptions& options,
      ActsExamples::SeedFilter&      seed_filter,
      std::shared_ptr<Acts::VectorMultiTrajectory> trajectory) const override {
    return trackFinder.findTracks(initialParameters, options,
                                  seed_filter,
                                  std::move(trajectory));
  };
};

}  // namespace

std::shared_ptr<ActsExamples::TrackFindingAlgorithm::TrackFinderFunction>
ActsExamples::TrackFindingAlgorithm::makeTrackFinderFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField) {
  Stepper stepper(std::move(magneticField));
  Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Navigator navigator(cfg);
  Propagator propagator(std::move(stepper), std::move(navigator));
  CKF trackFinder(std::move(propagator));

  // build the track finder functions. owns the track finder object.
  return std::make_shared<TrackFinderFunctionImpl>(std::move(trackFinder));
}
