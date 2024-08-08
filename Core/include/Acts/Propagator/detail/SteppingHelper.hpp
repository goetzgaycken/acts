// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <sstream>
#include <iostream>

#include <limits>

namespace Acts::detail {

/// Update surface status - Single component
///
/// This method intersect the provided surface and update the navigation
/// step estimation accordingly (hence it changes the state). It also
/// returns the status of the intersection to trigger onSurface in case
/// the surface is reached.
///
/// @param state [in,out] The stepping state (thread-local cache)
/// @param surface [in] The surface provided
/// @param boundaryTolerance [in] The boundary check for this status update
template <typename stepper_t>
Acts::Intersection3D::Status updateSingleSurfaceStatus(
    const stepper_t& stepper, typename stepper_t::State& state,
    const Surface& surface, std::uint8_t index, Direction navDir,
    const BoundaryTolerance& boundaryTolerance, ActsScalar surfaceTolerance,
    const Logger& logger) {
  ACTS_VERBOSE("Update single surface status for surface: "
               << surface.geometryId() << " index " << static_cast<int>(index));
  static const std::array<std::size_t,23> geoID = {
     1585267756029461248,
        1585267756029461504,
        1585267618590573312,
        1585267618590573568,
        1585267481151619840,
        1585267481151620096,
        1585267343712715520,
        1585267343712715776,
        1585267206273582336,
        1585267206273582592,
        576464325716219904,
        576464188277266432,
        576464050838308352,
        576463913399354880,
        576464325716220160,
        576464188277266688,
        576464050838308608,
        576463913399355136,
        936750646638415872,
        1008807553481577984,
        1080865010080549632,
        1080865010080549888,
        648518483780306176
  };

  auto sIntersection =
      surface.intersect(state.geoContext, stepper.position(state),
                        navDir * stepper.direction(state), boundaryTolerance,
                        surfaceTolerance)[index];
  if (std::find(geoID.begin(),geoID.end(), surface.geometryId().value())!=geoID.end()) {
     std::size_t geo_id = surface.geometryId().value();
     (void) geo_id;
     std::cout <<"DEBUG SteppingHelper " << __LINE__ << " geo " << surface.geometryId() << " [" << surface.geometryId().value() << "]" << std::endl;
  }

  std::stringstream out;
  out << "DEBUG updateSingleSurfaceStatus "  << surface.geometryId() << " [" << surface.geometryId().value() << "] ";
  // The intersection is on surface already
  if (sIntersection.status() == Intersection3D::Status::onSurface) {
    // Release navigation step size
    state.stepSize.release(ConstrainedStep::actor);
    ACTS_VERBOSE("Intersection: state is ON SURFACE");
    out << " on-surface" << std::endl;
    std::cout << out.str() << std::flush;
    return Intersection3D::Status::onSurface;
  }

  const double nearLimit = std::numeric_limits<double>::lowest();
  const double farLimit = state.stepSize.value(ConstrainedStep::aborter);

  if (sIntersection.isValid() &&
      detail::checkPathLength(sIntersection.pathLength(), nearLimit, farLimit,
                              logger)) {
    ACTS_VERBOSE("Surface is reachable");
    stepper.updateStepSize(state, sIntersection.pathLength(),
                           ConstrainedStep::actor);
    out << " reachabele " << nearLimit << " " << sIntersection.pathLength() <<  " < " << farLimit << std::endl;
    std::cout << out.str() << std::flush;
    return Intersection3D::Status::reachable;
  }

  ACTS_VERBOSE("Surface is NOT reachable");
  out << " not reachable" << std::endl;
  std::cout << out.str() << std::flush;
  return Intersection3D::Status::unreachable;
}

/// Update the Step size - single component
///
/// It takes a (valid) object intersection from the compatibleX(...)
/// calls in the geometry and updates the step size
///
/// @param state [in,out] The stepping state (thread-local cache)
/// @param oIntersection [in] The object that yielded this step size
/// @param release [in] A release flag
template <typename stepper_t, typename object_intersection_t>
void updateSingleStepSize(typename stepper_t::State& state,
                          const object_intersection_t& oIntersection,
                          bool release = true) {
  double stepSize = oIntersection.pathLength();
  state.stepSize.update(stepSize, ConstrainedStep::actor, release);
}

}  // namespace Acts::detail
