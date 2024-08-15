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
#include "Acts/Geometry/GeoIdRegistry.hpp"

#include <sstream>
#include <iostream>
#include <span>

#include <limits>

#include "SamplingHelper.hpp"

namespace Acts::detail {

namespace {

  Acts::SurfaceMultiIntersection intersect(const GeometryContext &gctx,
                                           const Surface& surface,
                                           std::uint8_t index,
                                           std::span<std::pair<Vector3,Vector3> > trajectory_samples,
                                           const BoundaryTolerance& boundaryTolerance,
                                           ActsScalar surfaceTolerance,
                                           double nearLimit,
                                           double farLimit,
                                           const Logger& logger) {

      std::stringstream msg;
      unsigned int n_max_samples = surface.associatedDetectorElement() != nullptr ? trajectory_samples.size() : 1;

      Dbg::GeoIdHelper  geo_ids;
      bool is_known = geo_ids.isKnown(surface.geometryId().value());
      
      for (unsigned int sample_i=0; sample_i<n_max_samples; ++sample_i) {

         auto sMultiIntersection = surface.intersect(gctx, trajectory_samples[sample_i].first,
                                                     trajectory_samples[sample_i].second,
                                                     boundaryTolerance,
                                                     surfaceTolerance);
         msg << "DEBUG updateSingleSurfaceStatus "  << surface.geometryId()
             << " [" << std::hex << surface.geometryId().value() << std::dec << (is_known ? "*" : "") << "] "
             << " pos: #" << sample_i << "/" << n_max_samples << " " << trajectory_samples[sample_i].first[0]
             << " " << trajectory_samples[sample_i].first[1] << " " << trajectory_samples[sample_i].first[2]
             << (sMultiIntersection.size() > index ? "" : " index-out-of-range")
             << " status : " << (sMultiIntersection.size() > index ? static_cast<int>( sMultiIntersection[index].status()) : -1)
             << " " << nearLimit << " < " << (sMultiIntersection.size() > index
                                              ?  sMultiIntersection[index].pathLength()
                                              : std::numeric_limits<double>::max())
             << " < " << farLimit 
             << std::endl;
         if (sMultiIntersection.size() > index) {

            if ((sMultiIntersection[index].status() == Intersection3D::Status::onSurface)
                || (sMultiIntersection[index].isValid() &&
                    detail::checkPathLength(sMultiIntersection[index].pathLength(), nearLimit, farLimit,
                                            logger))
                || sample_i+1==n_max_samples) {
               std::cout << msg.str() << std::flush;
               return sMultiIntersection;
            }
         }
      }
      std::cout << msg.str() << std::flush;
      throw std::runtime_error("No multi intersection with enough elements");
      
   }


template <typename stepper_t, typename state_t>
inline Acts::SurfaceMultiIntersection intersectionHelper (const stepper_t &stepper,
                                                   const state_t &state_stepping,
                                                   const Surface& surface,
                                                   std::uint8_t index,
                                                   Direction navDir,
                                                   const BoundaryTolerance& boundaryTolerance,
                                                   ActsScalar surfaceTolerance,
                                                   double nearLimit,
                                                   double farLimit,
                                                   const Logger& logger) {
   
   auto principalIntersection  = surface.intersect(state_stepping.geoContext, stepper.position(state_stepping),
                                                   navDir * stepper.direction(state_stepping),
                                                   (surface.associatedDetectorElement() != nullptr
                                                   ? BoundaryTolerance::Infinite()
                                                   : boundaryTolerance),
                                                   surfaceTolerance);
   //if constexpr(has_stepping<decltype(state)>::value) {
   //          state.gwer();
   if constexpr(Acts::SamplingHelper::has_cov<decltype(state_stepping)>::value)  {
       if ( Acts::SamplingHelper::getNSamplers(state_stepping)>0) {

      auto samplingIntersection =  std::visit( [&stepper,
                          &state_stepping,
                          navDir,
                          &surface,
                          index,
                          &boundaryTolerance,
                          &surfaceTolerance,
                          &nearLimit,
                          &farLimit,
                          &logger
                          ](const auto &elm) {
         constexpr unsigned int N_SAMPLES = elm.getSamples();
            std::array<std::pair<Vector3,Vector3>, N_SAMPLES>
               trajectory_samples = Acts::SamplingHelper::makeTrajectorySamples<N_SAMPLES>(stepper.position(state_stepping),
                                                                                           stepper.direction(state_stepping),
                                                                                           state_stepping.cov,
                                                                                           navDir,
                                                                                           Acts::SamplingHelper::getNSigmas(state_stepping));
            return intersect(state_stepping.geoContext,
                             surface,
                             index,
                             std::span(trajectory_samples.begin(),trajectory_samples.end()),
                             boundaryTolerance,
                             surfaceTolerance,
                             nearLimit,
                             farLimit,
                             logger);

      }, Acts::SamplingHelper::makeSampleValue(Acts::SamplingHelper::getNSamplers(state_stepping)) );
      if (principalIntersection[index].isValid() != samplingIntersection[index].isValid()
          || detail::checkPathLength(principalIntersection[index].pathLength(), nearLimit, farLimit,
                                     logger)
             != detail::checkPathLength(samplingIntersection[index].pathLength(), nearLimit, farLimit,
                                        logger)
          || principalIntersection[index].status() != samplingIntersection[index].status()) {
         Vector3 pos = stepper.position(state_stepping);

         Acts::Vector3 loc3Dframe = (surface.transform(state_stepping.geoContext).inverse()) * pos;
         Dbg::GeoIdHelper  geo_ids;
         bool is_known = geo_ids.isKnown(surface.geometryId().value());
         
         std::stringstream out;
         out << "DEBUG updateSingleSurfaceStatus sampling vs principal"  << surface.geometryId()
             << " [" << std::hex << surface.geometryId().value() << std::dec << (is_known ? "*" : "") << "] "
             << " " << pos[0]
             << " " << pos[1] << " " << pos[2]
             << " (d: " << loc3Dframe.z() << ")"
             << (samplingIntersection.size() > index ? "" : " index-out-of-range")
             << (principalIntersection.size() > index ? "" : " index-out-of-range")
             << " status :"
             << " " << (samplingIntersection.size() > index ? static_cast<int>( samplingIntersection[index].status()) : -1)
             << " " << (principalIntersection.size() > index ? static_cast<int>( principalIntersection[index].status()) : -1)
             << " " << nearLimit << " < "
             << (samplingIntersection.size() > index ?  samplingIntersection[index].pathLength() : std::numeric_limits<double>::max())
             << ", " << (principalIntersection.size() > index ?  principalIntersection[index].pathLength() : std::numeric_limits<double>::max())
             << " < " << farLimit 
             << std::endl;
         std::cout << out.str() << std::flush;
      }
   }
      // else {
      //    return surface.intersect(state_stepping.geoContext, stepper.position(state_stepping),
      //                             navDir * stepper.direction(state_stepping), boundaryTolerance,
      //                             surfaceTolerance);
      // }
   }
   return principalIntersection;

}
}
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

  Dbg::GeoIdHelper  geo_ids;

  const double nearLimit = std::numeric_limits<double>::lowest();
  const double farLimit = state.stepSize.value(ConstrainedStep::aborter);

auto sMultiIntersection  
   = intersectionHelper(stepper,
                                          state,
                                          surface,
                                          index,
                                          navDir,
                                          boundaryTolerance,
                                          surfaceTolerance,
                                          nearLimit,
                                          farLimit,
                        logger);
 
 if ( index>=sMultiIntersection.size() ) {
      throw std::runtime_error("No multi intersection with enough elements");
 }
 auto sIntersection = sMultiIntersection[index];
 if (sIntersection.status() == Intersection3D::Status::onSurface
     && sIntersection.pathLength()>0.001) {
    throw std::runtime_error("Not on surface.");
 }
      // surface.intersect(state.geoContext, stepper.position(state),
      //                   navDir * stepper.direction(state), boundaryTolerance,
      //                   surfaceTolerance)[index];

  bool is_known = geo_ids.isKnown(surface.geometryId().value());
  if (is_known) {

     std::size_t geo_id = surface.geometryId().value();
     (void) geo_id;
     std::cout <<"DEBUG SteppingHelper " << __LINE__ << " geo " << surface.geometryId() << " [" << std::hex << surface.geometryId().value() << std::dec << "*]" << std::endl;
  }

  unsigned int n_samples =0;
  if constexpr(Acts::SamplingHelper::has_cov<decltype(state)>::value)  {
     n_samples = Acts::SamplingHelper::getNSamplers(state);
  }

  std::stringstream out;
  Vector3 pos = stepper.position(state);
  out << "DEBUG updateSingleSurfaceStatus "  << surface.geometryId() << " ["
      << std::hex << surface.geometryId().value() << std::dec << (is_known ? "*" : "") << "] "
      << " pos: " << pos[0] << " " << pos[1] << " " << pos[2];
  // The intersection is on surface already
  if (sIntersection.status() == Intersection3D::Status::onSurface) {
    // Release navigation step size

  if  (surface.associatedDetectorElement() != nullptr) {
     Acts::Vector3 position = stepper.position(state);
     Acts::Vector3 loc3Dframe = (surface.transform(state.geoContext).inverse()) * position;
     if (std::abs(loc3Dframe.z())< s_onSurfaceTolerance) {
        state.stepSize.release(ConstrainedStep::actor);
        ACTS_VERBOSE("Intersection: state is ON SURFACE");
        out << " on-surface ("  << loc3Dframe.z() << "; index " << static_cast<unsigned int>(index) << "; pathLength  " <<  sIntersection.pathLength() << " )" << std::endl;
        std::cout << out.str() << std::flush;
        return Intersection3D::Status::onSurface;
     }
     else {
        if (detail::checkPathLength(loc3Dframe.z(), nearLimit, farLimit,
                                    logger)) {
           ACTS_VERBOSE("Surface is reachable");
           stepper.updateStepSize(state, sIntersection.pathLength(),
                                  ConstrainedStep::actor);
           out << " reachabele " << nearLimit << " " << loc3Dframe.z() << "; " << sIntersection.pathLength()<<  " < " << farLimit << std::endl;
           std::cout << out.str() << std::flush;
           return Intersection3D::Status::reachable;
        }
        else {
           ACTS_VERBOSE("Surface is NOT reachable");
           out << " not reachable"
               <<  (sIntersection.isValid()  ? " intersecting" : " not-intersecting" )
               << nearLimit << " " << loc3Dframe.z() << "; " <<  sIntersection.pathLength() <<  " < " << farLimit
               << " (samples " << n_samples << ")" 
               << std::endl;
           std::cout << out.str() << std::flush;
           return Intersection3D::Status::unreachable;
        }
     }
  }
  else {
        state.stepSize.release(ConstrainedStep::actor);
        ACTS_VERBOSE("Intersection: state is ON SURFACE");
        out << " on-surface (index " << static_cast<unsigned int>(index) << "; pathLength  " <<  sIntersection.pathLength() << " )" << std::endl;
        std::cout << out.str() << std::flush;
        return Intersection3D::Status::onSurface;
  }
  }
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
  out << " not reachable"
      <<  (sIntersection.isValid()  ? " intersecting" : " not-intersecting" )
      << nearLimit << " " << sIntersection.pathLength() <<  " < " << farLimit
      << " (samples " << n_samples << ")" 
      << std::endl;
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
