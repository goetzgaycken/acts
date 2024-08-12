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

   namespace {
   template <typename Object>
   using cov_type = decltype(std::declval<Object>().cov);

   template <typename  Object, typename = std::void_t<> >
   struct has_cov : std::false_type{};

   template <typename Object>
   struct has_cov<Object, std::void_t<cov_type<Object> > > : std::true_type{};
      
   template <typename Object>
   using stepping_type = decltype(std::declval<Object>().stepping);

   template <typename  Object, typename = std::void_t<> >
   struct has_stepping : std::false_type{};

   template <typename Object>
   struct has_stepping<Object, std::void_t<stepping_type<Object> > > : std::true_type{};

   std::array<std::pair<Vector3,Vector3>, 4> makeTrajectorySamples(const Acts::Vector3       &position,
                                                                   const Acts::Vector3       &dir,
                                                                   const BoundSquareMatrix   &curvi_cov,
                                                                   Acts::Direction           &navDir)
   {
      std::array<std::pair<Acts::Vector3,Acts::Vector3>, 4> trajectory_samples;
      trajectory_samples[0] = std::make_pair(position,
                                             dir);
    auto delta_phi = std::sqrt(curvi_cov(eBoundPhi,eBoundPhi));
    auto delta_theta = std::sqrt(curvi_cov(eBoundTheta,eBoundTheta));
    // @TODO switch axis for large eta ...
    Vector3 dir_u = Vector3::UnitZ().cross(trajectory_samples[0].second).normalized();
    Vector3 dir_v = trajectory_samples[0].second.cross(dir_u);

    constexpr double multiplier = 4;

    Vector3 delta_u = dir_u * std::sqrt(curvi_cov(eBoundLoc0,eBoundLoc0)) * multiplier;
    Vector3 delta_v = dir_v * std::sqrt(curvi_cov(eBoundLoc1,eBoundLoc1)) * multiplier;

    Vector3 delta_dirphi;
    delta_dirphi    << /* -sin(phi) * sin(theta) */ -trajectory_samples[0].second[1] * delta_phi * multiplier,
                       /*  cos(phi) * sin(theta) */  trajectory_samples[0].second[0] * delta_phi * multiplier,
                        0;
    // @TODO use jacobi ?
    double phi = std::atan2(trajectory_samples[0].second[1], trajectory_samples[0].second[0]);
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);
    double sin_theta = sin( std::acos(trajectory_samples[0].second[2]) );
    double norm_dir = trajectory_samples[0].second.norm();
    
    double sin_theta_alt = std::sqrt( std::max(0., 1 - trajectory_samples[0].second[2] *trajectory_samples[0].second[2]));
    double inv_sin_theta = 1./sin_theta_alt;
    double cos_phi_alt = trajectory_samples[0].second[0]*inv_sin_theta;
    double sin_phi_alt = trajectory_samples[0].second[1]*inv_sin_theta;
    
    Vector3 delta_dirtheta;
    delta_dirtheta <<   /*cos(phi) * cos(theta) */ cos_phi_alt* trajectory_samples[0].second[2] * delta_theta * multiplier,
                        /*sin(phi) * cos(theta) */ sin_phi_alt* trajectory_samples[0].second[2] * delta_theta * multiplier,
                        sin_theta_alt * delta_theta * multiplier;

    //    1., 1.    - 1 * dx  1 * dy   45
    //    1.,-1  ->   0 * dx -2 * dy  -45
    //    -1,-1  ->  -2 * dx  0 * dy  -45-90
    //    -1, 1  ->   0 * dx -2 * dy  -45-180

    static constexpr unsigned int N_SAMPLES = trajectory_samples.size()-1;
    static_assert( N_SAMPLES==3);
    static constexpr double phase0=M_PI/4.;
    static constexpr double phase_step = 2*M_PI/N_SAMPLES;
    static constexpr std::array<std::pair<double,double>,N_SAMPLES >  step_12 {
       std::make_pair(cos(phase0 + phase_step*0), sin(phase0 + phase_step*0 )),
       std::make_pair(cos(phase0 + phase_step*1), sin(phase0 + phase_step*1 )),
       std::make_pair(cos(phase0 + phase_step*2), sin(phase0 + phase_step*2 ))
    };
    // choose direction delta_dirphi, delta_dirtheta to maximise || pos + dir - pos0 ||
    double a1=delta_dirphi.dot(delta_u);
    double a2=delta_dirphi.dot(delta_v);
    if (std::abs(a2)>std::abs(a1)) {
       Vector3 tmp=delta_dirphi;
       double b1=delta_dirtheta.dot(delta_u);
       if (b1<0) {
          delta_dirtheta*=-1.;
       }
       
       delta_dirphi=delta_dirtheta;
       delta_dirtheta=tmp;
       if (a2<0) {
          delta_dirtheta*=-1.;
       }
    }
    else {
       if (a1<0) {
          delta_dirphi*=-1.;
       }
       double b2=delta_dirtheta.dot(delta_v);
       if (b2<0) {
          delta_dirtheta*=-1.;
       }
    }
    
    for (unsigned int i=0; i<step_12.size(); ++i) { 
       trajectory_samples[i+1] = std::make_pair( trajectory_samples[0].first
                                                   + delta_u * step_12[i].first
                                                   + delta_v * step_12[i].second,
                                                 navDir * (trajectory_samples[0].second
                                                                   + delta_dirphi   * step_12[i].first
                                                                   + delta_dirtheta * step_12[i].second));
    }
    return trajectory_samples;
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
   //if constexpr(has_stepping<decltype(state)>::value) {
   //          state.gwer();
       if constexpr(has_cov<decltype(state_stepping)>::value) {
          std::array<std::pair<Vector3,Vector3>, 4>
             trajectory_samples = makeTrajectorySamples(stepper.position(state_stepping),
                                                        stepper.direction(state_stepping),
                                                        state_stepping.cov,
                                                        navDir);
          std::stringstream msg;
          unsigned int n_max_samples = surface.associatedDetectorElement() != nullptr ? trajectory_samples.size() : 1;
          for (unsigned int sample_i=0; sample_i<n_max_samples; ++sample_i) {

             auto sMultiIntersection = surface.intersect(state_stepping.geoContext, trajectory_samples[sample_i].first,
                                                         trajectory_samples[sample_i].second,
                                                         boundaryTolerance,
                                                         surfaceTolerance);
             msg << "DEBUG updateSingleSurfaceStatus "  << surface.geometryId() << " [" << surface.geometryId().value() << "] "
                 << " pos: #" << sample_i << " " << trajectory_samples[sample_i].first[0]
                 << " " << trajectory_samples[sample_i].first[1] << " " << trajectory_samples[sample_i].first[2]
                 << (sMultiIntersection.size() > index ? "" : " index-out-of-range")
                 << " status : " << (sMultiIntersection.size() > index ? static_cast<int>( sMultiIntersection[index].status()) : -1)
                 << " " << nearLimit << " < " << (sMultiIntersection.size() > index ?  sMultiIntersection[index].pathLength() : std::numeric_limits<double>::max())
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
       else {
          return surface.intersect(state_stepping.geoContext, stepper.position(state_stepping),
                                   navDir * stepper.direction(state_stepping), boundaryTolerance,
                                   surfaceTolerance);
       }
       // }
       // else  {
       //    return surface.intersect(state.geoContext, stepper.position(state),
       //                             navDir * stepper.direction(state), boundaryTolerance,
       //                             surfaceTolerance);
       // }
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

  const double nearLimit = std::numeric_limits<double>::lowest();
  const double farLimit = state.stepSize.value(ConstrainedStep::aborter);

  auto sIntersection = intersectionHelper(stepper,
                                          state,
                                          surface,
                                          index,
                                          navDir,
                                          boundaryTolerance,
                                          surfaceTolerance,
                                          nearLimit,
                                          farLimit,
                                          logger)[index];
      // surface.intersect(state.geoContext, stepper.position(state),
      //                   navDir * stepper.direction(state), boundaryTolerance,
      //                   surfaceTolerance)[index];
  if (std::find(geoID.begin(),geoID.end(), surface.geometryId().value())!=geoID.end()) {
     std::size_t geo_id = surface.geometryId().value();
     (void) geo_id;
     std::cout <<"DEBUG SteppingHelper " << __LINE__ << " geo " << surface.geometryId() << " [" << surface.geometryId().value() << "]" << std::endl;
  }

  std::stringstream out;
  Vector3 pos = stepper.position(state);
  out << "DEBUG updateSingleSurfaceStatus "  << surface.geometryId() << " [" << surface.geometryId().value() << "] "
      << " pos: " << pos[0] << " " << pos[1] << " " << pos[2];
  // The intersection is on surface already
  if (sIntersection.status() == Intersection3D::Status::onSurface) {
    // Release navigation step size
    state.stepSize.release(ConstrainedStep::actor);
    ACTS_VERBOSE("Intersection: state is ON SURFACE");
    out << " on-surface" << std::endl;
    std::cout << out.str() << std::flush;
    return Intersection3D::Status::onSurface;
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
