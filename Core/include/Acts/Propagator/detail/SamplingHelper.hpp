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
#include <array>

#include <limits>

#include <variant>
#include <iostream>
#include <type_traits>


namespace Acts::SamplingHelper {
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
   
   template <typename Object>
   using options_type = decltype(std::declval<Object>().options);

   template <typename  Object, typename = std::void_t<> >
   struct has_options : std::false_type{};

   template <typename Object>
   struct has_options<Object, std::void_t<options_type<Object> > > : std::true_type{};




   template <typename Object>
   using extension_type = decltype(std::declval<Object>().extension);

   template <typename  Object, typename = std::void_t<> >
   struct has_extension : std::false_type{};

   template <typename Object>
   struct has_extension<Object, std::void_t<extension_type<Object> > > : std::true_type{};
   

   
   template <typename stepper_state_t>
   unsigned int getNSamplers(const stepper_state_t &state) {
      if constexpr(has_extension<stepper_state_t>::value) {
         return state.extension.n_samples;
      }
      else {
         return 0u;
      }
   }

   

   template <std::size_t N_SAMPLES>
   inline
   constexpr auto makeOffsets()  -> std::array<std::pair<double,double>,N_SAMPLES-1 >
      requires(N_SAMPLES>1)
   {
      constexpr double phase0=M_PI/4.;
      constexpr double phase_step = 2*M_PI/(N_SAMPLES-1);
      std::array<std::pair<double,double>,N_SAMPLES-1 >  step_12;
      for (unsigned int i=0; i<N_SAMPLES-1; ++i) {
         step_12[i]=std::make_pair(cos(phase0 + phase_step*i), sin(phase0 + phase_step*i ));
      }
      return step_12;
   }

   template <std::size_t N_SAMPLES>
   inline
   std::array<std::pair<Vector3,Vector3>, N_SAMPLES> makeTrajectorySamples(const Acts::Vector3       &position,
                                                                           const Acts::Vector3       &dir,
                                                                           const BoundSquareMatrix   &curvi_cov,
                                                                           Acts::Direction           &step_direction)
      requires( N_SAMPLES>0)
   {
      std::array<std::pair<Vector3,Vector3>, N_SAMPLES> trajectory_samples;
      trajectory_samples[0] = std::make_pair(position,
                                             step_direction * dir);
      if constexpr(N_SAMPLES>1) {
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

         static constexpr std::array<std::pair<double,double>,N_SAMPLES-1 >  step_12=makeOffsets<N_SAMPLES>();
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
                                                      step_direction * (trajectory_samples[0].second
                                                                        + delta_dirphi   * step_12[i].first
                                                                        + delta_dirtheta * step_12[i].second));
         }
      }
      return trajectory_samples;
   }

   template <unsigned int N_SAMPLES>
   struct NSamples {
      static constexpr unsigned int getSamples() { return n_samples; }
      static constexpr unsigned int n_samples = N_SAMPLES;
   };

   using AllowedSamples = std::variant<NSamples<1>, NSamples<4>, NSamples<8> > ;

   template <std::size_t N = std::variant_size_v<AllowedSamples> > 
   AllowedSamples makeSampleValue(unsigned int value) {
      if constexpr(N == 0) {
         throw std::runtime_error("Not allowed");
         return decltype(std::get<0>(AllowedSamples{})){};
      }
      else {
         if (value == std::remove_reference<decltype(std::get<N-1>(AllowedSamples{}))>::type ::n_samples) {
            return decltype(std::get<N-1>(AllowedSamples{})){};
         }
         return makeSampleValue<N-1>(value);
      }
   }
}
