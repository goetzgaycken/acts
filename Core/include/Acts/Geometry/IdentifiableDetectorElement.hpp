// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace  Acts {
   /// @brief Detector element with an ID.
   /// Simple detector element which is associated to surfaces and just provides a detector element ID.
   class IdentifiableDetectorElement : public Acts::DetectorElementBase
   {
   public:
      /// constructor
      /// @param the surface which represents the detector element.
      /// @param detector_element_id the ID of the detector element
      /// @param thickness the thickness of the detector element.
      /// Create an identifiable detector element.
      IdentifiableDetectorElement(const Acts::Surface &surface, uint64_t detector_element_id, float thickness)
        : m_surface(&surface),
          m_transform( surface.transform(Acts::GeometryContext()) ),
          m_id(detector_element_id),
          m_thickness(thickness)
      {}

      /// @param gctx The current geometry context object, e.g. alignment
      const Acts::Transform3& transform([[maybe_unused]] const Acts::GeometryContext& gctx) const override
         { return m_transform; }

      /// Return surface representation
      const Acts::Surface& surface() const override
        { return *m_surface; }

      double thickness() const override { return m_thickness; }

      /// The ID of this detector element.
      uint64_t detectorId() const { return m_id; }
   private:
      const Acts::Surface *m_surface;
      Acts::Transform3 m_transform;
      uint64_t m_id;
      float m_thickness;
   };
}

