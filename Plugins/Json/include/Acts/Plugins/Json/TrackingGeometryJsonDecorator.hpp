// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Plugins/Json/ITrackingGeometryJsonDecorator.hpp"
#include "Acts/Geometry/IdentifiableDetectorElement.hpp"
#include "Acts/Plugins/Json/VolumeBoundsJsonConverter.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"
#include "Acts/Plugins/Json/MaterialJsonConverter.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include <cstring>

namespace  Acts {

  class TrackingGeometryJsonDecorator : public Acts::ITrackingGeometryJsonDecorator {
  public:
     TrackingGeometryJsonDecorator(const Acts::GeometryContext &geometry_context) : m_geometryContext(&geometry_context) {}

    /// @brief Add extra elements to the json object already filled for the given
    /// surface
    ///
    /// @param surface the surface which was used to fill the json object
    /// @pram json the json object that is enhanced
    /// Will add:
    ///   - detector Ids for all surfaces which have an associated
    ///     detector element and which provide Ids.
    ///   - will convert the surface material to json.
    virtual void decorate(
      const Acts::Surface &surface,
      [[maybe_unused]] nlohmann::json &json) const override {
       if (surface.associatedDetectorElement()) {
           const Acts::IdentifiableDetectorElement *
              det_element = dynamic_cast<const Acts::IdentifiableDetectorElement *>(surface.associatedDetectorElement());
          if (det_element) {
            json["detID"]=det_element->detectorId();
          }
       }
       const Acts::ISurfaceMaterial *surface_material = surface.surfaceMaterialSharedPtr().get();
       if (surface_material) {
          Acts::to_json(json["material"],surface_material);
       }
    }

    /// @brief Add extra elements to the json object already filled for the given
    /// volume
    ///
    /// @param volume the volume which was used to fill the json object
    /// @param json the json object that is enhanced
    /// will add :
    ///   - volume bounds and transform
    ///   - layer information : properties of the surface representation, thickness
    virtual void decorate(
      const Acts::TrackingVolume &volume,
      nlohmann::json &json) const override {
       Acts::to_json(json["bounds"],volume.volumeBounds());
       Acts::to_json(json["transform"],volume.transform());
       if (volume.confinedLayers()) {
          Acts::Surface::SurfaceType surface_type_all = Acts::Surface::Other;
          bool mixed=false;
          bool first=true;
          bool has_layers=false;

          const std::vector<Acts::LayerPtr> &confined_layers = volume.confinedLayers()->arrayObjects();
          for (const Acts::LayerPtr &layer_ptr : confined_layers) {
             if (layer_ptr && layer_ptr->thickness()>0.) {
                has_layers=true;
             }
             const Acts::Surface& layer_surface = layer_ptr->surfaceRepresentation();
             Acts::Surface::SurfaceType surface_type = layer_surface.type();
             if (first) {
                surface_type_all = surface_type;
                first=false;
             }
             else if (surface_type != surface_type_all) {
                mixed=true;
             }
             if (mixed && has_layers) break;
          }

          if (has_layers) {
             nlohmann::json layer_container = nlohmann::json::array();

             //layer_container.reserve( confined_layers.size());
             for (const Acts::LayerPtr &layer_ptr : confined_layers) {
                if (layer_ptr && layer_ptr->thickness()>0.) {
                   has_layers=true;
                   const Acts::Surface& layer_surface = layer_ptr->surfaceRepresentation();
                   nlohmann::json layer_entry;
                   layer_entry["layer"] = layer_ptr->geometryId().layer();
                   Acts::toJson(layer_entry, layer_surface, *m_geometryContext);
                   layer_entry["thickness"]=layer_ptr->thickness();
                   layer_container.push_back(std::move(layer_entry) );
                }
             }
             json["layers"] = std::move(layer_container);
          }
          json["layerType"]=mixed ? std::string("mixed") : removeSurfaceSuffix(Acts::surfaceTypes[surface_type_all]);
       }
    }
   protected:
     const Acts::GeometryContext *m_geometryContext = nullptr;
      // copied from Acts/Plugins/Json/GeometryHierarchyMapJsonConverter.hpp
      static nlohmann::json encodeIdentifier(Acts::GeometryIdentifier id) {
         nlohmann::json encoded;
         // only store non-zero identifiers
         if (id.volume() != 0u) {
            encoded["volume"] = id.volume();
         }
         if (id.boundary() != 0u) {
            encoded["boundary"] = id.boundary();
         }
         if (id.layer() != 0u) {
            encoded["layer"] = id.layer();
         }
         if (id.approach() != 0u) {
            encoded["approach"] = id.approach();
         }
         if (id.sensitive() != 0u) {
            encoded["sensitive"] = id.sensitive();
         }
         return encoded;
      }

      static std::string removeSurfaceSuffix(std::string in) {
         constexpr std::string::size_type str_length = strlen("Surface");
         if  (in.size() > str_length && strncmp(in.data()+(in.size() - str_length), "Surface", str_length)==0) {
            return in.substr(0,in.size()-str_length);
         }
         else {
            return in;
         }
      }
   };
}
