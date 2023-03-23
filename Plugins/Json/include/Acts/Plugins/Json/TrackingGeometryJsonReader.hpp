// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"

#include <memory>
#include <array>
#include <limits>
#include <stdexcept>

namespace Acts {
   class IMaterialDecorator;

/// Helper to read a tracking geometry from a json object or file.
class TrackingGeometryJsonReader
{
public:

   struct Config {
      /// Name of the detector being build
      std::string detectorName = "";
      /// Logging level of the child tools
      Acts::Logging::Level toolLogLevel = Acts::Logging::INFO;
      /// Logging level
      Acts::Logging::Level logLevel = Acts::Logging::INFO;
   };

   /// Constructor of the tracking geometry json reader.
   /// @param cfg configuration object to set logging levels for the builders and the detector name.
   /// @param logger an optional external logger to replace the default logger.
   TrackingGeometryJsonReader(
       const Config& cfg,
       std::unique_ptr<const Acts::Logger> logger = std::unique_ptr<const Acts::Logger>{}) :
          m_cfg(cfg),
          m_logger(std::move(logger))
   {
      if (!m_logger) {
         m_logger = std::move( Acts::getDefaultLogger("TrackingGeometryJsonReader", m_cfg.logLevel) );
      }
   }

   /// Read a tracking geometry described by the given json file.
   /// @param file_name name of the json file which describes the tracking geometry
   /// @param mdecorator an invalid pointer or a decorator to override the material description inside the json file
   /// @return the tracking geometry.
   std::shared_ptr<const Acts::TrackingGeometry> read(
       const std::string &file_name,
       std::shared_ptr<const Acts::IMaterialDecorator> mdecorator = std::shared_ptr<const Acts::IMaterialDecorator>()) const;

   /// Create a tracking geometry described by the given json object.
   /// @param tracking_geometry_description the json object which describes the tracking geometry
   /// @param mdecorator an invalid pointer or a decorator to override the material description inside the json file
   /// @return the tracking geometry.
   std::shared_ptr<const Acts::TrackingGeometry> createTrackingGeometry(const nlohmann::json& tracking_geometry_description,
                                                                        std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) const;

protected:
   /// Private access to the logger
   const Acts::Logger& logger() const { return *m_logger; }

   /// Create a proto detector using the list of volumes described by the json object
   /// @param tracking_geometry_description the json object created by e.g. the tracking geometry writer
   /// @return the proto detector
   Acts::ProtoDetector createProtoDetector(const nlohmann::json& tracking_geometry_description) const;

   /// Create surfaces using the list of surfaces described by the json object
   /// @param tracking_geometry_description the json object created by e.g. the tracking geometry writer
   /// @param ignore_material if true the material description in the json object is ignored
   /// @return surface list, and list of associated detector elements
   std::tuple<std::vector<std::shared_ptr<Acts::Surface>>,
              std::vector<std::unique_ptr<DetectorElementBase>> >
      createSurfaces(const nlohmann::json& tracking_geometry_description,
                     bool ignore_material) const;

   Config m_cfg;

   /// Logging instance
   std::unique_ptr<const Acts::Logger> m_logger;

   static float scale(float a, float b) {
      float tmp=std::abs(a+b);
      return tmp == 0. ? 2. : tmp;
   }
   /// Test whether two floats agree to very high precision
   static bool floatsAgree(float a, float b) {
      return std::abs(a-b) < std::numeric_limits<float>::epsilon() * scale(a,b);
   }

   /// Helper class which describes a volume
   struct VolumeInfo {
      VolumeInfo(const std::string a_name, unsigned int an_id=std::numeric_limits<unsigned int>::max())
         : id(an_id) {proto_volume.name=a_name; }
      /// the name given to the volume
      Acts::ProtoVolume proto_volume;
      /// the external volume id
      unsigned int id;
      /// child indexes or empty if no child volumes
      std::vector<unsigned int> childs;
      /// parent index or UINT_MAX
      unsigned int parent = std::numeric_limits<unsigned int>::max();

      std::unique_ptr<Acts::VolumeBounds> volumeBounds;
      Acts::Transform3                    volumeTransform;

      struct Domain {
         enum EType {kOpen, kClosed, kNRangeTypes};
         Domain() = default;
         Domain(float a_min, float a_max)
            : min(std::min(a_min,a_max)), max(std::max(a_min,a_max)) {}
         float min = -std::numeric_limits<float>::max();
         float max = std::numeric_limits<float>::max();

         bool operator==(const Domain &a_domain) const {
            return floatsAgree(max,a_domain.max) && floatsAgree(min,a_domain.min);
         }
         bool operator!=(const Domain &a_domain) const {
            return not operator==(a_domain);
         }
         bool overlaps(const Domain &a_domain) const {
            return min < a_domain.max && max > a_domain.min;
         }
      };

      enum EDomainTypes { eR, eZ, eNDomains};
      std::array<Domain, eNDomains> domains;

      /// helper class to hold binning description
      struct Binning {
         enum EType {kOpen, kClosed, kNRangeTypes};
         Binning() = default;
         Binning(float a_min, float a_max, EType a_type)
            : min(std::min(a_min,a_max)), max(std::max(a_min,a_max)),type(a_type) {}
         float min = -std::numeric_limits<float>::max();
         float max = std::numeric_limits<float>::max();
         EType type = kNRangeTypes;
         bool isClosed() const {
            float a_min = min;
            float a_max = max;
            if (a_min<=-M_PI) a_min+=2*M_PI;
            if (a_max>M_PI)   a_max-=2*M_PI;
            return  (a_max-a_min>M_PI && std::abs(a_max-a_min)< 2*std::numeric_limits<float>::epsilon());
         }
         bool operator==(const Binning &a_binning) const {
            return type == a_binning.type && floatsAgree(max,a_binning.max) && floatsAgree(min,a_binning.min);
         }
         bool operator!=(const Binning &a_binning) const {
            return not operator==(a_binning);
         }
         bool overlaps(const Binning &a_binning) const {
            return min < a_binning.max && max > a_binning.min;
         }
      };

      unsigned short setDomainMask = 0; // specifies for a volume which of the domains (R,Z) are specified.
      unsigned short noOverlapMask = 0; // specifies in which domain the child volumes do not overlap.

      bool processed = false;
   };

   /// helper method to register a new volume and link its child volumes.
   /// @param volume_id the volume id
   /// @param name the volume name
   /// @param volume_name_map associative container for the volume information where the key is the volume id
   /// @param name_map inverse map which maps a volume name to its id.
   static VolumeInfo &registerVolume(unsigned int volume_id,
                                     const std::string &name,
                                     std::unordered_map<unsigned int, VolumeInfo > &volume_name_map,
                                     std::unordered_map<std::string, unsigned int > &name_map);

};
}

