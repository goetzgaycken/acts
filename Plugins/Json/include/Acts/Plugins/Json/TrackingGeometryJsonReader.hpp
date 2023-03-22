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
   TrackingGeometryJsonReader(
       const Config& cfg,
       std::unique_ptr<const Acts::Logger> logger = std::unique_ptr<const Acts::Logger>{}) : m_cfg(cfg) {
       if (!logger) {
          logger = Acts::getDefaultLogger("TrackingGeometryJsonReader", m_cfg.logLevel);
       }
   }

   std::shared_ptr<const Acts::TrackingGeometry> read(const std::string &file_name,
                                                      std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) const;
   std::shared_ptr<const Acts::TrackingGeometry> createTrackingGeometry(const nlohmann::json& tracking_geometry_description,
                                                                        std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) const;

protected:
   /// Private access to the logger
   const Acts::Logger& logger() const { return *m_logger; }

   /// Create a proto detector using the list of volumes described by the json object
   /// @param tracking_geometry_description the json object created by e.g. the tracking geometry writer
   /// @param surfaces the list of all surfaces described by the given json
   /// @param surfaces_per_volume surface indices per volume
   /// @return the proto detector
   Acts::ProtoDetector createProtoDetector(const nlohmann::json& tracking_geometry_description,
                                           const std::vector<std::shared_ptr<Acts::Surface>> &surfaces,
                                           const std::vector< std::vector< unsigned int >> &surfaces_per_volume) const;

   /// Create surfaces using the list of surfaces described by the json object
   /// @param tracking_geometry_description the json object created by e.g. the tracking geometry writer
   /// @return surface list, surfaces indices per volume, and list of associated detector elements
   std::tuple<std::vector<std::shared_ptr<Acts::Surface>>,
              std::vector< std::vector< unsigned int >>,
              std::vector<std::unique_ptr<DetectorElementBase>> >
      createSurfaces(const nlohmann::json& tracking_geometry_description) const;

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
      /// child indices or empty if no child volumes
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

