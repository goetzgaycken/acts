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

      std::string geantinoInputFileName;
      uint64_t maxGeantinoEntries;
      bool stop = false;
   };
   TrackingGeometryJsonReader(const Config& cfg)  : m_cfg(cfg) {}

   std::shared_ptr<const Acts::TrackingGeometry> read(const std::string &file_name,
                                                      std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) const;
   std::shared_ptr<const Acts::TrackingGeometry> createTrackingGeometry(const nlohmann::json& tracking_geometry_description,
                                                                        std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) const;

protected:
   Config m_cfg;

   Acts::ProtoDetector createProtoDetector(const nlohmann::json& tracking_geometry_description,
                                           const std::vector<std::shared_ptr<Acts::Surface>> &surfaces,
                                           const std::vector< std::vector< unsigned int >> &surfaces_per_volume) const;
   std::tuple<std::vector<std::shared_ptr<Acts::Surface>>,
              std::vector< std::vector< unsigned int >>,
              std::vector<std::unique_ptr<DetectorElementBase>> >
      createSurfaces(const nlohmann::json& tracking_geometry_description) const;

   template <typename enum_t, class array_t>
   static enum_t getNamedEnum(const array_t &enum_names, const std::string &the_name, const std::string &error_head="Unknwon enum: ") {
      assert( enum_names.size() < std::numeric_limits<unsigned short>::max());
      typename array_t::const_iterator
         enum_iter = std::find(enum_names.begin(), enum_names.end(), the_name);
      if ( enum_iter == enum_names.end()) {
         throw std::runtime_error(error_head+the_name);
      }
      return static_cast<enum_t>(enum_iter - enum_names.begin());
   }

   static float scale(float a, float b) {
      float tmp=std::abs(a+b);
      return tmp == 0. ? 2. : tmp;
   }
   static bool floatsAgree(float a, float b) {
      return std::abs(a-b) < std::numeric_limits<float>::epsilon() * scale(a,b);
   }

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

      // @TODO remove :
      enum EBinningVariable {kR, kPhi, kZ, kNBinningVariables};
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
         static float scale(float a, float b) {
            float tmp=std::abs(a+b);
            return tmp == 0. ? 2. : tmp;
         }
         bool overlaps(const Binning &a_binning) const {
            return min < a_binning.max && max > a_binning.min;
         }
      };
      std::array<Binning, kNBinningVariables> binning {};
      //  remove until here

      unsigned short setDomainMask = 0;
      unsigned short noOverlapMask = 0;

      bool processed = false;
   };

   static VolumeInfo &registerVolume(unsigned int volume_id,
                                     const std::string &name,
                                     std::unordered_map<unsigned int, VolumeInfo > &volume_name_map,
                                     std::unordered_map<std::string, unsigned int > &name_map);

};
}

