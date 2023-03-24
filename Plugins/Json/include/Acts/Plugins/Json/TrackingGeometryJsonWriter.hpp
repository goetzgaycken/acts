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
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <memory>

namespace Acts {
/// Helper to read a tracking geometry from a json object or file.
class TrackingGeometryJsonWriter
{
public:

   struct Config : public Acts::MaterialMapJsonConverter::Config {
      // /// Name of the detector being build
      // std::string detectorName = "";
      // /// Logging level of the child tools
      // Acts::Logging::Level toolLogLevel = Acts::Logging::INFO;
      // /// Logging level
      Acts::Logging::Level logLevel = Acts::Logging::INFO;
   };

   /// Constructor of the tracking geometry json reader.
   /// @param cfg configuration object to set logging levels for the builders and the detector name.
   /// @param logger an optional external logger to replace the default logger.
   TrackingGeometryJsonWriter(
       const Config& cfg,
       std::unique_ptr<const Acts::Logger> logger = std::unique_ptr<const Acts::Logger>{}) :
          m_cfg(cfg),
          m_logger(std::move(logger))
   {
      if (!m_logger) {
         m_logger = std::move( Acts::getDefaultLogger("TrackingGeometryJsonWriter", m_cfg.logLevel) );
      }
   }

   /// Read a tracking geometry described by the given json file.
   /// @param file_name name of the json file which describes the tracking geometry
   /// @param mdecorator an invalid pointer or a decorator to override the material description inside the json file
   /// @return the tracking geometry.
   void write(
       const std::shared_ptr<const Acts::TrackingGeometry> &tracking_geometry,
       const std::string &file_name,
       const Acts::GeometryContext& gctx) const;

   void writeNominal(
       const std::shared_ptr<const Acts::TrackingGeometry> &tracking_geometry,
       const std::string &file_name) const;

   /// Create a json object which describes the tracking geometry
   /// @return json object describing the tracking geometry
   const nlohmann::json trackingGeometryToJson(const std::shared_ptr<const Acts::TrackingGeometry> &tracking_geometry,
                                               const Acts::GeometryContext& gctx=Acts::GeometryContext()) const;

protected:
   /// Private access to the logger
   const Acts::Logger& logger() const { return *m_logger; }

   Config m_cfg;

   /// Logging instance
   std::unique_ptr<const Acts::Logger> m_logger;

};
}

