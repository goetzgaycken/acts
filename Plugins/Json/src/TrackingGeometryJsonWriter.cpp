#include "Acts/Plugins/Json/TrackingGeometryJsonWriter.hpp"
#include "Acts/Plugins/Json/TrackingGeometryJsonDecorator.hpp"
#include <fstream>
#include <signal.h>
#include <unistd.h>

namespace Acts {
   void TrackingGeometryJsonWriter::write(
       const std::shared_ptr<const Acts::TrackingGeometry> &tracking_geometry,
       const std::string &file_name,
       const Acts::GeometryContext& gctx) const {
      
      std::cout << "DEBUG " << getpid() << " wait for signal CONT " << std::endl;
      kill(getpid(), SIGSTOP);

      const nlohmann::json jout = trackingGeometryToJson(tracking_geometry,gctx);
      // And write the file
      std::ofstream ofj(file_name);
      ofj << jout << std::endl;
   }

   /// Create a json object which describes the tracking geometry
   /// @return json object describing the tracking geometry
   const nlohmann::json TrackingGeometryJsonWriter::trackingGeometryToJson(
      const std::shared_ptr<const Acts::TrackingGeometry> &tracking_geometry,
      const Acts::GeometryContext& gctx) const {

      TrackingGeometryJsonDecorator decorator(gctx);
      // Setup the converter config
      Acts::MaterialMapJsonConverter::Config cfg { m_cfg } ;
      cfg.context = gctx;

      // Evoke the converter
      Acts::MaterialMapJsonConverter jmConverter(cfg, m_cfg.logLevel);
      return jmConverter.trackingGeometryToJson(*tracking_geometry, &decorator);
   }
}

