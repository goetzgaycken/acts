#include "Acts/Plugins/Json/TrackingGeometryJsonWriter.hpp"
#include "Acts/Plugins/Json/TrackingGeometryJsonDecorator.hpp"
#include <fstream>

namespace Acts {


   void TrackingGeometryJsonWriter::writeNominal(
       const std::shared_ptr<const Acts::TrackingGeometry> &tracking_geometry,
       const std::string &file_name) const {
      write(tracking_geometry, file_name, Acts::GeometryContext());
   }

   void TrackingGeometryJsonWriter::write(
       const std::shared_ptr<const Acts::TrackingGeometry> &tracking_geometry,
       const std::string &file_name,
       const Acts::GeometryContext& gctx) const {

      const nlohmann::json jout = trackingGeometryToJson(tracking_geometry,gctx);
      // And write the file
      std::ofstream ofj(file_name);
      ofj << std::setw(4) << jout << std::endl;
   }

   class Counter {
   public:
      Counter &operator++() { ++m_counts; return *this; }
      unsigned int counts() const { return m_counts; }
      unsigned int m_counts = 0;
   };
   /// Create a json object which describes the tracking geometry
   /// @return json object describing the tracking geometry
   const nlohmann::json TrackingGeometryJsonWriter::trackingGeometryToJson(
      const std::shared_ptr<const Acts::TrackingGeometry> &tracking_geometry,
      const Acts::GeometryContext& gctx) const {

      std::map<std::string,Counter> volume_names;
      unsigned int n_surfaces=0;
      tracking_geometry->visitSurfaces([&volume_names,&n_surfaces](const Acts::Surface* surface) {
         ++n_surfaces;
         if (surface->associatedLayer()
             && surface->associatedLayer()->trackingVolume()) {            
            ++volume_names[surface->associatedLayer()->trackingVolume()->volumeName()];
         }
      });
      std::cout << "DEBUG tracking geometry with " << n_surfaces << " highest volume : "
                << tracking_geometry->highestTrackingVolume()->volumeName() << std::endl;
      for (const auto &elm : volume_names ) {
         std::cout << "DEBUG tracking volume " << elm.first << " has " << elm.second.counts() << " surfaces." << std::endl;
      }
      TrackingGeometryJsonDecorator decorator(gctx);
      // Setup the converter config
      Acts::MaterialMapJsonConverter::Config cfg { m_cfg } ;
      cfg.context = gctx;

      // Evoke the converter
      Acts::MaterialMapJsonConverter jmConverter(cfg, m_cfg.logLevel);
      return jmConverter.trackingGeometryToJson(*tracking_geometry, &decorator);
   }
}

