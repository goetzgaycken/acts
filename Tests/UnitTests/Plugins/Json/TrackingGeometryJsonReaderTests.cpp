#include <boost/test/unit_test.hpp>
#include "Acts/Plugins/Json/TrackingGeometryJsonReader.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <fstream>
#include <iostream>

BOOST_AUTO_TEST_SUITE(TrackingGeometryJsonReader)

BOOST_AUTO_TEST_CASE(ReadJson) {

 std::ifstream f("/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/geometry-maps-volbounds.json");
 nlohmann::json tracking_geometry_description = nlohmann::json::parse(f);
 nlohmann::json &o=tracking_geometry_description;
 // for (nlohmann::json::iterator it = o.begin(); it != o.end(); ++it) {
 //    std::cout << it.key() << " : " << it.value() << "\n";
 // }

 Acts::TrackingGeometryJsonReader::Config  reader_cfg;
 reader_cfg.detectorName="ITK";
 reader_cfg.toolLogLevel = Acts::Logging::VERBOSE;
 reader_cfg.logLevel = Acts::Logging::VERBOSE;

 Acts::TrackingGeometryJsonReader reader(reader_cfg);
 std::shared_ptr<const Acts::TrackingGeometry>
    tracking_geometry=reader.createTrackingGeometry(tracking_geometry_description);
 if (tracking_geometry) {
    Acts::GeometryContext geo_ctx;
    tracking_geometry->visitSurfaces([&geo_ctx](const Acts::Surface *a_surface) {
          if (a_surface) {
             std::cout <<"DEBUG TGsurface ";
             a_surface->toStream(geo_ctx,std::cout);
             std::cout << std::endl;
          }
       });
 }

}

BOOST_AUTO_TEST_SUITE_END()



