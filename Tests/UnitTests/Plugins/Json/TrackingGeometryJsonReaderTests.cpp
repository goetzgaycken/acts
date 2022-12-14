#include <boost/test/unit_test.hpp>
#include "Acts/Plugins/Json/TrackingGeometryJsonReader.hpp"

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

 std::unique_ptr<Acts::TrackingGeometry> 
    tracking_geometry=Acts::TrackingGeometryJsonReader::trackingGeometry(tracking_geometry_description);

}

BOOST_AUTO_TEST_SUITE_END()



