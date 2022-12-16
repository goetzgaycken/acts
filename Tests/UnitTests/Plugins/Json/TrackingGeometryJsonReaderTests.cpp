#include <boost/test/unit_test.hpp>
#include "Acts/Plugins/Json/TrackingGeometryJsonReader.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"

#include <fstream>
#include <iostream>

BOOST_AUTO_TEST_SUITE(TrackingGeometryJsonReader)

BOOST_AUTO_TEST_CASE(ProtoDetectorTest) {
   Acts::ProtoVolume vol_beam_pipe;
   vol_beam_pipe.name="BeamPipe::Barrel";
   vol_beam_pipe.layerType = Acts::Surface::SurfaceType::Cylinder;
   vol_beam_pipe.extent.set(Acts::binR,0., 25.434);
   Acts::ProtoVolume vol_pixel_inner_ecn;
   vol_pixel_inner_ecn.name="ITkPixelInner::NegativeEndcap";
   vol_pixel_inner_ecn.layerType = Acts::Surface::SurfaceType::Disc;
   vol_pixel_inner_ecn.extent.set(Acts::binZ,-3001.,-252.574);
   Acts::ProtoVolume vol_pixel_inner_barrel;
   vol_pixel_inner_barrel.name="ITkPixelInner::Barrel";
   vol_pixel_inner_barrel.layerType = Acts::Surface::SurfaceType::Cylinder;
   vol_pixel_inner_barrel.extent.set(Acts::binZ,-252.574, 252.574);
   Acts::ProtoVolume vol_pixel_inner_ecp;
   vol_pixel_inner_ecp.name="ITkPixelInner::PositiveEndcap";
   vol_pixel_inner_ecp.layerType = Acts::Surface::SurfaceType::Disc;
   vol_pixel_inner_barrel.extent.set(Acts::binZ,252.574,3001.);
   
   Acts::ProtoVolume vol_inner_pixel;
   vol_inner_pixel.name="{ ITkPixelInner::NegativeEndcap | ITkPixelInner::Barrel | ITkPixelInner::PositiveEndcap }";
   vol_inner_pixel.extent.set(Acts::binR, 25.434,128.641);
   vol_inner_pixel.extent.set(Acts::binZ,-3001.,3001.);
   vol_inner_pixel.constituentBinning = {
      Acts::BinningData(Acts::open,
                        Acts::binZ,
                        {0., 1.})};
   vol_inner_pixel.constituentVolumes.reserve(3);
   vol_inner_pixel.constituentVolumes.push_back( std::move( vol_pixel_inner_ecn) );
   vol_inner_pixel.constituentVolumes.push_back( std::move( vol_pixel_inner_barrel) );
   vol_inner_pixel.constituentVolumes.push_back( std::move( vol_pixel_inner_ecp) );

   
   Acts::ProtoVolume vol;
   vol.name="{ BeamPipe::Barrel | { ITkPixelInner::NegativeEndcap | ITkPixelInner::Barrel | ITkPixelInner::PositiveEndcap } }";
   vol.extent.set(Acts::binR,0. , 128.641);
   vol.constituentBinning = {
      Acts::BinningData(Acts::open,
                        Acts::binR,
                        {0., 1.})};
   vol.extent.set(Acts::binZ,-3001., 3001.);
   vol.constituentVolumes.reserve(2);
   vol.constituentVolumes.push_back( std::move( vol_beam_pipe) );
   vol.constituentVolumes.push_back( std::move( vol_inner_pixel) );

   Acts::ProtoDetector detector;
   detector.worldVolume = std::move(vol);
   detector.harmonize(true);
}
   
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



