#include <boost/test/unit_test.hpp>
#include "Acts/Plugins/Json/TrackingGeometryJsonReader.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"

#include <fstream>
#include <iostream>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <set>

namespace {

   inline double sqr(double a) { return a*a; }
   class Stat {
   public:
      Stat() = default;
      Stat(Stat &&) = default;
      Stat(const Stat &) = default;
      Stat(unsigned int nbins, double lower_edge, double upper_edge) {
         setBinning(nbins, lower_edge, upper_edge);
      }
      void add(double val) {
         m_n += 1.;
         m_sum += val;
         m_sum2 += val*val;
         m_min=std::min(val,m_min);
         m_max=std::max(val,m_max);
         if (!m_bins.empty()) {fill(val); }
      }
      double mean() const { return m_n>0. ? m_sum/m_n : 0.; }
      double rms()  const { return m_n>1. ? sqrt( (m_sum2 - m_sum*m_sum/m_n)/(m_n-1.)) : 0.; }
      double n()    const { return m_n; }
      double min()  const { return m_n>0 ? m_min : 0.; }
      double max()  const { return m_n>0 ? m_max : 0.; }

      Stat &operator =(const Stat &a) = default;
      Stat &operator =(Stat &&a) = default;
      Stat &operator +=(const Stat &a) {
         m_n += a.m_n;
         m_sum += a.m_sum;
         m_sum2 += a.m_sum2;
         m_min=std::min(a.m_min,m_min);
         m_max=std::max(a.m_max,m_max);
         if (m_bins.size() == a.m_bins.size()) {
            for (unsigned int bin_i=0; bin_i < m_bins.size(); ++bin_i) {
               m_bins[bin_i] += a.m_bins[bin_i];
            }
         }
         return *this;
      }
      double m_n   = 0.;
      double m_sum = 0.;
      double m_sum2 = 0.;
      double m_min = std::numeric_limits<double>::max();
      double m_max = -std::numeric_limits<double>::max();

      void setBinning(unsigned int n_bins, double lower_edge, double upper_edge) {
         m_bins.clear();
         m_bins.resize(n_bins+2,0);
         m_scale = n_bins / (upper_edge - lower_edge);
         m_lowerEdge=lower_edge-binWidth();
      }
      unsigned int bin( double value) {
         return value < m_lowerEdge ? 0 : std::min(m_bins.size()-1,
                                                   static_cast<std::size_t>( (value - m_lowerEdge) * m_scale));
      }
      double lowerEdge() const { return m_lowerEdge + binWidth(); }
      double upperEdge() const { return m_lowerEdge + (m_bins.size()-1)/m_scale; }
      double binWidth()  const { return 1/m_scale; }
      const std::vector<unsigned short> &bins() const { return m_bins; }
      void fill(double value) {
         ++m_bins.at(bin(value));
      }
      double m_lowerEdge;
      double m_scale;
      std::vector<unsigned short> m_bins;
   };
   std::ostream &operator<<(std::ostream &out, const Stat &a ) {
      out <<           std::setw(14) << a.min() << " < "
          <<           std::setw(14) << a.mean()
          << " +- " << std::setw(14) << a.rms()
          << " < "  << std::setw(14) << a.max()
          << " / "  << std::setw(9)  << a.n();
      if (!a.bins().empty()) {
         out << std::endl;
         out << a.lowerEdge() << ", " << a.lowerEdge()+a.binWidth() << ".." << a.upperEdge() << ": ";
         out << std::setw(5) << a.bins()[0] << " |";
         for (unsigned int bin_i=1; bin_i<a.bins().size()-1; ++bin_i) {
            out << " " << std::setw(5) << a.bins()[bin_i];
         }
         out << " | " << std::setw(5) << a.bins().back() ;
      }
      return out;
   }

   struct Counter {
      Counter &operator++() {
         ++m_value;
         return *this;
      }
      unsigned int m_value=0;
      std::vector<const Acts::Layer *> m_layer;
   };
}

namespace {
   void  dumpVolumes(const Acts::TrackingVolume* tracking_volume, std::string margin,
                     std::map<const Acts::Surface *,Counter> &duplicates,
                     std::map<const Acts::TrackingVolume *,Counter> &duplicate_volumes) {
      if (!tracking_volume) return;

      Stat n_surfaces;
      if (tracking_volume->confinedLayers()) {

         const std::vector<Acts::LayerPtr> &confined_layers = tracking_volume->confinedLayers()->arrayObjects();

         for (const Acts::LayerPtr &layer_ptr : confined_layers) {
            if (layer_ptr->surfaceArray()) {
               for (const Acts::Surface *a_surface : layer_ptr->surfaceArray()->surfaces()) {
                  std::pair<std::map<const Acts::Surface *,Counter>::iterator, bool>
                     ret_insert = duplicates.insert( std::make_pair(a_surface, Counter()));
                  ++ret_insert.first->second;
                  ret_insert.first->second.m_layer.push_back(layer_ptr.get());
                  if (ret_insert.first->second.m_value>1) {
                  std::pair<std::map<const Acts::TrackingVolume *,Counter>::iterator, bool>
                     ret_vol_insert = duplicate_volumes.insert( std::make_pair(tracking_volume, Counter()));
                  ++ret_vol_insert.first->second;
                  }
               }
               n_surfaces.add(layer_ptr->surfaceArray()->surfaces().size());
            }
         }
      }
      else {
         if (false) {
         unsigned int surface_counter=0;
         tracking_volume->visitSurfaces([&surface_counter,&duplicates,&duplicate_volumes, tracking_volume](const Acts::Surface *a_surface) {
               if (a_surface) {
                  ++surface_counter;
                  std::pair<std::map<const Acts::Surface *,Counter>::iterator, bool>
                     ret_insert = duplicates.insert( std::make_pair(a_surface, Counter()));
                  ++ret_insert.first->second;
                  if (ret_insert.first->second.m_value>1) {
                  std::pair<std::map<const Acts::TrackingVolume *,Counter>::iterator, bool>
                     ret_vol_insert = duplicate_volumes.insert( std::make_pair(tracking_volume, Counter()));
                  ++ret_vol_insert.first->second;
                  }
               }
            });
         n_surfaces.add(surface_counter);
         }
      }
      std::cout << "DEBUG " << margin << tracking_volume->volumeName()
                << std::endl
                << "DEBUG " << margin << "   " << n_surfaces << " surfaces / layer "
                << std::endl;
      if (tracking_volume->volumeName() == "ITkPixelInner::Barrel") {
         Acts::GeometryContext gctx;
         std::map<std::tuple<double, double, double>, std::pair<const Acts::Surface *,Counter> > sorted;
         tracking_volume->visitSurfaces([&gctx,tracking_volume, &sorted](const Acts::Surface *a_surface) {
               if (a_surface) {
                  Acts::Vector3 gpos = a_surface->center(gctx);
                  auto key = std::make_tuple(
                                                                 Acts::VectorHelpers::perp(gpos),
                                                                 Acts::VectorHelpers::phi(gpos),
                                                                 gpos[2]);
                  std::map<std::tuple<double, double, double>, std::pair<const Acts::Surface *,Counter> >::const_iterator
                     iter = sorted.find(key);
                  if (iter != sorted.end()) {
                     if (iter->second.first != a_surface) {
                        std::cout
                           << std::tuple<const Acts::Surface &, const Acts::GeometryContext&>(*iter->second.first,gctx) << std::endl
                           << std::tuple<const Acts::Surface &, const Acts::GeometryContext&>(*a_surface,gctx) << std::endl;

                     }
                  }

                  std::pair<std::map<std::tuple<double, double, double>, std::pair<const Acts::Surface *,Counter> >::iterator, bool>
                     ret =
                  sorted.insert( std::make_pair( std::make_tuple(
                                                                 Acts::VectorHelpers::perp(gpos),
                                                                 Acts::VectorHelpers::phi(gpos),
                                                                 gpos[2] ), std::make_pair(a_surface, Counter())));
                  ++ret.first->second.second;
               }
            });
         for (const std::pair< const std::tuple<double, double, double>, std::pair<const Acts::Surface *, Counter> > &elm : sorted ) {
         const Acts::Surface *a_surface = elm.second.first;
         Acts::Vector3 gpos = a_surface->center(gctx);
         std::cout  << tracking_volume->volumeName()
                    << " id: "
                    << std::setw(4) << a_surface->geometryId().volume()
                    << std::setw(4) << a_surface->geometryId().layer()
                    << " center: "
                    << std::setw(14) << Acts::VectorHelpers::perp(gpos)
                    << std::setw(14) << Acts::VectorHelpers::phi(gpos)
                    << std::setw(14) << gpos[2]
                    << " mult: " << std::setw(3) << elm.second.second.m_value
                    << " boundType: " <<  a_surface->bounds().type() << std::endl;
      }
      }
      margin+="  ";
      if (tracking_volume->confinedVolumes().get()) {
         for( const auto &tg  : tracking_volume->confinedVolumes()->arrayObjects() ) {
            dumpVolumes(tg.get(),margin,duplicates, duplicate_volumes);
         }
      }
      for( const auto &tg  : tracking_volume->denseVolumes() ) {
         dumpVolumes(tg.get(),margin, duplicates, duplicate_volumes);
      }
   }
}

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

 std::ifstream f("/data/goetz/ws/IDPVM/run/ITK_ttbar_mu200/geometry-maps-volbounds_2023.json");
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
 std::map<const Acts::Surface *,Counter> duplicates;
 std::map<const Acts::TrackingVolume *,Counter> duplicate_volumes;
 dumpVolumes(tracking_geometry->highestTrackingVolume(), "", duplicates, duplicate_volumes);
 Stat duplicate_stat(10,-0.5,9.5);
 Acts::GeometryContext gctx;
 std::cout << "Duplicates:" << std::endl;
 for (const std::pair<const Acts::Surface *,Counter> a_surface_stat : duplicates) {
    if (a_surface_stat.second.m_value > 1 ) {
       std::cout << "duplicated: " << a_surface_stat.second.m_value
                 << std::tuple<const Acts::Surface &, const Acts::GeometryContext&>(*a_surface_stat.first,gctx) << std::endl;
       unsigned int idx=0;
       std::set<const Acts::Layer *> used;
       for (const Acts::Layer *layer : a_surface_stat.second.m_layer ) {
          if (used.insert( layer ).second) {
             std::cout << std::setw(4) << idx << " "
                       << std::setw(4) << layer->geometryId().volume()
                       << std::setw(4) << layer->geometryId().layer()
                       << layer->trackingVolume()->volumeName() << ": " << std::endl
                       << std::tuple<const Acts::Surface &, const Acts::GeometryContext&>(layer->surfaceRepresentation(),gctx) << std::endl;
          }
       }
    }
    duplicate_stat.add( a_surface_stat.second.m_value);
 }
 std::cout << "duplicates " <<  duplicate_stat << std::endl;
 for (const std::pair<const Acts::TrackingVolume *,Counter> a_vol : duplicate_volumes) {
    std::cout << std::setw(6) << a_vol.second.m_value << " " << a_vol.first->volumeName() << std::endl;
 }

}

BOOST_AUTO_TEST_SUITE_END()

