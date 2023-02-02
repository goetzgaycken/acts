#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

#include "Acts/Plugins/Json/TrackingGeometryJsonReader.hpp"
#include "Acts/Plugins/Json/VolumeBoundsJsonConverter.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"

#include "Acts/Geometry/KDTreeTrackingGeometryBuilder.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"

#include "GeantinoReader.hpp"

#include <memory>
#include <array>
#include <deque>
#include <limits>
#include <stdexcept>

// Stat
#include <cmath>
#include <algorithm>
#include <vector>
#include <ostream>

// dump geom
#include <deque>

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
      double min()  const { return m_min; }
      double max()  const { return m_max; }

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

   bool overlaps(const Acts::Extent &a,  const Acts::Extent &b, Acts::BinningValue bValue) {
      return a.min(bValue) < b.max(bValue) && a.max(bValue) > b.min(bValue);
   }

   void dumpProtoVolumes(const Acts::ProtoDetector &detector) {
      // dump volume hierarchy
   std::deque<const Acts::ProtoVolume *> childs {&detector.worldVolume};
   while (!childs.empty()) {
      const Acts::ProtoVolume *a_child = childs.front();
      childs.pop_front();
      if (a_child and a_child->container.has_value()) {
         unsigned int container=0;
         unsigned int layer_contaier=0;
         for (const Acts::ProtoVolume &elm : a_child->container.value().constituentVolumes) {
            if (elm.container.has_value() and elm.container.value().layerContainer) {
               ++layer_contaier;
            }
            else {
               ++container;
            }
         }
         std::cout << "/---" << a_child->name
                   << ( /*a_child->layerContainer*/ true ?
                        ((a_child->internal.has_value()
                          and a_child->internal.value().layerType == Acts::Surface::SurfaceType::Disc) ? " (disc layers)"
                         : (a_child->internal.has_value()
                            and a_child->internal.value().layerType == Acts::Surface::SurfaceType::Cylinder) ? " (cylinder layers)"
                         : " (unknown layer type)")
                        : "")
                   << ( a_child->container.has_value() and  a_child->container.value().layerContainer ? " (layerContainer)" : "");
         if (container>0) {
            std::cout << " tracking volumes: " << container;
         }
         if (layer_contaier) {
            std::cout << " layer container: " << layer_contaier;
         }
         std::cout << std::endl;
         std::cout << a_child->extent.toString("    ") << std::endl;
         if (a_child->container.has_value()) {
            for (const Acts::ProtoVolume &elm : a_child->container.value().constituentVolumes) {
               std::cout << "| " << elm.name
                         << ( /*elm.layerContainer*/true ?
                              ((elm.internal.has_value()
                                and elm.internal.value().layerType == Acts::Surface::SurfaceType::Disc) ? " (disc layers)"
                               : (elm.internal.has_value()
                                  and elm.internal.value().layerType == Acts::Surface::SurfaceType::Cylinder) ? " (cylinder layers)"
                               : " (unknown layer type)")
                              : "")
                         << ( (elm.container.has_value()
                               and elm.container.value().layerContainer) ? " (layerContainer)" : "");
               if (!a_child->container.value().constituentBinning.empty()) {
                  std::cout << " (binning:";
                  for (const Acts::BinningData &binning : a_child->container.value().constituentBinning ) {
                     std::cout << " " << Acts::binningValueNames().at(binning.binvalue);
                  }
                  std::cout << ") ";
               }
                              std::cout << " type:" << (elm.internal.has_value()
                                                        ? a_child->internal.value().layerType : Acts::Surface::Other ) << std::endl;
               std::cout << "| " << elm.extent.toString("|   ") << std::endl;
               childs.push_back( &elm );
            }
         }
         std::cout << "\\---" << a_child->name << std::endl << std::endl;
      }
   }
   }

   void dumpTrackingVolume(const Acts::TrackingVolume *volume, const std::string &margin) {
      if (volume) {

         std::cout << margin << volume->volumeName() << std::endl;
         if (volume->confinedLayers()) {
         Stat n_surfaces;
         Stat n_surfaces_with_det;
         Stat n_sensitive_surfaces;
         Stat n_boundary_surfaces;
         Stat n_approach_surfaces;
         for (const std::shared_ptr<const Acts::Layer> &a_layer : volume->confinedLayers()->arrayObjects()) {
            if (a_layer->surfaceArray()) {
               unsigned short n_surfaces_counter = 0;
               unsigned short n_surfaces_with_det_counter = 0;
               unsigned short n_sensitive_surfaces_counter = 0;
               unsigned short n_boundary_surfaces_counter = 0;
               unsigned short n_approach_surfaces_counter = 0;
               for ( const Acts::Surface *surface : a_layer->surfaceArray()->surfaces() ) {
                  ++n_surfaces_counter;
                  if (surface->geometryId().boundary()!= 0u) {
                     ++n_boundary_surfaces_counter;
                  }
                  if (surface->geometryId().approach()!= 0u) {
                     ++n_approach_surfaces_counter;
                  }
                  if (surface->geometryId().sensitive()!= 0u) {
                     ++n_sensitive_surfaces_counter;
                  }
               }
               if (n_surfaces_counter>0) {
                  n_surfaces.add(n_surfaces_counter);
                  n_surfaces_with_det.add(n_surfaces_with_det_counter);
                  n_sensitive_surfaces.add(n_sensitive_surfaces_counter);
                  n_boundary_surfaces.add(n_boundary_surfaces_counter);
                  n_approach_surfaces.add(n_approach_surfaces_counter);
               }
            }
         }
         std::vector<std::pair<std::string, Stat *> > data;
         if (n_surfaces.mean() > 0. ) {
            data.reserve(5);
            data.push_back(std::make_pair("surfaces", &n_surfaces) );
            data.push_back(std::make_pair("surfaces with associated detectors", &n_surfaces_with_det) );
            data.push_back(std::make_pair("sensitve surfaces", &n_sensitive_surfaces) );
            data.push_back(std::make_pair("approach surfaces", &n_approach_surfaces) );
            data.push_back(std::make_pair("boundary surfaces", &n_boundary_surfaces) );
            for (const std::pair<std::string, Stat *> &a_data : data ) {
               if (a_data.second) {
                  std::cout << margin << " " << *a_data.second <<  " | " << a_data.first << std::endl;
               }
            }
         }
         }
      }
   }
   void dumpTrackingVolumeHierarchy(const Acts::TrackingVolume *volume) {
      if (volume) {
         std::vector<std::pair<const Acts::TrackingVolume *, std::string> > stack;
         stack.push_back(std::make_pair(volume,std::string("")));
         while (!stack.empty()) {
            auto [a_volume,margin] = stack.back();
            stack.pop_back();
            if (a_volume) {
               dumpTrackingVolume(a_volume, margin);
               if (a_volume->confinedVolumes()) {
                  for (const std::shared_ptr<const Acts::TrackingVolume> &child_volume : a_volume->confinedVolumes()->arrayObjects()) {
                     stack.push_back(std::make_pair(child_volume.get(), margin + "  "));
                  }
               }
            }
         }
      }
   }
   void dumpTrackingGeometry(const Acts::TrackingGeometry *geom) {
      if (geom) {
         dumpTrackingVolumeHierarchy(geom->highestTrackingVolume());
      }
   }
}


namespace Acts {
TrackingGeometryJsonReader::VolumeInfo &
TrackingGeometryJsonReader::registerVolume(unsigned int volume_id,
                                           const std::string &name,
                                           std::unordered_map<unsigned int, VolumeInfo > &volume_map,
                                           std::unordered_map<std::string, unsigned int > &name_map) {
   if (name.empty()) {
      throw std::runtime_error("volume without name");
   }

   std::pair<std::unordered_map<std::string, unsigned int>::iterator, bool>
      insert_key = name_map.insert( std::make_pair(name, name_map.size()));
   unsigned int new_key = insert_key.first->second;
   std::pair<std::unordered_map<unsigned int, VolumeInfo>::iterator, bool>
      insert_result = volume_map.insert( std::make_pair(new_key, VolumeInfo(name, volume_id)));
   VolumeInfo &new_element = insert_result.first->second;

   if (name[0]=='{') {
      int counter=0;
      std::string::size_type pos=2;
      std::string::size_type end_pos=name.size();
      std::string::size_type start_pos=pos;
      for (;pos<end_pos; ++pos) {
         int next_counter=counter;
         switch (name[pos]) {
         case '{': {
            ++pos; // skip whitespace following an opening bracket;
            ++next_counter;
            break;
         }
         case '}': {
            --next_counter;
            [[fallthrough]];
         }
         case '|': {
            if (counter<=0) {
               std::string::size_type name_end = ( pos > 0 ? pos - 1 : pos );
               std::string child_name = name.substr(start_pos,name_end-start_pos);
               std::pair<std::unordered_map<std::string, unsigned int>::iterator, bool>
                  insert_child_key = name_map.insert( std::make_pair(child_name, name_map.size()));
               unsigned int child_key = insert_child_key.first->second;
               //               std::cout << "DEBUG child " << child_key << " <" << child_name << ">" << std::endl;
               std::pair<std::unordered_map<unsigned int, VolumeInfo>::iterator, bool>
                  child_insert_result = volume_map.insert( std::make_pair(child_key,VolumeInfo(child_name)));

               child_insert_result.first->second.parent = new_key;
               new_element.childs.push_back(child_insert_result.first->first);
               start_pos=pos+2;
            }
            ++pos; // skip whitespace following a closing bracket;
            break;
         }
         }
         counter=next_counter;
      }
   }
   if (!insert_result.second) {
      insert_result.first->second.id     = volume_id;
   }
   return insert_result.first->second;
}

std::shared_ptr<const Acts::TrackingGeometry>
 TrackingGeometryJsonReader::read(const std::string &file_name, std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) const {

 std::ifstream f(file_name);
 nlohmann::json tracking_geometry_description = nlohmann::json::parse(f);
 return createTrackingGeometry(tracking_geometry_description, std::move(mdecorator));
}

Acts::ProtoDetector TrackingGeometryJsonReader::createProtoDetector(
 const nlohmann::json& tracking_geometry_description,
 const std::vector<std::shared_ptr<Acts::Surface>> &surfaces,
 const std::vector< std::vector< unsigned int >> &surfaces_per_volume) const {

   Acts::GeometryContext gctx;

   Acts::ProtoDetector detector;

   nlohmann::json volumes = tracking_geometry_description["Volumes"];
   int version = volumes["acts-geometry-hierarchy-map"]["format-version"].get<int>();
   static std::array<int,1> supported_versions{0};
   if (std::find(supported_versions.begin(),supported_versions.end(), version) == supported_versions.end()) {
      throw std::runtime_error("Unsupported input format.");
   }
   static const std::array<std::string, VolumeInfo::kNBinningVariables> binning_variables {
      std::string("binR"),
      std::string("binPhi"),std::string("binZ")};
   static const std::array<std::string, VolumeInfo::Binning::kNRangeTypes> binning_range_types {
      std::string("open"),
      std::string("closed")};
   static const std::map<std::string, Acts::Surface::SurfaceType> supported_layer_types {
      { std::string("Disc"), Acts::Surface::SurfaceType::Disc },
      { std::string("Cylinder"), Acts::Surface::SurfaceType::Cylinder }
   };

   static std::array<Acts::BinningValue, VolumeInfo::eNDomains> acts_domain_type
      {Acts::binR, Acts::binZ}; // must match order in VolumeInfos::EDomainTypes

   std::unordered_map<std::string, unsigned int > name_map;
   std::unordered_map<unsigned int, VolumeInfo> hierarchy;
   nlohmann::json volume_list  = volumes.at("entries");
   bool is_cylinder=true;
   bool is_z_axis_aligned=true;

   // create proto volumes for each volume and each layer of volume
   for( nlohmann::json a_volume : volume_list) {
      unsigned int volume_id = a_volume["volume"].get<unsigned int>();
      nlohmann::json value=a_volume["value"];

      VolumeInfo &element = registerVolume(volume_id, value["NAME"].get<std::string>(), hierarchy, name_map);
      element.volumeBounds = Acts::unqiueVolumeBoundsFromJson(value["bounds"]);

      nlohmann::json volume_transform = value["transform"];
      if (volume_transform.find("rotation") != volume_transform.end() and not volume_transform["rotation"].empty()) {
         // @TODO check that matrix is not identity
         is_z_axis_aligned = false;
      }
      Acts::from_json(value["transform"], element.volumeTransform);
      if (   not floatsAgree(element.volumeTransform.translation()[0],0.)
          or not floatsAgree(element.volumeTransform.translation()[1],0.)) {
         std::cout << "ERROR not a z-axis aligned tracking volume " << element.proto_volume.name << std::endl;
         is_z_axis_aligned = false;
      }

      const Acts::CylinderVolumeBounds *cylinder_bounds=dynamic_cast<const Acts::CylinderVolumeBounds *>(element.volumeBounds.get());

      if (not element.volumeBounds or element.volumeBounds->type() != Acts::VolumeBounds::eCylinder or not cylinder_bounds) {
         std::cout << "ERROR not a cylinder tracking volume " << element.proto_volume.name << std::endl;
         is_cylinder=false;
      }
      if (is_cylinder and is_z_axis_aligned) {
         element.domains.at(0) = VolumeInfo::Domain(cylinder_bounds->get(Acts::CylinderVolumeBounds::eMinR),
                                                    cylinder_bounds->get(Acts::CylinderVolumeBounds::eMaxR));
         element.domains.at(1) = VolumeInfo::Domain(  element.volumeTransform.translation()[2]
                                                    - cylinder_bounds->get(Acts::CylinderVolumeBounds::eHalfLengthZ),
                                                      element.volumeTransform.translation()[2]
                                                    + cylinder_bounds->get(Acts::CylinderVolumeBounds::eHalfLengthZ));
      }
      if (value.contains("layers")) {
         nlohmann::json layer_list  = value.at("layers");
         if (not layer_list.empty()) {
            element.proto_volume.container = Acts::ProtoVolume::ContainerStructure{};
         for( nlohmann::json a_layer : layer_list) {
            unsigned int layer_id = a_layer["layer"].get<unsigned int>();
            std::shared_ptr<Acts::Surface> layer_surface = Acts::surfaceFromJson(a_layer);
            if (!layer_surface || (layer_surface->type() != Acts::Surface::Disc && layer_surface->type() != Acts::Surface::Cylinder)) {
                  std::cout << "ERROR layer " << layer_id << " of volume " << element.proto_volume.name
                            << " is not a disc or cylinder but " << a_layer["type"].get<std::string>()
                            << std::endl;
            }
            else {
               nlohmann::json layer_transform_json = a_layer["transform"];
               bool layer_is_z_axis_aligned=true;
               if (layer_transform_json.find("rotation") != layer_transform_json.end()
                   and not layer_transform_json["rotation"].empty()) {
                  // @TODO check that matrix is not identity
                  layer_is_z_axis_aligned = false;
               }
               const Acts::Transform3& layer_transform = layer_surface->transform(gctx);
               if (   not floatsAgree(layer_transform.translation()[0],0.)
                   or not floatsAgree(layer_transform.translation()[1],0.)) {
                  layer_is_z_axis_aligned = false;
               }
               if (!layer_is_z_axis_aligned) {
                  std::cout << "ERROR layer " << layer_id << " of volume " << element.proto_volume.name
                            << " is not z-axis aligned." << std::endl;
               }
               else {
                  double layer_thickness = a_layer["thickness"].get<double>();
                  ProtoVolume layer_volume;
                  layer_volume.container=Acts::ProtoVolume::ContainerStructure{};
                  layer_volume.container.value().layerContainer=true;
                  layer_volume.internal = ProtoVolume::InternalStructure{};
                  {
                     std::stringstream tmp;
                     tmp << element.proto_volume.name << "_" << layer_id;
                     layer_volume.name = tmp.str();
                  }

                  // determine lauer type from the bounds of the surface representation
                  switch (layer_surface->bounds().type()) {
                  case Acts::SurfaceBounds::eCylinder: {
                     const Acts::CylinderBounds &layer_bounds = static_cast<const Acts::CylinderBounds &>(layer_surface->bounds());
                     layer_volume.extent.set(Acts::binR,
                                             layer_bounds.values().at(Acts::CylinderBounds::eR) - layer_thickness *.5,
                                             layer_bounds.values().at(Acts::CylinderBounds::eR) + layer_thickness *.5);
                     layer_volume.extent.set(Acts::binZ,
                                             layer_transform.translation()[2]
                                             - layer_bounds.values().at(Acts::CylinderBounds::eHalfLengthZ),
                                             layer_transform.translation()[2]
                                             + layer_bounds.values().at(Acts::CylinderBounds::eHalfLengthZ));
                     layer_volume.internal.value().layerType = Acts::Surface::SurfaceType::Cylinder;
                     break;
                  }
                  case Acts::SurfaceBounds::eDisc: {
                     const Acts::RadialBounds &layer_bounds = static_cast<const Acts::RadialBounds &>(layer_surface->bounds());
                     layer_volume.extent.set(Acts::binR,layer_bounds.rMin(), layer_bounds.rMax());
                     layer_volume.extent.set(Acts::binZ,
                                             layer_transform.translation()[2] - layer_thickness *.5,
                                             layer_transform.translation()[2] + layer_thickness *.5);
                     layer_volume.internal.value().layerType = Acts::Surface::SurfaceType::Disc;
                     break;
                  }
                  default: {
                     std::cout << "ERROR layer " << layer_id << " of volume " << element.proto_volume.name
                               << " has unsupported bounds :  "  << typeid(layer_surface->bounds() ).name()
                               << std::endl;
                     break;
                  }
                  }
                  element.proto_volume.container.value().constituentVolumes.push_back(std::move(layer_volume));
               }
            }
         }
         element.proto_volume.container.value().layerContainer = true;
         // identify non overlapping domains of the layers
         unsigned int child_domain_overlaps=0;
         for (unsigned int child_i=0;
              child_i<element.proto_volume.container.value().constituentVolumes.size();
              ++child_i) {
            for (unsigned int child_j=child_i+1;
                 child_j<element.proto_volume.container.value().constituentVolumes.size();
                 ++child_j) {
               unsigned int domain_i;
               for (Acts::BinningValue bValue : acts_domain_type) {
                  child_domain_overlaps |= (overlaps(element.proto_volume.container.value().constituentVolumes.at(child_i).extent,
                                                     element.proto_volume.container.value().constituentVolumes.at(child_j).extent,
                                                     bValue)
                                            << domain_i);
                  ++domain_i;
               }
            }
         }
         unsigned int non_overlapping_domains = 0;
         for(unsigned int domain_i =0 ; domain_i < acts_domain_type.size(); ++domain_i) {
            non_overlapping_domains |= (1<<domain_i);
         }
         element.setDomainMask = non_overlapping_domains;
         non_overlapping_domains &= ~child_domain_overlaps;
         element.noOverlapMask  |= non_overlapping_domains;

         // enable binning for all domains in which the layers do not overlap
         for(unsigned int domain_i =0 ; domain_i < acts_domain_type.size(); ++domain_i) {
            if (    (element.noOverlapMask & (1<<domain_i) ) ) {
               element.proto_volume.container.value().constituentBinning = {
                  Acts::BinningData(Acts::closed,
                                    acts_domain_type.at(domain_i),
                  {0., 1.})};
            }
         }
         }


      }
   }
   // the tracking volumes must by cylindrical and must be z-axis aligned.
   // everything else is not supported
   if (not is_cylinder or not is_z_axis_aligned) return detector;

   std::deque<unsigned int> element_queue;
   VolumeInfo *world=nullptr;

   // find proto volumes without constituents. These are the volumes to process first.
   // also set domain mask i.e. a mask indicating which domains are identical for all constituents
   // and a the noOverlapMask i.e. a mask indicating the domains in which constituents do not overlap
   for(std::pair<const unsigned int, VolumeInfo> &element : hierarchy) {
      if (element.second.parent == std::numeric_limits<unsigned int>::max()) {
         world=&element.second;
      }
      if (element.second.childs.empty()) {
         element_queue.push_back(element.first);
      }
      else {
         // identify common and overlapping domains of the constituents
         unsigned int non_identical_domain_mask=0;
         unsigned int child_domain_overlaps=0;
         for(unsigned int child_i=0; child_i<element.second.childs.size(); ++child_i) {
            unsigned int child_idx = element.second.childs.at(child_i);
            VolumeInfo &a_child_element = hierarchy.at(child_idx);
            for(unsigned int domain_i =0 ;
                domain_i < a_child_element.domains.size();
                ++domain_i) {
               assert( a_child_element.domain.size() == element.second.domains.size());
               if ( a_child_element.domains.at(domain_i)
                   != element.second.domains.at(domain_i)) {
                  non_identical_domain_mask |= (1<<domain_i);
               }
               for(unsigned int other_child_i=0;
                   other_child_i<child_i; ++other_child_i) {
                  VolumeInfo &other_child_element
                     = hierarchy.at( element.second.childs.at(other_child_i) );
                  if (a_child_element.domains.at(domain_i).overlaps(
                         other_child_element.domains.at(domain_i) )) {
                     child_domain_overlaps |= (1<<domain_i);
                     break;
                  }
               }
            }
         }
         unsigned int non_overlapping_domains = 0;
         for(unsigned int domain_i =0 ; domain_i < element.second.domains.size(); ++domain_i) {
            non_overlapping_domains |= (1<<domain_i);
         }
         unsigned int same_domain_mask = non_overlapping_domains & (~non_identical_domain_mask);
         non_overlapping_domains &= ~child_domain_overlaps;

         for(unsigned int child_i=0; child_i<element.second.childs.size(); ++child_i) {
            VolumeInfo &a_child_element = hierarchy.at( element.second.childs.at(child_i) );
            if (a_child_element.childs.empty()) {
               a_child_element.setDomainMask |= non_overlapping_domains;
            }
         }
         // set domain bit if domain is identical to the domain definition of all childs;
         element.second.setDomainMask |= same_domain_mask;
         // in the parent also set the domain in which the constituents do not overlap
         // since this is a domain which must be specified by the parent.
         element.second.setDomainMask |= non_overlapping_domains;
         element.second.noOverlapMask  |= non_overlapping_domains;

      }
   }

   // \/ loop over all elements
   // \/ loop over childs
   // \/ identify binning data of childs which is not overlapping with the binning of any of the childs
   //    set as extends
   // \/ set identical binning as binning for parent
   //    unless closed phi binning and min == max
   // \/ put elements without childs to queue.

   // process element queue
   // take element from front deque if childs have been processed if yes
   // add childs as constituents, mark as processed, put parent to end of queue
   // if childs are not processed put to end of queue
   while (!element_queue.empty()) {
      unsigned int element_idx = element_queue.front();
      element_queue.pop_front();
      VolumeInfo &an_element = hierarchy.at(element_idx);
      bool childs_processed=true;
      unsigned short common_child_domains=0;
      for(unsigned int domain_i =0 ; domain_i < an_element.domains.size(); ++domain_i) {
         common_child_domains |= (1<<domain_i);
      }

      for(unsigned int child_i=0; child_i<an_element.childs.size(); ++child_i) {
         VolumeInfo &a_child_element = hierarchy.at(an_element.childs.at(child_i) );
         childs_processed &= a_child_element.processed;
         if (not childs_processed) { break; }
            common_child_domains &= a_child_element.setDomainMask;
      }
      if (childs_processed ) {
         // @TODO sort childs by z or r
         if (!an_element.childs.empty()) {
            an_element.proto_volume.container = Acts::ProtoVolume::ContainerStructure{};
            an_element.proto_volume.container.value().constituentVolumes.reserve(an_element.childs.size());
            bool layer_container = hierarchy.at(an_element.childs.at(0)).proto_volume.container.value().layerContainer;
            bool mixed=false;

            for(unsigned int child_i=1; child_i<an_element.childs.size(); ++child_i) {
               VolumeInfo &a_child_element = hierarchy.at(an_element.childs.at(child_i) );
               if (a_child_element.proto_volume.container.value().layerContainer != layer_container) {
                  mixed=true;
                  break;
               }
            }
            if (mixed) {
               std::cout << "WARNING constituent volumes of " << an_element.proto_volume.name
                         << " are partially layer container and volume container." << std::endl;
            }
         }
         Acts::BinningValue last_domain_type {};
         unsigned int n_domain_types=0;
         bool propagate_domains=false;
         // set extent of proto volume
         for(unsigned int domain_i =0 ; domain_i < an_element.domains.size(); ++domain_i) {
            if ( (an_element.setDomainMask & (1<<domain_i)) ) {
               ++n_domain_types;
               last_domain_type = acts_domain_type.at(domain_i);
               an_element.proto_volume.extent.set(acts_domain_type.at(domain_i),
                                                   an_element.domains.at(domain_i).min,
                                                   an_element.domains.at(domain_i).max);

               if ( !(common_child_domains & (1<<domain_i)) ) {
                  propagate_domains=true;
               }
            }
            if (!an_element.childs.empty()) {
               if (    (an_element.noOverlapMask & (1<<domain_i) ) ) {
                  an_element.proto_volume.container.value().constituentBinning = {
                     Acts::BinningData(Acts::closed,
                                      acts_domain_type.at(domain_i),
                                       {0., 1.})};

               }
            }
         }
         // if the domain is not defined for any of the childs propagate the domain down to
         // all childs for which it is not defined and which have childs themselves
         if (propagate_domains) {
            std::vector<VolumeInfo *> parent_stack { &an_element };
            while (!parent_stack.empty()) {
               VolumeInfo &parent_element = *parent_stack.back();
               parent_stack.pop_back();
               for(unsigned int child_i=0; child_i<parent_element.childs.size(); ++child_i) {
                  VolumeInfo &a_child_element = hierarchy.at(parent_element.childs.at(child_i) );
                  //                  if (!a_child_element.childs.empty())
                  {
                     for(unsigned int domain_i =0 ; domain_i < an_element.domains.size(); ++domain_i) {
                        if (!(a_child_element.setDomainMask & (1<<domain_i))
                            && (an_element.setDomainMask & (1<<domain_i))) {
                           a_child_element.proto_volume.extent.set(acts_domain_type.at(domain_i),
                                                                   an_element.domains.at(domain_i).min,
                                                                   an_element.domains.at(domain_i).max);
                           a_child_element.setDomainMask |=  (1<<domain_i);
                           parent_stack.push_back(&a_child_element);
                        }
                     }
                  }
               }
            }
         }

         if (n_domain_types==1) {
            switch (last_domain_type){
            case Acts::binR: {
               if (an_element.proto_volume.container.value().layerContainer) {
                  an_element.proto_volume.internal.value().layerType = Acts::Surface::SurfaceType::Cylinder;
               }
               break;
            }
            case Acts::binZ: {
               if (an_element.proto_volume.container.value().layerContainer) {
                  an_element.proto_volume.internal.value().layerType = Acts::Surface::SurfaceType::Disc;
               }
               break;
            }
            default:
               break;
            }
         }
         an_element.processed = childs_processed;

         an_element.proto_volume.container.value().constituentVolumes.reserve( an_element.childs.size() );
         bool childs_have_constituents=false;
         Acts::Surface::SurfaceType child_surface_type=Acts::Surface::Other;
         for(unsigned int child_i=0; child_i<an_element.childs.size(); ++child_i) {
            an_element.proto_volume.container.value().constituentVolumes.push_back( std::move(hierarchy.at(an_element.childs.at(child_i) ).proto_volume) );
            child_surface_type = ((child_i==0 && an_element.proto_volume.container.value().constituentVolumes.back().internal.has_value())
                                  || (an_element.proto_volume.container.value().constituentVolumes.back().internal.has_value()
                                      && child_surface_type
                                      == an_element.proto_volume.container.value().constituentVolumes.back().internal.value().layerType))
               ? an_element.proto_volume.container.value().constituentVolumes.back().internal.value().layerType
               : Acts::Surface::Other;
            if (an_element.proto_volume.container.value().constituentVolumes.back().container.has_value()
                and not an_element.proto_volume.container.value().constituentVolumes.back().container.value().constituentVolumes.empty()) {
               childs_have_constituents=true;
            }
         }

         if (an_element.parent < hierarchy.size()) {
            if (std::find(element_queue.begin(),element_queue.end(),an_element.parent) == element_queue.end()) {
               element_queue.push_back(an_element.parent);
            }
         }
      }
      else {
         element_queue.push_back(element_idx);
      }
   }

   for (const std::pair<const unsigned int, VolumeInfo> &element : hierarchy) {
      if (not element.second.processed) {
         throw std::logic_error("Not all volumes have been processed.");
      }
   }
   if (world) {
      detector.worldVolume = std::move(world->proto_volume);
   }
   //dumpProtoVolumes(detector);

   return detector;
}
namespace {
   struct LayerId : std::pair<unsigned int, unsigned int> {
      LayerId (unsigned int volume_id, unsigned int layer_id )
         : std::pair<unsigned int, unsigned int>( std::make_pair(volume_id, layer_id)) {}
      unsigned int volumeId() const { return this->first; }
      unsigned int layerId() const { return this->second; }
   };
   struct LayerSurfaceInfo {
      unsigned int nSensitiveSurfaces = 0;
      std::vector< std::shared_ptr<Acts::Surface> > nonSensitiveSurfaces;
   };
}

std::pair<std::vector<std::shared_ptr<Acts::Surface>>, std::vector< std::vector< unsigned int >> >
TrackingGeometryJsonReader::createSurfaces(const nlohmann::json& tracking_geometry_description) const {
   nlohmann::json surfaces_in = tracking_geometry_description["Surfaces"];
   int version = surfaces_in["acts-geometry-hierarchy-map"]["format-version"].get<int>();
   static std::array<int,1> supported_versions{0};
   if (std::find(supported_versions.begin(),supported_versions.end(), version) == supported_versions.end()) {
      throw std::runtime_error("Unsupported input format.");
   }
   nlohmann::json surface_list  = surfaces_in.at("entries");


   std::pair<std::vector<std::shared_ptr<Acts::Surface>>, std::vector< std::vector< unsigned int >> >
      surfaces_out;
   surfaces_out.first.reserve( surface_list.size() );
   std::map<LayerId,LayerSurfaceInfo> layer_surface_map;

   unsigned int counter_i=0;
   for( nlohmann::json a_surface : surface_list) {
      // skip boundary surfaces
      // @TODO check for sensitive surfaces
      if (a_surface.contains("boundary") || a_surface.contains("approach")) {
         continue;
      }
      //      unsigned int volume_id = a_volume["volume"].get<unsigned int>();
      /// @return a shared_ptr to a surface object for type polymorphism
      std::shared_ptr<Acts::Surface> surface_ptr( surfaceFromJson(a_surface["value"]) );
      if (not surface_ptr.get()) {
         uint64_t geo_id = std::numeric_limits<uint64_t>::max();
         uint64_t det_id = std::numeric_limits<uint64_t>::max();
         if (a_surface.contains("value")) {
            if (a_surface["value"].contains("geo_id") ) {
               geo_id = a_surface["value"]["geo_id"].get<uint64_t>();
            }

            if (a_surface["value"].contains("detID") ) {
               det_id = a_surface["value"]["detID"].get<uint64_t>();
            }
         }
         std::cout << "WARNING failed to create surface " << counter_i << " IDs geo: " << geo_id
                   << " detector " << det_id
                   << std::endl;
      }
      else {

         LayerId  layer_id(a_surface["volume"].get<unsigned int>(),
                           a_surface.contains("layer") ? a_surface["layer"].get<unsigned int>() : 0);

         std::pair<std::map< LayerId, LayerSurfaceInfo>::iterator, bool >
            insert_ret = layer_surface_map.insert(std::make_pair( layer_id, LayerSurfaceInfo()));

         if (a_surface.contains("sensitive")) {
            surfaces_out.first.push_back( std::move(surface_ptr) );
            if (layer_id.volumeId() >= surfaces_out.second.size()) {
               surfaces_out.second.resize( layer_id.volumeId() + 1 );
            }
            surfaces_out.second.at(layer_id.volumeId()).push_back(surfaces_out.first.size()-1);
            ++insert_ret.first->second.nSensitiveSurfaces;
         }
         else {
            insert_ret.first->second.nonSensitiveSurfaces.push_back( std::move(surface_ptr));
         }
      }
      ++counter_i;
   }
   // for layers without sensitive surfaces add
   unsigned int sensitive=0;
   unsigned int non_sensitive=0;

   for (std::pair<const LayerId, LayerSurfaceInfo> &surface_info : layer_surface_map ) {
      std::cout << "DEBUG volume " << surface_info.first.volumeId()
                << " layer " << surface_info.first.layerId()
                << " sensitive " << surface_info.second.nSensitiveSurfaces
                << " non-sensitive " << surface_info.second.nonSensitiveSurfaces.size()
                << std::endl;
      sensitive+=surface_info.second.nSensitiveSurfaces;
      non_sensitive += surface_info.second.nonSensitiveSurfaces.size();
      if (surface_info.second.nSensitiveSurfaces<1 && surface_info.first.layerId() % 2 == 0) {
         unsigned int idx=0;
         if (surface_info.second.nonSensitiveSurfaces.size() != 1) {
            std::cout << "WARNING volume " << surface_info.first.volumeId()
                      << " layer " << surface_info.first.layerId()
                      << " has not any sensitive surfaces, but "
                      << surface_info.second.nonSensitiveSurfaces.size() << " non sensitive surfaces."
                      << " Expected exactly one." << std::endl;
         }

         for ( std::shared_ptr<Acts::Surface> &a_surface_ptr : surface_info.second.nonSensitiveSurfaces ) {
            surfaces_out.first.push_back( std::move(a_surface_ptr) );
            if (surface_info.first.volumeId() >= surfaces_out.second.size()) {
               surfaces_out.second.resize( surface_info.first.volumeId() + 1 );
            }
            surfaces_out.second.at(surface_info.first.volumeId()).push_back(surfaces_out.first.size()-1);
            std::cout << "DEBUG add non sensitive surface " << idx << " for volume " << surface_info.first.volumeId()
                      << " layer " << surface_info.first.layerId()
                      << std::endl;
            ++idx;
         }

      }
      else {
      }
   }
   std::cout << "DEBUG surfaces-out total " <<  surfaces_out.first.size() << " sensitive: " << sensitive
             << " non-sensitive " << non_sensitive << std::endl;

   return surfaces_out;

}

std::shared_ptr<const Acts::TrackingGeometry>
TrackingGeometryJsonReader::createTrackingGeometry(const nlohmann::json& tracking_geometry_description,
                                                   std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) const {
   auto [surfaces, surfaces_per_volume] = createSurfaces( tracking_geometry_description);

   Acts::ProtoDetector detector(  createProtoDetector(tracking_geometry_description, surfaces, surfaces_per_volume) );
   std::cout << detector.toString("") << std::endl;
   detector.harmonize(true);
   std::cout << "DEBUG harmonize." << std::endl;
   //   dumpProtoVolumes(detector);


   // @TODO copied from Examples/Detectors/Geant4Detector/src/Geant4DetectorService.cpp.
   //    Move to a helper function ?
   Acts::GeometryContext tContext;

   // Surface array creatorr
   auto surfaceArrayCreator =
      std::make_shared<const Acts::SurfaceArrayCreator>(
          Acts::SurfaceArrayCreator::Config(),
          Acts::getDefaultLogger("SurfaceArrayCreator", m_cfg.toolLogLevel));
   // Layer Creator
   Acts::LayerCreator::Config lcConfig;
   lcConfig.surfaceArrayCreator = surfaceArrayCreator;
   auto layerCreator = std::make_shared<Acts::LayerCreator>(
       lcConfig, Acts::getDefaultLogger("LayerCreator", m_cfg.toolLogLevel));
   // Layer array creator
   Acts::LayerArrayCreator::Config lacConfig;
   auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
       lacConfig,
       Acts::getDefaultLogger("LayerArrayCreator", m_cfg.toolLogLevel));
   // Tracking volume array creator
   Acts::TrackingVolumeArrayCreator::Config tvacConfig;
   auto tVolumeArrayCreator =
      std::make_shared<const Acts::TrackingVolumeArrayCreator>(
          tvacConfig, Acts::getDefaultLogger("TrackingVolumeArrayCreator",
                                             m_cfg.toolLogLevel));
   // configure the cylinder volume helper
   Acts::CylinderVolumeHelper::Config cvhConfig;
   cvhConfig.layerArrayCreator = layerArrayCreator;
   cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
   auto cylinderVolumeHelper =
      std::make_shared<const Acts::CylinderVolumeHelper>(
          cvhConfig,
          Acts::getDefaultLogger("CylinderVolumeHelper", m_cfg.toolLogLevel));

   // The KDT tracking geometry builder
   Acts::KDTreeTrackingGeometryBuilder::Config kdtgConfig;
   kdtgConfig.layerCreator = layerCreator;
   kdtgConfig.trackingVolumeHelper = cylinderVolumeHelper;

   kdtgConfig.surfaces = std::move(surfaces);
   kdtgConfig.protoDetector = std::move(detector);
   kdtgConfig.materialDecorator = std::move(mdecorator);

   // Make the builder
   auto kdtTrackingGeometryBuilder = Acts::KDTreeTrackingGeometryBuilder(
      kdtgConfig, Acts::getDefaultLogger("KDTreeTrackingGeometryBuilder",
                                         m_cfg.toolLogLevel));

   std::shared_ptr<const Acts::TrackingGeometry> tmp(
      kdtTrackingGeometryBuilder.trackingGeometry(tContext).release());
   dumpTrackingGeometry( tmp.get() );

   if (! m_cfg.geantinoInputFileName.empty()) {
      std::map<uint64_t,uint64_t> id_map = GeantinoReader::createIdRemap(m_cfg.geantinoInputFileName, *tmp,
                                                                         static_cast<Long64_t>(m_cfg.maxGeantinoEntries) );
      for (const std::pair<const uint64_t, uint64_t> elm : id_map) {
         std::cout << "GeoId-map " << elm.first << "  " << elm.second << std::endl;
      }
   }
   return tmp;

}
}
