#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"

#include "Acts/Plugins/Json/TrackingGeometryJsonReader.hpp"
#include "Acts/Plugins/Json/VolumeBoundsJsonConverter.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"

#include "Acts/Geometry/KDTreeTrackingGeometryBuilder.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"

#include <memory>
#include <array>
#include <deque>
#include <limits>
#include <stdexcept>

void Acts::TrackingGeometryJsonReader::dumpMap(std::unordered_map<unsigned int, Acts::TrackingGeometryJsonReader::VolumeInfo > &volume_name_map) {
   for (const std::pair<const unsigned int, Acts::TrackingGeometryJsonReader::VolumeInfo > &elm : volume_name_map) {
      std::cout << elm.first << " -> " << elm.second.proto_volume.name << " " <<elm.second.id
                << " childs  ";
      for (unsigned int a_child :  elm.second.childs) {
         std::cout << " " << a_child;
      }
      std::cout  << std::endl;
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

Acts::ProtoDetector TrackingGeometryJsonReader::createProtoDetector(const nlohmann::json& tracking_geometry_description) {
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


   std::unordered_map<std::string, unsigned int > name_map;
   std::unordered_map<unsigned int, VolumeInfo> hierarchy;
   nlohmann::json volume_list  = volumes.at("entries");
   bool is_cylinder=true;
   bool is_z_axis_aligned=true;

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
      if (value.contains("layerType")) {
         std::map<std::string, Acts::Surface::SurfaceType>::const_iterator
            layer_type = supported_layer_types.find(value["layerType"].get<std::string>());
         if (layer_type != supported_layer_types.end()) {
            element.proto_volume.layerType = layer_type->second;
            std::cout << "DEBUG layer type for "<< element.proto_volume.name << " "
                      << layer_type->first << " -> "  << static_cast<unsigned int>(layer_type->second) 
                      << std::endl;
         }
         else {
            std::cout << "ERROR unsupported layer type for  " << element.proto_volume.name
                      << " (layer type " << value["layerType"] << ")"<< std::endl;
         }
      }

      // for (nlohmann::json &a_binning : value["material"]["binUtility"]["binningdata"]) {
      //    element.binning.at(getNamedEnum<unsigned short>(binning_variables,
      //                                                    a_binning["value"].get<std::string>(),
      //                                                    "Unknwown binning: "))
      //       = VolumeInfo::Binning(a_binning["min"].get<float>(),
      //                          a_binning["max"].get<float>(),
      //                          getNamedEnum<VolumeInfo::Binning::EType>(binning_range_types,
      //                                                                a_binning["option"].get<std::string>(),
      //                                                                "Unknown binning option: " ));
      // }
   }
   if (not is_cylinder or not is_z_axis_aligned) return detector;

   // for( nlohmann::json a_volume : volume_list) {
   //    unsigned int volume_id = a_volume["volume"].get<unsigned int>();
   //    nlohmann::json value=a_volume["value"];
   //    (void) registerVolumeInfo(volume_id, value["NAME"].get<std::string>(), hierarchy, true);
   // }

   std::deque<unsigned int> element_queue;
   VolumeInfo *world=nullptr;
   for(std::pair<const unsigned int, VolumeInfo> &element : hierarchy) {
      if (element.second.parent == std::numeric_limits<unsigned int>::max()) {
         world=&element.second;
      }
      if (element.second.childs.empty()) {
         element_queue.push_back(element.first);
      }
      else {
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
         std::cout << "DEBUG "
                   << element.second.proto_volume.name << " no-overlaps:";
         for(unsigned int domain_i =0 ; domain_i < element.second.domains.size();
             ++domain_i) {
            if ((1<<domain_i) & non_overlapping_domains) {
               std::cout << " " << domain_i;
            }
         }
         std::cout << " same:";
         for(unsigned int domain_i =0 ; domain_i < element.second.domains.size();
             ++domain_i) {
            if ((1<<domain_i) & same_domain_mask) {
               std::cout << " " << domain_i;
            }
         }
         std::cout << std::endl;

         for(unsigned int child_i=0; child_i<element.second.childs.size(); ++child_i) {
            VolumeInfo &a_child_element = hierarchy.at( element.second.childs.at(child_i) );
            if (a_child_element.childs.empty()) {
               a_child_element.setDomainMask |= non_overlapping_domains;
            }
         }
         // set domain if domain is identical to the domain definition of all childs;
         element.second.setDomainMask |= same_domain_mask;
         element.second.setDomainMask |= non_overlapping_domains;
         element.second.noOverlapMask  |= non_overlapping_domains;
      }
      //         if (element.second.parent < hierarchy.size()) {
      //         }
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
   for(unsigned int child_i=0; child_i<hierarchy.size(); ++child_i) {
      VolumeInfo &a_child_element = hierarchy.at(child_i);
      std::cout << "#" << child_i << " " << a_child_element.proto_volume.name << std::endl;
   }
   std::cout << std::endl;

   static std::array<Acts::BinningValue, VolumeInfo::eNDomains> acts_domain_type
      {Acts::binR, Acts::binZ}; // must match order in VolumeInfos::EDomainTypes
   while (!element_queue.empty()) {
      unsigned int element_idx = element_queue.front();
      element_queue.pop_front();
      VolumeInfo &an_element = hierarchy.at(element_idx);
      {
         std::cout << "DEBUG processing ";
         const VolumeInfo &tmp = an_element;
         if (tmp.moved) {
            std::cout << "(";
         }
         std::cout << " #" << element_idx;
         if (tmp.processed) {
            std::cout << "*";
         }
         if (tmp.moved) {
            std::cout << ")";
         }
         std::cout << std::endl;
      }
      if (an_element.moved) {
         std::cout << "ERROR processing moved element " << element_idx << std::endl;
      }
      bool childs_processed=true;
      for(unsigned int child_i=0; child_i<an_element.childs.size(); ++child_i) {
         VolumeInfo &a_child_element = hierarchy.at(an_element.childs.at(child_i) );
         childs_processed &= a_child_element.processed;
         if (not childs_processed) { break; }
      }
      if (childs_processed ) {
         // @TODO sort childs by z or r
         an_element.proto_volume.constituentVolumes.reserve(an_element.childs.size());
         std::cout << "DEBUG " << an_element.proto_volume.name << " childs:";
         for(unsigned int child_i=0; child_i<an_element.childs.size(); ++child_i) {
            VolumeInfo &a_child_element = hierarchy.at(an_element.childs.at(child_i) );
            std::cout << " " << a_child_element.proto_volume.name;
         }
         std::cout << std::endl;
         for(unsigned int domain_i =0 ; domain_i < an_element.domains.size(); ++domain_i) {
            if ( (an_element.setDomainMask & (1<<domain_i)) ) {
               an_element.proto_volume.extent.set(acts_domain_type.at(domain_i),
                                                   an_element.domains.at(domain_i).min,
                                                   an_element.domains.at(domain_i).max);
               std::cout << "DEBUG " << an_element.proto_volume.name << " set domain "
                         << domain_i << " (" << acts_domain_type.at(domain_i) << ")"
                          << " " << an_element.domains.at(domain_i).min
                          << " .. " << an_element.domains.at(domain_i).max
                         << std::endl;
            }
            if (!an_element.childs.empty()) {
               if (    (an_element.noOverlapMask & (1<<domain_i) ) ) {
                  std::cout << "DEBUG " << an_element.proto_volume.name << " set constituentDomain "
                            << " " << domain_i << "(" << acts_domain_type.at(domain_i) << ")"
                            << std::endl;

                  an_element.proto_volume.constituentBinning = {
                     Acts::BinningData(Acts::open,
                                      acts_domain_type.at(domain_i),
                                       {0., 1.})};

               }
            }
         }
         an_element.processed = childs_processed;

         an_element.proto_volume.constituentVolumes.reserve( an_element.childs.size() );
         for(unsigned int child_i=0; child_i<an_element.childs.size(); ++child_i) {
            hierarchy.at(an_element.childs.at(child_i)).moved=true;
            std::cout << "DEBUG move #" << an_element.childs.at(child_i) << " to #" << element_idx << std::endl;
            an_element.proto_volume.constituentVolumes.push_back( std::move(hierarchy.at(an_element.childs.at(child_i) ).proto_volume) );
            //            std::cout << " " << hierarchy.at(an_element.childs.at(child_i) ).proto_volume.name;
         }

         if (an_element.parent < hierarchy.size()) {
            if (std::find(element_queue.begin(),element_queue.end(),an_element.parent) == element_queue.end()) {
               element_queue.push_back(an_element.parent);
            }
         }
         std::cout << "DEBUG queue :";
         for (auto elm : element_queue) {
            const VolumeInfo &tmp = hierarchy.at(elm);
            if (tmp.moved) {
               std::cout << "(";
            }
            std::cout << " #" << elm;
            if (tmp.processed) {
               std::cout << "*";
            }
            if (tmp.moved) {
               std::cout << ")";
            }
         }
         std::cout << std::endl;
      }
      else {
         element_queue.push_back(element_idx);
      }
   }
   for (const std::pair<const unsigned int, VolumeInfo> &element : hierarchy) {
      if (not element.second.processed) {
         std::cout << "DEBUG unprocessed volume " << element.second.proto_volume.name << " childs:";
         for(unsigned int child_i=0; child_i<element.second.childs.size(); ++child_i) {
            std::cout << " " << hierarchy.at(element.second.childs.at(child_i) ).proto_volume.name;
         }
         std::cout << std::endl;
      }
   }
   if (world) {
      detector.worldVolume = std::move(world->proto_volume);
   }
   std::deque<const Acts::ProtoVolume *> childs {&detector.worldVolume};
   while (!childs.empty()) {
      const Acts::ProtoVolume *a_child = childs.front();
      childs.pop_front();
      if (a_child) {
         std::cout << "/---" << a_child->name << std::endl;
         std::cout << a_child->extent.toString("    ") << std::endl;
         for (const ProtoVolume &elm : a_child->constituentVolumes) {
            std::cout << "| " << elm.name  << std::endl;
            std::cout << "| " << elm.extent.toString("|   ") << std::endl;
            childs.push_back( &elm );
         }
         std::cout << "\\---" << a_child->name << std::endl << std::endl;
      }
   }
   return detector;
}


std::vector<std::shared_ptr<Acts::Surface>>
TrackingGeometryJsonReader::createSurfaces(const nlohmann::json& tracking_geometry_description) {
   nlohmann::json surfaces_in = tracking_geometry_description["Surfaces"];
   int version = surfaces_in["acts-geometry-hierarchy-map"]["format-version"].get<int>();
   static std::array<int,1> supported_versions{0};
   if (std::find(supported_versions.begin(),supported_versions.end(), version) == supported_versions.end()) {
      throw std::runtime_error("Unsupported input format.");
   }
   nlohmann::json surface_list  = surfaces_in.at("entries");

   std::vector<std::shared_ptr<Acts::Surface>> surfaces_out;
   surfaces_out.reserve( surface_list.size() );
   unsigned int counter_i=0;
   for( nlohmann::json a_surface : surface_list) {
      //      unsigned int volume_id = a_volume["volume"].get<unsigned int>();
      /// @return a shared_ptr to a surface object for type polymorphism
      surfaces_out.push_back( surfaceFromJson(a_surface["value"]) );
      if (not surfaces_out.back().get()) {
         uint64_t geo_id = std::numeric_limits<uint64_t>::max();
         uint64_t det_id = std::numeric_limits<uint64_t>::max();
         try {
            geo_id = a_surface["value"]["geo_id"].get<uint64_t>();
         }
         catch (...) {
         }
         try {
            det_id = a_surface["value"]["detID"].get<uint64_t>();
         }
         catch (...) {
         }
         std::cout << "WARNING failed to create surface " << counter_i << " IDs geo: " << geo_id
                   << " detector " << det_id
                   << std::endl;
         surfaces_out.pop_back();
      }
      ++counter_i;
   }
   return surfaces_out;

}

std::shared_ptr<const Acts::TrackingGeometry>
TrackingGeometryJsonReader::createTrackingGeometry(const nlohmann::json& tracking_geometry_description) {
   Acts::ProtoDetector detector(  createProtoDetector(tracking_geometry_description) );
   std::cout << detector.toString("") << std::endl;
   detector.harmonize(true);
   std::vector<std::shared_ptr<Acts::Surface>> surfaces( createSurfaces( tracking_geometry_description) );

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

   // Make the builder
   auto kdtTrackingGeometryBuilder = Acts::KDTreeTrackingGeometryBuilder(
      kdtgConfig, Acts::getDefaultLogger("KDTreeTrackingGeometryBuilder",
                                         m_cfg.toolLogLevel));

   return std::shared_ptr<const Acts::TrackingGeometry>(
      kdtTrackingGeometryBuilder.trackingGeometry(tContext).release());

}
}
