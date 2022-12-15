#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Plugins/Json/TrackingGeometryJsonReader.hpp"
#include "Acts/Plugins/Json/VolumeBoundsJsonConverter.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"

#include <memory>
#include <array>
#include <deque>
#include <limits>
#include <stdexcept>

namespace dbg {
   void dumpMap(std::unordered_map<unsigned int, Acts::TrackingGeometryJsonReader::VolumeInfo > &volume_name_map) {
      for (const std::pair<const unsigned int, Acts::TrackingGeometryJsonReader::VolumeInfo > &elm : volume_name_map) {
         std::cout << elm.first << " -> " << elm.second.proto_volume.name << " " <<elm.second.id
                   << " childs  ";
         for (unsigned int a_child :  elm.second.childs) {
            std::cout << " " << a_child;
         }

         std::cout  << std::endl;
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

std::unique_ptr<Acts::TrackingGeometry> TrackingGeometryJsonReader::trackingGeometry(const nlohmann::json& tracking_geometry_description) {
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
         //         std::cout << "ERROR not a z-axis aligned tracking volume " << element.proto_volume.name << std::endl;
         is_z_axis_aligned = false;
      }

      const Acts::CylinderVolumeBounds *cylinder_bounds=dynamic_cast<const Acts::CylinderVolumeBounds *>(element.volumeBounds.get());

      if (not element.volumeBounds or element.volumeBounds->type() != Acts::VolumeBounds::eCylinder or not cylinder_bounds) {
         //         std::cout << "ERROR not a cylinder tracking volume " << element.proto_volume.name << std::endl;
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
   if (not is_cylinder or not is_z_axis_aligned) return std::unique_ptr<Acts::TrackingGeometry>{};

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
            for(unsigned int domain_i =0 ; domain_i < a_child_element.domains.size(); ++domain_i) {
               assert( a_child_element.domain.size() == element.second.domains.size());
               if (a_child_element.domains.at(domain_i) != element.second.domains.at(domain_i)) {
                  non_identical_domain_mask |= (1<<domain_i);
               }
               for(unsigned int other_child_i=0; other_child_i<child_i; ++other_child_i) {
                  VolumeInfo &other_child_element = hierarchy.at( element.second.childs.at(other_child_i) );
                  if (a_child_element.domains.at(domain_i).overlaps( other_child_element.domains.at(domain_i) )) {
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
         // set domain if domain is identical to the domain definition of all childs;
         element.second.setDomainMask |= same_domain_mask;
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

   static std::array<Acts::BinningValue, VolumeInfo::eNDomains> acts_domain_type
      {Acts::binR, Acts::binZ}; // must match order in VolumeInfos::EDomainTypes
   while (!element_queue.empty()) {
      unsigned int element_idx = element_queue.front();
      element_queue.pop_front();
      VolumeInfo &an_element = hierarchy.at(element_idx);

      bool childs_processed=true;
      for(unsigned int child_i=0; child_i<an_element.childs.size(); ++child_i) {
         VolumeInfo &a_child_element = hierarchy.at(an_element.childs.at(child_i) );
         childs_processed &= a_child_element.processed;
         if (not an_element.processed) { break; }
      }
      if (childs_processed ) {
         an_element.proto_volume.constituentVolumes.reserve(an_element.childs.size());
         //         std::cout << "DEBUG " << an_element.proto_volume.name << " childs:";
         // @TODO sort childs by z or r
         for(unsigned int child_i=0; child_i<an_element.childs.size(); ++child_i) {
            an_element.proto_volume.constituentVolumes.push_back( hierarchy.at(an_element.childs.at(child_i) ).proto_volume);
            //            std::cout << " " << hierarchy.at(an_element.childs.at(child_i) ).proto_volume.name;
         }
         //         std::cout << std::endl;
         for(unsigned int domain_i =0 ; domain_i < an_element.domains.size(); ++domain_i) {
            if ( (an_element.setDomainMask & (1<<domain_i)) ) {
               an_element.proto_volume.extent.set(acts_domain_type.at(domain_i),
                                                   an_element.domains.at(domain_i).min,
                                                   an_element.domains.at(domain_i).max);
               // std::cout << "DEBUG " << an_element.proto_volume.name << " set domain "
               //           << domain_i << " (" << acts_domain_type.at(domain_i) << ")"
               //           << " " << an_element.domains.at(domain_i).min
               //           << " .. " << an_element.domains.at(domain_i).max
               //           << std::endl;
            }
            if (!an_element.childs.empty()) {
               if (    (an_element.noOverlapMask & (1<<domain_i) ) ) {
                  // std::cout << "DEBUG " << an_element.proto_volume.name << " set constituentDomain "
                  //           << " " << domain_i << "(" << acts_domain_type.at(domain_i) << ")"
                  //           << std::endl;

                  an_element.proto_volume.constituentBinning = {
                     Acts::BinningData(Acts::open,
                                      acts_domain_type.at(domain_i),
                                       {0., 1.})};

               }
            }
            an_element.processed = childs_processed;
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
         std::cout << "DEBUG unprocessed volume " << element.second.proto_volume.name << " childs:";
         for(unsigned int child_i=0; child_i<element.second.childs.size(); ++child_i) {
            std::cout << " " << hierarchy.at(element.second.childs.at(child_i) ).proto_volume.name;
         }
         std::cout << std::endl;
      }
   }

   return std::unique_ptr<Acts::TrackingGeometry>{};
}
}
