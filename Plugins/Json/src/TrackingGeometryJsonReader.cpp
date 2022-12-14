#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Plugins/Json/TrackingGeometryJsonReader.hpp"


#include <memory>
#include <array>
#include <deque>
#include <limits>
#include <stdexcept>

namespace dbg {
   void dumpMap(std::unordered_map<unsigned int, Acts::TrackingGeometryJsonReader::Element > &volume_name_map) {
      for (const std::pair<const unsigned int, Acts::TrackingGeometryJsonReader::Element > &elm : volume_name_map) {
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
TrackingGeometryJsonReader::Element &
TrackingGeometryJsonReader::registerElement(unsigned int volume_id,
                                            const std::string &name,
                                            std::unordered_map<unsigned int, Element > &volume_map,
                                            std::unordered_map<std::string, unsigned int > &name_map) {
   if (name.empty()) {
      throw std::runtime_error("volume without name");
   }
   
   std::pair<std::unordered_map<std::string, unsigned int>::iterator, bool>
      insert_key = name_map.insert( std::make_pair(name, name_map.size()));
   unsigned int new_key = insert_key.first->second;
   std::pair<std::unordered_map<unsigned int, Element>::iterator, bool>
      insert_result = volume_map.insert( std::make_pair(new_key, Element(name, volume_id)));
   Element &new_element = insert_result.first->second;

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
               std::cout << "DEBUG child " << child_key << " <" << child_name << ">" << std::endl;
               std::pair<std::unordered_map<unsigned int, Element>::iterator, bool>
                  child_insert_result = volume_map.insert( std::make_pair(child_key,Element(child_name)));

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
   static const std::array<std::string, Element::kNBinningVariables> binning_variables {
      std::string("binR"),
      std::string("binPhi"),std::string("binZ")};
   static const std::array<std::string, Element::Binning::kNRangeTypes> binning_range_types {
      std::string("open"),
      std::string("closed")};
   std::unordered_map<std::string, unsigned int > name_map;
   std::unordered_map<unsigned int, Element> hierarchy;
   nlohmann::json volume_list  = volumes.at("entries");;
   for( nlohmann::json a_volume : volume_list) {
      unsigned int volume_id = a_volume["volume"].get<unsigned int>();
      nlohmann::json value=a_volume["value"];
      Element &element = registerElement(volume_id, value["NAME"].get<std::string>(), hierarchy, name_map);
      for (nlohmann::json &a_binning : value["material"]["binUtility"]["binningdata"]) {
         element.binning.at(getNamedEnum<unsigned short>(binning_variables,
                                                         a_binning["value"].get<std::string>(),
                                                         "Unknwown binning: "))
            = Element::Binning(a_binning["min"].get<float>(),
                               a_binning["max"].get<float>(),
                               getNamedEnum<Element::Binning::EType>(binning_range_types,
                                                                     a_binning["option"].get<std::string>(),
                                                                     "Unknown binning option: " ));
      }
   }
   // for( nlohmann::json a_volume : volume_list) {
   //    unsigned int volume_id = a_volume["volume"].get<unsigned int>();
   //    nlohmann::json value=a_volume["value"];
   //    (void) registerElement(volume_id, value["NAME"].get<std::string>(), hierarchy, true);
   // }

   std::deque<unsigned int> element_queue;
   Element *world=nullptr;
   for(std::pair<const unsigned int, Element> element : hierarchy) {
      if (element.second.parent == std::numeric_limits<unsigned int>::max()) {
         world=&element.second;
      }
      if (element.second.childs.empty()) {
         element_queue.push_back(element.first);
      }
      else {
         unsigned int non_identical_binning_mask=0;
         unsigned int child_binning_overlaps=0;
         for(unsigned int child_i=0; child_i<element.second.childs.size(); ++child_i) {
            unsigned int child_idx = element.second.childs.at(child_i);
            Element &a_child_element = hierarchy.at(child_idx);
            for(unsigned int binning_i =0 ; binning_i < a_child_element.binning.size(); ++binning_i) {
               assert( a_child_element.binning.size() == element.second.binning.size());
               if (a_child_element.binning.at(binning_i) != element.second.binning.at(binning_i)) {
                  non_identical_binning_mask |= (1<<binning_i);
               }
               for(unsigned int other_child_i=0; other_child_i<child_i; ++other_child_i) {
                  Element &other_child_element = hierarchy.at( element.second.childs.at(other_child_i) );
                  if (a_child_element.binning.at(binning_i).overlaps( other_child_element.binning.at(binning_i) )) {
                     child_binning_overlaps |= (1<<binning_i);
                     break;
                  }
               }
            }
         }
         unsigned int non_overlapping_binnings = 0;
         for(unsigned int binning_i =0 ; binning_i < element.second.binning.size(); ++binning_i) {
            non_overlapping_binnings |= (1<<binning_i);
         }
         unsigned int same_binning_mask = non_overlapping_binnings & (~non_identical_binning_mask);
         non_overlapping_binnings &= ~child_binning_overlaps;
         
         for(unsigned int child_i=0; child_i<element.second.childs.size(); ++child_i) {
            Element &a_child_element = hierarchy.at( element.second.childs.at(child_i) );
            if (a_child_element.childs.empty()) {
               a_child_element.setBinningMask |= non_overlapping_binnings;
            }
         }
         // set binning if binning is identical to the binning definition of all childs;
         element.second.setBinningMask |= same_binning_mask;
         element.second.noOverlapMask  |= non_overlapping_binnings;
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

   static std::array<Acts::BinningValue, Element::kNBinningVariables> acts_binning_type
      {Acts::binR, Acts::binPhi,Acts::binZ}; // must match order in Elements::EBinnningVariable
   static std::array<Acts::BinningOption, Element::Binning::kNRangeTypes> acts_range_type
   {Acts::open, Acts::closed}; // must match order in Elements::EBinnningVariable
   while (!element_queue.empty()) {
      unsigned int element_idx = element_queue.front();
      element_queue.pop_front();
      Element &an_element = hierarchy.at(element_idx);

      bool processed=true;
      for(unsigned int child_i=0; child_i<an_element.childs.size(); ++child_i) {
         Element &a_child_element = hierarchy.at(an_element.childs.at(child_i) );
         processed &= a_child_element.processed;
         if (not an_element.processed) { break; }
      }
      an_element.processed = processed;
      if (processed ) {
         an_element.proto_volume.constituentVolumes.reserve(an_element.childs.size());
         std::cout << "DEBUG " << an_element.proto_volume.name << " childs:";
         for(unsigned int child_i=0; child_i<an_element.childs.size(); ++child_i) {
            an_element.proto_volume.constituentVolumes.push_back( hierarchy.at(an_element.childs.at(child_i) ).proto_volume);
            std::cout << " " << hierarchy.at(an_element.childs.at(child_i) ).proto_volume.name;
         }
         std::cout << std::endl;
         for(unsigned int binning_i =0 ; binning_i < an_element.binning.size(); ++binning_i) {
            if ((an_element.setBinningMask & (1<<binning_i))
                && (binning_i != Element::kPhi || !an_element.binning.at(binning_i).isClosed())) {
               an_element.proto_volume.extent.set(acts_binning_type.at(binning_i),
                                                  an_element.binning.at(binning_i).min,
                                                  an_element.binning.at(binning_i).max);
               std::cout << "DEBUG " << an_element.proto_volume.name << " set extend "
                         << binning_i << " (" << acts_binning_type.at(binning_i) << ")"
                         << " " << an_element.binning.at(binning_i).min
                         << " .. " << an_element.binning.at(binning_i).max
                         << std::endl;
            }
            if (!an_element.childs.empty()) {
               if ((an_element.noOverlapMask & (1<<binning_i))
                   && (binning_i != Element::kPhi || !an_element.binning.at(binning_i).isClosed())) {
                  std::cout << "DEBUG " << an_element.proto_volume.name << " set constituentBinning "
                            << static_cast<unsigned int>(acts_range_type.at( an_element.binning.at(binning_i).type))
                            << " " << binning_i << "(" << acts_binning_type.at(binning_i) << ")"
                            << std::endl;

                  an_element.proto_volume.constituentBinning = {
                     Acts::BinningData(acts_range_type.at( an_element.binning.at(binning_i).type),
                                       acts_binning_type.at(binning_i),
                                       {0., 1.})};

               }
            }
         }
         element_queue.push_back(an_element.parent);
      }
      else {
         element_queue.push_back(element_idx);
      }
   }

   return std::unique_ptr<Acts::TrackingGeometry>{};
}
}
