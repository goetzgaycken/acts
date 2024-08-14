#include "Acts/Geometry/GeoIdRegistry.hpp"

namespace Dbg {
   std::mutex GeoIdRegistry::s_instanceMutex;
   std::unique_ptr<GeoIdRegistry> GeoIdRegistry::s_instance;
   std::unordered_map<unsigned int, std::unordered_set< std::size_t> > GeoIdRegistry::s_dummyMap;
   std::unordered_set< std::size_t> GeoIdRegistry::s_dummySet;
   thread_local unsigned int GeoIdRegistry::s_currentEventId;
   thread_local unsigned int GeoIdRegistry::s_currentObjectId;
   
}
   
