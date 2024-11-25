#include "Acts/Geometry/GeoIdRegistry.hpp"
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>

namespace Dbg {

   std::mutex Statistics::s_mutex;
   std::unique_ptr<Statistics> Statistics::s_instance;
   Statistics::~Statistics() {
      const std::array<std::string, kNCounter> counter_name {
         std::string("BoundCovChanged"),
         std::string("FullTransportJacobianChanged"),
         std::string("FreeTransportJacobianChanged"),
         std::string("FreeToPathDeritivesChanged"),
         std::string("BoundToFreeJacobianChanged"),
         std::string("Sampling"),
         std::string("NotSampling"),
         std::string("SamplingSteppingHelper"),
         std::string("NotSamplingSteppingHelper"),
         std::string("InvalidCovPos"),
         std::string("InvalidFreeCovPos"),
         std::string("InvalidCovNeg"),
         std::string("InvalidFreeCovNeg"),
         std::string("InvalidCovPosSteppingHelper"),
         std::string("InvalidFreeCovPosSteppingHelper"),
         std::string("InvalidCovNegSteppingHelper"),
         std::string("InvalidFreeCovNegSteppingHelper"),
         std::string("CovPos"),
         std::string("CovNeg"),
         std::string("CovPosSteppingHelper"),
         std::string("CovNegSteppingHelper"),
         std::string("CovPosCKF"),
         std::string("CovNegCKF"),
         std::string("InvalidCovPosCKF"),
         std::string("InvalidCovNegCKF"),
         std::string("Updates"),
         std::string("UpdatesValidInput"),
         std::string("UpdatesValidInputAndOutput")
      };
      unsigned int idx=0;
      std::stringstream out;
      for (const auto &elm : m_counter) {
         out << "DEBUG Statistics counter " << std::setw(9) << elm << " " << counter_name.at(idx) << std::endl;
         ++idx;
      }
      std::cout << out.str() << std::flush;
   }
   
   std::mutex GeoIdRegistry::s_instanceMutex;
   std::unique_ptr<GeoIdRegistry> GeoIdRegistry::s_instance;
   std::unordered_map<unsigned int, std::unordered_set< std::size_t> > GeoIdRegistry::s_dummyMap;
   std::unordered_set< std::size_t> GeoIdRegistry::s_dummySet;
   thread_local unsigned int GeoIdRegistry::s_currentEventId;
   thread_local unsigned int GeoIdRegistry::s_currentObjectId;
   
}
   
