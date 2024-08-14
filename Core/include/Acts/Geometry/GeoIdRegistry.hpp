#pragma once

#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <memory>
#include <limits>
#include <vector>
namespace Dbg {
class GeoIdRegistry {
public:
   void record(unsigned int event_id, unsigned int id, std::size_t geometry_id_value) {
      std::lock_guard<std::mutex> lock(m_mutex);
      m_geoIds[event_id][id].insert(geometry_id_value);
   }
   void record(unsigned int event_id, unsigned int id, const std::vector<std::size_t> &geometry_ids) {
      std::lock_guard<std::mutex> lock(m_mutex);
      for (auto geometry_id_value : geometry_ids) {
         m_geoIds[event_id][id].insert(geometry_id_value);
      }
   }
   void releaseEvent(unsigned int event_id) {
      std::lock_guard<std::mutex> lock(m_mutex);
      auto iter = m_geoIds.find(event_id);
      if (iter != m_geoIds.end()) {
         m_geoIds.erase(iter);
      }
   }
   const std::unordered_map<unsigned int, std::unordered_set< std::size_t> > &
   getEventGeoIds(unsigned int event_id=GeoIdRegistry::currentEventId()) const {
      std::lock_guard<std::mutex> lock(m_mutex);
      auto iter = m_geoIds.find(event_id);
      if (iter != m_geoIds.end()) {
         return iter->second;
      }
      return s_dummyMap;
   }

   const std::unordered_set< std::size_t> &
   getObjectGeoIds(unsigned int object_id=GeoIdRegistry::currentObjectId(),
                   unsigned int event_id=GeoIdRegistry::currentEventId()
                   ) const {
      std::lock_guard<std::mutex> lock(m_mutex);
      auto iter = m_geoIds.find(event_id);
      if (iter != m_geoIds.end()) {
         auto obj_iter = iter->second.find(object_id);
         if (obj_iter !=iter->second.end()) {
            return obj_iter->second;
         }
      }
      return s_dummySet;
   }
   const std::unordered_map <unsigned int, std::unordered_map<unsigned int, std::unordered_set< std::size_t> > > &
   getGeoIds() const { return m_geoIds; }

   static void setCurrentEventId( unsigned int event_id ) {
      s_currentEventId = event_id;
   }
   static unsigned int currentEventId() {
      return s_currentEventId;
   }
   static void setCurrentObjectId( unsigned int object_id ) {
      s_currentObjectId = object_id;
   }
   static unsigned int currentObjectId() {
      return s_currentObjectId;
   }

   static constexpr unsigned int invalidId() {
      return std::numeric_limits<unsigned int>::max();
   }
   static constexpr bool isValid(unsigned int id) {
      return id != invalidId();
   }

   static GeoIdRegistry &instance() {
      if (!s_instance) {
         std::lock_guard<std::mutex> lock(s_instanceMutex);
         if (!s_instance) {
            s_instance = std::make_unique<GeoIdRegistry>();
            // s_currentEventId = invalidId();
            // s_currentObjectId = invalidId();
         }
      }
      return *s_instance;
   }

   
private:
   thread_local static unsigned int s_currentEventId;
   thread_local static unsigned int s_currentObjectId;

   mutable std::mutex m_mutex;
   std::unordered_map <unsigned int, std::unordered_map<unsigned int, std::unordered_set< std::size_t> > > m_geoIds;
   static std::mutex s_instanceMutex;
   static std::unique_ptr<GeoIdRegistry> s_instance;
   static std::unordered_map<unsigned int, std::unordered_set< std::size_t> > s_dummyMap;
   static std::unordered_set< std::size_t> s_dummySet;
};

class EventIdGuard {
public:
   EventIdGuard(unsigned int event_id)
      : m_orig(GeoIdRegistry::currentEventId()),
        m_origObjId(GeoIdRegistry::currentObjectId())
   {
      GeoIdRegistry::setCurrentEventId( event_id);
      GeoIdRegistry::setCurrentObjectId( GeoIdRegistry::invalidId());
   }
   ~EventIdGuard() {
      GeoIdRegistry::setCurrentObjectId( m_origObjId);
      GeoIdRegistry::setCurrentEventId( m_orig);
   }
   unsigned int m_orig;
   unsigned int m_origObjId;
};
class ObjectIdGuard {
public:
   ObjectIdGuard(unsigned int object_id) : m_orig(GeoIdRegistry::currentObjectId() ){
      GeoIdRegistry::setCurrentObjectId( object_id);
   }
   ~ObjectIdGuard() {
      GeoIdRegistry::setCurrentObjectId( m_orig);
   }
   unsigned int m_orig;
};

   class GeoIdHelper {
   public:
      GeoIdHelper(unsigned int object_id=GeoIdRegistry::currentObjectId(),
                  unsigned int event_id=GeoIdRegistry::currentEventId()
                   )
      : m_eventId(event_id),
        m_objectId(object_id)
      {
         if (GeoIdRegistry::isValid(event_id)) {
            if (GeoIdRegistry::isValid(object_id)) {
               m_objectGeoIds = &GeoIdRegistry::instance().getObjectGeoIds(object_id, event_id);
            }
            else {
               m_eventGeoIds = &GeoIdRegistry::instance().getEventGeoIds(event_id);
            }
         }
      }

      bool isKnown(std::size_t geo_id) const {
         if (m_objectGeoIds) {
            return m_objectGeoIds->find(geo_id) != m_objectGeoIds->end();
         }
         else if (m_eventGeoIds) {
            for (const std::pair<const unsigned int,  std::unordered_set< std::size_t> >  &object_geo_ids : *m_eventGeoIds) {
               if (object_geo_ids.second.find(geo_id) != object_geo_ids.second.end()) {
                  return true;
               }
            }
         }
         return false;
      }

      unsigned int matchingObjectId(std::size_t geo_id) const {
         if (m_objectGeoIds) {
            if (m_objectGeoIds->find(geo_id) != m_objectGeoIds->end()) {
               return m_objectId;
            }
         }
         else if (m_eventGeoIds) {
            for (const std::pair<const unsigned int,  std::unordered_set< std::size_t> >  &object_geo_ids : *m_eventGeoIds) {
               if (object_geo_ids.second.find(geo_id) != object_geo_ids.second.end()) {
                  return object_geo_ids.first;
               }
            }
         }
         return GeoIdRegistry::invalidId();
      }

      unsigned int eventId() const {return m_eventId; }
      unsigned int objectId() const { return m_objectId; }
   private:
      unsigned int m_eventId = GeoIdRegistry::invalidId();
      unsigned int m_objectId = GeoIdRegistry::invalidId();
      const std::unordered_map<unsigned int, std::unordered_set< std::size_t> > *m_eventGeoIds = nullptr;
      const std::unordered_set< std::size_t> *m_objectGeoIds=nullptr;
   };
}
