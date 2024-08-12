// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Layer.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <atomic>

namespace Dbg {
   struct InstanceCounter {
      InstanceCounter(std::string &&name) :m_name(std::move(name)) {}
      ~InstanceCounter() {
         std::stringstream msg;
         msg << "DEBUG Layer Surface/Bounds intance counter " << m_name << " :: max surface types: " << m_nSurfaceTypes
             << " max bound types: " << m_nBoundTypes
             << " max bound instances: " << m_nBoundInstances
             << " max bound values: " << m_nBoundValues
             << " n-surfaces: " << (m_nSurfacesN > 0 ? m_nSurfacesSum / m_nSurfacesN : 0)
             << " +- " << (m_nSurfacesN > 1
                           ? sqrt((m_nSurfacesSum2 - m_nSurfacesSum * m_nSurfacesSum / m_nSurfacesN)  / (m_nSurfacesN-1))
                           : 0)
             << " < " << m_nSurfacesMax
             << std::endl;
         std::cout << msg.str() << std::flush;
      }
      void update( std::atomic<unsigned int> &max_val, unsigned int val) {
         for (;;) {
            unsigned int current = max_val;
            if (current >= val || max_val.compare_exchange_weak(current, val)) break;
         }
      }
      void updateNSurfaces(unsigned int n_surfaces) {
         update(m_nSurfacesMax, n_surfaces);
         m_nSurfacesSum += n_surfaces;
         m_nSurfacesSum2 += n_surfaces*n_surfaces;
         ++m_nSurfacesN;
      }
      void updateNSurfaceTypes(unsigned int n_surface_types) { update(m_nSurfaceTypes, n_surface_types); }
      void updateNBoundTypes(unsigned int n_bound_types) { update(m_nBoundTypes, n_bound_types); }
      void updateNBoundInstances(unsigned int n_bound_instances) { update(m_nBoundInstances, n_bound_instances); }
      void updateNBoundValues(unsigned int n_bound_values) { update(m_nBoundValues, n_bound_values); }
      std::atomic<unsigned int> m_nSurfacesMax{};
      std::atomic<unsigned int> m_nSurfacesSum{};
      std::atomic<unsigned int> m_nSurfacesSum2{};
      std::atomic<unsigned int> m_nSurfacesN{};
      std::atomic<unsigned int> m_nSurfaceTypes{};
      std::atomic<unsigned int> m_nBoundTypes{};
      std::atomic<unsigned int> m_nBoundInstances{};
      std::atomic<unsigned int> m_nBoundValues{};
      std::string m_name;
   };

   enum {
      kApproach,
      kSensitive,
      kLayerSurfaces,
      kNSurfaceCategories
   };
   std::array<std::unique_ptr<InstanceCounter>,kNSurfaceCategories> g_InstanceCounter {
      std::make_unique<InstanceCounter>("approach"),
      std::make_unique<InstanceCounter>("sensitive"),
      std::make_unique<InstanceCounter>("layer")
   };
   std::atomic<unsigned int> g_counter{};
   std::atomic<unsigned int> g_callCounter{};
}

Acts::Layer::Layer(std::unique_ptr<SurfaceArray> surfaceArray, double thickness,
                   std::unique_ptr<ApproachDescriptor> ades, LayerType laytyp)
    : m_nextLayers(NextLayers(nullptr, nullptr)),
      m_surfaceArray(surfaceArray.release()),
      m_layerThickness(thickness),
      m_approachDescriptor(nullptr),
      m_representingVolume(nullptr),
      m_layerType(laytyp),
      m_ssRepresentingSurface(1) {
  if (ades) {
    ades->registerLayer(*this);
    m_approachDescriptor = std::move(ades);
    m_ssApproachSurfaces = 1;  // indicates existence
  }
  // indicates existence of sensitive surfaces
  if (m_surfaceArray) {
    m_ssSensitiveSurfaces = 1;
  }
}

const Acts::ApproachDescriptor* Acts::Layer::approachDescriptor() const {
  return m_approachDescriptor.get();
}

Acts::ApproachDescriptor* Acts::Layer::approachDescriptor() {
  return const_cast<ApproachDescriptor*>(m_approachDescriptor.get());
}

void Acts::Layer::closeGeometry(const IMaterialDecorator* materialDecorator,
                                const GeometryIdentifier& layerID,
                                const GeometryIdentifierHook& hook,
                                const Logger& logger) {
  // set the volumeID of this
  assignGeometryId(layerID);
  // assign to the representing surface
  Surface* rSurface = const_cast<Surface*>(&surfaceRepresentation());
  if (materialDecorator != nullptr) {
    materialDecorator->decorate(*rSurface);
  }
  ACTS_DEBUG("layerID: " << layerID);

  rSurface->assignGeometryId(layerID);

  // also find out how the sub structure is defined
  if (surfaceRepresentation().surfaceMaterial() != nullptr) {
    m_ssRepresentingSurface = 2;
  }
  // loop over the approach surfaces
  if (m_approachDescriptor) {
    // indicates the existence of approach surfaces
    m_ssApproachSurfaces = 1;
    // loop through the approachSurfaces and assign unique GeomeryID
    GeometryIdentifier::Value iasurface = 0;
    for (auto& aSurface : m_approachDescriptor->containedSurfaces()) {
      auto asurfaceID = GeometryIdentifier(layerID).setApproach(++iasurface);
      auto mutableASurface = const_cast<Surface*>(aSurface);
      mutableASurface->assignGeometryId(asurfaceID);
      if (materialDecorator != nullptr) {
        materialDecorator->decorate(*mutableASurface);
      }
      // if any of the approach surfaces has material
      if (aSurface->surfaceMaterial() != nullptr) {
        m_ssApproachSurfaces = 2;
      }
    }
  }
  // check if you have sensitive surfaces
  if (m_surfaceArray) {
    // indicates the existence of sensitive surfaces
    m_ssSensitiveSurfaces = 1;
    // loop sensitive surfaces and assign unique GeometryIdentifier
    GeometryIdentifier::Value issurface = 0;
    for (auto& sSurface : m_surfaceArray->surfaces()) {
      auto ssurfaceID = GeometryIdentifier(layerID).setSensitive(++issurface);
      ssurfaceID = hook.decorateIdentifier(ssurfaceID, *sSurface);
      auto mutableSSurface = const_cast<Surface*>(sSurface);
      mutableSSurface->assignGeometryId(ssurfaceID);
      if (materialDecorator != nullptr) {
        materialDecorator->decorate(*mutableSSurface);
      }
      // if any of the sensitive surfaces has material
      if (sSurface->surfaceMaterial() != nullptr) {
        m_ssSensitiveSurfaces = 2;
      }
    }
  }
}

boost::container::small_vector<Acts::SurfaceIntersection, 10>
Acts::Layer::compatibleSurfaces(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const NavigationOptions<Surface>& options) const {
  // the list of valid intersection
  boost::container::small_vector<SurfaceIntersection, 10> sIntersections;

  std::stringstream out;
  out << "DEBUG compatibleSurfaces pos " << position[0] << " " << position[1] << " " << position[2] << " "
      << ", " << direction[0] << " " << direction[1] << " " << direction[2] << " ";
  
  // fast exit - there is nothing to
  if (!m_surfaceArray || !m_approachDescriptor) {
     out << std::endl;
     std::cout << out.str() << std::flush;
    return sIntersections;
  }

  double nearLimit = options.nearLimit;
  double farLimit = options.farLimit;

  // TODO this looks like a major hack; is this really needed?
  // (0) End surface check
  // @todo: - we might be able to skip this by use of options.pathLimit
  // check if you have to stop at the endSurface
  if (options.endObject != nullptr) {
    // intersect the end surface
    // - it is the final one don't use the boundary check at all
    SurfaceIntersection endInter =
        options.endObject
            ->intersect(gctx, position, direction, BoundaryTolerance::None())
            .closest();
    // non-valid intersection with the end surface provided at this layer
    // indicates wrong direction or faulty setup
    // -> do not return compatible surfaces since they may lead you on a wrong
    // navigation path
    out << " intersect w. endObj " << (endInter.isValid() ? " valid  " : " not-intersecting");
    if (endInter.isValid()) {
      farLimit = endInter.pathLength();
    } else {
       out << std::endl;
       std::cout << out.str() << std::flush;
      return sIntersections;
    }
  } else {
    // compatibleSurfaces() should only be called when on the layer,
    // i.e. the maximum path limit is given by the layer thickness times
    // path correction, we take a safety factor of 1.5
    // -> this avoids punch through for cylinders
    double pCorrection =
        surfaceRepresentation().pathCorrection(gctx, position, direction);
    farLimit = 1.5 * thickness() * pCorrection;
  }

  auto isUnique = [&](const SurfaceIntersection& b) {
    auto find_it = std::find_if(
        sIntersections.begin(), sIntersections.end(), [&b](const auto& a) {
          return a.object() == b.object() && a.index() == b.index();
        });
    return find_it == sIntersections.end();
  };

  // lemma 0 : accept the surface
  auto acceptSurface = [&options](const Surface& sf,
                                  bool sensitive = false) -> bool {
    // surface is sensitive and you're asked to resolve
    if (sensitive && options.resolveSensitive) {
      return true;
    }
    // next option: it's a material surface and you want to have it
    if (options.resolveMaterial && sf.surfaceMaterial() != nullptr) {
      return true;
    }
    // last option: resolve all
    return options.resolvePassive;
  };

  out << std::endl;
  // lemma 1 : check and fill the surface
  unsigned int n_surfaces{};
  unsigned int n_surface_types{};
  unsigned int n_bound_types{};
  unsigned int n_bound_instances{};
  unsigned int n_bound_values{};
  std::array<std::size_t,8> surface_types;
  std::array<std::size_t,8> bound_types;
  std::array<const void *,16> bound_instances;
  std::array<std::vector<double>,16> bound_values;
  std::size_t a_geo_id{};
  
  // [&sIntersections, &options, &parameters
  auto processSurface = [&](const Surface& sf, bool sensitive = false) {
    a_geo_id = sf.geometryId().value();
    std::size_t surface_type_hash =  typeid(sf).hash_code();
    ++n_surfaces;
    if (n_surface_types<surface_types.size()) {
    if (std::find(surface_types.begin(),surface_types.begin()+n_surface_types, surface_type_hash)==surface_types.begin()+n_surface_types) {
       surface_types[n_surface_types]=surface_type_hash;
       ++n_surface_types;
    }
    }
    if (n_bound_types<bound_types.size()) {
    std::size_t bound_type_hash =  typeid(sf.bounds()).hash_code();
    if (std::find(bound_types.begin(),bound_types.begin()+n_bound_types, bound_type_hash)==bound_types.begin()+n_bound_types) {
       bound_types[n_bound_types]=bound_type_hash;
       ++n_bound_types;
    }
    }
    if (n_bound_instances<bound_instances.size()) {
    const void *bound_instance =  &sf.bounds();
    if (std::find(bound_instances.begin(),bound_instances.begin()+n_bound_instances, bound_instance)==bound_instances.begin()+n_bound_instances) {
       bound_instances[n_bound_instances]=bound_instance;
       ++n_bound_instances;
    }
    }
    if (n_bound_values<bound_values.size()) {
       std::vector<double> this_bound_values =  sf.bounds().values();
       bool have=false;
       for (unsigned int idx=0; idx < n_bound_values; ++idx) {
          if (this_bound_values.size() == bound_values[idx].size()) {
             bool agree=true;
             for (unsigned int val_i=0; val_i<bound_values[idx].size(); ++val_i) {
                if (bound_values[idx][val_i]!=this_bound_values[val_i]) {
                   agree=false;
                   break;
                }
             }
             if (agree) {
                have=true;
                break;
             }
          }
       }
       if (!have) {
          bound_values[n_bound_values]=std::move(this_bound_values);
          if (Dbg::g_counter < 6000) {
             std::stringstream msg;
             msg << "DEBUG bound_values call " << Dbg::g_callCounter << " n=" << bound_values[n_bound_values].size() << ":";
             for (const double &elm : bound_values[n_bound_values]) {
                msg << std::setw(12) << elm;
             }
             msg<<std::endl;
             std::cout << msg.str() << std::flush;
             ++Dbg::g_counter;
          }
          ++n_bound_values;
       }
    }
    
    // veto if it's start surface
    if (options.startObject == &sf) {
      return;
    }
    // veto if it doesn't fit the prescription
    if (!acceptSurface(sf, sensitive)) {
      return;
    }
    BoundaryTolerance boundaryTolerance = options.boundaryTolerance;
    if (std::find(options.externalSurfaces.begin(),
                  options.externalSurfaces.end(),
                  sf.geometryId()) != options.externalSurfaces.end()) {
      boundaryTolerance = BoundaryTolerance::Infinite();
    }
    // the surface intersection
    //try {
    SurfaceIntersection sfi =
        sf.intersect(gctx, position, direction, boundaryTolerance).closest();
    out << "DEBUG processSurface pos " << position[0] << " " << position[1] << " " << position[2] << " intersect with "
        <<  sf.geometryId().value() << " " << ( sf.associatedDetectorElement() != nullptr ? "D" : "")
       << ( sf.surfaceMaterial() != nullptr ? "M" : "")
        << ( sfi.isValid() ? " intersecting " : " not-intersecting")
        << (detail::checkPathLength(sfi.pathLength(), nearLimit, farLimit) ? "" : " pathlegnth-check-failed")
        << (isUnique(sfi) ? "" : " not unique") << std::endl;
    
    if (sfi.isValid() &&
        detail::checkPathLength(sfi.pathLength(), nearLimit, farLimit) &&
        isUnique(sfi)) {
      sIntersections.push_back(sfi);
    }
    // }
    // catch (...) {
    // SurfaceIntersection sfi =
    //     sf.intersect(gctx, position, direction, BoundaryTolerance::None()).closest();
    // if (sfi.isValid() &&
    //     detail::checkPathLength(sfi.pathLength(), nearLimit, farLimit) &&
    //     isUnique(sfi)) {
    //   sIntersections.push_back(sfi);
    // }
    // }
  };

  // (A) approach descriptor section
  //
  // the approach surfaces are in principle always testSurfaces
  // - the surface on approach is excluded via the veto
  // - the surfaces are only collected if needed
  if (m_approachDescriptor &&
      (options.resolveMaterial || options.resolvePassive)) {
    // the approach surfaces
    const std::vector<const Surface*>& approachSurfaces =
        m_approachDescriptor->containedSurfaces();
    // we loop through and veto
    // - if the approach surface is the parameter surface
    // - if the surface is not compatible with the collect
    for (auto& aSurface : approachSurfaces) {
      processSurface(*aSurface);
    }
    Dbg::g_InstanceCounter[Dbg::kApproach]->updateNSurfaces(n_surfaces);
    Dbg::g_InstanceCounter[Dbg::kApproach]->updateNSurfaceTypes(n_surface_types);
    Dbg::g_InstanceCounter[Dbg::kApproach]->updateNBoundTypes(n_bound_types);
    Dbg::g_InstanceCounter[Dbg::kApproach]->updateNBoundInstances(n_bound_instances);
    Dbg::g_InstanceCounter[Dbg::kApproach]->updateNBoundValues(n_bound_values);

    a_geo_id=0u;
    n_surfaces=0;
    n_surface_types=0;
    n_bound_types=0;
    n_bound_instances=0;
    n_bound_values=0;
  }

  // (B) sensitive surface section
  //
  // check the sensitive surfaces if you have some
  if (m_surfaceArray && (options.resolveMaterial || options.resolvePassive ||
                         options.resolveSensitive)) {
    // get the candidates
    const std::vector<const Surface*>& sensitiveSurfaces =
        m_surfaceArray->neighbors(position);
    // loop through and veto
    // - if the approach surface is the parameter surface
    // - if the surface is not compatible with the type(s) that are collected
    for (auto& sSurface : sensitiveSurfaces) {
      processSurface(*sSurface, true);
    }
    Dbg::g_InstanceCounter[Dbg::kSensitive]->updateNSurfaces(n_surfaces);
    Dbg::g_InstanceCounter[Dbg::kSensitive]->updateNSurfaceTypes(n_surface_types);
    Dbg::g_InstanceCounter[Dbg::kSensitive]->updateNBoundTypes(n_bound_types);
    Dbg::g_InstanceCounter[Dbg::kSensitive]->updateNBoundInstances(n_bound_instances);
    Dbg::g_InstanceCounter[Dbg::kSensitive]->updateNBoundValues(n_bound_values);
    if (n_surfaces>100) {
       std::stringstream msg;
       msg << "DEBUG many surfaces " << n_surfaces << " : " << std::hex << a_geo_id << std::dec << std::endl;
       std::cout << msg.str() << std::flush;
    }

    a_geo_id=0u;
    n_surfaces=0;
    n_surface_types=0;
    n_bound_types=0;
    n_bound_instances=0;
    n_bound_values=0;
  }

  // (C) representing surface section
  //
  // the layer surface itself is a testSurface
  const Surface* layerSurface = &surfaceRepresentation();
  processSurface(*layerSurface);
  Dbg::g_InstanceCounter[Dbg::kLayerSurfaces]->updateNSurfaces(n_surfaces);
  Dbg::g_InstanceCounter[Dbg::kLayerSurfaces]->updateNSurfaceTypes(n_surface_types);
  Dbg::g_InstanceCounter[Dbg::kLayerSurfaces]->updateNBoundTypes(n_bound_types);
  Dbg::g_InstanceCounter[Dbg::kLayerSurfaces]->updateNBoundInstances(n_bound_instances);
  Dbg::g_InstanceCounter[Dbg::kLayerSurfaces]->updateNBoundValues(n_bound_values);
  
  a_geo_id=0u;
  n_surfaces=0;
  n_surface_types=0;
  n_bound_types=0;
  n_bound_instances=0;
  n_bound_values=0;
  ++Dbg::g_callCounter;
 
  out << "DEBUG processSurface pos " << position[0] << " " << position[1] << " " << position[2] << " intersections: ";
  for (const auto &elm : sIntersections) {
     out << ", " << elm.object()->geometryId().value();
  }
  out << std::endl;
  std::cout << out.str() << std::flush;

  return sIntersections;
}


boost::container::small_vector<Acts::SurfaceIntersection, 10>
Acts::Layer::compatibleSurfaces(
    const GeometryContext& gctx,
    std::span<std::pair<Vector3,Vector3> > trajectory_samples,
    const NavigationOptions<Surface>& options) const
{
  // the list of valid intersection
  boost::container::small_vector<SurfaceIntersection, 10> sIntersections;
  if (trajectory_samples.empty()) return sIntersections;

  std::stringstream out;
  out << "DEBUG compatibleSurfaces pos " << trajectory_samples[0].first[0] << " "
      << trajectory_samples[0].first[1] << " " << trajectory_samples[0].first[2] << ", "
      << trajectory_samples[0].second[0] << " "
      << trajectory_samples[0].second[1] << " " << trajectory_samples[0].second[2] << " ";
  // fast exit - there is nothing to
  if (!m_surfaceArray || !m_approachDescriptor) {
     out << std::endl;
     std::cout << out.str() << std::flush;
    return sIntersections;
  }

  double nearLimit = options.nearLimit;
  double farLimit = options.farLimit;

  // TODO this looks like a major hack; is this really needed?
  // (0) End surface check
  // @todo: - we might be able to skip this by use of options.pathLimit
  // check if you have to stop at the endSurface
  if (options.endObject != nullptr) {
    // intersect the end surface
    // - it is the final one don't use the boundary check at all
    SurfaceIntersection endInter =
        options.endObject
            ->intersect(gctx, trajectory_samples[0].first, trajectory_samples[0].second, BoundaryTolerance::None())
            .closest();
    // non-valid intersection with the end surface provided at this layer
    // indicates wrong direction or faulty setup
    // -> do not return compatible surfaces since they may lead you on a wrong
    // navigation path
    out << " intersect w. endObj " << (endInter.isValid() ? " valid  " : " not-intersecting");
    if (endInter.isValid()) {
      farLimit = endInter.pathLength();
    } else {
       out << std::endl;
       std::cout << out.str() << std::flush;
      return sIntersections;
    }
  } else {
    // compatibleSurfaces() should only be called when on the layer,
    // i.e. the maximum path limit is given by the layer thickness times
    // path correction, we take a safety factor of 1.5
    // -> this avoids punch through for cylinders
    double pCorrection =
        surfaceRepresentation().pathCorrection(gctx, trajectory_samples[0].first, trajectory_samples[0].second);
    farLimit = 1.5 * thickness() * pCorrection;
  }

  auto isUnique = [&](const SurfaceIntersection& b) {
    auto find_it = std::find_if(
        sIntersections.begin(), sIntersections.end(), [&b](const auto& a) {
          return a.object() == b.object() && a.index() == b.index();
        });
    return find_it == sIntersections.end();
  };

  // lemma 0 : accept the surface
  auto acceptSurface = [&options](const Surface& sf,
                                  bool sensitive = false) -> bool {
    // surface is sensitive and you're asked to resolve
    if (sensitive && options.resolveSensitive) {
      return true;
    }
    // next option: it's a material surface and you want to have it
    if (options.resolveMaterial && sf.surfaceMaterial() != nullptr) {
      return true;
    }
    // last option: resolve all
    return options.resolvePassive;
  };

  out << std::endl;
  // lemma 1 : check and fill the surface
  unsigned int n_surfaces{};
  unsigned int n_surface_types{};
  unsigned int n_bound_types{};
  unsigned int n_bound_instances{};
  unsigned int n_bound_values{};
  std::array<std::size_t,8> surface_types;
  std::array<std::size_t,8> bound_types;
  std::array<const void *,16> bound_instances;
  std::array<std::vector<double>,16> bound_values;
  std::size_t a_geo_id{};
  
  // [&sIntersections, &options, &parameters
  auto processSurface = [&](const Surface& sf, bool sensitive = false) {
    a_geo_id = sf.geometryId().value();
    std::size_t surface_type_hash =  typeid(sf).hash_code();
    ++n_surfaces;
    if (n_surface_types<surface_types.size()) {
    if (std::find(surface_types.begin(),surface_types.begin()+n_surface_types, surface_type_hash)==surface_types.begin()+n_surface_types) {
       surface_types[n_surface_types]=surface_type_hash;
       ++n_surface_types;
    }
    }
    if (n_bound_types<bound_types.size()) {
    std::size_t bound_type_hash =  typeid(sf.bounds()).hash_code();
    if (std::find(bound_types.begin(),bound_types.begin()+n_bound_types, bound_type_hash)==bound_types.begin()+n_bound_types) {
       bound_types[n_bound_types]=bound_type_hash;
       ++n_bound_types;
    }
    }
    if (n_bound_instances<bound_instances.size()) {
    const void *bound_instance =  &sf.bounds();
    if (std::find(bound_instances.begin(),bound_instances.begin()+n_bound_instances, bound_instance)==bound_instances.begin()+n_bound_instances) {
       bound_instances[n_bound_instances]=bound_instance;
       ++n_bound_instances;
    }
    }
    if (n_bound_values<bound_values.size()) {
       std::vector<double> this_bound_values =  sf.bounds().values();
       bool have=false;
       for (unsigned int idx=0; idx < n_bound_values; ++idx) {
          if (this_bound_values.size() == bound_values[idx].size()) {
             bool agree=true;
             for (unsigned int val_i=0; val_i<bound_values[idx].size(); ++val_i) {
                if (bound_values[idx][val_i]!=this_bound_values[val_i]) {
                   agree=false;
                   break;
                }
             }
             if (agree) {
                have=true;
                break;
             }
          }
       }
       if (!have) {
          bound_values[n_bound_values]=std::move(this_bound_values);
          if (Dbg::g_counter < 6000) {
             std::stringstream msg;
             msg << "DEBUG bound_values call " << Dbg::g_callCounter << " n=" << bound_values[n_bound_values].size() << ":";
             for (const double &elm : bound_values[n_bound_values]) {
                msg << std::setw(12) << elm;
             }
             msg<<std::endl;
             std::cout << msg.str() << std::flush;
             ++Dbg::g_counter;
          }
          ++n_bound_values;
       }
    }
    
    // veto if it's start surface
    if (options.startObject == &sf) {
      return;
    }
    // veto if it doesn't fit the prescription
    if (!acceptSurface(sf, sensitive)) {
      return;
    }
    BoundaryTolerance boundaryTolerance = options.boundaryTolerance;
    if (std::find(options.externalSurfaces.begin(),
                  options.externalSurfaces.end(),
                  sf.geometryId()) != options.externalSurfaces.end()) {
      boundaryTolerance = BoundaryTolerance::Infinite();
    }
    // the surface intersection
    //try {
    unsigned int n_samples_max = sensitive ? trajectory_samples.size() : 1;
    for (unsigned int sample_i=0; sample_i<n_samples_max; ++sample_i) {
    SurfaceIntersection sfi =
       sf.intersect(gctx, trajectory_samples[sample_i].first, trajectory_samples[sample_i].second, boundaryTolerance).closest();
    out << "DEBUG processSurface pos #" << sample_i << " " << trajectory_samples[sample_i].first[0]
        << " " << trajectory_samples[sample_i].first[1] << " " << trajectory_samples[sample_i].first[2] << " intersect with "
        <<  sf.geometryId().value()
        << " " << ( sf.associatedDetectorElement() != nullptr ? "D" : "")
        << ( sf.surfaceMaterial() != nullptr ? "M" : "")
        << ( sfi.isValid() ? " intersecting " : " not-intersecting")
        << (detail::checkPathLength(sfi.pathLength(), nearLimit, farLimit) ? " ok:" : " failed: ")
        << nearLimit << " < " << sfi.pathLength() << " < " << farLimit
        << (isUnique(sfi) ? "" : " not unique") << std::endl;
    
    if (sfi.isValid() &&
        detail::checkPathLength(sfi.pathLength(), nearLimit, farLimit) &&
        isUnique(sfi)) {
      sIntersections.push_back(sfi);
      break;
    }
    }
    // }
    // catch (...) {
    // SurfaceIntersection sfi =
    //     sf.intersect(gctx, position, direction, BoundaryTolerance::None()).closest();
    // if (sfi.isValid() &&
    //     detail::checkPathLength(sfi.pathLength(), nearLimit, farLimit) &&
    //     isUnique(sfi)) {
    //   sIntersections.push_back(sfi);
    // }
    // }
  };

  // (A) approach descriptor section
  //
  // the approach surfaces are in principle always testSurfaces
  // - the surface on approach is excluded via the veto
  // - the surfaces are only collected if needed
  if (m_approachDescriptor &&
      (options.resolveMaterial || options.resolvePassive)) {
    // the approach surfaces
    const std::vector<const Surface*>& approachSurfaces =
        m_approachDescriptor->containedSurfaces();
    // we loop through and veto
    // - if the approach surface is the parameter surface
    // - if the surface is not compatible with the collect
    for (auto& aSurface : approachSurfaces) {
      processSurface(*aSurface);
    }
    Dbg::g_InstanceCounter[Dbg::kApproach]->updateNSurfaces(n_surfaces);
    Dbg::g_InstanceCounter[Dbg::kApproach]->updateNSurfaceTypes(n_surface_types);
    Dbg::g_InstanceCounter[Dbg::kApproach]->updateNBoundTypes(n_bound_types);
    Dbg::g_InstanceCounter[Dbg::kApproach]->updateNBoundInstances(n_bound_instances);
    Dbg::g_InstanceCounter[Dbg::kApproach]->updateNBoundValues(n_bound_values);

    a_geo_id=0u;
    n_surfaces=0;
    n_surface_types=0;
    n_bound_types=0;
    n_bound_instances=0;
    n_bound_values=0;
  }

  // (B) sensitive surface section
  //
  // check the sensitive surfaces if you have some
  if (m_surfaceArray && (options.resolveMaterial || options.resolvePassive ||
                         options.resolveSensitive)) {
    // get the candidates
     //     trajectory_samples[sample_i].first, trajectory_samples[sample_i].second
   const std::vector<const Surface*>& sensitiveSurfaces =
        m_surfaceArray->neighbors(trajectory_samples[0].first);
    // loop through and veto
    // - if the approach surface is the parameter surface
    // - if the surface is not compatible with the type(s) that are collected
    for (auto& sSurface : sensitiveSurfaces) {
      processSurface(*sSurface, true);
    }
    Dbg::g_InstanceCounter[Dbg::kSensitive]->updateNSurfaces(n_surfaces);
    Dbg::g_InstanceCounter[Dbg::kSensitive]->updateNSurfaceTypes(n_surface_types);
    Dbg::g_InstanceCounter[Dbg::kSensitive]->updateNBoundTypes(n_bound_types);
    Dbg::g_InstanceCounter[Dbg::kSensitive]->updateNBoundInstances(n_bound_instances);
    Dbg::g_InstanceCounter[Dbg::kSensitive]->updateNBoundValues(n_bound_values);
    if (n_surfaces>100) {
       std::stringstream msg;
       msg << "DEBUG many surfaces " << n_surfaces << " : " << std::hex << a_geo_id << std::dec << std::endl;
       std::cout << msg.str() << std::flush;
    }

    a_geo_id=0u;
    n_surfaces=0;
    n_surface_types=0;
    n_bound_types=0;
    n_bound_instances=0;
    n_bound_values=0;
  }

  // (C) representing surface section
  //
  // the layer surface itself is a testSurface
  const Surface* layerSurface = &surfaceRepresentation();
  processSurface(*layerSurface);
  Dbg::g_InstanceCounter[Dbg::kLayerSurfaces]->updateNSurfaces(n_surfaces);
  Dbg::g_InstanceCounter[Dbg::kLayerSurfaces]->updateNSurfaceTypes(n_surface_types);
  Dbg::g_InstanceCounter[Dbg::kLayerSurfaces]->updateNBoundTypes(n_bound_types);
  Dbg::g_InstanceCounter[Dbg::kLayerSurfaces]->updateNBoundInstances(n_bound_instances);
  Dbg::g_InstanceCounter[Dbg::kLayerSurfaces]->updateNBoundValues(n_bound_values);
  
  a_geo_id=0u;
  n_surfaces=0;
  n_surface_types=0;
  n_bound_types=0;
  n_bound_instances=0;
  n_bound_values=0;
  ++Dbg::g_callCounter;
 
  out << "DEBUG processSurface pos " << trajectory_samples[0].first[0] << " "
      << trajectory_samples[0].first[1] << " " << trajectory_samples[0].first[2] << " intersections: ";
  for (const auto &elm : sIntersections) {
     out << ", " << elm.object()->geometryId().value()
        << " " << ( elm.object()->associatedDetectorElement() != nullptr ? "D" : "")
        << ( elm.object()->surfaceMaterial() != nullptr ? "M" : "")
        ;
  }
  out << std::endl;
  std::cout << out.str() << std::flush;

  return sIntersections;
}


Acts::SurfaceIntersection Acts::Layer::surfaceOnApproach(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const NavigationOptions<Layer>& options) const {
  // resolve directive based by options
  // - options.resolvePassive is on -> always
  // - options.resolveSensitive is on -> always
  // - options.resolveMaterial is on
  //   && either sensitive or approach surfaces have material
  bool resolvePS = options.resolveSensitive || options.resolvePassive;
  bool resolveMS = options.resolveMaterial &&
                   (m_ssSensitiveSurfaces > 1 || m_ssApproachSurfaces > 1 ||
                    (surfaceRepresentation().surfaceMaterial() != nullptr));

  // The Limits
  double nearLimit = options.nearLimit;
  double farLimit = options.farLimit;

  // Helper function to find valid intersection
  auto findValidIntersection =
      [&](const SurfaceMultiIntersection& sfmi) -> SurfaceIntersection {
    for (const auto& sfi : sfmi.split()) {
      if (sfi.isValid() &&
          detail::checkPathLength(sfi.pathLength(), nearLimit, farLimit)) {
        return sfi;
      }
    }

    // Return an invalid one
    return SurfaceIntersection::invalid();
  };

  // Approach descriptor present and resolving is necessary
  if (m_approachDescriptor && (resolvePS || resolveMS)) {
    SurfaceIntersection aSurface = m_approachDescriptor->approachSurface(
        gctx, position, direction, options.boundaryTolerance, nearLimit,
        farLimit);
    return aSurface;
  }

  // Intersect and check the representing surface
  const Surface& rSurface = surfaceRepresentation();
  auto sIntersection =
      rSurface.intersect(gctx, position, direction, options.boundaryTolerance);
  return findValidIntersection(sIntersection);
}
