// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/IWriter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

#include <cstdint>
#include <mutex>
#include <string>
#include <cmath>
#include <vector>
#include <type_traits>

class TFile;
//class TTree;
//class TBranch;

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TClass.h"
#include "TVirtualStreamerInfo.h"

namespace ActsExamples {

/// Write out particles as a flat TTree.
///
/// Each entry in the TTree corresponds to one particle for optimum writing
/// speed. The event number is part of the written data.
///
/// Safe to use from multiple writer threads. To avoid thread-saftey issues,
/// the writer must be the sole owner of the underlying file. Thus, the
/// output file pointer can not be given from the outside.
class RootDigiBase {
 public:
  struct ConfigBase {
    /// truth particle collection to write.
    std::string particles = "particles_final";
    /// name of the measurements collection
    std::string measurements = "measurements";
    /// name for the map between measurements and truth particles
    std::string measurementParticlesMap = "measurement_particles_map";
    /// Path to the output file.
    std::string filePath;
    /// Output file access mode.
    std::string fileMode;
    /// Name of the tree within the output file.
    std::string treeName = "Digi";

    /// truth particle branch prefix
    std::string truthParticlePrefix;
    /// truth vertex branch prefix
    std::string truthVertexPrefix;
    /// truth vertex branch prefix
    std::string eventInfoPrefix;
    /// truth  branch prefix
    std::vector<std::string> measurementPrefix;
    /// for each measurement prefix the type of the measurement : 0 pixel, 1 strips 
    std::vector<int> measurementType;

    /// the tracking geometry necessary to match measurements to surfaces and/or perform local to global transformations
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };
  RootDigiBase() = default;
  /// Construct the particle writer.
  ///
  /// @params cfg is the configuration object
  /// @params lvl is the loggin\g level
  RootDigiBase(ActsExamples::RootDigiBase::ConfigBase cfg,
               const std::string &name,
               bool write,
               Acts::Logging::Level lvl);

  /// Ensure underlying file is closed.
  virtual ~RootDigiBase();

  /// End-of-run hook
  ProcessCode endRunBase();

  /// Get readonly access to the config parameters
  const ConfigBase& config() const { return m_cfg; }

 protected:
  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] particles are the particle to be written

   
  ConfigBase m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  std::mutex m_mutex;
  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;
  std::vector<TBranch *> m_activeBranches;
  /// Event identifier.
  struct EventInfo {
     //     UInt_t          runNumber;
     uint64_t       eventNumber;
     //     UInt_t          lumiBlock;
     //     UInt_t          timeStamp;
     //     UInt_t          timeStampNSOffset;
     //     UInt_t          bcid;
  };
  // struct TruthEvent {
  //    std::vector<vector<float> > weights;
  //    std::vector<float>          crossSection;
  //    std::vector<float>          crossSectionError;
  // };
  using Int_t = int32_t;
  using UInt_t = uint32_t;
   
  struct ParticleContainer {
     static constexpr unsigned int kNParticlesMax = 15000;
     std::vector<int>     pdgId;
     std::vector<int>     barcode; // @TODO should be uint64 (schema evolution? )
     std::vector<int>     status;
     Int_t           n_prodVtxLink = 0;
     // UInt_t          TruthParticlesAux_prodVtxLink_m_persKey[1];   //[TruthParticlesAux.prodVtxLink_]
     UInt_t          prodVtxLink_m_persIndex[kNParticlesMax];   //[TruthParticlesAux.prodVtxLink_]
     Int_t           n_decayVtxLink = 0;
     //UInt_t          decayVtxLink_m_persKey[kMaxTruthParticlesAux_decayVtxLink];   //[TruthParticlesAux.decayVtxLink_]
     UInt_t          decayVtxLink_m_persIndex[kNParticlesMax];   //[TruthParticlesAux.decayVtxLink_]
     std::vector<float>   px;
     std::vector<float>   py;
     std::vector<float>   pz;
     std::vector<float>   e;
     std::vector<float>   m;
     void reserve(std::size_t new_capacity) {
        pdgId.reserve(new_capacity);
        barcode.reserve(new_capacity);
        status.reserve(new_capacity);
        if (new_capacity>kNParticlesMax) throw std::runtime_error("Exceeded Maximum number of truth particles");
        px.reserve(new_capacity);
        py.reserve(new_capacity);
        pz.reserve(new_capacity);
        e.reserve(new_capacity);
        m.reserve(new_capacity);
     }
     void push_back( int a_pdg_id, int a_barcode, int a_status, int a_prod_vtx_idx, int a_decay_vtx_idx,
                     float a_px, float a_py, float a_pz, float a_e, float a_m) {
        pdgId.push_back(a_pdg_id);
        barcode.push_back(a_barcode);
        status.push_back(a_status);
        assert( n_prodVtxLink < kNParticlesMax);
        prodVtxLink_m_persIndex[++n_prodVtxLink]=a_prod_vtx_idx;
        assert( n_decayVtxLink < kNParticlesMax);
        decayVtxLink_m_persIndex[++n_decayVtxLink]=a_decay_vtx_idx;
        px.push_back(a_px);
        py.push_back(a_py);
        pz.push_back(a_pz);
        e.push_back(a_e);
        m.push_back(a_m);
     }
     void push_back( int a_pdg_id, int a_barcode, int a_status, int a_prod_vtx_idx, int a_decay_vtx_idx,
                     float a_px, float a_py, float a_pz, float a_m) {
        push_back(a_pdg_id, a_barcode, a_status, a_prod_vtx_idx, a_decay_vtx_idx,a_px,a_py,a_pz,a_m,
                  static_cast<float>(std::sqrt( std::pow(a_px,2) + std::pow(a_py,2) + std::pow(a_pz,2) + std::pow(a_m,2) )));
     }
     void clear() {
        pdgId.clear();
        barcode.clear();
        status.clear();
        n_prodVtxLink=0;
        n_decayVtxLink=0;
        px.clear();
        py.clear();
        pz.clear();
        e.clear();
        m.clear();
     }
  };
  struct ParticleVertexContainer {
     //     std::vector<int>     id;
     std::vector<int>     barcode; 
     //     std::vector<std::vector<ElementLink<DataStd::Vector<xAOD::TruthParticle_v1> > > > TruthVerticesAux_incomingParticleLinks;
     //     std::vector<std::vector<ElementLink<DataStd::Vector<xAOD::TruthParticle_v1> > > > TruthVerticesAux_outgoingParticleLinks;
     std::vector<float>   x;
     std::vector<float>   y;
     std::vector<float>   z;
     std::vector<float>   t;
     void reserve(std::size_t new_capacity) {
        barcode.reserve(new_capacity);
        x.reserve(new_capacity);
        y.reserve(new_capacity);
        z.reserve(new_capacity);
        t.reserve(new_capacity);
     }
     void push_back(int a_barcode, float a_x, float a_y, float a_z, float a_t) {
        barcode.push_back(a_barcode);
        x.push_back(a_x);
        y.push_back(a_y);
        z.push_back(a_z);
        t.push_back(a_t);
     }
     void clear() {
        barcode.clear();
        x.clear();
        y.clear();
        z.clear();
     }
  };
  struct ClusterContainer {
     std::vector<float>    localX;
     std::vector<float>    localY;
     std::vector<float>    localXError;
     std::vector<float>    localYError;
     std::vector<float>    localXYCorrelation;
     std::vector<float>    *globalX = nullptr;
     std::vector<float>    globalY;
     std::vector<float>    globalZ;
     std::vector<unsigned long> detectorElementID;
     void reserve(std::size_t new_capacity) {
        localX.reserve(new_capacity);
        localY.reserve(new_capacity);
        localXError.reserve(new_capacity);
        localYError.reserve(new_capacity);
        localXYCorrelation.reserve(new_capacity);
        globalX.reserve(new_capacity);
        globalY.reserve(new_capacity);
        globalZ.reserve(new_capacity);
        detectorElementID.reserve(new_capacity);
     }
     void push_back(float a_localX, float a_localY, float a_localXError, float a_localYError, float a_localXYCorrelation,
                    float a_globalX, float a_globalY, float a_globalZ, unsigned long a_detectorElementID) {
        localX.push_back(a_localX);
        localY.push_back(a_localY);
        localXError.push_back(a_localXError);
        localYError.push_back(a_localYError);
        localXYCorrelation.push_back(a_localXYCorrelation);
        globalX.push_back(a_globalX);
        globalY.push_back(a_globalY);
        globalZ.push_back(a_globalZ);
        detectorElementID.push_back(a_detectorElementID);
     }
     void clear() {
        localX.clear();
        localY.clear();
        localXError.clear();
        localYError.clear();
        localXYCorrelation.clear();
        globalX.clear();
        globalY.clear();
        globalZ.clear();
        detectorElementID.clear();
     }
  };
  struct SDOInfoContainer {
   std::vector<std::vector<int> >                 *sdo_words = nullptr;
   std::vector<std::vector<std::vector<int> > >   *sim_depositsBarcode = nullptr;  // @TODO should be uint64 (schema evolution? )
   std::vector<std::vector<std::vector<float> > > *sim_depositsEnergy = nullptr;
     void reserve(std::size_t new_capacity) {
        if (sdo_words) sdo_words->reserve(new_capacity);
        if (sim_depositsBarcode) sim_depositsBarcode->reserve(new_capacity);
        if (sim_depositsEnergy)  sim_depositsEnergy->reserve(new_capacity);
     }
     void push_back(std::vector< int> &&a_sdo_words,
                    std::vector<std::vector< int> > &&a_sim_depositsBarcode,
                    std::vector<std::vector<float> > &&a_sim_depositsEnergy ) {
        assert( a_sdo_words.size() == a_sim_depositsEnergy.size() && a_sim_depositsEnergy.size() == a_sim_depositsBarcode.size());
        assert(    std::accumulate(a_sim_depositsBarcode.begin(),a_sim_depositsBarcode.begin(),0u,[](auto &elm) { return elm.size(); })
                == std::accumulate(a_sim_depositsEnergy.begin(), a_sim_depositsEnergy.begin(), 0u,[](auto &elm) { return elm.size(); }) );
        sdo_words->push_back(std::move(a_sdo_words));
        sim_depositsBarcode->push_back(std::move(a_sim_depositsBarcode));
        sim_depositsEnergy->push_back(std::move(a_sim_depositsEnergy));
     }
     void clear() {
        if (sdo_words) sdo_words->clear();
        if (sim_depositsBarcode) sim_depositsBarcode->clear();
        if (sim_depositsEnergy)  sim_depositsEnergy->clear();
     }
  };
  struct PixelContainer {
   std::vector<std::vector<int> >   *loc0idx = nullptr;
   std::vector<std::vector<int> >   *loc1idx = nullptr;
   std::vector<std::vector<float> > *charge = nullptr;
   std::vector<int>                 *waferID = nullptr;
     void reserve(std::size_t new_capacity) {
        if (loc0idx) loc0idx->reserve(new_capacity);
        if (loc1idx) loc1idx->reserve(new_capacity);
        if (charge) charge->reserve(new_capacity);
        if (waferID) waferID->reserve(new_capacity);
     }
     void push_back(std::vector<int> &&a_loc0idx, std::vector<int> &&a_loc1idx, std::vector<float> &&a_charge, int a_waferID) {
        assert( a_loc1idx.size() == a_loc1idx.size() && a_loc0idx.size() == a_charge.size());
        loc0idx->push_back(std::move(a_loc0idx));
        loc1idx->push_back(std::move(a_loc1idx));
        charge->push_back(std::move(a_charge));
        waferID->push_back(a_waferID);
     }
     void clear() {
        if (loc0idx) loc0idx->clear();
        if (loc1idx) loc1idx->clear();
        if (charge) charge->clear();
        if (waferID) waferID->clear();
     }
  };

  struct StripContainer {
   std::vector<int>               *SiWidth = nullptr;
   std::vector<int>               *hitsInThirdTimeBin = nullptr;
   std::vector<int>               *side = nullptr;
   std::vector<std::vector<int> > *strip = nullptr;
   std::vector<std::vector<int> > *timebin = nullptr;
   std::vector<std::vector<int> > *groupsize = nullptr;
     void reserve(std::size_t new_capacity) {
        if (SiWidth) SiWidth->reserve(new_capacity);
        if (hitsInThirdTimeBin) hitsInThirdTimeBin->reserve(new_capacity);
        if (side) side->reserve(new_capacity);
        if (strip) strip->reserve(new_capacity);
        if (timebin) timebin->reserve(new_capacity);
        if (groupsize) groupsize->reserve(new_capacity);
     }
     void push_back(int a_SiWidth, int a_hitsInThirdTimeBin, int a_side, std::vector<int> &&a_strip, std::vector<int> &&a_timebin, std::vector<int> &&a_groupsize) {
        assert( a_strip.size() == a_timebin.size() && a_strip.size() == a_groupsize.size());
        SiWidth->push_back(a_SiWidth);
        hitsInThirdTimeBin->push_back(a_hitsInThirdTimeBin);
        side->push_back(a_side);
        strip->push_back(std::move(a_strip));
        timebin->push_back(std::move(a_timebin));
        groupsize->push_back(std::move(a_groupsize));
     }
     void clear() {
        if (SiWidth) SiWidth->clear();
        if (hitsInThirdTimeBin) hitsInThirdTimeBin->clear();
        if (side) side->clear();
        if (strip) strip->clear();
        if (timebin) timebin->clear();
        if (groupsize) groupsize->clear();
     }
  };

  EventInfo               m_eventInfo;
  ParticleContainer       m_truthParticles;
  ParticleVertexContainer m_truthVertices;
  enum EClusterTypes {
     kPixel,
     kStrip,
     kNClusterTypes
  };
  std::array< ClusterContainer,kNClusterTypes > m_clusters;
  std::array< SDOInfoContainer,kNClusterTypes >     m_sdoInfo;
  PixelContainer m_pixelRDOs;
  StripContainer m_stripRDOs;
  bool           m_write = false;
protected:
  TLeaf *getLeaf(std::string leaf_name ) const;
  template <class T>
  TBranch *setObjectBranchAddress(const char *branch_name, T *adr);
  template <class T>
  void connectBranch(std::string branch_name, T *adr);

  void connectEventInfo(const std::string &prefix, ActsExamples::RootDigiBase::EventInfo &event_info);
  void connectParticles(const std::string &prefix, ActsExamples::RootDigiBase::ParticleContainer &particles);
  void connectVertices(const std::string &prefix, ActsExamples::RootDigiBase::ParticleVertexContainer &vertices);
  void connectClusters(const std::string &prefix, ActsExamples::RootDigiBase::ClusterContainer &clusters);
  void connectPixelRDOs(const std::string &prefix, ActsExamples::RootDigiBase::PixelContainer &pixel_rdos);
  void connectStripRDOs(const std::string &prefix, ActsExamples::RootDigiBase::StripContainer &strip_rdos);
  void connectSDOs(const std::string &prefix, ActsExamples::RootDigiBase::SDOInfoContainer &sdos);

  std::vector<std::function<void() > > m_cleanup;
};

class RootDigiWriter : public RootDigiBase, public virtual IWriter {
public:

   struct Config : public ConfigBase {};
     RootDigiWriter();
     RootDigiWriter(Config cfg, Acts::Logging::Level lvl);
     std::string name() const override { return "RootDigiWriter"; }

     ProcessCode write(const AlgorithmContext& context)  override;

     ProcessCode endRun() override { return endRunBase();}
};

class RootDigiReader : public RootDigiBase, public virtual IReader {
public:
     struct Config : public ConfigBase {};
     RootDigiReader();
     RootDigiReader(Config cfg, Acts::Logging::Level lvl);
     ~RootDigiReader() { endRunBase(); }

     std::string name() const override { return "RootDigiReader"; };

     std::pair<size_t, size_t> availableEvents() const override;

     ProcessCode read(const AlgorithmContext& context)  override;
};

inline TLeaf *ActsExamples::RootDigiBase::getLeaf(std::string leaf_name ) const {
   TLeaf *the_leaf = nullptr;
   if (!m_write) {
      the_leaf = m_tree->GetLeaf(leaf_name.c_str());
      if (!the_leaf) {
         std::string::size_type pos = leaf_name.find("Aux.");
         if (pos!= std::string::npos) {
            std::string tmp( leaf_name.substr(0,pos+3)+"Dyn."+leaf_name.substr(pos+4,leaf_name.size()-pos-4));
            std::cout << "DEBUG try alternative leave name " << leaf_name << " -> " << tmp << std::endl;
            the_leaf = m_tree->GetLeaf(tmp.c_str());
         }
      }
   }
   return the_leaf;
}

template <class T>
inline TBranch *ActsExamples::RootDigiBase::setObjectBranchAddress(const char *branch_name, T *adr) {
   static_assert( std::is_pointer<T>() );
   if (!*adr) {
      auto tmp = new typename std::remove_pointer<T>::type;
      *adr=tmp;
      const T expected_value = *adr;
      std::cout << "DEBUG register cleanup for : " << typeid(*expected_value).name() << " adr = " << *adr << std::endl;
      m_cleanup.emplace_back( [adr, expected_value]() {
            if (*adr == expected_value) {
               delete *adr;
               *adr=nullptr;
            }
         });
   }
   TBranch *branch_addr=nullptr;
   m_tree->SetBranchAddress(branch_name, *adr, &branch_addr);
   return branch_addr;
}

template <class T>
inline void ActsExamples::RootDigiBase::connectBranch(std::string branch_name, T *adr) {
   TBranch *branch_addr=nullptr;
   if (m_write) {
      branch_addr = m_tree->Branch(branch_name.c_str(), adr);
   }
   else {
      TLeaf *the_leaf = getLeaf(branch_name);
      TBranch *the_branch = (the_leaf ? the_leaf->GetBranch() : nullptr);
      if (the_branch) {
         EDataType the_type;
         TClass *the_class;
         if constexpr(std::is_pointer<T>()) {
               the_type = TDataType::GetType(typeid(typename std::remove_pointer<T>::type ));
               the_class = TClass::GetClass<typename std::remove_pointer<T>::type >();
         }
         else {
               the_type = TDataType::GetType(typeid(T));
               the_class = TClass::GetClass<T>();
         }
         std::cout << "DEBUG connect branch " << branch_name
                   << " -> " << (the_branch ? the_branch->GetName() : "")
                   << " destination type: " << (the_class ? the_class->GetName() : "")
                   << " " << the_type
                   << " branch " << (the_branch ? the_branch->GetClassName() : "")
                   << " leaf " << (the_leaf ? the_leaf->GetTypeName() : "")
                   << std::endl;

         if (the_class) {
            TClass *leaf_class=TClass::GetClass(the_leaf->GetTypeName());
            std::cout << "DEBUG branch type " << the_branch->GetClassName() << " ["
                      << the_branch->GetNleaves()
                      << " ]"
                      << " leaf " << the_leaf->GetTypeName() << " " << (leaf_class ? leaf_class->GetCheckSum() : 0u)
                      << " desired " <<  the_class->GetName() << std::endl;
            if (the_class->GetName() != the_leaf->GetTypeName() ) {
               TVirtualStreamerInfo *converter = the_class->GetConversionStreamerInfo(the_leaf->GetTypeName(),
                                                                                      the_class->GetClassVersion());
               if (!converter) {
                  ACTS_ERROR("No converter for " << branch_name << " from "
                             << the_leaf->GetTypeName() << " to " << the_class->GetName());
                  return;
               }
            }

            m_tree->SetBranchStatus(the_branch->GetName(),1);
            // check whether a pointer to a pointer to the object or just a pointer to the object is expected
            // (note T is a pointer to the object)
            if constexpr(std::is_pointer<T>()) {
                  if (the_branch->GetClassName() == the_leaf->GetTypeName()) {
                     m_tree->SetBranchAddress(the_branch->GetName(), adr, &branch_addr);
                  }
                  else {
                     branch_addr = setObjectBranchAddress(the_branch->GetName(),adr);
                  }
            }
            else {
               m_tree->SetBranchAddress(the_branch->GetName(), adr, &branch_addr);
            }
            branch_addr=the_branch;
         }
         else {
            if (a_branch->GetNleaves()==1 && a_branch->IsA() == TBranchElement::Class()) {
               TBranchElement *branch_element = static_cast<TBranchElement*>(a_branch);
               if (branch_element->GetType() == the_type()) {
                  m_tree->SetBranchAddress(the_branch->GetName(), adr, &branch_addr);
               }
            }
         }
      }
   }

   if (branch_addr) {
      ACTS_INFO("Connected branch " << branch_addr->GetName() );
      m_activeBranches.push_back(branch_addr);
   }
   else {
      ACTS_ERROR("Failed to " << (m_write ? "create" : "get") << " branch " << branch_name );
   }
}

}  // namespace ActsExamples
