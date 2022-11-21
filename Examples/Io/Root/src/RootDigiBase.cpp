// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootDigiBase.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
namespace ActsExamples {
using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
}

#include <ios>
#include <stdexcept>
#include <signal.h>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

ActsExamples::RootDigiBase::RootDigiBase(
    ActsExamples::RootDigiBase::ConfigBase cfg,
    const std::string &name,
    bool write,
    Acts::Logging::Level lvl)
   : m_cfg(cfg),
     m_logger(Acts::getDefaultLogger(name, lvl)),
     m_write(write)
 {
  // inputParticles is already checked by base constructor
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file path");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // open root file and create the tree
  if (m_write) {
     static const std::array<const char *,4> write_modes{"NEW","RECREATE","CREATE","UPDATE"};
     bool supported=false;
     for (const char *a_mode : write_modes ) {
        if (m_cfg.fileMode == a_mode) supported=true;
     }
     if (!supported) throw std::runtime_error("unsupported file mode for writing");
  }
  else {
     static const std::array<const char *,2> read_modes{"","READ"};
     bool supported=false;
     for (const char *a_mode : read_modes ) {
        if (m_cfg.fileMode == a_mode) supported=true;
     }
     if (!supported) throw std::runtime_error("unsupported file mode for reading");
  }
  m_file = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_file == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  m_file->cd();
  m_tree = static_cast<TTree *>(m_file->Get(m_cfg.treeName.c_str()));
  if (m_write) {
     if (m_tree) std::runtime_error("Tree already exists.");
     m_tree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
     if (m_tree == nullptr) {
        throw std::bad_alloc();
     }
  }
  else {
     // initially disable all branches, then enable all used branches
     m_tree->SetBranchStatus("*",0);
     m_tree->SetMakeClass(1);
  }

  std::cout << "DEBUG " << getpid() << " wait for signal CONT " << std::endl;
  kill(getpid(), SIGSTOP);

  // setup the branches
  connectEventInfo(m_cfg.eventInfoPrefix, m_eventInfo);
  connectParticles(m_cfg.truthParticlePrefix, m_truthParticles);
  connectVertices(m_cfg.truthVertexPrefix, m_truthVertices);
  if (m_cfg.measurementPrefix.size() != m_cfg.measurementType.size()) {
     throw std::range_error("Number of measurement prefixes and measurement types do not match.");
  }
  assert (m_clusters.size() == m_sdoInfos.size());
  int mask=0;
  for (unsigned int measurement_i=0; measurement_i<m_cfg.measurementPrefix.size(); ++measurement_i) {
     unsigned int type_i = m_cfg.measurementType[measurement_i];
     if (type_i >= kNClusterTypes) {
        throw std::runtime_error("Invalid measurement type. Only 0 : pixel, 1 : strips are currently understood.");
     }
     if (mask & (1<<type_i)) {
        throw std::runtime_error("MeasurementType used more than once.");
     }
     mask |= (1<<type_i);
     
     connectClusters(m_cfg.measurementPrefix[measurement_i], m_clusters[measurement_i]);     
     connectSDOs(m_cfg.measurementPrefix[measurement_i], m_sdoInfo[measurement_i]);
     switch (type_i) {
     case kPixel:
        connectPixelRDOs(m_cfg.measurementPrefix[measurement_i], m_pixelRDOs);
        break;
     case kStrip:
        connectStripRDOs(m_cfg.measurementPrefix[measurement_i], m_stripRDOs);
        break;
     default:
        break;
     }
  }
}


void ActsExamples::RootDigiBase::connectEventInfo(const std::string &prefix, ActsExamples::RootDigiBase::EventInfo &event_info) {
   //   connectBranch((prefix+"runNumber").c_str(),         &m_eventInfo.runNumber);
   connectBranch((prefix+"eventNumber").c_str(),       &event_info.eventNumber);
   //   connectBranch((prefix+"lumiBlock").c_str(),         &m_eventInfo.lumiBlock);
   //   connectBranch((prefix+"timeStamp").c_str(),         &m_eventInfo.timeStamp);
   //   connectBranch((prefix+"timeStampNSOffset").c_str(), &m_eventInfo.timeStampNSOffset);
   //   connectBranch((prefix+"bcid").c_str(),              &m_eventInfo.bcid);
}

void ActsExamples::RootDigiBase::connectParticles(const std::string &prefix, ActsExamples::RootDigiBase::ParticleContainer &particles) {
   connectBranch((prefix+"pdgId").c_str(),                    &particles.pdgId);
   connectBranch((prefix+"barcode").c_str(),                  &particles.barcode);
   connectBranch((prefix+"status").c_str(),                   &particles.status);
   connectBranch((prefix+"prodVtxLink").c_str(),              &particles.n_prodVtxLink);
   //   connectBranch((prefix+"prodVtxLink.m_persKey").c_str(),     particles.prodVtxLink_m_persKey);
   connectBranch((prefix+"prodVtxLink.m_persIndex").c_str(),   particles.prodVtxLink_m_persIndex);
   connectBranch((prefix+"decayVtxLink").c_str(),             &particles.n_decayVtxLink);
   //   connectBranch((prefix+"decayVtxLink.m_persKey").c_str(),    particles.decayVtxLink_m_persKey);
   connectBranch((prefix+"decayVtxLink.m_persIndex").c_str(),  particles.decayVtxLink_m_persIndex);
   connectBranch((prefix+"px").c_str(),                       &particles.px);
   connectBranch((prefix+"py").c_str(),                       &particles.py);
   connectBranch((prefix+"pz").c_str(),                       &particles.pz);
   connectBranch((prefix+"e").c_str(),                        &particles.e);
   connectBranch((prefix+"m").c_str(),                        &particles.m);
}

void ActsExamples::RootDigiBase::connectVertices(const std::string &prefix, ActsExamples::RootDigiBase::ParticleVertexContainer &vertices) {
   //   connectBranch((prefix+"id").c_str(),                    &vertices.id);
   connectBranch((prefix+"barcode").c_str(),               &vertices.barcode);
   //   connectBranch((prefix+"incomingParticleLinks").c_str(), &vertices.incomingParticleLinks);
   //   connectBranch((prefix+"outgoingParticleLinks").c_str(), &vertices.outgoingParticleLinks);
   connectBranch((prefix+"x").c_str(),                     &vertices.x);
   connectBranch((prefix+"y").c_str(),                     &vertices.y);
   connectBranch((prefix+"z").c_str(),                     &vertices.z);
   connectBranch((prefix+"t").c_str(),                     &vertices.t);
}
void ActsExamples::RootDigiBase::connectClusters(const std::string &prefix, ActsExamples::RootDigiBase::ClusterContainer &clusters) {
   //   connectBranch((prefix+"identifier").c_str(),         &clusters.identifier);
   //   connectBranch((prefix+"rdoIdentifierList").c_str(),  &clusters.rdoIdentifierList);
   connectBranch((prefix+"localX").c_str(),             &clusters.localX);
   connectBranch((prefix+"localY").c_str(),             &clusters.localY);
   connectBranch((prefix+"localXError").c_str(),        &clusters.localXError);
   connectBranch((prefix+"localYError").c_str(),        &clusters.localYError);
   connectBranch((prefix+"localXYCorrelation").c_str(), &clusters.localXYCorrelation);
   connectBranch((prefix+"globalX").c_str(),            &clusters.globalX);
   connectBranch((prefix+"globalY").c_str(),            &clusters.globalY);
   connectBranch((prefix+"globalZ").c_str(),            &clusters.globalZ);
   connectBranch((prefix+"detectorElementID").c_str(),  &clusters.detectorElementID);
}
void ActsExamples::RootDigiBase::connectPixelRDOs(const std::string &prefix, ActsExamples::RootDigiBase::PixelContainer &pixel_rdos) {
   connectBranch((prefix+"rdo_phi_pixel_index").c_str(), &pixel_rdos.loc0idx);
   connectBranch((prefix+"rdo_eta_pixel_index").c_str(), &pixel_rdos.loc0idx);
   connectBranch((prefix+"rdo_charge").c_str(),          &pixel_rdos.charge);
   connectBranch((prefix+"waferID").c_str(),             &pixel_rdos.waferID);
}
void ActsExamples::RootDigiBase::connectStripRDOs(const std::string &prefix, ActsExamples::RootDigiBase::StripContainer &strip_rdos) {
   connectBranch((prefix+"SiWidth").c_str(),            &strip_rdos.SiWidth);
   connectBranch((prefix+"hitsInThirdTimeBin").c_str(), &strip_rdos.hitsInThirdTimeBin);
   connectBranch((prefix+"side").c_str(),               &strip_rdos.side);
   connectBranch((prefix+"rdo_strip").c_str(),          &strip_rdos.strip);
   connectBranch((prefix+"rdo_timebin").c_str(),        &strip_rdos.timebin);
   connectBranch((prefix+"rdo_groupsize").c_str(),      &strip_rdos.groupsize);
}

void ActsExamples::RootDigiBase::connectSDOs(const std::string &prefix, ActsExamples::RootDigiBase::SDOInfoContainer &sdos) {
   connectBranch((prefix+"sdo_words").c_str(),           &sdos.sdo_words);
   connectBranch((prefix+"sdo_depositsBarcode").c_str(), &sdos.sim_depositsBarcode);
   connectBranch((prefix+"sdo_depositsEnergy").c_str(),  &sdos.sim_depositsEnergy);
}

ActsExamples::RootDigiBase::~RootDigiBase() {
  if (m_file != nullptr) {
    m_file->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootDigiBase::endRunBase() {
  if (m_write && m_tree && m_file) {
    m_file->cd();
    m_file->Write();
    ACTS_INFO("Wrote tree '" << m_tree->GetName() << "' in '"
              << m_file->GetName()  << "'");
  }

  if (m_file != nullptr) {
    m_file->Close();
  }

  return ProcessCode::SUCCESS;
}

ActsExamples::RootDigiWriter::RootDigiWriter() {
   m_cfg.fileMode="RECREATE";
   m_write=true;
}
ActsExamples::RootDigiWriter::RootDigiWriter(ActsExamples::RootDigiWriter::Config cfg, Acts::Logging::Level lvl)
   :ActsExamples::RootDigiBase::RootDigiBase(cfg, "RootDigiWriter", true, lvl)
{
}
  
ActsExamples::ProcessCode ActsExamples::RootDigiWriter::write(const AlgorithmContext& ctx) {
  if (m_file == nullptr) {
    ACTS_ERROR("Missing output file");
    return ProcessCode::ABORT;
  }
  const ActsExamples::SimParticleContainer &particles = ctx.eventStore.get<ActsExamples::SimParticleContainer>(m_cfg.particles);
   
  // ensure exclusive access to tree/file while writing
  std::lock_guard<std::mutex> lock(m_mutex);

  m_eventInfo.eventNumber = ctx.eventNumber;
  m_truthParticles.clear();
  std::vector< std::array<float,4 > > vertex_position;
  
  m_truthVertices.clear();
  m_truthParticles.clear();
  m_truthParticles.reserve( particles.size());
  for (const auto& particle : particles) {
     //     void push_back( int a_pdg_id, int a_barcode, int a_status, int a_prod_vtx_idx, int a_decay_vtx_idx,
     //              float a_px, float a_py, float a_pz, float a_m) {
     const auto p = particle.absoluteMomentum() / Acts::UnitConstants::MeV;
     const auto particle_direction = particle.unitDirection();
     m_truthParticles.push_back(particle.pdg(),particle.particleId().particle(),2,-particle.particleId().vertexPrimary(),0,
                                p * particle_direction.x(),
                                p * particle_direction.y(),
                                p * particle_direction.z(),
                                particle.mass() / Acts::UnitConstants::MeV);
     if (particle.particleId().vertexPrimary()>= vertex_position.size()) {
        vertex_position.resize(particle.particleId().vertexPrimary()+1);
     }
     vertex_position[particle.particleId().vertexPrimary()] = std::array<float,4> {
        static_cast<float>(particle.fourPosition().x() / Acts::UnitConstants::mm),
        static_cast<float>(particle.fourPosition().y() / Acts::UnitConstants::mm),
        static_cast<float>(particle.fourPosition().z() / Acts::UnitConstants::mm),
        static_cast<float>(particle.fourPosition().w() / Acts::UnitConstants::ns)
     }; 
  }
  int idx=0;
  for (const std::array<float, 4> &vertex : vertex_position) {
     m_truthVertices.push_back(-idx, vertex[0],vertex[1],vertex[2],vertex[3]);
     ++idx;
  }

  const auto& measurements =
     ctx.eventStore.get<ActsExamples::MeasurementContainer>(m_cfg.measurements);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.measurementParticlesMap);

  //m_cfg.trackingGeometry), ctx.geoContext
  for (unsigned int meas_i=0; meas_i<measurements.size(); ++meas_i) {

     const auto &meas=measurements[meas_i];

     const auto& slink =
        std::visit([](const auto& x) { return &x.sourceLink(); }, meas);

     const auto geoId = slink->geometryId();

     const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
     auto [localPos, localCov, meas_dim] = std::visit(
                                                      [](const auto& measurement) {
                                                         auto expander = measurement.expander();
                                                         Acts::BoundVector par = expander * measurement.parameters();
                                                         Acts::BoundSymMatrix cov =
                                                         expander * measurement.covariance() * expander.transpose();
                                                         // extract local position
                                                         Acts::Vector2 lpar(par[Acts::eBoundLoc0], par[Acts::eBoundLoc1]);
                                                         // extract local position covariance.
                                                         Acts::SymMatrix2 lcov = cov.block<2, 2>(Acts::eBoundLoc0, Acts::eBoundLoc0);
                                                         return std::make_tuple(lpar, lcov, measurement.size());
                                                      },
                                                      meas);
     Acts::Vector3 globalPos = surface->localToGlobal(ctx.geoContext, localPos, Acts::Vector3());

     EClusterTypes cluster_type = meas_dim == 2 ? kPixel : kStrip; // @TODO should use geoId

     //     void push_back(float a_localX, float a_localY, float a_localXError, float a_localYError, float a_localXYCorrelation,
     //               float a_globalX, float a_globalY, float a_globalZ, unsigned long a_detectorElementID) {
     m_clusters[ cluster_type ].push_back( localPos[0],localPos[1],localCov(0,0),localCov(0,0),localCov(0,1),
                                                             globalPos[0],globalPos[1],globalPos[2], geoId.value());

     std::vector<int> sdo_word;
     std::vector< std::vector<int> > sdo_barcode;
     std::vector< std::vector<float> > sdo_depositsEnergy;
     auto matchingParticles = makeRange(hitParticlesMap.equal_range(meas_i));
     sdo_word.reserve(sdo_word.size());
     sdo_barcode.reserve(sdo_barcode.size());
     sdo_depositsEnergy.reserve(sdo_depositsEnergy.size());

     for (auto hitParticle : matchingParticles) {
        sdo_word.push_back(0);
        sdo_barcode.emplace_back();
        sdo_barcode.back().push_back(hitParticle.second.particle());
        sdo_depositsEnergy.emplace_back();
        sdo_depositsEnergy.back().push_back(std::numeric_limits<float>::epsilon());
     }
     m_sdoInfo[cluster_type].push_back(std::move(sdo_word), std::move(sdo_barcode), std::move(sdo_depositsEnergy));
  }

  m_tree->Fill();


  return ProcessCode::SUCCESS;
}

ActsExamples::RootDigiReader::RootDigiReader() {
   m_cfg.fileMode="READ";
   m_write=false;
}

ActsExamples::RootDigiReader::RootDigiReader(ActsExamples::RootDigiReader::Config cfg, Acts::Logging::Level lvl)
   : ActsExamples::RootDigiBase::RootDigiBase(cfg, "RootDigiWriter", false, lvl)
{
}

std::pair<size_t, size_t> ActsExamples::RootDigiReader::availableEvents() const {
   return std::make_pair(0u, m_tree->GetEntries() );
}

ActsExamples::ProcessCode ActsExamples::RootDigiReader::read(const AlgorithmContext& ctx) {
    auto entry = ctx.eventNumber;
    // if (not m_cfg.orderedEvents and entry < m_entryNumbers.size()) {
    //   entry = m_entryNumbers[entry];
    // }
    // @TODO in case the input is a chain, detect tree changes and re-collect the active branches.
    //    m_tree->GetEntry(entry);
    Long64_t current_tree_entry = m_tree->LoadTree(entry);
    if (current_tree_entry < 0) return ProcessCode::ABORT;

    for (TBranch *branch : m_activeBranches) {
       Long_t read_bytes = branch->GetEntry(current_tree_entry);
       if (read_bytes <=0 ) {
          ACTS_ERROR("Failed to read branch " << branch->GetName() << " from file " << (m_tree->GetDirectory() ? m_tree->GetDirectory()->GetFile()->GetName() : "") );
       }
       else {
          TClass *the_class;
          EDataType the_type;
          branch->GetExpectedType(the_class,the_type);
          std::cout << "DEBUG " << name() << " read " << read_bytes << " for " << branch->GetName()
                    << " adr= " << static_cast<const void *>(branch->GetAddress())
                    << " type: " << (the_class && the_class->GetTypeInfo() ? the_class->GetTypeInfo()->name() : "")
                    << " " << the_type
                    << std::endl;
       }
       
    }

    std::cout << "DEBUG " << name() << " " << m_eventInfo.eventNumber
              << " barcode=" << m_truthParticles.barcode.size() << " (" << static_cast<const void *>(&m_truthParticles.barcode) << ") "
              << " globalX=" << m_clusters[kPixel].globalX.size()  << " (" << static_cast<const void *>(&m_clusters[kPixel].globalX) << ") "
              << std::endl;
    unsigned int counter=0;
    for ( const auto &elm : m_truthParticles.barcode) {
       if (counter %20==0) {
          if (counter>0) std::cout << std::endl;
          std::cout << "DEBUG " << name();
       }
       std::cout << " " << elm << std::endl;
       ++counter;
    }
    if (counter>0) std::cout << std::endl;
    counter=0;
    if (m_clusters[kPixel].globalX) {
    for ( const auto &elm : *(m_clusters[kPixel].globalX)) {
       if (counter %20==0) {
          if (counter>0) std::cout << std::endl;
          std::cout << "DEBUG " << name();
       }
       std::cout << " " << elm << std::endl;
       ++counter;
    }
    }
    if (counter>0) std::cout << std::endl;

    // std::cout <<"DEBUG " << name() << " wait for CONT " << getpid() << std::endl;
    // kill (getpid(), SIGSTOP);
  return ProcessCode::SUCCESS;
}
