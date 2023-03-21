// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"

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

#include "Acts/EventData/SourceLink.hpp"


#include <ios>
#include <stdexcept>
#include <list>
#include <signal.h>
#include <unistd.h>
#include <regex>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"

//#define DEBUG_READ 1


namespace Test {
   class IdentifiableDetectorElement : public Acts::DetectorElementBase
   {
   public:
      IdentifiableDetectorElement(const Acts::Surface *surface, uint64_t detector_element_id, float thickness)
        : m_surface(surface),
          m_transform( surface->transform(Acts::GeometryContext()) ),
          m_id(detector_element_id),
          m_thickness(thickness)
      {}

      ~IdentifiableDetectorElement() {
         //         std::cout << "IdentifiableDetectorElement::dtor " << m_id << std::endl;
      }

      /// @param gctx The current geometry context object, e.g. alignment
      const Acts::Transform3& transform([[maybe_unused]] const Acts::GeometryContext& gctx) const override
        { return m_transform; }
      /// Return surface representation
      const Acts::Surface& surface() const override
        { return *m_surface; }

      double thickness() const override { return m_thickness; }

      uint64_t detectorId() const { return m_id; }
   private:
      const Acts::Surface *m_surface;
      Acts::Transform3 m_transform;
      uint64_t m_id;
      float m_thickness;
   };
}

namespace {
   inline double sqr(double a) { return a*a; }
   unsigned int gCounter=0;
}

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
  setupTree();
  setupValidationHists();
  readSurfaceCenters();
 }
void ActsExamples::RootDigiBase::setupTree() {
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

void ActsExamples::RootDigiBase::readSurfaceCenters() {
   if (!m_cfg.surfaceCenterFile.empty()) {
      std::ifstream in(m_cfg.surfaceCenterFile);
      while (in) {
         uint64_t id;
         std::array<double,3> center;
         in >> id >> center[0] >> center[1];
         if (in) {
            in >> center[2];
            // std::cout << id << " " << center[0] << " " << center[1] << " " << center[2] << std::endl;
            m_surfaceCenters.insert( std::make_pair(id, center) );
         }
      }
      std::cout << "DEBUG read " << m_surfaceCenters.size() << " surface centers." << std::endl;
   }
}

void ActsExamples::RootDigiBase::setupValidationHists() {
   if (!m_cfg.outputFilePath.empty()) {
      static const std::array<std::string,kNHistCategories> cat_names{ std::string("_geoMatched"),std::string("_noGeoMatch")};
      static const std::array<std::string,kNHistTypes> type_names{ std::string("MeasPos_RZ"), std::string("MeasPos_RPhi"),std::string("MeasPos_ZPhi")};

      //      const std::array<float,16> z_slices{-2800,-2500,-2200,-1954,-1900,-1700,-1600,-1490,
      //                                   1490, 1600, 1700, 1900, 1954, 2200, 2500, 2800 };
      //      const float r=350;
      std::array<const char *,19> zslice_names{"Pixel",
                                               "StripOuter",
                                               "StripEM8","StripEM7","StripEM6","StripEM5","StripEM4","StripEM3","StripEM2", "StripEM1",
                                               "StripB",
                                               "StripEP1","StripEP2","StripEP3","StripEP4","StripEP5","StripEP6", "StripEP7","StripEP8" };

      for (unsigned int cat_i=0; cat_i < kNHistCategories; ++cat_i ) {
         //enum EHistType {kPosREta, kPosRPhi, kPosEtaPhi, kNHistTypes};
         m_validationeHist.at(kPosRZ).at(cat_i)=new TH2F((type_names[kPosRZ]+cat_names[cat_i]).c_str(),
                                                           (type_names[kPosRZ]+cat_names[cat_i]).c_str(),
                                                           1000,350,1200,2000,14000,3000);
         m_validationeHist.at(kPosRPhi).at(cat_i)=new TH2F((type_names[kPosRPhi]+cat_names[cat_i]).c_str(),
                                                           (type_names[kPosRPhi]+cat_names[cat_i]).c_str(),
                                                           1000,350,1200,500,0, M_PI);
         m_validationeHist.at(kPosZPhi).at(cat_i)=new TH2F((type_names[kPosZPhi]+cat_names[cat_i]).c_str(),
                                                           (type_names[kPosZPhi]+cat_names[cat_i]).c_str(),
                                                           2000,1400,3000,500,0, M_PI);
         unsigned int z_i=0;
         m_validationeHist.at(kPosPixel+z_i).at(cat_i)=new TH2F((std::string(zslice_names.at(z_i))+cat_names[cat_i]).c_str(),
                                                                (std::string(zslice_names.at(z_i))+cat_names[cat_i]).c_str(),
                                                                350,0,350,100,-M_PI, M_PI);
         for (; ++z_i<kNHistTypes-kPosPixel;) {
            m_validationeHist.at(kPosPixel+z_i).at(cat_i)=new TH2F((std::string(zslice_names.at(z_i))+cat_names[cat_i]).c_str(),
                                                                    (std::string(zslice_names.at(z_i))+cat_names[cat_i]).c_str(),
                                                                    500,350,1050,100,-M_PI, M_PI);
         }

      }
   }
}

void ActsExamples::RootDigiBase::fillValidationHists(const Acts::Vector3 &global_pos, EHistCategory category) {
   if (!m_cfg.outputFilePath.empty()) {
      double phi = Acts::VectorHelpers::phi(global_pos);
      //      double eta = Acts::VectorHelpers::eta(global_pos);
      double z = global_pos[2];
      double r = Acts::VectorHelpers::perp(global_pos);

      static const std::array<float,18> z_slices{-5000,-2800,-2500,-2200,-1954,-1900,-1700,-1600,-1490,
                                                  1490, 1600, 1700, 1900, 1954, 2200, 2500, 2800,5000 };
      static const float r_pixel=350;

      if (r>350 && z>1400 && phi>=0) {
         m_validationeHist.at(kPosRZ).at(category)->Fill(r, z);
         m_validationeHist.at(kPosRPhi).at(category)->Fill(r,phi);
         m_validationeHist.at(kPosZPhi).at(category)->Fill(z,phi);
         m_validationeHist.at(kPosZPhi).at(category)->Fill(z,phi);
      }
      if (r<r_pixel) {
         m_validationeHist.at(kPosPixel).at(category)->Fill(r,phi);
      }
      else {
         unsigned int z_i=0;
         if (std::abs(z)<3000) {
            for (; ++z_i<z_slices.size();) {
               if (z<z_slices.at(z_i)) break;
            }
         }
         m_validationeHist.at(kPosStripOuter+z_i).at(category)->Fill(r,phi);
      }
   }
}

namespace {
   class GlobalRootRestore {
   public:
      GlobalRootRestore() : m_lastFile(gFile), m_lastDir(gDirectory) {}
      ~GlobalRootRestore() { gFile= m_lastFile; gDirectory = m_lastDir; if (gDirectory) {gDirectory->cd(); } }
   private:
      TFile      *m_lastFile;
      TDirectory *m_lastDir;
   };
}
void ActsExamples::RootDigiBase::writeValidationHists() const {
  if (!m_cfg.outputFilePath.empty()) {
     GlobalRootRestore restore;

     std::unique_ptr<TFile> output( TFile::Open(m_cfg.outputFilePath.c_str(),"RECREATE"));
     if (output && output->IsWritable()) {
        output->cd();
        for (unsigned int type_i=0; type_i < kNHistTypes; ++type_i) {
           for (unsigned int cat_i=0; cat_i < kNHistCategories; ++cat_i ) {
              m_validationeHist.at(type_i).at(cat_i)->Write();
           }
        }
     }
     else {
        ACTS_ERROR("Failed to create output file " << m_cfg.outputFilePath);
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
   connectBranch((prefix+"prodVtxLink_").c_str(),             &particles.n_prodVtxLink);
   //   connectBranch((prefix+"prodVtxLink.m_persKey").c_str(),     particles.prodVtxLink_m_persKey);
   connectBranch((prefix+"prodVtxLink.m_persIndex").c_str(),   particles.prodVtxLink_m_persIndex);
   connectBranch((prefix+"decayVtxLink_").c_str(),             &particles.n_decayVtxLink);
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
   connectBranch((prefix+"rdo_eta_pixel_index").c_str(), &pixel_rdos.loc1idx);
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
   if (m_tree) {
      m_events=m_tree->GetEntries();
      if (m_cfg.singleParticle) {
         for (TBranch *branch : m_activeBranches) {
            std::string branch_name (branch->GetName());
            if (branch_name.find(cfg.truthParticlePrefix) != std::string::npos &&  branch_name.find("status")!= std::string::npos) {

               uint64_t  entries = m_tree->GetEntries();
               m_events=0;
               for (uint64_t entry_i=0; entry_i< entries; ++ entry_i) {
                  Long64_t current_tree_entry = m_tree->LoadTree(entry_i);
                  if (current_tree_entry < 0) {
                     break;
                  }
                  branch->GetEntry(current_tree_entry);
                  std::size_t end_particle =m_truthParticles.status.size();
                  for (unsigned int particle_i = 0; particle_i < end_particle; ++particle_i) {
                     if (m_truthParticles.status.at(particle_i)==1) {
                        ++m_events;
                     }
                  }
               }
               std::cout << "DEBUG RootDigiReader particles " << m_events << std::endl;
               break;
            }
         }
          //          if (m_truthParticles.status.at(particle_i)==1) {

      }
   }
}

std::pair<size_t, size_t> ActsExamples::RootDigiReader::availableEvents() const {
   // if (m_cfg.singleParticle) {
   //    // count particles
   //    return std::make_pair(0u, 1000 /*m_tree->GetEntries()*/ );
   // }
   // else {
   //    return std::make_pair(0u, m_tree->GetEntries() );
   // }
   return std::make_pair(0u, m_events );
}


namespace {
   bool altBoundsCheck(uint64_t detElementId,
                       const Acts::Vector2 &locpos,
                       const std::vector<Acts::Vector2> &corners, bool verbose=false) {
      const double p0=locpos[0];
      const double p1=locpos[1];
      const double x0=corners[0][0];
      const double y0=corners[0][1];
      const double x1=corners[1][0];
      const double y1=corners[1][1];
      const double x2=corners[2][0];
      const double y2=corners[2][1];
      const double x3=corners[3][0];
      const double y3=corners[3][1];

      // express the coordiantes of the measurement in some wierd wierd coordinate system
      // for which s=0,t=0 is the center of the polygon and the 4 corners
      // are reached for (s=-1,+1 , t=-1,+1)
      // s goes in direction : from  edge corner 0, corner 1 to edge corner 3, corner 2
      // t goes in direction : from  edge corner 0, corner 3 to edge corner 1, corner 2
      // locX = -(s*t*x3-t*x3+s*x3-x3-s*t*x2-t*x2-s*x2-x2+s*t*x1+t*x1-s*x1-x1-s*t*x0+t*x0+s*x0-x0)/4
      // locY = -(s*t*y3-t*y3+s*y3-y3-s*t*y2-t*y2-s*y2-y2+s*t*y1+t*y1-s*y1-y1-s*t*y0+t*y0+s*y0-y0)/4
      //
      // solution for s,t:
      // s=  ((t+1)*x3+(t+1)*x2+(1-t)*x1+(1-t)*x0-4*p0)
      //    /((t+1)*x3+(-t-1)*x2+(t-1)*x1+(1-t)*x0);
      //
      // t*(2*x2*y3-2*p0*y3-2*x3*y2+2*p0*y2+2*x0*y1-2*p0*y1-2*x1*y0+2*p0*y0+2*p1*x3-2*p1*x2+2*p1*x1-2*p1*x0)
      // +t^2*(x2*y3-2*x1*y3+x0*y3-x3*y2+x1*y2+2*x3*y1-x2*y1-x0*y1-x3*y0+x1*y0)
      // +x2*y3+2*x1*y3-x0*y3-2*p0*y3-x3*y2-x1*y2+2*p0*y2-2*x3*y1+x2*y1-x0*y1+2*p0*y1+x3*y0+x1*y0-2*p0*y0+2*p1*x3-2*p1*x2-2*p1*x1+2*p1*x0
      // =0

      // a*t^2 + b*t + c = 0
      // (%o206) b = 2*(x2*y3+p0*(-y3+y2-y1+y0)-x3*y2+x0*y1-x1*y0+p1*(x3-x2+x1-x0))
      // (%o207) a = x2*y3-x1*y3-x3*y2+x0*y2+x3*y1-x0*y1-x2*y0+x1*y0
      // (%o208) c = x2*y3+x1*y3+p0*(-2*y3+2*y2+2*y1-2*y0)-x3*y2-x0*y2-x3*y1-x0*y1+x2*y0+x1*y0+p1*(2*x3-2*x2-2*x1+2*x0)
      const double b = 2*(x2*y3-x3*y2+x0*y1-x1*y0 + p0 * (-y3+y2-y1+y0) + p1 * (x3-x2+x1-x0));
      const double a = x2*y3-x1*y3-x3*y2+x0*y2+x3*y1-x0*y1-x2*y0+x1*y0;
      const double c = x2*y3+x1*y3-x3*y2-x0*y2-x3*y1-x0*y1+x2*y0+x1*y0 + p0 * (-2*y3+2*y2+2*y1-2*y0) + p1 * (2*x3-2*x2-2*x1+2*x0);
                   
      // => t =  (-b +- sqrt ( b^2 - 4*a*c))/(2*a)
      const  double q = b*b - 4 * a * c;
      if (q<0) {
         std::cout << "DEBUG failed to compute skewed coordinates for detector element " << detElementId
                   << " localPos " << p0 << ", " << p1
                   << std::endl;
         return false;
      }
      else {
         double t1 = (-b - sqrt(q))/(2*a);
         double t2 = (-b + sqrt(q))/(2*a);
         auto solution_s = [&](double t) {
            return ((t+1)*x3+(t+1)*x2+(1-t)*x1+(1-t)*x0-4*p0)
                  /((t+1)*x3+(-t-1)*x2+(t-1)*x1+(1-t)*x0);
         };
         auto the_test = [&](double s, double t) {
            return std::make_pair(-(s*t*x3-t*x3+s*x3-x3-s*t*x2-t*x2-s*x2-x2+s*t*x1+t*x1-s*x1-x1-s*t*x0+t*x0+s*x0-x0)/4,
                                  -(s*t*y3-t*y3+s*y3-y3-s*t*y2-t*y2-s*y2-y2+s*t*y1+t*y1-s*y1-y1-s*t*y0+t*y0+s*y0-y0)/4);
         };
         bool is_inside = ( (std::abs(solution_s(t1))<1 && std::abs(t1)<1) || (std::abs(solution_s(t2))<1 && std::abs(t2)<1));
         constexpr auto sqr = [](double arg) {return arg*arg;};
         if (verbose && !is_inside) {
            std::cout << "DEBUG skewed coordinates for " << detElementId << " " << p0 << ", " << p1 << " -> "
                      << solution_s(t1) << " , " << t1
                      << ( (std::abs(solution_s(t1))<1 && std::abs(t1)<1) ? " * " : "")
                      << " or " << solution_s(t2) << " , " << t2
                      << ( (std::abs(solution_s(t2))<1 && std::abs(t2)<1) ? " * " : "")
                      << ( (std::abs(solution_s(t1))<1 && std::abs(t1)<1) || (std::abs(solution_s(t2))<1 && std::abs(t2)<1) ? " inside" : " OUTSIDE")
                      << std::endl;
         }
         auto test1 = the_test ( solution_s(t1),t1);
         auto test2 = the_test ( solution_s(t2),t2);
         double chi2_1= sqrt( sqr(test1.first-p0)/1e-5 + sqr(test1.second-p1)/1e-5 );
         double chi2_2= sqrt( sqr(test2.first-p0)/1e-5 + sqr(test2.second-p1)/1e-5 );
         if (chi2_1>1 || chi2_2>1) {
            std::cout << "DEBUG skewed coordinates test " << detElementId << " " << p0 << ", " << p1 << " " 
                      << test1.first << " , "  << test1.second
                      <<  " (" << chi2_1 << ") "
                      << " or "
                      << test2.first << " , "  << test2.second
                      <<  " (" << chi2_2 << ") "
                      << std::endl;
         }
         return is_inside;
      }
   }
}

ActsExamples::ProcessCode ActsExamples::RootDigiReader::read(const AlgorithmContext& ctx) {
    auto entry = ctx.eventNumber;

    // if (not m_cfg.orderedEvents and entry < m_entryNumbers.size()) {
    //   entry = m_entryNumbers[entry];
    // }
    // @TODO in case the input is a chain, detect tree changes and re-collect the active branches.
    //    m_tree->GetEntry(entry);
    std::lock_guard<std::mutex> lock(m_mutex);

    
    if (m_cfg.stop) {
       std::cout << "DEBUG " << getpid() << " wait for signal CONT " << std::endl;
       kill(getpid(), SIGSTOP);
    }

  // std::regex ignore(".*HGTD_NOTHING.*");

  // m_cfg.trackingGeometry->visitSurfaces([&ctx, &ignore](const Acts::Surface* surface) {
  //       std::cout  << "DEBUG visit surfaces id " <<  ( surface->associatedLayer() ? surface->associatedLayer()->geometryId() : 0u)
  //                  << " "
  //                  <<  ( surface->associatedLayer()  && surface->associatedLayer()->trackingVolume()
  //                        ? surface->associatedLayer()->trackingVolume()->volumeName()
  //                        : "(unknown)" )
  //                  << std::endl;
  //       if (surface->associatedLayer()
  //           && surface->associatedLayer()->trackingVolume()
  //           && std::regex_match(surface->associatedLayer()->trackingVolume()->volumeName(), ignore)) return;
  //       const auto &trans = surface->transform(ctx.geoContext);
  //       for (unsigned int i=0; i< trans.rows(); ++i) {
  //          std::cout << i << ":";
  //          for (unsigned int j=0; j< trans.cols(); ++j) {
  //             std::cout << " " << std::setw(14) << trans(i,j);
  //          }
  //          std::cout << std::endl;
  //       }

  //    });

    // if (!m_file) {
    //    setupTree();
    // }
    bool read_next_entry=true;
    if (m_cfg.singleParticle) {
       if (m_nextEntry > 0 ) {
          m_currentParticle = entry-m_entryOffset;
          if (m_currentParticle<m_stableParticles.size()) {
             std::cout << "DEBUG \"read\" particle " << m_currentParticle << std::endl;
             read_next_entry=false;
          }
          else {
             m_entryOffset+=m_currentParticle;
             m_currentParticle=0;
          }
       }
       entry = m_nextEntry;
    }

    if (read_next_entry) {
       std::cout << "DEBUG read entry " << entry << std::endl;
    Long64_t current_tree_entry = m_tree->LoadTree(entry);
    if (current_tree_entry < 0) return ProcessCode::ABORT;
    m_nextEntry=entry+1;

    for (TBranch *branch : m_activeBranches) {
       if (!branch) continue;
       Long_t read_bytes = branch->GetEntry(current_tree_entry);
       if (read_bytes <=0 ) {
          ACTS_ERROR("Failed to read branch " << branch->GetName() << " from file " << (m_tree->GetDirectory() ? m_tree->GetDirectory()->GetFile()->GetName() : "") );
       }
       else {
          std::cout << "DEBUG " << name() << " read " << read_bytes << " for " << branch->GetName()
                    << " adr= " << static_cast<const void *>(branch->GetAddress())
                    << std::endl;
          TClass *the_class=nullptr;
          EDataType the_type;
          branch->GetExpectedType(the_class,the_type);
          std::cout << "DEBUG " << name() << " read " << read_bytes << " for " << branch->GetName()
                    << " adr= " << static_cast<const void *>(branch->GetAddress())
                    << " type: " << (the_class && the_class->GetTypeInfo() ? the_class->GetTypeInfo()->name() : "")
                    << " " << the_type
                    << std::endl;
       }

    }
    if (m_cfg.singleParticle) {

       std::set<int> particles_with_measurements;
       for (const SDOInfoContainer &sdo_info : m_sdoInfo ) {
          for (unsigned int meas_i=0; meas_i < sdo_info.sim_depositsBarcode->size(); ++meas_i) {
             for ( const std::vector<int>  &sim_hit_list : sdo_info.sim_depositsBarcode->at(meas_i) ) {
                for ( const int  &sim_barcode : sim_hit_list ) {
                   particles_with_measurements.insert( sim_barcode);
                }
             }
          }
       }
       m_stableParticles.clear();
       unsigned int end_particle=m_truthParticles.barcode.size();
       m_stableParticles.reserve( m_truthParticles.barcode.size() );
       for (unsigned int particle_i = 0; particle_i < end_particle; ++particle_i) {
          if (m_truthParticles.status.at(particle_i)==1) {
             //if (   particles_with_measurements.find( m_truthParticles.barcode.at(particle_i) )
             // != particles_with_measurements.end()) {
             float pt = sqrt( sqr( m_truthParticles.px.at(particle_i))
                              + sqr(  m_truthParticles.py.at(particle_i) ));

             m_stableParticles.push_back( std::make_pair( pt, particle_i ) );
          }
       }
       std::sort( m_stableParticles.begin(), m_stableParticles.end(), [](std::pair<float, unsigned int> &a,
                                                                         std::pair<float, unsigned int> &b) {
                     return a.first > b.first;
                  });
       unsigned int part_i=0;
       for (const std::pair<float, unsigned int> &particle_idx : m_stableParticles) {
          std::cout << "DEBUG stable particle "
                    << part_i
                    << " idx " << particle_idx.second
                    << " barcode " << m_truthParticles.barcode.at(particle_idx.second)
                    << " pdgId "  << m_truthParticles.pdgId.at(particle_idx.second)
                    << " pt "  << particle_idx.first
                    << std::endl;
          ++part_i;
       }
    }

    if ( m_stat.empty() ) {
       m_stat.emplace_back( m_cfg.eventInfoPrefix,std::vector<std::pair<std::string,Stat> >() );
       m_stat.emplace_back( m_cfg.truthParticlePrefix, std::vector<std::pair<std::string,Stat> >() );
       m_stat.emplace_back( m_cfg.truthVertexPrefix, std::vector<std::pair<std::string,Stat> >() );
       for (unsigned int measurement_i=0; measurement_i<m_cfg.measurementPrefix.size(); ++measurement_i) {
          m_stat.emplace_back( m_cfg.measurementPrefix[measurement_i], std::vector<std::pair<std::string,Stat> >() );
          m_stat.emplace_back( m_cfg.measurementPrefix[measurement_i]+" SDO", std::vector<std::pair<std::string,Stat> >() );
          switch (m_cfg.measurementType[measurement_i] ) {
          case kPixel:
             m_stat.emplace_back( m_cfg.measurementPrefix[measurement_i]+" pixelRDOs", std::vector<std::pair<std::string,Stat> >() );
             break;
          case kStrip:
             m_stat.emplace_back( m_cfg.measurementPrefix[measurement_i]+" stripRDOs", std::vector<std::pair<std::string,Stat> >() );
             break;
          }
       }
    }


    unsigned int stat_i=0;
    m_eventInfo.fillStat(m_stat.at(stat_i++).second);
    m_truthParticles.fillStat(m_stat.at(stat_i++).second);
    m_truthVertices.fillStat(m_stat.at(stat_i++).second);
    for (unsigned int measurement_i=0; measurement_i<m_cfg.measurementPrefix.size(); ++measurement_i) {
       m_clusters[measurement_i].fillStat(m_stat.at(stat_i++).second);
       m_sdoInfo[measurement_i].fillStat(m_stat.at(stat_i++).second);
       switch (m_cfg.measurementType[measurement_i] ) {
       case kPixel:
          m_pixelRDOs.fillStat(m_stat.at(stat_i++).second);
          break;
       case kStrip:
          m_stripRDOs.fillStat(m_stat.at(stat_i++).second);
          break;
       }
    }

# ifdef DEBUG_READ
    std::cout << "DEBUG " << name() << " " << m_eventInfo.eventNumber
              << " barcode=" << m_truthParticles.barcode.size() << " (" << static_cast<const void *>(&m_truthParticles.barcode) << ") "
              << " globalX=" << m_clusters[kPixel].globalX->size()  << " (" << static_cast<const void *>(m_clusters[kPixel].globalX) << ") "
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
# endif
    }

  std::map<int, ActsFatras::Barcode>  barcode_map = convertParticles(ctx);
  convertMeasurements(ctx, barcode_map);
  return ProcessCode::SUCCESS;
}

std::map<int, ActsFatras::Barcode> ActsExamples::RootDigiReader::convertParticles(const AlgorithmContext& ctx) {
  SimParticleContainer particles;
  m_truthParticles.checkDimensions();
  m_truthVertices.checkDimensions();
  std::map<int, ActsFatras::Barcode> barcode_map;
  unsigned int start_particle=(m_cfg.singleParticle ? m_stableParticles.at(m_currentParticle).second : 0 );
  unsigned int end_particle=(m_cfg.singleParticle ? start_particle+1 : m_truthParticles.pdgId.size() );
  std::cout << "DEBUG convert particles  " << start_particle << " .. " << end_particle;
  if (m_cfg.singleParticle) {
     std::cout << " idx " << m_stableParticles.at(m_currentParticle).second
               << " barcode " << m_truthParticles.barcode.at(start_particle)
               << " pdgId "  << m_truthParticles.pdgId.at(start_particle)
               << " pt "  << m_stableParticles.at(m_currentParticle).first;
  }
  std::cout << std::endl;

  for (unsigned int particle_i = start_particle; particle_i < end_particle; ++particle_i) {
     // reject all non stable truth particles
     if (!m_cfg.singleParticle) {
     if (m_truthParticles.pdgId.at(particle_i) ==21 || m_truthParticles.pdgId.at(particle_i)==22)  {
        continue;
     }
     if (m_truthParticles.status.at(particle_i) != 1) {
        std::cout << "DEBUG reject non stable truth particle " << m_truthParticles.barcode.at(particle_i)
                  << " "  << m_truthParticles.pdgId.at(particle_i)
                  <<  " " << m_truthParticles.status.at(particle_i)
                  << std::endl;
        continue;
     }
     double pt = sqrt( sqr(m_truthParticles.px.at(particle_i)+m_truthParticles.py.at(particle_i)) );
     if (pt<100) continue;
     }

     // somehow map the input barcode and vertex index to the ActsFatras::Barcode
     unsigned int vertex_i = (m_truthParticles.n_prodVtxLink>0 && particle_i < static_cast<unsigned int>(m_truthParticles.n_prodVtxLink)
                              ? m_truthParticles.prodVtxLink_m_persIndex[ particle_i ]
                              : std::numeric_limits<unsigned int>::max() );
     if  (vertex_i & 0xff000000) {
        if (m_truthParticles.n_prodVtxLink<=0 || particle_i >= static_cast<unsigned int>(m_truthParticles.n_prodVtxLink))  {
           ACTS_ERROR("No vertex for  " << particle_i << " !< " << m_truthParticles.n_prodVtxLink << " -> " << (vertex_i & 0xffffff) );
        }
        else {
           ACTS_ERROR("Precision loss for vertex index of particle " << particle_i << " : " << vertex_i
                      << " -> " << (vertex_i & 0xffffff) );
        }
     }

     const auto pid = ActsFatras::Barcode(0u)
        .setParticle(m_truthParticles.barcode.at(particle_i) & 0xffff)
        .setSubParticle((m_truthParticles.barcode.at(particle_i)>>16) & 0xffff)
        .setVertexPrimary( (vertex_i & 0xfff) )
        .setVertexSecondary( ((vertex_i>>12) & 0xfff) );

     std::pair<std::map<int, ActsFatras::Barcode>::const_iterator, bool>
        insert_result =  barcode_map.insert( std::make_pair(m_truthParticles.barcode.at(particle_i), pid) );
     if (!insert_result.second) {
        ACTS_ERROR("Barcode " << m_truthParticles.barcode.at(particle_i) << " not unique. Old dest : " << insert_result.first->second
                   << " New dest " << pid);
     }
     std::cout << "DEBUG " << insert_result.first->first  << " -> " << insert_result.first->second
               << " ( " << m_truthParticles.barcode.at(particle_i) << " -> " << pid << ")" << std::endl;

     //  set vertix part
     ActsFatras::Particle particle(pid,
                                   Acts::PdgParticle(m_truthParticles.pdgId.at(particle_i)),
                                   m_truthParticles.pdgId.at(particle_i) >= 0 ? 1. : -1, // charge
                                   m_truthParticles.m.at(particle_i));
     Acts::Vector3 dir(m_truthParticles.px.at(particle_i),
                       m_truthParticles.py.at(particle_i),
                       m_truthParticles.pz.at(particle_i));
     double abs_p = dir.norm();
     particle.setAbsoluteMomentum(abs_p);
     particle.setDirection(dir * (1/abs_p));
     if (vertex_i< m_truthVertices.x.size()) {
        particle.setPosition4(m_truthVertices.x.at(vertex_i),
                              m_truthVertices.y.at(vertex_i),
                              m_truthVertices.z.at(vertex_i),
                              m_truthVertices.t.at(vertex_i));
     }



     particles.insert( std::move(particle) );
  }
  ctx.eventStore.add(m_cfg.particles, std::move(particles));
  return barcode_map;
}
   std::ostream &operator<<(std::ostream &out, const Acts::Transform3 &trans) {
      for (unsigned int i=0; i< trans.rows(); ++i) {
         out << i << ":";
         for (unsigned int j=0; j< trans.cols(); ++j) {
            out << " " << std::setw(14) << trans(i,j);
         }
         out << std::endl;
      }
      return out;
   }

namespace {

   struct Counter {
      Counter() {
         for (unsigned int i=0; i<3; ++i) {
            m_box[i][0]= std::numeric_limits<float>::max();
            m_box[i][1]=-std::numeric_limits<float>::max();
         }
      }
      //      Counter &operator++() { ++m_counter; return *this; }
      void updateBox(const Acts::Vector3 &global_pos) {
         for (unsigned int i=0; i<3; ++i) {
            m_box[i][0]=std::min(m_box[i][0],static_cast<float>(global_pos[i]));
            m_box[i][1]=std::max(m_box[i][1],static_cast<float>(global_pos[i]));
         }
         ++m_counter;
      }
      unsigned int m_counter =0;
      std::array<std::array<float,2> , 3> m_box;
      // {
      //    {std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()},
      //    {std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()},
      //    {std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()}
      // } ;
   };
   struct DistanceRange {
      void update(float distance) {
         m_min=std::min(m_min,distance);
         m_max=std::max(m_max,distance);
         ++m_n;
       }
      float m_min=std::numeric_limits<float>::max();
      float m_max=-std::numeric_limits<float>::max();
      unsigned int m_n=0;
   };
   ActsExamples::Stat log_distance( 10, -6,4);
   std::unordered_map<  Acts::GeometryIdentifier, Counter > m_destinationMultiplicity;
   std::unordered_map< unsigned long, DistanceRange > match_distance_range;
   void checkDestinationMultiplicity() {
      for (const std::pair<  const Acts::GeometryIdentifier, Counter > & elm : m_destinationMultiplicity) {
         if (elm.second.m_counter != 1) {
            std::cout << "DEBUG multiplicity  " << elm.first << " : " << std::setw(3) << elm.second.m_counter << " |";
            for (unsigned int i=0; i<3; ++i) {
               std::cout << " " << std::setw(14) << elm.second.m_box[i][0] << "..." << elm.second.m_box[i][1];
                  //<< std::endl;
            }
            std::cout << std::endl;
         }
      }
      DistanceRange total;
      ActsExamples::Stat total_stat_count;
      ActsExamples::Stat total_stat_min;
      ActsExamples::Stat total_stat_max;
      unsigned long max_distance=0;
      for (const std::pair< const unsigned long, DistanceRange > & distance_range : match_distance_range ) {
         if (distance_range.second.m_max > total.m_max) { 
            max_distance = distance_range.first;
         }
         total_stat_max.add(distance_range.second.m_max);
         total_stat_min.add(distance_range.second.m_min);
         total_stat_count.add(distance_range.second.m_n);
         total.update(distance_range.second.m_min);
         total.update(distance_range.second.m_max);
      }
      std::cout << "DEBUG distance range " << total.m_min << " ... " << total.m_max
                << " worst : " << max_distance << std::endl;
      std::cout << "DEBUG min distance: " << total_stat_min << std::endl;
      std::cout << "DEBUG max distance: " << total_stat_max << std::endl;
      std::cout << "DEBUG count: "        << total_stat_count << std::endl;
      std::cout << "DEBUG log distance: " << log_distance << std::endl;
   }

   static void checkSurface(unsigned int layer_i,
                            unsigned int surface_i,
                            const Acts::Surface *a_surface,
                            const Acts::GeometryContext &geo_ctx,
                            const Acts::Vector3 &global_pos,
                            Acts::BoundaryCheck bcheck,
                            const double max_error) {
      if (a_surface) {
         Acts::Vector3 dummy_momentum{1.,1.,1.};
         Acts::Result<Acts::Vector2> result = a_surface->globalToLocal(geo_ctx,
                                                                       global_pos,
                                                                       dummy_momentum,
                                                                       max_error);
         if (result .ok()) {
            Acts::Vector3 on_surface = a_surface->localToGlobal(geo_ctx,
                                                                result.value(),
                                                                dummy_momentum);
            double distance2 = (global_pos - on_surface).squaredNorm();

            const auto &trans = a_surface->transform(geo_ctx);
            std::array<Acts::BoundaryCheck,3> tests{ Acts::BoundaryCheck(true,true, 1e-1,1e-1),
                  Acts::BoundaryCheck(true,true, 1e-2,1e-2),
                  Acts::BoundaryCheck(true,true, 1e-3,1e-3) };

            unsigned int tol_i;
            for (tol_i=tests.size(); tol_i-->0; ) {
               if (a_surface->bounds().inside(result.value(), tests.at(tol_i))) break;
            }
            ++tol_i;

            bool on_layer = a_surface->associatedLayer()
               ? a_surface->associatedLayer()->isOnLayer(geo_ctx,global_pos, bcheck)
               : false;
            Acts::BoundaryCheck bcheck2(true,true, 1,1);
            bool on_layer2 = (on_layer || a_surface->associatedLayer()->isOnLayer(geo_ctx,global_pos, bcheck2));

            Acts::Vector3 center_pos = a_surface->center(geo_ctx);
            std::cout << "DEBUG " <<gCounter << " " << layer_i << " " << surface_i
                      << " : " << sqrt(sqr(center_pos[0])+sqr(center_pos[1])) << ",  " << center_pos[2]
                      << " -> " << sqrt(sqr(on_surface[0])+sqr(on_surface[1])) << ",  " << on_surface[2]
                      << " | " << sqrt(distance2) << " [ " << result.value()[0]  << ", "  << result.value()[1] << " ] "
                      << (tol_i>0 ? " (in-side bounds) " : " ") << (tol_i) << " "
                      << (on_layer ? " (on layer)" : (a_surface->associatedLayer()
                                                      ? ( on_layer2 ? "(~on layer)" : "")
                                                      : "(NO LAYER)")) << std::endl
               //                      << "  " << a_surface->bounds()
                      << std::endl;
            if (!on_layer && tol_i>0 && a_surface->associatedLayer()) {
               const Acts::Layer *a_layer = a_surface->associatedLayer();
               Acts::Result<Acts::Vector2> result2 = a_layer->surfaceRepresentation().globalToLocal(geo_ctx,
                                                                              global_pos,
                                                                              dummy_momentum,
                                                                              3);
               double distance3=std::numeric_limits<double>::max();
               if (result2.ok()) {
                  Acts::Vector3 on_surface3 = a_layer->surfaceRepresentation().localToGlobal(geo_ctx,
                                                                                             result2.value(),
                                                                                             dummy_momentum);
                  distance3 = (global_pos - on_surface3).squaredNorm();
               }
               else {
                  distance3=-1;
               }

               bool on_layer3 = a_surface->associatedLayer()->isOnLayer(geo_ctx,global_pos, bcheck);


               std::cout << global_pos << " "
                         << a_layer->geometryId()
                         << (a_layer->representingVolume() ? " (vol) " : "(surf)")
                         << distance3
                         << " thickness : " << a_layer->thickness() << " : "
                         << std::tuple<const Acts::Surface&,
                                       const Acts::GeometryContext&>(a_layer->surfaceRepresentation(),
                                                                     geo_ctx)
                         << std::endl;
               if (a_layer->representingVolume()) {
                  std::cout << "   " << typeid(a_layer->representingVolume()).name() << std::endl;
               }
            }
         }
      }
   }

   static void dumpVolumeInfo(const Acts::TrackingVolume *final_volume,
                              const Acts::GeometryContext &geo_ctx,
                              const Acts::Vector3 &global_pos,
                              double min_distance2) {
      if (final_volume && std::abs(min_distance2)>1) {
         if (gCounter<100) {
            std::cerr << "ERROR no match  for " << sqrt(sqr(global_pos[0])+sqr(global_pos[1]))
                      << ", " <<  global_pos[2] << " -> " << final_volume->volumeName()
                      << " min surf. dist. " << sqrt(min_distance2)
                      << std::endl;
            for (const auto &boundary_surface :  final_volume->boundarySurfaces()) {
                  std::cout << "DEBUG " << gCounter << " " << final_volume->volumeName()
                            << std::tuple<const Acts::Surface&,
                                          const Acts::GeometryContext&>(boundary_surface->surfaceRepresentation(),
                                                                        geo_ctx)
                            << std::endl;

            }
            bool on_layer=false;
            Acts::BoundaryCheck bcheck(true,true, 1e-1,1e-1);
            const double max_error = 1e-1;
            {
               unsigned layer_i=0;
               Acts::Vector3 dummy_momentum{1.,1.,1.};
               for (const auto &a_layer : final_volume->confinedLayers()->arrayObjects() ) {
                  if (a_layer->isOnLayer(geo_ctx, global_pos, bcheck)) {
                     std::cout << "DEBUG " << gCounter << " " << layer_i << " " << final_volume->volumeName()
                               << " thickness : " << a_layer->thickness() << " "
                               << std::tuple<const Acts::Surface&,
                                             const Acts::GeometryContext&>(a_layer->surfaceRepresentation(),
                                                                           geo_ctx)
                               << std::endl;
                     unsigned int surface_i=0;
                     for (const Acts::Surface *a_surface: a_layer->surfaceArray()->surfaces()) {
                        checkSurface(layer_i, surface_i, a_surface, geo_ctx,global_pos, bcheck, max_error);
                        ++surface_i;
                     }
                     on_layer=true;
                  }
                  ++layer_i;
               }
            }

            if (!on_layer) {
               unsigned layer_i=0;
               for (const auto &a_layer : final_volume->confinedLayers()->arrayObjects() ) {
                  unsigned int n_surfaces=0;
                  if (a_layer->surfaceArray()) {
                     n_surfaces += a_layer->surfaceArray()->surfaces().size();
                  }
                  std::cout << "DEBUG " << gCounter << " " << layer_i << " " << final_volume->volumeName()
                            << " thickness : " << a_layer->thickness() << " surfaces : " << n_surfaces << " "
                            << std::tuple<const Acts::Surface&,
                                          const Acts::GeometryContext&>(a_layer->surfaceRepresentation(),
                                                                        geo_ctx)
                            << std::endl;
                  ++layer_i;
               }

               unsigned int idx=0;
               final_volume->visitSurfaces([&geo_ctx, &idx, &global_pos, &bcheck, &max_error](const Acts::Surface* surface) {
                     checkSurface(99, idx, surface, geo_ctx,global_pos, bcheck, max_error);

                     // Acts::Vector3 center_pos = surface->center(geo_ctx);
                     // std::cout << "DEBUG " << gCounter << " " << idx << " surface of " << final_volume->volumeName()
                     //           << " : " << sqrt(sqr(center_pos[0])+sqr(center_pos[1])) << ",  " << center_pos[2]
                     //           << std::endl;
                     ++idx;
                  });
            }
            ++gCounter;
         }
      }
   }

   static std::string cornersToStr(const Acts::Surface *a_surface,
                                   const Acts::GeometryContext &geo_ctx) {
      std::stringstream out;
      std::vector<Acts::Vector2> corners;
      if (a_surface) {
         const Acts::AnnulusBounds *annulus_bounds = dynamic_cast<const Acts::AnnulusBounds *>( &a_surface->bounds() );
         if (annulus_bounds) {
            corners = annulus_bounds->corners();
         }
         else {
            const Acts::PlanarBounds *bounds = dynamic_cast<const Acts::PlanarBounds *>( &a_surface->bounds() );
            if (bounds) {
               corners = bounds->vertices(1);
            }
         }
      }

      // const Acts::AnnulusBounds * annulus_bounds = dynamic_cast<const Acts::AnnulusBounds *>( &a_surface->bounds() );
      // if (annulus_bounds) {
      Acts::Vector3 dummy_momentum{1.,1.,1.};


      //    std::vector<Acts::Vector2> corners = annulus_bounds->corners();
      unsigned int corner_i=0;
      for (const Acts::Vector2 &a_corner  : corners) {
         Acts::Vector3 glob = a_surface->localToGlobal(geo_ctx,
                                                       a_corner,
                                                       dummy_momentum);
         out << " DEBUG surfaceBound corners " << a_surface->geometryId().value() << " "
             << corner_i << " : " << glob[0] << " , " << glob[1] << " , " << glob[2] << std::endl;
         ++corner_i;
      }
      return out.str();
   }

   template <typename visitor_t>
   void visitLayerSurfaces(const Acts::TrackingVolume *volume,
                           const Acts::GeometryContext&geo_ctx,
                           const Acts::Vector3 &global_pos,
                           const Acts::BoundaryCheck &bcheck,
                           visitor_t&& visitor) {
      if (volume) {
         for (const auto &a_layer : volume->confinedLayers()->arrayObjects() ) {
            if (a_layer->isOnLayer(geo_ctx, global_pos, bcheck)) {
               for (const Acts::Surface *a_surface: a_layer->surfaceArray()->surfaces()) {
                  visitor(a_surface);
               }
               return;
            }
         }
         std::cout << "DEBUG no matching layer in " << volume->volumeName() << " for "
                   << sqrt(sqr(global_pos[0])+sqr(global_pos[1])) << " , "
                   << global_pos[2]
                   << " search all surfaces."
                   << std::endl;
         dumpVolumeInfo(volume, geo_ctx, global_pos, std::numeric_limits<double>::max());
         volume->visitSurfaces(visitor);
      }
   }

   static void testSurface(const ActsExamples::AlgorithmContext& ctx,
                           const Acts::Surface *surface, const Acts::Vector3 &global_pos, double max_error) {
      Acts::BoundaryCheck bcheck(true,true, 1e-2,1e-2);
      Acts::Vector3 dummy_momentum{};
      Acts::Result<Acts::Vector2> result = surface->globalToLocal(ctx.geoContext,
                                                                  global_pos,
                                                                  dummy_momentum,
                                                                  max_error);
      if (result .ok()) {
         //    //                            << " trans " <<surface->transform(ctx.geoContext)
         //           << std::endl;
         const auto &trans = surface->transform(ctx.geoContext);
                  //                  Acts::BoundaryCheck bcheck(true,true, 1e-1,1e-1);
         bool inside_bounds = surface->bounds().inside(result.value(), bcheck);

         // if (!inside_bounds ) {
         //    const Acts::AnnulusBounds * annulus_bounds = dynamic_cast<const Acts::AnnulusBounds *>( &surface->bounds() );
         //    if (annulus_bounds) {
         //       std::vector<Acts::Vector2> corners = annulus_bounds->corners();
         //       bool alt_inside_bounds = altBoundsCheck(detector_element_id, result.value(),corners);
         //       if (alt_inside_bounds != inside_bounds) {
         //          std::cout << "DEBUG bounds check disagree " << inside_bounds << " != " << alt_inside_bounds << std::endl;
         //       }
         //    }
         // }
         double distance2=-1.;
         if (inside_bounds) {
            Acts::Vector3 on_surface = surface->localToGlobal(ctx.geoContext,
                                                              result.value(),
                                                              dummy_momentum);
            double distance2 = (global_pos - on_surface).squaredNorm();
            double distance=sqrt(distance2);
            // std::cout << "DEBUG " << global_pos << " -> " << result.value() 
            //           << " id " <<  ( surface->associatedLayer() ? surface->associatedLayer()->geometryId() : 0u)
            //           << " distance = " << sqrt(distance2)
            //           << std::endl;
            // if (distance2<min_distance2) {
            //    min_distance2 = distance2;
            //    min_surface = surface;
            // }
         }
      }
   }

   std::pair<Acts::GeometryIdentifier,bool> findGeoId(const ActsExamples::AlgorithmContext& ctx,
                                                      const Acts::TrackingGeometry &trackingGeometry,
                                                      std::unordered_map<unsigned long, Acts::GeometryIdentifier> &geo_map,
                                                      unsigned long detector_element_id,
                                                      const Acts::Vector3 &global_pos,
                                                      const double &max_error,
                                                      bool dump_surface=false,
                                                      bool debug=false) {

      std::pair<Acts::GeometryIdentifier,bool> ret{};
      //std::regex ignore(".*HGTD.*");
      std::unordered_map<unsigned long, Acts::GeometryIdentifier>::const_iterator
         geomap_iter = geo_map.find(detector_element_id);
      if (debug || geomap_iter == geo_map.end()) {
         //{

         const Acts::TrackingVolume* volume = trackingGeometry.lowestTrackingVolume(ctx.geoContext,
                                                                              global_pos);
         const Acts::TrackingVolume* final_volume = volume;
         // while (volume  != final_volume) {
         //    final_volume = volume->lowestTrackingVolume(ctx.geoContext,
         //                                               global_pos,
         //                                               max_error);
         // }
         double min_distance2=std::numeric_limits<double>::max();
         if (final_volume) {
            const Acts::Surface* min_surface=nullptr;

            //            final_volume->visitSurfaces()
            Acts::BoundaryCheck bcheck(true,true, 1e-2,1e-2);
            visitLayerSurfaces(final_volume,
                               ctx.geoContext,
                               global_pos,
                               bcheck,
                               [&ctx,
                                &bcheck,
                                &max_error,
                                &global_pos,
                                &min_distance2,
                                &min_surface,
                                detector_element_id](const Acts::Surface* surface) {
                  // if (surface->associatedLayer()
                  //     && surface->associatedLayer()->trackingVolume()
                  //     && std::regex_match(surface->associatedLayer()->trackingVolume()->volumeName(), ignore)) return;
                                  uint64_t sur_geo_id = surface->geometryId().value();
                                  bool dump = (detector_element_id == 576815344803381248); //(sur_geo_id == 1441152705392670976 || sur_geo_id == 1441152705392672512);
               Acts::Vector3 dummy_momentum{};
               Acts::Result<Acts::Vector2> result = surface->globalToLocal(ctx.geoContext,
                                                                           global_pos,
                                                                           dummy_momentum,
                                                                           max_error);
               if (result .ok()) {
                  //    //                            << " trans " <<surface->transform(ctx.geoContext)
                  //           << std::endl;
                  const auto &trans = surface->transform(ctx.geoContext);
                  //                  Acts::BoundaryCheck bcheck(true,true, 1e-1,1e-1);
                  bool inside_bounds = surface->bounds().inside(result.value(), bcheck);


                  // if (!inside_bounds ) {
                  //    const Acts::AnnulusBounds * annulus_bounds = dynamic_cast<const Acts::AnnulusBounds *>( &surface->bounds() );
                  //    if (annulus_bounds) {
                  //       std::vector<Acts::Vector2> corners = annulus_bounds->corners();
                  //       bool alt_inside_bounds = altBoundsCheck(detector_element_id, result.value(),corners);
                  //       if (alt_inside_bounds != inside_bounds) {
                  //          std::cout << "DEBUG bounds check disagree " << inside_bounds << " != " << alt_inside_bounds << std::endl;
                  //       }
                  //    }
                  // }
                  double distance2=-1.;
                  if (inside_bounds) {
                     Acts::Vector3 on_surface = surface->localToGlobal(ctx.geoContext,
                                                                       result.value(),
                                                                       dummy_momentum);
                     distance2 = (global_pos - on_surface).squaredNorm();
                     // std::cout << "DEBUG " << global_pos << " -> " << result.value() 
                     //           << " id " <<  ( surface->associatedLayer() ? surface->associatedLayer()->geometryId() : 0u)
                     //           << " distance = " << sqrt(distance2)
                     //           << std::endl;
                     if (distance2<min_distance2) {
                        min_distance2 = distance2;
                        min_surface = surface;
                     }

                  }

                  if (dump) {
                     std::cout << "DEBUG " << global_pos[0] << " , " << global_pos[1] << ", " << global_pos[2]
                               << " -> " << result.value()[0]  << ", " << result.value()[1]
                               << " det " << detector_element_id
                               << " geoId " <<  sur_geo_id 
                               << " " << distance2
                               << (inside_bounds ? " inside " : "")
                               << " " << (min_surface ? min_surface->geometryId().value() : 0u)
                               << std::endl;
                     // std::cout << "DEBUG " << global_pos << " -> " << result.value() 
                     //           << " id " <<  sur_geo_id
                     //           << " det " << static_cast<const void *>(surface->associatedDetectorElement())
                  }

               }
            });
            if (min_surface) {
               const Acts::Layer* layer = min_surface->associatedLayer();
               if (dump_surface) {
                     std::cout << "DEBUG surface for detector element " << detector_element_id << " mapped to "
                               << min_surface->geometryId();
                     min_surface->toStream(ctx.geoContext, std::cout);
                     std::cout << std::endl;
                     const Acts::AnnulusBounds * annulus_bounds = dynamic_cast<const Acts::AnnulusBounds *>( &min_surface->bounds() );
                     if (annulus_bounds) {
                        Acts::Vector3 dummy_momentum{1.,1.,1.};
                        std::vector<Acts::Vector2> corners = annulus_bounds->corners();
                        unsigned int corner_i=0;
                        for (const Acts::Vector2 &a_corner  : corners) {
                           Acts::Vector3 glob = min_surface->localToGlobal(ctx.geoContext,
                                                                           a_corner,
                                                                           dummy_momentum);
                           std::cout << " DEBUG surfaceBound corners " << detector_element_id << " " 
                                     << corner_i << " : " << glob[0] << " , " << glob[1] << " , " << glob[2] << std::endl;
                           ++corner_i;
                        }
                     }

               }
               if (layer) {
                  if (min_distance2>1e-4*1e-4) {
                     Acts::Vector3 dummy_momentum{1.,1.,1.};
                     Acts::Vector3 on_surface{0.,0.,0.};
                     Acts::Result<Acts::Vector2> result = min_surface->globalToLocal(ctx.geoContext,
                                                                                 global_pos,
                                                                                 dummy_momentum,
                                                                                 max_error);
                     if (result .ok()) {
                        on_surface = min_surface->localToGlobal(ctx.geoContext,
                                                                result.value(),
                                                                dummy_momentum);
                     }

                     std::cout << "WARNING large matching distance for detector element " << detector_element_id << " mapped to "
                               << min_surface->geometryId()
                               << " pos: "
                               << global_pos[0] << ", " << global_pos[1] << ", " << global_pos[2]
                               << " -> "
                               << on_surface[0] << ", " << on_surface[1] << ", " << on_surface[2]
                               << " distance=" << sqrt(min_distance2) << std::endl;
                     const Acts::AnnulusBounds * annulus_bounds = dynamic_cast<const Acts::AnnulusBounds *>( &min_surface->bounds() );
                     if (annulus_bounds) {
                        std::vector<Acts::Vector2> corners = annulus_bounds->corners();
                        unsigned int corner_i=0;
                        for (const Acts::Vector2 &a_corner  : corners) {
                           Acts::Vector3 glob = min_surface->localToGlobal(ctx.geoContext,
                                                                           a_corner,
                                                                           dummy_momentum);
                           std::cout << " DEBUG surfaceBound corners " << detector_element_id << " "
                                     << corner_i << " : " << glob[0] << " , " << glob[1] << " , " << glob[2] << std::endl;
                           ++corner_i;
                        }
                     }
                     dumpVolumeInfo(final_volume, ctx.geoContext, global_pos, min_distance2);


                  }

                  match_distance_range[detector_element_id].update( sqrt(min_distance2) );
                  log_distance.add( min_distance2>0 ? log( sqrt(min_distance2) ) / log(10.) : -308 );
                  if (geomap_iter != geo_map.end()) {

                     uint64_t sur_geo_id = min_surface->geometryId().value();
                     bool dump = detector_element_id == 576815344803381248; // (sur_geo_id == 1441152705392670976 || sur_geo_id == 1441152705392672512);
                     if (dump) {
                        std::cout << "DEBUG detector element " << detector_element_id << " mapped to "
                                  << geomap_iter->second.value() << " " << geomap_iter->second
                                  << " distance=" << sqrt(min_distance2)
                                  << " pos " << global_pos[0] << ", " << global_pos[1] << ", " << global_pos[2]
                                  << std::endl;
                     }

                     if (min_surface->geometryId() == geomap_iter->second) return std::make_pair(geomap_iter->second,
                                                                                           sqrt(min_distance2) < 1e-3);
                     std::cout << "ERROR Geometry element for " << detector_element_id  <<  " and pos: "
                               << global_pos[0] << ", " << global_pos[1] << ", " << global_pos[2] << " ."
                               << min_surface->geometryId()
                               << " Previously detector element " << detector_element_id
                               << " was matched to " << geomap_iter->second
                               << std::endl;

                     auto surface = trackingGeometry.findSurface(geomap_iter->second);
                     if (surface) {
                        std::cout << "DEBUG old surface " << surface->geometryId().value() << " "
                                  << surface->geometryId() << " surface " << static_cast<const void *>(surface)
                                  << " " << std::tuple<const Acts::Surface&,
                                                       const Acts::GeometryContext&>(*surface,ctx.geoContext)
                                  << cornersToStr(surface, ctx.geoContext)
                                  << std::endl;
                     }
                     std::cout << "DEBUG new surface " << min_surface->geometryId().value() << " "
                               << min_surface->geometryId() << " surface " <<  static_cast<const void *>(min_surface)
                               << " " << std::tuple<const Acts::Surface&,
                                                    const Acts::GeometryContext&>(*min_surface,ctx.geoContext)
                               << cornersToStr(min_surface, ctx.geoContext)
                               << std::endl;
                     if (dump) {
                           dumpVolumeInfo(final_volume, ctx.geoContext, global_pos, std::numeric_limits<double>::max());
                     }

                  }
                  else {

                  std::pair<std::unordered_map<unsigned long, Acts::GeometryIdentifier>::const_iterator, bool>
                     insert_result = geo_map.insert( std::make_pair(detector_element_id, min_surface->geometryId()));
                  if (insert_result.second) {
                     m_destinationMultiplicity[insert_result.first->second].updateBox(global_pos);
                     if (debug) {
                        std::cout << "DEBUG detector element " << detector_element_id << " mapped to "
                                  << insert_result.first->second.value() << " " << insert_result.first->second
                                  << " distance=" << sqrt(min_distance2)
                                  << " pos " << global_pos[0] << ", " << global_pos[1] << ", " << global_pos[2]
                                  << std::endl;
                        if (min_surface->geometryId().value() == 576460889742376960) {
                           dumpVolumeInfo(final_volume, ctx.geoContext, global_pos, std::numeric_limits<double>::max());
                        }

                     }
                     return std::make_pair(insert_result.first->second,
                                           sqrt(min_distance2) < 1e-3);
                  }
                  else {
                     std::cout << "WARNING detector element " << detector_element_id << " is already mapped to "
                               << insert_result.first->second << " new: " << min_surface->geometryId() << std::endl;
                     return std::make_pair(min_surface->geometryId(),
                                           sqrt(min_distance2) < 1e-3);
                  }
                  }
               }
               else {
                  if (geomap_iter != geo_map.end()) {
                     std::cout << "ERROR No associated layer for matching geometry element input:" << detector_element_id  <<  " pos: "
                               << global_pos[0] << ", " << global_pos[1] << ", " << global_pos[2] << " ."
                               << " Although detector element " << detector_element_id
                               << " was matched to " << geomap_iter->second
                               << std::endl;
                     dumpVolumeInfo(final_volume, ctx.geoContext, global_pos, std::numeric_limits<double>::max());
                  }
               }

            }
         }

         std::stringstream msg;
         msg << "No matching geometry element found for " << detector_element_id  <<  " pos: "
             << global_pos[0] << ", " << global_pos[1] << ", " << global_pos[2] << " .";
         //         throw std::runtime_error( msg.str() );
         std::cerr << msg.str() << std::endl;
         dumpVolumeInfo(final_volume, ctx.geoContext, global_pos, min_distance2);

         return std::make_pair(0u,false);

      }
      return std::make_pair(geomap_iter->second,false);
   }
   ActsFatras::Barcode findParticleBarCode(const std::map<int, ActsFatras::Barcode> &barcode_map, int source_barcode) {
      std::map<int, ActsFatras::Barcode>::const_iterator  barcode_iter = barcode_map.find(source_barcode);
      if (barcode_iter == barcode_map.end()) {
         std::stringstream msg;
         msg << "Unknown particle barcode " << source_barcode <<  " .";
         throw std::runtime_error( msg.str() );
      }
      return barcode_iter->second;
   }
}

void ActsExamples::RootDigiReader::dumpParticleHits(const Acts::GeometryContext& gctx,
                                                    const Acts::TrackingGeometry &trackingGeometry,
                                                    std::unordered_map<unsigned long, Acts::GeometryIdentifier> &geo_map) const {
   std::map<int, unsigned int> barcodeToParticle;
   std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > particleToHit;
   std::vector< std::pair<float, unsigned int> > particle_order;

  unsigned int start_particle=(m_cfg.singleParticle ? m_stableParticles.at(m_currentParticle).second : 0 );
  unsigned int end_particle=(m_cfg.singleParticle ? start_particle+1 : m_truthParticles.pdgId.size() );

   for (unsigned int particle_i = start_particle; particle_i < end_particle; ++particle_i) {
      if (m_truthParticles.status.at(particle_i)==1) {
         if (!barcodeToParticle.insert( std::make_pair( m_truthParticles.barcode.at(particle_i), particle_i ) ).second) {
            std::cout << "WARNING barcode " <<  m_truthParticles.barcode.at(particle_i) << " of " << particle_i << " not unique." 
                      << std::endl;
         }
         particle_order.push_back( std::pair<float,unsigned int>(  sqr(m_truthParticles.px.at(particle_i))
                                                                 + sqr(m_truthParticles.py.at(particle_i)),
                                                                 particle_i));
      }
   }
   std::sort(particle_order.begin(), particle_order.end(),[](std::pair<float, unsigned int> &a, std::pair<float, unsigned int> &b) {
      return a.first > b.first;
   });

   unsigned int container_i=0;
   for (const SDOInfoContainer &sdo_info : m_sdoInfo ) {
      for( unsigned int meas_i=0; meas_i<sdo_info.sim_depositsBarcode->size(); ++meas_i) {
         for ( const std::vector<int>  &sim_hit_list : sdo_info.sim_depositsBarcode->at(meas_i)) {
            unsigned int idx1=0;
            for ( const int  &sim_barcode : sim_hit_list ) {
               std::map<int, unsigned int>::const_iterator particle_iter = barcodeToParticle.find(sim_barcode);
               if (particle_iter == barcodeToParticle.end()) continue;
               std::vector<std::pair<unsigned int, unsigned int> > &particle_hits = particleToHit[ particle_iter->second ];
               if (std::find(particle_hits.begin(),particle_hits.end(), std::make_pair(container_i, meas_i)) == particle_hits.end()) {
                  particle_hits.push_back(std::make_pair(container_i, meas_i));
               }
            }
         }
      }
      ++container_i;
   }
   Acts::Vector3 dummy_momentum{1.,1.,1.};
   for (std::pair<float, unsigned int>  &part_idx : particle_order ) {
      unsigned int part_i = part_idx.second;
      double norm_scale = 1./sqrt(part_idx.first + sqr(m_truthParticles.pz.at(part_i)));
      const std::vector<std::pair<unsigned int, unsigned int> > &particle_hits = particleToHit[ part_i ];
      if (particle_hits.size()>0) {
         std::cout << "ORIGPART " << m_truthParticles.pdgId.at(part_i)
                   << " " << m_truthParticles.barcode.at(part_i)
                   << " " << sqrt(part_idx.first)
                   << " " << m_truthParticles.px.at(part_i)*norm_scale
                   << " " << m_truthParticles.py.at(part_i)*norm_scale
                   << " " << m_truthParticles.pz.at(part_i)*norm_scale
                   << " " << particle_hits.size();




         for (const std::pair<unsigned int, unsigned int> &container_meas : particle_hits) {
            std::set<int> contributions;
            for (const std::vector<int> &sim_hit_list:
                    m_sdoInfo.at(container_meas.first).sim_depositsBarcode->at(container_meas.second) ) {
               for (int barcode : sim_hit_list) {

                  std::map<int, unsigned int>::const_iterator iter =  barcodeToParticle.find(barcode);
                  if (iter != barcodeToParticle.end() && m_truthParticles.status.at(iter->second)==1) {
                     int pdg_id = m_truthParticles.pdgId.at(iter->second);
                     if (pdg_id != 21 && pdg_id != 22) {
                        contributions.insert(barcode);
                     }
                  }
               }
            }
            if (contributions.size()==0) {
               std::cerr << "WARNING particle track without conributing stable, charged particles:" << std::endl;
               for (const std::vector<int> &sim_hit_list:
                       m_sdoInfo.at(container_meas.first).sim_depositsBarcode->at(container_meas.second) ) {
                  for (int barcode : sim_hit_list) {
                     
                     std::cerr << "WARNING contribution to " <<container_meas.first << ", " << container_meas.second
                               << " -> " << barcode << ":";
                     std::map<int, unsigned int>::const_iterator iter =  barcodeToParticle.find(barcode);
                     if (iter != barcodeToParticle.end()) {
                        std::cerr << " status=" << m_truthParticles.status.at(iter->second)
                                  << " pdgID=" << m_truthParticles.pdgId.at(iter->second);
                     }
                     std::cerr << std::endl;
                  }
               }
            }
            std::cout << " " << contributions.size();
            std::cout << " " << m_clusters.at(container_meas.first).globalX->at(container_meas.second)
                      << " " << m_clusters.at(container_meas.first).globalY->at(container_meas.second)
                      << " " << m_clusters.at(container_meas.first).globalZ->at(container_meas.second);

            std::vector<Acts::Vector2> corners;
            const Acts::Surface *surface=nullptr;
            std::unordered_map<unsigned long, Acts::GeometryIdentifier>::const_iterator
               geomap_iter = geo_map.find(m_clusters.at(container_meas.first).detectorElementID->at(container_meas.second));
            if (geomap_iter != geo_map.end()) {
               surface = trackingGeometry.findSurface(geomap_iter->second);
               if (surface) {
                  const Acts::AnnulusBounds *annulus_bounds = dynamic_cast<const Acts::AnnulusBounds *>( &surface->bounds() );
                  if (annulus_bounds) {
                     corners = annulus_bounds->corners();
                  }
                  else {
                     const Acts::PlanarBounds *bounds = dynamic_cast<const Acts::PlanarBounds *>( &surface->bounds() );
                     if (bounds) {
                        corners = bounds->vertices(1);
                     }
                  }
               }
            }
            std::cout << " " << corners.size();
            for (const Acts::Vector2 &a_corner  : corners) {
               assert( surface != nullptr);
               Acts::Vector3 glob = surface->localToGlobal(gctx,
                                                           a_corner,
                                                           dummy_momentum);
               std::cout << " " << glob[0] << " " << glob[1] << " " << glob[2];
            }
         }
         std::cout << std::endl;
      }
   }
}


void ActsExamples::RootDigiReader::convertMeasurements(const AlgorithmContext& ctx,
                                                       const std::map<int, ActsFatras::Barcode> &barcode_map) {
  std::list<IndexSourceLink> sourceLinkStorage;
  IndexSourceLinkContainer sourceLinks;
  MeasurementContainer measurements;
  IndexMultimap<ActsFatras::Barcode> measurementParticlesMap;
  std::unordered_map<unsigned long, Acts::GeometryIdentifier> geo_map;
  //  IndexMultimap<Index> measurementSimHitsMap;
  std::regex ignore(".*HGTD.*");

  const Acts::TrackingGeometry *tg=m_cfg.trackingGeometry.get();
  m_cfg.trackingGeometry->visitSurfaces([&geo_map,&ctx, tg](const Acts::Surface* surface) {
        if (surface && surface->associatedDetectorElement()) {
           const Test::IdentifiableDetectorElement *
              det_element=dynamic_cast<const Test::IdentifiableDetectorElement *>(surface->associatedDetectorElement());
           if (det_element) {
              std::pair<std::unordered_map<unsigned long, Acts::GeometryIdentifier>::iterator, bool>
                 insert_result =  geo_map.insert( std::make_pair( det_element->detectorId(), surface->geometryId()) );
              if (!insert_result.second) {
                 const Acts::Surface *old_surf = tg->findSurface(insert_result.first->second);
                 std::cout << "ERROR failed to insert ID mapping " << det_element->detectorId() << " -> " << surface->geometryId().value()
                           << " previously mapped to " << insert_result.first->second << std::endl
                           << " old : " << old_surf->geometryId().value() << " : "
                           << std::tuple<const Acts::Surface&,
                                         const Acts::GeometryContext&>(*old_surf,ctx.geoContext)
                           << " new : "<< surface->geometryId().value() << " : "
                           << std::tuple<const Acts::Surface&,
                                         const Acts::GeometryContext&>(*surface,ctx.geoContext)
                           <<std::endl;
              }
           }
           else {
              std::cout << "ERROR unsupported detector element " << typeid(*surface->associatedDetectorElement()).name()
                        << std::endl;
           }

        }
      });
  std::cout << "DEBUG RootDigiReader::convertMeasurements geo_map size from tg : " <<  geo_map.size() << std::endl;

  // m_cfg.trackingGeometry->visitSurfaces([&ctx, &ignore](const Acts::Surface* surface) {
  //       std::cout  << "DEBUG visit surfaces id " <<  ( surface->associatedLayer() ? surface->associatedLayer()->geometryId() : 0u)
  //                  << " "
  //                  <<  ( surface->associatedLayer()  && surface->associatedLayer()->trackingVolume()
  //                        ? surface->associatedLayer()->trackingVolume()->volumeName()
  //                        : "(unknown)" )
  //                  << std::endl;
  //       if (surface->associatedLayer()
  //           && surface->associatedLayer()->trackingVolume()
  //           && std::regex_match(surface->associatedLayer()->trackingVolume()->volumeName(), ignore)) return;
  //       const auto &trans = surface->transform(ctx.geoContext);
  //       for (unsigned int i=0; i< trans.rows(); ++i) {
  //          std::cout << i << ":";
  //          for (unsigned int j=0; j< trans.cols(); ++j) {
  //             std::cout << " " << std::setw(14) << trans(i,j);
  //          }
  //          std::cout << std::endl;
  //       }

  //    });

  if (false)  {
  for (const std::pair<const uint64_t, std::array<double,3 > > &elm : m_surfaceCenters) {
     std::cout <<"Find geoID for surface " << elm.first << " : " << elm.second[0] << ", " << elm.second[1] << ", " << elm.second[2]
               << std::endl;
           //        Acts::GeometryIdentifier geoId
           std::pair<Acts::GeometryIdentifier,bool> ret =  findGeoId(ctx,
                                                                     *(m_cfg.trackingGeometry),
                                                                     geo_map,
                                                                     elm.first,
                                                                     Acts::Vector3(elm.second[0],elm.second[1],elm.second[2]),
                                                                     1e-3,
                                                                     false);
  }
  }

  // {
  //    std::array<Acts::Vector3, 2>  pos {
  //          Acts::Vector3{642.409302, 126.368721, -1522.74902},
  //          Acts::Vector3{642.409302, 126.368721, -1522.74902},
  //          };
  //    std::array<uint64_t,2> geo_ids {1441152705392445696, 1441152705392443648 };
  //    for(const Acts::Vector3 &a_pos : pos) {
  //       for (uint64_t geo_id : geo_ids) {
  //          testSurface(ctx, m_cfg.trackingGeometry->findSurface(geo_id), a_pos, 1e-2);
  //       }
  //    }
  // }
  for (const std::pair<const int,  ActsFatras::Barcode> &elm : barcode_map ) {
     std::cout << "DEBUG convert measruements with hits from " << elm.first << " -> " << elm.second << std::endl;
  }
  unsigned int container_i=0;
  std::vector<unsigned int> measurements_per_container;
  measurements_per_container.reserve(m_clusters.size());
  for (const ClusterContainer &cluster_container : m_clusters ) {
     cluster_container.checkDimensions();
     const SDOInfoContainer &sdo_info = m_sdoInfo.at(container_i);
     sdo_info.checkDimensions(cluster_container.localX->size());
     std::vector< unsigned int > mismatch;
     unsigned int n_measurements_per_container=0;
     for (unsigned int meas_i=0; meas_i < cluster_container.localX->size(); ++meas_i) {
        if (m_cfg.singleParticle) {
           bool accept=false;
           for ( const std::vector<int>  &sim_hit_list : sdo_info.sim_depositsBarcode->at(meas_i) ) {
              for ( const int  &sim_barcode : sim_hit_list ) {
                 if (barcode_map.find(sim_barcode) != barcode_map.end()) {
                    accept=true;
                    break;
                 }
              }
              if (accept) break;
           }
           if (!accept) continue;
           std::cout << "DEBUG convert measurement " << meas_i << std::endl;
        }
        Index measurementIdx = measurements.size();
        //        Acts::GeometryIdentifier geoId
        std::pair<Acts::GeometryIdentifier,bool> ret =  findGeoId(ctx,
                                                    *(m_cfg.trackingGeometry),
                                                    geo_map,
                                                    cluster_container.detectorElementID->at(meas_i),
                                                    Acts::Vector3(cluster_container.globalX->at(meas_i),
                                                                  cluster_container.globalY->at(meas_i),
                                                                  cluster_container.globalZ->at(meas_i)),
                                                    std::max(cluster_container.localXError->at(meas_i),
                                                             cluster_container.localYError->at(meas_i)));
        Acts::GeometryIdentifier geoId = ret.first;
        if (geoId.value() == 576460889742376960) {
           auto surface = m_cfg.trackingGeometry->findSurface(geoId);
           if (surface) {
              std::cout << "DEBUG search surface matching " << geoId.value() << " " << geoId << " surface " << surface->geometryId().value()
                        <<  surface->geometryId() << " " << std::tuple<const Acts::Surface&,
                                                                       const Acts::GeometryContext&>(*surface,ctx.geoContext)
                        << std::endl;
           }
           else {
              std::cout << "DEBUG search surface matching " << geoId.value() << " " << geoId << " NO surface." << std::endl;
           }
        }

        if (!m_cfg.trackingGeometry->findSurface(geoId)) {
           std::stringstream msg;
           msg << "ERROR invalid geoID for "
               << cluster_container.globalX->at(meas_i)
               << cluster_container.globalY->at(meas_i)
               << cluster_container.globalZ->at(meas_i)
               << " det=" << cluster_container.detectorElementID->at(meas_i)
               << " -> " << geoId.value() << " " << geoId;
           std::runtime_error(msg.str());
        }
        fillValidationHists(Acts::Vector3(cluster_container.globalX->at(meas_i),
                                          cluster_container.globalY->at(meas_i),
                                          cluster_container.globalZ->at(meas_i)),
                            ret.second ? kMatched : kUnmatched);
        if (!ret.second) mismatch.push_back(meas_i);

        sourceLinkStorage.emplace_back(geoId, measurementIdx);
        IndexSourceLink& sourceLink = sourceLinkStorage.back();

        std::cout << "DEBUG add measurement " << geoId << " -> " << measurementIdx << std::endl;

        sourceLinks.insert(sourceLink);

        // @TODO currently only 2D measurements are supported.
        std::array<Acts::BoundIndices, 2> indices {Acts::BoundIndices(0u),Acts::BoundIndices(1u)};
        Acts::ActsVector<2>    par { cluster_container.localX->at(meas_i),
                                     cluster_container.localY->at(meas_i)} ;
        Acts::ActsSymMatrix<2> cov;
        cov(0,0) = cluster_container.localXError->at(meas_i);
        cov(1,1) = cluster_container.localYError->at(meas_i);
        cov(0,1) = cluster_container.localXYCorrelation->at(meas_i);
        cov(1,0) = cluster_container.localXYCorrelation->at(meas_i);

        {
           auto surface = m_cfg.trackingGeometry->findSurface(geoId);
           Acts::Vector3 globalPos_converted = surface->localToGlobal(ctx.geoContext, par, Acts::Vector3());
           Acts::Vector3 globalPos_orig( cluster_container.globalX->at(meas_i),
                                         cluster_container.globalY->at(meas_i),
                                         cluster_container.globalZ->at(meas_i));
           Acts::Vector3 dummy_momentum{1.,1.,1.};
           Acts::Result<Acts::Vector2> localPos_converted = surface->globalToLocal(ctx.geoContext,
                                                                                   globalPos_orig,
                                                                                   dummy_momentum,
                                                                                   1e-3);

           std::cout << "DEBUG convert measurement " << meas_i << " coord : "
                     << par[0] << ", " << par[1] << "[ " << cov(0,0) << ", " << cov(1,1) << "; " << cov(0,1) << "] "
                     << " | " << cluster_container.globalX->at(meas_i)
                     << ", " << cluster_container.globalY->at(meas_i)
                     << ", " << cluster_container.globalZ->at(meas_i)
                     << std::endl
                     << " -> " << localPos_converted.value()[0] << ", " << localPos_converted.value()[1]
                     << (localPos_converted.ok() ? " OK  " : " ERR ")
                     << " -> " << globalPos_converted[0]
                     << ", " << globalPos_converted[1]
                     << ", " << globalPos_converted[2]
                     << std::endl
                     << cornersToStr(surface, ctx.geoContext)
                     << std::endl;
           if (   std::abs(cluster_container.globalX->at(meas_i)-globalPos_converted[0])>1e-3
               || std::abs(cluster_container.globalY->at(meas_i)-globalPos_converted[1])>1e-3
               || std::abs(cluster_container.globalZ->at(meas_i)-globalPos_converted[2])>1e-3) {
              std::cout << "ERROR global coordinate mismatch for measurement " << meas_i << " " << container_i
                        << " : deltas "
                        << (cluster_container.globalX->at(meas_i)-globalPos_converted[0])
                        << ", " << (cluster_container.globalY->at(meas_i)-globalPos_converted[1])
                        << ", " <<(cluster_container.globalZ->at(meas_i)-globalPos_converted[2])
                        << std::endl;
              double par0=par[0];
              par[0]=par[1];
              par[1]=par0;
              globalPos_converted = surface->localToGlobal(ctx.geoContext, par, Acts::Vector3());
              std::cout << "DEBUG try swapped local coordinates: "
                     << " -> " << globalPos_converted[0]
                        << ", " << globalPos_converted[1]
                        << ", " << globalPos_converted[2]
                        << " deltas : "
                        << (cluster_container.globalX->at(meas_i)-globalPos_converted[0])
                        << ", " << (cluster_container.globalY->at(meas_i)-globalPos_converted[1])
                        << ", " <<(cluster_container.globalZ->at(meas_i)-globalPos_converted[2])
                        << std::endl;

              par[0] = localPos_converted.value()[0];
              par[1] = localPos_converted.value()[1];
              double sig_r = cov(0,0);
              cov(0,0) = cov(1,1);
              cov(1,1) = sig_r;
           }
        }

        switch (container_i) {
           case 0: /* pixel */
              measurements.emplace_back( Acts::Measurement<Acts::BoundIndices, 2>(Acts::SourceLink{geoId, sourceLink},
                                                                                  indices, par, cov) );
              break;
           case 1: /* strips */
              {

                  std::array<Acts::BoundIndices, 1> strip_indices{indices[0]};
                  Acts::ActsVector<1> strip_par {par[0]};
                  Acts::ActsSymMatrix<1> strip_cov {cov(0,0) };

                 measurements.emplace_back( Acts::Measurement<Acts::BoundIndices, 1>(Acts::SourceLink{geoId, sourceLink},
                                                                                     strip_indices, strip_par, strip_cov) );
                 break;
              }
           default :
              throw std::runtime_error("Currently assuming container_i==0 -> pixel; container_i==1 -> strips. Nothing else supported.");
           }
        ++n_measurements_per_container;
        unsigned int n_contributions=0;
        unsigned int n_contributions_validBarcode=0;
        double       energy=0.;
        double       energy_validBarcode=0.;
        unsigned int idx0=0;
        const std::vector<std::vector<float> > &sim_hit_energy_list =sdo_info.sim_depositsEnergy->at(meas_i);

        for ( const std::vector<int>  &sim_hit_list : sdo_info.sim_depositsBarcode->at(meas_i) ) {
           unsigned int idx1=0;
           for ( const int  &sim_barcode : sim_hit_list ) {
              bool has_barcode=false;
              try {
                 ActsFatras::Barcode particle_id = findParticleBarCode(barcode_map, sim_barcode);
                 std::pair<HitParticlesMap::const_iterator, HitParticlesMap::const_iterator>
                    particle_range = measurementParticlesMap.equal_range(meas_i);
                 for (auto hitParticle : makeRange(particle_range)) {
                    if (hitParticle.second == particle_id) {
                       has_barcode = true;
                       break;
                    }
                 }
                 std::cout << "DEBUG sim barcode " << sim_barcode << " -> " << particle_id
                           << " | " << (has_barcode ? " has barcode" : "")
                           << std::endl;
                 if (!has_barcode) {
                    measurementParticlesMap.emplace_hint( (particle_range.first == particle_range.second
                                                           ? measurementParticlesMap.end()
                                                           : particle_range.second),
                                                         measurementIdx,
                                                         particle_id);
                 }
              }
              catch (std::runtime_error &err) {
                 std::cout << "DEBUG missing target barcode for " << measurementIdx << " source=" << sim_barcode
                           << " : " << err.what()
                           << std::endl;
              }

              {
                 if (idx0<sim_hit_energy_list.size() && idx1<sim_hit_energy_list[idx0].size()) {
                    if (has_barcode &&  sim_barcode>0 && sim_barcode<std::numeric_limits<int>::max()-1) {
                       ++n_contributions_validBarcode;
                       energy_validBarcode += sim_hit_energy_list[idx0][idx1];
                    }
                    else {
                       energy += sim_hit_energy_list[idx0][idx1];
                       ++n_contributions;
                    }
                 }
                 else {
                    std::cout << "ERROR depostits: barcode and energy dimensions mismatch: "
                              << idx0 << ", " << idx1
                              << std::endl;
                 }
              }
              ++idx1;
           }
           ++idx0;
        }
        energy += energy_validBarcode;
        m_assoStat.at(kEnergyFraction).add(energy > 0 ? energy_validBarcode / energy : 0.);
        m_assoStat.at(kEnergy).add(energy);
        m_assoStat.at(kContributionsValidBarcode).add(n_contributions_validBarcode);
        m_assoStat.at(kContributions).add(n_contributions_validBarcode + n_contributions);

     }

     for (unsigned int meas_i  : mismatch) {
        Index measurementIdx = measurements.size();
        std::cout << "DEBUG search again " << meas_i << " : " << cluster_container.detectorElementID->at(meas_i)
                  << std::endl;
        std::map<uint64_t, std::array<double,3 > >::const_iterator 
           iter = m_surfaceCenters.find(cluster_container.detectorElementID->at(meas_i));
        if (iter != m_surfaceCenters.end()) {
           //        Acts::GeometryIdentifier geoId
           std::pair<Acts::GeometryIdentifier,bool> ret =  findGeoId(ctx,
                                                                     *(m_cfg.trackingGeometry),
                                                                     geo_map,
                                                                     cluster_container.detectorElementID->at(meas_i),
                                                                     Acts::Vector3(iter->second[0],iter->second[1],iter->second[2]),
                                                                     std::max(cluster_container.localXError->at(meas_i),
                                                                              cluster_container.localYError->at(meas_i)));
        }

     }
     measurements_per_container.push_back( n_measurements_per_container);
     
     ++container_i;
  }
  {
  container_i=0;
  for (unsigned int value : measurements_per_container) {
     std::cout << "DEBUG measurements per container " << container_i << " : " << value << std::endl;
     ++container_i;
  }
  }
  checkDestinationMultiplicity();

  if (m_cfg.stop) {
     //     std::cout << "DEBUG " << getpid() << " wait for signal CONT " << std::endl;
     //     kill(getpid(), SIGSTOP);
  }


  dumpParticleHits(ctx.geoContext,
                   *(m_cfg.trackingGeometry),
                   geo_map);

  ctx.eventStore.add(m_cfg.sourceLinks, std::move(sourceLinks));
  ctx.eventStore.add(m_cfg.sourceLinks + "__storage",
                     std::move(sourceLinkStorage));
  ctx.eventStore.add(m_cfg.measurements, std::move(measurements));
  ctx.eventStore.add(m_cfg.measurementParticlesMap,
                     std::move(measurementParticlesMap));
}


void ActsExamples::RootDigiReader::showStat() const {
   for (const std::pair< std::string, std::vector<std::pair<std::string, Stat>  > > &stat_list :  m_stat ) {
      for (const std::pair<std::string, Stat> &a_stat :  stat_list.second ) {
         std::cout << stat_list.first << " " << a_stat.first << " : " << a_stat.second << std::endl;
      }
   }
   static const std::array<const char *,kNAssocStat>
      label{"Energy", "Energy fraction w barcode","contributions", "Contributions w barcode"};
   unsigned idx=0;
   for (const Stat &a_stat :  m_assoStat) {
      std::cout << label.at(idx) << " : " << a_stat << std::endl;
      ++idx;
   }

}
