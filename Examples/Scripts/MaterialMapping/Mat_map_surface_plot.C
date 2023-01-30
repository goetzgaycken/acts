// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <TROOT.h>
#include "TH2F.h"
#include "TText.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TSystem.h"

#include "materialPlotHelper.cpp"

#include <fstream>
#include <iostream>
#include <sstream>

#include <map>
#include <array>
#include <limits>
#include <cmath>
/// Draw and save the histograms.
class MyStat {
public:
   void add(double value) {
      ++m_n;
      m_sum += value;
      m_sum2 += value*value;
      m_min=std::min(m_min,value);
      m_max=std::max(m_max,value);
   }
   unsigned int n() const { return m_n; } 
   double max() const { return m_max; }
   double min() const { return m_min; }
   double mean() const { return (m_n>0 ? m_sum/m_n : 0 ); }
   double rms() const { return (m_n>1 ? sqrt((m_sum2 - m_sum*m_sum/m_n)/(m_n-1))+1 : 0 ); }
   
   unsigned int m_n=0;
   double m_sum =0.;
   double m_sum2 = 0;
   double m_min = std::numeric_limits<double>::max();
   double m_max = -std::numeric_limits<double>::max();
};
class RangeKey {
public:
   RangeKey(float a_min, float a_max) : m_min(a_min), m_max(a_max) {}
   bool operator<(const RangeKey &a) const {
      return m_min < a.m_min;
   }
   float min() const { return m_min; }
   float max() const { return m_max; }
   bool contained( float val ) const {
      return m_min <= val && val <= m_max;
   }
   float m_min;
   float m_max;
};

using Map3D = std::multimap<RangeKey, std::multimap<RangeKey, std::multimap< RangeKey , std::pair< std::array<float,3>, uint64_t> > > >;

class SurfaceIDMap {
public:
   std::unordered_map<uint64_t, std::array<MyStat, 3>  > m_geoIdToCenter;
   void addSurface(uint64_t id, float x, float y, float z) {
      std::pair<std::unordered_map<uint64_t, std::array<MyStat, 3>  >::iterator, bool>
         ret = m_geoIdToCenter.insert( std::make_pair( id, std::array<MyStat, 3>{} ) );
      ret.first->second.at(0).add(x);
      ret.first->second.at(1).add(y);
      ret.first->second.at(2).add(z);
   }
   Map3D getCenterToGeoIDMap() {
      Map3D  center_to_geo_id; 
      for(const auto &elm : m_geoIdToCenter ) {
         float lower_x = static_cast<float>(elm.second.at(0).min()-2*elm.second.at(0).rms());
         float lower_y = static_cast<float>(elm.second.at(1).min()-2*elm.second.at(1).rms());
         float lower_z = static_cast<float>(elm.second.at(2).min()-2*elm.second.at(2).rms());
         float mean_x = static_cast<float>(elm.second.at(0).mean());
         float mean_y = static_cast<float>(elm.second.at(1).mean());
         float mean_z = static_cast<float>(elm.second.at(2).mean());
         bool next=false;
         std::multimap<RangeKey, std::multimap<RangeKey, std::multimap< RangeKey , std::pair< std::array<float,3>, uint64_t> > > >::iterator
            x_iter = center_to_geo_id.lower_bound(RangeKey(lower_x,lower_x));
         for ( ;x_iter != center_to_geo_id.end() && x_iter->first.min() <= elm.second.at(0).min(); ++x_iter) {
            if (x_iter->first.contained(mean_x)) {
               std::multimap<RangeKey, std::multimap< RangeKey , std::pair< std::array<float,3>, uint64_t> > >::iterator
                  y_iter = x_iter->second.lower_bound(RangeKey(lower_y,lower_y));
               for (;y_iter != x_iter->second.end() && y_iter->first.min() <= elm.second.at(1).min(); ++y_iter) {
                  if (y_iter->first.contained(mean_y)) {
                     y_iter->second.insert( std::make_pair(RangeKey(lower_z, static_cast<float>(elm.second.at(2).max()+2*elm.second.at(2).rms())),
                                                           std::make_pair( std::array<float,3>{mean_x, mean_y, mean_z}, elm.first )
                                                           ) );
                     std::cout << "Insert: "
                               << x_iter->first.min() << " .. "
                               << x_iter->first.max() << " , "
                               << y_iter->first.min() << " .. "
                               << y_iter->first.max() << " : "
                               << lower_x << ", " << lower_y << ", " << lower_z << " [ " << mean_x << " +- " << elm.second.at(0).rms()
                               << ", " << mean_y << " +- " << elm.second.at(1).rms()
                               << ", " << mean_z << " +- " << elm.second.at(2).rms()
                               << " ]"
                               << elm.first << std::endl;
                     std::cout << "       "
                               << elm.second.at(0).min() << " .. " << elm.second.at(0).max()
                               << elm.second.at(1).min() << " .. " << elm.second.at(1).max()
                               << elm.second.at(2).min() << " .. " << elm.second.at(2).max()
                               << std::endl;
                     
                     next=true;
                     break;
               // std::multimap< RangeKey , std::pair< std::array<float,3>, uint64_t> >::iterator
               //    z_iter = center_to_geo_id.lower_bound(lower_y);
               // for (;z_iter != center_to_geo_id.end() && z_iter->first.min() <= elm.second.at(2).min(); ++y_iter) {
               //    if (z_iter->first.contained( mean_z)) {
               //       z_iter->second
               //    }
               // }
                  }
               }
               if (next) break;
               x_iter->second.insert(std::make_pair(RangeKey(lower_y, static_cast<float>(elm.second.at(1).max()+2*elm.second.at(1).rms())),
                                                    std::multimap< RangeKey , std::pair< std::array<float,3>, uint64_t> >()))
                  ->second.insert(std::make_pair(RangeKey(lower_z, static_cast<float>(elm.second.at(2).max()+2*elm.second.at(2).rms())),
                                                       std::make_pair( std::array<float,3>{mean_x, mean_y, mean_z}, elm.first )
                                                       ) );
               std::cout << "Insert: "
                         << x_iter->first.min() << " .. "
                         << x_iter->first.max() << " : "
                         << lower_x << ", " << lower_y << ", " << lower_z << " [ " << mean_x << " +- " << elm.second.at(0).rms()
                         << ", " << mean_y << " +- " << elm.second.at(1).rms()
                         << ", " << mean_z << " +- " << elm.second.at(2).rms()
                         << " ]"
                         << elm.first << std::endl;
               std::cout << "       "
                         << elm.second.at(0).min() << " .. " << elm.second.at(0).max()
                         << elm.second.at(1).min() << " .. " << elm.second.at(1).max()
                         << elm.second.at(2).min() << " .. " << elm.second.at(2).max()
                         << std::endl;
               next=true;
            
            }
         }
         if (next) continue;
         center_to_geo_id.insert(std::make_pair(RangeKey(lower_x, static_cast<float>(elm.second.at(0).max()+2*elm.second.at(0).rms())),
                                                std::multimap<RangeKey, std::multimap< RangeKey , std::pair< std::array<float,3>, uint64_t> > >()))
            ->second.insert(std::make_pair(RangeKey(lower_y, static_cast<float>(elm.second.at(1).max()+2*elm.second.at(1).rms())),
                                                 std::multimap< RangeKey , std::pair< std::array<float,3>, uint64_t> >()))
            ->second.insert(std::make_pair(RangeKey(lower_z, static_cast<float>(elm.second.at(2).max()+2*elm.second.at(2).rms())),
                                          std::make_pair( std::array<float,3>{mean_x, mean_y, mean_z}, elm.first )
                                          ) );
         std::cout << "Insert: "
                   << lower_x << ", " << lower_y << ", " << lower_z << " [ " << mean_x << " +- " << elm.second.at(0).rms()
                   << ", " << mean_y << " +- " << elm.second.at(1).rms()
                   << ", " << mean_z << " +- " << elm.second.at(2).rms()
                   << " ]"
                   << elm.first << std::endl;
         std::cout << "       "
                   << elm.second.at(0).min() << " .. " << elm.second.at(0).max()
                   << elm.second.at(1).min() << " .. " << elm.second.at(1).max()
                   << elm.second.at(2).min() << " .. " << elm.second.at(2).max()
                   << std::endl;
      }
      return center_to_geo_id;
   }
};

inline double sqr( double a) { return a*a; }
inline double distance2( float x, float y, float z, const std::array<float, 3> &point) {
   return sqr( x - point.at(0) ) + sqr( y - point.at(1) ) + sqr( z - point.at(2) );
}
uint64_t findGeoID(float x_val, float y_val, float z_val,
                   Map3D &center_to_geo_id) {
   double min_dist = std::numeric_limits<double>::max();
   uint64_t best_id = 0u;
   unsigned int n_x=0u;
   unsigned int n_y=0u;
   unsigned int n_z=0u;
   unsigned int n_tests=0u;
   std::multimap<RangeKey, std::multimap<RangeKey, std::multimap< RangeKey , std::pair< std::array<float,3>, uint64_t> > > >::iterator
      x_iter = center_to_geo_id.lower_bound(RangeKey(x_val,x_val));
   for ( ;x_iter != center_to_geo_id.end() && x_iter->first.min() <= x_val; ++x_iter) {
      ++n_x;
      if (x_iter->first.contained(x_val)) {
         std::multimap<RangeKey, std::multimap< RangeKey , std::pair< std::array<float,3>, uint64_t> > >::iterator
            y_iter = x_iter->second.lower_bound(RangeKey(y_val,y_val));
         for (;y_iter != x_iter->second.end() && y_iter->first.min() <= y_val; ++y_iter) {
            ++n_y;
            if (y_iter->first.contained(y_val)) {
               std::multimap< RangeKey , std::pair< std::array<float,3>, uint64_t> >::iterator
                  z_iter = y_iter->second.lower_bound(RangeKey(z_val,z_val));
               for (;z_iter != y_iter->second.end() && y_iter->first.min() <= y_val; ++y_iter) {
                  ++n_z;
                  if (z_iter->first.contained(z_val)) {
                     ++n_tests;
                     double dist = distance2(x_val, y_val, z_val, z_iter->second.first);
                     if (dist < min_dist) {
                        min_dist=dist;
                        best_id = z_iter->second.second;
                     }
                  }
               }
            }
         }
      }
   }
   if (best_id == 0u) {
      std::cerr << "ERROR failed to find geo ID for " << x_val << ", " << y_val  << ", " << z_val << " tests: "
                << n_x << ", " << n_y << ", " << n_z << " -> "
                << n_tests << std::endl;
   }
   return best_id;
}
   

void plot(std::vector<TH2F*> Map, const sinfo& surface_info, const std::string& name){

  std::string out_name = name+"/"+surface_info.name+"/"+surface_info.name+"_"+surface_info.idname;
  gSystem->Exec( Form("mkdir %s", (name+"/"+surface_info.name).c_str()) );

  // Disk
  if(surface_info.type == 2  || surface_info.type == 4){

    TText *vol = new TText(.1,.95,surface_info.name.c_str());
    vol->SetNDC();
    TText *surface = new TText(.1,.9,surface_info.id.c_str());
    surface->SetNDC();
    TText *surface_z = new TText(.1,.85,("Z = " + to_string(surface_info.pos)).c_str() );
    surface_z->SetNDC();

    TCanvas *c1 = new TCanvas("c1","mat_X0",1200,1200);
    c1->SetRightMargin(0.14);
    c1->SetTopMargin(0.14);
    c1->SetLeftMargin(0.14);
    c1->SetBottomMargin(0.14);
    Map[0]->Draw("COLZ");
    vol->Draw();
    surface->Draw();
    surface_z->Draw();
    c1->Print( (out_name+"_X0.pdf").c_str());
    //c1->Print( (out_name+"_X0.root").c_str());

    delete c1;

    delete vol;
    delete surface;
    delete surface_z;
  }

  // Cylinder
  if(surface_info.type == 1){

    TText *vol = new TText(.1,.95,surface_info.name.c_str());
    vol->SetNDC();
    TText *surface = new TText(.1,.9,surface_info.id.c_str());
    surface->SetNDC();
    TText *surface_r = new TText(.1,.85,("R = " + to_string(surface_info.pos)).c_str() );
    surface_r->SetNDC();
    TCanvas *c1 = new TCanvas("c1","mat_X0",1200,1200);
    c1->SetRightMargin(0.14);
    c1->SetTopMargin(0.14);
    c1->SetLeftMargin(0.14);
    c1->SetBottomMargin(0.14);
    Map[0]->Draw("COLZ");
    vol->Draw();
    surface->Draw();
    surface_r->Draw();
    c1->Print( (out_name+"_X0.pdf").c_str());
    //c1->Print( (out_name+"_X0.root").c_str());

    delete c1;

    delete vol;
    delete surface;
    delete surface_r;
  }
  return;
}

/// Initialise the histograms for each surface.

void Initialise_hist(std::vector<TH2F*>& surface_hist,
  const sinfo& surface_info){

  TH2F * Map_X0;
  TH2F * Map_L0;

  TH2F * Map_scale;

  if(surface_info.type == 1){
    Map_X0    = new TH2F(("Map_X0_"+surface_info.idname).c_str(),("Map_X0_"+surface_info.idname).c_str(),
                         50,-6,6,50,-3.2,3.2);
    Map_L0    = new TH2F(("Map_L0_"+surface_info.idname).c_str(),("Map_L0_"+surface_info.idname).c_str(),
                         50,-6,6,50,-3.2,3.2);
    Map_scale = new TH2F(("Map_scale_"+surface_info.idname).c_str(),("Map_scale_"+surface_info.idname).c_str(),
                          50,-6,6,50,-3.2,3.2);
    Map_X0->GetXaxis()->SetTitle("Eta");
    Map_X0->GetYaxis()->SetTitle("Phi");
    Map_X0->GetZaxis()->SetTitle("X0");
    Map_L0->GetXaxis()->SetTitle("Eta");
    Map_L0->GetYaxis()->SetTitle("Phi");
    Map_L0->GetZaxis()->SetTitle("L0");
  }

  if(surface_info.type == 2 || surface_info.type == 4){
    Map_X0    = new TH2F(("Map_X0_"+surface_info.idname).c_str(),("Map_X0_"+surface_info.idname).c_str(),
                          50,-1*surface_info.range_max,surface_info.range_max,50,-1*surface_info.range_max, surface_info.range_max);
    Map_L0    = new TH2F(("Map_L0_"+surface_info.idname).c_str(),("Map_L0_"+surface_info.idname).c_str(),
                          50,-1*surface_info.range_max,surface_info.range_max,50,-1*surface_info.range_max, surface_info.range_max);
    Map_scale = new TH2F(("Map_scale_"+surface_info.idname).c_str(),("Map_scale_"+surface_info.idname).c_str(),
                          50,-1*surface_info.range_max,surface_info.range_max,50,-1*surface_info.range_max, surface_info.range_max);
    Map_X0->GetXaxis()->SetTitle("X [mm]");
    Map_X0->GetYaxis()->SetTitle("Y [mm]");
    Map_X0->GetZaxis()->SetTitle("X0");
    Map_L0->GetXaxis()->SetTitle("X [mm]");
    Map_L0->GetYaxis()->SetTitle("Y [mm]");
    Map_L0->GetZaxis()->SetTitle("L0");
  }
  std::vector<TH2F*> v_hist;
  v_hist.push_back(Map_X0);
  v_hist.push_back(Map_L0);
  v_hist.push_back(Map_scale);
  surface_hist = v_hist;
}

/// Fill the histograms for each surfaces.

void Fill(std::map<uint64_t,std::vector<TH2F*>>& surface_hist,  std::map<uint64_t,sinfo>& surface_info,
          const std::string& input_file, const int& nbprocess, SurfaceIDMap *id_map=nullptr, Map3D *id_remap=nullptr){

  std::map<std::string,std::string> surface_name;

  std::map<uint64_t,float> surface_weight;

  //Get file, tree and set top branch address
  TFile *tfile = new TFile(input_file.c_str());
  TTree *tree = (TTree*)tfile->Get("material-tracks");

  float v_phi   = 0;
  float v_eta   = 0;
  std::vector<float> *mat_X0   = 0;
  std::vector<float> *mat_L0   = 0;
  std::vector<float> *mat_step_length = 0;

  std::vector<uint64_t> *sur_id = 0;
  std::vector<int32_t> *sur_type = 0;
  std::vector<float> *sur_x = 0;
  std::vector<float> *sur_y = 0;
  std::vector<float> *sur_z = 0;
  std::vector<float> *sur_range_min = 0;
  std::vector<float> *sur_range_max = 0;

  tree->SetBranchAddress("v_phi",&v_phi);
  tree->SetBranchAddress("v_eta",&v_eta);
  tree->SetBranchAddress("mat_X0",&mat_X0);
  tree->SetBranchAddress("mat_L0",&mat_L0);
  tree->SetBranchAddress("mat_step_length",&mat_step_length);

  tree->SetBranchAddress("sur_id",&sur_id);
  tree->SetBranchAddress("sur_type",&sur_type);
  tree->SetBranchAddress("sur_x",&sur_x);
  tree->SetBranchAddress("sur_y",&sur_y);
  tree->SetBranchAddress("sur_z",&sur_z);
  tree->SetBranchAddress("sur_range_min",&sur_range_min);
  tree->SetBranchAddress("sur_range_max",&sur_range_max);

  int nentries = tree->GetEntries();
  if(nentries > nbprocess && nbprocess != -1) nentries = nbprocess;
  // Loop over all the material tracks.
  for (Long64_t i=0;i<nentries; i++) {
    if(i%10000==0) std::cout << "processed " << i << " events out of " << nentries << std::endl;
    tree->GetEntry(i);

    // Reset the weight
    for (auto weight_it = surface_weight.begin(); weight_it != surface_weight.end(); weight_it++){
      weight_it->second = 0;
    }
    // loop over all the material hits to do initialisation and compute weight
    for(int j=0; j<mat_X0->size(); j++ ){

      // Ignore surface of incorrect type
      if(sur_type->at(j) == -1) continue;

      if (id_map) {
         id_map->addSurface(sur_id->at(j),sur_x->at(j),sur_y->at(j),sur_z->at(j));
      }
      uint64_t the_sur_id= sur_id->at(j);
      if (id_remap) {
         the_sur_id =  findGeoID(sur_x->at(j),sur_y->at(j),sur_z->at(j), *id_remap);
      }

      // If a surface was never encountered initialise the hist, info and weight
      if(surface_hist.find(the_sur_id)==surface_hist.end()){
        float pos;
        if(sur_type->at(j) == 1){
          pos = sqrt(sur_x->at(j)*sur_x->at(j)+sur_y->at(j)*sur_y->at(j));
        }
        if(sur_type->at(j) == 2 || sur_type->at(j) == 4){
          pos = sur_z->at(j);
        }
        surface_weight[the_sur_id] = 0;
        Initialise_info(surface_info[the_sur_id], surface_name, the_sur_id, sur_type->at(j), pos, sur_range_min->at(j), sur_range_max->at(j));
        Initialise_hist(surface_hist[the_sur_id], surface_info[the_sur_id]);
      }
      // Weight for each surface = number of hit associated to it.
      surface_weight[the_sur_id]++;
    }

    // loop over all the material hit to fill the histogram
    for(int j=0; j<mat_X0->size(); j++ ){

      // Ignore surface of incorrect type
      if(sur_type->at(j) == -1) continue;
      uint64_t the_sur_id= sur_id->at(j);
      if (id_remap) {
         the_sur_id =  findGeoID(sur_x->at(j),sur_y->at(j),sur_z->at(j), *id_remap);
      }

      if(sur_type->at(j) == 1){
        surface_hist[the_sur_id][0]->Fill(v_eta, v_phi, (mat_step_length->at(j)/mat_X0->at(j)));
        surface_hist[the_sur_id][1]->Fill(v_eta, v_phi, (mat_step_length->at(j)/mat_L0->at(j)));
        surface_hist[the_sur_id][2]->Fill(v_eta, v_phi, (1/surface_weight[the_sur_id]));
      }
      if(sur_type->at(j) == 2 || sur_type->at(j) == 4){
        surface_hist[the_sur_id][0]->Fill(sur_x->at(j), sur_y->at(j), (mat_step_length->at(j)/mat_X0->at(j)));
        surface_hist[the_sur_id][1]->Fill(sur_x->at(j), sur_y->at(j), (mat_step_length->at(j)/mat_L0->at(j)));
        surface_hist[the_sur_id][2]->Fill(sur_x->at(j), sur_y->at(j), (1/surface_weight[the_sur_id]));
      }
    }
  }
  // Normalise the histograms
  for (auto hist_it = surface_hist.begin(); hist_it != surface_hist.end(); hist_it++){
    hist_it->second[0]->Divide(hist_it->second[2]);
    hist_it->second[1]->Divide(hist_it->second[2]);
  }
}

/// Plot the material on each surface.
/// nbprocess : number of parameter to be processed.
/// name : name of the output directory.

void Mat_map_surface_plot(std::string input_file = "", int nbprocess = -1, std::string name = ""){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  std::map<uint64_t,std::vector<TH2F*>> surface_hist;
  std::map<uint64_t,sinfo> surface_info;

  Fill(surface_hist, surface_info, input_file, nbprocess);
  for (auto hist_it = surface_hist.begin(); hist_it != surface_hist.end(); hist_it++){
    plot(hist_it->second, surface_info[hist_it->first], name);
    for (auto hist : hist_it->second){
      delete hist;
    }
    hist_it->second.clear();
  }
}
