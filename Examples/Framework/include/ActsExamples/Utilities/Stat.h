#pragma once

#include <cmath>
#include <algorithm>
#include <vector>
#include <ostream>

namespace ActsExamples {

   inline double sqr(double a) { return a*a; }
   class Stat {
   public:
      Stat() = default;
      Stat(Stat &&) = default;
      Stat(const Stat &) = default;
      Stat(unsigned int nbins, double lower_edge, double upper_edge) {
         setBinning(nbins, lower_edge, upper_edge);
      }
      void add(double val) {
         m_n += 1.;
         m_sum += val;
         m_sum2 += val*val;
         m_min=std::min(val,m_min);
         m_max=std::max(val,m_max);
         if (!m_bins.empty()) {fill(val); }
      }
      double mean() const { return m_n>0. ? m_sum/m_n : 0.; }
      double rms()  const { return m_n>1. ? sqrt( (m_sum2 - m_sum*m_sum/m_n)/(m_n-1.)) : 0.; }
      double n()    const { return m_n; }
      double min()  const { return m_min; }
      double max()  const { return m_max; }

      Stat &operator =(const Stat &a) = default;
      Stat &operator =(Stat &&a) = default;
      Stat &operator +=(const Stat &a) {
         m_n += a.m_n;
         m_sum += a.m_sum;
         m_sum2 += a.m_sum2;
         m_min=std::min(a.m_min,m_min);
         m_max=std::max(a.m_max,m_max);
         if (m_bins.size() == a.m_bins.size()) {
            for (unsigned int bin_i=0; bin_i < m_bins.size(); ++bin_i) {
               m_bins[bin_i] += a.m_bins[bin_i];
            }
         }
         return *this;
      }
      double m_n   = 0.;
      double m_sum = 0.;
      double m_sum2 = 0.;
      double m_min = std::numeric_limits<double>::max();
      double m_max = -std::numeric_limits<double>::max();

      void setBinning(unsigned int n_bins, double lower_edge, double upper_edge) {
         m_bins.clear();
         m_bins.resize(n_bins+2,0);
         m_scale = n_bins / (upper_edge - lower_edge);
         m_lowerEdge=lower_edge-binWidth();
      }
      unsigned int bin( double value) {
         return value < m_lowerEdge ? 0 : std::min(m_bins.size()-1,
                                                   static_cast<std::size_t>( (value - m_lowerEdge) * m_scale));
      }
      double lowerEdge() const { return m_lowerEdge + binWidth(); }
      double upperEdge() const { return m_lowerEdge + (m_bins.size()-1)/m_scale; }
      double binWidth()  const { return 1/m_scale; }
      const std::vector<unsigned short> &bins() const { return m_bins; }
      void fill(double value) {
         ++m_bins.at(bin(value));
      }
      double m_lowerEdge;
      double m_scale;
      std::vector<unsigned short> m_bins;
   };
   std::ostream &operator<<(std::ostream &out, const Stat &a );
}
