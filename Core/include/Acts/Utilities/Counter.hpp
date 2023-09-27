#pragma once

#include <cmath>
#include <iomanip>
#include <ostream>
#include <string>
#include <sstream>
#include <vector>
#include <mutex>
#include <array>
#include <map>
#include <cassert>

#include <stdexcept>
#include <iostream>

namespace DbgUtils {

   //
   // USAGE:
   // StatusCode TrackFindingAlg::execute(const EventContext &ctx) const
   // {
   //    DbgUtils::CounterGuard set_counter(name());
   //    ...
   //    ...count...
   //    ...
   //    if (msgLvl(MSG::INFO)) {
   //        const DbgUtils::Counter &counter=DbgUtils::Counter::getThreadCounter();
   //        std::stringstream msg_out;
   //        counter.dump(msg_out);
   //        msg(MSG::INFO) << msg_out.str() << endmsg;
   // }

   
class Stat {
public:
   void add(double val) {
      ++m_n;
      m_sum += val;
      m_sum2 += val*val;
      m_min=std::min(m_min,val);
      m_max=std::max(m_max,val);
   }
   unsigned int n() const { return m_n; }
   double min() const { return m_min; }
   double max() const { return m_max; }
   double mean() const { return m_n>0 ? m_sum/m_n : 0.; }
   double rms2() const { return m_n>1 ? (m_sum2 - m_sum *m_sum/m_n)/(m_n-1) : 0.; }
   double rms() const { return std::sqrt( rms2() ); }
   Stat &operator -=(const Stat &b) {
      if (m_n != b.m_n) {
         double scale = 1. * m_n / b.m_n;
         double corr = 2*m_sum * b.m_sum * scale/m_n;
         m_sum  -= b.m_sum * scale;
         m_sum2 += b.m_sum2 * scale;
         m_sum2 -= corr;
      }
      else {
         double corr = 2*m_sum * b.m_sum / m_n;
         m_sum -= b.m_sum;
         m_sum2 += b.m_sum2;
         m_sum2 -= corr;
      }
      m_min = m_min - b.m_max;
      m_max = m_max - b.m_min;
      return *this;
   }
   void reset() {
      m_n=0;
      m_sum=0.;
      m_sum2=0.;
      m_min=std::numeric_limits<double>::max();
      m_max=-std::numeric_limits<double>::max();
   }

   unsigned int m_n=0;
   double m_sum=0.;
   double m_sum2=0.;
   double m_min=std::numeric_limits<double>::max();
   double m_max=-std::numeric_limits<double>::max();
};

inline std::ostream &operator<<(std::ostream &out, const Stat &stat) {
   if (stat.n() > 1) {
   out << std::setw(14) << stat.min() << " < "
       << std::setw(14) << stat.mean() << " +- " << std::setw(14) <<  stat.rms() << " < "
       << std::setw(14) << stat.max()
       << " / " << std::setw(9) << stat.n();
   }
   else {
      out << std::setw(14*4+9+3*3+4) << stat.mean();
   }
   return out;
}

class StatHist : public Stat {
public:
   StatHist() : StatHist(20,0,40) {}

   StatHist(unsigned int n_bins, float xmin, float xmax) {
      setBinning(n_bins,xmin,xmax);
   }

   void setBinning(unsigned int n_bins, float xmin, float xmax)
   {
      m_xmin=xmin;
      m_scale = ( n_bins / (xmax-xmin) );
      m_xmin -= 1./m_scale;
      m_histogram.resize(n_bins+2,0u);
   }

   void add(double val) {
      Stat::add(val);
      if (!m_histogram.empty()) {
         unsigned int bin = std::min( static_cast<unsigned int>(m_histogram.size()-1), static_cast<unsigned int>( std::max(0.,(val - m_xmin)*m_scale)) );
         ++m_histogram.at(bin) ;
      }
   }
   void reset() {
      Stat::reset();
      for (unsigned int &bin : m_histogram) {
         bin = 0u;
      }
   }

   std::string histogramToString() const {
      std::stringstream msg;
      if (m_histogram.size()>2) {
         unsigned int max_val = 0;
         for (const auto &count : m_histogram) {
            max_val = std::max(max_val, count);
         }
         double bin_width=1./m_scale;
         unsigned int w = static_cast<unsigned int>(log(1.*max_val) / log(10.))+1;
         msg << (m_xmin+bin_width) << " .. " << ((m_histogram.size()-2)/m_scale + m_xmin+bin_width) << " : "
            << std::setw(w) << m_histogram[0] << " |";
         for (unsigned int i=1; i<m_histogram.size()-1; ++i) {
            msg << " " << std::setw(w) << m_histogram[i];
         }
         msg << " | " << std::setw(w) << m_histogram.back();
      }
      return msg.str();
   }

   double m_xmin;
   double m_scale;
   std::vector<unsigned int> m_histogram;
};

   inline double sqr(double a) { return a*a; }
   std::string ratioWithError( const Stat &a, const Stat &b);

  class Counter {
  public:
     class LabelSetter;
     friend class LabelSetter;

     std::vector< unsigned int> m_detailedCounter {};

     std::vector< StatHist> m_stat {};

     static std::mutex s_counterLabelMutex ;

     struct CounterInfo {
        std::string m_title;
        std::string m_file;
        unsigned int m_line = 0;
     };

     struct HistoInfo {
        unsigned int m_nBins=25;
        float        m_xmin =-.5;
        float        m_xmax = 25-.5;
     };

     static std::map< std::string, unsigned int> s_counterId;
     static std::map< std::string, unsigned int> s_statId;
     enum ECounterType {
        kCounter,
        kHistogram,
        kNCounterTypes
     };
     static std::array<std::map< unsigned int, HistoInfo >,  kNCounterTypes> s_counterBinning;
     static std::array<std::map< unsigned int, CounterInfo >,kNCounterTypes> s_counterLabel;
     unsigned int m_idMask = std::numeric_limits<unsigned int>::max();

  protected:
     static void equalOrThrowCounterId(unsigned int counter_id, const char *name, const char *file, unsigned int line);
     static void equalOrThrowHistoId(unsigned int histo_id, const char *name, const char *file, unsigned int line);
     static void equalBinningOrThrow(unsigned int histo_id, unsigned int bins, float xmin, float xmax);
  public:     
     class CounterId {
     private:
        unsigned int m_id = std::numeric_limits<unsigned int>::max();
     public:
        CounterId(const char *name, const char *file, unsigned int line);
        CounterId(  const char *name,
                    const char *file,
                    unsigned int line,
                    unsigned int bins,
                    float xmin,
                    float xmax);

        unsigned int id() const { return m_id; }
     private:
        static unsigned int getId(std::map< std::string, unsigned int> &id_map, const char *name) {
           std::string tmp(name);
           std::map< std::string, unsigned int>::iterator iter = id_map.find( tmp );
           if (iter == id_map.end()) {
              std::pair<std::map< std::string, unsigned int>::iterator, bool >
                 ret = id_map.insert( std::make_pair( std::move(tmp), id_map.size()));
              iter = ret.first;
           }
           return iter->second;
        }


     };

     Counter() {}
     Counter(bool dummy);
     void dump(std::ostream &out) const;
     void resetCounter();
     void updateBinning(unsigned int start=0);

     static std::map<std::string, Counter > s_algCounter ;
     static void setThreadCounter(const std::string &name);
     static void resetThreadCounter();
     static Counter &getThreadCounter();
     static const Counter *getCounter(const std::string &name);
     static unsigned int &getThreadCounter(unsigned int counter_id) {
        Counter &counter = getThreadCounter();
        // std::cout << "DEBUG getThreadCounter " << counter_id  << " counter "
        //           << static_cast< const void *>(&counter)
        //           << std::endl;
        // assert( counter_id < counter.m_detailedCounter.size() );
        counter_id  &= counter.m_idMask;
        if (counter_id >= counter.m_detailedCounter.size()) {
           counter.m_detailedCounter.resize( std::max(counter_id+1, static_cast<unsigned int>(Counter::s_counterId.size())), 0u);
        }
        // if ( counter_id < counter.m_detailedCounter.size() ) {
        //    throw std::range_error("Counter index out of bounrds.");
        // }
        return counter.m_detailedCounter[counter_id ];
     }
     template <class T>
     static void gatherStat(unsigned int  stat_i, T value) {
        Counter *counter=statCounter(stat_i);
        counter->m_stat[stat_i & counter->m_idMask].add(static_cast<double>(value));
     }


  private:
     static Counter *statCounter(unsigned int  stat_i);
     static Counter s_dummy;
  };

   class CounterGuard {
   public:
      CounterGuard(const std::string &a_name) { DbgUtils::Counter::setThreadCounter(a_name ) ; }
      ~CounterGuard() { DbgUtils::Counter::resetThreadCounter( ) ; }
   };

}
#define ENABLE_DEBUG_COUNTER
#ifdef ENABLE_DEBUG_COUNTER
#define DEBUG_ADD_TO_COUNTER(name, file, line, increment)   \
   { static const DbgUtils::Counter::CounterId counter_id(name, file, line );       \
     DbgUtils::Counter::getThreadCounter(counter_id.id()) += increment; } \
   do {} while (0)

#define DEBUG_INCREMENT_COUNTER(name, file, line)  DEBUG_ADD_TO_COUNTER(name, file, line, 1)

#define DEBUG_CREATE_ID(name, file, line, var_name, nbins, xmin, xmax) \
   const DbgUtils::Counter::CounterId var_name(name,file, line, nbins, xmin,xmax);

#define DEBUG_HISTOGRAM_COUNTER(name, file, line, increment, nbins, xmin, xmax) \
   { static const DbgUtils::Counter::CounterId counter_id(name,file, line, nbins, xmin,xmax); \
     DbgUtils::Counter::gatherStat(counter_id.id(), increment); }                \
   do {} while (0)

#define DEBUG_HISTOGRAM_COUNTER_USING_ID(counter_id, increment) \
     DbgUtils::Counter::gatherStat(counter_id.id(), increment); \
   do {} while (0)
#else
#define DEBUG_ADD_TO_COUNTER(name, file, line, increment)   \
   do {} while (0)

#define DEBUG_INCREMENT_COUNTER(name, file, line)  DEBUG_ADD_TO_COUNTER(name, file, line, 1)

#define DEBUG_CREATE_ID(name, file, line, var_name)     \
   constexpr int var_name = 0

#define DEBUG_HISTOGRAM_COUNTER(name, file, line, increment, nbins, xmin, xmax) \
   do {} while (0)
#define DEBUG_HISTOGRAM_COUNTER_USING_ID(counter_id, increment) \
   do {} while (0)
#endif
