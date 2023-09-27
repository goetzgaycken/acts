#include "Acts/Utilities/Counter.hpp"
#include <iostream>
#include <iomanip>

namespace DbgUtils {
   thread_local Counter *g_currentCounter;

   std::mutex Counter::s_counterLabelMutex;
   std::array<std::map< unsigned int, Counter::HistoInfo >,  Counter::kNCounterTypes> Counter::s_counterBinning;
   std::array<std::map< unsigned int, Counter::CounterInfo >,Counter::kNCounterTypes> Counter::s_counterLabel;

   std::map< std::string, unsigned int> Counter::s_counterId;
   std::map< std::string, unsigned int> Counter::s_statId;

   std::map<std::string, Counter > Counter::s_algCounter;

   Counter Counter::s_dummy(true);

   Counter::Counter(bool dummy) {
      if (dummy) {
         m_idMask = 0;
         m_stat.resize(1);
         m_detailedCounter.resize(1);
         std::lock_guard<std::mutex> lock( Counter::s_counterLabelMutex);

         s_counterId.insert( std::make_pair( std::string(""), 0));
         s_statId.insert( std::make_pair( std::string(""), 0));
         Counter::s_counterLabel[kHistogram].insert( std::make_pair(0,
                                                                    CounterInfo { std::string(""),
                                                                       std::string(""),
                                                                       0u}));
         Counter::s_counterLabel[kCounter].insert( std::make_pair(0,
                                                                    CounterInfo { std::string(""),
                                                                       std::string(""),
                                                                       0u}));
         Counter::s_counterBinning[kHistogram].insert( std::make_pair(0,
                                                                      HistoInfo { 1, 0., 1.}));
         m_stat[0] .setBinning(1,0.,1.);
      }
   }

   
   void Counter::setThreadCounter(const std::string &name) {
      std::lock_guard<std::mutex> lock( Counter::s_counterLabelMutex);
      std::cout << "DEBUG Counter::setThreadCounter " << name << std::endl;

      std::map<std::string, Counter >::iterator
         counter_iter = s_algCounter.find(name);
      if (counter_iter == s_algCounter.end()) {
         std::pair<std::map<std::string, Counter >::iterator, bool>
            ret = s_algCounter.insert( std::make_pair(name, Counter() ) );
         if (!ret.second) {
            throw std::runtime_error("Failed to create counter for algorithm.");
         }
         counter_iter = ret.first;
      }
      if (   counter_iter->second.m_detailedCounter.size() + counter_iter->second.m_stat.size()
          <  s_counterId.size() + s_statId.size()) {

         if (counter_iter->second.m_detailedCounter.size() <  s_counterId.size()) {
            counter_iter->second.m_detailedCounter.resize( s_counterId.size() );
         }
         if (counter_iter->second.s_counterId.size() <  s_statId.size()) {
            unsigned int old_size =static_cast<unsigned int>(counter_iter->second.m_stat.size());
            counter_iter->second.m_stat.resize( s_statId.size() );
            counter_iter->second.updateBinning(old_size);
         }
      }
      g_currentCounter = &(counter_iter->second); //&(s_algCounter[name]);
   }


   void Counter::equalOrThrowCounterId(unsigned int counter_id, const char *name, const char *file, unsigned int line) {
      const CounterInfo &info = Counter::s_counterLabel[kCounter].at(counter_id);
      if (   info.m_file != file
             || info.m_line != line
             || info.m_title != name) {
         std::stringstream msg;
         msg << "Already registered counter " << counter_id << ". Failed to insert " << name << " ( " << file << " : " << line << " )"
             << " Registered: " << info.m_title << ", " << info.m_file << ":" << info.m_line;
         throw std::logic_error(msg.str());
      }
   }

   void Counter::equalOrThrowHistoId(unsigned int histo_id, const char *name, const char *file, unsigned int line) {
      const CounterInfo &info = Counter::s_counterLabel[kHistogram][histo_id];
      if (   info.m_file != file
             || info.m_line != line
             || info.m_title != name) {
         std::stringstream msg;
         msg << "Already registered histogram " << histo_id << ". Failed to insert " << name << " ( " << file << " : " << line << " )"
             << " Registered: " << info.m_title << ", " << info.m_file << ":" << info.m_line;
         throw std::logic_error(msg.str());
      }
   }

   void Counter::equalBinningOrThrow(unsigned int histo_id, unsigned int bins, float xmin, float xmax) {
      HistoInfo hist_info = Counter::s_counterBinning[kHistogram][histo_id];
      if (hist_info.m_nBins != bins
          || hist_info.m_xmin != xmin
          || hist_info.m_xmin != xmax) {                  
         std::stringstream msg;
         msg << "Already registered binning for " << histo_id << ". Reject " << bins << ", " << xmin << ", " << xmax
             << " Registered: " <<  hist_info.m_nBins << ", " << hist_info.m_xmin << ", " << hist_info.m_xmax;
         throw std::logic_error(msg.str());
      }
   }

   Counter::CounterId::CounterId(const char *name, const char *file, unsigned int line) {
      std::lock_guard<std::mutex> lock(Counter::s_counterLabelMutex);
      unsigned int counter_id  = getId( s_counterId, name);
      if (!Counter::s_counterLabel[kCounter].insert( std::make_pair(counter_id,
                                                                    CounterInfo { std::string(name),
                                                                       std::string(file),
                                                                       line}) ).second) {
         Counter::equalOrThrowCounterId(counter_id, name, file, line);
      }
      m_id=counter_id;
   }
   Counter::CounterId::CounterId(const char *name,
                                 const char *file,
                                 unsigned int line,
                                 unsigned int bins,
                                 float xmin,
                                 float xmax) {
      std::lock_guard<std::mutex> lock(Counter::s_counterLabelMutex);
      unsigned int counter_id  = getId( s_statId, name);
      if (!Counter::s_counterLabel[kHistogram].insert( std::make_pair(counter_id,
                                                                      CounterInfo { std::string(name),
                                                                         std::string(file),
                                                                         line}) ).second) {
         Counter::equalOrThrowHistoId(counter_id, name, file, line);
      }
      if (!Counter::s_counterBinning[kHistogram].insert( std::make_pair(counter_id,
                                                                        HistoInfo { bins, xmin, xmax}) ).second) {
         Counter::equalBinningOrThrow(counter_id, bins, xmin, xmax);
      }
      m_id=counter_id;
   }
      
   void Counter::updateBinning(unsigned int start) {
      for (unsigned int stat_i = start ; stat_i < m_stat.size(); ++stat_i) {
         std::map< unsigned int, HistoInfo >::const_iterator binning_iter =  s_counterBinning[kHistogram].find(stat_i);
         if (binning_iter != s_counterBinning[kHistogram].end() && m_stat.at(stat_i).n() == 0) {
            m_stat[stat_i] .setBinning(binning_iter->second.m_nBins,
                                       binning_iter->second.m_xmin,
                                       binning_iter->second.m_xmax);
         }
      }
   }

   const Counter *Counter::getCounter(const std::string &name) {
      std::lock_guard<std::mutex> lock( Counter::s_counterLabelMutex);
      std::map<std::string, Counter >::const_iterator
         counter_iter = s_algCounter.find(name);
      return counter_iter == s_algCounter.end() ? nullptr : &(counter_iter->second);
   }


   void Counter::resetThreadCounter() {
      {
      std::lock_guard<std::mutex> lock( Counter::s_counterLabelMutex);
      std::cout << "DEBUG Counter::resetThreadCounter " << std::endl;
      }
      g_currentCounter = nullptr;
   }
   Counter &Counter::getThreadCounter() {
      Counter *counter = g_currentCounter;
      if (!counter) {
         return s_dummy;
         //  throw std::logic_error("No counter set for this thread.");
      }
      return *counter;
   }

   Counter *Counter::statCounter(unsigned int  stat_i) {
      Counter *counter = g_currentCounter;
      if (!counter) {
         // throw std::logic_error("No counter set for this thread.");
         counter=&s_dummy;
      }
      stat_i &= counter->m_idMask;
      if  (stat_i >= counter->m_stat.size() ) {
         unsigned int old_size =static_cast<unsigned int>(counter->m_stat.size());
         counter->m_stat.resize( std::max(stat_i+1, static_cast<unsigned int>(Counter::s_statId.size())) );
         counter->updateBinning(old_size);
      }
      return counter;
   }
   

   void Counter::dump(std::ostream &msg_out) const {
     const CounterInfo *first_counter=nullptr;
     std::size_t max_equal=std::numeric_limits<std::size_t>::max();
     std::size_t max_length=0u;
     for (unsigned int type_i=0; type_i<kNCounterTypes; ++type_i) {
        for (const std::pair< const unsigned int, CounterInfo  > &counter_label : s_counterLabel[type_i] ) {
           if (!first_counter) {
              first_counter = &counter_label.second;
              max_equal= counter_label.second.m_file.size();
              max_length= counter_label.second.m_title.size();
              std::size_t pos = counter_label.second.m_file.rfind("/");
              if (pos != std::string::npos && pos+1 < counter_label.second.m_file.size()) {
                 max_equal = pos+1;
              }
              continue;
           }
           max_length=std::max(max_length, counter_label.second.m_title.size());
           std::size_t a_length = std::min(counter_label.second.m_file.size(), max_equal);
           std::size_t a_equal=0;
           for (;a_equal<a_length; ++a_equal) {
              if (counter_label.second.m_file[a_equal] != first_counter->m_file[a_equal]) break;
           }
           max_equal=std::min(max_equal, a_equal);
        }
     }
     std::string empty{};
     std::stringstream a_msg;
     a_msg << "---- Detailed counter:" << std::endl;
     for (unsigned int counter_i=0; counter_i < m_detailedCounter.size(); ++counter_i) {
        if (m_detailedCounter[counter_i]>0) {
           std::map< unsigned int, CounterInfo >::const_iterator label_iter =  s_counterLabel[kCounter].find(counter_i);
           a_msg << std::setw(12) << m_detailedCounter[counter_i]
                 << " " << std::left << std::setw(max_length)
                 << (label_iter != s_counterLabel[kCounter].end() ? label_iter->second.m_title : empty)
                 << std::right
                 << " | "
                 << (label_iter != s_counterLabel[kCounter].end()
                     ? (label_iter->second.m_file.size()  >= max_equal
                        ? label_iter->second.m_file.substr(max_equal, label_iter->second.m_file.size() - max_equal)
                        : label_iter->second.m_file.substr(0,label_iter->second.m_file.size() ) )
                     : empty.substr(0,empty.size()))
                 <<  " : "
                 << (label_iter != s_counterLabel[kCounter].end()
                     ? label_iter->second.m_line
                     : 0)
                 << std::endl;
        }
     }
     msg_out << a_msg.str(); // << std::endl;
     a_msg.str("");
     a_msg << "---- Detailed stat:" << std::endl;
     for (unsigned int stat_i=0; stat_i < m_stat.size(); ++stat_i) {
        if (m_stat[stat_i].n()>0) {
           std::map< unsigned int, CounterInfo >::const_iterator label_iter =  s_counterLabel[kHistogram].find(stat_i);
           a_msg << m_stat[stat_i]
                 << " " << std::left << std::setw(max_length)
                 << (label_iter != s_counterLabel[kHistogram].end() ? label_iter->second.m_title : empty)
                 << std::right
                 << " | "
                 << (label_iter != s_counterLabel[kHistogram].end()
                     ? (label_iter->second.m_file.size()  >= max_equal
                        ? label_iter->second.m_file.substr(max_equal, label_iter->second.m_file.size() - max_equal)
                        : label_iter->second.m_file.substr(0,label_iter->second.m_file.size() ) )
                     : empty.substr(0,empty.size()))
                 <<  " : "
                 << (label_iter != s_counterLabel[kHistogram].end()
                     ? label_iter->second.m_line
                     : 0)
                 << std::endl
                 << m_stat[stat_i].histogramToString()
                 << std::endl;
        }
     }
     msg_out << a_msg.str();
   }
   void Counter::resetCounter() {
     for (unsigned int counter_i=0; counter_i < m_detailedCounter.size(); ++counter_i) {
        m_detailedCounter[counter_i]=0;
     }
     for (unsigned int stat_i=0; stat_i < m_stat.size(); ++stat_i) {
        m_stat[stat_i].reset();
     }
   }

   std::string ratioWithError( const Stat &a, const Stat &b) {
      std::stringstream msg;
      double ratio = (a.mean() / b.mean());
      msg << std::setw(14) << ratio;
      if (a.n() > 1 || b.n()>1) {
         // sqrt( 
         //       (1/b * sigma_a)^2
         //     + (a/b^2 * signa_b)^2 )

         double uncertainty = std::sqrt( a.rms2() / sqr( b.mean() ) + b.rms2() * sqr( a.mean()/sqr(b.mean()) ) )  / ratio * 100;
         msg << "(+-" << std::setprecision(1) << std::fixed << std::setw(5) <<  uncertainty << "%)";
      }
      return msg.str();
   }
}

namespace {
   void dummy()  {
    DEBUG_HISTOGRAM_COUNTER("DummyHisto" ,__FILE__,__LINE__, 1, 25,-.5, 25-.5);
    DEBUG_INCREMENT_COUNTER("DummyCounter" ,__FILE__,__LINE__);
   }
}
