#include "ActsExamples/Utilities/Stat.h"
#include <iomanip>

namespace ActsExamples {
   std::ostream &operator<<(std::ostream &out, const Stat &a ) {
      out <<           std::setw(14) << a.min() << " < "
          <<           std::setw(14) << a.mean()
          << " +- " << std::setw(14) << a.rms()
          << " < "  << std::setw(14) << a.max()
          << " / "  << std::setw(9)  << a.n();
      if (!a.bins().empty()) {
         out << std::endl;
         out << a.lowerEdge() << ", " << a.lowerEdge()+a.binWidth() << ".." << a.upperEdge() << ": ";
         out << std::setw(5) << a.bins()[0] << " |";
         for (unsigned int bin_i=1; bin_i<a.bins().size()-1; ++bin_i) {
            out << " " << std::setw(5) << a.bins()[bin_i];
         }
         out << " | " << std::setw(5) << a.bins().back() ;
      }
      return out;
   }
}
