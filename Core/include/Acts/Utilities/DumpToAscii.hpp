#include <ostream>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

namespace Acts {
   void dumpParameters(std::ostream &out,
                       const Acts::BoundVector &curvi_param,
                       const Acts::BoundSquareMatrix &curvi_cov);
   void dumpPosition(std::ostream &out,
                     const Acts::Vector3 &position);
   void dumpMomentum(std::ostream &out,
                     const Acts::Vector3 &momentum);

   void dumpMagField(std::ostream &out,
                     const Acts::Vector3 &mag_field);
      
   void dumpSurface(std::ostream &out,
                    const Acts::GeometryContext &tgContext,
                    const Acts::Surface &surface,
                    const Acts::Vector3 &position,
                    const Acts::Vector3 &direction);

   Acts::Vector3 boundCenter(const Acts::GeometryContext &tgContext,
                             const Acts::Surface &surface,
                             const Acts::Vector3 &position,
                             const Acts::Vector3 &direction);

}
