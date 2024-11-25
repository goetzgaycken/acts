#include "Acts/Utilities/DumpToAscii.hpp"

#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

#include <sstream>
#include <stdexcept>

namespace Acts {
   void dumpParameters(std::ostream &out,
                       const Acts::BoundVector &curvi_param,
                       const Acts::BoundSquareMatrix &curvi_cov) {
      out << " parm";
      for(const auto &param : curvi_param ) {
         out << " " << param;
      }
      //      std::array<double, 6> scale { 1., 1., 1., 1., 1., 1. };
      out << " paramCov";
      for (unsigned int row_i=0; row_i < curvi_cov.rows(); ++row_i) {
         for (unsigned int col_i=0; col_i < curvi_cov.cols(); ++col_i) {
            out << " " << curvi_cov(row_i,col_i);
            //* scale.at(row_i)*scale.at(col_i);
         }
      }
   }
   void dumpPosition(std::ostream &out,
                     const Acts::Vector3 &position) {
      out << " pos " << position[0] << " " << position[1] << " " << position[2];
   }
   void dumpMomentum(std::ostream &out,
                     const Acts::Vector3 &momentum) {
      out << " mom " << momentum[0] << " " << momentum[1] << " " << momentum[2];      
   }

   void dumpMagField(std::ostream &out,
                     const Acts::Vector3 &mag_field)
   {
      out << " mag " << mag_field[0] << " " << mag_field[1] << " " << mag_field[2];
   }


   Acts::Vector3 boundCenter(const Acts::GeometryContext &tgContext,
                             const Acts::Surface &surface,
                             const Acts::Vector3 &position,
                             const Acts::Vector3 &direction)  {
      const auto *bounds = &(surface.bounds());
      const Acts::PlanarBounds *planar_bounds =  dynamic_cast<const Acts::PlanarBounds *>(bounds);
      std::vector<Acts::Vector2> vertices;
      if (planar_bounds) {
         vertices = planar_bounds->vertices(12);
      }
      else {
         const Acts::DiscBounds *disc_bounds =  dynamic_cast<const Acts::DiscBounds *>(bounds);
         if (disc_bounds) {
            vertices = disc_bounds->vertices(120);
         }
      }
      if (!vertices.empty()) {
         Acts::Vector2 center = Vector2::Zero();
         for (const Acts::Vector2 &vertex : vertices ) {
            center += vertex;
         }
         center /= (1. * vertices.size());

         Acts::Vector3 a_vertex_global( surface.center(tgContext)
                                  +   surface.referenceFrame(tgContext, position, direction)
                                    * Acts::Vector3{center(0,0), center(1,0), 0.} );
         return a_vertex_global;
      }
      else {
         return Acts::Vector3::Zero();
      }
   }
   
      
   void dumpSurface(std::ostream &out,
                    const Acts::GeometryContext &tgContext,
                    const Acts::Surface &surface,
                    const Acts::Vector3 &position,
                    const Acts::Vector3 &direction) {
      const Acts::RotationMatrix3 rframe = surface.referenceFrame(tgContext, position, direction);
      out << " surfaceRot";
      for (unsigned int row_i=0; row_i < rframe.rows(); ++row_i) {
         for (unsigned int col_i=0; col_i < rframe.cols(); ++col_i) {
            out << " " << rframe(row_i,col_i);
         }
      }
      out << " surfaceCenter";
      const auto &center = surface.center(tgContext);
      for (unsigned int row_i=0; row_i < center.rows(); ++row_i) {
         for (unsigned int col_i=0; col_i < center.cols(); ++col_i) {
            out << " " << center(row_i,col_i);
         }
      }
      out << " " << typeid(surface).name();
      double r = 0.;
      const Acts::CylinderSurface *cylinder = dynamic_cast<const Acts::CylinderSurface *>(&surface);
      if (cylinder) {
         r = cylinder->bounds().get(Acts::CylinderBounds::eR);
      }
      out << " "  << r;
      const auto *bounds = &(surface.bounds());
      const Acts::PlanarBounds *planar_bounds =  dynamic_cast<const Acts::PlanarBounds *>(bounds);
      std::vector<Acts::Vector2> vertices;
      if (planar_bounds) {
         vertices = planar_bounds->vertices(12);
      }
      else {
         const Acts::DiscBounds *disc_bounds =  dynamic_cast<const Acts::DiscBounds *>(bounds);
         if (disc_bounds) {
            vertices = disc_bounds->vertices(120);
         }
      }
      if (!vertices.empty()) {
         out << " surfaceBounds " << vertices.size();
         for (const Acts::Vector2 &vertex : vertices ) {
            out << " " << vertex(0,0)  << " " << vertex(1,0);
         }
      }
      out << " geoID " << surface.geometryId().value();
   }

}
