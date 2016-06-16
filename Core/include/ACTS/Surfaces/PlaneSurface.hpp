// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlaneSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_PLANESURFACE_H
#define ACTS_SURFACES_PLANESURFACE_H 1

#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Identifier.hpp"

namespace Acts {

class DetectorElementBase;

///
/// @class PlaneSurface
///
/// Class for a planaer in the TrackingGeometry.
///
/// The PlaneSurface extends the Surface class with the possibility to
/// convert local to global positions (vice versa).
///
/// @image html PlaneSurface.gif
///
class PlaneSurface : public Surface
{
public:
  /// Default Constructor - needed for persistency
  PlaneSurface();

  /// Copy Constructor
  PlaneSurface(const PlaneSurface& psf);

  /// Copy Constructor with shift
  PlaneSurface(const PlaneSurface& psf, const Transform3D& transf);

  /// Dedicated Constructor with normal vector 
  PlaneSurface(const Vector3D& position, const Vector3D& normal);

  /// Constructor from DetectorElementBase - potentially with identifier 
  PlaneSurface(const DetectorElementBase& detelement,
               const Identifier&          identifier = Identifier());

  /// Constructor for planar Surface without Bounds
  /// @param htrans transform in 3D that positions this surface              
  PlaneSurface(std::shared_ptr<Transform3D> htrans);

  /// Constructor for Planes with shared bounds object 
  /// @param htrans transform in 3D that positions this surface              
  /// @param pbounds bounds object to describe the actual surface area 
  /// @attention the pointer to pbounds must not be a nullptr 
  PlaneSurface(std::shared_ptr<Transform3D>        htrans,
               std::shared_ptr<const PlanarBounds> pbounds);

  /// Destructor
  virtual ~PlaneSurface();

  /// Assignment operator
  /// @param psf source PlaneSurface for assignment
  PlaneSurface&
  operator=(const PlaneSurface& psf);

  /// Comparison: equality operator
  /// @param sf source Surface for comparison
  virtual bool
  operator==(const Surface& sf) const override;

  /// Virtual constructor with optional shift 
  /// ownership of the shift transform is not given !!
  /// @copydoc Surface::clone
  virtual PlaneSurface*
  clone(const Transform3D* shift = nullptr) const override;

  /// Normal vector return
  /// @param lpos is the local position is ignored
  /// return a Vector3D by value
  const Vector3D normal(const Vector2D& lpos = Vector2D()) const;

  /// Return the surface type 
  virtual SurfaceType
  type() const override
  {
    return Surface::Plane;
  }

  /// Return method for bounds object of this surfrace
  virtual const PlanarBounds&
  bounds() const override;

  /// Geometrical on surface test
  /// This method returns true if the GlobalPosition is on the Surface for both,
  /// within or without check of whether the local position is inside boundaries or not
  /// @param gpos global position to be checked
  /// @param bchk gboundary check directive
  virtual bool
  isOnSurface(const Vector3D&      gpos,
              const BoundaryCheck& bchk = true) const override;

  /// @copydoc Surface::localToGlobal
  /// For planar surfaces the momentum is ignroed in the local to global transformation              
  virtual void
  localToGlobal(const Vector2D& lpos,
                const Vector3D& mom,
                Vector3D&       gpos) const override;

  /// @copydoc Surface::globalToLocal
  /// For planar surfaces the momentum is ignroed in the gloabl to l transformation              
  virtual bool
  globalToLocal(const Vector3D& gpos,
                const Vector3D& mom,
                Vector2D&       lpos) const override;

  ///  fast straight line intersection schema - standard: provides closest
  /// intersection and (signed) path length
  ///  forceDir is to provide the closest forward solution
  /// 
  ///  <b>mathematical motivation:</b>
  /// 
  ///  the equation of the plane is given by: <br>
  ///  @f$ \vec n \cdot \vec x = \vec n \cdot \vec p,@f$ <br>
  ///  where @f$ \vec n = (n_{x}, n_{y}, n_{z})@f$ denotes the normal vector of
  /// the plane,
  ///  @f$ \vec p = (p_{x}, p_{y}, p_{z})@f$ one specific point on the plane and
  /// @f$ \vec x = (x,y,z) @f$ all possible points
  ///  on the plane.<br>
  ///  Given a line with:<br>
  ///  @f$ \vec l(u) = \vec l_{1} + u \cdot \vec v @f$, <br>
  ///  the solution for @f$ u @f$ can be written:
  ///  @f$ u = \frac{\vec n (\vec p - \vec l_{1})}{\vec n \vec v}@f$ <br>
  ///  If the denominator is 0 then the line lies:
  ///  - either in the plane
  ///  - perpenticular to the normal of the plane
  /// 
  /// 
  virtual Intersection
  intersectionEstimate(const Vector3D&      pos,
                       const Vector3D&      dir,
                       bool                 forceDir,
                       const BoundaryCheck& bchk = true) const override;

  /// Return properly formatted class name for screen output 
  virtual std::string
  name() const override
  {
    return "Acts::PlaneSurface";
  }

protected:                                                     
  /// PlanarBounds - this can be nullptr if the Surface is a PROXY
  std::shared_ptr<const PlanarBounds>       m_bounds;
  Vector3D                                  m_normal;
};

inline PlaneSurface*
PlaneSurface::clone(const Transform3D* shift) const
{
  if (shift) new PlaneSurface(*this, *shift);
  return new PlaneSurface(*this);
}

inline const PlanarBounds&
PlaneSurface::bounds() const
{
  if (m_bounds) return (*m_bounds.get());
  if (Surface::m_associatedDetElement
      && Surface::m_associatedDetElementId.is_valid()) {
    return m_associatedDetElement->bounds(Surface::m_associatedDetElementId);
  }
  return m_associatedDetElement->bounds();
}

inline const Vector3D 
PlaneSurface::normal(const Vector2D& lpos = Vector2D()) const
{
    return m_normal;
}


inline Intersection
PlaneSurface::intersectionEstimate(const Vector3D&      pos,
                                   const Vector3D&      dir,
                                   bool                 forceDir,
                                   const BoundaryCheck& bchk) const
{
  double denom = dir.dot(normal());
  if (denom) {
    double   u = (normal().dot((center() - pos))) / (denom);
    Vector3D intersectPoint(pos + u * dir);
    // evaluate the intersection in terms of direction
    bool isValid = forceDir ? (u > 0.) : true;
    // evaluate (if necessary in terms of boundaries)
    isValid = bchk ? (isValid && isOnSurface(intersectPoint, bchk)) : isValid;
    // return the result
    return Intersection(intersectPoint, u, isValid);
  }
  return Intersection(pos, 0., false);
}

}  // end of namespace

#endif  // ACTS_SURFACES_PLANESURFACE_H
