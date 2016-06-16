// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/ConeSurface.hpp"
#include "ACTS/Surfaces/RealQuadraticEquation.hpp"
#include <assert.h>
#include <iomanip>
#include <iostream>

Acts::ConeSurface::ConeSurface(const ConeSurface& csf)
  : Acts::Surface(csf), m_bounds(csf.m_bounds)
{
}

Acts::ConeSurface::ConeSurface(const ConeSurface&       csf,
                               const Acts::Transform3D& transf)
  : Acts::Surface(csf, transf)
  , m_bounds(csf.m_bounds)
{
}

Acts::ConeSurface::ConeSurface(std::shared_ptr<Acts::Transform3D> htrans,
                               double                             alpha,
                               bool                               symmetric)
  : Acts::Surface(htrans)
  , m_bounds(std::make_shared<Acts::ConeBounds>(alpha, symmetric))
{
}

Acts::ConeSurface::ConeSurface(std::shared_ptr<Acts::Transform3D> htrans,
                               double                             alpha,
                               double                             zmin,
                               double                             zmax,
                               double                             halfPhi)
  : Acts::Surface(htrans)
  , m_bounds(std::make_shared<Acts::ConeBounds>(alpha, zmin, zmax, halfPhi))
{
}

Acts::ConeSurface::ConeSurface(std::shared_ptr<Acts::Transform3D>      htrans,
                               std::shared_ptr<const Acts::ConeBounds> cbounds)
  : Acts::Surface(htrans), m_bounds(cbounds)
{
  assert(cbounds);
}

Acts::ConeSurface::~ConeSurface()
{
}

Acts::Vector3D
Acts::ConeSurface::binningPosition(Acts::BinningValue bValue) const
{
  // special binning type for R-type methods
  if (bValue == Acts::binR || bValue == Acts::binRPhi)
    return Acts::Vector3D(
        center().x() + bounds().r(), center().y(), center().z());
  // give the center as default for all of these binning types
  // binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta
  return Acts::Surface::binningPosition(bValue);
}

Acts::ConeSurface&
Acts::ConeSurface::operator=(const ConeSurface& csf)
{
  if (this != &csf) {
    Acts::Surface::operator=(csf);
    m_bounds               = csf.m_bounds;
  }
  return *this;
}

const Acts::Vector3D
Acts::ConeSurface::rotSymmetryAxis() const
{
  return std::move(transform().rotation().col(2));
}

const Acts::RotationMatrix3D
Acts::ConeSurface::measurementFrame(const Acts::Vector3D& pos,
                                    const Acts::Vector3D&) const
{
  Acts::RotationMatrix3D mFrame;
  // construct the measurement frame
  Acts::Vector3D measY(
      transform().rotation().col(2));  // measured Y is the z axis
  Acts::Vector3D measDepth
      = Acts::Vector3D(pos.x(), pos.y(), 0.)
            .unit();  // measured z is the position transverse normalized
  Acts::Vector3D measX(
      measY.cross(measDepth).unit());  // measured X is what comoes out of it
  // the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  //!< @TODO fold in alpha
  // return it
  return std::move(mFrame);
}

void
Acts::ConeSurface::localToGlobal(const Acts::Vector2D& lpos,
                                 const Acts::Vector3D&,
                                 Acts::Vector3D& gpos) const
{
  // create the position in the local 3d frame
  double         r   = lpos[Acts::eLOC_Z] * bounds().tanAlpha();
  double         phi = lpos[Acts::eLOC_RPHI] / r;
  Acts::Vector3D loc3Dframe(r * cos(phi), r * sin(phi), lpos[Acts::eLOC_Z]);
  // transport it to the globalframe
  gpos = transform() * loc3Dframe;
}

bool
Acts::ConeSurface::globalToLocal(const Acts::Vector3D& gpos,
                                 const Acts::Vector3D&,
                                 Acts::Vector2D& lpos) const
{
  const Acts::Transform3D& surfaceTrans = transform();
  Acts::Transform3D        inverseTrans(surfaceTrans.inverse());
  Acts::Vector3D           loc3Dframe(inverseTrans * gpos);
  double                   r = loc3Dframe.z() * bounds().tanAlpha();
  lpos = Acts::Vector2D(r * atan2(loc3Dframe.y(), loc3Dframe.x()),
                          loc3Dframe.z());
  // now decide on the quility of the transformation
  double inttol = r * 0.0001;
  inttol        = (inttol < 0.01) ? 0.01 : 0.01;  // ?
  return ((fabs(loc3Dframe.perp() - r) > inttol) ? false : true);
}

Acts::Intersection
Acts::ConeSurface::intersectionEstimate(const Acts::Vector3D& pos,
                                        const Acts::Vector3D& dir,
                                        bool                  forceDir,
                                        const BoundaryCheck&  bchk) const
{
  // transform to a frame with the cone along z, with the tip at 0
  Acts::Vector3D tpos1 = transform().inverse() * pos;
  Acts::Vector3D tdir  = transform().inverse().linear() * dir;
  // see the header for the formula derivation
  double tan2Alpha = bounds().tanAlpha() * bounds().tanAlpha(),
         A         = tdir.x() * tdir.x() + tdir.y() * tdir.y()
      - tan2Alpha * tdir.z() * tdir.z(),
         B = 2 * (tdir.x() * tpos1.x() + tdir.y() * tpos1.y()
                  - tan2Alpha * dir.z() * tpos1.z()),
         C = tpos1.x() * tpos1.x() + tpos1.y() * tpos1.y()
      - tan2Alpha * tpos1.z() * tpos1.z();
  if (A == 0.) A += 1e-16;  // avoid div by zero

  // use Andreas' quad solver, much more stable than what I wrote
  Acts::RealQuadraticEquation solns(A, B, C);

  Acts::Vector3D solution(0., 0., 0.);
  double         path    = 0.;
  bool           isValid = false;
  if (solns.solutions != Acts::none) {
    double         t1 = solns.first;
    Acts::Vector3D soln1Loc(tpos1 + t1 * dir);
    isValid = forceDir ? (t1 > 0.) : true;
    // there's only one solution
    if (solns.solutions == Acts::one) {
      solution = soln1Loc;
      path     = t1;
    } else {
      double         t2 = solns.second;
      Acts::Vector3D soln2Loc(tpos1 + t2 * dir);
      // both solutions have the same sign
      if (t1 * t2 > 0. || !forceDir) {
        if (t1 * t1 < t2 * t2) {
          solution = soln1Loc;
          path     = t1;
        } else {
          solution = soln2Loc;
          path     = t2;
        }
      } else {
        if (t1 > 0.) {
          solution = soln1Loc;
          path     = t1;
        } else {
          solution = soln2Loc;
          path     = t2;
        }
      }
    }
  }
  solution = transform() * solution;

  isValid = bchk ? (isValid && isOnSurface(solution, bchk)) : isValid;
  return Acts::Intersection(solution, path, isValid);
}

double
Acts::ConeSurface::pathCorrection(const Acts::Vector3D& pos,
                                  const Acts::Vector3D& mom) const
{
  // (cos phi cos alpha, sin phi cos alpha, sgn z sin alpha)
  bool applyTransform = !(transform().isApprox(Acts::Transform3D::Identity()));
  Acts::Vector3D posLocal = applyTransform ? transform().inverse() * pos : pos;
  double         phi      = posLocal.phi();
  double         sgn      = posLocal.z() > 0. ? -1. : +1.;
  Acts::Vector3D normalC(cos(phi) * bounds().cosAlpha(),
                         sin(phi) * bounds().cosAlpha(),
                         sgn * bounds().sinAlpha());
  if (applyTransform) normalC = transform() * normalC;
  // back in global frame
  double cAlpha = normalC.dot(mom.unit());
  return fabs(1. / cAlpha);
}
