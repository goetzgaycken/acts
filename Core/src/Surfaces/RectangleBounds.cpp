// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/RectangleBounds.hpp"

#include "Acts/Surfaces/detail/BoundaryCheckHelper.hpp"

#include <iomanip>
#include <iostream>

#include "Acts/Utilities/Counter.hpp"

bool Acts::RectangleBounds::inside(
    const Acts::Vector2& lposition,
    const Acts::BoundaryTolerance& boundaryTolerance) const {
  if (boundaryTolerance.hasAbsoluteCartesian()) {
     
     DEBUG_INCREMENT_COUNTER("RectangularBound",__FILE__,__LINE__);
     const BoundaryTolerance::AbsoluteCartesian& tolerance = boundaryTolerance.asAbsoluteCartesian();
     bool ret=    lposition(0,0) > (m_min.x()-tolerance.tolerance0)
            && lposition(0,0) < (m_max.x()+tolerance.tolerance0)
            && lposition(1,0) > (m_min.y()-tolerance.tolerance1)
            && lposition(1,0) < (m_max.y()+tolerance.tolerance1);
     bool ret2= (   lposition(0,0) > (m_min.x())
                    && lposition(0,0) < (m_max.x())
                    && lposition(1,0) > (m_min.y())
                    && lposition(1,0) < (m_max.y()) );
     if (ret != ret2) {
        if (ret) {
           if (lposition(0,0) < m_min.x() /*&& lposition(0,0) > m_min.x()-tolerance.tolerance0*/) {
              DEBUG_HISTOGRAM_COUNTER("RectangularBound_posMinX",__FILE__,__LINE__, lposition(0,0) - m_min.x() , 20,-5.,0.);
           }
           if (lposition(0,0) > m_max.x() /*&& lposition(0,0) < m_min.x()-tolerance.tolerance0*/) {
              DEBUG_HISTOGRAM_COUNTER("RectangularBound_posMaxX",__FILE__,__LINE__, lposition(0,0) - m_max.x() , 20,0.,5.);
           }
           if (lposition(1,0) < m_min.y() /*&& lposition(0,0) > m_min.x()-tolerance.tolerance0*/) {
              DEBUG_HISTOGRAM_COUNTER("RectangularBound_posMinY",__FILE__,__LINE__, lposition(1,0) - m_min.y() , 20,-5.,0.);
           }
           if (lposition(1,0) > m_max.y() /*&& lposition(0,0) < m_min.x()-tolerance.tolerance0*/) {
              DEBUG_HISTOGRAM_COUNTER("RectangularBound_posMaxY",__FILE__,__LINE__, lposition(1,0) - m_max.y() , 20,0.,5.);
           }
           DEBUG_INCREMENT_COUNTER("RectangularBound_posTol",__FILE__,__LINE__);
           //           ++g_rectangleBorderCount;
        }
        if (ret2) {
           if (lposition(0,0) < m_min.x()-tolerance.tolerance0 /*&& lposition(0,0) > m_min.x()-tolerance.tolerance0*/) {
              DEBUG_HISTOGRAM_COUNTER("RectangularBound_negMinX",__FILE__,__LINE__, lposition(0,0) - m_min.x(), 20,0.,5.);
           }
           if (lposition(0,0) > m_max.x()+tolerance.tolerance0 /*&& lposition(0,0) < m_min.x()-tolerance.tolerance0*/) {
              DEBUG_HISTOGRAM_COUNTER("RectangularBound_negMaxX",__FILE__,__LINE__, lposition(0,0) - m_max.x() , 20,-5.,0.);
           }
           if (lposition(1,0) < m_min.y()-tolerance.tolerance1 /*&& lposition(0,0) > m_min.x()-tolerance.tolerance0*/) {
              DEBUG_HISTOGRAM_COUNTER("RectangularBound_negMinY",__FILE__,__LINE__, lposition(1,0) - m_min.y() , 20,0,5.);
           }
           if (lposition(1,0) > m_max.y()+tolerance.tolerance1 /*&& lposition(0,0) < m_min.x()-tolerance.tolerance0*/) {
              DEBUG_HISTOGRAM_COUNTER("RectangularBound_negMaxY",__FILE__,__LINE__, lposition(1,0) - m_max.y() , 20,5.,0.);
           }
           DEBUG_INCREMENT_COUNTER("RectangularBound_negTol",__FILE__,__LINE__);
           //++g_invRectangleBorderCount;
        }
        
     }
     return ret ;
  }
   
  return detail::insideAlignedBox(m_min, m_max, boundaryTolerance, lposition,
                                  std::nullopt);
}

std::vector<Acts::Vector2> Acts::RectangleBounds::vertices(
    unsigned int /*lseg*/) const {
  // counter-clockwise starting from bottom-left corner
  return {m_min, {m_max.x(), m_min.y()}, m_max, {m_min.x(), m_max.y()}};
}

const Acts::RectangleBounds& Acts::RectangleBounds::boundingBox() const {
  return (*this);
}

// ostream operator overload
std::ostream& Acts::RectangleBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::RectangleBounds:  (hlX, hlY) = "
     << "(" << 0.5 * (get(eMaxX) - get(eMinX)) << ", "
     << 0.5 * (get(eMaxY) - get(eMinY)) << ")";
  sl << "\n(lower left, upper right):\n";
  sl << min().transpose() << "\n" << max().transpose();
  sl << std::setprecision(-1);
  return sl;
}
