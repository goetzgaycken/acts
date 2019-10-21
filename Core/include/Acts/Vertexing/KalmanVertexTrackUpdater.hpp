// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/KalmanVertexUpdater.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

/// @class KalmanVertexTrackUpdater
///
/// @brief Refits a single track with the knowledge of
/// the vertex it has originated from
/// Based on R. Frühwirth et al.
/// Vertex reconstruction and track bundling at the lep collider using
/// robust Algorithms Computer Physics Comm.: 96 (1996) 189, chapter 2.1
///
/// @tparam input_track_t Track object type
template <typename input_track_t>
class KalmanVertexTrackUpdater {
 public:
  /// @struct Configuration struct
  struct Config {
    /// Kalman vertex updater
    KalmanVertexUpdater<input_track_t> vtx_updater;
  };

  /// Constructor
  KalmanVertexTrackUpdater(const Config& config = Config()) : m_cfg(config) {}

  /// @brief Refits a single track with the knowledge of
  /// the vertex it has originated from
  ///
  /// @param gctx The Geometry Context
  /// @param track Track to update
  /// @param vtx Vertex `track` belongs to
  Result<void> update(const GeometryContext& gctx,
                      TrackAtVertex<input_track_t>& track,
                      const Vertex<input_track_t>* vtx) const;

 private:
  /// Configuration object
  const Config m_cfg;
};

}  // Namespace Acts

#include "Acts/Vertexing/KalmanVertexTrackUpdater.ipp"
