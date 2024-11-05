// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// for definitions of Calibrator, MeasurementSelector
#include "Acts/TrackFinding/CombinatorialKalmanFilterExtensions.hpp"

namespace Acts {
/// to be processed by the CKF
template <typename source_link_iterator_t>
using SourceLinkAccessorDelegate =
    Delegate<std::pair<source_link_iterator_t, source_link_iterator_t>(
        const Surface&)>;

/// Delegate type that retrieves a range of source links to for a given surface

/// expected max number of track states that are expected to be added by
/// stateCandidateCreator
/// @note if the number of states exceeds this number dynamic memory allocation will occur.
///       the number is chosen to yield a container size of 64 bytes.
static constexpr std::size_t s_maxBranchesPerSurface = 10;

namespace CkfTypes {

template <typename T>
using BranchVector = boost::container::small_vector<T, s_maxBranchesPerSurface>;

}  // namespace CkfTypes



/// @brief Get source link range for the given surface and pass to derived class for
///      the actual creation of track states.
///
/// - First get a source link range covering relevant measurements for the given surface
///   This task is delegated to a SourceLinkAccessor.
/// - Then, pass the source link range to the derived class which will do everything else.
template <typename track_state_creator_derived_t,
          typename source_link_iterator_t,
          typename track_container_t>
struct ComposableTrackStateCreator  {
  using SourceLinkAccessor = SourceLinkAccessorDelegate<source_link_iterator_t>;
  using TrackStatesResult =
     Acts::Result<CkfTypes::BranchVector<TrackIndexType>>;
  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  using TrackProxy = typename track_container_t::TrackProxy;
  using TrackStateProxy = typename track_container_t::TrackStateProxy;
  using BoundState = std::tuple<BoundTrackParameters, BoundMatrix, double>;

  SourceLinkAccessor sourceLinkAccessor;
   
  template <typename T_SourceLinkAccessor>
  static SourceLinkAccessor sourceLinkAccessorDelegate(const T_SourceLinkAccessor &sourceLinkAccessor) {
     SourceLinkAccessor delegate;
     delegate.template connect<&T_SourceLinkAccessor::range>(&sourceLinkAccessor);
     return delegate;
  }

  const track_state_creator_derived_t &derived() const
   { return *static_cast<const track_state_creator_derived_t *>(this); }

  /// @brief extend the trajectory onto the given surface.
  ///
  /// @param gctx The geometry context to be used for this task
  /// @param calibrationContext The calibration context used to fill the calibrated data
  /// @param surface The surface onto which the trajectory is extended
  /// @param boundState the predicted bound state on the given surface
  /// @param prevTip the tip of the trajectory which is to be extended
  /// @param bufferTrajectory a temporary buffer which can be used to create
  ///        temporary track states before the selection.
  /// @param trackStateCandidates a temporary buffer which can be used to
  ///        to keep track of newly created temporary track states.
  /// @param trajectory the trajectory to be extended.
  /// @param logger a logger for messages.
  /// 
  /// @return a list of indices of newly created track states which extend the
  ///    trajectory onto the given surface and match the bound state, or an error.
  ///
  /// Extend or branch onto the given surface. This may create new track states
  /// using measurements which match the predicted bound state. This may create
  /// multiple branches. The new track states still miss the "filtered" data.
  Result<CkfTypes::BranchVector<TrackIndexType>> extendTrajectoryOntoSurface (
      const GeometryContext& gctx,
      const CalibrationContext& calibrationContext,
      [[maybe_unused]] const Surface& surface, const BoundState& boundState,
      TrackIndexType prevTip, TrackStateContainerBackend& bufferTrajectory,
      std::vector<TrackStateProxy>& trackStateCandidates,
      TrackStateContainerBackend& trajectory, const Logger& logger) const {

    TrackStatesResult tsRes = TrackStatesResult::success({});
    using SourceLinkRange = decltype(sourceLinkAccessor(surface));
    std::optional<SourceLinkRange>
       slRange = sourceLinkAccessor(surface);
    if (slRange.has_value() && slRange->first != slRange->second) {
       auto [slBegin, slEnd] = *slRange;
       tsRes = derived().createSourceLinkTrackStates(
          gctx, calibrationContext, surface, boundState, slBegin, slEnd,
          prevTip, bufferTrajectory, trackStateCandidates, trajectory, logger);
    }
    return tsRes;
  }
};

/// @brief Create track states for selected measurements from a source link range.
///
/// - First create temporary track states for all measurements defined
///   by a source link range, calibrate the measurements and fill the
///   the calibrated data of these track states.
/// - The measurement selection is delegated to a dedicated measurement selector.
/// - Finally add branches to the given trajectory for the selected, temporary
///   track states. The track states of these branches still lack the filtered
//    data which is to be filled by the next stage e.g. the CombinatorialKalmanFilter.
template <typename source_link_iterator_t, typename track_container_t>
struct TrackStateCandidatorCreatorImpl
   : ComposableTrackStateCreator<TrackStateCandidatorCreatorImpl<source_link_iterator_t, track_container_t>,
                                 source_link_iterator_t,
                                 track_container_t>
{

  using candidate_container_t =
      typename std::vector<typename track_container_t::TrackStateProxy>;
  using MeasurementSelector =
      Delegate<Result<std::pair<typename candidate_container_t::iterator,
                                typename candidate_container_t::iterator>>(
          candidate_container_t& trackStates, bool&, const Logger&)>;

  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  using TrackProxy = typename track_container_t::TrackProxy;
  using TrackStateProxy = typename track_container_t::TrackStateProxy;
  using BoundState = std::tuple<BoundTrackParameters, BoundMatrix, double>;


  typename CombinatorialKalmanFilterExtensions<track_container_t>::Calibrator
      calibrator;
  MeasurementSelector measurementSelector{
      DelegateFuncTag<voidMeasurementSelector>{}};
   
  template < typename measurement_selector_t>
  static MeasurementSelector
  measurementSelectorDelegate(const measurement_selector_t &measurementSelector)
   {
      MeasurementSelector delegate;
      delegate.template connect<&measurement_selector_t::template select<TrackStateContainerBackend>>(&measurementSelector);
      return delegate;
   }
   
  template < typename calibrator_t>
  static typename CombinatorialKalmanFilterExtensions<track_container_t>::Calibrator
  calibratorDelegate(const calibrator_t &calibrator)
   {
      typename CombinatorialKalmanFilterExtensions<track_container_t>::Calibrator
         delegate;
      delegate.template connect<&calibrator_t::calibrate>(&calibrator);
      return delegate;
   }
   
  /// Create track states for selected measurements given by the source links
  ///
  /// @param gctx The current geometry context
  /// @param calibrationContext pointer to the current calibration context
  /// @param surface the surface the sourceLinks are associated to
  /// @param boundState Bound state from the propagation on this surface
  /// @param slBegin Begin iterator for sourceLinks
  /// @param slEnd End iterator for sourceLinks
  /// @param prevTip Index pointing at previous trajectory state (i.e. tip)
  /// @param bufferTrajectory a buffer for temporary candidate track states
  /// @param trackStateCandidates a buffer for temporary track state proxies for candidates
  /// @param trajectory the trajectory to which new track states for selected measurements will be added
  /// @param logger the logger for messages.
  Result<CkfTypes::BranchVector<TrackIndexType>> createSourceLinkTrackStates(
      const GeometryContext& gctx,
      const CalibrationContext& calibrationContext,
      [[maybe_unused]] const Surface& surface, const BoundState& boundState,
      source_link_iterator_t slBegin, source_link_iterator_t slEnd,
      TrackIndexType prevTip, TrackStateContainerBackend& bufferTrajectory,
      std::vector<TrackStateProxy>& trackStateCandidates,
      TrackStateContainerBackend& trajectory, const Logger& logger) const {
    using PM = TrackStatePropMask;

    using ResultTrackStateList =
        Acts::Result<CkfTypes::BranchVector<TrackIndexType>>;
    ResultTrackStateList resultTrackStateList{
        CkfTypes::BranchVector<TrackIndexType>()};
    const auto& [boundParams, jacobian, pathLength] = boundState;

    trackStateCandidates.clear();
    if constexpr (std::ranges::random_access_range<source_link_iterator_t>) {
      trackStateCandidates.reserve(std::distance(slBegin, slEnd));
    }

    // Calibrate all the source links on the surface since the selection has
    // to be done based on calibrated measurement
    for (auto it = slBegin; it != slEnd; ++it) {
      // get the source link
      const auto sourceLink = *it;

      // prepare the track state
      PM mask = PM::Predicted | PM::Jacobian | PM::Calibrated;
      if (it != slBegin) {
        // not the first TrackState, only need uncalibrated and calibrated
        mask = PM::Calibrated;
      }

      ACTS_VERBOSE("Create temp track state with mask: " << mask);
      // CAREFUL! This trackstate has a previous index that is not in this
      // MultiTrajectory Visiting brackwards from this track state will
      // fail!
      auto ts = bufferTrajectory.makeTrackState(mask, prevTip);

      if (it == slBegin) {
        // only set these for first
        ts.predicted() = boundParams.parameters();
        if (boundParams.covariance()) {
          ts.predictedCovariance() = *boundParams.covariance();
        }
        ts.jacobian() = jacobian;
      } else {
        // subsequent track states can reuse
        auto& first = trackStateCandidates.front();
        ts.shareFrom(first, PM::Predicted);
        ts.shareFrom(first, PM::Jacobian);
      }

      ts.pathLength() = pathLength;
      ts.setReferenceSurface(boundParams.referenceSurface().getSharedPtr());

      // now calibrate the track state
      calibrator(gctx, calibrationContext, sourceLink, ts);

      trackStateCandidates.push_back(ts);
    }

    bool isOutlier = false;
    Result<std::pair<typename std::vector<TrackStateProxy>::iterator,
                     typename std::vector<TrackStateProxy>::iterator>>
        selectorResult =
            measurementSelector(trackStateCandidates, isOutlier, logger);
    if (!selectorResult.ok()) {
      ACTS_ERROR("Selection of calibrated measurements failed: "
                 << selectorResult.error());
      resultTrackStateList =
          ResultTrackStateList::failure(selectorResult.error());
    } else {
      auto selectedTrackStateRange = *selectorResult;
      resultTrackStateList = processSelectedTrackStates(
          selectedTrackStateRange.first, selectedTrackStateRange.second,
          trajectory, isOutlier, logger);
    }

    return resultTrackStateList;
  }

  /// Create track states for the given trajectory from candidate track states
  ///
  /// @param begin begin iterator of the list of candidate track states
  /// @param end end iterator of the list of candidate track states
  /// @param trackStates the trajectory to which the new track states are added
  /// @param isOutlier true if the candidate(s) is(are) an outlier(s).
  /// @param logger the logger for messages
  Result<CkfTypes::BranchVector<TrackIndexType>> processSelectedTrackStates(
      typename std::vector<TrackStateProxy>::const_iterator begin,
      typename std::vector<TrackStateProxy>::const_iterator end,
      TrackStateContainerBackend& trackStates, bool isOutlier,
      const Logger& logger) const {
    using PM = TrackStatePropMask;

    using ResultTrackStateList =
        Acts::Result<CkfTypes::BranchVector<TrackIndexType>>;
    ResultTrackStateList resultTrackStateList{
        CkfTypes::BranchVector<TrackIndexType>()};
    CkfTypes::BranchVector<TrackIndexType>& trackStateList =
        *resultTrackStateList;
    trackStateList.reserve(end - begin);

    std::optional<TrackStateProxy> firstTrackState{std::nullopt};
    for (auto it = begin; it != end; ++it) {
      auto& candidateTrackState = *it;

      PM mask = PM::Predicted | PM::Filtered | PM::Jacobian | PM::Calibrated;
      if (it != begin) {
        // subsequent track states don't need storage for these as they will
        // be shared
        mask &= ~PM::Predicted & ~PM::Jacobian;
      }
      if (isOutlier) {
        // outlier won't have separate filtered parameters
        mask &= ~PM::Filtered;
      }

      // copy this trackstate into fitted states MultiTrajectory
      auto trackState =
          trackStates.makeTrackState(mask, candidateTrackState.previous());
      ACTS_VERBOSE("Create SourceLink output track state #"
                   << trackState.index() << " with mask: " << mask);

      if (it != begin) {
        // assign indices pointing to first track state
        trackState.shareFrom(*firstTrackState, PM::Predicted);
        trackState.shareFrom(*firstTrackState, PM::Jacobian);
      } else {
        firstTrackState = trackState;
      }

      // either copy ALL or everything except for predicted and jacobian
      trackState.allocateCalibrated(candidateTrackState.calibratedSize());
      trackState.copyFrom(candidateTrackState, mask, false);

      auto typeFlags = trackState.typeFlags();
      typeFlags.set(TrackStateFlag::ParameterFlag);
      typeFlags.set(TrackStateFlag::MeasurementFlag);
      if (trackState.referenceSurface().surfaceMaterial() != nullptr) {
        typeFlags.set(TrackStateFlag::MaterialFlag);
      }
      if (isOutlier) {
        // propagate information that this is an outlier state
        ACTS_VERBOSE(
            "Creating outlier track state with tip = " << trackState.index());
        typeFlags.set(TrackStateFlag::OutlierFlag);
      }

      trackStateList.push_back(trackState.index());
    }
    return resultTrackStateList;
  }

  /// Default measurement selector which will return all measurements
  /// @param candidates Measurement track state candidates
  static Result<std::pair<typename std::vector<TrackStateProxy>::iterator,
                          typename std::vector<TrackStateProxy>::iterator>>
  voidMeasurementSelector(typename std::vector<TrackStateProxy>& candidates,
                          bool& /*isOutlier*/, const Logger& /*logger*/) {
    return std::pair{candidates.begin(), candidates.end()};
  };
   
};

/// @brief Base class for the ComposableTrackStateCreator which simply delegates its tasks
///
/// A base class for the ComposableTrackStateCreator which simply delegates
/// the track state creation and measurement selection to a dedicated
/// implementation.
template <typename source_link_iterator_t, typename track_container_t>
struct TrackStateCreatorDelegate
   : ComposableTrackStateCreator<TrackStateCreatorDelegate<source_link_iterator_t, track_container_t>,
                                 source_link_iterator_t,
                                 track_container_t>
{
  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  using TrackProxy = typename track_container_t::TrackProxy;
  using TrackStateProxy = typename track_container_t::TrackStateProxy;
  using BoundState = std::tuple<BoundTrackParameters, BoundMatrix, double>;

  /// Delegate definition to create track states for selected measurements
  ///
  /// @note expected to iterator over the given sourceLink range,
  ///       select measurements, and create track states for
  ///       which new tips are to be created, more over the outlier
  ///       flag should be set for states that are outlier.
  ///
  /// @param geoContext The current geometry context
  /// @param calibrationContext pointer to the current calibration context
  /// @param surface the surface at which new track states are to be created
  /// @param boundState the current bound state of the trajectory
  /// @param slBegin Begin iterator for sourceLinks
  /// @param slEnd End iterator for sourceLinks
  /// @param prevTip Index pointing at previous trajectory state (i.e. tip)
  /// @param bufferTrajectory a temporary trajectory which can be used to create temporary track states
  /// @param trackStateCandidates a temporary buffer that can be used to collect track states
  /// @param trajectory the trajectory to which the new states are to be added
  /// @param logger a logger for messages
  using TrackStateCandidateCreator =
      Delegate<Result<CkfTypes::BranchVector<TrackIndexType>>(
          const GeometryContext& geoContext,
          const CalibrationContext& calibrationContext, const Surface& surface,
          const BoundState& boundState, source_link_iterator_t slBegin,
          source_link_iterator_t slEnd, TrackIndexType prevTip,
          TrackStateContainerBackend& bufferTrajectory,
          std::vector<TrackStateProxy>& trackStateCandidates,
          TrackStateContainerBackend& trajectory, const Logger& logger)>;

  TrackStateCandidateCreator trackStateCreator;

  template <typename track_state_candidate_creator_t>
  static TrackStateCandidateCreator
  makeTrackStateCandidateCreator(const track_state_candidate_creator_t &track_state_candidate_creator) {
     TrackStateCandidateCreator delegate;
     delegate.template connect<&track_state_candidate_creator_t
                  ::template createSourceLinkTrackStates<source_link_iterator_t> >(&track_state_candidate_creator);
     return delegate;
  }

  Result<CkfTypes::BranchVector<TrackIndexType>> createSourceLinkTrackStates(
      const GeometryContext& gctx,
      const CalibrationContext& calibrationContext,
      [[maybe_unused]] const Surface& surface, const BoundState& boundState,
      source_link_iterator_t slBegin, source_link_iterator_t slEnd,
      TrackIndexType prevTip, TrackStateContainerBackend& bufferTrajectory,
      std::vector<TrackStateProxy>& trackStateCandidates,
      TrackStateContainerBackend& trajectory, const Logger& logger) const {
     return trackStateCreator(gctx, calibrationContext,
          surface, boundState, slBegin, slEnd, prevTip, bufferTrajectory, trackStateCandidates,
          trajectory, logger);
  }

};

}
