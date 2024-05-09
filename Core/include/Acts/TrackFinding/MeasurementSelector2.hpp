#pragma once

// Alternative measurement selector
//
// This measurement selector is assuming the following
// - the number of selected measurements is small (~<10)
// - the number of measurement candidates is typically >> the
//   number of selected measuremebts
// - the total number of candidate measurements can be large
// - there is a simple one-to-one relation between bound state
//   parameters and measurement coordinates.

#include "Acts/Utilities/MultiDimVariant.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"

// for BaseTypes
#include "Acts/EventData/TrackStateProxy.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/EventData/TrackParameters.hpp"

// for MeasurementSizeMax
#include "Acts/EventData/MultiTrajectory.hpp"

template <typename traj_t>
struct BaseTypes {
   template <std::size_t N>
   using Measurement = typename Acts::detail_lt::Types<N>::Coefficients;

   template <std::size_t N>
   using MeasurementCovariance = typename Acts::detail_lt::Types<N>::Covariance;

   using trajectory_t = traj_t;
   using TrackStateProxy = typename traj_t::TrackStateProxy;
   using MatrixFloatType = Acts::ActsScalar;
   //   using source_link_iterator_t = typename traj_t::source_link_iterator_t;
};

template <class T_BaseTypes>
struct MeasurementData {
   // template <std::size_t M>
   // using PairType = std::pair< typename T_BaseTypes::template Measurement<M>, typename T_BaseTypes::template MeasurementCovariance<M> >;

   // helper struct to aid template parameter deduction (@TODO better solution ?)
   template <std::size_t M>
   struct MeasCovPair : std::pair< typename T_BaseTypes::template Measurement<M>, typename T_BaseTypes::template MeasurementCovariance<M> >
   {
      //      using std::pair< typename T_BaseTypes::template Measurement<M>, typename T_BaseTypes::template MeasurementCovariance<M> >::pair;
      MeasCovPair &operator=(const std::pair< typename T_BaseTypes::template Measurement<M>, typename T_BaseTypes::template MeasurementCovariance<M> > &a) {
         this->PairType.operator=(a);
         return *this;
      }
      MeasCovPair &operator=(std::pair< typename T_BaseTypes::template Measurement<M>, typename T_BaseTypes::template MeasurementCovariance<M> > &&a) {
         this->PairType.operator=(std::move(a));
         return *this;
      }
      MeasCovPair() = default;
      MeasCovPair(const std::pair< typename T_BaseTypes::template Measurement<M>, typename T_BaseTypes::template MeasurementCovariance<M> > &a)
         : std::pair< typename T_BaseTypes::template Measurement<M>, typename T_BaseTypes::template MeasurementCovariance<M> >(a)
      {}
      MeasCovPair(std::pair< typename T_BaseTypes::template Measurement<M>, typename T_BaseTypes::template MeasurementCovariance<M> > &&a)
         : std::pair< typename T_BaseTypes::template Measurement<M>, typename T_BaseTypes::template MeasurementCovariance<M> >(std::move(a))
      {}
   };

   template <std::size_t M>
   using type = MeasCovPair<M> ;
};

// the variant allowing measurements :
// std::variant < (Meas<1>,MeasCov<1>), (Meas<2>,MeasCov<2>), ... , (Meas<N>,MeasCov<N>)  >
// with ( , )      == pair
//      Meas<N>    == T_MultiTrajectory::TrackStateProxy:: Measurement<N>
//      MeasCov<N> == T_MultiTrajectory::TrackStateProxy:: MeasurementCovariance<N>
template <std::size_t N, class T_BaseTypes>
using MultiDimMeasurementData = typename MultiDimVariant::MakeMultiDimVariant< MeasurementData<T_BaseTypes>, N>::variant_type;

// Map from the measurement to bound state domain
// it is assumed that there is a simple unambiguous association
// between coordinates of the measurement domain and the bound state domain
// e.g. measurement coordinate 1 maps to loc0 of the bound state.
struct ParameterMapping {
   template <std::size_t M>
   using type = std::array<unsigned char, M>;
};

// variant to hold maps for mapping measurement coordinates to bound parameters for measurements with 1,2, ... N dimensions
// std::variant< array<uchar_t, 1> , array<uchar_t, 2>, ... array<uchar_t, N> >
template <std::size_t N>
using MultiDimParameterMap = typename MultiDimVariant::MakeMultiDimVariant<ParameterMapping, N>::variant_type;

/// matrix adapter for Eigen
template <class T_Matrix>
constexpr std::size_t matrixColumns() { return T_Matrix::ColsAtCompileTime;}
template <class T_Matrix>
constexpr std::size_t matrixRows() { return T_Matrix::RowsAtCompileTime;}

template <typename T_Float, class T_Matrix>
auto matrixTypeCast(const T_Matrix &matrix) { return matrix.template cast<T_Float>(); }

template <class T_Matrix>
auto transpose(const T_Matrix &matrix) { return matrix.transpose(); }

template <class T_Matrix>
auto invert(const T_Matrix &matrix) { return matrix.inverse(); }


// template <std::size_t N, class T_Matrix>
// constexpr std::size_t destinationMatrixColumns() { return T_Matrix::ColsAtCompileTime==1 ? 1 : N;}

// utility to "project" a bound state parameter vector or covariance matrix onto the 1,2, ... N dimensional measurement domain
// #TODO allow to influence resulting matrix type ?
template <std::size_t N,class T_ResultType,class T_Matrix>
T_ResultType project(ParameterMapping::type<N> parameter_map, const T_Matrix &matrix)
{
   using MatrixIndexMapType = unsigned char; // "char" to reduce the size of the map, and if not wide enough this entire
                                             //        concept is likely inefficient.
   using MatrixIndexType = unsigned int;     // @TODO or std::size_t ? does not matter

   // ensure that index types are wide enough
   static_assert( matrixRows<T_Matrix>() < std::numeric_limits<MatrixIndexMapType>::max());
   static_assert( N*matrixRows<T_Matrix>() < std::numeric_limits<MatrixIndexType>::max());

   T_ResultType ret;
   if constexpr(matrixColumns<T_Matrix>() == 1) {
      // handle projection of paramteter vector
      for (MatrixIndexType meas_i=0; meas_i<N; ++meas_i) {
         assert( meas_i < parameter_map.size() );
         ret(meas_i,0) = matrix( parameter_map[meas_i], 0);
      }
   }
   else {
      // handle projection of covariance matrix
      // "project" matrix
      for (MatrixIndexType meas_i=0; meas_i<N; ++meas_i) {
         assert( meas_i < parameter_map.size());
         MatrixIndexType param_i = parameter_map[meas_i];
         for (MatrixIndexType meas_j=0; meas_j<N; ++meas_j) {
            assert( meas_j < parameter_map.size());
            ret(meas_i,meas_j) = matrix(param_i, parameter_map[meas_j]);
         }
      }
   }
   return ret;
}

// Helper struct for mapping a bound state to one of the measurement variants.
// The struct holds a bound state and will produce measurement, covariance pairs
// for a given map from the measurement to the bound state domain.
template <class T_BaseTypes, class T_BoundParam, std::size_t DIMMAX>
struct ProjectionHelper {
   // function to produce a measurement and covariance pair from a bound state for a certain mapping
   // the function can be applied to a 1,2, ... or n-dimensional parameter map and will create a measurement
   // covariance pair matching the dimension.
   template <std::size_t N>
   MultiDimMeasurementData<DIMMAX, T_BaseTypes> operator()(const ParameterMapping::type<N> &parameter_mapping) const {
      MultiDimMeasurementData<DIMMAX, T_BaseTypes>
         ret = std::make_pair(project<N,typename T_BaseTypes::template Measurement<N> >(parameter_mapping, m_boundParam->parameters()),
                              project<N, typename T_BaseTypes::template MeasurementCovariance<N> >(parameter_mapping, m_boundParam->covariance().value()));
      return ret;
   }
   const T_BoundParam *m_boundParam = nullptr;
};

// helper method to create a measurement, covariance pair for a given map from measurement to parameter space
// the dimensionality of the returned pair will correspond to the dimensionality of the given parameter map
template <std::size_t DIMMAX,
          class T_BaseTypes,
          class T_BoundParam>
inline  MultiDimMeasurementData<DIMMAX, T_BaseTypes> makePredicted( const MultiDimParameterMap<DIMMAX> &parameter_map,
                                                                          const T_BoundParam &bound_param) {
   // call projection from bound to N-dim measurement, where N depends on the dimension of the projection matrix
   ProjectionHelper<T_BaseTypes, T_BoundParam, DIMMAX> projection_helper{ &bound_param };
   MultiDimMeasurementData<DIMMAX, T_BaseTypes>
      ret = std::visit< MultiDimMeasurementData<DIMMAX, T_BaseTypes>  > ( projection_helper, parameter_map);
   return ret;
}

// helper method to compute a chi2 for the difference of two "measurement" and covariance pairs
template <class T_BaseTypes, std::size_t N>
double computeChi2(const typename T_BaseTypes::template Measurement<N> &a,
                   const typename T_BaseTypes::template MeasurementCovariance<N> &a_cov,
                   const typename T_BaseTypes::template Measurement<N> &b,
                   const typename T_BaseTypes::template MeasurementCovariance<N> &b_cov) {
   typename T_BaseTypes::template MeasurementCovariance<N> inv_ab_cov( invert(a_cov+b_cov) );
   typename T_BaseTypes::template Measurement<N> diff( a-b);
   return (transpose(diff) * inv_ab_cov * diff)(0,0);
}

// helper struct to compute the chi2 for measurements of 1,2,... or N dimensions.
// It is assumed that the dimensionality of the predicted position is known and
// matched by the provided measurement covariance pairs.
template <class T_BaseTypes, std::size_t N_Pred>
struct Chi2Functor {
   template <std::size_t N>
   double operator()(const typename MeasurementData<T_BaseTypes>::template type<N> &val) {
      if constexpr( N == N_Pred) {
         return computeChi2<T_BaseTypes, N_Pred>(val.first, val.second, m_predicted->first, m_predicted->second);
      }
      else {
         assert( N == N_Pred );
         return 0.;
      }
   }
   const typename MeasurementData<T_BaseTypes>::template type<N_Pred> *m_predicted;
};

// Collection to hold the n-"best" candidates
// The objects of type PayloadType must support assignment operation, and must
// be default constructible. Moreover it must be possible to provide a
// "comparison" operator to order the payload objects.
template <std::size_t N, class PayloadType >
struct TopCollection {
   using IndexType = unsigned short; // @TODO or char ? If N>>10 this concept is likely
                                     //       inefficient
   //   using PayloadType = Payload<DIM>;
   TopCollection(std::size_t max_n) {
      init(max_n);
   }

   // @param max_n the maximum number of top-candidates is fixed by the template parameter
   //    N but can be reduced further to this number
   void init(std::size_t max_n) {
      assert( max_n < N);
      m_nextSlot=0;
      m_maxSlots=max_n;
      m_order[0]=0;
   }
   // @param get a slot to hold a new candidate which is not necessarily accepted in the list
   //     the n-top candidates
   PayloadType &slot() {
      return m_slots[m_order[m_nextSlot] ];
   }
   // @param idx get the specified filled slot (read only) indicated by the index, where the index does not
   //    indicate the order in the top-candidate list
   const PayloadType &getSlot(IndexType idx) const {
      return m_slots[idx];
   }
   // @param idx get the specified filled slot indicated by the index, where the index does not
   //    indicate the order in the top-candidate list
   PayloadType &getSlot(IndexType idx) {
      return m_slots[idx];
   }
   // @param test whether the given index points to one of the accepted top candidates
   bool isValid(IndexType idx) const {
      return idx < m_nextSlot;
   }
   // Accept the element of the latest slot provided there is still a free slot or it is better than
   // the worst element.
   // @param comparison operator to compute e.g. a<b where "smaller" means "better"
   void acceptAndSort(std::function<bool(const PayloadType &a, const PayloadType &b)> comparison) {
      // bubble best element to top
      for (unsigned int slot_i = m_nextSlot;
           slot_i-- > 0
           && !comparison(m_slots[ m_order[slot_i] ],m_slots[ m_order[slot_i+1] ]);) {
         std::swap(m_order[slot_i],m_order[slot_i+1]);
      }
      // if there are still free slot increase the number of used slots
      if (m_nextSlot < m_maxSlots) {
         ++m_nextSlot;
         m_order[m_nextSlot]=m_nextSlot;
      }
   }

   bool empty() const {
      return m_nextSlot==0;
   }

   // helper to iterate over the slot indices in sorting order from best to worst
   typename std::array<IndexType, N+1>::const_iterator begin() { return m_order.begin(); }
   typename std::array<IndexType, N+1>::const_iterator end()   { return m_order.begin()+m_nextSlot; }

   std::array<PayloadType, N+1> m_slots;         // storage for the slots
   std::array<IndexType,   N+1> m_order;         // order of the filled slots
   IndexType                    m_nextSlot = 0;  // the index of the next free slot
   IndexType                    m_maxSlots = 0;  // maximum number of top-slots
};

// decorated measurement
// decorate measurement and its covariance with a chi2 and an outlier flag
template <std::size_t DIMMAX, class T_BaseTypes>
struct MatchingMeasurement {
   MultiDimMeasurementData<DIMMAX, T_BaseTypes> m_measurement;
   std::optional<Acts::SourceLink>              m_sourceLink;
   float                                        m_chi2;
   bool m_isOutLier;
};

// Extensions for the MeasurementSelector
template <std::size_t DIMMAX, class T_BaseTypes>
struct MeasurementSelectorExtensions {
   using Projector =
      Acts::Delegate<MultiDimParameterMap<DIMMAX>(const Acts::GeometryContext&,
                                                  const Acts::CalibrationContext *,
                                                  const Acts::Surface &surface )>;

   using MeasurementGetter =
      Acts::Delegate<MultiDimMeasurementData<DIMMAX,T_BaseTypes>(const Acts::GeometryContext&,
                                                                 const Acts::CalibrationContext *,
                                                                 const Acts::SourceLink&)>;

   Projector projector;
   MeasurementGetter measurementGetter;
};

//loop over given measurements, and select the n-best measurements compatible with the prediction
template <std::size_t NMeasMax, std::size_t DIMMAX, class T_BaseTypes, typename source_link_iterator_t>
struct MeasurementLoop {
   using  MeasurementGetter=
      Acts::Delegate<MultiDimMeasurementData<DIMMAX,T_BaseTypes>(const Acts::GeometryContext&,
                                                                 const Acts::CalibrationContext *,
                                                                 const Acts::SourceLink&)>;

   source_link_iterator_t soureLinksBegin;
   source_link_iterator_t sourceLinksEnd;
   double chi2Cut = std::numeric_limits<double>::max();
   const Acts::GeometryContext* geometryContext;
   const Acts::CalibrationContext* calibrationContext;
   typename MeasurementSelectorExtensions<DIMMAX, T_BaseTypes>::MeasurementGetter measurementGetter;

   mutable TopCollection<NMeasMax, MatchingMeasurement<DIMMAX, T_BaseTypes> > *selectedMeasurements = nullptr;

   template <std::size_t N>
   void operator()(const  typename MeasurementData<T_BaseTypes>::template type<N> &predicted) const {
      Chi2Functor<T_BaseTypes, N> chi2_func  { &predicted };
      for (auto sourceLinkIter = soureLinksBegin; sourceLinkIter != sourceLinksEnd; ++sourceLinkIter) {
        // get the source link
        //MultiDimMeasurementData<DIMMAX,T_BaseTypes>
        Acts::SourceLink source_link = *sourceLinkIter;
        MatchingMeasurement<DIMMAX, T_BaseTypes>
           &matching_measurement =selectedMeasurements->slot();
        matching_measurement.m_measurement= measurementGetter(*geometryContext, calibrationContext,source_link);
        matching_measurement.m_chi2 = std::visit<double>( chi2_func, matching_measurement.m_measurement);
        if (matching_measurement.m_chi2<chi2Cut) {
           matching_measurement.m_sourceLink=source_link;
           selectedMeasurements->acceptAndSort([](const MatchingMeasurement<DIMMAX,T_BaseTypes> &a,
                                                  const MatchingMeasurement<DIMMAX,T_BaseTypes> &b) {
              return a.m_chi2 < b.m_chi2;
           });
        }
      }
   }
};

struct MeasurementSelector2Cuts {
  /// bins in |eta| to specify variable selections
  std::vector<float> etaBins{};
  /// Maximum local chi2 contribution.
  std::vector<std::pair<float, float> > chi2CutOff{ {15,25} };
  /// Maximum number of associated measurements on a single surface.
  std::vector<std::size_t> numMeasurementsCutOff{1};
};

template <std::size_t NMeasMax,
          std::size_t DIMMAX,
          class T_BaseTypes,
          class T_BoundState = std::tuple<Acts::BoundTrackParameters, Acts::BoundMatrix, double> >
struct MeasurementSelector2 {
   using Config = Acts::GeometryHierarchyMap<MeasurementSelector2Cuts>;
   Config m_config;
   MeasurementSelectorExtensions<DIMMAX, T_BaseTypes> m_extensions;

   struct ProjectorBitSetMaker {
      template <std::size_t N>
      Acts::ProjectorBitset operator()(const ParameterMapping::type<N> &parameter_map) const {
         constexpr std::size_t nrows = Acts::MultiTrajectoryTraits::MeasurementSizeMax;
         constexpr std::size_t ncols = Acts::eBoundSize;

         std::bitset<nrows * ncols> proj_bitset {};

         for (unsigned int col_i=0; col_i<N; ++col_i) {
            unsigned int row_i = parameter_map[col_i];
            unsigned int idx = col_i *nrows + row_i;      // @TODO row major or column major ?
            proj_bitset[ (nrows * ncols - 1) - idx ] = 1;
         }
         return proj_bitset.to_ullong();
      }

      static Acts::ProjectorBitset create(const MultiDimParameterMap<DIMMAX> &parameter_map) {
          ProjectorBitSetMaker projectorBitSetMaker;
          return std::visit<Acts::ProjectorBitset>(projectorBitSetMaker, parameter_map);
       }

   };
   // Measurement selector which returns the n-most compatible measurements
   // Loop over the given measurements, calibrate them and compute a chi2
   // from distance to the prediction given by the bound state.
   // Keep the n-best measurements and flag them as outliers if the chi2
   // is too large.
   template <typename source_link_iterator_t>
   Acts::Result<TopCollection<NMeasMax, MatchingMeasurement<DIMMAX, T_BaseTypes> > >
   compatibleMeasurements(const Acts::GeometryContext& geometryContext,
                          const Acts::CalibrationContext* calibrationContext,
                          const Acts::Surface& surface,
                          const T_BoundState& boundState,
                          source_link_iterator_t sourceLinkBegin,
                          source_link_iterator_t sourceLinkEnd,
                          const Acts::Logger& logger,
                          MultiDimParameterMap<DIMMAX> &parameter_map) const {

      Acts::Result< TopCollection<NMeasMax, MatchingMeasurement<DIMMAX, T_BaseTypes> > >
         result{ TopCollection<NMeasMax, MatchingMeasurement<DIMMAX, T_BaseTypes> >(0) };
      if (sourceLinkBegin != sourceLinkEnd) { // @TODO already tested by the caller ...

         // Get geoID of this surface
         auto geoID = surface.geometryId();
         // Find the appropriate cuts
         auto cuts = m_config.find(geoID);
         if (cuts == m_config.end()) {
            // for now we consider missing cuts an unrecoverable error
            // TODO consider other options e.g. do not add measurements at all (not
            // even as outliers)
            return Acts::CombinatorialKalmanFilterError::MeasurementSelectionFailed;
         }

         //         const std::vector<std::pair<float, float> >& chi2CutOff = cuts->chi2CutOff;
         std::size_t eta_bin = getEtaBin(std::get<0>(boundState), cuts->etaBins);
         const std::pair<float,float> maxChi2Cut = ! cuts->chi2CutOff.empty()
            ? cuts->chi2CutOff[ std::min(cuts->chi2CutOff.size()-1, eta_bin) ]
            : std::make_pair<float,float>(std::numeric_limits<float>::max(),
                                          std::numeric_limits<float>::max());
         const std::size_t numMeasurementsCut = (!cuts->numMeasurementsCutOff.empty()
                                                 ? cuts->numMeasurementsCutOff[ std::min(cuts->numMeasurementsCutOff.size()-1, eta_bin) ]
                                                 : NMeasMax);
         ACTS_VERBOSE("Get cut for eta-bin="
                      << (eta_bin < cuts->etaBins.size() ? cuts->etaBins[eta_bin] : std::numeric_limits<float>::max())
                      << ": chi2 (max,max-outlier) " << maxChi2Cut.first << ", " << maxChi2Cut.second
                      << " max.meas." << numMeasurementsCut);

         parameter_map = m_extensions.projector(geometryContext,
                                                calibrationContext,
                                                surface);

         MultiDimMeasurementData<DIMMAX,T_BaseTypes>
            predicted( makePredicted<DIMMAX,T_BaseTypes>( parameter_map, std::get<0>(boundState)));

         TopCollection<NMeasMax, MatchingMeasurement<DIMMAX, T_BaseTypes> > &compatible_measurements = *result;
         compatible_measurements.init(numMeasurementsCut);

         std::visit<void>( MeasurementLoop<NMeasMax,DIMMAX,T_BaseTypes, source_link_iterator_t>
                                    {sourceLinkBegin,
                                     sourceLinkEnd,
                                     maxChi2Cut.second, // also accept outliers
                                     &geometryContext,
                                     calibrationContext,
                                     m_extensions.measurementGetter,
                                     &compatible_measurements },
                           predicted);

         // flag outlier
         for (auto idx : compatible_measurements) {
            compatible_measurements.getSlot(idx).m_isOutLier
               = compatible_measurements.getSlot(idx).m_chi2 >maxChi2Cut.first;
         }

      }
      return result;
   }

   template <typename parameters_t>
   static std::size_t getEtaBin(const parameters_t& boundState,
                                const std::vector<float> &etaBins) {
      if (etaBins.empty()) {
         return 0u;  // shortcut if no etaBins
      }
      const float eta = std::abs(std::atanh(std::cos(boundState.parameters()[Acts::eBoundTheta])));
      std::size_t bin = 0;
      for (auto etaBin : etaBins) {
         if (etaBin >= eta) {
            break;
         }
         bin++;
      }
      return bin;
   }


   struct MeasurementSetter {
      mutable typename T_BaseTypes::TrackStateProxy *trackState;
      template <std::size_t DIM>
      void operator()(const typename MeasurementData<T_BaseTypes>::template type<DIM> &measurement_data) const {
         trackState->allocateCalibrated(DIM);
         trackState->template calibrated<DIM>()           = matrixTypeCast<typename T_BaseTypes::MatrixFloatType>(measurement_data.first);
         trackState->template calibratedCovariance<DIM>() = matrixTypeCast<typename T_BaseTypes::MatrixFloatType>(measurement_data.second);
      }
   };

   template <typename source_link_iterator_t>
   Acts::Result<std::pair<
                   typename std::vector<typename T_BaseTypes::TrackStateProxy>::iterator,
                   typename std::vector<typename T_BaseTypes::TrackStateProxy>::iterator>>
   createSourceLinkTrackStates(const Acts::GeometryContext& geoContext,
                               const Acts::CalibrationContext* calibrationContext,
                               const Acts::Surface& surface,
                               const T_BoundState& boundState,
                               source_link_iterator_t sourceLinkBegin,
                               source_link_iterator_t sourceLinkEnd,
                               std::size_t prevTip,
                               typename T_BaseTypes::trajectory_t& trajectory,
                               std::vector<typename T_BaseTypes::TrackStateProxy> &trackStateCandidates,
                               bool &isOutlier,
                               const Acts::Logger& logger) const {
      MultiDimParameterMap<DIMMAX> parameter_map;
      Acts::Result<TopCollection<NMeasMax, MatchingMeasurement<DIMMAX, T_BaseTypes> > >
          compatible_measurements = compatibleMeasurements(geoContext,
                                                           calibrationContext,
                                                           surface,
                                                           boundState,
                                                           sourceLinkBegin,
                                                           sourceLinkEnd,
                                                           logger,
                                                           parameter_map);
       trackStateCandidates.clear();
       if (!compatible_measurements.ok()) {
         return Acts::CombinatorialKalmanFilterError::MeasurementSelectionFailed;
       }
       else if (compatible_measurements.ok() && !compatible_measurements->empty()) {
          if constexpr (std::is_same_v<
                        typename std::iterator_traits<
                        source_link_iterator_t>::iterator_category,
                        std::random_access_iterator_tag>) {
             unsigned int n_measurements
                = std::max(std::count_if ( compatible_measurements->begin(),
                                           compatible_measurements->end(),
                                           [&compatible_measurements]
                                           (typename TopCollection<NMeasMax, MatchingMeasurement<DIMMAX, T_BaseTypes> >::IndexType idx) {
                                              return !compatible_measurements->getSlot(idx).m_isOutLier;
                                           }),
                           1u);
             trackStateCandidates.reserve(n_measurements);
          }
       }
       trajectory.clear();
       using PM = Acts::TrackStatePropMask;
       Acts::ProjectorBitset projector_bitset = ProjectorBitSetMaker::create(parameter_map);

       // create track states for all compatible measurement candidates
       const auto &boundParams = std::get<0>(boundState);
       const auto &pathLength = std::get<2>(boundState);

       bool first_state=true;
       for (typename TopCollection<NMeasMax, MatchingMeasurement<DIMMAX, T_BaseTypes> >::IndexType
               idx: *compatible_measurements) {
          // prepare the track state
          PM mask = PM::Calibrated;

          if (first_state) {

             // not the first TrackState, only need uncalibrated and calibrated
             mask = PM::Predicted | PM::Jacobian | PM::Calibrated;
             isOutlier = compatible_measurements->getSlot(idx).m_isOutLier;
          }
          else if (compatible_measurements->getSlot(idx).m_isOutLier) {
             break;
          }
          ACTS_VERBOSE("Create temp track state with mask: " << mask);
          // CAREFUL! This trackstate has a previous index that is not in this
          // MultiTrajectory Visiting brackwards from this track state will
          // fail!
          typename T_BaseTypes::TrackStateProxy ts = trajectory.makeTrackState(mask, prevTip);
          if (first_state) {
             // only set these for first
             ts.predicted() = boundParams.parameters();
             if (boundParams.covariance()) {
                ts.predictedCovariance() = *boundParams.covariance();
             }
             ts.jacobian() = std::get<1>(boundState);
             first_state=false;
          } else {
             // subsequent track states can reuse
             auto& first = trackStateCandidates.front();
             ts.shareFrom(first, PM::Predicted);
             ts.shareFrom(first, PM::Jacobian);
          }
          ts.pathLength() = pathLength;

          ts.setReferenceSurface(boundParams.referenceSurface().getSharedPtr());
          ts.setProjectorBitset(projector_bitset); // @TODO needed ?
          ts.setUncalibratedSourceLink(compatible_measurements->getSlot(idx).m_sourceLink.value());

          // @TODO set sourcelink
          std::visit<void>( MeasurementSetter{&ts}, compatible_measurements->getSlot(idx).m_measurement);
          trackStateCandidates.push_back(ts);
       }
       return
       (trackStateCandidates.begin() != trackStateCandidates.end())
        ? Acts::Result<std::pair<
                          typename std::vector<typename T_BaseTypes::TrackStateProxy>::iterator,
                          typename std::vector<typename T_BaseTypes::TrackStateProxy>::iterator> >
             ::success(std::make_pair(trackStateCandidates.begin(), trackStateCandidates.end()))
       : Acts::Result<std::pair<
                          typename std::vector<typename T_BaseTypes::TrackStateProxy>::iterator,
                          typename std::vector<typename T_BaseTypes::TrackStateProxy>::iterator> >
             ::failure(Acts::CombinatorialKalmanFilterError::MeasurementSelectionFailed); // @TODO no compatible candidates == error ?
    }


};
