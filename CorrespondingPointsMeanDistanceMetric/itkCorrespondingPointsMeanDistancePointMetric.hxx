/*=========================================================================
 *
 *  Copyright UMC Utrecht and contributors
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#pragma once

#ifndef __itkCorrespondingPointsMeanDistancePointMetric_hxx
#define __itkCorrespondingPointsMeanDistancePointMetric_hxx

#include "itkCorrespondingPointsMeanDistancePointMetric.h"

namespace itk
{

/**
 * ******************* Constructor *******************
 */

template< class TFixedPointSet, class TMovingPointSet >
CorrespondingPointsMeanDistancePointMetric< TFixedPointSet, TMovingPointSet >
::CorrespondingPointsMeanDistancePointMetric()
{} // end Constructor


/**
 * ******************* Initialize *******************
 */

template< class TFixedPointSet, class TMovingPointSet >
void
CorrespondingPointsMeanDistancePointMetric< TFixedPointSet, TMovingPointSet >
::Initialize( void ) throw ( ExceptionObject )
{
  /** Clean state from previous resolution */
  m_LastDimensionValues.clear();
  m_Landmarks.clear();

  /** Initialize transform, interpolator, etc. */
  Superclass::Initialize();

  MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();
  if( !movingPointSet )
  {
    itkExceptionMacro( << "Moving point set has not been assigned" );
  }

  /** Retrieve slowest varying dimension */
  m_LastDimension     = movingPointSet->PointDimension - 1;

  /** Retrieve observed size of the slowest varying dimension */
  /** Loop over the corresponding points. */
  for( auto p = movingPointSet->GetPoints()->Begin(); p != movingPointSet->GetPoints()->End(); ++p )
  {
    m_LastDimensionValues.insert( p.Value()[ m_LastDimension ] );
  }

  /** Create a map between the last dimension channel and the set of
   * points on that channel */
  std::map<unsigned int, std::vector< MovingPointType > > pointsPerImage;
  for ( auto const v : m_LastDimensionValues ) {
    pointsPerImage[ v ] = {};
  }

  /** Populate the point set for each image */
  for( auto p = movingPointSet->GetPoints()->Begin(); p != movingPointSet->GetPoints()->End(); ++p )
  {
    auto t = p.Value()[ m_LastDimension ];
    pointsPerImage[ t ].push_back( p.Value() );
  }

  m_NumberOfPointsPerImage = 0;
  for ( auto const & pair : pointsPerImage ) {
    /** Validate the number of points */
    if (m_NumberOfPointsPerImage == 0) {
      m_NumberOfPointsPerImage = pair.second.size();
      m_Landmarks.resize(m_NumberOfPointsPerImage);
    }
    else {
      if (m_NumberOfPointsPerImage != pair.second.size()) {
        itkExceptionMacro( << "Different number of corresponding points in the stacked images" );
      }
    }

    /** Fill one array with all the correspondences of the same point
     * across all images */
    for ( size_t i = 0; i < pair.second.size(); ++i ) {
      auto p = pair.second[ i ];
      m_Landmarks[ i ].occurrences.push_back( p );
    }
  }

  /** Compute the mean for the occurrences of each landmark */
  for ( auto & landmark : m_Landmarks )
  {
    landmark.mean.Fill( NumericTraits< MovingPointValueType >::ZeroValue() );
    for ( const auto & p : landmark.occurrences )
    {
      landmark.mean.GetVnlVector() += p.GetVnlVector();
    }
    landmark.mean.GetVnlVector() /= landmark.occurrences.size();
  }

  m_NumberOfLandmarks = m_Landmarks.size();

} // end Initialize()

/**
 * ******************* GetValue *******************
 */

template< class TFixedPointSet, class TMovingPointSet >
typename CorrespondingPointsMeanDistancePointMetric< TFixedPointSet, TMovingPointSet >::MeasureType
CorrespondingPointsMeanDistancePointMetric< TFixedPointSet, TMovingPointSet >
::GetValue( const TransformParametersType & parameters ) const
{
  /** TODO
   */
  MeasureType measure = NumericTraits< MeasureType >::Zero;
  DerivativeType dummyvalue = DerivativeType( this->GetNumberOfParameters() );
  dummyvalue.Fill( NumericTraits< DerivativeValueType >::ZeroValue() );

  this->GetValueAndDerivative( parameters, measure, dummyvalue );

  return measure;

} // end GetValue()


/**
 * ******************* GetDerivative *******************
 */

template< class TFixedPointSet, class TMovingPointSet >
void
CorrespondingPointsMeanDistancePointMetric< TFixedPointSet, TMovingPointSet >
::GetDerivative( const TransformParametersType & parameters,
  DerivativeType & derivative ) const
{
  /** When the derivative is calculated, all information for calculating
   * the metric value is available. It does not cost anything to calculate
   * the metric value now. Therefore, we have chosen to only implement the
   * GetValueAndDerivative(), supplying it with a dummy value variable.
   */
  MeasureType dummyvalue = NumericTraits< MeasureType >::Zero;
  this->GetValueAndDerivative( parameters, dummyvalue, derivative );

} // end GetDerivative()


/**
 * ******************* GetValueAndDerivative *******************
 */


template< class TFixedPointSet, class TMovingPointSet >
void
CorrespondingPointsMeanDistancePointMetric< TFixedPointSet, TMovingPointSet >
::GetValueAndDerivative( const TransformParametersType & parameters,
  MeasureType & value, DerivativeType & derivative ) const
{
  /***************************************************************************
   *                                                                         *
   * METRIC DESCRIPTION                                                      *
   * ------------------                                                      *
   *                                                                         *
   * Align the landmarks on their mean position.  Minimise the mean          *
   * Euclidean distance between the mapped mean position of each             *
   * landmark and all its occurrences in the moving images                   *
   *                                                                         *
   *   E(mu) = 1/n Σ_i Σ_j [ ( p_ij - T( 1/m Σ_k p_ik ) ]                    *
   *                                                                         *
   * where p_i is the i-th landmark, p_ij is the position of landmark p_i    *
   * in the j-th image, n is the number of landmarks, m the number of        *
   * subjects, T the transform with parameters mu.                           *
   *                                                                         *
   * STRUCTURE OF THE MOVING POINTS FILE                                     *
   * -----------------------------------                                     *
   *                                                                         *
   * Only the moving points file is used, the fixed must be provided only    *
   * as a dummy to prevent errors, but it is not used in this metric.        *
   *                                                                         *
   * It is important to store the moving points in the correct order.        *
   * Given `m` moving images (stored in `m` channels in the time axis),      *
   * the points must be sorted in the file as                                *
   *                                                                         *
   *   x00 y00 z00 0  // image 0, landmark 0                                 *
   *   x01 y01 z01 0  // image 0, landmark 1                                 *
   *        ...                                                              *
   *   x0n y0n z0n 0  // image 0, landmark n                                 *
   *                                                                         *
   *   x10 y10 z10 1  // image 1, landmark 0                                 *
   *   x11 y11 z11 1  // image 1, landmark 1                                 *
   *        ...                                                              *
   *   x1n y1n z1n 1  // image 1, landmark n                                 *
   *     .........                                                           *
   *   xm0 ym0 zm0 m  // image m, landmark 0                                 *
   *   xm1 ym1 zm1 m  // image m, landmark 1                                 *
   *       ...                                                               *
   *   xmn ymn zmn m  // image m, landmark n                                 *
   *                                                                         *
   ***************************************************************************/

  /** Sanity checks. */
  MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();
  if( !movingPointSet )
  {
    itkExceptionMacro( << "Moving point set has not been assigned" );
  }

  /** Initialize some variables */
  this->m_NumberOfPointsCounted = 0;
  MeasureType measure = NumericTraits< MeasureType >::Zero;
  derivative = DerivativeType( this->GetNumberOfParameters() );
  derivative.Fill( NumericTraits< DerivativeValueType >::ZeroValue() );
  NonZeroJacobianIndicesType nzji( this->m_Transform->GetNumberOfNonZeroJacobianIndices() );
  TransformJacobianType jacobian;

  InputPointType  movingPoint, meanPoint;
  OutputPointType mappedPoint;

  /** Call non-thread-safe stuff, such as:
   *   this->SetTransformParameters( parameters );
   *   this->GetImageSampler()->Update();
   * Because of these calls GetValueAndDerivative itself is not thread-safe,
   * so cannot be called multiple times simultaneously.
   * This is however needed in the CombinationImageToImageMetric.
   * In that case, you need to:
   * - switch the use of this function to on, using m_UseMetricSingleThreaded = true
   * - call BeforeThreadedGetValueAndDerivative once (single-threaded) before
   *   calling GetValueAndDerivative
   * - switch the use of this function to off, using m_UseMetricSingleThreaded = false
   * - Now you can call GetValueAndDerivative multi-threaded.
   */
  this->BeforeThreadedGetValueAndDerivative( parameters );
  
  /***************************************************************************
   * Pass over the points to compute the actual metric and derivative *
   ***************************************************************************/
  PointIterator pointItMoving = movingPointSet->GetPoints()->Begin();
  PointIterator pointEnd      = movingPointSet->GetPoints()->End();
  int point_no = 0; // position of the current point in the list
  while( pointItMoving != pointEnd )
  {
    /** Index of the current landmark */
    auto const landmark_no = point_no % m_NumberOfLandmarks;

    /** Get the current corresponding point. */
    movingPoint = pointItMoving.Value();
    meanPoint = this->m_Landmarks[ landmark_no ].mean;
    meanPoint[ m_LastDimension ] = movingPoint[ this->m_LastDimension ];

    /** Transform the point. */
    mappedPoint = this->m_Transform->TransformPoint( meanPoint );

    /** Check if it is inside the mask */
    if( this->m_MovingImageMask.IsNotNull() && ! this->m_MovingImageMask->IsInside( mappedPoint ) )
    {
      ++pointItMoving;
      ++point_no;
      continue;
    }

    this->m_NumberOfPointsCounted++;

    /** Get the TransformJacobian dT/dmu. */
    this->m_Transform->GetJacobian( meanPoint, jacobian, nzji );

    /** Update the mean
     * The landmarks are sorted such that `point_no % m_NumberOfLandmarks`
     * is equal to the index of the current landmark.  */
    VnlVectorType diffPoint = ( movingPoint - mappedPoint ).GetVnlVector();
    diffPoint[ m_LastDimension ] = 0.0;
    MeasureType   distance  = diffPoint.magnitude();
    measure += distance;

    /** Calculate the contributions to the derivatives with respect to each parameter. */
    if( distance > vcl_numeric_limits< MeasureType >::epsilon() )
    {
      VnlVectorType diff_2 = diffPoint / distance;
      if( nzji.size() == this->GetNumberOfParameters() )
      {
        /** Loop over all Jacobians. */
        derivative -= diff_2 * jacobian;
      }
      else
      {
        /** Only pick the nonzero Jacobians. */
        for( unsigned int i = 0; i < nzji.size(); ++i )
        {
          const unsigned int index  = nzji[ i ];
          VnlVectorType      column = jacobian.get_column( i );
          derivative[ index ] -= dot_product( diff_2, column );
        }
      }
    } // end if distance != 0

    ++pointItMoving;
    ++point_no;

  } // end loop over all moving points

  /** Copy the measure to value. */
  value = measure;
  if( this->m_NumberOfPointsCounted > 0 )
  {
    derivative /= this->m_NumberOfPointsCounted;
    value       = measure / this->m_NumberOfPointsCounted;
  }

} // end GetValueAndDerivative()


} // end namespace itk

#endif // end #ifndef __itkCorrespondingPointsMeanDistancePointMetric_hxx
