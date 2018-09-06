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
#ifndef __itkCorrespondingPointsMeanDistancePointMetric_h
#define __itkCorrespondingPointsMeanDistancePointMetric_h

#include "itkSingleValuedPointSetToPointSetMetric.h"
#include "itkPoint.h"
#include "itkPointSet.h"
#include "itkImage.h"

#include <set>
#include <vector>

namespace itk
{

/** \class CorrespondingPointsMeanDistancePointMetric
 * \brief Computes the Euclidean distance between corresponding points
 * over the last dimension in a moving point-set.
 *
 * The 4D points must be specified in the usual format, with the additional
 * constraints that they have to be sorted by dimension, and the order
 * of points in each dimension must be consistent.
 *
 * The fixed point set is ignored.
 *
 * Correspondence is needed.
 *
 *
 * \ingroup RegistrationMetrics
 */

template< class TFixedPointSet, class TMovingPointSet >
class CorrespondingPointsMeanDistancePointMetric :
  public SingleValuedPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
{
public:

  /** Standard class typedefs. */
  typedef CorrespondingPointsMeanDistancePointMetric Self;
  typedef SingleValuedPointSetToPointSetMetric<
    TFixedPointSet, TMovingPointSet >               Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CorrespondingPointsMeanDistancePointMetric,
    SingleValuedPointSetToPointSetMetric );

  /** Types transferred from the base class */
  typedef typename Superclass::TransformType           TransformType;
  typedef typename Superclass::TransformPointer        TransformPointer;
  typedef typename Superclass::TransformParametersType TransformParametersType;
  typedef typename Superclass::TransformJacobianType   TransformJacobianType;

  typedef typename Superclass::MeasureType                MeasureType;
  typedef typename Superclass::DerivativeType             DerivativeType;
  typedef typename Superclass::DerivativeValueType        DerivativeValueType;
  typedef typename Superclass::FixedPointSetType          FixedPointSetType;
  typedef typename Superclass::MovingPointSetType         MovingPointSetType;
  typedef typename Superclass::FixedPointSetConstPointer  FixedPointSetConstPointer;
  typedef typename Superclass::MovingPointSetConstPointer MovingPointSetConstPointer;

  typedef typename Superclass::PointIterator     PointIterator;
  typedef typename Superclass::PointDataIterator PointDataIterator;

  typedef typename Superclass::InputPointType    InputPointType;
  typedef typename Superclass::OutputPointType   OutputPointType;
  typedef typename OutputPointType::CoordRepType CoordRepType;
  typedef vnl_vector< CoordRepType >             VnlVectorType;

  typedef typename Superclass::NonZeroJacobianIndicesType NonZeroJacobianIndicesType;

  /** Initialise */
  void Initialize( void ) throw ( ExceptionObject );

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & parameters ) const;

  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType & parameters,
    DerivativeType & Derivative ) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & parameters,
    MeasureType & Value, DerivativeType & Derivative ) const;

protected:

  typedef typename TMovingPointSet::PointType MovingPointType;
  typedef typename MovingPointType::ValueType MovingPointValueType;

  struct Landmark {
      std::vector < MovingPointType > occurrences;
      MovingPointType mean;
  };

  CorrespondingPointsMeanDistancePointMetric();
  virtual ~CorrespondingPointsMeanDistancePointMetric() {}

  unsigned int m_LastDimension;
  unsigned int m_NumberOfPointsPerImage;
  unsigned int m_NumberOfLandmarks;
  std::set< unsigned int > m_LastDimensionValues;
  std::vector< Landmark > m_Landmarks;

private:

  CorrespondingPointsMeanDistancePointMetric( const Self & ); // purposely not implemented
  void operator=( const Self & );                             // purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCorrespondingPointsMeanDistancePointMetric.hxx"
#endif

#endif
