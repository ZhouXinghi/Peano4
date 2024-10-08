// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org

#include "Measurement.h"

#include <cmath>
#include <limits>
#include <sstream>

#include "tarch/Assertions.h"
#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"


tarch::logging::Log tarch::timing::Measurement::_log("tarch::timing::Measurement");


tarch::timing::Measurement::Measurement(const double& accuracy):
  _accuracy(accuracy),
  _accumulatedValue(0.0),
  _lastValue(0.0),
  _accumulatedSquares(0.0),
  _numberOfMeasurements(0.0),
  _isAccurateValue(false),
  _min(std::numeric_limits<double>::max()),
  _max(-std::numeric_limits<double>::max()),
  _minMeasurement(-1),
  _maxMeasurement(-1) {

  assertion1(accuracy >= 0.0, accuracy);
}


void tarch::timing::Measurement::erase() {
  _accumulatedValue     = 0.0;
  _accumulatedSquares   = 0.0;
  _numberOfMeasurements = 0.0;
  _isAccurateValue      = false;
  _min                  = std::numeric_limits<double>::max();
  _max                  = -std::numeric_limits<double>::max();
  _minMeasurement       = -1;
  _maxMeasurement       = -1;
}


double tarch::timing::Measurement::getMeanValue() const {
  assertion(_numberOfMeasurements > 0.0);
  return _accumulatedValue / _numberOfMeasurements;
}


double tarch::timing::Measurement::getMeanValueOfNextStep(double newValue) const {
  assertion(_numberOfMeasurements > 0.0);
  return (_accumulatedValue + newValue) / (_numberOfMeasurements + 1.0);
}


double tarch::timing::Measurement::getAccumulatedValue() const { return _accumulatedValue; }


double tarch::timing::Measurement::getMostRecentValue() const {
  return _lastValue;
}


double tarch::timing::Measurement::getValue() const {
  if (tarch::la::smaller(_numberOfMeasurements, 1.0)) {
    assertionEquals(_accumulatedValue, 0.0);
    return 0.0;
  } else {
    return getMeanValue();
  }
}


double tarch::timing::Measurement::getStandardDeviation() const {
  if (tarch::la::smaller(_numberOfMeasurements, 1.0)) {
    assertionEquals(_accumulatedValue, 0.0);
    return 0.0;
  } else {
    return std::sqrt(std::abs(_accumulatedSquares / _numberOfMeasurements - getMeanValue() * getMeanValue()));
  }
}


bool tarch::timing::Measurement::isAccurateValue() const { return _isAccurateValue; }


void tarch::timing::Measurement::setAccuracy(const double& value) {
  assertion(value >= 0.0);
  _accuracy = value;
}


void tarch::timing::Measurement::increaseAccuracy(const double& factor) {
  assertion(factor > 0.0);
  _accuracy /= factor;
}


void tarch::timing::Measurement::setValue(const double& value) {
  if (tarch::la::smallerEquals(_numberOfMeasurements, 1.0)) {
    assertion(!_isAccurateValue);
  } else {
    const double differenceDueToNewValue = tarch::la::smallerEquals(getMeanValue(), 1.0)
                                             ? std::abs(getMeanValueOfNextStep(value) - getMeanValue())
                                             : std::abs(getMeanValueOfNextStep(value) / getMeanValue() - 1.0);
    _isAccurateValue                     = _numberOfMeasurements > 0 && differenceDueToNewValue < _accuracy;
  }

  if (_min > value) {
    _min            = value;
    _minMeasurement = _numberOfMeasurements;
  }
  if (_max < value) {
    _max            = value;
    _maxMeasurement = _numberOfMeasurements;
  }

  _accumulatedValue     += value;
  _accumulatedSquares   += value * value;
  _numberOfMeasurements += 1.0;
  _lastValue             = value;

  logDebug(
    "setValue",
    "set "
      << static_cast<int>(_numberOfMeasurements) << "th value = " << value << ", current deviation="
      << ((_accumulatedValue + value) / (_numberOfMeasurements + 1) - _accumulatedValue / _numberOfMeasurements)
      << ", accumulatedValue=" << _accumulatedValue << ", _numberOfMeasurements" << _numberOfMeasurements
      << ", isAccurateValue=" << _isAccurateValue << ", _accuracy=" << _accuracy
  );
}


int tarch::timing::Measurement::getNumberOfMeasurements() const { return static_cast<int>(_numberOfMeasurements); }


std::string tarch::timing::Measurement::toString() const {
  std::ostringstream msg;
  msg << "(avg=";
  if (_numberOfMeasurements > 0) {
    msg << getValue();
  } else {
    msg << "nan";
  }
  msg
    << ",#measurements=" << _numberOfMeasurements << ",max=" << _max << "(value #"
    << (static_cast<int>(_maxMeasurement)) << ")"
    << ",min=" << _min << "(value #" << (static_cast<int>(_minMeasurement)) << ")";

  if (_numberOfMeasurements > 0) {
    msg
      << ",+" << (_max - getValue()) / getValue() * 100.0 << "%"
      << ",-" << (getValue() - _min) / getValue() * 100.0 << "%"
      << ",std-deviation=" << getStandardDeviation();
  }

  if (_accuracy > 0.0) {
    msg << ",eps=" << _accuracy;
    if (!_isAccurateValue) {
      msg << " [not a valid averaged value yet]";
    } else {
      msg << " [converged]";
    }
  }

  msg << ")";

  return msg.str();
}


double tarch::timing::Measurement::max() const {
  assertion1(isAccurateValue(), toString());
  return _max;
}


double tarch::timing::Measurement::min() const {
  assertion1(isAccurateValue(), toString());
  return _min;
}


double tarch::timing::Measurement::getAccuracy() const { return _accuracy; }
