// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once


#include <chrono>
#include <ctime>


#include "tarch/logging/Log.h"


namespace tarch::timing {
  class Watch;
} // namespace tarch::timing


/**
 * A simple class that has to be included to measure the clock ticks required
 * for an operation. To use it you've to create an instance of this class at
 * the beginning of the operation block you want to measure the time spent
 * within it:
 *
 * <pre>
 *   void anyOperation() {
 *     tarch::timing::Watch watch("MyClass", "anyOperation()");
 *     ...
 *   }
 * </pre>
 *
 * The result the of the measurement is written to log.info level. Note that
 * this operation works within operation blocks (e.g., a for-loop), too.
 *
 * For an extended measurement of calender (i.e., user) time and a saving of the
 * data and access via getters, we implemented parts of what Markus Langlotz
 * gave us for the ancient stokes project. Additionally, the counter overflow
 * and the non-access to the measured time were reasons for introducing this
 * extended version. Furthermore, we now may choose specific computation
 * intervals to be counted in between the start and the stop command.
 *
 * @version $Revision: 1.9 $
 * @author  Tobias Weinzierl, Tobias Neckel
 */
class tarch::timing::Watch {
private:
  /**
   * Log device the result is written to.
   */
  tarch::logging::Log _log;

  /**
   * Flag to distinguish the (original) standard watch. Is used in the
   * destructor.
   */
  bool _plotResultInDestructor;

  /**
   * Stores the name of the operation the watch is used within.
   */
  std::string _operationName;

  /**
   * Holds the clock ticks at the beginning of the time measurement.
   */
  std::clock_t _startClockTicks;

  /**
   * Holds the time at the beginning of the time measurement.
   */
  std::chrono::steady_clock::time_point _startTime; // steady_clock replaced
                                                    // the
                                                    // high_resolution_clock
                                                    // as the latter can be
                                                    // jumpy.

  /**
   * Holds the elapsed processor time.
   */
  std::clock_t _elapsedClockTicks;

  /**
   * Holds the elapsed calendar time.
   */
  double _elapsedTime;

  /**
   * Has stopTimer() been called before.
   */
  bool _isRunning;

public:
  /**
   * Construct a watch
   *
   * @param className     Name of the class the watch is used within (for log).
   * @param operationName Name of the operation the watch is used within.
   * @param extendedWatchType Bool for extended.
   */
  Watch(
    const std::string& className,
    const std::string& operationName,
    bool               plotResultInDestructor,
    bool               startToTickImmediately = true
  );

  /**
   * For standard version (): Stops the watch and plots the time spent. For
   * the extended version, this is the usual destructor without any special
   * actions.
   */
  virtual ~Watch();

  /**
   * (Re)Start the Timer
   *
   * This method starts the timer. Actually, it restarts the time as the
   * constructor also invokes this operation. You don't have to call start()
   * after stop(). If you omit another call to start(), you just won't reset
   * the baseline, i.e., all getXXXTime() calls will still refer to the same
   * time zero when the Watch has been constructed.
   *
   * @see stopTimer()
   */
  void start();


  /**
   * Stop timer
   *
   * This method stops the timer. You may stop a timer multiple time because
   * the stop does not reset the start time. Please note that all the getters
   * return meaningful results if and only if you call stopTimer().
   */
  void stop();


  /**
   * Return CPU Time in Seconds
   *
   * This method returns the elapsed CPU time between the start and stop
   * command of the timer, i.e., the clock ticks actually spent by the process.
   * Take care: This counter might overflow, especially on a 32 bit
   * architecture.
   *
   * !!! Multithreading
   *
   * If you use multithreading, CPU time sums up the clock ticks invested by
   * the individual cores, i.e., if you switch from a one core to a dual core
   * machine, the CPU time should roughly double.
   *
   * @return CPU time in seconds
   */
  double getCPUTime();

  /**
   * Equals getCPUTime() but returns the clock ticks instead of the time in
   * seconds.
   */
  std::clock_t getCPUTicks();

  /**
   * This method returns the elapsed calendar time between the start and stop
   * command of the timer, i.e., the real world time. The result is given in
   * seconds.
   */
  double getCalendarTime();

  /**
   * This operation returns whether the watch is currently on. If the routine
   * returns false, it is currently stopped. The concept of on/off is a
   * logical concept. Internally, a watch simply administers system time
   * snapshots, while the get routines return differences of these snapshots.
   */
  bool isOn() const;
};
