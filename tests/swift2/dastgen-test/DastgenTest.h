#pragma once

#include "globaldata/dummyPart.h"

/**
 * @file DastgenTest.h
 * @brief Contains routines intended to test that the classes and attributes
 * DaStGen2 generates compile and run properly, as used in the `dastgen-test`
 * example and generated by DastgenTestDummyParticle.py.
 */


namespace swift2 {
  namespace dastgenTest {

    /**
     * Main testing function.
     * Should call getters and setters for all available variables,
     * where applicable. Also try constructor versions of particles.
     **/
    template <typename ParticleContainer>
    void checkDastgen(ParticleContainer& assignedParticles);


    /**
     * Perform all the checks on a particle's boolean members
     */
    void checkBooleans(
      tests::swift2::dastgenTest::globaldata::dummyPart& particle
    );


    /**
     * Perform all the checks on a particle's double members
     */
    void checkDoubles(
      tests::swift2::dastgenTest::globaldata::dummyPart& particle
    );


    /**
     * Perform all the checks on a particle's enum members
     */
    void checkEnums(tests::swift2::dastgenTest::globaldata::dummyPart& particle
    );


    /**
     * Perform all the checks on a particle's integer members
     */
    void checkIntegers(
      tests::swift2::dastgenTest::globaldata::dummyPart& particle
    );


    /**
     * Perform all the checks on a particle's string members
     */
    void checkStrings(
      tests::swift2::dastgenTest::globaldata::dummyPart& particle
    );


    /**
     * Perform all the checks on a particle's boolean array members
     */
    void checkBooleanArrays(
      tests::swift2::dastgenTest::globaldata::dummyPart& particle
    );


    /**
     * Perform all the checks on a particle's double array members
     */
    void checkDoubleArrays(
      tests::swift2::dastgenTest::globaldata::dummyPart& particle
    );


    /**
     * Perform all the checks on a particle's integer array members
     */
    void checkIntegerArrays(
      tests::swift2::dastgenTest::globaldata::dummyPart& particle
    );


    /**
     * Perform all the checks on a particle's double array members
     */
    void checkPeanoDoubleArrays(
      tests::swift2::dastgenTest::globaldata::dummyPart& particle
    );


    /**
     * Perform all the checks on a particle's integer array members
     */
    void checkPeanoIntegerArrays(
      tests::swift2::dastgenTest::globaldata::dummyPart& particle
    );


    /**
     * Perform all the checks on a particle's user defined types.
     */
    void checkUserDefinedType(
      tests::swift2::dastgenTest::globaldata::dummyPart& particle
    );


    /**
     * Check the constructor methods.
     */
    void checkConstructors();


    /**
     * A dummy function to move the simulation forward in time.
     * The time step size is hard-coded such that the sim runs
     * for 2 steps. */
    template <typename Particle>
    void dummyMoveForwardInTime();


    /**
     * Report time-step size and current time at the end of a step.
     **/
    template <typename Particle>
    void reportStep(const std::string& particleName);

  } // namespace dastgenTest
} // namespace swift2


#include "DastgenTest.cpph"
