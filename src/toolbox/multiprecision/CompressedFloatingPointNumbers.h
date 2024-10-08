// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once


namespace toolbox {
  /**
   * The heap namespace holds all routines that we use in connection with
   * Peano's heap data storage. This means
   *
   * - the actual heap classes
   * - allocators used for heaps and architectures requiring concurrent and aligned layouts
   * - data conversion/compression routines
   */
  namespace multiprecision {
    /**
     * Takes a double and returns the exponent and the mantissa.
     *
     * Is just a wrapper around C++'s std::ldexp.
     */
    template <class T>
    void decompose(
      const double& value,
      char&         exponent,
      T&            mantissa
    );

    /**
     * Decompose floating point value
     *
     * Decomposes the input value into eight different variants with
     * 1,2,3,... bytes (chars) per mantissa. It also gives the errors.
     *
     * You may want to cast one of the integers into a sequence of chars
     * or other data types. While the array mantissa comprises long ints,
     * only the first kth bytes are used in the k-1th entry.
     *
     * Negative mantissa pose a challenge to the byte compression. Due to the
     * integer codes, we may not just hold negative values in the long ints.
     * -1 has all bits set, e.g. So if someone later only stores the minor k
     * bytes of mantissa, this will fail. So we restrict to unsigned values
     * first as we extract the sign and then multiply the biggest byte with
     * the sign in the end. This implies that the biggest byte of each
     * mantissa entry actually holds only seven digits.
     *
     * Please note that there are two compose() operations, i.e. two
     * counterparts of decompose(). Ensure you use the right one.
     *
     * @return exponent Sequence of exponents belonging to the eight decompositions.
     * @return mantissa Sequence of mantissas belonging to the eight decompositions.
     * @return error    Errors associated to the decompositions. They are all positive values.
     *
     * <h2> Bug/interoperability </h2>
     *
     * This operation works if and only if your system maps a long int onto a 64 bit value.
     *
     * We noticed that this operation yields invalid results with GCC if we
     * switch on optimisation. It works fine with ICC and any optimisation
     * level. We thus make the code translate with -O0 locally if ICC is not
     * used.
     *
     * <h2> Continue to work with results </h2>
     *
     * See decompose() how to continue to work with one of the result entries.
     *
     * @see compose( char, long int, int )
     *
     */
    void decomposeIntoEightVariants(
      double        value,
      char          exponent[8],
      long int      mantissa[8],
      double        error[8]
    );

    /**
     * @see decomposeIntoEightVariants
     */
    void decomposeIntoFourVariants(
      double        value,
      char          exponent[4],
      int           mantissa[4],
      double        error[4]
    );

    /**
     * Analyses the handed data and determines the most aggressive compression.
     *
     * The operation returns the number of bytes used in the mantissa. If you
     * use this index and reduce it by one, you find the right entry in
     * decomposeIntoFourVariants().
     *
     * @return Value between 1 and 7 or 0 if data does not hold any information
     */
    int findMostAgressiveCompression(
      double        value,
	  double        maxAbsoluteError
    );


    /**
     * Wrapper around findMostAgressiveCompression() that works for an array of
     * count entries.
     *
     * I run over the array of length count and check for every single entry
     * whether we could store it with fewer than 8 bytes without harming the
     * maxError constraint. This check relies on
     * findMostAgressiveCompression(double,double,bool). The routine then
     * returns the minimal number of bytes you have to invest to encode this
     * whole array such that the result still remains in the error bounds.
     *
     * <h2> Usage pattern </h2>
     *
     * This is the standard workflow:
     *
     * - Convert an array into its hierarchical representation. Usually, I
     *   determine the mean value and then store all remaining values within the
     *   array as actual value minus mean.
     * - Call this routine on the remaining array.
     * - If the result is bigger or equal 8, I store/send away the original
     *   array. Otherwise:
     * - I store the mean as 8 bytes into a char stream. I then loop over the
     *   array and per array entry call decompose. I know how many bytes of
     *   exponent and mantissa are valid from the result of this routine
     *   and thus enqueue only those guys in an output stream/MPI message.
     */
    int findMostAgressiveCompression(
      const double  values[],
      int           count,
	  double        maxError
    );


    /**
     * If you wanna the decompose values into compressed char series, use
     *
     * <pre>
     char* pMantissa = reinterpret_cast<char*>( &(mantissa) );
     for (int j=0; j<requiredBytesPerMantissa; j++) {
       do something with mantissa[j]
     }
       </pre>
     *
     * Don't forget to invoke shrink_to_fit() on your std::vector containers if
     * you use them. Otherwise, C++ most likely won't free your memory.
     */
    void decompose(
      double        value,
      char&         exponent,
      long int&     mantissa,
      int           bytesUsedForMantissa
    );

    /**
     * Counterpart of decompose( const double&, char& , T& )
     *
     * As it is the counterpart, code typically looks similar to
     * <pre>
  char exponent  = 0;
  long int mantissa = 0;
  char* pMantissa = reinterpret_cast char* ( &(mantissa) );

  for (int j=bytesForMantissa-1; j>=0; j--) {
    pMantissa[j] = CompressedDataHeap::getInstance().getData( heapIndex )[compressedDataHeapIndex]._persistentRecords._u;
    compressedDataHeapIndex--;
  }
  exponent = CompressedDataHeap::getInstance().getData( heapIndex )[compressedDataHeapIndex]._persistentRecords._u;
</pre>
     *
     * There are two evergreens of pitfalls:
     *
     * - If you use push and pop to add stuff to the char heaps when you
     *   decompose, you have to invert all the accesses: you push the
     *   exponent first but you pop it last, and you pop the mantissas
     *   bytes in inverted order. I personally prefer to use a resize()
     *   command on the vector to avoid re-allocation and to access the
     *   vector entries directly (see snippet above). In this case, I
     *   could do without an inversion of access order.
     * - Please ensure that you set the long int mantissa above to 0!
     *   As you overwrite only parts of it (with compressed data), you
     *   might otherwise end up with garbage in the higher bytes.
     */
    template <class T>
    double compose(
      const char&  exponent,
      const T&     mantissa
    );

    /**
     * Counterpart of decompose( double, char[8], long int[8], double[8] ).
     *
     * Please note that bytesUsed is basically which entry of the decompose
     * return arrays you used. However, it is increased by one, i.e. if you
     * pass 1, the operation assumes that you've taken the entry exponent[0]
     * of decompose, e.g.
     */
    double compose(
      char         exponent,
      long int     mantissa,
      int          bytesUsedForMantissa
    );
  }
}


#include "CompressedFloatingPointNumbers.cpph"

