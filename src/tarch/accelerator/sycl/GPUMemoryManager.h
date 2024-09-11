// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#include "tarch/multicore/multicore.h"


#if !defined(_TARCH_ACCELERATOR_SYCL_GPU_MEMORY_MANAGER) and defined(GPUOffloadingSYCL)
#define _TARCH_ACCELERATOR_SYCL_GPU_MEMORY_MANAGER


#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/logging/Log.h"

#include <map>
#include <list>
#include <vector>

#include "tarch/accelerator/accelerator.h"


namespace tarch {
  namespace accelerator {
    namespace sycl {
      class GPUMemoryManager;
    }
  }
}

/**
 * Singleton that manages persistent memory on the GPUs
 *
 * The idea is that this is basically a heap rep for the GPU and that code
 * requiring GPU memory does not allocate this memory but instead grabs it
 * from here. As it is not freed afterwards, we avoid overheads that arise
 * from frequent re-allocation.
 *
 * It seems that memory allocation on GPUs is very slow with most OMP runtimes.
 * Therefore, it is reasonable to manage the GPU memory in the user space.
 */
class tarch::accelerator::sycl::GPUMemoryManager {
  private:
    static tarch::logging::Log          _log;

    tarch::multicore::BooleanSemaphore  _semaphore;

    /**
     * Map of pages (memory blocks) per size.
     *
     */
    typedef std::map< int, std::vector<void*> >  FreePages;

    /**
     * Counterpart of FreePages, i.e. map of memory blocks (addresses) onto
     * sizes.
     */
    typedef std::map< void*, int >               UsedPages;

    struct DeviceMemory {
      FreePages freePages;
      UsedPages usedPages;
    };

    /**
     * Map from device number onto memory
     *
     * I don't want to use the copy here as key, so I instead use a pointer.
     */
    std::map<::sycl::queue*, DeviceMemory>  _deviceMemory;

    GPUMemoryManager() = default;

  public:
    GPUMemoryManager( const GPUMemoryManager& copy ) = delete;

    /**
     * Invoke clear().
     */
    virtual ~GPUMemoryManager();

    static GPUMemoryManager& getInstance();

    /**
     * @param size Count of block, i.e. number of doubles
     *
     * @return This pointer is a device pointer, i.e. do not use it on the host.
     */
    template <class T>
    T* allocate( int size, ::sycl::queue& device ) {
      return reinterpret_cast<T*>( allocateRawData(size * static_cast<int>(sizeof(T)), device ) );
    }

    /**
     * Return a pointer to new device data of size size bytes on device.
     *
     * This operation is thread-safe. If tries to recycle data that has been
     * allocated before. If this is not successful, it acquires new data on
     * device and then hands out this new data.
     *
     * Before we return the address of the new data, we memorise it as booked
     * data and store the size corresponding to this device. This is important:
     * When we return the data through free(), we have to look up what the
     * underlying size is and bucket sort the pointer for reusage
     * appropriately.
     *
     * The routine could, in theory, hand out a pointer to a memory block which
     * points to a chunk of memory that is bigger than size. In the context of
     * Peano this would however be stupid: Most Peano codes work with rather
     * regular chunks of data and if we use big chunks to host smaller data
     * footprint, we actually run risk to fragment the device memory.
     *
     * To make the allocation efficient, i.e. to avoid too many allocations,
     * we do not acquire one more entry per request, but we immediately
     * acquire several entries. The number of entries that we grab in one
     * rush is given by growth. Many codes acquire data that depends on the
     * Dimensions. Therefore, I use this magic default here.
     *
     * @param growth Number greater than zero which tells the code how many entries
     * (pages) to book whenever it runs out of memory.
     */
    void* allocateRawData( int size, ::sycl::queue&  device );

    /**
     * Free does not actually free data, but instead hands the memory block
     * back to the singleton, such that the next allocateRawData() call can
     * reuse the memory region.
     */
    void free(void* devicePointer, ::sycl::queue&  device);

    /**
     * Idempotent function
     *
     * Run over _deviceMemory, i.e. handle device by device. If a device still
     * has used pages, the user code has not properly freed data. We issue an
     * error.
     *
     * This operation is thread-safe.
     */
    void clear();
};


#endif
