// **********************************************************************************************
// ***                                     !!!WARNING!!!                                      ***
// *** WARNING: AUTO GENERATED FILE! DO NOT MODIFY BY HAND! YOUR CHANGES WILL BE OVERWRITTEN! ***
// ***                                     !!!WARNING!!!                                      ***
// ***                  Generated by Peano's Python API: www.peano-framework.org              ***
// **********************************************************************************************
#pragma once

#include "peano4/datamanagement/CellMarker.h"
#include "tarch/multicore/Tasks.h"

#include "exahype2/EnclaveBookkeeping.h"
#include "exahype2/EnclaveTask.h"
#include "repositories/SolverRepository.h"

{{SOLVER_INCLUDES}}

#include <vector>
#include <memory>


#ifdef UseSmartMPI
#include "smartmpi.h"
#endif

{% for item in NAMESPACE -%}
  namespace {{ item }} {

{%- endfor %}
  class {{CLASSNAME}};

{% for item in NAMESPACE -%}
  }
{%- endfor %}

/**
 * Single task that can also take multiple tasks and deploy them to the GPU
 *
 * @author Tobias Weinzierl
 */
class {{NAMESPACE | join("::")}}::{{CLASSNAME}}: public ::exahype2::EnclaveTask
#ifdef UseSmartMPI
, public smartmpi::Task
#endif
{
  private:
    static tarch::logging::Log  _log;

    #ifdef UseSmartMPI
    /*
     * This uniquely identifies an enclave task. Unlike _enclaveTaskTypeId it
     * differentiates between tasks within a task type, not between task types.
     *
     * Crucially, the value of this ID is not necessarily the same as that of
     * tarch::multicore::Task::_id . That is because _id is reset each time
     * we move a task to another rank and reconstitute it there (see the constructor
     * of tarch::multicore::Task). Instead, _remoteTaskId tracks the ID of
     * the original task object (i.e. the task originally spawned and
     * only then moved). In smartmpi we need to keep track of the task's
     * ID so that it can be bookmarked correctly after being moved around.
     *
     * As such moveTask(...), sends _id if the task has not yet been
     * moved and _remoteTaskId if it has already been moved. Similarly,
     * _remoteTaskId is always sent when forwarding task outcomes since the
     * task will already have been moved.
     */
    int          _remoteTaskId = -1;
    #endif

    /**
     * This is a unique static id which we use both for task fusion and within
     * SmartMPI to distinguish different task types.
     */
    static int                                _enclaveTaskTypeId;

  public:
    static int getEnclaveTaskTypeId();

    /**
     * Simple call to compute kernels.
     *
     * The finite volume solver assumes that _inputValues are reconstructed
     * values, i.e., a temporary data field. So we don't have to allocate
     * anything but we have to free stuff after the run. This is done here in
     * the Finite Volume solver. We assume that the input data has been
     * allocated through tarch::allocateMemory() with the flag
     * ManagedSharedAcceleratorDeviceMemory.
     *
     * @return Max eigenvalue over cell
     */
    static double applyKernelToCell(
      const ::peano4::datamanagement::CellMarker& marker,
      double                                      t,
      double                                      dt,
      double* __restrict__                        reconstructedPatch,
      double* __restrict__                        targetPatch
    );

    /**
     * Task constructor
     *
     * The task construction injects a call to applyKernelToCell() into
     * the generic enclave superclass. We assume that reconstructedPatch
     * is created on the heap before - typically in the actual action
     * set. As applyKernelToCell() erases the
     * reconstructedPatch argument, we don't have to call a free() in the
     * task constructor of the functor that we pass to the enclave
     * task (superclass). Likewise, the superclass EnclaveTask tasks
     * care of the allocation of the output data if no output (nullptr)
     * is specified.
     */
    {{CLASSNAME}}(
      const ::peano4::datamanagement::CellMarker& marker,
      double                                      t,
      double                                      dt,
      double* __restrict__                        reconstructedPatch,
      double* __restrict__                        output
    );

    #ifdef UseSmartMPI
    bool isSmartMPITask() const;

    virtual void runLocally() override;
    virtual void moveTask(int rank, int tag, MPI_Comm communicator) override;
    virtual void runLocallyAndSendTaskOutputToRank(int rank, int tag, MPI_Comm communicator) override;
    virtual void forwardTaskOutputToRank(int rank, int tag, MPI_Comm communicator) override;

    static smartmpi::Task* receiveTask(int rank, int tag, MPI_Comm communicator);
    static smartmpi::Task* receiveOutcome(int rank, int tag, MPI_Comm communicator, const bool intentionToForward);
    #endif

    {% if STATELESS_PDE_TERMS %}
    virtual bool fuse(const std::list<Task*>& otherTasks, int targetDevice = Host) override;
    virtual bool canFuse() const override;
    {% endif %}
};

{# Empty line here #}