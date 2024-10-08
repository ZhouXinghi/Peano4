// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once

#include "TraversalObserver.h"

#include "tarch/logging/Log.h"
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/mpi/BooleanSemaphore.h"

#include "tarch/plotter/griddata/unstructured/vtk/VTUTextFileWriter.h"

#include "config.h"

namespace peano4 {
  namespace grid {
    class TraversalVTKPlotter;
  }
}

/**
 * Observer which pipes the automaton transitions into a VTK file
 *
 * While we use the up-to-date vtk format, the observer plots the whole thing
 * as a discontinuous unstructured mesh. It is not particular sophisticated.
 *
 * The plotter can write whole time series. For this, you have to invoke
 * startNewSnapshot() prior to each plot. It is the latter which also ensures
 * that parallel plots in an MPI environment do work.
 *
 * <h2> Parallel plotting </h2>
 *
 * Each tree dumps its own vtk file. That is, each thread and each rank in
 * theory might write its file parallel to the other guys. VTK/VTU offers us
 * to define a metafile (pvtu) which collocates various dumps. As we create
 * one observer per thread through clone(), every thread on every rank
 * has its instance and pipes its data. getFilename() ensures that no file
 * is overwritten. It combines the tree number with a counter, and _counter,
 * which is static, is incremented through endTraversalOnRank() which I expect
 * the user to call once after each traversal.
 *
 * <h2> Known bugs </h2>
 *
 * As the MPI domain decomposition creates fake observers for the master of
 * a local rank when it is created, we'll have multiple entries for forking
 * ranks in the meta file.
 */
class peano4::grid::TraversalVTKPlotter: public peano4::grid::TraversalObserver {
  protected:
    static tarch::logging::Log                 _log;

    const std::string                                                                _filename;
    const int                                                                        _spacetreeId;

    /**
     * Does the actual plotting, i.e. all checks/decision making is already done before
     */
    void plotCell(
        const GridTraversalEvent&  event
    );

  private:
    tarch::plotter::griddata::unstructured::vtk::VTUTextFileWriter*                  _writer;
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::VertexWriter*    _vertexWriter;
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellWriter*      _cellWriter;
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellDataWriter*  _spacetreeIdWriter;
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellDataWriter*  _coreWriter;

    static tarch::mpi::BooleanSemaphore                                              _sempahore;

  public:
    /**
     * You have to invoke startNewSnapshot() if you wanna have a pvd file
     * immediately after you've created this observer in the main code.
     *
     * If this guy is ran on the global master,
     */
    TraversalVTKPlotter( const std::string& filename, int treeId=-1 );
    virtual ~TraversalVTKPlotter();

    virtual void beginTraversal(
      const tarch::la::Vector<Dimensions,double>&  x,
      const tarch::la::Vector<Dimensions,double>&  h
    ) override;

    virtual void endTraversal(
      const tarch::la::Vector<Dimensions,double>&  x,
      const tarch::la::Vector<Dimensions,double>&  h
    ) override;

    virtual void loadCell(
      const GridTraversalEvent&  event
    ) override;

    virtual void storeCell(
      const GridTraversalEvent&  event
    ) override;

    virtual void enterCell(
      const GridTraversalEvent&  event
    ) override;

    virtual void leaveCell(
      const GridTraversalEvent&  event
    ) override;

    virtual TraversalObserver* clone(int spacetreeId) override;

    /**
     * Obviously empty for this particular observer.
     */
    virtual std::vector< GridControlEvent > getGridControlEvents() const override;
};
