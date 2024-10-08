// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once

#include "PatchWriter.h"
#include "config.h"

#include "tarch/logging/Log.h"

#include <vector>

#ifdef UseHDF5
#include "hdf5.h"
#endif

namespace tarch {
  namespace plotter {
    namespace griddata {
      namespace blockstructured {
        class PeanoHDF5PatchFileWriter;
      }
    }
  }
}

/**
 * HDF 5 writer
 *
 * This class works if and only if you have compiled Peano with --with-hdf5 (autotools).
 */
class tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter: public tarch::plotter::griddata::blockstructured::PatchWriter {
  protected:
    static tarch::logging::Log _log;
    static const std::string HEADER;

    const int  _dimensions;
    const bool _compress;

    int _vertexCounter;
    int _cellCounter;

    #ifdef UseHDF5
    hid_t       _file;

    hid_t  createDataTableProperties(int lineWidth, int rowCount) const;
    #endif

    bool        _isOpen;

    std::vector<double>  _geometryData;

    /**
     * See the cookbook. At any time, the writer pipes data only into one
     * dataset (subdirectory) which is identified through a unique number.
     * Yet, HDF5 works with identifiers (string) instead of numbers, so you
     * have to convert it through getNameOfCurrentDataset().
     */
    int         _numberOfActiveDataset;

    /**
     * @see _numberOfActiveDataset
     */
    std::string  getNameOfCurrentDataset() const;

  public:
    class CellDataWriter: public tarch::plotter::griddata::blockstructured::PatchWriter::CellDataWriter {
      protected:
        tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter& _writer;

        const std::string    _identifier;
        const int            _numberOfUnknowns;
        std::vector<double>  _data;

      public:
        CellDataWriter(
          const std::string& identifier,
          int                unknownsPerAxis,
          int                numberOfUnknowns, const std::string& description,
          const std::string& metaData,
          double*            mapping,
          tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter& writer
        );

        virtual ~CellDataWriter();

         /**
          * Write data for one cell.
          *
          * @param index Index of the cell. This index has to equal the index
          *              used for the cell within the VTKWriter class
          *              interface.
          * @param value Value for the cell.
          */
        virtual void plotCell( int index, double value ) override;
        virtual void plotCell( int index, double* values ) override;

        virtual void close() override;

        virtual void assignRemainingCellsDefaultValues() override;

        virtual int getFirstCellWithinPatch(int index) const override;
     };

     class VertexDataWriter: public tarch::plotter::griddata::blockstructured::PatchWriter::VertexDataWriter {
       protected:
         tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter& _writer;

         const std::string    _identifier;
         const int            _numberOfUnknowns;
         std::vector<double>  _data;

       public:
         VertexDataWriter(
           const std::string& identifier,
           int                unknownsPerAxis,
           int                numberOfUnknowns, const std::string& description,
           const std::string& metaData,
           double*            mapping,
           tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter& writer
         );

         ~VertexDataWriter();

         /**
          * Write data for one cell.
          *
          * @param index Index of the vertex. This index has to equal the index
          *              used for the cell within the VTKWriter class
          *              interface.
          * @param value Value for the cell.
          */
         virtual void plotVertex( int index, double value ) override;
         virtual void plotVertex( int index, double* values ) override;

         virtual void close() override;

         /**
          * @see close()
          */
         virtual void assignRemainingVerticesDefaultValues() override;

         virtual int getFirstVertexWithinPatch(int index) const override;
     };

    PeanoHDF5PatchFileWriter(
      int                  dimension,
      const std::string&   filename,
      bool                 append,
      bool                 compress
    );

    virtual ~PeanoHDF5PatchFileWriter();

    /**
     * Caller has to destroy this instance manually.
     */
    virtual CellDataWriter*    createCellDataWriter( const std::string& identifier, int unknownsPerAxis, int recordsPerCell, const std::string& description ) override;
    virtual CellDataWriter*    createCellDataWriter( const std::string& identifier, int unknownsPerAxis, int recordsPerCell, const std::string& description, const std::string& metaData ) override;

    /**
     * The mapping is an additional field that has d * (n+1)^d doubles that
     * describe how the vertices within a unit cube are distributed. d is the
     * dimension of the plotter, n is the number of cells per axis.
     */
    virtual CellDataWriter*    createCellDataWriter( const std::string& identifier, int unknownsPerAxis, int recordsPerCell, const std::string& description, const std::string& metaData, double* mapping ) override;

    /**
     * Caller has to destroy this instance manually.
     */
    virtual VertexDataWriter*  createVertexDataWriter( const std::string& identifier, int unknownsPerAxis, int recordsPerVertex, const std::string& description ) override;
    virtual VertexDataWriter*  createVertexDataWriter( const std::string& identifier, int unknownsPerAxis, int recordsPerVertex, const std::string& description, const std::string& metaData  ) override;
    virtual VertexDataWriter*  createVertexDataWriter( const std::string& identifier, int unknownsPerAxis, int recordsPerVertex, const std::string& description, const std::string& metaData, double* mapping ) override;

    virtual int plotPatch(
      const tarch::la::Vector<2,double>& offset,
      const tarch::la::Vector<2,double>& size
    ) override;

    virtual int plotPatch(
      const tarch::la::Vector<3,double>& offset,
      const tarch::la::Vector<3,double>& size
    ) override;

    /**
     * @return Write has been successful
     */
    virtual bool writeToFile() override;

    /**
     * @return Whether writer is ready to accept data.
     */
    virtual bool isOpen() override;

    /**
     * Clear the writer, i.e. erase all the data. However, as the writer does
     * not track how many vertex and cell writers you've created, it's up to
     * you to ensure that none of these instances is left.
     */
    virtual void clear() override;

    void addMetaData(const std::string& metaData);
};
