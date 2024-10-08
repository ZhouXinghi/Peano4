#include "tarch/plotter/griddata/unstructured/vtk/VTKBinaryFileWriter.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdio.h>

tarch::logging::Log tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::_log("tarch::plotter::griddata::"
                                                                                           "unstructured::vtk::"
                                                                                           "VTKBinaryFileWriter");

const std::string tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::HEADER
  = "# vtk DataFile Version 2.0\n "
    "Generated by Peano3 output component $Revision: 1.2 $ Author: Tobias Weinzierl\n "
    "BINARY\n ";

tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::VTKBinaryFileWriter(
  const std::string&                                 fileName,
  const std::string&                                 indexFileName,
  tarch::plotter::PVDTimeSeriesWriter::IndexFileMode mode,
  double                                             timeStamp,
  const int                                          precision
):
  _writtenToFile(false),
  _precision(precision),
  _doubleOrFloat(setDoubleOrFloatString(precision)),
  _numberOfVertices(0),
  _numberOfCells(0),
  _numberOfCellEntries(0),
  _fileName(fileName) {
  if (fileName.rfind(".vtk") != std::string::npos) {
    logWarning(
      "writeToFile()",
      "filename should not end with .vtk as routine adds extension automatically. Chosen filename prefix=" << fileName
    );
  }
  if (mode != tarch::plotter::PVDTimeSeriesWriter::IndexFileMode::NoIndexFile and indexFileName.rfind(".pvd") != std::string::npos) {
    logWarning(
      "writeToFile()",
      "index filename should not end with .pvd as routine adds extension automatically. Chosen filename prefix="
        << indexFileName
    );
  }

  // VTK does not support more precise values
  const double DefaultTimeStampPrecision = 1e-5;

  switch (mode) {
  case tarch::plotter::PVDTimeSeriesWriter::IndexFileMode::CreateNew:
    tarch::plotter::PVDTimeSeriesWriter::createEmptyIndexFile(indexFileName);
    tarch::plotter::PVDTimeSeriesWriter::appendNewData(indexFileName, fileName + ".vtk", timeStamp);
    break;
  case tarch::plotter::PVDTimeSeriesWriter::IndexFileMode::AppendNewData:
    if (not std::filesystem::exists(indexFileName + ".pvd")) {
      logInfo("PeanoTextPatchFileWriter(...)", "no index file " << indexFileName << " found. Create new one");
      tarch::plotter::PVDTimeSeriesWriter::createEmptyIndexFile(indexFileName);
    } else if (tarch::la::smaller(
                 timeStamp,
                 tarch::plotter::PVDTimeSeriesWriter::getLatestTimeStepInIndexFile(indexFileName),
                 tarch::la::relativeEpsNormaledAgainstValueGreaterOne(
                   timeStamp,
                   tarch::plotter::PVDTimeSeriesWriter::getLatestTimeStepInIndexFile(indexFileName),
                   DefaultTimeStampPrecision
                 )
               )) {
      logWarning(
        "PeanoTextPatchFileWriter(...)",
        "there is an index file "
          << indexFileName << " with data for time stamp "
          << tarch::plotter::PVDTimeSeriesWriter::getLatestTimeStepInIndexFile(indexFileName)
          << ". Will be overwritten as we dump data for time " << timeStamp
      );
      tarch::plotter::PVDTimeSeriesWriter::createEmptyIndexFile(indexFileName);
    }

    tarch::plotter::PVDTimeSeriesWriter::appendNewData(indexFileName, fileName + ".vtk", timeStamp);
    break;
  case tarch::plotter::PVDTimeSeriesWriter::IndexFileMode::NoIndexFile:
    break;
  }
}

tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::~VTKBinaryFileWriter() {
  if (!_writtenToFile) {
    assertionEqualsMsg(
      _numberOfVertices,
      0,
      "Still vertices in vtk writer pipeline. Maybe you forgot to call writeToFile() on a data vtk writer?"
    );
    assertionEqualsMsg(
      _numberOfCells,
      0,
      "Still cells in vtk writer pipeline. Maybe you forgot to call writeToFile() on a data vtk writer?"
    );
    assertionEqualsMsg(
      _numberOfCellEntries,
      0,
      "Still cell entries in vtk writer pipeline. Maybe you forgot to call writeToFile() on a data vtk writer?"
    );
  }
}

void tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::clear() {
  _writtenToFile       = false;
  _numberOfVertices    = 0;
  _numberOfCells       = 0;
  _numberOfCellEntries = 0;
  _vertexDescription.clear();
  _cellDescription.clear();
  _cellTypeDescription.clear();
  _vertexDataDescription.clear();
  _cellDataDescription.clear();
}

bool tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::writeToFile() {
  assertion(!_writtenToFile);

  std::ostringstream filenameStream;
  filenameStream << _fileName << ".vtk";
  const std::string filename = filenameStream.str();

  std::ofstream out;
  out.open(filename.c_str(), std::ios::binary);
  if ((!out.fail()) && out.is_open()) {
    logDebug("close()", "opened data file " + filename);

    out << HEADER << std::endl << std::endl;

    out << "DATASET UNSTRUCTURED_GRID" << std::endl
        << "POINTS " << _numberOfVertices << " " << _doubleOrFloat << std::endl;
    out << _vertexDescription.rdbuf() << std::endl << std::endl;

    out << "CELLS " << _numberOfCells << " " << _numberOfCellEntries << std::endl;
    out << _cellDescription.rdbuf() << std::endl << std::endl;

    out << "CELL_TYPES " << _numberOfCells << std::endl;
    out << _cellTypeDescription.rdbuf() << std::endl << std::endl;

    if (_numberOfVertices > 0 && !_vertexDataDescription.str().empty()) {
      out << "POINT_DATA " << _numberOfVertices << std::endl;
      out << _vertexDataDescription.rdbuf() << std::endl << std::endl;
    }

    if (_numberOfCells > 0 && !_cellDataDescription.str().empty()) {
      out << "CELL_DATA " << _numberOfCells << std::endl;
      out << _cellDataDescription.rdbuf() << std::endl << std::endl;
    }

    logDebug("close()", "data written to " + filename);

    _writtenToFile = true;
    return true;
  } else {
    logError("close()", "unable to write output file " + filename);
    return false;
  }
}

bool tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::isOpen() { return !_writtenToFile; }

tarch::plotter::griddata::unstructured::UnstructuredGridWriter::VertexWriter* tarch::plotter::griddata::unstructured::
  vtk::VTKBinaryFileWriter::createVertexWriter() {
  return new tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::VertexWriter(*this);
}

tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellWriter* tarch::plotter::griddata::unstructured::
  vtk::VTKBinaryFileWriter::createCellWriter() {
  return new tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::CellWriter(*this);
}

void tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::validateDataWriterIdentifier(
  const std::string& identifier
) const {
  if (identifier.empty()) {
    logWarning(
      "validateDataWriterIdentifier(string)",
      "identifier for vtk file is empty. Spaces are not allowed for vtk data field identifiers and some vtk "
      "visualisers might crash."
    );
  }
  if (identifier.find(' ') != std::string::npos) {
    logWarning(
      "validateDataWriterIdentifier(string)",
      "identifier \""
        << identifier
        << "\" contains spaces. Spaces are not allowed for vtk data field identifiers and some vtk visualisers might "
           "crash."
    );
  }
}

tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellDataWriter* tarch::plotter::griddata::unstructured::
  vtk::VTKBinaryFileWriter::createCellDataWriter(const std::string& identifier, int recordsPerCell) {
  validateDataWriterIdentifier(identifier);
  return new tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::CellDataWriter(
    identifier, *this, recordsPerCell
  );
}

tarch::plotter::griddata::unstructured::UnstructuredGridWriter::VertexDataWriter* tarch::plotter::griddata::
  unstructured::vtk::VTKBinaryFileWriter::createVertexDataWriter(const std::string& identifier, int recordsPerVertex) {
  validateDataWriterIdentifier(identifier);
  return new tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter::VertexDataWriter(
    identifier, *this, recordsPerVertex
  );
}
