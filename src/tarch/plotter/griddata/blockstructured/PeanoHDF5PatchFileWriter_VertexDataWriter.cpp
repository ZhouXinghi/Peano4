#include "tarch/plotter/griddata/blockstructured/PeanoHDF5PatchFileWriter.h"


tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter::VertexDataWriter::VertexDataWriter(
  [[maybe_unused]] const std::string&                                                   identifier,
  [[maybe_unused]] int                                                                  unknownsPerAxis,
  [[maybe_unused]] int                                                                  numberOfUnknowns,
  [[maybe_unused]] const std::string&                                                   description,
  [[maybe_unused]] const std::string&                                                   metaData,
  [[maybe_unused]] double*                                                              mapping,
  [[maybe_unused]] tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter& writer
):
  _writer(writer),
  _identifier(identifier),
  _numberOfUnknowns(numberOfUnknowns) {
#ifdef UseHDF5
  logDebug("VertexDataWriter(...)", "create numberofunknowns entry");

  /**
   * Create scalar attribute.
   */
  hid_t numberOfUnknownsDataSpace = H5Screate(H5S_SCALAR);
  hid_t numberOfUnknownsAttribute = H5Acreate(
    _writer._file,
    (_writer.getNameOfCurrentDataset() + "/vertexdata/numberofunknowns/" + _identifier).c_str(),
    H5T_NATIVE_INT,
    numberOfUnknownsDataSpace,
    H5P_DEFAULT,
    H5P_DEFAULT
  );

  /**
   * Write scalar attribute.
   */
  H5Awrite(numberOfUnknownsAttribute, H5T_NATIVE_INT, &numberOfUnknowns);

  H5Aclose(numberOfUnknownsAttribute);
  H5Sclose(numberOfUnknownsDataSpace);
#endif

  if (!metaData.empty()) {
#ifdef UseHDF5
    hid_t metaDataAttribute = H5Screate(H5S_SCALAR);
    hid_t metaDataType      = H5Tcopy(H5T_C_S1);

    H5Tset_size(metaDataType, metaData.size());
    H5Tset_strpad(metaDataType, H5T_STR_NULLTERM);
    hid_t attribute = H5Acreate2(
      _writer._file,
      (_writer.getNameOfCurrentDataset() + "/vertexdata/metadata/" + _identifier).c_str(),
      metaDataType,
      metaDataAttribute,
      H5P_DEFAULT,
      H5P_DEFAULT
    );

    /*
     * Write string attribute.
     */
    H5Awrite(attribute, metaDataType, metaData.c_str());

    /*
     * Close attribute and file dataspaces, and datatype.
     */
    H5Aclose(attribute);
    H5Sclose(metaDataAttribute);
#endif
  }

  if (mapping != nullptr) {
#ifdef UseHDF5
    //
    // Create the data space with unlimited dimensions.
    //
    hsize_t tableDimensions[] = {
      static_cast<hsize_t>(_writer._dimensions),
      static_cast<hsize_t>(std::pow(_numberOfUnknowns, _writer._dimensions))};

    //
    // Set up handles/tables
    //
    hid_t dataTable = H5Screate_simple(2, tableDimensions, NULL);
    hid_t dataset   = H5Dcreate(
      _writer._file,
      (_writer.getNameOfCurrentDataset() + "/vertexdata/mapping/" + _identifier).c_str(),
      H5T_NATIVE_DOUBLE,
      dataTable,
      H5P_DEFAULT,
      H5P_DEFAULT,
      H5P_DEFAULT
    );

    //
    // Write data
    //
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mapping);

    //
    // Close/release handles
    //
    H5Dclose(dataset);
    H5Sclose(dataTable);
#endif
  }
}


tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter::VertexDataWriter::~VertexDataWriter() {}


void tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter::VertexDataWriter::plotVertex(
  int index, double value
) {
  while (static_cast<int>(_data.size()) < (index + 1) * _numberOfUnknowns) {
    _data.resize((index + 1) * _numberOfUnknowns);
  }

  _data[index * _numberOfUnknowns] = value;

  for (int i = 1; i < _numberOfUnknowns; i++) {
    _data[index * _numberOfUnknowns + i] = 0.0;
  }
}


void tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter::VertexDataWriter::plotVertex(
  [[maybe_unused]] int index, [[maybe_unused]] double* values
) {
  /*
    while (static_cast<int>(_data.size())<(index+1)*_numberOfUnknowns) {
      _data.resize( (index+1)*_numberOfUnknowns );
    }

    for (int i=0; i<numberOfValues; i++) {
      _data[index*_numberOfUnknowns+i] = values[i];
    }

    for (int i=numberOfValues; i<_numberOfUnknowns; i++) {
      _data[index*_numberOfUnknowns+i] = 0.0;
    }
  */
}

void tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter::VertexDataWriter::close() {
  assignRemainingVerticesDefaultValues();

#ifdef UseHDF5
  logDebug("close()", "create data table of " << _identifier);

  const int lineLenght = std::pow(_numberOfUnknowns, _writer._dimensions);
  assertion1(_data.size() % lineLenght == 0, _identifier);

  //
  // Create the data space with unlimited dimensions.
  //
  hsize_t tableDimensions[] = {static_cast<hsize_t>(lineLenght), static_cast<hsize_t>(_data.size()) / lineLenght};

  //
  // Set up handles/tables
  //
  hid_t dataTable = H5Screate_simple(2, tableDimensions, NULL);
  hid_t dataset   = H5Dcreate(
    _writer._file,
    (_writer.getNameOfCurrentDataset() + "/vertexdata/data/" + _identifier).c_str(),
    H5T_NATIVE_DOUBLE,
    dataTable,
    H5P_DEFAULT,
    _writer.createDataTableProperties(lineLenght, static_cast<int>(_data.size()) / lineLenght),
    H5P_DEFAULT
  );

  logDebug("close()", "write data table of " << _identifier);

  //
  // Write data
  //
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, _data.data());

  //
  // Close/release handles
  //
  H5Dclose(dataset);
  H5Sclose(dataTable);
#endif
}


void tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter::VertexDataWriter::
  assignRemainingVerticesDefaultValues() {
  while (static_cast<int>(_data.size()) < _writer._vertexCounter * _numberOfUnknowns) {
    _data.push_back(0.0);
  }
}


int tarch::plotter::griddata::blockstructured::PeanoHDF5PatchFileWriter::VertexDataWriter::getFirstVertexWithinPatch(
  [[maybe_unused]] int index
) const {
  return -1;
}
