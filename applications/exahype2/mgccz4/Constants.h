//
// Peano4 data file
// Generated by Peano's Python API
// www.peano-framework.org
// This is generated. Be careful with adding your own stuff
//
#ifndef __EXAMPLES_EXAHYPE2_MGCCZ4_CONSTANTS__
#define __EXAMPLES_EXAHYPE2_MGCCZ4_CONSTANTS__


#include <string>


#include <bitset>



namespace examples{
namespace exahype2{
namespace mgccz4{

  const std::initializer_list<double> DomainOffset = {-0.5,-0.5,-0.5};
  const std::initializer_list<double> DomainSize = {1.0,1.0,1.0};
  constexpr auto TerminalTime = 1.0;
  constexpr auto FirstPlotTimeStamp = 0.0;
  constexpr auto TimeInBetweenPlots = 0.04;
  constexpr auto PlotterPrecision = 5;
  const std::bitset<3> PeriodicBC = 0+1+2+4;

}
}
}


#endif

