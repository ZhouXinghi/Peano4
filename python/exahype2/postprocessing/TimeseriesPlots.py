# This file is part of the ExaHyPE2 project. For conditions of distribution and 
# use, please see the copyright notice at www.peano-framework.org
from .PerformanceData import PerformanceData
from .utils           import next_symbol
from .utils           import next_markevery

import matplotlib.pyplot as plt


from enum import Enum
class XAxis(Enum):
  Iterations    = 1
  RealTime      = 2
  SimulatedTime = 3



def plot_runtime_per_step(performance_data,label=None,xaxis=XAxis.RealTime):
  """
  
  Label: String or None
  
  """    
  x_data = []
  y_data = []
  if xaxis.value>1:
    i=0
    while i<len(performance_data.get_time_step_real_time_stamps()) and i<len(performance_data.get_time_per_time_step()):
      if 1:#i%xaxis.value==0:
        x_data.append(performance_data.get_time_step_real_time_stamps()[i])
        y_data.append(performance_data.get_time_per_time_step()[i])
      else:
        #x_data[-1]  = performance_data.get_time_step_time_stamps()[i]
        #y_data[-1] += performance_data.get_time_per_time_step()[i]
        pass
      i+=1
      
    x_data = x_data[0:-1]
    y_data = y_data[0:-1]
  else:
    x_data = performance_data.get_time_step_real_time_stamps()
    y_data = performance_data.get_time_per_time_step()
  plt.plot( x_data, y_data, next_symbol(), markevery=next_markevery(performance_data.timesteps()), label=label  )
  
  if len(performance_data.plotting_time_stamp)>0:
    max_time_per_time_step = max(performance_data.get_time_per_time_step())
    min_time_per_time_step = min(performance_data.get_time_per_time_step())
    for i in performance_data.plotting_time_stamp:
      plt.plot( [i,i], [min_time_per_time_step,max_time_per_time_step], "--", color="#000000"  )
    plt.plot( [performance_data.plotting_time_stamp[0],performance_data.plotting_time_stamp[0]], [min_time_per_time_step,max_time_per_time_step], "--", color="#000000", label="plot"  )
    
  plt.xlabel( "Runtime [t]=s" )
  plt.ylabel( "Time per step [t]=s" )


def plot_time_step_size_per_step(performance_data,label=None,xaxis=XAxis.RealTime):
  """
       
  """
  if xaxis==XAxis.RealTime and not performance_data.uses_local_timestepping():
    if len(performance_data.get_time_step_real_time_stamps()) != len(performance_data.get_time_step_time_step_size()):
      raise Exception( "Size of fields do not match: " + str(len(performance_data.get_time_step_real_time_stamps())) + " vs. " + str(len(performance_data.get_time_step_time_step_size())))
    plt.plot( performance_data.get_time_step_real_time_stamps(), performance_data.get_time_step_time_step_size(), next_symbol(), markevery=next_markevery(performance_data.timesteps()), label=label  )
    plt.xlabel( "Real time [t]=s" )
  elif xaxis==XAxis.RealTime and performance_data.uses_local_timestepping():
    symbol = next_symbol()
    x_data = performance_data.get_time_step_real_time_stamps()
    y_data = performance_data.get_time_step_time_step_size()[0]
    if len(x_data)>len(y_data):
      print( "WARNING: too many time stamps (does not fit to data points). Remove first entry")
      x_data = x_data[1:]
    plt.plot( x_data, y_data, symbol, markevery=next_markevery(len(y_data)), label=label  )
    y_data = performance_data.get_time_step_time_step_size()[1]
    if len(x_data)>len(y_data):
      print( "WARNING: too many time stamps (does not fit to data points). Remove first entry")
      x_data = x_data[1:]
    plt.plot( x_data, y_data, "-"+symbol, markevery=next_markevery(len(y_data)), label=label  )
    plt.xlabel( "Real time [t]=s" )
  elif xaxis==XAxis.Iterations and not performance_data.uses_local_timestepping():
    y_data = performance_data.get_time_step_time_step_size()
    plt.plot( range(0,len(y_data)), y_data, next_symbol(), markevery=next_markevery(len(y_data)), label=label  )
    plt.xlabel( "Simulation step" )
  elif xaxis==XAxis.Iterations and performance_data.uses_local_timestepping():
    y_data = performance_data.get_time_step_time_step_size()[0]
    symbol = next_symbol()
    plt.plot( range(0,len(y_data)), y_data, symbol, markevery=next_markevery(len(y_data)), label=label  )
    y_data = performance_data.get_time_step_time_step_size()[1]
    plt.plot( range(0,len(y_data)), y_data, "-" + symbol, markevery=next_markevery(len(y_data)), label=label  )
    plt.xlabel( "Simulation step" )
  elif xaxis==XAxis.SimulatedTime and not performance_data.uses_local_timestepping():
    if len(performance_data.get_time_step_simulated_time_stamps()) != len(performance_data.get_time_step_time_step_size()):
      raise Exception( "Size of fields do not match: " + str(len(performance_data.get_time_step_simulated_time_stamps())) + " vs. " + str(len(performance_data.get_time_step_time_step_size())))
    plt.plot( performance_data.get_time_step_simulated_time_stamps(), performance_data.get_time_step_time_step_size(), next_symbol(), markevery=next_markevery(performance_data.timesteps()), label=label  )
    plt.xlabel( "Simulated time" )
  elif xaxis==XAxis.SimulatedTime and performance_data.uses_local_timestepping():
    if len(performance_data.get_time_step_simulated_time_stamps()) != len(performance_data.get_time_step_time_step_size()):
      raise Exception( "Size of fields do not match: " + str(len(performance_data.get_time_step_simulated_time_stamps())) + " vs. " + str(len(performance_data.get_time_step_time_step_size())))
    symbol = next_symbol()
    x_data = performance_data.get_time_step_simulated_time_stamps()[0]
    y_data = performance_data.get_time_step_time_step_size()[0]
    x_data = x_data[len(x_data)-len(y_data):]
    plt.plot( x_data, y_data, symbol, markevery=next_markevery(len(y_data)), label=label  )
    y_data = performance_data.get_time_step_time_step_size()[1]
    plt.plot( x_data, y_data, "-"+symbol, markevery=next_markevery(len(y_data)), label=label  )
    plt.xlabel( "Simulated time $T$" )
  else:
    raise Exception( "enum value not supported" )
  plt.ylabel( "Time step size $\Delta T$" )



def plot_updates_per_step(performance_data,label=None,xaxis=XAxis.RealTime):
  """
  
  """
  if xaxis==XAxis.RealTime:
    plt.plot( performance_data.get_time_step_real_time_stamps(), performance_data.get_upates(), next_symbol(), markevery=next_markevery(performance_data.get_updates()), label=label  )
    plt.xlabel( "Real time [t]=s" )
  elif xaxis==XAxis.Iterations:
    y_data = performance_data.get_updates()
    plt.plot( range(0,len(y_data)), y_data, next_symbol(), markevery=next_markevery(len(y_data)), label=label  )
    plt.xlabel( "Simulation step" )
  elif xaxis==XAxis.SimulatedTime:
    if performance_data.uses_local_timestepping():
      x_data = performance_data.get_time_step_simulated_time_stamps()[0]
    else:
      x_data = performance_data.get_time_step_simulated_time_stamps()
    plt.plot( x_data, performance_data.get_updates(), next_symbol(), markevery=next_markevery(performance_data.timesteps()), label=label  )
    plt.xlabel( "Simulated time" )
  else:
    raise Exception( "enum value not supported" )
  plt.ylabel( "Updates" )
