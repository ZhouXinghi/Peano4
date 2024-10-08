find_program(SLURM_SBATCH_COMMAND sbatch DOC "Path to the SLURM sbatch executable")
find_program(SLURM_SRUN_COMMAND srun DOC "Path to the SLURM srun executable")
mark_as_advanced(SLURM_SRUN_COMMAND SLURM_SBATCH_COMMAND)

if(SLURM_SRUN_COMMAND AND SLURM_SBATCH_COMMAND)
  set(SLURM_FOUND TRUE)
  if(NOT SLURM_FIND_QUIETLY)
    message(STATUS "Found Slurm")
  endif ()
else ()
  set(SLURM_FOUND FALSE)
  if(SLURM_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find Slurm")
  elseif(NOT SLURM_FIND_QUIETLY)
    message(WARNING "Could not find Slurm")
  endif()
endif()
