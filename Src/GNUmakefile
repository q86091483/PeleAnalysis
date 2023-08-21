AMREX_HOME          ?= /home/b/bsavard/zisen347/codes/exascale/amrex
PELE_PHYSICS_HOME   ?= /home/b/bsavard/zisen347/codes/exascale/PelePhysics
PELE_ANALYSIS_HOME  ?= /scratch/b/bsavard/zisen347/PeleAnalysis

# AMReX
DEBUG	      = FALSE
DIM	        = 3
PRECISION   = DOUBLE

# Compilation
COMP          = intel
USE_MPI       = TRUE
USE_OMP       = FALSE
USE_CUDA      = FALSE

Chemistry_Model = DME_30atm_new
#Reactions_dir   = Fuego
Eos_Model       = Fuego
Transport_Model = Simple
Reactions_dir := Fuego

# Executable
# Uncomment the desired exec
#EBASE         = grad
#EBASE         = combinePlts
#EBASE         = isosurface
#EBASE         = flattenAMRFile
#EBASE         = stream
#EBASE         = partStream
#EBASE         = blowOut
#EBASE         = stream2plt
#EBASE         = sampleStreamlines
#EBASE         = streamScatter
#EBASE         = streamSub
#EBASE         = streamTubeStats
#EBASE         = curvature
#EBASE         = curvature_old
#EBASE         = surfMEFtoDAT
#EBASE         = surfDATtoMEF
#EBASE         = buildPMF
#EBASE         = subPlt
#EBASE         = conditionalMean
#EBASE         = jpdf
#EBASE         = avgToPlane
#EBASE          = avgPlotfiles
#EBASE          = amrToFE
#EBASE         = conditionalISRN
EBASE         = gradISRN

NEEDS_f90_SRC = FALSE
#NEEDS_f90_SRC = TRUE

Pdirs   := Base Boundary AmrCore Extern/amrdata #LinearSolvers/MLMG
Ppack   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

PPdirs  := Utility/PltFileManager
Ppack   += $(foreach dir, $(PPdirs), $(PELE_PHYSICS_HOME)/$(dir)/Make.package)

ifeq ($(EBASE),partStream)
  USE_PARTICLES = TRUE
  Pdirs += Particle
  CEXE_headers += StreamData.H    StreamPC.H
  CEXE_sources += StreamData.cpp  StreamPC.cpp
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp
ifeq ($(NEEDS_f90_SRC),TRUE)
  f90EXE_sources += $(EBASE)_nd.f90
endif
INCLUDE_LOCATIONS += .
VPATH_LOCATIONS   += .

ifeq ($(EBASE), conditionalISRN)
  #INCLUDE_LOCATIONS += $(HDF5_HOME)/include
  LIBRARIES         += -lhdf5_cpp -lhdf5
endif

ifeq ($(EBASE), gradISRN)
  #INCLUDE_LOCATIONS += $(HDF5_HOME)/include
  LIBRARIES         += -lhdf5_cpp -lhdf5
endif

include $(Ppack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

SDF_LOC = ../Tools/SDFGen
include $(SDF_LOC)/Make.package
INCLUDE_LOCATIONS += $(SDF_LOC)
VPATH_LOCATIONS += $(SDF_LOC)

vpath %.c   : $(VPATH_LOCATIONS)
vpath %.h   : $(VPATH_LOCATIONS)
vpath %.cpp : $(VPATH_LOCATIONS)
vpath %.H   : $(VPATH_LOCATIONS)
vpath %.F   : $(VPATH_LOCATIONS)
vpath %.f   : $(VPATH_LOCATIONS)
vpath %.f90 : $(VPATH_LOCATIONS)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.Analysis
