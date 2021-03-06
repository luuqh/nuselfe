################################################################################
# Parallel SELFE Makefile
#
# User makes environment settings for particular OS / PLATFORM / COMPILER / MPI
# below. These settings are based on the following options.
#
# Compiler name:
#   FCS: Serial compiler (for utilities)
#   FCP: Parallel compiler
#   FLD: Linker (in general same as parallel compiler)
#
# Compilation flags
#   FCSFLAGS: Flags for serial compilation
#   FCPFLAGS: Flags for parallel compilation (including all pre-processing flags)
#   FLDFLAGS: Flags for linker (e.g., -O2)
#
# Preprocessor flags:
#   DEBUG: Enable debugging code
#   ORDERED_SUM: Enable globally ordered sums & dot-products for bit reproducibility
#     of state quantities independent of number of processors (note: this can
#     significantly degrade performance);
#   INCLUDE_TIMING: Enable wallclock timing of code (note: this can have slight
#     effect on performance);
#   MPI_VERSION = 1 or 2: Version of MPI (try 2 first, if compile fails due to mpi
#     related errors then switch to version 1;
#
# Libraries (needed for parallel code)
#   MTSLIBS: Flags for linking ParMeTiS/MeTiS libaries
#   ALTLIBS: Flags for linking alternate solver libraries (LAPACK or ITPACK,
#            these are just for testing)
################################################################################


################################################################################
# Comment out the following line if you don't want the build dependent
# upon the Makefile, otherwise use the name of the Makefile
################################################################################
#MAKEFILE = Makefile

#ENV         = XEON
#ENV         = CANOPUS
#ENV         = SIRIUS
#ENV          = FEDORA
ENV          = CETO
################################################################################
# Define pre-processing flags.
# Comment out unwanted flags.
################################################################################
# -DSELFE always on
  PPFLAGS = -DSELFE 

# Wave models need a coupler
# SWAN
#  USE_SWAN = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_SWAN
  #MPI version
#  SWAN_DIR = /disk/ambcs03/0/yinglong/SWAN/swan4072_MPI/

# Comment out if timer report is not needed.
#  PPFLAGS := $(PPFLAGS) -DINCLUDE_TIMING

# Precip/evaporation model
#  PPFLAGS := $(PPFLAGS) -DPREC_EVAP 

# MM5 in heat exchange model
#  PPFLAGS := $(PPFLAGS) -DMM5

# GOTM turbulence closure model
#  USE_GOTM = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_GOTM 

# Ecological model
#  USE_ECO = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_ECO

# Sediment model
#  USE_SED = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_SED

# Oil spill model
#  USE_OIL = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_OIL

# Debug mode (more time consuming)
#  PPFLAGS := $(PPFLAGS) -DDEBUG 



################################################################################
# Environment for Linux / 32 bit /  Intel Compiler / MPICH
################################################################################
ifeq ($(ENV),XEON)
  FCP = mpif90 #MPI compiler
  FLD = $(FCP)
# MPI vserion (1 or 2) 
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=1 
  FCPFLAGS = $(PPFLAGS) -O2 -Bstatic -assume byterecl #MPI code
  FLDFLAGS = -O2 #for final linking of object files
  #####Libraries
  MTSLIBS = -L/home/users/yinglong/SELFE/ParMetis-3.1-32bit -lparmetis -lmetis
  CDFLIBS = -L/usr/local/netcdf/lib -lnetcdf 
  CDFMOD = -I/usr/local/netcdf/include # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =  -I/home/users/yinglong/GOTM/gotm-3.2.5/32bit/gotm-3.2.5/modules/IFORT/ #modules
    GTMLIBS = -L/home/users/yinglong/GOTM/gotm-3.2.5/32bit/gotm-3.2.5/lib/IFORT/ -lturbulence_prod  -lutil_prod 
  else
    GTMMOD =
    GTMLIBS =
  endif
endif

################################################################################
# Environment for Linux / 64 bit /  Intel Compiler / MPICH (CETO)
################################################################################
ifeq ($(ENV),CETO)
  #FCP = /home/yinglong/mpich2/bin/mpif90
  FCP = /opt/mpich2/mpich2-1.0.8-intel/bin/mpif90
  FLD = $(FCP)
# MPI vserion (1 or 2) 
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=2 
  FCPFLAGS = $(PPFLAGS) -O3 -Bstatic -assume byterecl #MPI code
  FLDFLAGS = -O2 #for final linking of object files
  #####Libraries
  MTSLIBS = -L/data1/tmsxh/ParMetis-3.1-64bit/ -lparmetis -lmetis
  CDFLIBS = -L/usr/local/netcdf4-intel/lib -lnetcdff  
  CDFMOD = -I/usr/local/netcdf4-intel/include/ # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =  -I/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/modules/IFORT/ #modules
    GTMLIBS = -L/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/lib/IFORT/ -lturbulence_prod  -lutil_prod 
  else
    GTMMOD  =
    GTMLIBS =
  endif
endif
################################################################################
# Environment for FEDORA Desktop
################################################################################
ifeq ($(ENV),FEDORA)
  USE_WRAP = yes
  FCP = mpif90 -f90=ifort
  FLD = $(FCP)
# MPI vserion (1 or 2) 
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=2 -DUSE_WRAP
  FCPFLAGS = $(PPFLAGS) -O2 -Bstatic -assume byterecl #MPI code
  FLDFLAGS = -O2 #for final linking of object files
  #####Libraries
  MTSLIBS = -L/opt/ParMetis-3.1.1/ -lparmetis -lmetis
  MTSMOD = -I/opt/ParMetis-3.1.1/     #wrapc.c uses parmetis.h
  CDFLIBS = -L/home/trim/intel/lib/ -lnetcdff 
  CDFMOD = -I/home/trim/intel/include/ # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =  -I/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/modules/IFORT/ #modules
    GTMLIBS = -L/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/lib/IFORT/ -lturbulence_prod  -lutil_prod
  else
    GTMMOD =
    GTMLIBS =
  endif
endif

################################################################################
# Environment for SIRIUS cluster
################################################################################
ifeq ($(ENV),SIRIUS)
  USE_WRAP = yes
  FCP = mpif90
  FLD = $(FCP)
# MPI vserion (1 or 2) 
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=2 -DUSE_WRAP
  FCPFLAGS = $(PPFLAGS) -O2 -Bstatic -assume byterecl #MPI code
  FLDFLAGS = -O2 #for final linking of object files
  #####Libraries
  MTSLIBS = -L/share/apps/ParMetis/ -lparmetis -lmetis
  MTSMOD = -I/share/apps/ParMetis/ #wrapc.c uses parmetis.h
  CDFLIBS = -L/share/apps/netcdf/lib/ -lnetcdf -L/opt/intel/fce/10.1.015/lib/ -lirc
  CDFMOD = -I/share/apps/netcdf/include/ # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =  -I/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/modules/IFORT/ #modules
    GTMLIBS = -L/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/lib/IFORT/ -lturbulence_prod  -lutil_prod
  else
    GTMMOD =
    GTMLIBS =
  endif
endif

################################################################################
# Executable and source code names (without extensions)
################################################################################
EXEC   = pelfe_3.0b_jet_netcdf_v4.1
MAIN   = elfe_main
GLBL   = elfe_glbl
MSGP   = elfe_msgp
GRID   = grid_subs
BKTRK  = bktrk_subs
SOLVER = solver_subs
IO     = io_subs
MISC   = misc_subs
READ_PARAM = read_param
TRANS  = transport
SFLUX  = sflux_9c
LAPACK = lap
WRAPC  = wrapc

# SWAN
COUPLER_SWAN = coupler_elfe_swan

# Ecology
BIOINI     = bio_init
BIOL       = biology
BIOPAR     = bio_param
ECLIG      = eclight
ECOSIM     = ecosim
ECORD      = read_eco_input
SPEC       = specir

#Netcdf 
NSF_VAR   = nsf_self_vars
NSF_CON   = nsf_container
CHECK_ERR = check_error
DATA_TYPE = data_type_mod
################################################################################
# Object files; modules compile first; modules that depend on other modules compile last.
################################################################################
OBJS = $(GLBL).o \
       $(MSGP).o 

ifdef USE_ECO
  OBJS := $(OBJS) $(BIOPAR).o $(ECLIG).o $(BIOL).o $(BIOINI).o $(ECOSIM).o $(ECORD).o $(SPEC).o
endif

OBJS := $(OBJS) \
       $(MAIN).o \
       $(GRID).o \
       $(IO).o \
       $(MISC).o \
       $(READ_PARAM).o \
       $(SOLVER).o \
       $(BKTRK).o \
       $(TRANS).o \
       $(SFLUX).o \
       $(LAPACK).o

ifdef USE_WRAP
  OBJS := $(OBJS) $(WRAPC).o
endif

ifdef USE_SWAN
  SWAN_OBJS = \
  $(SWAN_DIR)/swmod1.o \
  $(SWAN_DIR)/swmod2.o \
  $(SWAN_DIR)/m_constants.o \
  $(SWAN_DIR)/m_fileio.o \
  $(SWAN_DIR)/serv_xnl4v5.o \
  $(SWAN_DIR)/mod_xnl4v5.o \
  $(SWAN_DIR)/SwanGriddata.o \
  $(SWAN_DIR)/SwanGridobjects.o \
  $(SWAN_DIR)/SwanCompdata.o \
  $(SWAN_DIR)/swanmain.o \
  $(SWAN_DIR)/swanpre1.o \
  $(SWAN_DIR)/swanpre2.o \
  $(SWAN_DIR)/swancom1.o \
  $(SWAN_DIR)/swancom2.o \
  $(SWAN_DIR)/swancom3.o \
  $(SWAN_DIR)/swancom4.o \
  $(SWAN_DIR)/swancom5.o \
  $(SWAN_DIR)/swanout1.o \
  $(SWAN_DIR)/swanout2.o \
  $(SWAN_DIR)/swanser.o \
  $(SWAN_DIR)/swanparll.o \
  $(SWAN_DIR)/SwanReadGrid.o \
  $(SWAN_DIR)/SwanReadADCGrid.o \
  $(SWAN_DIR)/SwanReadTriangleGrid.o \
  $(SWAN_DIR)/SwanReadEasymeshGrid.o \
  $(SWAN_DIR)/SwanInitCompGrid.o \
  $(SWAN_DIR)/SwanCheckGrid.o \
  $(SWAN_DIR)/SwanCreateEdges.o \
  $(SWAN_DIR)/SwanGridTopology.o \
  $(SWAN_DIR)/SwanGridVert.o \
  $(SWAN_DIR)/SwanGridCell.o \
  $(SWAN_DIR)/SwanGridFace.o \
  $(SWAN_DIR)/SwanPrintGridInfo.o \
  $(SWAN_DIR)/SwanFindPoint.o \
  $(SWAN_DIR)/SwanPointinMesh.o \
  $(SWAN_DIR)/SwanBpntlist.o \
  $(SWAN_DIR)/SwanPrepComp.o \
  $(SWAN_DIR)/SwanVertlist.o \
  $(SWAN_DIR)/SwanCompUnstruc.o \
  $(SWAN_DIR)/SwanDispParm.o \
  $(SWAN_DIR)/SwanPropvelX.o \
  $(SWAN_DIR)/SwanSweepSel.o \
  $(SWAN_DIR)/SwanPropvelS.o \
  $(SWAN_DIR)/SwanTranspAc.o \
  $(SWAN_DIR)/SwanTranspX.o \
  $(SWAN_DIR)/SwanInterpolatePoint.o \
  $(SWAN_DIR)/SwanInterpolateAc.o \
  $(SWAN_DIR)/SwanInterpolateOutput.o \
  $(SWAN_DIR)/SwanConvAccur.o \
  $(SWAN_DIR)/SwanConvStopc.o \
  $(SWAN_DIR)/SwanFindObstacles.o \
  $(SWAN_DIR)/SwanCrossObstacle.o \
  $(SWAN_DIR)/SwanComputeForce.o \
  $(SWAN_DIR)/SwanIntgratSpc.o \
  $(SWAN_DIR)/ocpids.o \
  $(SWAN_DIR)/ocpcre.o \
  $(SWAN_DIR)/ocpmix.o \

  OBJS := $(OBJS) $(COUPLER_SWAN).o $(SWAN_OBJS)
endif

#NETCDF
OBJS := $(OBJS) nsf_self_vars.o nsf_container.o check_error.o data_type_mod.o jet.o
#OBJS := $(OBJS) nsf_self_vars.o nsf_container.o check_error.o data_type_mod.o 
################################################################################
# Library list
################################################################################
LIBS = $(MTSLIBS) $(CDFLIBS) $(GTMLIBS)


################################################################################
# Build commands for SELFE
################################################################################
#all: $(EXEC) convert_param read_output5 ptrack2

ifdef USE_SWAN
  FLDFLAGS := $(FLDFLAGS) -ipo -xW -mp1 -W0 -assume byterecl -traceback
endif

$(EXEC): $(OBJS) $(MAKEFILE)
	$(FLD) $(FLDFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

$(GLBL).o: $(GLBL).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(GLBL).F90

$(MSGP).o: $(MSGP).F90 $(GLBL).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(MSGP).F90

#NETCDF
data_type_mod.o: data_type_mod.F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(CDFMOD) data_type_mod.F90

check_error.o: check_error.F90 data_type_mod.o $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(CDFMOD) check_error.F90

nsf_container.o: nsf_container.F90 check_error.o $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(CDFMOD) nsf_container.F90

nsf_self_vars.o: nsf_self_vars.F90 nsf_container.o $(GLBL).o $(MSGP).o $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(CDFMOD) nsf_self_vars.F90

jet.o: jet.F90  $(GLBL).o $(MSGP).o $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(CDFMOD) jet.F90

$(MAIN).o: $(MAIN).F90 $(GLBL).F90 $(MSGP).F90 nsf_self_vars.o $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(GTMMOD) $(CDFMOD) $(MAIN).F90

$(GRID).o: $(GRID).F90 $(GLBL).F90 $(MSGP).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(GRID).F90

$(BKTRK).o: $(BKTRK).F90 $(GLBL).F90 $(MSGP).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(BKTRK).F90

$(SOLVER).o: $(SOLVER).F90 $(GLBL).F90 $(MSGP).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(SOLVER).F90

$(IO).o: $(IO).F90 $(GLBL).F90 $(MSGP).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(IO).F90

$(MISC).o: $(MISC).F90 $(GLBL).F90 $(MSGP).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(MISC).F90

$(READ_PARAM).o: $(READ_PARAM).F90 $(GLBL).F90 $(MSGP).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(READ_PARAM).F90

$(TRANS).o: $(TRANS).F90 $(GLBL).F90 $(MSGP).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(TRANS).F90

$(SFLUX).o: $(SFLUX).F90 $(GLBL).F90 $(MSGP).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(CDFMOD) $(SFLUX).F90
$(LAPACK).o: $(LAPACK).F90 $(GLBL).F90 $(MSGP).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(LAPACK).F90

$(WRAPC).o: $(WRAPC).c $(MAKEFILE)
	mpicc -c $(WRAPC).c $(MTSMOD)

# SWAN-SELFE coupler
$(COUPLER_SWAN).o: $(COUPLER_SWAN).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(COUPLER_SWAN).F90
	#$(FCP) -c -O2 $(COUPLER_SWAN).F90

# Ecology
$(BIOINI).o: $(BIOINI).F90 $(GLBL).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(BIOINI).F90

$(BIOL).o: $(BIOL).F90 $(GLBL).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(BIOL).F90

$(BIOPAR).o: $(BIOPAR).F90 $(GLBL).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(BIOPAR).F90

$(ECLIG).o: $(ECLIG).F90 $(GLBL).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(ECLIG).F90

$(ECOSIM).o: $(ECOSIM).F90 $(GLBL).F90 $(MSGP).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(ECOSIM).F90

$(ECORD).o: $(ECORD).F90 $(GLBL).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(ECORD).F90

$(SPEC).o: $(SPEC).F90 $(GLBL).F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(SPEC).F90

#NETCDF

#$(NSF_VAR).o: $(NSF_VAR).F90 $(NSF_CON).o $(GLBL).o $(MSGP).o $(MAKEFILE)
#	$(FCP) -c $(FCPFLAGS) $(CDFMOD) $(NSF_VAR).F90

#$(NSF_CON).o: $(NSF_CON).F90 $(CHECK_ERR).o $(MAKEFILE)
#	$(FCP) -c $(FCPFLAGS) $(CDFMOD) $(NSF_CON).F90

#$(CHECK_ERR).o: $(CHECK_ERR).F90 $(DATA_TYPE).o $(MAKEFILE)
#	$(FCP) -c $(FCPFLAGS) $(CDFMOD) $(CHECK_ERR).F90

#$(DATA_TYPE).o: $(DATA_TYPE).F90
#	$(FCP) -c $(FCPFLAGS) $(CDFMOD) $(DATA_TYPE).F90



################################################################################
# Build commands for utilities
################################################################################
#convert_param: convert_param.f90 $(MAKEFILE)
#	$(FCS) $(FCSFLAGS) convert_param.f90 -o convert_param
#
#read_output5: read_output5.f90 $(MAKEFILE)
#	$(FCS) $(FCSFLAGS) read_output5.f90 -o read_output5
#
#ptrack2: ptrack2.f90 $(MAKEFILE)
#	$(FCS) $(FCSFLAGS) ptrack2.f90 -o ptrack2


################################################################################
# Clean up or Clobber
################################################################################
# Clean up objects and modules
clean:
	rm -f *.o *.mod $(EXEC)
