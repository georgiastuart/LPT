ifeq ($(PLATFORM),lonestar-intel)
  ifeq ($(BUILDTYPE),lpt_serial)
    COMPILER := ifort
    FLAGS := -DVERBOSE=4 -DDEBUG -check all -traceback
  endif
  ifeq ($(BUILDTYPE),lpt_parallel)
    COMPILER := mpif90
    FLAGS := -DVERBOSE=1 -DMPI # -traceback # -check all -DDEBUG # -DDEBUG_MPI
  endif
  ifeq ($(NETCDF),enable)
    NETCDF_FLAGS := -DNETCDF=4 -I${TACC_NETCDF_INC} -L${TACC_NETCDF_LIB} -lnetcdf -lnetcdff -L${TACC_HDF5_LIB} -Wl,-rpath,${TACC_HDF5_LIB} -L${TACC_HDF5_LIB} -lhdf5 -lz
    HDF5_FLAGS :=
  endif
endif

SOURCE_FILES = lpt_comm.F90 lpt_data.F90 lpt_kdtree.F90 lpt_netcdf.F90 lpt_read.F90 lpt_write.F90 lpt_main.F90 lpt_print.F90 lpt_drog.F90 lpt_oil.F90
SOURCE_OBJECTS := $(patsubst %.F, %.o, $(SOURCE_FILES) )

%.o : %.F90
	$(COMPILER) $(FLAGS) -c $< 

all : lpt_serial lpt_parallel

ifeq ($(MAKELEVEL),0)
  lpt_serial :
	$(MAKE) BUILDTYPE=lpt_serial $@
  lpt_parallel :
	$(MAKE) BUILDTYPE=lpt_parallel $@
else
  lpt_serial :: $(SOURCE_OBJECTS)
	$(COMPILER) $(FLAGS) $^ $(NETCDF_FLAGS) $(HDF5_FLAGS) -o $@
  lpt_parallel :: $(SOURCE_OBJECTS)
	$(COMPILER) $(FLAGS) $^ $(NETCDF_FLAGS) $(HDF5_FLAGS) -o $@
endif

clean :
	rm -f *.mod lpt_serial lpt_parallel
