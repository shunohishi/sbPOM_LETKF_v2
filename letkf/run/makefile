#
# sbPOM makefile
#

#-----------------------------------------------------------------------
# Settings that depend on the system and the compiler
#-----------------------------------------------------------------------
# Set macros
CPP = cpp -P

### Fugaku ###
FC=mpifrtpx
LD=mpifrtpx

CLEAN = rm

# Set libraries and include files
#----Fugaku
#NETCDFINC = -I/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-c-4.9.2-cczsjh6lmalmjqhr72yatvxzmc3iwdl2/include
#NETCDFFINC = -I/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-fortran-4.6.1-kjm7jl5naxulm7neknsgxnsbkf7iet5j/include
#NETCDFLIB = -L/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-c-4.9.2-cczsjh6lmalmjqhr72yatvxzmc3iwdl2/lib -lnetcdf
#NETCDFFLIB = -L/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-fortran-4.6.1-kjm7jl5naxulm7neknsgxnsbkf7iet5j/lib -lnetcdff  -lnetcdf -lnetcdf -lm
#HDF5INC = 
#HDF5LIB =
#STATIC = -lhdf5_hl -lhdf5 -lm -lsz -lbz2 -lzstd -lblosc -lmpi_cxx -lfjprofmpi

#FFLAGS = -Kfast -KSVE -Kparallel -Kopenmp -Kcmodel=large -Nfjomplib -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
#LIBS = -Kfast -KSVE -Kparallel -Kopenmp -Kcmodel=large -Nfjomplib -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFLIB) $(NETCDFLIB) $(HDF5LIB) $(STATIC)

#---JSS3
NETCDFINC = -I/opt/JX/oss/aarch64/netcdf-parallel/4.9.2/include
NETCDFFINC = -I/opt/JX/oss/aarch64/netcdf-fortran-parallel/4.6.1/include
NETCDFLIB = -L/opt/JX/oss/aarch64/netcdf-parallel/4.9.2/lib
NETCDFFLIB = -L/opt/JX/oss/aarch64/netcdf-fortran-parallel/4.6.1/lib
HDF5INC = -I/opt/JX/oss/aarch64/phdf5/1.12.0/include
HDF5LIB = -L/opt/JX/oss/aarch64/phdf5/1.12.0/lib

#MPI+openmp/Parallel NETCDF
FFLAGS = -Kfast -KSVE -Kparallel -Kopenmp -Kcmodel=large -Nfjomplib -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
LIBS = -Kfast -KSVE -Kparallel -Kopenmp -Kcmodel=large -Nfjomplib -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFLIB) -lnetcdff ${HDF5LIB} $(NETCDFLIB) -lnetcdf -lnetcdf -lm -lhdf5_hl -lhdf5 -lm -lbz2 -lzstd -lxml2 -lcurl -lfjprofmpi

#-----------------------------------------------------------------------
# Set the executable
#-----------------------------------------------------------------------
BIN = letkf.exe

#-----------------------------------------------------------------------
# Define source directory
#-----------------------------------------------------------------------
SRCDIR = src

#-----------------------------------------------------------------------
# Define objects
#-----------------------------------------------------------------------
OBJS = 	common_setting.o  \
	common_mpi.o \
	common.o \
	common_io.o \
	common_pom.o \
	common_letkf.o  \
	letkf_obs.o \
	letkf_tools.o \
	letkf.o

VPATH = $(SRCDIR)

#-----------------------------------------------------------------------
# Set implicit rules for compilation
#-----------------------------------------------------------------------
%.o: %.f
	@echo
	$(FC) -c $(FFLAGS) $<

%.o: %.f90
	@echo
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90 %.o
	@echo

#-----------------------------------------------------------------------
# Set implicit rules for dependencies
#-----------------------------------------------------------------------
%.f: %.F
	@echo
	$(CPP) $(FFLAGS) $< > $*.f

#-----------------------------------------------------------------------
# Create the executable
#-----------------------------------------------------------------------
$(BIN): $(OBJS)
	@echo
	$(LD) -o $(BIN) $(FFLAGS) $(OBJS) $(LIBS)

#-----------------------------------------------------------------------
# Cleaning target
#-----------------------------------------------------------------------
clean:
#	cd $(SRCDIR)
	@rm -f *.o *.mod
