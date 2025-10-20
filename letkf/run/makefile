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
#---Fugaku (As of Dec. 2024)
#NETCDFINC = -I/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-c-4.9.2-cczsjh6lmalmjqhr72yatvxzmc3iwdl2/include
#NETCDFFINC = -I/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-fortran-4.6.1-kjm7jl5naxulm7neknsgxnsbkf7iet5j/include
#NETCDFLIB = -L/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-c-4.9.2-cczsjh6lmalmjqhr72yatvxzmc3iwdl2/lib -lnetcdf
#NETCDFFLIB = -L/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-fortran-4.6.1-kjm7jl5naxulm7neknsgxnsbkf7iet5j/lib -lnetcdff  -lnetcdf -lnetcdf -lm
#HDF5INC = 
#HDF5LIB =
#STATIC = -lhdf5_hl -lhdf5 -lm -lsz -lbz2 -lzstd -lblosc -lmpi_cxx -lfjprofmpi

#FFLAGS = -Kfast -KSVE -Kparallel -Kopenmp -Kcmodel=large -Nfjomplib -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
#LIBS = -Kfast -KSVE -Kparallel -Kopenmp -Kcmodel=large -Nfjomplib -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFLIB) $(NETCDFLIB) $(HDF5LIB) $(STATIC)

#---Fugaku (As of Oct. 2025)
NETCDFINC = -I/vol0004/apps/oss/spack-v1.0.1/opt/spack/linux-a64fx/netcdf-c-4.9.2-gxiznee6h7r35ke6tivigi54gvzefl6v/include
NETCDFFINC = -I/vol0004/apps/oss/spack-v1.0.1/opt/spack/linux-a64fx/netcdf-fortran-4.6.1-4jlnxiqcepafk4y3whu4763jzuyjvxr7/include -I/vol0004/apps/oss/spack-v1.0.1/opt/spack/linux-a64fx/netcdf-fortran-4.6.1-4jlnxiqcepafk4y3whu4763jzuyjvxr7/include 
NETCDFLIB = -L/vol0004/apps/oss/spack-v1.0.1/opt/spack/linux-a64fx/netcdf-c-4.9.2-gxiznee6h7r35ke6tivigi54gvzefl6v/lib -lnetcdf
NETCDFFLIB = -L/vol0004/apps/oss/spack-v1.0.1/opt/spack/linux-a64fx/netcdf-fortran-4.6.1-4jlnxiqcepafk4y3whu4763jzuyjvxr7/lib -lnetcdff  -lnetcdf -lnetcdf -lm
HDF5INC = 
HDF5LIB =
STATIC = -lhdf5_hl -lhdf5 -lm -lsz -lbz2 -lzstd -lblosc -lmpi_cxx -lfjprofmpi

FFLAGS = -Kfast -KSVE -Kparallel -Kopenmp -Kcmodel=large -Nfjomplib -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
LIBS = -Kfast -KSVE -Kparallel -Kopenmp -Kcmodel=large -Nfjomplib -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFLIB) $(NETCDFLIB) $(HDF5LIB) $(STATIC)

#---JSS3
#NETCDFINC = -I/opt/JX/oss/aarch64/netcdf-parallel/4.9.2/include
#NETCDFFINC = -I/opt/JX/oss/aarch64/netcdf-fortran-parallel/4.6.1/include
#NETCDFLIB = -L/opt/JX/oss/aarch64/netcdf-parallel/4.9.2/lib
#NETCDFFLIB = -L/opt/JX/oss/aarch64/netcdf-fortran-parallel/4.6.1/lib
#HDF5INC = -I/opt/JX/oss/aarch64/phdf5/1.12.0/include
#HDF5LIB = -L/opt/JX/oss/aarch64/phdf5/1.12.0/lib

#FFLAGS = -Kfast -KSVE -Kparallel -Kopenmp -Kcmodel=large -Nfjomplib -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
#LIBS = -Kfast -KSVE -Kparallel -Kopenmp -Kcmodel=large -Nfjomplib -Nalloc_assign -SCALAPACK -SSL2 $(NETCDFFLIB) -lnetcdff ${HDF5LIB} $(NETCDFLIB) -lnetcdf -lnetcdf -lm -lhdf5_hl -lhdf5 -lm -lbz2 -lzstd -lxml2 -lcurl -lfjprofmpi

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
	common_mt19937ar.o \
	common.o \
	common_io.o \
	common_pom.o \
	common_letkf.o  \
	common_lpfgm.o  \
	common_obs.o \
	letkf_tools.o \
	lpfgm_tools.o \
	main.o

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
