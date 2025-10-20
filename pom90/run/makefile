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
#----Fugaku (As of Dec. 2024)
#NETCDFINC = -I/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-c-4.9.2-cczsjh6lmalmjqhr72yatvxzmc3iwdl2/include
#NETCDFFINC = -I/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-fortran-4.6.1-kjm7jl5naxulm7neknsgxnsbkf7iet5j/include
#NETCDFLIB = -L/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-c-4.9.2-cczsjh6lmalmjqhr72yatvxzmc3iwdl2/lib -lnetcdf
#NETCDFFLIB = -L/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-fortran-4.6.1-kjm7jl5naxulm7neknsgxnsbkf7iet5j/lib -lnetcdff  -lnetcdf -lnetcdf -lm
#HDF5INC = 
#HDF5LIB =
#STATIC = -lhdf5_hl -lhdf5 -lm -lsz -lbz2 -lzstd -lblosc -lmpi_cxx -lfjprofmpi

#FFLAGS = -Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign -Nfjomplib $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
#LIBS = -Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign -Nfjomplib $(NETCDFFLIB) $(NETCDFLIB) $(HDF5LIB) $(STATIC)

#---Fugaku (As of Oct. 2025)
NETCDFINC = -I/vol0004/apps/oss/spack-v1.0.1/opt/spack/linux-a64fx/netcdf-c-4.9.2-gxiznee6h7r35ke6tivigi54gvzefl6v/include
NETCDFFINC = -I/vol0004/apps/oss/spack-v1.0.1/opt/spack/linux-a64fx/netcdf-fortran-4.6.1-4jlnxiqcepafk4y3whu4763jzuyjvxr7/include -I/vol0004/apps/oss/spack-v1.0.1/opt/spack/linux-a64fx/netcdf-fortran-4.6.1-4jlnxiqcepafk4y3whu4763jzuyjvxr7/include 
NETCDFLIB = -L/vol0004/apps/oss/spack-v1.0.1/opt/spack/linux-a64fx/netcdf-c-4.9.2-gxiznee6h7r35ke6tivigi54gvzefl6v/lib -lnetcdf
NETCDFFLIB = -L/vol0004/apps/oss/spack-v1.0.1/opt/spack/linux-a64fx/netcdf-fortran-4.6.1-4jlnxiqcepafk4y3whu4763jzuyjvxr7/lib -lnetcdff  -lnetcdf -lnetcdf -lm
HDF5INC = 
HDF5LIB =
STATIC = -lhdf5_hl -lhdf5 -lm -lsz -lbz2 -lzstd -lblosc -lmpi_cxx -lfjprofmpi

FFLAGS = -Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign -Nfjomplib $(NETCDFFINC) $(NETCDFINC) $(HDF5INC)
LIBS = -Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign -Nfjomplib $(NETCDFFLIB) $(NETCDFLIB) $(HDF5LIB) $(STATIC)

#---JSS3
#NETCDFINC = -I/opt/JX/oss/aarch64/netcdf-parallel/4.9.2/include
#NETCDFFINC = -I/opt/JX/oss/aarch64/netcdf-fortran-parallel/4.6.1/include
#NETCDFLIB = -L/opt/JX/oss/aarch64/netcdf-parallel/4.9.2/lib -lnetcdf -lnetcdf -lm
#NETCDFFLIB = -L/opt/JX/oss/aarch64/netcdf-fortran-parallel/4.6.1/lib -lnetcdff
#HDF5INC =  -I/opt/JX/oss/aarch64/phdf5/1.12.0/include
#HDF5LIB = -L/opt/JX/oss/aarch64/phdf5/1.12.0/lib
#STATIC = -lhdf5_hl -lhdf5 -lm -lbz2 -lzstd -lxml2 -lcurl -lfjprofmpi

#FFLAGS = -Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign $(NETCDFFINC) $(HDF5INC) $(NETCDFINC)
#LIBS = -Kfast -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign $(NETCDFFLIB) $(HDF5LIB) $(NETCDFLIB) $(STATIC)

#-Debug
#FFLAGS = -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign $(NETCDFFINC) $(HDF5INC) $(NETCDFINC)
#LIBS = -Kopenmp -Kparallel -Kcmodel=large -Nalloc_assign $(NETCDFFLIB) $(HDF5LIB) $(NETCDFLIB) $(STATIC)

#-----------------------------------------------------------------------
# Executable file
#-----------------------------------------------------------------------

BIN = pom.exe

#-----------------------------------------------------------------------
# Source directory
#-----------------------------------------------------------------------

SRCDIR = src

#-----------------------------------------------------------------------
# Objects
#-----------------------------------------------------------------------

OBJS =	common_pom_var.o	\
	mod_const.o		\
	mod_thermo.o		\
	mod_blk_coare.o		\
	mod_blk_ecmwf.o		\
	mod_blk_ncar.o		\
	mod_aerobulk_compute.o	\
	mod_aerobulk.o		\
	parallel_mpi.o          \
	merge_scatter.o         \
	timelib.o		\
	random.o\
	read_netcdf.o           \
	fluxlib.o		\
	solver.o		\
	turb_MYNNF_lev2.5.o	\
	profq_mynnf.o		\
	bounds_forcing.o	\
	daily_average.o		\
	whole_average.o		\
	write_netcdf.o          \
	advance.o		\
	initialize.o		\
	pom.o

VPATH = $(SRCDIR)

#-----------------------------------------------------------------------
# Rule for compile
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
# Rule for dependencies
#-----------------------------------------------------------------------

%.f: %.F
	@echo
	$(CPP) $(FFLAGS) $< > $*.f

#-----------------------------------------------------------------------
# Create an executable file
#-----------------------------------------------------------------------

$(BIN): $(OBJS)
	@echo
	$(LD) -o $(BIN) $(FFLAGS) $(OBJS) $(LIBS)

#-----------------------------------------------------------------------
# Clean
#-----------------------------------------------------------------------

clean:
	@rm -f *.o *.mod *~ *exe
