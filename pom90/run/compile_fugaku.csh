#!/bin/csh

setenv SPACK_ROOT /vol0004/apps/oss/spack
source /vol0004/apps/oss/spack/share/spack/setup-env.csh
spack load ${netcdf_fj}

make
