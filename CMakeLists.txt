cmake_minimum_required(VERSION 2.6)

INCLUDE (CheckLibraryExists)
INCLUDE (FindMPI)
# http://www.cmake.org/Wiki/CMake:How_To_Write_Platform_Checks
# http://www.cmake.org/cmake/help/cmake2.6docs.html

############################################################
# FORTRAN support                                          #
############################################################

set(CMAKE_Fortran_COMPILER_INIT gfortran)

project(3DEX)
enable_language(Fortran)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_LIST_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_LIST_DIR}/lib)

# LAPACK
CHECK_LIBRARY_EXISTS(liblapack.a DGETRF $ENV{LAPACK} foundLapack) 
if (foundLapack)
   link_directories($ENV{LAPACK})
endif (foundLapack)

# CFITSIO
link_directories($ENV{CFITSIO})

# HEALPIX 
include_directories($ENV{HEALPIX}/include)
link_directories($ENV{HEALPIX}/lib)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  set(CMAKE_Fortran_FLAGS "-O3 -DGFORTRAN -fno-second-underscore -frecursive -fopenmp")
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  set(CMAKE_Fortran_FLAGS "-TODO")
elseif (Fortran_COMPILER_NAME STREQUAL "g77") 
  set(CMAKE_Fortran_FLAGS "-TODO")
else (Fortran_COMPILER_NAME STREQUAL "gfortran")
  set(CMAKE_Fortran_FLAGS "-TODO")
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

set(LIBS "-lcfitsio -lhealpix -llapack -lblas -fopenmp")
set(3DEXTAG "-L${LIBRARY_OUTPUT_PATH} -l3dex")

# build executables
set(MODULES "f3dex_stats.mod" "f3dex_cosmotools.mod" "f3dex_fitstools.mod"
"f3dex_utils.mod" "f3dex_transforms.mod")
set(MODULES_SRC "src/f90/mod/f3dex_fitstools.f90"
"src/f90/mod/f3dex_cosmotools.f90" "src/f90/mod/f3dex_utils.f90"
"src/f90/mod/f3dex_stats.f90" "src/f90/mod/f3dex_transforms.f90"
"src/f90/external/rfsfns.f"  "src/f90/external/quadpack.f90"  
"src/f90/external/quadrule.f90")
set(EXECUTABLES survey2almn survey2almn_lm survey2almn_interactive
 almnfile2rmap compute_legen_plm compute_qln compute_qml_pl)

add_library(3dex ${MODULES_SRC})   

foreach (p ${EXECUTABLES})
  add_executable (${p} "src/f90/${p}/${p}.f90")
  target_link_libraries(${p} ${3DEXTAG} ${LIBS})
endforeach (p)

############################################################
# CXX support                                              #
############################################################

find_package(MPI REQUIRED)
if (MPI_FOUND)

  if(APPLE)
    include_directories($ENV{HEALPIX}/src/cxx/osx/include)
    link_directories($ENV{HEALPIX}/src/cxx/osx/lib)
  elseif(UNIX)
    include_directories($ENV{HEALPIX}/src/cxx/generic_gcc/include)
    link_directories($ENV{HEALPIX}/src/cxx/generic_gcc/lib)
  endif(APPLE)

  set(CMAKE_C_COMPILER mpicc)
  set(CMAKE_CXX_COMPILER mpicxx)

  set(CMAKE_CXX_FLAGS "-O3 -fomit-frame-pointer -fno-common -fPIC -fopenmp")
  set(CXX_LIBS "-lhealpix_cxx -lcfitsio -lcxxsupport -lc_utils -lfftpack -lhealpix -lhpxgif -lgfortran")

  add_executable(survey2almn_cxx src/cxx/sbtools.cpp src/cxx/survey.cpp src/cxx/survey2almn.cpp)
  target_link_libraries(survey2almn_cxx ${3DEXTAG} ${CXX_LIBS})

  add_executable(survey2almn_mpi src/cxx/sbtools.cpp src/cxx/survey.cpp src/cxx/survey2almn_mpi.cpp)
  target_link_libraries(survey2almn_mpi ${3DEXTAG} ${CXX_LIBS})

endif (MPI_FOUND)

############################################################
# Test framework                                           #
############################################################

set(TESTEXECUTABLES utest_inverse utest_trace utest_legendre utest_traceprod utest_legendre_angmatrices)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/tests)
foreach (p ${TESTEXECUTABLES})
  add_executable (${p} "src/f90/utests/${p}.f90")
  target_link_libraries(${p} ${3DEXTAG} ${LIBS})
endforeach (p)

enable_testing()

foreach (p ${TESTEXECUTABLES})
  add_test (${p} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${p})
endforeach (p)
