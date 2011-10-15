cmake_minimum_required(VERSION 2.6)

set(CMAKE_Fortran_COMPILER_INIT gfortran)

project(3DEX)
enable_language(Fortran)


# CFITSIO
link_directories($ENV{CFITSIO_DIR})

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

set(LIBS "-lcfitsio -lhealpix")
set(3DEXTAG "-l3dex")

# build executables
set(MODULES f3dex.mod" "cosmotools.mod)
set(MODULES_SRC "src/f90/mod/f3dex.f90" "src/f90/mod/cosmotools.f90" "src/f90/mod/rfsfns.f"  "src/f90/mod/quadpack.f90")
set(EXECUTABLES survey2almn survey2almn_lm survey2almn_interactive almnfile2rmap almn2cln)

add_library(3dex ${MODULES_SRC})   

foreach (p ${EXECUTABLES})
  add_executable (${p} "src/f90/${p}/${p}.f90" ${MODULES_SRC})
  target_link_libraries(${p} ${LIBS} ${3DEX})
endforeach (p)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib3dex.a
	DESTINATION ${PROJECT_SOURCE_DIR}/lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${MODULES}
	DESTINATION ${PROJECT_SOURCE_DIR}/include)
foreach (p ${EXECUTABLES})
   install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${p}
	DESTINATION ${PROJECT_SOURCE_DIR}/bin)
endforeach (p)

