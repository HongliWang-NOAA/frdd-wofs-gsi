# ===============================================================
# Minimal FindNetCDF.cmake for Intel + MPI + MKL setup
# Assumes you provide:
#   NetCDF_INCLUDE_DIR
#   NetCDF_LIBRARY
#   NetCDF_Fortran_INCLUDE_DIR
#   NetCDF_Fortran_LIBRARY
# ===============================================================

if(NOT DEFINED NetCDF_INCLUDE_DIR OR NOT DEFINED NetCDF_LIBRARY)
  message(FATAL_ERROR "Please set NetCDF_INCLUDE_DIR and NetCDF_LIBRARY")
endif()

if(NOT DEFINED NetCDF_Fortran_INCLUDE_DIR OR NOT DEFINED NetCDF_Fortran_LIBRARY)
  message(FATAL_ERROR "Please set NetCDF_Fortran_INCLUDE_DIR and NetCDF_Fortran_LIBRARY")
endif()

# Include paths
set(NETCDF_INCLUDES "${NetCDF_INCLUDE_DIR}" CACHE PATH "NetCDF include directory")
set(NETCDF_INCLUDES_F90 "${NetCDF_Fortran_INCLUDE_DIR}" CACHE PATH "NetCDF Fortran include directory")

# Library paths
set(NETCDF_LIBRARIES "${NetCDF_LIBRARY}" CACHE FILEPATH "NetCDF C library")
set(NETCDF_LIBRARIES_F90 "${NetCDF_Fortran_LIBRARY}" CACHE FILEPATH "NetCDF Fortran library")

# Set interface flags
set(NETCDF_F90 YES CACHE STRING "Force Fortran 90 interface to YES")
set(NetCDF_has_interfaces YES CACHE STRING "NetCDF interfaces found")

# Message for confirmation
message(STATUS "Using NetCDF C library: ${NETCDF_LIBRARIES}")
message(STATUS "Using NetCDF Fortran library: ${NETCDF_LIBRARIES_F90}")
message(STATUS "Forcing NetCDF F90 interface to YES")
