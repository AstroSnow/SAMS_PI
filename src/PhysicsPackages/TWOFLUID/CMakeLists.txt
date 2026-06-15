cmake_minimum_required(VERSION 3.16)
project(nlevel_T_calc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. Find Boost (Using CONFIG fixes the CMP0167 warning)
find_package(Boost REQUIRED CONFIG)

# 2. Attempt to find NetCDF using modern package config files smoothly
find_package(NetCDF COMPONENTS C QUIET)

# 3. Add the executable
add_executable(nlevel_T_calc nlevel_T_calc.cpp)

# 4. Determine the correct Target or Library variable for NetCDF
if(TARGET NetCDF::NetCDF_C)
    set(NETCDF_LINK_TARGET NetCDF::NetCDF_C)
elseif(TARGET NetCDF::NetCDF)
    set(NETCDF_LINK_TARGET NetCDF::NetCDF)
else()
    # Fallback to the legacy variable if no modern target exists
    set(NETCDF_LINK_TARGET ${NetCDF_LIBRARIES})
endif()

# 4b. Fallback find for library/include if standard variables are empty
if(NOT NetCDF_INCLUDE_DIRS)
    find_path(NetCDF_INCLUDE_DIRS netcdf.h)
endif()
if(NOT NetCDF_LIBRARIES)
    find_library(NetCDF_LIBRARIES netcdf)
endif()

# 5. Apply the settings to your target
target_include_directories(nlevel_T_calc PRIVATE ${NetCDF_INCLUDE_DIRS})
target_link_libraries(nlevel_T_calc PRIVATE Boost::headers)

if(NETCDF_LINK_TARGET)
    target_link_libraries(nlevel_T_calc PRIVATE ${NETCDF_LINK_TARGET})
else()
    target_link_libraries(nlevel_T_calc PRIVATE ${NetCDF_LIBRARIES})
endif()

# --- DEBUGGING OUTPUT ---
message(STATUS "--------------------------------------------")
message(STATUS "Boost Found:             ${Boost_FOUND}")
message(STATUS "Boost Include Dirs:      ${Boost_INCLUDE_DIRS}")
message(STATUS "NetCDF Found:            ${NetCDF_FOUND}")
message(STATUS "NetCDF Include Dirs:     ${NetCDF_INCLUDE_DIRS}")
message(STATUS "NetCDF Libraries:        ${NetCDF_LIBRARIES}")
message(STATUS "NetCDF Link Target:      ${NETCDF_LINK_TARGET}")
message(STATUS "C++ Standard Used:       ${CMAKE_CXX_STANDARD}")
message(STATUS "--------------------------------------------")
