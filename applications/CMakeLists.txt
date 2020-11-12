cmake_minimum_required( VERSION 3.1 )

project(ThermalAnalysis-Applications)

if( EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake" )
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup(TARGETS)
endif()

# these packages don't provide conan_basic_setup() support yet
find_package(UnitConvert)
find_package(BoostUnitDefinitions)

add_subdirectory( .. ThermalAnalysis)

add_executable( tempBuilder tempBuilder.cpp )
target_link_libraries( tempBuilder ThermalAnalysis::ThermalAnalysis CONAN_PKG::boost CONAN_PKG::yaml-cpp UnitConvert::UnitConvert BoostUnitDefinitions::BoostUnitDefinitions )

add_executable( thermal-damage-trender thermal-damage-trender.cpp )
target_link_libraries( thermal-damage-trender ThermalAnalysis::ThermalAnalysis CONAN_PKG::boost UnitConvert::UnitConvert BoostUnitDefinitions::BoostUnitDefinitions )

