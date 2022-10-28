cmake_minimum_required(VERSION 3.18)

# If not set by the user, make a release build.
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()
message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

project(primme
  VERSION 3.2
  DESCRIPTION "PRIMME: PReconditioned Iterative MultiMethod Eigensolver"
  LANGUAGES C CXX
  HOMEPAGE_URL https://github.com/primme/primme)

option(WITH_FLOAT "Build PRIMME with support for floats" OFF)
option(WITH_HALF "Build PRIMME with support for half floats" OFF)
option(WITH_FORTRAN "Build PRIMME Fortran 90 tests" OFF)

# set C99 standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

if (POLICY CMP0042)
  # enable MACOSX_RPATH by default
  cmake_policy (SET CMP0042 NEW)
endif ()

# Find our extra modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(PrimmeTests)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

set(PRIMME_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/primme")
set(PRIMME_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/primme")

if (WITH_HALF)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPRIMME_WITH_HALF")
endif()
if (WITH_FLOAT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPRIMME_WITH_FLOAT")
endif()

find_package(BLAS REQUIRED)
if (NOT TARGET LAPACK::LAPACK)
  find_package(LAPACK REQUIRED)
endif()
primme_test_fortran_underscore()

CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)                                                                                    
if (HAVE_LIB_M)
  set(EXTRA_LIBS ${EXTRA_LIBS} m) 
endif()

add_library(primme)
set(primme_MAJOR_VERSION 3)
set(primme_MINOR_VERSION 9)
set(primme_VERSION ${primme_MAJOR_VERSION}.${primme_MINOR_VERSION})
target_include_directories(primme
  PUBLIC
    $<INSTALL_INTERFACE:${PRIMME_INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/include>
)
target_link_libraries(primme PUBLIC LAPACK::LAPACK BLAS::BLAS ${EXTRA_LIBS})

add_subdirectory(include)
add_subdirectory(src)

enable_testing()
add_subdirectory(examples)

install(TARGETS primme
  EXPORT primmeTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

configure_file(cmake/primme-config.cmake.in "${PROJECT_BINARY_DIR}/primme-config.cmake" @ONLY)
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/primme-config-version.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY SameMajorVersion
  )
install(EXPORT primmeTargets DESTINATION "${PRIMME_INSTALL_CMAKEDIR}")
install(
  FILES
    "${PROJECT_BINARY_DIR}/primme-config.cmake"
    "${PROJECT_BINARY_DIR}/primme-config-version.cmake"
  DESTINATION "${PRIMME_INSTALL_CMAKEDIR}"
)