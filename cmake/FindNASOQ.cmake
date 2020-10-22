# - Try to find NASOQ based on the NASOQ CMake File

if (TARGET NASOQ::NASOQ)
  # already found, do nothing
else()

#  set(CMAKE_CXX_STANDARD 11)
  find_path( NASOQ_BASE_DIR
    NAMES nasoq_driver.cpp
    PATHS $ENV{NASOQ_DIR}
      /usr/local/include
      /usr/local/include/nasoq/
    REQUIRED)

  if (NASOQ_FIND_REQUIRED AND NOT NASOQ_BASE_DIR)
    message(SEND_ERROR "NASOQ_BASE_DIR not found")
  endif()

  message("Jo NASOQ_BASE_DIR: ${NASOQ_BASE_DIR}")

  # TODO: remove MKL dependency using OpenBLAS submodule
  find_package(MKL REQUIRED)

  # Adding hints for suitesparse
  set(SUITESPARSE_INCLUDE_DIR_HINTS ${SUITESPARSE_INCLUDE_DIR_HINTS} ${SUITE_ROOT_PATH}/include)
  set(SUITESPARSE_LIBRARY_DIR_HINTS ${SUITESPARSE_LIBRARY_DIR_HINTS} ${SUITE_ROOT_PATH}/lib)
  set(BLA_STATIC TRUE)
  find_package(SUITESPARSE OPTIONAL_COMPONENTS)

  #TODO adding METIS submodule
  find_package(METIS REQUIRED)

  find_package (Eigen3 REQUIRED NO_MODULE)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(NASOQ  DEFAULT_MSG
                                    MKL_INCLUDE_DIR METIS_INCLUDE_DIR)


  add_library (NASOQ::NASOQ INTERFACE IMPORTED)
  if(UNIX )
    target_compile_definitions(NASOQ::NASOQ INTERFACE "m64")
  endif()

  target_include_directories(NASOQ::NASOQ INTERFACE
    ${MKL_INCLUDE_DIR}
    ${SUITESPARSE_INCLUDE_DIRS}
    ${METIS_INCLUDE_DIRS}
    "${NASOQ_BASE_DIR}/symbolic/"
    "${NASOQ_BASE_DIR}/common/"
    "${NASOQ_BASE_DIR}/ldl/"
    "${NASOQ_BASE_DIR}/matrixMatrix/"
    "${NASOQ_BASE_DIR}/matrixVector/"
    "${NASOQ_BASE_DIR}/linear_solver/"
    "${NASOQ_BASE_DIR}/gmres/"
    "${NASOQ_BASE_DIR}/QP/"
    "${NASOQ_BASE_DIR}/triangularSolve/"
    "${NASOQ_BASE_DIR}/smp-format/"
    "${NASOQ_BASE_DIR}/eigen_interface/")

  target_link_libraries(NASOQ::NASOQ INTERFACE
    ${MKL_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${METIS_LIBRARIES})

endif(TARGET NASOQ::NASOQ)


