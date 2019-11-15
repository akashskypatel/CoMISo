cmake_minimum_required (VERSION 3.9)
if (WIN32)
  cmake_minimum_required (VERSION 3.13)
endif()

project(CoMISo VERSION 1.1.0)

# add our macro directory to cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set (CMAKE_DEBUG_POSTFIX "d")

include (ACGCommon)

# Flag to enable or disable QT widgets in CoMISo
set( WANT_COMISO_QT ON CACHE BOOL "Enable Qt support in CoMISo" )

if( WANT_COMISO_QT )

  acg_qt5 ()

  if( QT5_FOUND)
    add_definitions (-DQT5_FOUND)
    set (COMISO_QT_CONFIG_FILE_SETTINGS "#define COMISO_QT_AVAILABLE 1" )
  else()
    set (COMISO_QT_CONFIG_FILE_SETTINGS "#define COMISO_QT_AVAILABLE 0" )
  endif()

endif( WANT_COMISO_QT )

acg_get_version ()

include (ACGOutput)

set(COMISO_INCLUDE_DIRECTORIES "")
set(COMISO_LINK_DIRECTORIES "")
set(COMISO_LINK_LIBRARIES "")
set(NEED_LAPACK "")

option(COMISO_BUILD_WITHOUT_BOOST "Build CoMISo without boost" OFF)

FIND_PACKAGE( Boost 1.59.0 COMPONENTS system filesystem regex QUIET)
if(Boost_FOUND AND NOT COMISO_BUILD_WITHOUT_BOOST)
  set (COMISO_BOOST_CONFIG_FILE_SETTINGS "#define COMISO_BOOST_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_DIRECTORIES    ${Boost_LIBRARY_DIR} )
  list( APPEND COMISO_LINK_LIBRARIES  ${Boost_LIBRARIES} )

  option(Boost_USE_STATIC_LIBS "Link with static Boost libraries" OFF)
  mark_as_advanced(Boost_USE_STATIC_LIBS)
  if(Boost_USE_STATIC_LIBS)
    set(Boost_USE_STATIC_LIBS ON)
  else()
      # One must add -DBOOST_ALL_DYN_LINK to DEFINITIONS to use Boost
      # auto-link with shared libraries.
    add_definitions( -DBOOST_ALL_DYN_LINK )
  endif()

else()
  set (COMISO_BOOST_CONFIG_FILE_SETTINGS "#define COMISO_BOOST_AVAILABLE 0" )
  if(COMISO_BUILD_WITHOUT_BOOST)
    message (STATUS "Building without boost (COMISO_BUILD_WITHOUT_BOOST flag is set)")
  else()
    message (STATUS "Boost not found!")
  endif()
endif ()


find_package (GMM REQUIRED)
if (GMM_FOUND)
  set (COMISO_GMM_CONFIG_FILE_SETTINGS "#define COMISO_GMM_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${GMM_INCLUDE_DIR} )
  list (APPEND NEED_LAPACK "GMM")
else()
  set (COMISO_GMM_CONFIG_FILE_SETTINGS "#define COMISO_GMM_AVAILABLE 0" )
endif ()

set(TMP_CMAKE_FIND_LIBRARY_PREFIXES "${CMAKE_FIND_LIBRARY_PREFIXES}")
set(CMAKE_FIND_LIBRARY_PREFIXES lib "") #Our blas librariy is called libblas.lib. Is there a better way to find it than this?
find_package (BLAS)
if (BLAS_FOUND )
  set (COMISO_BLAS_CONFIG_FILE_SETTINGS "#define COMISO_BLAS_AVAILABLE 1" )

  list( APPEND COMISO_INCLUDE_DIRECTORIES ${BLAS_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_DIRECTORIES   ${BLAS_LIBRARY_DIRS} ${BLAS_LIBRARY_DIR} )
  list( APPEND COMISO_LINK_LIBRARIES  ${BLAS_LIBRARIES} )
else()
  set (COMISO_BLAS_CONFIG_FILE_SETTINGS "#define COMISO_BLAS_AVAILABLE 0" )
  message (STATUS "BLAS not found!")
endif ()
set(CMAKE_FIND_LIBRARY_PREFIXES "${TMP_CMAKE_FIND_LIBRARY_PREFIXES}")

find_package (ADOLC)
if (ADOLC_FOUND)
  set (COMISO_ADOLC_CONFIG_FILE_SETTINGS "#define COMISO_ADOLC_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${ADOLC_INCLUDE_DIR} )
  list( APPEND COMISO_LINK_DIRECTORIES   ${ADOLC_LIBRARY_DIR} )
  list( APPEND COMISO_LINK_LIBRARIES  ${ADOLC_LIBRARIES} )
else ()
  set (COMISO_ADOLC_CONFIG_FILE_SETTINGS "#define COMISO_ADOLC_AVAILABLE 0" )
  message (STATUS "ADOLC not found!")
endif ()

find_package (SUITESPARSE)
if (SUITESPARSE_FOUND )
  list (APPEND NEED_LAPACK "SUITESPARSE")

  set (COMISO_SUITESPARSE_CONFIG_FILE_SETTINGS "#define COMISO_SUITESPARSE_AVAILABLE 1" )

  list( APPEND COMISO_INCLUDE_DIRECTORIES ${SUITESPARSE_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_DIRECTORIES   ${SUITESPARSE_LIBRARY_DIRS} )
  list( APPEND COMISO_LINK_LIBRARIES  ${SUITESPARSE_LIBRARIES} )

  # special handling, since spqr is incorrect in several distributions
  if(SUITESPARSE_SPQR_VALID)
    set (COMISO_SUITESPARSE_SPQR_CONFIG_FILE_SETTINGS "#define COMISO_SUITESPARSE_SPQR_AVAILABLE 1" )
  else()
    message (STATUS "SUITESPARSE SPQR seems to be invalid!")
    set (COMISO_SUITESPARSE_SPQR_CONFIG_FILE_SETTINGS "#define COMISO_SUITESPARSE_SPQR_AVAILABLE 0" )
  endif()

else ()
  message (STATUS "SUITESPARSE not found!")
  set (COMISO_SUITESPARSE_CONFIG_FILE_SETTINGS "#define COMISO_SUITESPARSE_AVAILABLE 0" )
  set (COMISO_SUITESPARSE_SPQR_CONFIG_FILE_SETTINGS "#define COMISO_SUITESPARSE_SPQR_AVAILABLE 0" )
endif()

if(NOT DISABLE_MPI)
  find_package (MPI)
  if (MPI_FOUND )
    set (COMISO_MPI_CONFIG_FILE_SETTINGS "#define COMISO_MPI_AVAILABLE 1" )
    list( APPEND COMISO_INCLUDE_DIRECTORIES ${MPI_INCLUDE_PATH} )
    list( APPEND COMISO_LINK_LIBRARIES  ${MPI_CXX_LIBRARIES} )
  else ()
    message (STATUS "MPI not found!")
    set (COMISO_MPI_CONFIG_FILE_SETTINGS "#define COMISO_MPI_AVAILABLE 0" )
  endif ()
else ()
  message (STATUS "MPI disabled!")
  set (COMISO_MPI_CONFIG_FILE_SETTINGS "#define COMISO_MPI_AVAILABLE 0" )
endif ()

if(NOT DISABLE_PETSC)
  find_package (PETSC)
  if (PETSC_FOUND AND MPI_FOUND)
    set (COMISO_PETSC_CONFIG_FILE_SETTINGS "#define COMISO_PETSC_AVAILABLE 1" )
    list( APPEND COMISO_INCLUDE_DIRECTORIES ${PETSC_INCLUDE_DIRS} )
    list( APPEND COMISO_LINK_LIBRARIES  ${PETSC_LIBRARY} )
  else ()
    message (STATUS "PETSC or MPI not found!")
    set (COMISO_PETSC_CONFIG_FILE_SETTINGS "#define COMISO_PETSC_AVAILABLE 0" )
  endif ()
else ()
  message (STATUS "PETSC disabled!")
  set (COMISO_PETSC_CONFIG_FILE_SETTINGS "#define COMISO_PETSC_AVAILABLE 0" )
endif()


if(NOT DISABLE_TAO)
  find_package (TAO)
  if (TAO_FOUND AND PETSC_FOUND AND MPI_FOUND)
    set (COMISO_TAO_CONFIG_FILE_SETTINGS "#define COMISO_TAO_AVAILABLE 1" )
    list( APPEND COMISO_INCLUDE_DIRECTORIES ${TAO_INCLUDE_DIRS} )
    list( APPEND COMISO_LINK_DIRECTORIES ${TAO_LIBRARY_DIR} )
    list( APPEND COMISO_LINK_LIBRARIES  ${TAO_LIBRARY} )
  else ()
    message (STATUS "TAO or PETSC or MPI not found!")
    set (COMISO_TAO_CONFIG_FILE_SETTINGS "#define COMISO_TAO_AVAILABLE 0" )
  endif ()
else ()
  message (STATUS "TAO disabled!")
  set (COMISO_TAO_CONFIG_FILE_SETTINGS "#define COMISO_TAO_AVAILABLE 0" )
endif ()

find_package (METIS)
if (METIS_FOUND )
  set (COMISO_METIS_CONFIG_FILE_SETTINGS "#define COMISO_METIS_AVAILABLE 1" )

  list( APPEND COMISO_INCLUDE_DIRECTORIES ${METIS_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_DIRECTORIES   ${METIS_LIBRARY_DIRS} )
  list( APPEND COMISO_LINK_LIBRARIES  ${METIS_LIBRARIES} )
else()
  set (COMISO_METIS_CONFIG_FILE_SETTINGS "#define COMISO_METIS_AVAILABLE 0" )
  message (STATUS "METIS not found!")
endif ()

find_package (MUMPS)
if (MUMPS_FOUND )
  set (COMISO_MUMPS_CONFIG_FILE_SETTINGS "#define COMISO_MUMPS_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${MUMPS_INCLUDE_DIR} )
  list( APPEND COMISO_LINK_LIBRARIES ${MUMPS_LIBRARY} )
else ()
  message (STATUS "MUMPS not found!")
  set (COMISO_MUMPS_CONFIG_FILE_SETTINGS "#define COMISO_MUMPS_AVAILABLE 0" )
endif ()

find_package (IPOPT)
if (IPOPT_FOUND)
  set (COMISO_IPOPT_CONFIG_FILE_SETTINGS "#define COMISO_IPOPT_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${IPOPT_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_DIRECTORIES ${IPOPT_LIBRARY_DIRS} )
  list( APPEND COMISO_LINK_LIBRARIES ${IPOPT_LIBRARIES} )
  if ( IPOPT_HSL_LIBRARY_DIR )
    set (COMISO_HSL_CONFIG_FILE_SETTINGS "#define COMISO_HSL_AVAILABLE 1" )
  else ()
    set (COMISO_HSL_CONFIG_FILE_SETTINGS "#define COMISO_HSL_AVAILABLE 0" )
  endif()
else ()
  message (STATUS "IPOPT or dependency not found!")
  set (COMISO_IPOPT_CONFIG_FILE_SETTINGS "#define COMISO_IPOPT_AVAILABLE 0" )
  set (COMISO_HSL_CONFIG_FILE_SETTINGS "#define COMISO_HSL_AVAILABLE 0" )
endif ()

find_package (EIGEN3)
if (EIGEN3_FOUND )
  set (COMISO_EIGEN3_CONFIG_FILE_SETTINGS "#define COMISO_EIGEN3_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIR} )
else ()
  message (FATAL_ERROR "EIGEN3 not found!")
  set (COMISO_EIGEN3_CONFIG_FILE_SETTINGS "#define COMISO_EIGEN3_AVAILABLE 0" )
endif ()

find_package (DCO)
if (DCO_FOUND )
  set (COMISO_DCO_CONFIG_FILE_SETTINGS "#define COMISO_DCO_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${DCO_INCLUDE_DIR} )
else ()
  message (STATUS "DCO not found!")
  set (COMISO_DCO_CONFIG_FILE_SETTINGS "#define COMISO_DCO_AVAILABLE 0" )
endif ()

find_package (CBC)
if (CBC_FOUND )
  set (COMISO_CBC_CONFIG_FILE_SETTINGS "#define COMISO_CBC_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${CBC_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_LIBRARIES  ${CBC_LIBRARIES} )
else ()
  message (STATUS "CBC not found!")
  set (COMISO_CBC_CONFIG_FILE_SETTINGS "#define COMISO_CBC_AVAILABLE 0" )
endif ()

find_package (CLP)
if (CLP_FOUND )
  set (COMISO_CLP_CONFIG_FILE_SETTINGS "#define COMISO_CLP_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${CLP_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_LIBRARIES  ${CLP_LIBRARIES} )
else ()
  message (STATUS "CLP not found!")
  set (COMISO_CLP_CONFIG_FILE_SETTINGS "#define COMISO_CLP_AVAILABLE 0" )
endif ()

find_package (CGL)
if (CGL_FOUND )
  set (COMISO_CGL_CONFIG_FILE_SETTINGS "#define COMISO_CGL_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${CGL_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_LIBRARIES  ${CGL_LIBRARIES} )
else ()
  message (STATUS "CGL not found!")
  set (COMISO_CGL_CONFIG_FILE_SETTINGS "#define COMISO_CGL_AVAILABLE 0" )
endif ()

find_package (COINUTILS)
if (COINUTILS_FOUND )
  set (COMISO_COINUTILS_CONFIG_FILE_SETTINGS "#define COMISO_COINUTILS_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${COINUTILS_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_LIBRARIES  ${COINUTILS_LIBRARIES} )
else ()
  message (STATUS "COINUTILS not found!")
  set (COMISO_COINUTILS_CONFIG_FILE_SETTINGS "#define COMISO_COINUTILS_AVAILABLE 0" )
endif ()

find_package (OSI)
if (OSI_FOUND )
  set (COMISO_OSI_CONFIG_FILE_SETTINGS "#define COMISO_OSI_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${OSI_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_LIBRARIES  ${OSI_LIBRARIES} )
else ()
  message (STATUS "OSI not found!")
  set (COMISO_OSI_CONFIG_FILE_SETTINGS "#define COMISO_OSI_AVAILABLE 0" )
endif ()

find_package (Taucs QUIET)
set (SUPRESS_TAUCS OFF CACHE BOOL "Set to ON if CoMISo should be built without TAUCS.")
if (TAUCS_FOUND AND NOT SUPRESS_TAUCS)
  list (APPEND NEED_LAPACK "Taucs")
  set (COMISO_TAUCS_CONFIG_FILE_SETTINGS "#define COMISO_TAUCS_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${TAUCS_INCLUDE_DIR} )
  list( APPEND COMISO_LINK_LIBRARIES  ${TAUCS_LIBRARY} )

else ()
  message (STATUS "TAUCS not found!")
  set (COMISO_TAUCS_CONFIG_FILE_SETTINGS "#define COMISO_TAUCS_AVAILABLE 0" )
endif ()

find_package (GUROBI)
if (GUROBI_FOUND )
  set (COMISO_GUROBI_CONFIG_FILE_SETTINGS "#define COMISO_GUROBI_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${GUROBI_INCLUDE_DIRS} )
#  list( APPEND COMISO_LINK_DIRECTORIES    ${GUROBI_LIBRARY_DIR} )
  list( APPEND COMISO_LINK_LIBRARIES  ${GUROBI_LIBRARIES} )
else ()
  message (STATUS "GUROBI not found!")
  set (COMISO_GUROBI_CONFIG_FILE_SETTINGS "#define COMISO_GUROBI_AVAILABLE 0" )
endif ()

find_package (MOSEK)
if (MOSEK_FOUND )
  set (COMISO_MOSEK_CONFIG_FILE_SETTINGS "#define COMISO_MOSEK_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${MOSEK_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_DIRECTORIES    ${MOSEK_LIBRARY_DIR} )
  list( APPEND COMISO_LINK_LIBRARIES  ${MOSEK_LIBRARIES} )
else ()
  message (STATUS "MOSEK not found!")
  set (COMISO_MOSEK_CONFIG_FILE_SETTINGS "#define COMISO_MOSEK_AVAILABLE 0" )
endif ()


find_package (ARPACK)
if (ARPACK_FOUND )
  set (COMISO_ARPACK_CONFIG_FILE_SETTINGS "#define COMISO_ARPACK_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${ARPACK_INCLUDE_DIR} )
#  list( APPEND COMISO_LINK_DIRECTORIES    ${ARPACK_LIBRARY_DIR} )
  list( APPEND COMISO_LINK_LIBRARIES  ${ARPACK_LIBRARY} )
else ()
  message (STATUS "ARPACK not found!")
  set (COMISO_ARPACK_CONFIG_FILE_SETTINGS "#define COMISO_ARPACK_AVAILABLE 0" )
endif ()

find_package (CPLEX)
if (CPLEX_FOUND )
  set (COMISO_CPLEX_CONFIG_FILE_SETTINGS "#define COMISO_CPLEX_AVAILABLE 1" )
  list( APPEND COMISO_INCLUDE_DIRECTORIES ${CPLEX_INCLUDE_DIRS} )
  list( APPEND COMISO_LINK_LIBRARIES  ${CPLEX_LIBRARIES} )

  #enable c++ support
   add_definitions(-DIL_STD)
else ()
  message (STATUS "CPLEX not found!")
  set (COMISO_CPLEX_CONFIG_FILE_SETTINGS "#define COMISO_CPLEX_AVAILABLE 0" )
endif ()

set(TMP_CMAKE_FIND_LIBRARY_PREFIXES "${CMAKE_FIND_LIBRARY_PREFIXES}")
set(CMAKE_FIND_LIBRARY_PREFIXES lib "") #Our lapack library is called liblapack.lib. Is there a better way to find it than this?
if (NEED_LAPACK AND NOT SUITESPARSE_FOUND)
  FIND_PACKAGE(LAPACK)
  if (LAPACK_FOUND)
    list( APPEND COMISO_INCLUDE_DIRECTORIES ${LAPACK_INCLUDE_DIR} )
    list( APPEND COMISO_LINK_DIRECTORIES    ${LAPACK_LIBRARY_DIR} )
    list( APPEND COMISO_LINK_LIBRARIES  ${LAPACK_LIBRARIES} )
  else()
    message(STATUS "LAPACK not found!")
  endif()
endif(NEED_LAPACK AND NOT SUITESPARSE_FOUND)
set(CMAKE_FIND_LIBRARY_PREFIXES "${TMP_CMAKE_FIND_LIBRARY_PREFIXES}")

# source code directories
set (directories
  .
  Base/Code
  Base/Debug
  Base/Security
  Base/Utils
  Solver
  NSolver
  EigenSolver
  Config
  Utils
  QtWidgets
)

# collect all header,source and ui files
acg_append_files (headers "*.hh" ${directories})
acg_append_files (sources "*.cc" ${directories})

if( QT5_FOUND )
acg_append_files (ui "*.ui" ${directories})
endif()

macro (of_list_filter _list)
  if (WIN32)
   foreach (_element ${${_list}})
     if (_element MATCHES "gnuplot_i\\.(cc|hh)$")
       list (REMOVE_ITEM ${_list} ${_element})
     endif ()
   endforeach ()
  endif ()
endmacro ()

of_list_filter  ( headers )
of_list_filter  ( sources )

# remove template cc files from source file list
acg_drop_templates (sources)

if ( QT5_FOUND )
  # Enable automoc
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)

  list (APPEND COMISO_LINK_LIBRARIES ${QT_LIBRARIES})
endif()

# TODO: remove this line when we can assume people have a cmake 3.13.
# Then link directories is set as target property below.
link_directories(${COMISO_LINK_DIRECTORIES})

acg_add_library (CoMISo SHARED ${ui} ${sources} ${headers} )
add_library(CoMISo::CoMISo ALIAS CoMISo)

target_include_directories(CoMISo
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> # TODO: this is dangerous!
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/>
)

target_include_directories(CoMISo PUBLIC "${COMISO_INCLUDE_DIRECTORIES}")

#get_target_property(INC_DIRS CoMISo INCLUDE_DIRECTORIES)
#get_target_property(INT_INC_DIRS CoMISo INTERFACE_INCLUDE_DIRECTORIES)
#message("Target CoMISo include dirs: ${INC_DIRS}")
#message("Target CoMISo int include dirs: ${INT_INC_DIRS}")
#message("COMISO_INCLUDE_DIRECTORIES: ${COMISO_INCLUDE_DIRECTORIES}")
#message("COMISO_LINK_DIRECTORIES: ${COMISO_LINK_DIRECTORIES}")

if (COMISO_LINK_DIRECTORIES)
  set_target_properties(CoMISo PROPERTIES LINK_DIRECTORIES "${COMISO_LINK_DIRECTORIES}")
endif()
if (WIN32)
#TODO: do this also for non windows platforms as soon as we can require cmake 3.13
  target_link_directories(CoMISo PUBLIC ${COMISO_LINK_DIRECTORIES})
endif()

# generate dllexport macros on windows
if (WIN32)
  add_definitions(-DCOMISODLL)
  add_definitions(-DBASEDLL)
  add_definitions(-D_SCL_SECURE_NO_DEPRECATE)
  # ignore "... needs to have dll-interface to be used by clients"
  #add_definitions(/wd4251)
  target_compile_definitions(CoMISo PRIVATE /wd4251)
endif ()

# enable debug stuff
set (COMISO_ENABLE_DEBUG_OUTPUT TRUE CACHE BOOL "Enable CoMISo Debug Output")
if (COMISO_ENABLE_DEBUG_OUTPUT)
  #add_definitions(-DDEB_ON)
  target_compile_definitions(CoMISo PRIVATE -DDEB_ON)
endif ()

if (WIN32)
  #add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX)
  target_compile_definitions(CoMISo PRIVATE -D_USE_MATH_DEFINES -DNOMINMAX)
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "[cC][lL][aA][nN][gG]")
  # disable unused parameter warning
  #add_definitions(-Wno-unused-parameter -Wno-deprecated-register)
  target_compile_options(CoMISo PRIVATE -Wno-unused-parameter -Wno-deprecated-register)
endif()
if ( QT5_FOUND )

  if ( NOT WIN32 )
    # Allow targets depending on ACG find the autogenerated header files
    #target_include_directories( CoMISo PUBLIC
    #                           ${CMAKE_CURRENT_BINARY_DIR}/CoMISo_autogen/include )
  else()
     # Allow targets depending on ACG find the autogenerated header files
     #target_include_directories( CoMISo PUBLIC
     #                          $<$<CONFIG:Debug>:${CMAKE_CURRENT_BINARY_DIR}/CoMISo_autogen/include_Debug/>
     #                          $<$<CONFIG:Release>:${CMAKE_CURRENT_BINARY_DIR}/CoMISo_autogen/include_Release/> )
  endif()

endif()

target_link_libraries (CoMISo
  ${COMISO_LINK_LIBRARIES}
)


# display results
acg_print_configure_header (COMISO "CoMISo")

# write config file
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/Config/config.hh.in"
     "CoMISo/Config/config.hh" @ONLY IMMEDIATE)



#######################################################################
# Configure the examples last to be sure, that all configure files
# of the library are already there
#######################################################################

set (COMISO_BUILD_EXAMPLES TRUE CACHE BOOL "Build CoMISo Examples")

if (COMISO_BUILD_EXAMPLES )
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/factored_solver/CMakeLists.txt" )
    add_subdirectory (Examples/factored_solver)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/quadratic_solver/CMakeLists.txt" )
    add_subdirectory (Examples/quadratic_solver)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/test2/CMakeLists.txt" )
    add_subdirectory (Examples/test2)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_quadratic_example/CMakeLists.txt" )
    add_subdirectory (Examples/small_quadratic_example)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_factored_example/CMakeLists.txt" )
    add_subdirectory (Examples/small_factored_example)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/super_sparse_matrix/CMakeLists.txt" )
    add_subdirectory (Examples/super_sparse_matrix)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/eigen_solver/CMakeLists.txt" )
    add_subdirectory (Examples/eigen_solver)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_nsolver/CMakeLists.txt" )
    add_subdirectory (Examples/small_nsolver)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_eigenproblem/CMakeLists.txt" )
    add_subdirectory (Examples/small_eigenproblem)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_miqp/CMakeLists.txt" )
    add_subdirectory (Examples/small_miqp)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_nleast_squares/CMakeLists.txt" )
    add_subdirectory (Examples/small_nleast_squares)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_sparseqr/CMakeLists.txt" )
    add_subdirectory (Examples/small_sparseqr)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_quadratic_resolve_example/CMakeLists.txt" )
    add_subdirectory (Examples/small_quadratic_resolve_example)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_cplex_soc/CMakeLists.txt" )
    add_subdirectory (Examples/small_cplex_soc)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_adolc/CMakeLists.txt" )
    add_subdirectory (Examples/small_adolc)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_dco/CMakeLists.txt" )
    add_subdirectory (Examples/small_dco)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/vector1_adolc/CMakeLists.txt" )
    add_subdirectory (Examples/vector1_adolc)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_linear_problem/CMakeLists.txt" )
    add_subdirectory (Examples/small_linear_problem)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/crossfield3d/CMakeLists.txt" )
    add_subdirectory (Examples/crossfield3d)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/crossfield3d/CMakeLists.txt" )
    add_subdirectory (Examples/crossfield3d_dco)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_finite_element/CMakeLists.txt" )
    add_subdirectory (Examples/small_finite_element)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_AQP/CMakeLists.txt" )
    add_subdirectory (Examples/small_AQP)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/finite_element_integrability_problem/CMakeLists.txt" )
    add_subdirectory (Examples/finite_element_integrability_problem)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_mosek_native/CMakeLists.txt" )
    add_subdirectory (Examples/small_mosek_native)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_mosek_fusion_sdp/CMakeLists.txt" )
    add_subdirectory (Examples/small_mosek_fusion_sdp)
  endif()
  if( EXISTS "${CMAKE_SOURCE_DIR}/Examples/small_symmetric_dirichlet/CMakeLists.txt" )
  add_subdirectory (Examples/small_symmetric_dirichlet)
  endif()

endif (COMISO_BUILD_EXAMPLES )

# Only create install target, when we are building CoMISo standalone
if(NOT COMISO_NO_INSTALL)

  include(CMakePackageConfigHelpers)
  include(GNUInstallDirs)
  
  target_include_directories(CoMISo
    PUBLIC
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/CoMISo>
  )
  
  set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/CoMISo)

  write_basic_package_version_file(
      CoMISoConfigVersion.cmake
      VERSION ${PROJECT_VERSION}
      COMPATIBILITY SameMajorVersion)

  install(FILES
      "${CMAKE_CURRENT_BINARY_DIR}/CoMISoConfigVersion.cmake"
      DESTINATION ${INSTALL_CONFIGDIR})

  # Install library
  install(TARGETS CoMISo EXPORT CoMISoTargets
          RUNTIME DESTINATION lib
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)

  # Install Header Files
  install(DIRECTORY .
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/CoMISo # include/CoMISo
          FILES_MATCHING
          PATTERN "*.hh"
          PATTERN "*T.cc"
          PATTERN "CVS" EXCLUDE
          PATTERN ".svn" EXCLUDE
          PATTERN ".git" EXCLUDE
          PATTERN "Examples" EXCLUDE
          PATTERN "CI" EXCLUDE
          PATTERN "tmp" EXCLUDE
          PATTERN "Templates" EXCLUDE
          PATTERN "Debian*" EXCLUDE)

  # Install Config File
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CoMISo/Config/config.hh
          DESTINATION include/CoMISo/Config)

  # Create export target to reference CoMISo installation
  install(EXPORT CoMISoTargets
          FILE CoMISoConfig.cmake
          NAMESPACE CoMISo::
          DESTINATION ${INSTALL_CONFIGDIR})

  export(EXPORT CoMISoTargets
         NAMESPACE CoMISo::
         FILE CoMISoConfig.cmake)
  export(PACKAGE CoMISo)
endif()

list(REMOVE_AT CMAKE_MODULE_PATH -1)
