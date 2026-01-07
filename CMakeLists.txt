cmake_minimum_required (VERSION 3.18)

project(CoMISo VERSION 1.1.0 LANGUAGES C CXX)

if(POLICY CMP0100)
  cmake_policy(SET CMP0100 NEW)
endif()

if(POLICY CMP0167) # boost finder removal
    cmake_policy(SET CMP0167 NEW)
endif()

if(CMAKE_VERSION VERSION_LESS 3.21.0)
  if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(PROJECT_IS_TOP_LEVEL ON)
  else()
    set(PROJECT_IS_TOP_LEVEL OFF)
  endif()
endif()

# add our macro directory to cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake-library/finders
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake-library/VCI)

set (CMAKE_DEBUG_POSTFIX "d")

include(FindPkgConfig)
include(VCICommon)
include(VCIQt)
include(VCIOutput)


if(PROJECT_IS_TOP_LEVEL)
# Select appropriate C++ version.
    set(COMISO_CXX_STANDARD 17 CACHE STRING "C++ standard version to use")
    set_property(CACHE COMISO_CXX_STANDARD PROPERTY STRINGS 11 14 17 20 23)
    set(CMAKE_CXX_STANDARD ${COMISO_CXX_STANDARD})
    set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
    set(CMAKE_CXX_EXTENSIONS FALSE)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/bin")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
    set(CMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY ON) # deprecated since CMake 3.16
    set(CMAKE_FIND_USE_PACKAGE_REGISTRY OFF) # since CMake 3.16
endif()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build shared instead of static libraries" ON)
endif()


option(COMISO_ENABLE_DEFAULT "Default value for COMISO_ENABLE_* - only has an effect on first cmake run. Not recommended, default will change to OFF in the future." ON)

set(COMISO_ENABLE_GMM     ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with GMM (if found)" )
set(COMISO_ENABLE_ARPACK  ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with ARPACK (if found)" )
set(COMISO_ENABLE_BOOST   ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with Boost (if found)" )
set(COMISO_ENABLE_CBC     ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with CBC (if found)" )
set(COMISO_ENABLE_CGL     ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with CGL (if found)" )
set(COMISO_ENABLE_CLP     ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with CLP (if found)" )
set(COMISO_ENABLE_CPLEX   ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with CPLEX (if found)" )
set(COMISO_ENABLE_DCO     ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with DCO (if found)" )
set(COMISO_ENABLE_GUROBI  ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with Gurobi (if found)" )
set(COMISO_ENABLE_IPOPT   ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with IPOPT (if found)" )
set(COMISO_ENABLE_LPSOLVE ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with lpsolve (if found)" )
set(COMISO_ENABLE_METIS   ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with METIS (if found)" )
set(COMISO_ENABLE_MOSEK   ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with Mosek (if found)" )
set(COMISO_ENABLE_MPI     ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with MPI (if found). Only used for TAO and PETSc." )
set(COMISO_ENABLE_MUMPS   ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with MUMPS (if found)" )
set(COMISO_ENABLE_OSI     ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with OSI (if found)" )
set(COMISO_ENABLE_OSQP    ${COMISO_ENABLE_DEFAULT}   CACHE BOOL "Build CoMISo with OSQP (if found)" )
set(COMISO_ENABLE_PETSC   ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with PETSc (if found)." )
set(COMISO_ENABLE_QT      ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with Qt widgets (if found)")
set(COMISO_ENABLE_SUITESPARSE_DEFAULT   ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Default for COMISO_ENABLE_SUITESPARSE_* - only has an effect on first cmake run!" )
set(COMISO_ENABLE_TAO     ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with TAO (if found)." )
set(COMISO_ENABLE_TAUCS   ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with TAUCS (if found)" )
set(COMISO_ENABLE_TINYAD  ${COMISO_ENABLE_DEFAULT}  CACHE BOOL "Build CoMISo with TinyAD (if found)" )

set(COMISO_OSQP_PRERELEASE_API FALSE CACHE BOOL "Use pre-1.0 API of OSQP (DEPRECATED)" )

set (COMISO_ENABLE_DEBUG_OUTPUT TRUE CACHE BOOL "Enable CoMISo Debug Output")

### Deprecated aliases:
if(DEFINED WANT_COMISO_QT)
    message(WARNING "CoMISo: Legacy option WANT_COMISO_QT was passed, use COMISO_ENABLE_QT instead.")
    set(COMISO_ENABLE_QT ${WANT_COMISO_QT})
endif()

if (DISABLE_MPI)
    message(WARNING "CoMISo: Legacy option DISABLE_MPI was passed, use COMISO_ENABLE_MPI=FALSE instead.")
    set(COMISO_ENABLE_MPI FALSE)
endif()

if (DISABLE_PETSC)
    message(WARNING "CoMISo: Legacy option DISABLE_PETSC was passed, use COMISO_ENABLE_PETSC=FALSE instead.")
    set(COMISO_ENABLE_PETSC FALSE)
endif()

if (DISABLE_TAO)
    message(WARNING "CoMISo: Legacy option DISABLE_TAO was passed, use COMISO_ENABLE_TAO=FALSE instead.")
    set(COMISO_ENABLE_TAO FALSE)
endif()

if(DEFINED ENABLE_IPOPT)
    message(WARNING "CoMISo: Legacy option ENABLE_IPOPT was passed, use COMISO_ENABLE_IPOPT instead.")
    set(COMISO_ENABLE_IPOPT ${ENABLE_IPOPT})
endif()

if(DEFINED COMISO_BUILD_WITHOUT_BOOST AND COMISO_BUILD_WITHOUT_BOOST)
    message(WARNING "CoMISo: Legacy option COMISO_BUILD_WITHOUT_BOOST was passed, use COMISO_ENABLE_BOOST=FALSE instead.")
    set(COMISO_ENABLE_BOOST OFF)
endif()


add_library(CoMISo
    ./Base/Debug/DebCallStack.cc
    ./Base/Debug/DebConfig.cc
    ./Base/Debug/DebFile.cc
    ./Base/Debug/DebIndexMeshOut.cc
    ./Base/Debug/DebStream.cc
    ./Base/Utils/BaseError.cc
    ./Base/Utils/Environment.cc
    ./Base/Utils/FileOutput.cc
    ./Base/Utils/IOutputStream.cc
    ./Base/Utils/NullOutputStream.cc
    ./Base/Utils/OStringStream.cc
    ./Base/Utils/Option.cc
    ./Base/Utils/RedirectStream.cc
    ./Base/Utils/StopWatch.cc
    ./EigenSolver/ArpackSolver.cc
    ./NSolver/BoundConstraint.cc
    ./NSolver/CBCSolver.cc
    ./NSolver/COMISOSolver.cc
    ./NSolver/CPLEXSolver.cc
    ./NSolver/CombinedProblem.cc
    ./NSolver/ConeConstraint.cc
    ./NSolver/ConstraintTools.cc
    #./NSolver/DOCloudCache.cc
    #./NSolver/DOCloudJob.cc
    #./NSolver/DOCloudSolver.cc
    ./NSolver/FiniteElementProblem.cc
    ./NSolver/GUROBISolver.cc
    ./NSolver/GurobiHelper.cc
    ./NSolver/IPOPTProblemInstance.cc
    ./NSolver/IPOPTSolver.cc
    ./NSolver/LPSolveSolver.cc
    ./NSolver/LeastSquaresProblem.cc
    ./NSolver/LinearConstraint.cc
    ./NSolver/LinearConstraintHandlerElimination.cc
    ./NSolver/LinearConstraintHandlerPenalty.cc
    ./NSolver/LinearProblem.cc
    ./NSolver/NPDerivativeChecker.cc
    ./NSolver/NPLinearConstraints.cc
    ./NSolver/NPTiming.cc
    ./NSolver/NProblemInterface.cc
    ./NSolver/NewtonSolver.cc
    ./NSolver/OSQPSolver.cc
    ./NSolver/SymmetricDirichletProblem.cc
    ./NSolver/TAOSolver.cc
    ./NSolver/TapeIDSingleton.cc
    ./NSolver/TruncatedNewtonPCG.cc
    ./NSolver/cURLpp.cc
    ./Solver/CholmodSolver.cc
    ./Solver/ConstrainedSolver.cc
    ./Solver/EigenLDLTSolver.cc
    ./Solver/EigenLSQSolver.cc
    ./Solver/GMM_Tools.cc
    ./Solver/IterativeSolver.cc
    ./Solver/MISolver.cc
    ./Solver/MultiDimConstrainedSolver.cc
    ./Solver/SparseQRSolver.cc
    ./Solver/SubProblemSolver.cc
    ./Solver/TaucsSolver.cc
    ./Solver/UMFPACKSolver.cc
    ./Utils/CoMISoError.cc
    ./Utils/ExactConstraintSatisfaction.cc
    ./Utils/ExactConstraintProjection.cc
    ./Utils/MatrixDecomposition.cc
)


#target_sources(CoMISo PRIVATE ${sources})

set_target_properties (CoMISo PROPERTIES
    VERSION ${CoMISo_VERSION_MAJOR}.${CoMISo_VERSION_MINOR}
    SOVERSION ${CoMISo_VERSION_MAJOR}.${CoMISo_VERSION_MINOR}
    POSITION_INDEPENDENT_CODE ON
    AUTOMOC OFF
    AUTOUIC OFF
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN TRUE
    )

if (COMISO_ENABLE_DEBUG_OUTPUT)
  target_compile_definitions(CoMISo PUBLIC -DDEB_ON)
endif ()

add_library(CoMISo::CoMISo ALIAS CoMISo)



target_include_directories(CoMISo
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> # TODO: this is dangerous!
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/CoMISo>
)

include(GenerateExportHeader)

option(COMISO_EXCLUDE_DEPRECATED "Exclude deprecated parts of CoMISo" FALSE)
if (COMISO_EXCLUDE_DEPRECATED)
  set(NO_BUILD_DEPRECATED DEFINE_NO_DEPRECATED)
else()
  set(NO_BUILD_DEPRECATED "")
endif()


generate_export_header(CoMISo
    BASE_NAME COMISO
    EXPORT_FILE_NAME "src/CoMISo/Config/Export.hh"
    ${NO_BUILD_DEPRECATED}
)

generate_export_header(CoMISo
    BASE_NAME BASE
    EXPORT_FILE_NAME "src/Base/Config/Export.hh"
    ${NO_BUILD_DEPRECATED}
)


# Maybe CoMISo was added using add_subdirectory and Eigen is already available?
if (NOT TARGET Eigen3::Eigen)
    find_package (Eigen3 REQUIRED CONFIG)
endif()
target_link_libraries     (CoMISo PUBLIC Eigen3::Eigen)
set (COMISO_EIGEN3_CONFIG_FILE_SETTINGS "#define COMISO_EIGEN3_AVAILABLE 1" )

if( COMISO_ENABLE_QT )
  set(QT_REQUIRED_PACKAGES Core Widgets)
  vci_qt ()

  if( QT_FOUND)
    target_sources(CoMISo PRIVATE ./QtWidgets/MISolverDialogUI.cc)
    set_target_properties (CoMISo PROPERTIES
        AUTOMOC ON
        AUTOUIC ON)
    target_link_libraries(CoMISo PUBLIC ${QT_LIBRARIES})
    set (COMISO_QT_CONFIG_FILE_SETTINGS "#define COMISO_QT_AVAILABLE 1" )
  else()
    set (COMISO_QT_CONFIG_FILE_SETTINGS "#define COMISO_QT_AVAILABLE 0" )
  endif()
endif( COMISO_ENABLE_QT )



set(NEED_LAPACK "")

find_package(OpenMP)
if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(CoMISo PUBLIC OpenMP::OpenMP_CXX)
endif()


if(COMISO_ENABLE_BOOST)
    find_package( Boost 1.59.0 COMPONENTS system regex QUIET MODULE)
    if(Boost_FOUND)
        set (COMISO_BOOST_CONFIG_FILE_SETTINGS "#define COMISO_BOOST_AVAILABLE 1" )
        target_link_libraries     (CoMISo PRIVATE ${Boost_LIBRARIES})

        # TODO: is this something a consuming project should set instead?
        option(Boost_USE_STATIC_LIBS "Link with static Boost libraries" OFF)
        mark_as_advanced(Boost_USE_STATIC_LIBS)
        if(Boost_USE_STATIC_LIBS)
            set(Boost_USE_STATIC_LIBS ON)
        else()
            # One must add -DBOOST_ALL_DYN_LINK to DEFINITIONS to use Boost
            # auto-link with shared libraries.
            add_definitions( -DBOOST_ALL_DYN_LINK )
        endif()
        if (APPLE)
            # cf https://github.com/boostorg/boost_install/issues/47#issuecomment-1548879217
            # if we can avoid regex, we don't need this:
            find_package(ICU COMPONENTS uc i18n data REQUIRED)
            add_library(icuuc ALIAS ICU::uc)
            add_library(icui18n ALIAS ICU::i18n)
            add_library(icudata ALIAS ICU::data)
        endif()
    else()
        set (COMISO_BOOST_CONFIG_FILE_SETTINGS "#define COMISO_BOOST_AVAILABLE 0" )
        message (STATUS "CoMISo: optional dependency Boost not found.")
    endif()
else()
    set (COMISO_BOOST_CONFIG_FILE_SETTINGS "#define COMISO_BOOST_AVAILABLE 0" )
    message (STATUS "CoMISo: Building without boost (COMISO_ENABLE_BOOST=OFF)")
endif()


if(COMISO_ENABLE_GMM)
    if (TARGET GMM::GMM)
        message (STATUS "CoMISo: Building with supplied GMM target")
    else()
        find_package(GMM)
        if (TARGET GMM::GMM)
            message (STATUS "CoMISo: Building with GMM found via find_package()")
        endif()
    endif()
    if(TARGET GMM::GMM)
        set (COMISO_GMM_CONFIG_FILE_SETTINGS "#define COMISO_GMM_AVAILABLE 1" )
        set(COMISO_GMM_AVAILABLE 1)
        list (APPEND NEED_LAPACK "GMM")
        target_link_libraries     (CoMISo PUBLIC GMM::GMM)
    else()
        set(COMISO_GMM_AVAILABLE 0)
        set (COMISO_GMM_CONFIG_FILE_SETTINGS "#define COMISO_GMM_AVAILABLE 0" )
        message (WARNING "CoMISo: COMISO_ENABLE_GMM set, but GMM not found.")
    endif()
else ()
    set(COMISO_GMM_AVAILABLE 0)
    set (COMISO_GMM_CONFIG_FILE_SETTINGS "#define COMISO_GMM_AVAILABLE 0" )
    message (STATUS "CoMISo: Building without GMM (set COMISO_ENABLE_GMM to enable)")
endif()

find_package (BLAS)
if (BLAS_FOUND )
  set (COMISO_BLAS_CONFIG_FILE_SETTINGS "#define COMISO_BLAS_AVAILABLE 1" )

  message (STATUS "CoMISo: BLAS found: ${BLAS_LIBRARIES}")
  target_link_libraries     (CoMISo PRIVATE BLAS::BLAS)
else()
  set (COMISO_BLAS_CONFIG_FILE_SETTINGS "#define COMISO_BLAS_AVAILABLE 0" )
  message (STATUS "CoMISo: optional dependency BLAS not found.")
endif ()

find_package (ADOLC QUIET)
if (ADOLC_FOUND)
  set (COMISO_ADOLC_CONFIG_FILE_SETTINGS "#define COMISO_ADOLC_AVAILABLE 1" )
  target_include_directories(CoMISo PUBLIC ${ADOLC_INCLUDE_DIR})
  target_link_directories   (CoMISo PUBLIC ${ADOLC_LIBRARY_DIR})
  target_link_libraries     (CoMISo PUBLIC ${ADOLC_LIBRARIES})
else ()
  set (COMISO_ADOLC_CONFIG_FILE_SETTINGS "#define COMISO_ADOLC_AVAILABLE 0" )
  message (STATUS "CoMISo: optional dependency ADOLC not found.")
endif ()

SET(SUITESPARSE_PACKAGES AMD CAMD CCOLAMD CHOLMOD COLAMD CXSparse GraphBLAS KLU KLU_CHOLMOD LAGraph LDL ParU RBio SPEX SPQR UMFPACK)

SET(SUITESPARSE_PACKAGES_FOUND "")
SET(SUITESPARSE_PACKAGES_NOT_FOUND "")
SET(SUITESPARSE_PACKAGES_DISABLED "")
foreach(_package ${SUITESPARSE_PACKAGES})
    option(COMISO_ENABLE_SUITESPARSE_${_package}   "Build CoMISo with SuiteSparse component ${_package} (if found)" ${COMISO_ENABLE_SUITESPARSE_DEFAULT})
    if(COMISO_ENABLE_SUITESPARSE_${_package})
        find_package(${_package} CONFIG QUIET)
        if (${${_package}_FOUND})
            target_link_libraries(CoMISo PUBLIC SuiteSparse::${_package})
            list(APPEND SUITESPARSE_PACKAGES_FOUND ${_package})
            set(COMISO_SUITESPARSE_${_package}_AVAILABLE TRUE)
        else()
            list(APPEND SUITESPARSE_PACKAGES_NOT_FOUND ${_package})
            set(COMISO_SUITESPARSE_${_package}_AVAILABLE FALSE)
        endif()
    else()
        list(APPEND SUITESPARSE_PACKAGES_DISABLED ${_package})
        set(COMISO_SUITESPARSE_${_package}_AVAILABLE FALSE)
    endif()
endforeach()

if(SUITESPARSE_PACKAGES_FOUND)
    message (STATUS "${PROJECT_NAME}: SuiteSparse packages found: ${SUITESPARSE_PACKAGES_FOUND}")
endif()
if(SUITESPARSE_PACKAGES_NOT_FOUND)
    message (STATUS "${PROJECT_NAME}: SuiteSparse packages not found: ${SUITESPARSE_PACKAGES_NOT_FOUND}")
endif()
if(SUITESPARSE_PACKAGES_DISABLED)
    message (STATUS "${PROJECT_NAME}: SuiteSparse packages disabled: ${SUITESPARSE_PACKAGES_DISABLED}")
endif()

if(COMISO_ENABLE_MPI)
    find_package (MPI COMPONENTS CXX) # https://cmake.org/cmake/help/latest/module/FindMPI.html
    if (MPI_CXX_FOUND )

        set (COMISO_MPI_CONFIG_FILE_SETTINGS "#define COMISO_MPI_AVAILABLE 1" )
        message (WARNING "CoMISo: Building with MPI found via find_package(): ${MPI_CXX_LIBRARIES}")
        target_link_libraries(CoMISo PRIVATE MPI::MPI_CXX)
    else ()
        message (WARNING "CoMISo: COMISO_ENABLE_MPI is set, but MPI not found.")
        set (COMISO_MPI_CONFIG_FILE_SETTINGS "#define COMISO_MPI_AVAILABLE 0" )
    endif ()
else ()
    message (STATUS "CoMISo: building without MPI (enable with COMISO_ENABLE_MPI)!")
    set (COMISO_MPI_CONFIG_FILE_SETTINGS "#define COMISO_MPI_AVAILABLE 0" )
endif ()

if(COMISO_ENABLE_PETSC)
  find_package (PETSC)
  if (PETSC_FOUND AND MPI_FOUND)
    set (COMISO_PETSC_CONFIG_FILE_SETTINGS "#define COMISO_PETSC_AVAILABLE 1" )
    target_include_directories(CoMISo PRIVATE ${PETSC_INCLUDE_DIRS})
    target_link_libraries     (CoMISo PRIVATE ${PETSC_LIBRARY})
  else ()
    message (STATUS "PETSC or MPI not found!")
    set (COMISO_PETSC_CONFIG_FILE_SETTINGS "#define COMISO_PETSC_AVAILABLE 0" )
  endif ()
else ()
  message (STATUS "PETSC disabled!")
  set (COMISO_PETSC_CONFIG_FILE_SETTINGS "#define COMISO_PETSC_AVAILABLE 0" )
endif()


if(COMISO_ENABLE_TAO)
  find_package (TAO)
  if (TAO_FOUND AND PETSC_FOUND AND MPI_FOUND)
    set (COMISO_TAO_CONFIG_FILE_SETTINGS "#define COMISO_TAO_AVAILABLE 1" )
    target_include_directories(CoMISo PRIVATE ${TAO_INCLUDE_DIRS})
    target_link_directories   (CoMISo PRIVATE ${TAO_LIBRARY_DIR})
    target_link_libraries     (CoMISo PRIVATE ${TAO_LIBRARY})
  else ()
    message (STATUS "TAO or PETSC or MPI not found!")
    set (COMISO_TAO_CONFIG_FILE_SETTINGS "#define COMISO_TAO_AVAILABLE 0" )
  endif ()
else ()
  message (STATUS "TAO disabled!")
  set (COMISO_TAO_CONFIG_FILE_SETTINGS "#define COMISO_TAO_AVAILABLE 0" )
endif ()

if (COMISO_ENABLE_METIS)
    if (NOT TARGET METIS::METIS)
        find_package (METIS)
    endif()
    if (TARGET METIS::METIS)
        message (STATUS "CoMISo: Building with METIS")
        set (COMISO_METIS_CONFIG_FILE_SETTINGS "#define COMISO_METIS_AVAILABLE 1" )

        target_link_libraries(CoMISo PRIVATE METIS::METIS)
    else()
        set (COMISO_METIS_CONFIG_FILE_SETTINGS "#define COMISO_METIS_AVAILABLE 0" )
        message (STATUS "CoMISo: COMISO_ENABLE_METIS set, but METIS not found.")
    endif ()
else()
    message (STATUS "CoMISo: Building without METIS")
    set (COMISO_METIS_CONFIG_FILE_SETTINGS "#define COMISO_METIS_AVAILABLE 0" )
endif()

if (COMISO_ENABLE_MUMPS)
    find_package (MUMPS)
    if (MUMPS_FOUND )
        message (STATUS "CoMISo: Building with MUMPS")
        set (COMISO_MUMPS_CONFIG_FILE_SETTINGS "#define COMISO_MUMPS_AVAILABLE 1" )
        target_include_directories(CoMISo PRIVATE ${MUMPS_INCLUDE_DIR})
        target_link_libraries     (CoMISo PRIVATE ${MUMPS_LIBRARY})
    else ()
        message (STATUS "CoMISo: COMISO_ENABLE_MUMPS set, but MUMPS not found!")
        set (COMISO_MUMPS_CONFIG_FILE_SETTINGS "#define COMISO_MUMPS_AVAILABLE 0" )
    endif()
else()
    message (STATUS "CoMISo: Building without MUMPS")
    set (COMISO_MUMPS_CONFIG_FILE_SETTINGS "#define COMISO_MUMPS_AVAILABLE 0" )
endif ()


if(COMISO_ENABLE_IPOPT)
    if (TARGET IPOPT::IPOPT)
          target_link_libraries (CoMISo PUBLIC PkgConfig::IPOPT)
          set (COMISO_IPOPT_CONFIG_FILE_SETTINGS "#define COMISO_IPOPT_AVAILABLE 1" )
          set (COMISO_HSL_CONFIG_FILE_SETTINGS "#define COMISO_HSL_AVAILABLE 0" )
          message(STATUS "CoMISo: linking against provided target IPOPT::IPOPT.")
    elseif (NOT TARGET PkgConfig::IPOPT)
        pkg_check_modules(IPOPT IMPORTED_TARGET ipopt GLOBAL)
    endif()
    if (TARGET PkgConfig::IPOPT)
          target_link_libraries (CoMISo PUBLIC PkgConfig::IPOPT)
          set (COMISO_IPOPT_CONFIG_FILE_SETTINGS "#define COMISO_IPOPT_AVAILABLE 1" )
          set (COMISO_HSL_CONFIG_FILE_SETTINGS "#define COMISO_HSL_AVAILABLE 0" )
          message(STATUS "CoMISo: Linking against IPOPT found via pkgconfig")
    else()
        find_package (IPOPT)

        if (IPOPT_FOUND)
          set (COMISO_IPOPT_CONFIG_FILE_SETTINGS "#define COMISO_IPOPT_AVAILABLE 1" )
          target_include_directories(CoMISo PUBLIC ${IPOPT_INCLUDE_DIRS})
          target_link_directories   (CoMISo PUBLIC ${IPOPT_LIBRARY_DIRS})
          target_link_libraries     (CoMISo PUBLIC ${IPOPT_LIBRARIES})
          if ( IPOPT_HSL_LIBRARY_DIR )
            set (COMISO_HSL_CONFIG_FILE_SETTINGS "#define COMISO_HSL_AVAILABLE 1" )
          else ()
            set (COMISO_HSL_CONFIG_FILE_SETTINGS "#define COMISO_HSL_AVAILABLE 0" )
          endif()
          message(STATUS "CoMISo: Linking against IPOPT found via find_package()")
        else ()
          message (WARNING "CoMISo: IPOPT or dependency not found, but COMISO_ENABLE_IPOPT is set!")
          set (COMISO_IPOPT_CONFIG_FILE_SETTINGS "#define COMISO_IPOPT_AVAILABLE 0" )
          set (COMISO_HSL_CONFIG_FILE_SETTINGS "#define COMISO_HSL_AVAILABLE 0" )
        endif()
    endif()
else()
    set (COMISO_IPOPT_CONFIG_FILE_SETTINGS "#define COMISO_IPOPT_AVAILABLE 0" )
    set (COMISO_HSL_CONFIG_FILE_SETTINGS "#define COMISO_HSL_AVAILABLE 0" )
    message(STATUS "CoMISo: Building without IPOPT. Set COMISO_ENABLE_IPOPT to enable.")
endif()


if (COMISO_ENABLE_DCO)
    find_package (DCO)
    if (DCO_FOUND )
    set (COMISO_DCO_CONFIG_FILE_SETTINGS "#define COMISO_DCO_AVAILABLE 1" )
    target_include_directories   (CoMISo PUBLIC ${DCO_INCLUDE_DIR})
    else ()
    message (STATUS "DCO not found!")
    set (COMISO_DCO_CONFIG_FILE_SETTINGS "#define COMISO_DCO_AVAILABLE 0" )
    endif ()
endif()

if (COMISO_ENABLE_CBC)
    find_package (CBC)
    if (CBC_FOUND )
        set (COMISO_CBC_CONFIG_FILE_SETTINGS "#define COMISO_CBC_AVAILABLE 1" )
        target_include_directories(CoMISo PRIVATE ${CBC_INCLUDE_DIRS})
        target_link_libraries     (CoMISo PRIVATE ${CBC_LIBRARIES})
    else ()
        message (WARNING "CoMISo: Set COMISO_ENABLE_CBC is set, but CBC not found.")
        set (COMISO_CBC_CONFIG_FILE_SETTINGS "#define COMISO_CBC_AVAILABLE 0" )
    endif ()
else()
    message (STATUS "CoMISo: building without optional dependency CBC. Set COMISO_ENABLE_CBC to enable.")
    set (COMISO_CBC_CONFIG_FILE_SETTINGS "#define COMISO_CBC_AVAILABLE 0" )
endif()

if (COMISO_ENABLE_CLP)
    find_package (CLP)
    if (CLP_FOUND )
        set (COMISO_CLP_CONFIG_FILE_SETTINGS "#define COMISO_CLP_AVAILABLE 1" )
        target_include_directories(CoMISo PRIVATE ${CLP_INCLUDE_DIRS})
        target_link_libraries     (CoMISo PRIVATE ${CLP_LIBRARIES})
    else ()
        message (STATUS "CLP not found!")
        set (COMISO_CLP_CONFIG_FILE_SETTINGS "#define COMISO_CLP_AVAILABLE 0" )
    endif ()
endif()

if (COMISO_ENABLE_CGL)
    find_package (CGL)
    if (CGL_FOUND )
        set (COMISO_CGL_CONFIG_FILE_SETTINGS "#define COMISO_CGL_AVAILABLE 1" )
        target_include_directories(CoMISo PRIVATE ${CGL_INCLUDE_DIRS})
        target_link_libraries     (CoMISo PRIVATE ${CGL_LIBRARIES})
    else ()
        message (STATUS "CGL not found!")
        set (COMISO_CGL_CONFIG_FILE_SETTINGS "#define COMISO_CGL_AVAILABLE 0" )
    endif ()
endif()

# XXX TODO: clean this up: get coinutils automatically in finders for coin projects that depend on it
find_package (CoinUtils)
if (COINUTILS_FOUND )
  set (COMISO_COINUTILS_CONFIG_FILE_SETTINGS "#define COMISO_COINUTILS_AVAILABLE 1" )
  target_link_libraries     (CoMISo PRIVATE Coin::CoinUtils)
else ()
  message (STATUS "COINUTILS not found!")
  set (COMISO_COINUTILS_CONFIG_FILE_SETTINGS "#define COMISO_COINUTILS_AVAILABLE 0" )
endif ()

if (COMISO_ENABLE_OSI)
    find_package (OSI)
    if (OSI_FOUND )
        set (COMISO_OSI_CONFIG_FILE_SETTINGS "#define COMISO_OSI_AVAILABLE 1" )
        target_include_directories(CoMISo PRIVATE ${OSI_INCLUDE_DIRS})
        target_link_libraries     (CoMISo PRIVATE ${OSI_LIBRARIES})
    else ()
        message (STATUS "OSI not found!")
        set (COMISO_OSI_CONFIG_FILE_SETTINGS "#define COMISO_OSI_AVAILABLE 0" )
    endif ()
endif()

if (COMISO_ENABLE_TAUCS)
    find_package (Taucs QUIET)
    if (TAUCS_FOUND AND NOT SUPRESS_TAUCS)
    list (APPEND NEED_LAPACK "Taucs")
    set (COMISO_TAUCS_CONFIG_FILE_SETTINGS "#define COMISO_TAUCS_AVAILABLE 1" )
    target_include_directories(CoMISo PRIVATE ${TAUCS_INCLUDE_DIR})
    target_link_libraries     (CoMISo PRIVATE ${TAUCS_LIBRARY})

    else ()
    message (STATUS "TAUCS not found!")
    set (COMISO_TAUCS_CONFIG_FILE_SETTINGS "#define COMISO_TAUCS_AVAILABLE 0" )
    endif ()
else()
    set (COMISO_TAUCS_CONFIG_FILE_SETTINGS "#define COMISO_TAUCS_AVAILABLE 0" )
endif()

if (COMISO_ENABLE_GUROBI)
    if (NOT TARGET Gurobi::GurobiCXX )
        find_package (Gurobi)
    endif()
    if (TARGET Gurobi::GurobiCXX )
        message (STATUS "CoMISo: Building with Gurobi")
        set (COMISO_GUROBI_CONFIG_FILE_SETTINGS "#define COMISO_GUROBI_AVAILABLE 1" )
        target_link_libraries     (CoMISo PRIVATE Gurobi::GurobiC  Gurobi::GurobiCXX)
    else ()
        message (STATUS "CoMISo: COMISO_ENABLE_GUROBI set, but Gurobi not found!")
        set (COMISO_GUROBI_CONFIG_FILE_SETTINGS "#define COMISO_GUROBI_AVAILABLE 0" )
    endif()
else()
    set (COMISO_GUROBI_CONFIG_FILE_SETTINGS "#define COMISO_GUROBI_AVAILABLE 0" )
    message (STATUS "CoMISo: Building without Gurobi")
endif ()

if (COMISO_ENABLE_MOSEK)
    find_package (Mosek QUIET)
    if (MOSEK_FOUND )
        set (COMISO_MOSEK_CONFIG_FILE_SETTINGS "#define COMISO_MOSEK_AVAILABLE 1" )
        target_link_libraries (CoMISo PRIVATE Mosek::FusionCXX)
        message (STATUS "CoMISo: Building with MOSEK")
    else ()
        message (WARNING "CoMISo: COMISO_ENABLE_MOSEK is set, but MOSEK not found")
        set (COMISO_MOSEK_CONFIG_FILE_SETTINGS "#define COMISO_MOSEK_AVAILABLE 0" )
    endif ()
else()
    message (STATUS "CoMISo: Building without MOSEK (enable with COMISO_ENABLE_MOSEK)")
    set (COMISO_MOSEK_CONFIG_FILE_SETTINGS "#define COMISO_MOSEK_AVAILABLE 0" )
endif()


if (COMISO_ENABLE_ARPACK)
    find_package (ARPACK)
    if (ARPACK_FOUND )
        set (COMISO_ARPACK_CONFIG_FILE_SETTINGS "#define COMISO_ARPACK_AVAILABLE 1" )
        target_include_directories(CoMISo PRIVATE ${ARPACK_INCLUDE_DIR})
        target_link_libraries     (CoMISo PRIVATE ${ARPACK_LIBRARY})
        message (STATUS "CoMISo: Building with ARPACK found via find_package")
    else ()
        message (WARNING "CoMISo: COMISO_ENABLE_ARPACK set, but ARPACK not found.")
        set (COMISO_ARPACK_CONFIG_FILE_SETTINGS "#define COMISO_ARPACK_AVAILABLE 0" )
    endif ()
else()
    message (STATUS "CoMISo: Building without ARPACK (enable with COMISO_ENABLE_ARPACK).")
    set (COMISO_ARPACK_CONFIG_FILE_SETTINGS "#define COMISO_ARPACK_AVAILABLE 0" )
endif ()

if (COMISO_ENABLE_CPLEX)
    find_package (CPLEX)
    if (CPLEX_FOUND )
        set (COMISO_CPLEX_CONFIG_FILE_SETTINGS "#define COMISO_CPLEX_AVAILABLE 1" )
        target_include_directories(CoMISo PRIVATE ${CPLEX_INCLUDE_DIRS})
        target_link_libraries     (CoMISo PRIVATE ${CPLEX_LIBRARIES})

        #enable c++ support
        add_definitions(-DIL_STD)
    else ()
        message (WARNING "CoMISo: but COMISO_ENABLE_CPLEX is set, but CPLEX not found")
        set (COMISO_CPLEX_CONFIG_FILE_SETTINGS "#define COMISO_CPLEX_AVAILABLE 0" )
    endif ()
else()
    message (STATUS "CoMISo: Building without CPLEX (enable with COMISO_ENABLE_CPLEX)")
    set (COMISO_CPLEX_CONFIG_FILE_SETTINGS "#define COMISO_CPLEX_AVAILABLE 0" )
endif()

if (COMISO_ENABLE_TINYAD)
    if (TARGET TinyAD AND NOT TARGET tinyad::tinyad)
        add_library(tinyad::tinyad ALIAS TinyAD)
    endif()
    if (NOT TARGET tinyad::tinyad)
        find_package (TINYAD)
    endif()
    if (TARGET tinyad::tinyad)
        message(STATUS "CoMISo: Building CoMISo with TinyAD support")
        set (COMISO_TINYAD_CONFIG_FILE_SETTINGS "#define COMISO_TINYAD_AVAILABLE 1" )
        target_link_libraries(CoMISo PUBLIC tinyad::tinyad)
    else()
        message(STATUS "CoMISo: Building CoMISo without TinyAD support (not found)")
        set (COMISO_TINYAD_CONFIG_FILE_SETTINGS "#define COMISO_TINYAD_AVAILABLE 0" )
    endif ()
else()
    message(STATUS "CoMISo: Building without TinyAD support (enable with COMISO_ENABLE_TINYAD)")
    set (COMISO_TINYAD_CONFIG_FILE_SETTINGS "#define COMISO_TINYAD_AVAILABLE 0" )
endif ()

# XXX TODO: clean this up:
set(TMP_CMAKE_FIND_LIBRARY_PREFIXES "${CMAKE_FIND_LIBRARY_PREFIXES}")
set(CMAKE_FIND_LIBRARY_PREFIXES lib "") #Our lapack library is called liblapack.lib. Is there a better way to find it than this?
if (NEED_LAPACK AND NOT SUITESPARSE_FOUND) # TODO: why does this check for suitesparse? does suitesparse always bring lapack in?
  find_package(LAPACK)
  if (LAPACK_FOUND)
      message(STATUS "CoMISo: Found LAPACK: ${LAPACK_LIBRARIES}")
      target_link_libraries     (CoMISo PRIVATE LAPACK::LAPACK)
  else()
      message(STATUS "CoMISo: LAPACK not found!")
  endif()
endif(NEED_LAPACK AND NOT SUITESPARSE_FOUND)
set(CMAKE_FIND_LIBRARY_PREFIXES "${TMP_CMAKE_FIND_LIBRARY_PREFIXES}")


if(COMISO_ENABLE_LPSOLVE)
    if (NOT TARGET lpsolve::lpsolve )
        find_package (LPSolve)
    endif()
    if (TARGET lpsolve::lpsolve )
        message (STATUS "CoMISo: building with lpsolve")
        set (COMISO_LPSOLVE_CONFIG_FILE_SETTINGS "#define COMISO_LPSOLVE_AVAILABLE 1" )
        target_link_libraries     (CoMISo PUBLIC lpsolve::lpsolve)
    else ()
        message (WARNING "CoMISo: COMISO_ENABLE_LPSOLVE set, but lpsolve not found.")
        set (COMISO_LPSOLVE_CONFIG_FILE_SETTINGS "#define COMISO_LPSOLVE_AVAILABLE 0" )
    endif ()
else()
    message (STATUS "CoMISo: building without lpsolve (enable with COMISO_ENABLE_LPSOLVE)")
    set (COMISO_LPSOLVE_CONFIG_FILE_SETTINGS "#define COMISO_LPSOLVE_AVAILABLE 0" )
endif()

# generate dllexport macros on windows
if (WIN32)
  target_compile_definitions(CoMISo PRIVATE
      -D_SCL_SECURE_NO_DEPRECATE
      -D_USE_MATH_DEFINES
      -DNOMINMAX
      )
endif ()
# ignore "... needs to have dll-interface to be used by clients"
target_compile_options(CoMISo PRIVATE "$<$<CXX_COMPILER_ID:MSVC>:/wd4251>")

# enable debug stuff

if (CMAKE_CXX_COMPILER_ID MATCHES "[cC][lL][aA][nN][gG]")
  # disable unused parameter warning
  #add_definitions(-Wno-unused-parameter -Wno-deprecated-register)
  target_compile_options(CoMISo PRIVATE -Wno-unused-parameter -Wno-deprecated-register)
endif()

if(COMISO_ENABLE_OSQP)
    if (NOT TARGET osqp::osqp AND NOT TARGET osqp::osqpstatic)
        find_package(osqp CONFIG)
    endif()
    if (TARGET osqp::osqp OR TARGET osqp::osqpstatic)
        message("OSQP found")
        if (TARGET osqp::osqp)
            target_link_libraries (CoMISo PUBLIC osqp::osqp)
        else()
            target_link_libraries (CoMISo PUBLIC osqp::osqpstatic)
        endif()
        # I think osqp should link with CMAKE_DL_LIBS, but it doesn't.
        # So we need to do it here.
        #target_link_libraries(CoMISo PRIVATE ${CMAKE_DL_LIBS})
        set (COMISO_OSQP_CONFIG_FILE_SETTINGS "#define COMISO_OSQP_AVAILABLE 1" )
        if (COMISO_OSQP_PRERELEASE_API)
            string(APPEND COMISO_OSQP_CONFIG_FILE_SETTINGS "\n#define COMISO_OSQP_NEW_API 0")
        else()
            string(APPEND COMISO_OSQP_CONFIG_FILE_SETTINGS "\n#define COMISO_OSQP_NEW_API 1")
        endif()

    else()
        message (STATUS "CoMISo: optional dependency OSQP not found.")
        set (COMISO_OSQP_CONFIG_FILE_SETTINGS "#define COMISO_OSQP_AVAILABLE 0" )
    endif()

else()
    message (STATUS "CoMISo: Building without OSQP (enable with COMISO_ENABLE_OSQP)")
    set (COMISO_OSQP_CONFIG_FILE_SETTINGS "#define COMISO_OSQP_AVAILABLE 0" )
endif()

set_target_properties(CoMISo PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

# display results
vci_print_configure_header (COMISO "CoMISo")

# write config file
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/Config/config.hh.in"
     "CoMISo/Config/config.hh" @ONLY IMMEDIATE)

# TODO: split into multiple config files to minimize re-compilation when unrelated dependencies change:
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/Config/config_suitesparse.hh.in"
     "CoMISo/Config/config_suitesparse.hh" @ONLY IMMEDIATE)



#######################################################################
# Configure the examples last to be sure, that all configure files
# of the library are already there
#######################################################################

option (COMISO_BUILD_EXAMPLES "Build CoMISo examples" ${PROJECT_IS_TOP_LEVEL})

if (COMISO_BUILD_EXAMPLES )
    add_subdirectory (Examples/)
endif ()

option(COMISO_NO_INSTALL "avoid cmake install/export (hack)" OFF) # TODO: do we still need this? where?

# Only create install target, when we are building CoMISo standalone
if(NOT COMISO_NO_INSTALL)

  include(CMakePackageConfigHelpers)
  include(GNUInstallDirs)

  set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/CoMISo)

  write_basic_package_version_file(
      CoMISoConfigVersion.cmake
      VERSION ${PROJECT_VERSION}
      COMPATIBILITY SameMajorVersion)

  install(FILES
      "${CMAKE_CURRENT_BINARY_DIR}/CoMISoConfigVersion.cmake"
      DESTINATION ${INSTALL_CONFIGDIR})


  set(INSTALL_TARGETS "CoMISo")
  if(TARGET Mosek::FusionCXX)
      list(APPEND INSTALL_TARGETS "FusionCXX")
  endif()
  if(TARGET Gurobi::GurobiCXX)
      list(APPEND INSTALL_TARGETS "GurobiCXX")
  endif()
  # Install library
  install(TARGETS ${INSTALL_TARGETS}
      EXPORT CoMISoTargets
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      )

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

set(COMISO_ENABLE_UNITTESTS ${PROJECT_IS_TOP_LEVEL} CACHE BOOL "Build CoMISo unit tests.")
if (COMISO_ENABLE_UNITTESTS)
    enable_testing()
    find_package(GTest QUIET)

    if (NOT TARGET GTEST::gtest)
        message(STATUS "CoMISo tests enabled, but GTest not found - downloading using FetchContent.")
        include(FetchContent)
        FetchContent_Declare(googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG        v1.17.x
        )
        if(WIN32)
            # avoid linking errors, cf https://stackoverflow.com/questions/12540970/how-to-make-gtest-build-mdd-instead-of-mtd-by-default-using-cmake
            set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        endif()
        FetchContent_MakeAvailable(googletest)

        # gtest targets have their output directories explicitly set to something different
        # than what we use. Fix it so the tests run on windows (without rpath).
        foreach (_target gtest gtest_main)
            set_target_properties(${_target} PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
                ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
            )
        endforeach()
    endif()
    include(GoogleTest)
    add_subdirectory(Tests)
endif()
