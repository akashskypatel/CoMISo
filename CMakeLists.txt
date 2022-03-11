cmake_minimum_required(VERSION 3.13)

project(CoMISo)
add_library(CoMISo)
target_include_directories(CoMISo PUBLIC ${PROJECT_SOURCE_DIR}/..)


# Syntax: ensure_initialised(var1 var2 ...)
#
# Ensure that all variables in the argument list have been initialised. If a
# variable has not been defined, it is set to a blank string.
macro(ensure_initialised)
  foreach(var IN ITEMS ${ARGN})
    if(NOT DEFINED ${var})
      set(${var} "")
    endif()
  endforeach()
endmacro(ensure_initialised)

# Add a library to CoMISo.
# Syntax: comiso_add_lib(<name> <INTERFACE|PUBLIC|PRIVATE> [REQUIRED]  TARGETS [<target>...])
# If one of the targets specified after TARGETS is defined, CoMISo is linked against it
# as INTERFACE, PUBLIC, or PRIVATE library. Also, COMISO_$name_AVAILABLE is set to 1.
# If non of the targets is found, COMISO_$name_AVAILABLE is set to 0.
# If non of the targets is found and additionally REQUIRED is specified, the call will fail.
function(comiso_add_lib name visibility)
  set(options REQUIRED)
  set(singleValueArgs "")
  set(multiValueArgs TARGETS)
  cmake_parse_arguments(COMISO_ADD_LIB "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  # initialize variable to 0, later set to 1 if a target is found
  set(COMISO_${name}_AVAILABLE 0)

  # Check all provided targets
  foreach(target IN ITEMS ${COMISO_ADD_LIB_TARGETS})
    if (TARGET ${target})
      set(COMISO_${name}_AVAILABLE 1)
      target_link_libraries(CoMISo ${visibility} ${target})
      break() # dont try to find the other target names for the same target
    endif()
  endforeach()

  # Check if REQUIRED library was found
  if (COMISO_ADD_LIB_REQUIRED AND COMISO_${name}_AVAILABLE EQUAL 0)
    message(FATAL_ERROR "Required target ${name} not available")
  endif()

  # Make COMISO_${name}_AVAILABLE available in parent scope
  set(COMISO_${name}_AVAILABLE ${COMISO_${name}_AVAILABLE} PARENT_SCOPE)
endfunction(comiso_add_lib)

comiso_add_lib(IPOPT PRIVATE TARGETS Ipopt)
comiso_add_lib(OSQP PRIVATE TARGETS OSQP)
comiso_add_lib(EIGEN PUBLIC REQUIRED TARGETS Eigen3::Eigen Eigen eigen)
comiso_add_lib(BASE PUBLIC REQUIRED TARGETS Base)


# Function to add subdirectories and group source files
function(comiso_add_subdir subdir)
  set(full_subdir "${PROJECT_SOURCE_DIR}/${subdir}")
  add_subdirectory(${full_subdir})

  #various ways of organizing these files are possible, by default split them in different sub-groups cmk, hdr, imp, src
  source_group(${subdir}/cmk FILES ${subdir}/CMakeLists.txt)

  ensure_initialised(my_headers my_t_impls my_sources)
  source_group(${subdir}/hdr FILES ${my_headers})
  source_group(${subdir}/imp FILES ${my_t_impls})
  source_group(${subdir}/src FILES ${my_sources})

  # but each developer can customize this locally as they see fit, e.g., put all in a single subfolder (below)
  # source_group(${subdir} FILES ${subdir}/CMakeLists.txt ${my_headers} ${my_t_impls} ${my_sources})

  # exclude template implementaiton (T.cc) files from the build
  set_source_files_properties(${my_t_impls} PROPERTIES HEADER_FILE_ONLY TRUE)

  ensure_initialised(HEADERS SOURCES)
  set(HEADERS ${HEADERS} ${subdir}/CMakeLists.txt ${my_headers} ${my_t_impls} PARENT_SCOPE)
  set(SOURCES ${SOURCES} ${my_sources} PARENT_SCOPE)
endfunction(comiso_add_subdir)

comiso_add_subdir(Config) # for config.hh.in make sure to add all libs above before adding Config subdir
comiso_add_subdir(EigenSolver)
comiso_add_subdir(NSolver)
comiso_add_subdir(Solver)
comiso_add_subdir(Utils)

target_sources(CoMISo PRIVATE ${SOURCES} ${HEADERS})

# generate dllexport macros on windows
if (MSVC)
	
	#add_definitions(-DCOMISODLL) #build static for now
	#add_definitions(-DINCLUDE_TEMPLATES) (disabled, use explicit instantiation instead)
	add_definitions(-D_SCL_SECURE_NO_DEPRECATE)
	add_definitions(-D_USE_MATH_DEFINES)
	add_definitions(-DNOMINMAX)
	# add_definitions(-DIL_STD) #set CPLEX to use its STL interface, change this?!

	add_definitions(/W3) #add appropriate warnings flags for this project, the compilation raises to too may warnings
	add_definitions(/MP) #build on all cores 
	
endif (MSVC)

