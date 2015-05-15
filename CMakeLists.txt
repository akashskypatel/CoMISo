cmake_minimum_required(VERSION 3.0)

project(CoMISo)
	
my_add_subdir(EigenSolver)
my_add_subdir(Config)
my_add_subdir(NSolver)
my_add_subdir(Solver)
my_add_subdir(Utils)

include_directories(${PROJECT_SOURCE_DIR}/..)
include_directories(${3P_DIR}/Eigen/${INC_DIR})
include_directories(${3P_DIR}/Gmm/${INC_DIR})
# include_directories(${3P_DIR}/CPLEX/${INC_DIR}) # CPLEX is no longer used
include_directories(${3P_DIR}/Gurobi/${INC_DIR})
include_directories(${3P_DIR}/IPOPT/${INC_DIR})
include_directories(${3P_DIR}/Cbc/${INC_DIR})
include_directories(${3P_DIR}/libcurl/${INC_DIR})
include_directories(${3P_DIR}/Boost/${INC_DIR})

# message(STATUS "${SOURCES}" "${HEADES}")
add_library(CoMISo ${SOURCES} ${HEADERS})
target_link_libraries(CoMISo Base)

# generate dllexport macros on windows
if (MSVC)
	
	#add_definitions(-DCOMISODLL) #build static for now
	#add_definitions(-DINCLUDE_TEMPLATES) (disabled, use explicit instantiation instead)
	add_definitions(-D_SCL_SECURE_NO_DEPRECATE)
	add_definitions(-D_USE_MATH_DEFINES)
	add_definitions(-DNOMINMAX)
	# add_definitions(-DIL_STD) #set CPLEX to use its STL interface, change this?!

	add_definitions(/W1) #add appropriate warnings flags for this project, the compilation raises to too may warnings
	add_definitions(/MP) #build on all cores 
	
endif (MSVC)

