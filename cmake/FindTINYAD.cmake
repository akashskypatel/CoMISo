set(TINYAD_DIR $ENV{TINYAD_DIR} CACHE PATH "TinyAD include folder (contains TinyAD/Scalar.h).")
find_path(TINYAD_INCLUDE_DIR
           NAMES TinyAD/Scalar.hh
           PATHS ${TINYAD_DIR}
                 ${TINYAD_DIR}/include
                 /usr/include
                 /usr/local/include
          )
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TINYAD  DEFAULT_MSG
                                  TINYAD_INCLUDE_DIR)
mark_as_advanced(TINYAD_INCLUDE_DIR)

if (TINYAD_FOUND AND NOT TARGET tinyad::tinyad)
    add_library(tinyad::tinyad INTERFACE IMPORTED)
    target_include_directories(tinyad::tinyad INTERFACE ${TINYAD_INCLUDE_DIR})
endif()

