######################################################################
## CMakeLists.txt --- gsIncompressibleFlow
## This file is part of the G+Smo library.
##
######################################################################

## gsIncompressibleFlow module
project(gsIncompressibleFlow)

# Apply G+Smo config
include(gsConfig)
	
## Collect files
aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_H  )
aux_cpp_directory        (${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_CPP)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_HPP)

if( (NOT GISMO_BUILD_LIB) )
  aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_INS)
  if(${PROJECT_NAME}_INS)
    LIST( REMOVE_ITEM ${PROJECT_NAME}_CPP ${${PROJECT_NAME}_INS})
  endif()
endif()

set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} "${PROJECT_SOURCE_DIR}/src"
  CACHE INTERNAL "Gismo include directories")
  
include_directories(${GISMO_INCLUDE_DIRS})

# Add object library
add_library(${PROJECT_NAME} OBJECT
  ${${PROJECT_NAME}_H}
  ${${PROJECT_NAME}_HPP}
  ${${PROJECT_NAME}_CPP} )

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  FOLDER "G+Smo modules" )

set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "G+Smo modules" )

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/src"
  DESTINATION include/gismo
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

# add filedata folder
add_definitions(-DFLOW_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/filedata/")

# add example files
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    target_link_libraries(${tarname} gismo)
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})
