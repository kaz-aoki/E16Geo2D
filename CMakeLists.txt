cmake_minimum_required(VERSION 2.8)

project ("E16Geo2D")

set(CMAKE_VERBOSE_MAKEFILE TRUE)

file(GLOB libE16GeoSRCs RELATIVE ${CMAKE_SOURCE_DIR} src/*.C src/*.cc)

#file GLOB does not work for rootcint. The command thinks a single giant file name instead of a set of file names.
#I cannot avoid so I decided to just list.
#file(GLOB libE16GeoHs RELATIVE ${CMAKE_SOURCE_DIR} include/*.h )

set (libE16GeoHs E16Geo2DMan.h LinkDef.h )

message("DEBUG ===== SRC files")
message("${libE16GeoSRCs}")

#string (REPLACE ";" " " libE16GeoHs "${libE16GeoHs}")
#set(TMP ${libE16GeoHs})

add_custom_command (OUTPUT Dict.C
		  COMMAND rootcint -f ${CMAKE_CURRENT_BINARY_DIR}/Dict.C -c ${libE16GeoHs}
		  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/include"
)
#DEPENDS should be added to add_custom_command. However, this requires "include/xxx.h" while I have is "xxx.h"
#think about it later.

list(APPEND libE16GeoSRCs Dict.C)
message("${libE16GeoSRCs}")

include_directories(${ROOTSYS}/include ${CMAKE_SOURCE_DIR}/include)
link_directories(${ROOTSYS}/lib64)

execute_process(COMMAND root-config --cflags OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE cflags)
execute_process(COMMAND root-config --glibs OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE linkflags)
#execute_process(COMMAND root-config --libdir OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE linkflags)

add_library(E16Geo2D SHARED ${libE16GeoSRCs})
#target_link_options(E16Geo2D PRIVATE -L${linkflags} -lPhysics -lGeom -lHist)
target_link_options(E16Geo2D PRIVATE ${linkflags})

set ( CMAKE_INSTALL_PREFIX "../bin" )
set ( CMAKE_CXX_FLAGS "${cflags} -DE16_ROOT6 -I/e16/w/r01/aoki/E16/geant/geant4.10.06.p03-install/include/Geant4/")
set ( CMAKE_EXE_LINKER_FLAGS ${linkflags} )


