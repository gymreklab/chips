file(GLOB SOURCES *.cpp *.h *.hpp)
add_library(ChIPs ${SOURCES})
target_link_libraries(ChIPs ${htslib_static} ${zlib_static})
target_include_directories(ChIPs PUBLIC "${PROJECT_BINARY_DIR}")
