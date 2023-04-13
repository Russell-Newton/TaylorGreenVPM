add_library(vpm "")
target_sources(vpm PRIVATE vpm.cpp)
target_sources(vpm PRIVATE treecode.cpp)
target_include_directories(vpm PUBLIC ../include)
