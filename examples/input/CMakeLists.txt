# find_package(nlohmann_json 3.2.0 REQUIRED)

# add_library(input_ex ...)

add_executable(input_ex input_ex.cpp)
#target_link_libraries(input_ex PRIVATE nlohmann_json::nlohmann_json)