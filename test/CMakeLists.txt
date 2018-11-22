find_package(Boost REQUIRED)

file(GLOB test_pathes  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

foreach(test_path ${test_pathes})
    get_filename_component(test_name ${test_path} NAME_WE)
    add_executable(test_${test_name} ${test_path})
    add_test(NAME ${test_name} COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target test_${test_name} )
endforeach()
