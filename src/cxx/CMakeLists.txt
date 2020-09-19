find_package(Boost)
if(BOOST_FOUND)
    message("Boost Found")
    message(BOOST_VERSION)
endif()

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

find_package(blaze)
if(blaze_FOUND)
    add_library(blaze_target INTERFACE)
    target_link_libraries(blaze_target INTERFACE blaze::blaze)
endif()

find_package (Eigen3 3.3 NO_MODULE)

add_library(mpc_mosek_eigen SHARED mpc_mosek_eigen.h mpc_mosek_eigen.cpp mpc_mosek_utility_eigen.h)
add_library(mosek_eigen_blaze SHARED main.cpp)
target_include_directories(mpc_mosek_eigen PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
if(TARGET Eigen3::Eigen)
    target_link_libraries(mpc_mosek_eigen Eigen3::Eigen)
else()
    target_include_directories(mpc_mosek_eigen PUBLIC "${EIGEN3_INCLUDE_DIR}")
endif()
target_include_directories(mosek_eigen_blaze PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(mpc_mosek_eigen mosek)
target_link_libraries(mosek_eigen_blaze mpc_mosek_eigen blaze_target)
message("before exporting")
message("MEB_INCLUDE_DIR is ${MEB_INCLUDE_DIR}")
export(TARGETS mosek_eigen_blaze mpc_mosek_eigen blaze_target FILE mosek_eigen_blaze-targets.cmake)
message("MEB_INCLUDE_DIR is ${MEB_INCLUDE_DIR}")
install(TARGETS mosek_eigen_blaze mpc_mosek_eigen blaze_target
        EXPORT mosek_eigen_blaze-targets
        INCLUDES DESTINATION ${MEB_INCLUDE_DIR}
        LIBRARY DESTINATION ${MEB_LIBRARY_DIR}
        PUBLIC_HEADER DESTINATION ${MEB_INCLUDE_DIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
message("after exporting")
configure_file(mosek_eigen_blaze.pc.in mosek_eigen_blaze.pc @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/mosek_eigen_blaze.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
configure_package_config_file(mosek_eigen_blaze-config.cmake.in mosek_eigen_blaze-config.cmake
        INSTALL_DESTINATION ${MEB_DIR}
        PATH_VARS MEB_INCLUDE_DIR MEB_LIBRARY_DIR
        NO_CHECK_REQUIRED_COMPONENTS_MACRO)

install(EXPORT mosek_eigen_blaze-targets DESTINATION ${MEB_DIR})

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/mosek_eigen_blaze-config.cmake"
        DESTINATION ${MEB_DIR})


add_executable(test_meb main.cpp)

target_link_libraries(test_meb ${Boost_LIBRARIES} Eigen3::Eigen mpc_mosek_eigen blaze_target)
#target_link_libraries(test_blaze )
