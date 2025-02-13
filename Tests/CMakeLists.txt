add_executable(comiso_tests test.cpp)
target_link_libraries(comiso_tests PRIVATE Catch2::Catch2WithMain CoMISo::CoMISo)
catch_discover_tests(comiso_tests)
