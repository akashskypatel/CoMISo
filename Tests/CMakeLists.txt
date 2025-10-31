add_executable(comiso_tests
    test.cpp
    exact_constraint_projection.cpp
    )
target_link_libraries(comiso_tests PRIVATE
    CoMISo::CoMISo
    GTest::gtest GTest::gtest_main)

gtest_discover_tests(comiso_tests
    EXTRA_ARGS "--gtest_color=yes"
    DISCOVERY_MODE PRE_TEST
)

