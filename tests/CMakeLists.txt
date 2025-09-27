# Example test CMake config (adjust when using real test framework)
add_executable(test_dummy tests/test_core.cpp)
add_test(NAME DummyTest COMMAND test_dummy)
