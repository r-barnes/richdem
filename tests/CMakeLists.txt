cmake_minimum_required (VERSION 3.9)

project (richdem_tests)

add_executable(richdem_tests.exe
  fsm_tests.cpp
  test_main.cpp
  tests.cpp
)

target_link_libraries(richdem_tests.exe PRIVATE richdem)
target_compile_definitions(richdem_tests.exe PRIVATE RICHDEM_NO_PROGRESS)
